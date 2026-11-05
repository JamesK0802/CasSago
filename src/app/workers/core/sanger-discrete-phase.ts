import {materialize,PathEdit} from './sanger-path-common';
import {normalizeMoleculeEdits} from './sanger-edit-normalize';
import type{RawSequenceScorer}from './sanger-raw-interval-fit';
import {SangerSequenceFitCache,sequenceFitLoss,sequenceFitTolerance} from './sanger-sequence-fit-cache';
const bases='ACGT';
export interface SangerPersistentOptimizer{options:{set(values:Record<string,number|boolean|string>):unknown};passModel(model:unknown):unknown;setBasis(basis:unknown):unknown;getBasis():unknown;run():{modelStatus:number};getSolution():{colValue:ArrayLike<number>};dispose():void}
export interface SangerLinearOptimizer {solve(program:string,options:Record<string,number|boolean|string>):{Status:string;Columns:Record<string,{Primal:number}>};createModel?:()=>SangerPersistentOptimizer}
/** Variable projection alternates whole-molecule mixture fitting with legal
 * binary allele sequences. Validated integer-ADC data use measured rounding
 * intervals; other inputs retain centroid simplex fitting. A geometric search
 * over group frequencies escapes base/rate stationary points. Equal-length
 * transfers preserve invariant-flank marginals; near-coincident cross-length
 * transfers are bounded proposals requiring a complete free refit. Every proposal is scored
 * against the complete measured signal without truth-derived candidates. */
export function refineDiscretePhase(refs:{sequence:string;cuts:number[]}[],prefix:PathEdit[][][],initial:PathEdit[][][],positions:number[],observed:number[],optimizer:SangerLinearOptimizer,progress:(s:string)=>void=()=>{},target=1,exchange=true,rawScorer?:RawSequenceScorer){
 const profiles=initial.map(gs=>gs.map(es=>es.map(e=>({...e})))),sequences:string[]=[],keys=new Map<string,number>(),models:{h:number;first:number;second:number}[]=[],cache=new SangerSequenceFitCache(sequences,positions,observed,512,rawScorer);
 let previous:{ids:number[];weights:number[]}|undefined;
 const fit=()=>{const ids:number[]=[],where=new Map<number,{h:number;first:number;second:number}>();for(let h=0;h<refs.length;h++)for(let i=0;i<prefix[h].length;i++)for(let j=0;j<profiles[h].length;j++){const sequence=materialize(refs[h].sequence,[...prefix[h][i],...profiles[h][j]])!,key=h+':'+sequence;let id=keys.get(key);if(id===undefined){id=sequences.length;keys.set(key,id);sequences.push(sequence);models.push({h,first:i,second:j});}if(!where.has(id))where.set(id,{h,first:i,second:j});if(!ids.includes(id))ids.push(id);}const old=new Map(previous?.ids.map((id,i)=>[id,previous!.weights[i]])||[]),result=cache.fit(ids,ids.map(id=>old.get(id)||0));previous={ids,weights:result.weights};return{...result,ids,where};};
 let current=fit(),changes=0;const history=[current.rmse];
 const active=(h:number)=>{if(!current.rawActivePositions)return true;const shifts=target===0?[0]:prefix[h].map(es=>es.reduce((n,e)=>n+e.inserted.length-e.end+e.start,0)),cut=refs[h].cuts[target],lo=cut-3+Math.min(...shifts),hi=cut+13+Math.max(...shifts);return current.rawActivePositions.some(p=>p>=lo&&p<=hi);};
 const referenceOrder=()=>refs.map((_,h)=>{const shifts=target===0?[0]:prefix[h].map(es=>es.reduce((n,e)=>n+e.inserted.length-e.end+e.start,0)),cut=refs[h].cuts[target],lo=cut-3+Math.min(...shifts),hi=cut+13+Math.max(...shifts);return{h,error:(current.rawResiduals||[]).reduce((n,r)=>n+(r.position>=lo&&r.position<=hi?r.error:0),0)};}).sort((a,b)=>b.error-a.error||a.h-b.h).map(r=>r.h);
 for(let pass=0;pass<6;pass++){
  if(current.selectionLoss!==undefined&&current.selectionLoss<1e-7)break;
  let changed=false;
  for(const h of referenceOrder()){
   if(!active(h))continue;
   progress(`Joint base assignment · pass ${pass+1} · reference ${h+1} / ${refs.length}`);
   const root=profiles[h],groups=root.map((es,j)=>({events:es,size:es.reduce((s,e)=>s+e.inserted.length-e.end+e.start,0),sequence:materialize(refs[h].sequence,es)!,editable:j!==root.findIndex(es=>es.length===0),weights:new Map<number,number>()}));
   current.ids.forEach((id,k)=>{const m=current.where.get(id)!;if(m.h!==h)return;const shift=target===0?0:prefix[h][m.first].reduce((s,e)=>s+e.inserted.length-e.end+e.start,0),g=groups[m.second];g.weights.set(shift,(g.weights.get(shift)||0)+current.weights[k]);});
   const cut=refs[h].cuts[target],q=Array.from({length:16},(_,i)=>cut-3+i),shifts=groups.flatMap(g=>[...g.weights.keys()]),lo=q[0]+Math.min(...shifts),hi=q.at(-1)!+Math.max(...shifts),rows=positions.flatMap((p,i)=>p>=lo&&p<=hi?[i]:[]),outside=current.ids.flatMap((id,i)=>models[id].h!==h?[{id,w:current.weights[i]}]:[]),prediction=cache.predict(outside.map(x=>x.id),outside.map(x=>x.w)),y=rows.map(i=>observed.slice(4*i,4*i+4).map((v,b)=>v-prediction[4*i+b]));
   const totalBounds=current.rawResponse?rawScorer?.boundsFor?.(rows.map(i=>positions[i]),current.rawResponse):undefined,bounds=totalBounds?{scale:current.rawResponse!.amplitude,lower:totalBounds.lower.map((r,i)=>r.map((v,b)=>v-prediction[4*rows[i]+b])),upper:totalBounds.upper.map((r,i)=>r.map((v,b)=>v-prediction[4*rows[i]+b]))}:undefined;
   const candidate=assignSequenceBases(refs[h].sequence,cut,groups,rows.map(i=>positions[i]),y,optimizer,bounds);
   if(candidate){profiles[h]=candidate;const trial=fit();if(trial.converged&&sequenceFitLoss(trial)<sequenceFitLoss(current)-sequenceFitTolerance(current)){current=trial;history.push(current.rmse);changes++;changed=true;continue;}profiles[h]=root;}
   if(rawScorer&&bounds&&current.rawResponse){
    // Escape a frequency/base coordinate stationary point without changing any
    // measured length marginal. Transfer mass only between equal-length groups,
    // preserving each upstream shift total; then discover all bases jointly.
    progress(`Measured frequency continuation · target ${target+1} · reference ${h+1} / ${refs.length}`);
    const step=.125/current.rawResponse.amplitude,rate=(j:number)=>[...groups[j].weights.values()].reduce((s,v)=>s+v,0);
    frequencySearch:for(let a=0;a<groups.length;a++)for(let b=a+1;b<groups.length;b++)if(groups[a].size===groups[b].size||Math.abs(rate(a)-rate(b))<=4/current.rawResponse.amplitude){
     for(const [donor,receiver]of [[a,b],[b,a]]){const total=rate(donor);if(total<=0)continue;
      // Near-coincident weights of different length classes can exchange
      // base assignments too. Small cross-class perturbations are proposals;
      // only the complete free coefficient refit can accept them. No invariant
      // length marginal is claimed for this temporary cross-class proposal.
      const maxStep=groups[a].size===groups[b].size?total/2:Math.min(total/2,4/current.rawResponse.amplitude);
      for(let delta=step;delta<=maxStep;delta*=2){
       const altered=groups.map(g=>({...g,weights:new Map(g.weights)}));for(const [shift,w]of groups[donor].weights){const move=delta*w/total;altered[donor].weights.set(shift,w-move);altered[receiver].weights.set(shift,(altered[receiver].weights.get(shift)||0)+move);}
       const proposal=assignSequenceBases(refs[h].sequence,cut,altered,rows.map(i=>positions[i]),y,optimizer,bounds,1);if(!proposal)continue;profiles[h]=proposal;const trial=fit();
       if(trial.converged&&sequenceFitLoss(trial)<sequenceFitLoss(current)-sequenceFitTolerance(current)){current=trial;history.push(current.rmse);changes++;changed=true;break frequencySearch;}profiles[h]=root;
      }
     }
    }
   }
  }
  if(rawScorer&&changed)continue;
  if(exchange){
   // Deletions have no mutable bases. Move two legal deletion origins together
   // and refit the whole mixture, so a fixed-weight phase step cannot trap them.
   for(const h of referenceOrder()){
   if(!active(h))continue;
    const root=profiles[h],cut=refs[h].cuts[target],seq=refs[h].sequence;
    const alternatives=root.map(es=>es.length===1&&es[0].inserted===''?
     Array.from({length:7},(_,i)=>{const start=cut-3+i;return[{start,end:start+es[0].end-es[0].start,inserted:''}];})
      .filter(es=>es[0].start>=0&&es[0].end<=seq.length):[es]);
    let best=current,bestState=root;progress(`Testing deletion origins · target ${target+1} · reference ${h+1} / ${refs.length}`);
    for(let a=0;a<root.length;a++)if(alternatives[a].length>1)
     for(let b=a+1;b<root.length;b++)if(alternatives[b].length>1){
      for(const left of alternatives[a])for(const right of alternatives[b]){
       const state=[...root];state[a]=left;state[b]=right;profiles[h]=state;
       const trial=fit();
       if(trial.converged&&sequenceFitLoss(trial)<sequenceFitLoss(best)-sequenceFitTolerance(best)){best=trial;bestState=state;}
      }
     }
    profiles[h]=bestState;
    if(bestState!==root){current=best;history.push(current.rmse);changes++;changed=true;}
   }
   // A phase exchange is scored after refitting all mixture coefficients.
   // This crosses coordinate-stationary points without a multi-level beam.
   for(const h of referenceOrder()){
   if(!active(h))continue;
    const root=profiles[h],seq=refs[h].sequence,cut=refs[h].cuts[target],size=(es:PathEdit[])=>es.reduce((n,e)=>n+e.inserted.length-e.end+e.start,0),mutable=(es:PathEdit[])=>size(es)<0?[]:size(es)>0?Array.from({length:size(es)},(_,i)=>es[0].start+i):Array.from({length:16},(_,i)=>cut-3+i),words=root.map(es=>materialize(seq,es)!),bestRoot=current;
    let best=bestRoot,bestState=root;progress(`Testing coupled base tracts · target ${target+1} · reference ${h+1} / ${refs.length}`);
    const replace=(es:PathEdit[],qs:number[],other:string)=>{const n=size(es),word=materialize(seq,es)!;if(n>0){const start=es[0].start,inserted=[...word.slice(start,start+n)];for(const q of qs)inserted[q-start]=other[q];return[{start,end:start,inserted:inserted.join('')}];}const start=Math.min(...qs,...es.map(e=>e.start)),end=Math.max(...qs.map(q=>q+1),...es.map(e=>e.end)),inserted=[...word.slice(start,end)];for(const q of qs)inserted[q-start]=other[q];const out=normalizeMoleculeEdits(seq,[{start,end,inserted:inserted.join('')}]);if(out.length&&Math.max(cut-3,out.at(-1)!.end-10)>Math.min(cut+3,out[0].start))return null;return out;};
    for(let a=0;a<root.length;a++)if(root[a].length)for(let b=a+1;b<root.length;b++)if(root[b].length){const qs=mutable(root[a]).filter(q=>mutable(root[b]).includes(q)&&words[a][q]!==words[b][q]);for(let lo=0;lo<qs.length;lo++)for(let hi=lo;hi<qs.length;hi++){const tract=qs.slice(lo,hi+1),left=replace(root[a],tract,words[b]),right=replace(root[b],tract,words[a]);if(!left||!right)continue;const state=[...root];state[a]=left;state[b]=right;profiles[h]=state;const trial=fit();if(trial.converged&&sequenceFitLoss(trial)<sequenceFitLoss(best)-sequenceFitTolerance(best)){best=trial;bestState=state;}}}
    profiles[h]=bestState;if(bestState!==root){current=best;history.push(current.rmse);changes++;changed=true;}
   }
  }
  if(!changed)break;
 }
 return{profiles,changes,rmse:current.rmse,history};
}
interface Group {events:PathEdit[];size:number;sequence:string;weights:Map<number,number>;editable?:boolean}
export function assignSequenceBases(sequence:string,cut:number,groups:Group[],positions:number[],observed:number[][],optimizer:SangerLinearOptimizer,intervals?:{lower:number[][];upper:number[][];scale?:number},proposalTimeLimit=10){
 const q=Array.from({length:16},(_,i)=>cut-3+i),index=new Map(q.map((p,i)=>[p,i])),lines=['Minimize'],constraints:string[]=[],binary:string[]=[],cost:string[]=[],checks:{terms:[number,string][];op:string;rhs:number}[]=[];let row=0;
 const term=(coefficient:number,name:string)=>(coefficient<0?'- ':'+ ')+Math.abs(coefficient)+' '+name,add=(terms:[number,string][],op:string,rhs:number)=>{checks.push({terms,op,rhs});constraints.push(' r'+row+++': '+terms.filter(([c])=>c!==0).map(([c,n])=>term(c,n)).join(' ')+' '+op+' '+rhs);},bit=(g:number,k:number,b:number)=>`b${g}_${k}_${b}`,origin=(g:number,i:number)=>`z${g}_${i}`;
 groups.forEach((g,j)=>{
  if(!(g.editable??g.events.length>0))return;const origins=Array.from({length:7},(_,i)=>origin(j,i));binary.push(...origins);add(origins.map(n=>[1,n]),'=',1);
  q.forEach((p,k)=>{const bits=Array.from({length:4},(_,b)=>bit(j,k,b));binary.push(...bits);add(bits.map(n=>[1,n]),'=',1);
   for(let i=0;i<7;i++){const start=cut-3+i;if(start<=p&&p<start+(g.size>0?g.size:g.size===0?10:0))continue;const at=p<start?p:p-g.size,b=bases.indexOf(sequence[at]);if(b<0)throw Error('Incomplete reference coverage in discrete phase.');add([[1,origins[i]],[-1,bits[b]]],'<=',0);}
  });
 });
 positions.forEach((p,i)=>{for(let b=0;b<4;b++){
  const terms:[number,string][]=[],error=`e${i}_${b}`;let fixed=0;
  groups.forEach((g,j)=>{for(const [shift,w]of g.weights){const at=p-shift,k=index.get(at);if((g.editable??g.events.length>0)&&k!==undefined)terms.push([w,bit(j,k,b)]);else if((g.events.length?sequence[at<q[0]?at:at-g.size]:g.sequence[at])===bases[b])fixed+=w;}});
  const scale=intervals?.scale||1,scaled=terms.map(([c,n])=>[c*scale,n]as[number,string]);cost.push(error);add([...scaled,[-1,error]],'<=',((intervals?.upper[i][b]??observed[i][b])-fixed)*scale);add([...scaled,[1,error]],'>=',((intervals?.lower[i][b]??observed[i][b])-fixed)*scale);
 }});
 if(!binary.length)return null;lines.push(' obj: '+cost.join(' + '),'Subject To',...constraints,'Binary',...binary,'End');
 const answer=optimizer.solve(lines.join('\n'),{output_flag:false,time_limit:proposalTimeLimit,mip_rel_gap:1e-7,mip_abs_gap:1e-9,mip_feasibility_tolerance:1e-9,primal_feasibility_tolerance:1e-9});
 const read=(name:string)=>answer.Columns[name]?.Primal||0;
 if(answer.Status!=='Optimal'){
  if(!/time/i.test(answer.Status)||binary.some(name=>!Number.isFinite(answer.Columns[name]?.Primal)||Math.abs(read(name)-Math.round(read(name)))>1e-7))return null;
  for(const c of checks){const value=c.terms.reduce((sum,[a,n])=>sum+a*read(n),0),error=c.op==='='?Math.abs(value-c.rhs):c.op==='<='?value-c.rhs:c.rhs-value;if(error>1e-7)return null;}
 }

 return groups.map((g,j)=>{if(!(g.editable??g.events.length>0))return[];let choice=0;for(let i=1;i<7;i++)if(read(origin(j,i))>read(origin(j,choice)))choice=i;const start=cut-3+choice,length=g.size>0?g.size:g.size===0?10:0,word=Array.from({length},(_,i)=>{const k=index.get(start+i)!;let dye=0;for(let b=1;b<4;b++)if(read(bit(j,k,b))>read(bit(j,k,dye)))dye=b;return bases[dye];}).join('');return normalizeMoleculeEdits(sequence,[{start,end:start+(g.size<0?-g.size:g.size===0?10:0),inserted:word}]);});
}
