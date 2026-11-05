import {screenPathCandidates} from './sanger-path-screen';
import {JunctionInput,JunctionFit,JunctionModel,JunctionReference,junctionStarts,moleculeJunctionStarts} from './sanger-junction-mixture';
import {controlKnownContextEmission,controlContextProbability} from './sanger-control-response';
import {fitEmissionMixture,emissionVector} from './sanger-emission';
import {mixPhase,phaseNeighbor} from './sanger-signal';
export interface PathOrigin {start:number;weight:number;probabilities?:number[][]}
export interface PathEvent {target:number;size:number;origins:PathOrigin[];profiled?:boolean}
const bases='ACGT',copy=(v:number[][]|undefined)=>v?.map(a=>[...a]);
export const pathKey=(m:{size:number;sizes?:number[]})=>m.sizes?.join('/')??String(m.size);
export const isWildType=(m:{size:number;sizes?:number[]})=>m.sizes?m.sizes.every(v=>v===0):m.size===0;
/** Ordered physical events on the ORIGINAL reference. Every later measured
 * onset includes the accumulated sizes of all preceding events. */
export function pathBase(sequence:string,events:{start:number;size:number;probabilities?:number[][]}[],p:number):number[]{
 let delta=0;for(const e of events){const start=e.start+delta;if(p<start)break;if(e.size>0&&p<start+e.size)return e.probabilities?.[p-start]||[.25,.25,.25,.25];delta+=e.size;}
 const b=sequence[p-delta];return Array.from(bases,x=>b==='N'?.25:Number(x===b));
}
function combinations(events:PathEvent[]){
 let rows:{events:{target:number;start:number;size:number;probabilities?:number[][]}[];weight:number}[]=[{events:[],weight:1}];
 for(const e of events){const next:typeof rows=[];for(const a of rows)for(const o of e.origins){if(o.weight<=1e-9)continue;const previous=a.events.at(-1);if(previous&&o.start<previous.start+Math.max(0,-previous.size))continue;next.push({events:[...a.events,{target:e.target,start:o.start,size:e.size,probabilities:o.probabilities}],weight:a.weight*o.weight});}rows=next;}
 const total=rows.reduce((s,x)=>s+x.weight,0);return rows.map(x=>({...x,weight:x.weight/Math.max(1e-12,total)}));
}
export function pathColumn(ref:JunctionReference,m:JunctionModel,positions:number[],lag:number,theta:number[]):number[]{
 const choices=combinations(m.path||[]),out=Array(positions.length*4).fill(0);
 for(const choice of choices){
  if(positions.length<48){
   const one=(p:number)=>{const pi=[-1,0,1].map(d=>pathBase(ref.sequence,choice.events,p+d));if(!ref.contextResponse)return pi[1];const hard=pi.map(v=>v.findIndex(x=>x===1));return hard.every(v=>v>=0)?controlKnownContextEmission(bases[hard[0]],bases[hard[1]],bases[hard[2]],ref.contextResponse,p):controlContextProbability(pi[0],pi[1],pi[2],ref.contextResponse,p);};
   positions.forEach((p,j)=>{const v=lag?mixPhase(one(p),one(phaseNeighbor(p,lag)),lag):one(p);for(let b=0;b<4;b++)out[4*j+b]+=choice.weight*v[b];});continue;
  }
  const parts:string[]=[],soft=new Map<number,number[]>();let at=0,delta=0;
  for(const event of choice.events){
   parts.push(ref.sequence.slice(at,event.start));
   if(event.size>0){let inserted='';for(let j=0;j<event.size;j++){const pi=event.probabilities?.[j]||[.25,.25,.25,.25],hard=pi.findIndex(v=>v===1);inserted+=hard<0?'N':bases[hard];if(hard<0)soft.set(event.start+delta+j,pi);}parts.push(inserted);}
   at=event.start+Math.max(0,-event.size);delta+=event.size;
  }
  parts.push(ref.sequence.slice(at));const sequence=parts.join('');
  const probability=(p:number)=>soft.get(p)||Array.from(bases,b=>Number(sequence[p]===b));
  const one=(p:number)=>{
   if(!ref.contextResponse)return probability(p);
   const left=sequence[p-1]||'',current=sequence[p]||'',right=sequence[p+1]||'';
   return left!=='N'&&current!=='N'&&right!=='N'?controlKnownContextEmission(left,current,right,ref.contextResponse,p):controlContextProbability(probability(p-1),probability(p),probability(p+1),ref.contextResponse,p);
  };
  positions.forEach((p,j)=>{const a=one(p),v=lag?mixPhase(a,one(phaseNeighbor(p,lag)),lag):a;for(let b=0;b<4;b++)out[4*j+b]+=choice.weight*v[b];});
 }
 return out;
}
/** Multi-target composition of the calibrated junction model, using its same
 * normalized-mixture regression and conditional origin/base updates. Candidate
 * columns are generated in bounded working sets; no observations are masked. */

/** Map the SAME physical edited molecule to another primer/read axis. Sizes
 * retain their target labels; reverse reads reverse/complement inserted bases. */
export function mapPathModel(m:JunctionModel,map:{reverse:boolean;offset:number}):JunctionModel {
 const paths=(m.path||[]).map(e=>({...e,origins:e.origins.map(o=>({...o,start:map.reverse?map.offset-o.start-Math.max(0,-e.size):map.offset+o.start,probabilities:o.probabilities?(map.reverse?[...o.probabilities].reverse().map(v=>[...v].reverse()):copy(o.probabilities)):undefined}))}));
 return {...m,path:map.reverse?paths.reverse():paths};
}
const readRows=new WeakMap<JunctionInput,Map<number,{read:number;p:number}>>();
function axes(input:JunctionInput){let a=readRows.get(input);if(!a){a=new Map;input.reads!.forEach((r,read)=>r.rowIds.forEach((id,j)=>a!.set(id,{read,p:r.positions[j]})));readRows.set(input,a);}return a;}
export function jointPathColumn(input:JunctionInput,m:JunctionModel,rows=input.positions):number[]{
 if(!input.reads)return pathColumn(input.references[m.reference],m,rows,input.lag,input.theta);
 const axis=axes(input),out=Array(rows.length*4).fill(0);
 input.reads.forEach((read,index)=>{const ids=rows.flatMap((id,j)=>axis.get(id)?.read===index?[j]:[]),positions=ids.map(j=>axis.get(rows[j])!.p),mapped=mapPathModel(m,read.mappings[m.reference]),values=pathColumn(read.references[m.reference],mapped,positions,read.lag,read.theta);ids.forEach((row,k)=>{for(let b=0;b<4;b++)out[4*row+b]=values[4*k+b];});});return out;
}
function originRows(input:JunctionInput,m:JunctionModel,target:number,size:number,starts:number[],delta:number):number[]{
 if(!input.reads){const lo=starts[0]+delta-2,hi=starts.at(-1)!+delta+Math.max(0,size)+3;return input.positions.flatMap((p,j)=>p>=lo&&p<hi?[j]:[]);}
 const axis=axes(input);return input.positions.flatMap((id,j)=>{const {read:index,p}=axis.get(id)!,read=input.reads![index],map=read.mappings[m.reference],cuts=read.references[m.reference].cuts!,cut=cuts[target-1],shift=cuts.reduce((s,c,t)=>s+(c<cut?(m.sizes?.[t]||0):0),0),a=map.reverse?map.offset-starts.at(-1)!-Math.max(0,-size):map.offset+starts[0],b=map.reverse?map.offset-starts[0]-Math.max(0,-size):map.offset+starts.at(-1)!;return p>=a+shift-2&&p<b+shift+Math.max(0,size)+3?[j]:[];});
}
function refinePathSet(input:JunctionInput,seeds:JunctionModel[]=[]):JunctionFit {
 const {references:refs,positions,observed,positionWeights:precision}=input,n=observed.length;
 const cuts=refs.map(r=>r.cuts||[r.cut]),range=Array.from({length:input.maxDeletion+input.maxInsertion+1},(_,i)=>i-input.maxDeletion);
 const models:JunctionModel[]=[],columns:number[][]=[],known=new Set<string>();let currentWeights:number[]=[];const candidateCount=cuts.reduce((sum,c)=>sum+Math.pow(range.length,c.length),0);
 const make=(reference:number,sizes:number[]):JunctionModel=>({reference,size:sizes.reduce((s,v)=>s+v,0),sizes:[...sizes],start:cuts[reference][0],inserted:'',path:sizes.flatMap((size,t)=>{
  if(!size)return [];const starts=moleculeJunctionStarts(refs[reference].sequence,cuts[reference][t],size,input.radius),centre=cuts[reference][t]-Math.floor(Math.max(0,-size)/2),start=starts.reduce((a,b)=>Math.abs(b-centre)<Math.abs(a-centre)?b:a,starts[0]);
  const donor=[...seeds,...models].map((m,i)=>({m,w:i<seeds.length?1:(currentWeights[i-seeds.length]||0)})).filter(x=>x.m.reference===reference&&x.w>1e-8&&x.m.path?.some(e=>e.target===t+1&&e.size===size)).sort((a,b)=>b.w-a.w)[0]?.m.path?.find(e=>e.target===t+1&&e.size===size);
  return [{target:t+1,size,profiled:donor?.profiled,origins:donor?donor.origins.map(o=>({...o,probabilities:copy(o.probabilities)})):start===undefined?[]:[{start,weight:1,probabilities:size>0?Array.from({length:size},()=>[.25,.25,.25,.25]):undefined}]}];
 }).sort((a,b)=>cuts[reference][a.target-1]-cuts[reference][b.target-1])});
 const key=(m:JunctionModel)=>m.reference+':'+pathKey(m),column=(m:JunctionModel,rows=positions)=>jointPathColumn(input,m,rows);
 const add=(m:JunctionModel,c?:number[])=>{if(known.has(key(m)))return;known.add(key(m));models.push(m);columns.push(c||column(m));};
 for(const c of input.classes||[])add(make(c.reference,c.sizes||[c.size]));
 let fit=input.fit(columns,observed,precision),iterations=0,searchComplete=!!input.classes;
 const loss=(cols:number[][],y:number[],w:number[],q?:number[])=>{let error=0;for(let p=0;p<y.length/4;p++){const v=[0,0,0,0];cols.forEach((c,i)=>{if(w[i])for(let b=0;b<4;b++)v[b]+=w[i]*c[4*p+b];});const total=v.reduce((s,x)=>s+x,0);for(let b=0;b<4;b++)error+=(q?.[p]??1)*(y[4*p+b]-(total?v[b]/total:0))**2;}return error;};
 let objective=loss(columns,observed,fit.weights,precision);const objectiveHistory=[objective];
 // A working-set cap controls memory, not allele fraction. Incomplete searches
 // return numerical estimates and an explicit convergence diagnostic.
 for(let round=0;round<12;round++){
  iterations++;currentWeights=fit.weights;
  const sum=Array(n).fill(0);columns.forEach((c,i)=>{if(fit.weights[i])for(let k=0;k<n;k++)sum[k]+=fit.weights[i]*c[k];});
  for(let id=0;id<models.length;id++){
   const mass=fit.weights[id],m=models[id];if(mass<1e-8||!m.path?.length)continue;
   let delta=0;
   for(const target of cuts[m.reference].map((_,t)=>t).sort((a,b)=>cuts[m.reference][a]-cuts[m.reference][b])){
    const event=m.path.find(e=>e.target===target+1),size=m.sizes![target];if(!event){delta+=size;continue;}
    const starts=moleculeJunctionStarts(refs[m.reference].sequence,cuts[m.reference][target],size,input.radius);if(!starts.length)continue;
    const ids=originRows(input,m,target+1,size,starts,delta);if(!ids.length){delta+=size;continue;}
    const rows=ids.map(j=>positions[j]),y=ids.flatMap(j=>observed.slice(4*j,4*j+4)),q=ids.map(j=>precision?.[j]??1),old=ids.flatMap(j=>columns[id].slice(4*j,4*j+4)),other=ids.flatMap(j=>sum.slice(4*j,4*j+4)).map((v,k)=>v-mass*old[k]);
    const seed=event.origins.reduce((a,b)=>a.weight>=b.weight?a:b),origins=starts.map(start=>{const previous=event.origins.find(o=>o.start===start);return {start,weight:previous?.weight||0,probabilities:copy(previous?.probabilities||seed.probabilities)};});
    const conditional=(origin:PathOrigin)=>{const saved=event.origins;event.origins=[{...origin,weight:1}];const c=column(m,rows);event.origins=saved;return c.map((v,k)=>other[k]+mass*v);};
    const variants=origins.map(conditional);
    const refine=(o:number,current:number[],fraction:number,sweeps:number)=>{
     if(fraction<=1e-8||size<=0)return;
     for(let sweep=0;sweep<sweeps;sweep++){
      let change=0;
      for(let k=0;k<size;k++){
       const basis=Array.from({length:4},(_,base)=>conditional({...origins[o],probabilities:origins[o].probabilities!.map((v,j)=>j===k?Array.from(bases,(_,d)=>Number(base===d)):v)}).map((v,j)=>current[j]+fraction*(v-variants[o][j]))),previous=origins[o].probabilities![k],estimate=fitEmissionMixture(basis,y,q);
       if(loss(basis,y,estimate.weights,q)>loss(basis,y,previous,q)+1e-12)continue;
       change=Math.max(change,...estimate.weights.map((v,b)=>Math.abs(v-previous[b])));origins[o].probabilities![k]=estimate.weights;
       const next=conditional(origins[o]);for(let j=0;j<current.length;j++)current[j]+=fraction*(next[j]-variants[o][j]);variants[o]=next;
      }
      if(change<1e-6)break;
     }
    };
    // As in the single-target engine, profile all unknown insertion origins
    // once, then refine only origins carrying conditional molecular mass.
    if(size>0&&!event.profiled){for(let o=0;o<origins.length;o++)refine(o,[...variants[o]],1,5);event.profiled=true;}
    const candidate=fitEmissionMixture(variants,y,q),before=loss([sum],observed,[1],precision);
    if(size>0){const current=y.map((_,j)=>variants.reduce((s,v,o)=>s+candidate.weights[o]*v[j],0));for(let o=0;o<origins.length;o++)refine(o,current,candidate.weights[o],3);}
    const saved=event.origins;event.origins=origins.map((o,i)=>({...o,weight:candidate.weights[i]})).filter(o=>o.weight>1e-8);const next=column(m),after=loss([sum.map((v,k)=>v+mass*(next[k]-columns[id][k]))],observed,[1],precision);
    if(after<=before+1e-12){for(let k=0;k<n;k++)sum[k]+=mass*(next[k]-columns[id][k]);columns[id]=next;}else event.origins=saved;
    delta+=size;
   }
  }
  const next=input.fit(columns,observed,precision);if(loss(columns,observed,next.weights,precision)<=loss(columns,observed,fit.weights,precision)+1e-12)fit=next;
  const value=loss(columns,observed,fit.weights,precision);objectiveHistory.push(value);
  if(objective-value<1e-8){objective=value;break;}objective=value;
 }
 const predictColumns=(rows:number[])=>models.map(m=>column(m,rows));
 return {...fit,models,columns,iterations,duplicateLengths:0,objective,objectiveHistory,predictColumns,searchComplete,candidateCount};
}

/** Sparse working-set search over the same physical forward model. A proposed
 * compound edit is profiled while it can carry mass, before dense single-edit
 * columns absorb its two observable shift marginals. BIC chooses model size;
 * the retained molecular fractions are then unpenalized estimates. */
export function fitJunctionPaths(input:JunctionInput):JunctionFit {
 if(input.classes)return refinePathSet(input);
 const cuts=input.references.map(r=>r.cuts||[r.cut]),range=Array.from({length:input.maxDeletion+input.maxInsertion+1},(_,i)=>i-input.maxDeletion),count=cuts.reduce((n,c)=>n+range.length**c.length,0);
 const classes=(models:JunctionModel[])=>models.map(m=>({reference:m.reference,size:m.size,sizes:m.sizes}));
 const key=(m:JunctionModel)=>m.reference+':'+pathKey(m);
 const wt=input.references.map((r,h)=>({reference:h,size:0,sizes:cuts[h].map(()=>0)}));
 let current=refinePathSet({...input,classes:wt}),iterations=current.iterations,complete=false;
 const loss=(fit:JunctionFit)=>{let sse=0;for(let p=0;p<input.positions.length;p++){const v=[0,0,0,0];fit.columns.forEach((c,i)=>{if(fit.weights[i])for(let b=0;b<4;b++)v[b]+=fit.weights[i]*c[4*p+b];});const sum=v.reduce((s,v)=>s+v,0);for(let b=0;b<4;b++)sse+=(input.positionWeights?.[p]??1)*(input.observed[4*p+b]-(sum?v[b]/sum:0))**2;}return sse;};
 const bic=(fit:JunctionFit)=>{let df=Math.max(0,fit.weights.filter(w=>w>1e-8).length-1);fit.models.forEach((m,i)=>{if(fit.weights[i]<=1e-8)return;for(const e of m.path||[]){const active=e.origins.filter(o=>o.weight>1e-6);df+=Math.max(0,active.length-1)+active.length*3*Math.max(0,e.size);}});const n=3*input.positions.length;return n*Math.log(Math.max(1e-15,loss(fit)/n))+df*Math.log(n);};
 let cost=bic(current);const history=[loss(current)],startCache=new Map<string,number>();
 const initialStart=(h:number,t:number,size:number)=>{const key=h+':'+t+':'+size;if(startCache.has(key))return startCache.get(key)!;const starts=moleculeJunctionStarts(input.references[h].sequence,cuts[h][t],size,input.radius),centre=cuts[h][t]-Math.floor(Math.max(0,-size)/2),start=starts.reduce((a,b)=>Math.abs(b-centre)<Math.abs(a-centre)?b:a,starts[0]);startCache.set(key,start);return start;};
 const make=(h:number,sizes:number[]):JunctionModel=>({reference:h,size:sizes.reduce((s,v)=>s+v,0),sizes:[...sizes],start:Math.min(...cuts[h]),inserted:'',path:sizes.flatMap((size,t)=>{if(!size)return[];const start=initialStart(h,t,size);return[{target:t+1,size,origins:start===undefined?[]:[{start,weight:1,probabilities:size>0?Array.from({length:size},()=>[.25,.25,.25,.25]):undefined}]}];}).sort((a,b)=>cuts[h][a.target-1]-cuts[h][b.target-1])});
 for(let round=0;round<24;round++){
  const ranked=screenPathCandidates(input,current,range,make);
  let best=current,bestCost=cost;
  for(const candidate of ranked){
   const trial=refinePathSet({...input,classes:[...classes(current.models),{reference:candidate.m.reference,size:candidate.m.size,sizes:candidate.m.sizes}]},current.models),value=bic(trial);iterations+=trial.iterations;
   if(value<bestCost-1e-6){best=trial;bestCost=value;}
  }
  if(best===current){complete=true;break;}
  current=best;cost=bestCost;history.push(loss(current));
  // Remove redundant physical groups after every birth, refitting both origins
  // and fractions. WT remains available for every submitted reference.
  for(let sweep=0;sweep<3;sweep++){
   let reduced=current,reducedCost=cost;
   for(let i=0;i<current.models.length;i++){if(isWildType(current.models[i]))continue;const trial=refinePathSet({...input,classes:classes(current.models.filter((_,j)=>i!==j))},current.models),value=bic(trial);iterations+=trial.iterations;if(value<reducedCost-1e-6){reduced=trial;reducedCost=value;}}
   if(reduced===current)break;current=reduced;cost=reducedCost;
  }
 }
 return{...current,iterations,searchComplete:complete,candidateCount:count,objectiveHistory:history};
}

export function pathAnnotations(m:JunctionModel,fraction:number,origin:number,length:number,reverse:boolean){
 return combinations(m.path||[]).map(c=>({size:m.size,pathSizes:m.sizes,fraction:fraction*c.weight,events:c.events.map(e=>{let start=e.start+origin,end=start+Math.max(0,-e.size);if(reverse)[start,end]=[length-end,length-start];return {target:e.target,start,end,inserted:Math.max(0,e.size),insertedBaseProbabilities:e.probabilities?(reverse?[...e.probabilities].reverse().map(v=>[...v].reverse()):copy(e.probabilities)):undefined};})}));
}
