import {materialize,PathEdit} from './sanger-path-common';
import {SangerSequenceFitCache} from './sanger-sequence-fit-cache';
import {normalizeMoleculeEdits} from './sanger-edit-normalize';
const bases='ACGT';
/** Bounded two-cycle phase exchanges with variable projection. Candidate
 * hypotheses are generated against measured fluorescence with free target
 * pairing coefficients within one reference. Other references are held fixed
 * only during proposal generation; the caller must score returned profiles
 * with the full joint reference fit before accepting them. Only improvements exceeding the incumbent solver dual-gap
 * bound are committed. Search budgets are not identifiability certificates;
 * this function neither reads ground truth nor introduces ambiguous bases. */
export function refineConditionalPhase(refs:{sequence:string;cuts:number[]}[],prefix:PathEdit[][][],initial:PathEdit[][][],positions:number[],observed:number[],progress:(s:string)=>void=()=>{},targetIndex=1){
 const profiles=initial.map(gs=>{const nonWT=gs.filter(es=>es.length).map(es=>es.map(e=>({...e}))),all:PathEdit[][]=[[],...nonWT];while(all.length<6)all.push([]);return all;}),models:{h:number;first:number;second:number;sequence:string}[]=[],sequences:string[]=[],keys=new Map<string,number>(),positionIndex=new Map(positions.map((p,i)=>[p,i]));
 let cache:SangerSequenceFitCache|undefined;
 let previousWeights:number[]|undefined;
 let current:any,previousFit:{ids:number[];weights:number[]}|undefined;
 const blocks=new Map<number,{key:string;cache:SangerSequenceFitCache}>();
 const fit=(scope?:number)=>{const ids:number[]=[];for(let h=0;h<refs.length;h++)for(let i=0;i<prefix[h].length;i++)for(let j=0;j<profiles[h].length;j++){
  const sequence=materialize(refs[h].sequence,[...prefix[h][i],...profiles[h][j]])!,key=h+':'+sequence;let id=keys.get(key);if(id===undefined){id=models.length;keys.set(key,id);models.push({h,first:i,second:j,sequence});sequences.push(sequence);}if(!ids.includes(id))ids.push(id);
 }cache??=new SangerSequenceFitCache(sequences,positions,observed);
 if(scope===undefined){const old=new Map(previousFit?.ids.map((id,i)=>[id,previousFit!.weights[i]])||[]),result=cache.fit(ids,ids.map(id=>old.get(id)||0));previousFit={ids,weights:result.weights};return{...result,ids};}
 const outside=current.ids.flatMap((id:number,i:number)=>models[id].h!==scope?[{id,w:current.weights[i]}]:[]),key=outside.map((x:any)=>x.id+':'+x.w).join('/'),share=1-outside.reduce((s:number,x:any)=>s+x.w,0),prediction=cache.predict(outside.map((x:any)=>x.id),outside.map((x:any)=>x.w)),target=observed.map((v,i)=>(v-prediction[i])/share),own=ids.filter(id=>models[id].h===scope);let block=blocks.get(scope);if(!block||block.key!==key){block={key,cache:new SangerSequenceFitCache(sequences,positions,target)};blocks.set(scope,block);}const result=block.cache.fit(own),weights=ids.map(id=>{const at=own.indexOf(id);return at>=0?result.weights[at]*share:outside.find((x:any)=>x.id===id)?.w||0;});return{...result,weights,rmse:result.rmse*share,optimalityGap:(result.optimalityGap||0)*share*share,ids};};
 current=fit();let changes=0;const history=[current.rmse];
 for(let pass=0;pass<3;pass++){
  let changed=false;
  // Keep a bounded beam of coupled phase exchanges. Two individually worse
  // exchanges can jointly restore a consistent multi-base genotype.
  for(let h=0;h<refs.length;h++){
   const neighbors=(state:PathEdit[][])=>{
    const out:PathEdit[][][]=[],seen=new Set<string>(),cut=refs[h].cuts[targetIndex],size=(es:PathEdit[])=>es.reduce((n,e)=>n+e.inserted.length-e.end+e.start,0),mutable=(es:PathEdit[])=>size(es)<0?[]:size(es)>0?Array.from({length:size(es)},(_,i)=>es[0].start+i):Array.from({length:16},(_,i)=>cut-3+i);
    const add=(state:PathEdit[][])=>{const sorted=[[],...state.slice(1).sort((a,b)=>JSON.stringify(a).localeCompare(JSON.stringify(b)))],key=JSON.stringify(sorted);if(!seen.has(key)){seen.add(key);out.push(sorted);}};
    const change=(es:PathEdit[],q:number,b:string)=>{const n=size(es),seq=materialize(refs[h].sequence,es)!;if(n>0){const start=es[0].start,word=[...seq.slice(start,start+n)];word[q-start]=b;return[{start,end:start,inserted:word.join('')}];}const start=Math.min(q,...es.map(e=>e.start)),end=Math.max(q+1,...es.map(e=>e.end)),word=[...seq.slice(start,end)];word[q-start]=b;const events=normalizeMoleculeEdits(refs[h].sequence,[{start,end,inserted:word.join('')}]);return events.length&&(events[0].start<cut-3||events[0].start>cut+3||events.at(-1)!.end-events[0].start>10)?null:events;};
    for(let a=1;a<state.length;a++){
     const old=state[a],n=size(old),seq=materialize(refs[h].sequence,old)!;
     for(const q of mutable(old))for(const base of bases){if(base===seq[q])continue;const es=change(old,q,base);if(es){const next=[...state];next[a]=es;add(next);}}
     if(old.length)for(let start=cut-3;start<=cut+3;start++){const es=normalizeMoleculeEdits(refs[h].sequence,[{start,end:start+(n<0?-n:n===0?10:0),inserted:n<0?'':seq.slice(start,start+(n>0?n:10))}]),next=[...state];next[a]=es;add(next);}

    }
    // Exchange a whole phase tract, not just one base. Inserted words and
    // substitutions may require several coordinated base assignments before
    // the jointly fitted signal improves.
    const changeTract=(es:PathEdit[],qs:number[],other:string)=>{const n=size(es),seq=materialize(refs[h].sequence,es)!;
     if(n>0){const start=es[0].start,word=[...seq.slice(start,start+n)];for(const q of qs)word[q-start]=other[q];return[{start,end:start,inserted:word.join('')}];}
     const start=Math.min(...qs,...es.map(e=>e.start)),end=Math.max(...qs.map(q=>q+1),...es.map(e=>e.end)),word=[...seq.slice(start,end)];for(const q of qs)word[q-start]=other[q];const events=normalizeMoleculeEdits(refs[h].sequence,[{start,end,inserted:word.join('')}]);return events.length&&(events[0].start<cut-3||events[0].start>cut+3||events.at(-1)!.end-events[0].start>10)?null:events;
    };
    for(let a=1;a<state.length;a++)for(let b=a+1;b<state.length;b++){
     const sa=materialize(refs[h].sequence,state[a])!,sb=materialize(refs[h].sequence,state[b])!,different=mutable(state[a]).filter(q=>mutable(state[b]).includes(q)&&sa[q]!==sb[q]);
     for(let lo=0;lo<different.length;lo++)for(let hi=lo;hi<different.length;hi++){const qs=different.slice(lo,hi+1),na=changeTract(state[a],qs,sb),nb=changeTract(state[b],qs,sa);if(!na||!nb)continue;const next=[...state];next[a]=na;next[b]=nb;add(next);}
    }
    return out;
   };
   const root=profiles[h],first:{state:PathEdit[][];fit:any}[]=[];for(const state of neighbors(root)){profiles[h]=state;first.push({state,fit:fit(h)});}first.sort((a,b)=>a.fit.rmse-b.fit.rmse);let best={state:root,fit:current};for(const [candidate,parent]of first.slice(0,16).entries()){progress(`Sequence phase · pass ${pass+1} · reference ${h+1} / ${refs.length} · ${candidate+1} / ${Math.min(16,first.length)}`);for(const state of [parent.state,...neighbors(parent.state)]){profiles[h]=state;const trial=fit(h);if(trial.rmse*trial.rmse<best.fit.rmse*best.fit.rmse-Math.max(1e-15,2*(best.fit.optimalityGap||0)))best={state,fit:trial};}}
   profiles[h]=best.state;if(best.state!==root){current=best.fit;changes++;changed=true;history.push(current.rmse);}
  }
  progress('Joint sequence phasing · pass '+(pass+1));
  if(!changed)break;
 }
 return{profiles,rmse:current.rmse,changes,history};
}
