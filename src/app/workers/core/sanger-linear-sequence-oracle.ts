import {materialize,PathEdit} from './sanger-path-common';
import {normalizeMoleculeEdits} from './sanger-edit-normalize';
export interface LinearSequenceOracleInput {reference:string;cuts:number[];cost:number[][];radius:number;maxDeletion:number;maxInsertion:number;maxSubstitution?:number}
/** Exact linear minimization over one I/D/S tract per separated target.
 * Cost is additive per emitted base. Context-dependent yields require a
 * higher-order emission oracle and are deliberately outside this contract.
 * Dynamic programming retains one optimum per cumulative read shift. Inserted
 * strings are optimized base by base, so no 4^length dictionary is enumerated. */
export function linearSequenceOracle(input:LinearSequenceOracleInput):{sequence:string;events:PathEdit[];cost:number}{
 const {reference,cuts,cost,radius,maxDeletion,maxInsertion}=input,n=reference.length,sub=input.maxSubstitution??10,bases='ACGT',span=Math.max(maxDeletion,sub);
 if(!cuts.length||cuts.some((c,i)=>c<0||c>n||(i>0&&c-radius<=Math.min(n,cuts[i-1]+radius+span))))throw Error('Linear sequence oracle requires separated target windows.');
 const at=(p:number,b:string)=>cost[p]?.[bases.indexOf(b)]??0,minimum=cost.map(row=>Math.min(...row)),calls=cost.map(row=>bases[row.indexOf(Math.min(...row))]),prefixes=new Map<number,Float64Array>();
 const prefix=(shift:number)=>{let v=prefixes.get(shift);if(!v){v=new Float64Array(n+1);for(let r=0;r<n;r++)v[r+1]=v[r]+at(r+shift,reference[r]);prefixes.set(shift,v);}return v;};
 const word=(p:number,length:number)=>Array.from({length},(_,j)=>calls[p+j]??'A').join(''),wordCost=(p:number,length:number)=>{let sum=0;for(let j=0;j<length;j++)sum+=minimum[p+j]??0;return sum;};
 interface State {cost:number;events:PathEdit[]}
 let states=new Map<number,State>([[0,{cost:0,events:[]}]]),anchor=0;
 for(const cut of cuts){const nextAnchor=Math.min(n,cut+radius+span),next=new Map<number,State>();
  const keep=(shift:number,value:number,events:PathEdit[])=>{const old=next.get(shift);if(!old||value<old.cost||(value===old.cost&&events.length<old.events.length))next.set(shift,{cost:value,events});};
  for(const [shift,state]of states){const pre=prefix(shift);keep(shift,state.cost+pre[nextAnchor]-pre[anchor],state.events);
   for(let start=Math.max(anchor,cut-radius);start<=Math.min(n,cut+radius);start++){
    const before=state.cost+pre[start]-pre[anchor];
    const add=(end:number,length:number)=>{const newShift=shift+length-(end-start),post=prefix(newShift),inserted=word(start+shift,length),value=before+wordCost(start+shift,length)+post[nextAnchor]-post[end];keep(newShift,value,[...state.events,{start,end,inserted}]);};
    for(let length=1;length<=maxDeletion&&start+length<=nextAnchor;length++)add(start+length,0);
    for(let length=1;length<=maxInsertion;length++)add(start,length);
    for(let length=1;length<=sub&&start+length<=nextAnchor;length++)add(start+length,length);
   }
  }
  states=next;anchor=nextAnchor;
 }
 let best:{cost:number;events:PathEdit[]}|undefined;
 for(const [shift,state]of states){const pre=prefix(shift),value=state.cost+pre[n]-pre[anchor];if(!best||value<best.cost||(value===best.cost&&state.events.length<best.events.length))best={cost:value,events:state.events};}
 const events=normalizeMoleculeEdits(reference,best!.events);return {cost:best!.cost,events,sequence:materialize(reference,events)!};
}
