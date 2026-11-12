import {Prepared,PathEdit,materialize} from './sanger-path-common';
/** Enumerate the bounded allocation of allele slots across measured net-length
 * classes. Each allocation gets its own search budget; a cheap pseudo-trace
 * hypothesis cannot starve a different class allocation before raw scoring. */
export function lengthClassAllocations(sizes:number[],slots=6){
 const types=[...new Set(sizes)].sort((a,b)=>a-b),out:Record<string,number>[]=[];
 const visit=(at:number,left:number,counts:Record<string,number>)=>{if(at===types.length){if(!left)out.push({...counts});return;}const remaining=types.length-at-1;for(let n=1;n<=left-remaining;n++){counts[String(types[at])]=n;visit(at+1,left-n,counts);}};visit(0,slots,{});return out;
}
export {adaptivePaths} from './sanger-adaptive-paths';

import{fitSparseFrequencyProfile}from'./sanger-sparse-profile';
import type{SangerLinearOptimizer}from'./sanger-discrete-phase';
/** Every input and every allele count uses the same joint sparse profile fit.
 * No frequency-subset enumeration or fixed allele-count limit is used here. */
export function frequencyProfiles(ref:Prepared,families:{size:number;fraction:number}[],positions:number[],observed:number[],noise:number,budget:number,progress:(s:string)=>void,thresholdPercent:number,optimizer?:SangerLinearOptimizer){
 return fitSparseFrequencyProfile(ref,families,positions,observed,noise,budget,progress,thresholdPercent,optimizer);
}
