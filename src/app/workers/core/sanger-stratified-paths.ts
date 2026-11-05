import {decodeFrequencyPaths} from './sanger-frequency-path-decoder';
import type {Prepared,PathEdit} from './sanger-path-common';
/** Enumerate the bounded allocation of allele slots across measured net-length
 * classes. Each allocation gets its own search budget; a cheap pseudo-trace
 * hypothesis cannot starve a different class allocation before raw scoring. */
export function lengthClassAllocations(sizes:number[],slots=6){
 const types=[...new Set(sizes)].sort((a,b)=>a-b),out:Record<string,number>[]=[];
 const visit=(at:number,left:number,counts:Record<string,number>)=>{if(at===types.length){if(!left)out.push({...counts});return;}const remaining=types.length-at-1;for(let n=1;n<=left-remaining;n++){counts[String(types[at])]=n;visit(at+1,left-n,counts);}};visit(0,slots,{});return out;
}
export function proposeStratifiedPaths(ref:Prepared,families:{reference:number;size:number;fraction:number}[],positions:number[],observed:number[],noise:number,maxVariants=5,statesPerAllocation=200,profilesPerAllocation=4,progress:(stage:string)=>void=()=>{}){
 const allocations=lengthClassAllocations(families.filter(f=>f.fraction>0).map(f=>f.size),maxVariants+1),options:{events:PathEdit[];fraction:number}[][]=[],seen=new Set<string>(),details:any[]=[],scores:{rmse:number;maximum:number}[]=[];
 for(const [i,counts]of allocations.entries()){
  progress(`Comparing length-class allocations · ${i+1} / ${allocations.length}`);
  try{const decoded=decodeFrequencyPaths([ref],families,positions,observed,noise,3,10,maxVariants,statesPerAllocation,()=>{},true,counts),best=decoded.candidates.sort((a,b)=>a.rmse-b.rmse).slice(0,profilesPerAllocation);details.push({counts,...decoded.details[0],profiles:best.length});for(const c of best){const key=JSON.stringify(c.groups.map(g=>g.events).sort((a,b)=>JSON.stringify(a).localeCompare(JSON.stringify(b))));if(!seen.has(key)){seen.add(key);options.push(c.groups);scores.push({rmse:c.rmse,maximum:c.maximum});}}}catch(error){details.push({counts,profiles:0,error:String(error)});}
 }
 return{options,details,scores};
}
