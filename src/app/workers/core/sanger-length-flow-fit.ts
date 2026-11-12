import {sparseLengthFlowBasis} from './sanger-length-flow-basis';
import {fitGroupedIntensities} from './sanger-grouped-mixture';
import type {SimplexFit} from './sanger-signal';
export const lengthFlowDiagnostics:{observe?:(event:Record<string,unknown>)=>void}={};
export interface LengthFlowWindow {classes:{reference:number;size:number;column:number[]}[];observed:number[]}
/** Fits one molecule distribution across all invariant
 * flanks. Column pricing uses additive gradients and a shortest path through
 * cumulative length states; every transition obeys the per-target edit limits.
 * Edited reference shares are fitted jointly, independently of control brightness.
 * The same procedure handles one or multiple windows, without enumerating all
 * molecule paths. No genotype, reporting threshold or allele count is used. */
export function fitLengthFlows(windows:LengthFlowWindow[],shares:number[],deletion:number,insertion:number,robustDelta=Infinity){
 const observed=windows.flatMap(w=>w.observed),offsets=windows.map((_,t)=>windows.slice(0,t).reduce((n,w)=>n+w.observed.length,0)),paths:{reference:number;sizes:number[];column:number[]}[]=[],keys=new Set<string>();
 const add=(reference:number,sizes:number[])=>{const key=reference+':'+sizes.join(',');if(keys.has(key))return false;keys.add(key);paths.push({reference,sizes,column:windows.flatMap((w,t)=>w.classes.find(c=>c.reference===reference&&c.size===sizes[t])!.column)});return true;};
 shares.forEach((_,h)=>add(h,windows.map(()=>0)));
 let fit:SimplexFit,optimalityGap=Infinity,warm:number[]|undefined;
 for(;;){
  lengthFlowDiagnostics.observe?.({phase:'fit',paths:paths.length,delta:robustDelta,windows:windows.length});
  fit=fitGroupedIntensities(paths.map(p=>p.column),observed,paths.map(()=>0),[1],robustDelta,warm);if(!fit.converged)throw Error('Joint invariant-flank path fit did not converge: gap='+fit.optimalityGap+'; delta='+robustDelta);fit={...fit,weights:sparseLengthFlowBasis(paths,fit.weights)};
  const denominator=Array(observed.length/4).fill(0),prediction=Array(observed.length).fill(0);
  paths.forEach((p,j)=>{const weight=fit.weights[j];if(!weight)return;for(let k=0;k<p.column.length;k++){prediction[k]+=weight*p.column[k];denominator[k>>2]+=weight*p.column[k];}});
  for(let k=0;k<prediction.length;k++)prediction[k]/=Math.max(1e-30,denominator[k>>2]);
  const costs=windows.map((w,t)=>w.classes.map(c=>{let gradient=0;for(let k=0;k<c.column.length;k+=4){const at=offsets[t]+k,total=c.column[k]+c.column[k+1]+c.column[k+2]+c.column[k+3],den=Math.max(1e-30,denominator[at>>2]);for(let b=0;b<4;b++)gradient+=Math.max(-robustDelta,Math.min(robustDelta,prediction[at+b]-observed[at+b]))*(c.column[k+b]-prediction[at+b]*total)/den;}return gradient;})),cost=(h:number,t:number,size:number)=>costs[t][windows[t].classes.findIndex(c=>c.reference===h&&c.size===size)],proposals:{reference:number;sizes:number[];value:number}[]=[],terminalProposals:{reference:number;sizes:number[];value:number}[]=[];
  for(let h=0;h<shares.length;h++){
   let states=new Map<number,{value:number;sizes:number[]}>([[0,{value:0,sizes:[]}]]);
   for(let t=0;t<windows.length;t++){const next=new Map<number,{value:number;sizes:number[]}>();for(const c of windows[t].classes)if(c.reference===h){let best:{value:number;sizes:number[]}|undefined;for(const[a,previous]of states)if(c.size-a>=-deletion&&c.size-a<=insertion){const candidate={value:previous.value+cost(h,t,c.size),sizes:[...previous.sizes,c.size]};if(!best||candidate.value<best.value)best=candidate;}if(best)next.set(c.size,best);}states=next;}
   terminalProposals.push(...[...states.values()].map(p=>({reference:h,...p})));const best=[...states.values()].reduce((a,b)=>a.value<b.value?a:b);proposals.push({reference:h,...best});
  }
  const current=paths.reduce((n,p,j)=>n+fit.weights[j]*p.sizes.reduce((n,size,t)=>n+cost(p.reference,t,size),0),0);optimalityGap=Math.max(0,current-Math.min(...proposals.map(p=>p.value)))/Math.max(1,observed.length);
  lengthFlowDiagnostics.observe?.({phase:'priced',paths:paths.length,delta:robustDelta,gap:optimalityGap,rmse:fit.rmse});
  if(optimalityGap<=1e-9)break;
  const previous=paths.length;let added=false;for(const p of terminalProposals)if(p.value<current)added=add(p.reference,p.sizes)||added;warm=[...fit.weights,...Array(paths.length-previous).fill(0)];if(!added)throw Error('Joint length path pricing did not converge: '+optimalityGap);
 }
 return{...fit,optimalityGap,referenceShares:shares.map((_,h)=>paths.reduce((n,p,j)=>n+(p.reference===h?fit.weights[j]:0),0)),paths:paths.map((p,j)=>({reference:p.reference,sizes:p.sizes,fraction:fit.weights[j]})),lengths:windows.map((w,t)=>shares.map((_,h)=>w.classes.flatMap(c=>c.reference===h?[{size:c.size,fraction:paths.reduce((n,p,j)=>n+(p.reference===h&&p.sizes[t]===c.size?fit.weights[j]:0),0)}]:[])))};
}
