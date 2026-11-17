import type {LinearTerm} from './sanger-marginal-lp';

/** Tighten a conditional base fraction without treating numerator and group
 * abundance as independent. Each callback returns an ORIGINAL-domain global
 * dual outer bound. For mass m in [a,b], a bound on n-r*m yields a safe ratio
 * bound after choosing a or b according to the bound's sign.
 *
 * Marker-defined groups can supply distinct lower/upper numerator forms:
 * F+M-G <= child-base mass <= F. This retains unobserved linkage.
 * These are model-feasible outer envelopes, not statistical intervals.
 */
export function conditionalRatioEnvelope(
 numeratorLower:LinearTerm[],numeratorUpper:LinearTerm[],mass:LinearTerm[],
 massRange:[number,number],estimate:number,range:[number,number],
 bound:(terms:LinearTerm[],maximize:boolean)=>number|undefined,
):[number,number]{
 if(massRange[0]<=1e-8)return range;
 const centered=(numerator:LinearTerm[])=>{
  const terms=new Map<number,number>();
  for(const[i,w]of numerator)terms.set(i,(terms.get(i)||0)+w);
  for(const[i,w]of mass)terms.set(i,(terms.get(i)||0)-estimate*w);
  return [...terms].filter(([,w])=>w!==0)as LinearTerm[];
 };
 let[lower,upper]=range;
 const lo=bound(centered(numeratorLower),false),hi=bound(centered(numeratorUpper),true);
 // Retain the original finite-precision objective margin and the fitted
 // point. No new physical row, noise limit or hypothesis is introduced.
 if(lo!==undefined&&Number.isFinite(lo)){
  const delta=lo-1e-7,safe=estimate+delta/massRange[delta<0?0:1];
  lower=Math.max(lower,Math.min(estimate,Math.max(0,safe)));
 }
 if(hi!==undefined&&Number.isFinite(hi)){
  const delta=hi+1e-7,safe=estimate+delta/massRange[delta>0?0:1];
  upper=Math.min(upper,Math.max(estimate,Math.min(1,safe)));
 }
 return[lower,upper];
}
