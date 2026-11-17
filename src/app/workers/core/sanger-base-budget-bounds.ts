import type {LinearTerm} from './sanger-marginal-lp';
/** Intersect report envelopes with necessary consequences of the ORIGINAL
 * nonnegative departure budget. This does not fit, prune, reweight or phase a
 * sequence. Preserve the fitted percentages and every model-feasible profile.
 */
export function tightenBaseBudgetBounds(forms:LinearTerm[][],minimumMass:number,estimate:number[],lower:number[],upper:number[],budgetUpper:(terms:LinearTerm[])=>number|undefined){
 if(minimumMass>1e-8)for(let b=0;b<4;b++){
  const hi=budgetUpper(forms[b]),other=budgetUpper(forms.flatMap((f,k)=>k===b?[]:Array.from(f)));
  if(hi!==undefined)upper[b]=Math.min(upper[b],Math.max(estimate[b],Math.min(1,(hi+1e-7)/minimumMass)));
  if(other!==undefined)lower[b]=Math.max(lower[b],Math.min(estimate[b],Math.max(0,1-(other+1e-7)/minimumMass)));
 }
 for(let b=0;b<4;b++){
  lower[b]=Math.max(lower[b],Math.min(estimate[b],Math.max(0,1-upper.reduce((n,v,k)=>n+(k===b?0:v),0))));
  upper[b]=Math.min(upper[b],Math.max(estimate[b],Math.min(1,1-lower.reduce((n,v,k)=>n+(k===b?0:v),0))));
 }
}
