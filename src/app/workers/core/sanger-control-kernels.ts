import type {Prepared} from './sanger-path-common';
/** Position-specific stochastic dye responses calibrated ONLY on control.
 * Start with the global base/dye matrix. Remove surplus observed-colour mass
 * proportionally from its source rows, then send that mass to deficit colours.
 * Row sums stay one, columns stay nonnegative, and the known control base
 * mixture exactly reproduces the measured control. Unobserved bases retain
 * the global response. Equal contexts reuse the same matrix, for any input. */
export function measuredControlKernels(refs:Prepared[],shares:number[],means:number[][],control:number[][],nucleotideWeights?:number[][]):number[][][]{
 return control.map((measured,p)=>{const weights=nucleotideWeights?.[p]||Array.from({length:4},(_,b)=>refs.reduce((n,r,h)=>n+(r.sequence[p]==='ACGT'[b]?shares[h]:0),0)),sum=measured.reduce((n,v)=>n+v,0);if(sum<=0||Math.abs(weights.reduce((n,v)=>n+v,0)-1)>1e-9)return means;const observed=measured.map(v=>v/sum),predicted=Array.from({length:4},(_,j)=>weights.reduce((n,w,b)=>n+w*means[b][j],0));if(predicted.every((v,j)=>Math.abs(v-observed[j])<=1e-12))return means;
  const deficit=observed.map((v,j)=>Math.max(0,v-predicted[j])),total=deficit.reduce((n,v)=>n+v,0),shrink=predicted.map((v,j)=>v>0?Math.min(1,observed[j]/v):1);
  return means.map((row,b)=>{if(weights[b]<=1e-12)return row;const kept=row.map((v,j)=>v*shrink[j]),lost=Math.max(0,1-kept.reduce((n,v)=>n+v,0)),corrected=kept.map((v,j)=>v+(total?lost*deficit[j]/total:0)),normalizer=corrected.reduce((n,v)=>n+v,0);return corrected.map(v=>Math.max(0,v/normalizer));});
 });
}
