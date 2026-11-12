import {LP_INFINITY,type MarginalProgram} from './sanger-marginal-lp';
/** One pass of necessary interval implications for auxiliary master boxes.
 * Bounds come only from existing rows and nonnegative original boxes. A
 * roundoff margin enlarges every implied upper bound; no feasible point is
 * removed. Physical-model rows and original variable boxes are untouched. */
export function tightenAuxiliaryBoxes(p:MarginalProgram,firstAuxiliary:number){
 for(const row of p.rows)for(const sign of [1,-1]){const rhs=sign===1?row.upper:-row.lower;if(rhs>=LP_INFINITY)continue;let minimum=0,magnitude=Math.abs(rhs);for(const[i,w]of row.terms){const a=sign*w;if(a<0){const v=a*p.upper[i];minimum+=v;magnitude+=Math.abs(v);}}
  for(const[i,w]of row.terms){const a=sign*w;if(i<firstAuxiliary||a<=0)continue;const pad=8*(row.terms.length+4)*Number.EPSILON*(1+magnitude)/a,bound=(rhs-minimum)/a+pad;if(Number.isFinite(bound))p.upper[i]=Math.min(p.upper[i],Math.max(0,bound));}
 }
}
