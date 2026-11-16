import type {LinearTerm,MarginalProgram} from './sanger-marginal-lp';
import {LP_INFINITY} from './sanger-marginal-lp';
import {linearDualOuterLower} from './sanger-linear-dual-outer';
/** Optional read-only development audit; never supplies bounds or fit inputs. */
export const budgetEnvelopeDiagnostics:{observe?:(value:{group:string;position:number;beforeLower:number[];beforeUpper:number[];afterLower:number[];afterUpper:number[]})=>void}={};

/** Necessary upper envelope from ONE EXISTING nonnegative budget row.
 * If x>=0, a>=0 and a.x<=B, then c.x<=B max_i(c_i/a_i)
 * for nonnegative c supported by positive a. No rows/boxes/noise are changed;
 * this is not an attained optimum or a statistical confidence interval.
 * The original inference rows are append-only; invalidate on row replacement.
 */
export function nonnegativeBudgetEnvelope(p:MarginalProgram){
 const id=p.dualHint?.budgetRow,row=id===undefined?undefined:p.rows[id],terms=row?.terms,B=row?.upper;
 const cost=new Map<number,number>();let valid=!!row&&B!==undefined&&Number.isFinite(B)&&B>=0&&B<LP_INFINITY;
 if(valid)for(const[i,w]of terms!){if(!Number.isSafeInteger(i)||i<0||i>=p.upper.length||!Number.isFinite(w)||w<0||cost.has(i)){valid=false;break;}cost.set(i,w);}
 // The existing BLIND primary dual can prove that some departure cost is
 // mandatory outside the queried columns. Reuse that certificate; never
 // reserve the representative's arbitrary allocation or shrink the budget.
 let required=0,image:Float64Array|undefined,magnitude:Float64Array|undefined,operations:Uint32Array|undefined;
 const hint=p.dualHint;
 if(valid&&hint&&p.cost.every((c,i)=>Number.isFinite(c)&&c>=0&&c===(cost.get(i)||0))){
  const planes=hint.planes||[],linear:MarginalProgram={upper:p.upper,cost:p.cost,rows:p.rows.concat(planes.map(plane=>({terms:plane.terms,lower:-LP_INFINITY,upper:plane.upper})))},dual=new Float64Array(linear.rows.length);
  for(let k=0;k<Math.min(hint.shadow.length,p.rows.length);k++){const r=p.rows[k];let y=hint.shadow[k];if(!Number.isFinite(y)||y>0&&r.lower<=-LP_INFINITY||y<0&&r.upper>=LP_INFINITY)y=0;dual[k]=y;}
  planes.forEach((_,k)=>dual[p.rows.length+k]=-1);
  const lower=linearDualOuterLower(linear,p.cost,dual);
  if(lower!==undefined&&lower>0){
   required=lower;image=new Float64Array(p.upper.length);magnitude=new Float64Array(p.upper.length);operations=new Uint32Array(p.upper.length);
   for(let k=0;k<dual.length;k++)if(dual[k])for(const[i,w]of linear.rows[k].terms){const v=dual[k]*w;image[i]+=v;magnitude[i]+=Math.abs(v);operations[i]++;}
  }
 }
 return( objective:LinearTerm[]):number|undefined=>{
  if(!valid||p.rows[id!]!==row||row!.terms!==terms||row!.upper!==B)return;
  const coefficients=new Map<number,number>(),errors=new Map<number,number>();let count=0;
  for(const[i,w]of objective){
   if(!Number.isSafeInteger(i)||i<0||i>=p.upper.length||!Number.isFinite(w)||w<0)return;
   if(!w)continue;const a=cost.get(i);if(!(a!>0))return;
   const old=coefficients.get(i)||0,sum=old+w,part=sum-old,error=(old-(sum-part))+(w-part);if(!Number.isFinite(sum))return;
   coefficients.set(i,sum);errors.set(i,(errors.get(i)||0)+Math.abs(error));count++;
  }
  if(!count)return 0;
  const guard=(4*count+32)*Number.EPSILON;if(guard>=.01)return;
  let ratio=0;
  for(const[i,w]of coefficients){const upper=w+(errors.get(i)||0)+guard*(1+w),r=upper/cost.get(i)!;if(!Number.isFinite(r))return;ratio=Math.max(ratio,r+guard*(1+r));}
  let outside=0;
  if(required>0&&image&&magnitude&&operations){
   let loss=0;
   for(const i of coefficients.keys()){
    // Upper-bound the removed dual contribution despite cancellation. A
    // larger subtraction only weakens the necessary outside-cost bound.
    const n=2*operations[i]+8,gamma=n*Number.EPSILON,roundoff=gamma<.01?gamma/(1-gamma)*(1+magnitude[i]):Infinity;
    loss+=Math.min(cost.get(i)!,Math.max(0,image[i]+roundoff))*p.upper[i];
   }
   outside=Math.max(0,required-loss-guard*(1+Math.abs(required)+Math.abs(loss)));
  }
  const remaining=Math.max(0,B!-outside),value=remaining*ratio,bound=value+guard*(1+Math.abs(value));return Number.isFinite(bound)?bound:undefined;
 };
}
