import {cachedLinearDualOuter} from './sanger-linear-dual-outer';
import {quadraticCut} from './sanger-quadratic-envelope';
import {LP_INFINITY,type MarginalProgram,type LinearTerm} from './sanger-marginal-lp';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';

/** Report-only joint LP proposals retain EVERY original column/linear row.
 * Implied SOC support planes strengthen the relaxation; original cones remain
 * unchanged in the scientific model. Independently checked full-row/box duals
 * are outer envelopes, not attained molecular solutions or confidence limits.
 */
export function jointBoundCertificates(original:MarginalProgram,optimizer:SangerLinearOptimizer){
 const rows=original.rows.slice(),x=original.upper.map(()=>0);
 for(const plane of original.dualHint?.planes||[])rows.push({terms:plane.terms,lower:-LP_INFINITY,upper:plane.upper});
 for(const q of original.quadratics||[]){
  const columns=new Set(q.terms.flatMap(form=>form.map(([i])=>i)));
  for(const i of columns){x[i]=1;const cut=quadraticCut(q,x);x[i]=0;
   if(cut.value>0)rows.push({terms:cut.terms,lower:-LP_INFINITY,upper:cut.upper+1e-10*(1+Math.abs(cut.upper))});
  }
 }
 const relaxation:MarginalProgram={...original,immutableRows:true,rows,quadratics:undefined,quadraticCost:undefined,localBlocks:undefined,originPricing:undefined},linear=cachedLinearDualOuter(relaxation,optimizer,30,true),incomplete=new Set<string>();
 let unavailable=false;
 return{stats:linear.stats,bound(terms:LinearTerm[],maximum=false):number|undefined{
  const key=JSON.stringify([maximum,terms]);if(unavailable||incomplete.has(key))return;
  const partial=linear.stats.partial,failed=linear.stats.failures,invalid=linear.stats.invalidDuals;
  try{const lo=linear.lower(terms,maximum);
   // A stopped proposal must not silently count as completed detail. The
   // caller falls back to a certified solve of the full original model.
   if(linear.stats.partial!==partial||linear.stats.failures!==failed||linear.stats.invalidDuals!==invalid){incomplete.add(key);unavailable=true;return;}
   return maximum?-lo:lo;
  }catch{incomplete.add(key);unavailable=true;return;}
 },dispose(){linear.dispose();incomplete.clear();}};
}
