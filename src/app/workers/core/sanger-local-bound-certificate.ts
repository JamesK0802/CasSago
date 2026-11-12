import {quadraticCut} from './sanger-quadratic-envelope';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import type {MarginalInternalGroup} from './sanger-marginal-engine';
import {MarginalLinearSolver,LP_INFINITY,type MarginalProgram,type LinearTerm} from './sanger-marginal-lp';
/** Restricted coordinate solves are only proposals. Lift their duals into the
 * ORIGINAL program and compute a global Lagrangian bound (including every
 * omitted column's box contribution). Accept only a matching feasible primal
 * witness for exact bounds. outerBound exposes the independently valid
 * global dual envelope even when the local primal does not attain it. No row, geometry,
 * allele, noise allowance or measurement is removed from a reported bound. */
export function localBoundCertificates(p:MarginalProgram,anchor:number[],groups:MarginalInternalGroup[],optimizer:SangerLinearOptimizer){
 const allLocal=new Set(p.localBlocks?.flat()||[]);
 const adjacent=p.upper.map(()=>[]as number[]),values=p.rows.map((r,k)=>{for(const[i]of r.terms)adjacent[i].push(k);return r.terms.reduce((n,[i,w])=>n+w*anchor[i],0);}),equalities=p.rows.flatMap((r,k)=>r.lower===r.upper&&Math.abs(r.lower)<LP_INFINITY&&r.terms.length&&r.terms.every(([i])=>!allLocal.has(i))?[k]:[]).sort((a,b)=>p.rows[a].terms.length-p.rows[b].terms.length),cache=new Map<number,{solver:MarginalLinearSolver;columns:number[];indices:Map<number,number>;rows:number[];fixed:number[];globalRows:MarginalProgram['rows']}|null>(),stats={queries:0,certified:0,outer:0,maximumSlack:0,fallback:0,seconds:0};
 const geometryColumns=p.localBlocks?p.upper.flatMap((_,i)=>!allLocal.has(i)?[i]:[]):[];
 const geometryMap=new Map(geometryColumns.map((i,j)=>[i,j]));
 const geometryProgram:MarginalProgram|undefined=p.localBlocks?{
  upper:geometryColumns.map(i=>p.upper[i]),cost:geometryColumns.map(()=>0),
  rows:p.rows.filter(r=>r.terms.every(([i])=>geometryMap.has(i))).map(r=>({...r,terms:r.terms.map(([i,w])=>[geometryMap.get(i)!,w]as LinearTerm)})),
  quadratics:p.quadratics?.map(q=>({...q,terms:q.terms.map(r=>r.map(([i,w])=>[geometryMap.get(i)!,w]as LinearTerm))}))
 }:undefined;
 // Necessary cone support planes give a safe LP outer relaxation. The
 // original nonlinear model remains intact for inference and joint retrieval.
 if(geometryProgram){const x=geometryProgram.upper.map(()=>0);for(const q of geometryProgram.quadratics||[]){const columns=new Set(q.terms.flatMap(row=>row.map(([i])=>i)));for(const i of columns){x[i]=1;const cut=quadraticCut(q,x);x[i]=0;if(cut.value>0)geometryProgram.rows.push({terms:cut.terms,lower:-LP_INFINITY,upper:cut.upper+1e-10*(1+Math.abs(cut.upper))});}}geometryProgram.quadratics=undefined;}
 const geometry=geometryProgram?new MarginalLinearSolver(geometryProgram,optimizer):undefined;
 const hint={image:new Map<number,number>(),constant:0,magnitude:0};for(let k=0;k<(p.dualHint?.shadow.length||0);k++){let y=p.dualHint!.shadow[k];const r=p.rows[k];if(!Number.isFinite(y)||!r)continue;if(y>0&&r.lower<=-LP_INFINITY||y<0&&r.upper>=LP_INFINITY)y=0;if(!y)continue;const v=y*(y>0?r.lower:r.upper);hint.constant+=v;hint.magnitude+=Math.abs(v);for(const[i,w]of r.terms)hint.image.set(i,(hint.image.get(i)||0)+y*w);}
 for(const plane of p.dualHint?.planes||[]){hint.constant-=plane.upper;hint.magnitude+=Math.abs(plane.upper);for(const[i,w]of plane.terms)hint.image.set(i,(hint.image.get(i)||0)-w);}
 const budgetTerms=new Map(p.dualHint?p.rows[p.dualHint.budgetRow].terms:[]),sameDepartureCost=p.cost.every((c,i)=>c===(budgetTerms.get(i)||0));let required=hint.constant,requiredMagnitude=hint.magnitude;for(const[i,w]of hint.image){const v=Math.min(0,p.cost[i]-w)*p.upper[i];required+=v;requiredMagnitude+=Math.abs(v);}
 function create(q:number){
  const columns=[...new Set(groups.flatMap(g=>{const f=g.flow?.get(q);return f?[...f.unchanged.map(c=>c.column),...f.changed]:[];}))],indices=new Map(columns.map((v,i)=>[v,i]));if(!columns.length)return null;
  // The fitted representative may spend its entire optional departure budget
  // elsewhere. That choice is NOT a constraint on this coordinate's range.
  // Reserve only an independently certified mandatory departure elsewhere.
  let outside=required;for(const i of columns){const c=p.cost[i],image=hint.image.get(i)||0;outside+=(Math.min(0,-image)-Math.min(0,c-image))*p.upper[i];}outside=sameDepartureCost?Math.max(0,outside-128*Number.EPSILON*(1+requiredMagnitude)):0;
  const rows=[...new Set(columns.flatMap(i=>adjacent[i]))].sort((a,b)=>a-b),fixed:number[]=[],globalRows:MarginalProgram['rows']=[],program:MarginalProgram={upper:columns.map(i=>p.upper[i]),cost:columns.map(()=>0),rows:rows.map(k=>{const original=p.rows[k],budget=k===p.dualHint?.budgetRow,r=budget?{terms:original.terms.filter(([i])=>indices.has(i)),lower:0,upper:original.upper-outside}:original,terms=r.terms.flatMap(([i,w])=>indices.has(i)?[[indices.get(i)!,w]as LinearTerm]:[]),constant=budget?0:r.terms.reduce((s,[i,w])=>s+(indices.has(i)?0:w*anchor[i]),0);globalRows.push(r);fixed.push(original.terms.reduce((s,[i,w])=>s+(indices.has(i)?0:w*anchor[i]),0));const lower=r.lower<=-LP_INFINITY?-LP_INFINITY:r.lower-constant,upper=r.upper>=LP_INFINITY?LP_INFINITY:r.upper-constant,proposal=terms.reduce((n,[i,w])=>n+w*Math.max(0,Math.min(p.upper[columns[i]],anchor[columns[i]])),0);
   // This LP proposes dual multipliers, never a replacement physical model.
   // Fixing a verified finite-precision geometry can otherwise make a tiny
   // conditional LP infeasible. Include its boxed anchor in the proposal;
   // lift the dual using UNCHANGED original global rows below.
   return{terms,lower:lower<=-LP_INFINITY?-LP_INFINITY:Math.min(lower,proposal-1e-10),upper:upper>=LP_INFINITY?LP_INFINITY:Math.max(upper,proposal+1e-10)};})};
  return{solver:new MarginalLinearSolver(program,optimizer),columns,indices,rows,fixed,globalRows};
 }
 function calculate(q:number,terms:LinearTerm[],maximize=false,exact=true){const start=performance.now();stats.queries++;try{let s=cache.get(q);if(s===undefined){
  // Coordinate queries are consecutive. Retaining every native LP (and its
  // factorization) across all coordinates can exhaust the WASM heap. Release
  // the previous workspace; every original column and constraint is rebuilt
  // unchanged when that coordinate is requested again.
  for(const previous of cache.values())previous?.solver.dispose();cache.clear();s=create(q);cache.set(q,s);
 }if(!s){
  // Guard cycles can contain only phase/origin expressions and no local
  // substitution variables. Their existing geometry relaxation already
  // supplies a valid global outer bound; a full joint solve is unnecessary.
  if(!exact&&geometry&&terms.every(([i])=>geometryMap.has(i))){try{const answer=geometry.solve(terms.map(([i,w])=>[geometryMap.get(i)!,w]as LinearTerm),maximize),lower=geometry.objectiveLowerBound()??((maximize?-1:1)*answer.value-1e-7),bound=maximize?-lower:lower,at=terms.reduce((n,[i,w])=>n+w*anchor[i],0);if(Number.isFinite(bound)){stats.outer++;return maximize?Math.max(bound,at):Math.min(bound,at);}}catch{/* Preserve the full-model fallback if this certificate is incomplete. */}}
  stats.fallback++;return;
 }
  const local=terms.flatMap(([i,w])=>s!.indices.has(i)?[[s!.indices.get(i)!,w]as LinearTerm]:[]),fixed=terms.reduce((n,[i,w])=>n+(s!.indices.has(i)?0:w*anchor[i]),0),answer=s.solver.solve(local,maximize),primal=fixed+answer.value,dual=s.solver.rowDuals(),sign=maximize?-1:1,reduced=new Map<number,number>();for(const[i,w]of terms)reduced.set(i,(reduced.get(i)||0)+sign*w);let lower=0,magnitude=0,budgetMultiplier=0;
  for(let k=0;k<dual.length;k++){let y=dual[k];if(!Number.isFinite(y))return;const r=s.globalRows[k];if(y>0&&r.lower<=-LP_INFINITY||y<0&&r.upper>=LP_INFINITY)y=0;if(!y)continue;if(s.rows[k]===p.dualHint?.budgetRow)budgetMultiplier=Math.max(0,-y);const v=y*(y>0?r.lower:r.upper);lower+=v;magnitude+=Math.abs(v);for(const[i,w]of r.terms)reduced.set(i,(reduced.get(i)||0)-y*w);}
  const certificate=(alpha:number)=>{const coefficients=new Map(reduced);let result=lower+alpha*hint.constant,scale=magnitude+alpha*hint.magnitude;if(alpha)for(const[i,w]of hint.image)coefficients.set(i,(coefficients.get(i)||0)-alpha*w);
   if(geometry){try{const terms=[...coefficients].flatMap(([i,w])=>geometryMap.has(i)&&w?[[geometryMap.get(i)!,w]as LinearTerm]:[]),answer=geometry.solve(terms),bound=geometry.objectiveLowerBound()??answer.value-1e-7;result+=bound;scale+=Math.abs(bound);for(const[i,r]of coefficients)if(!geometryMap.has(i)&&r<0){const v=r*p.upper[i];result+=v;scale+=Math.abs(v);}return result-64*Number.EPSILON*(1+scale);}catch{/* The original box certificate remains valid if a geometry query is incomplete. */}}
   // Coordinate ascent in ORIGINAL equality multipliers. Signed origin
   // recurrences and phase-flow conservation connect an omitted context to
   // its measured state bounds, without fixing its origin allocation.
   for(let sweep=0;sweep<3;sweep++)for(const k of (sweep%2?equalities:equalities.slice().reverse())){const r=p.rows[k];if(r.terms.some(([i])=>s!.indices.has(i)))continue;const breaks=r.terms.map(([i,w])=>({at:(coefficients.get(i)||0)/w,weight:Math.abs(w)*p.upper[i]})).sort((a,b)=>a.at-b.at);let derivative=r.lower+r.terms.reduce((n,[i,w])=>n+(w<0?-w*p.upper[i]:0),0),y=0;for(const point of breaks){derivative-=point.weight;if(derivative<=1e-15){y=point.at;break;}}if(!y||!Number.isFinite(y))continue;const gain=y*r.lower+r.terms.reduce((n,[i,w])=>{const c=coefficients.get(i)||0;return n+(Math.min(0,c-y*w)-Math.min(0,c))*p.upper[i];},0);if(gain< -1e-12)continue;result+=y*r.lower;scale+=Math.abs(y*r.lower);for(const[i,w]of r.terms)coefficients.set(i,(coefficients.get(i)||0)-y*w);}
   for(const[i,r]of coefficients)if(r<0){const v=r*p.upper[i];result+=v;scale+=Math.abs(v);}return result-64*Number.EPSILON*(1+scale);
  };
  // Reuse the original whole-signal departure dual. This retains mandatory
  // substitutions at other cycles when bounding an individual coordinate.
  // Both choices are independently valid global bounds; take the tighter one.
  const lowerBound=Math.max(certificate(0),budgetMultiplier?certificate(budgetMultiplier):-Infinity);let global=maximize?-lowerBound:lowerBound;
  const anchorValue=terms.reduce((n,[i,w])=>n+w*anchor[i],0);if(exact){if(maximize?global<anchorValue-1e-7:global>anchorValue+1e-7){stats.fallback++;return;}}else{
   // The verified fitted allocation has finite-precision row residuals. An
   // outer envelope can safely include it as well as the independently valid
   // global bound, rather than launching a full-model solve to reconcile a
   // small numerical discrepancy. Exact endpoint claims remain unchanged.
   global=maximize?Math.max(global,anchorValue):Math.min(global,anchorValue);
  }
  const gap=maximize?global-primal:primal-global;
  if(!Number.isFinite(global)||exact&&(gap< -1e-8||gap>2e-8)){stats.fallback++;return;}
  // The local rows are the original rows with the fixed columns substituted.
  // All untouched rows still have the verified anchor value.
  if(exact)for(let k=0;k<s.rows.length;k++){const r=p.rows[s.rows[k]],v=s.fixed[k]+s.solver.program.rows[k].terms.reduce((n,[i,w])=>n+w*answer.x[i],0);if(v<r.lower-1e-7||v>r.upper+1e-7){stats.fallback++;return;}}
  if(exact)stats.certified++;else{stats.outer++;stats.maximumSlack=Math.max(stats.maximumSlack,gap);}return global;
 }catch{stats.fallback++;return;}finally{stats.seconds+=(performance.now()-start)/1000;}}
 return{bound:(q:number,t:LinearTerm[],maximum=false)=>calculate(q,t,maximum,true),outerBound:(q:number,t:LinearTerm[],maximum=false)=>calculate(q,t,maximum,false),stats,dispose(){geometry?.dispose();for(const s of cache.values())s?.solver.dispose();}};
}
