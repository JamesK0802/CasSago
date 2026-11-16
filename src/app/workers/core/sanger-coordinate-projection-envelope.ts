import {VerifiedConicSolver} from './sanger-conic-solver';
import {linearDualOuterLower} from './sanger-linear-dual-outer';
import {LP_INFINITY,type MarginalProgram,type LinearTerm} from './sanger-marginal-lp';
import {packLinearTerms} from './sanger-packed-terms';
import type {MarginalInternalGroup} from './sanger-marginal-engine';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';

/** Optional numerical OUTER envelope. Keep shared variables and all coordinates
 * whose unchanged bases/free insertions depend on edit geometry. Omitted local
 * variables retain their full original boxes, never their fitted values.
 * Projected solves supply multipliers ONLY: lift them to every original row,
 * original SOC support plane and original box for an independent global bound.
 * No projected primal is accepted as an original fit or genotype.
 */
export function coordinateProjectionEnvelope(p:MarginalProgram,anchor:number[],groups:MarginalInternalGroup[],optimizer:SangerLinearOptimizer){
 const sensitive=new Set<number>(),local=new Set(p.localBlocks?.flat()||[]),keep=new Set<number>();
 const stats={attempts:0,certifiedOuter:0,failures:0,seconds:0,preparationSeconds:0,originalColumns:p.upper.filter(v=>v>0).length,projectedColumns:0,skippedWithoutReduction:false,lastFailure:undefined as string|undefined};
 let solver:VerifiedConicSolver|undefined,columns:number[]=[],index=new Map<number,number>(),originalRows:number[]=[],ready=false,unavailable=false;
 const cache=new Map<string,number>();
 for(const g of groups)for(const[q,f]of g.flow||[])if(new Set(f.unchanged.map(c=>c.base)).size>1||g.edges.some(e=>e.free.has(q)))sensitive.add(q);
 for(const g of groups)for(const e of g.edges)for(const q of e.free.keys())sensitive.add(q);
 // Retain coordinates with resolved departure in the existing BLIND anchor.
 // This only chooses a numerical outer proposal workspace: omitted variables
 // remain free, and every bound is still lifted through all original rows.
 // Use the original primal certificate tolerance to ignore numerical dust.
 const departure=new Map<number,number>();
 for(const g of groups)for(const[q,f]of g.flow||[])departure.set(q,(departure.get(q)||0)+f.changed.reduce((sum,i)=>sum+Math.max(0,anchor[i]),0));
 for(const[q,mass]of departure)if(mass>1e-7)sensitive.add(q);
 const add=(i:number)=>{if(p.upper[i]>0)keep.add(i);};
 p.upper.forEach((_,i)=>{if(!local.has(i))add(i);});
 for(const q of p.quadratics||[])for(const form of q.terms)for(const[i]of form)add(i);
 for(const q of sensitive)for(const g of groups){const f=g.flow?.get(q);if(f){f.unchanged.forEach(c=>add(c.column));f.changed.forEach(add);}for(const e of g.edges)for(const i of e.free.get(q)||[])add(i);}
 // This chooses whether a numerical proposal is worth allocating; it does
 // not restrict any original hypothesis, constraint or inference pipeline.
 if(!optimizer.conic||p.quadraticCost?.length||!local.size||keep.size*2>stats.originalColumns){unavailable=true;stats.skippedWithoutReduction=true;}
 const prepare=()=>{
  if(ready||unavailable)return;const start=performance.now();originalRows=[];columns=[...keep].sort((a,b)=>a-b);index=new Map(columns.map((i,j)=>[i,j]));stats.projectedColumns=columns.length;
  const projection:MarginalProgram={upper:columns.map(i=>p.upper[i]),cost:columns.map(()=>0),rows:[],quadratics:(p.quadratics||[]).map(q=>({...q,terms:q.terms.map(f=>f.filter(([i])=>p.upper[i]>0).map(([i,w])=>[index.get(i)!,w]as LinearTerm))}))};
  for(let k=0;k<p.rows.length;k++){
   const row=p.rows[k],terms:LinearTerm[]=[];let lo=0,hi=0,magnitude=0,at=0;
   for(const[i,w]of row.terms){const j=index.get(i);if(j!==undefined){terms.push([j,w]);at+=w*anchor[i];}else{const value=w*p.upper[i];lo+=Math.min(0,value);hi+=Math.max(0,value);magnitude+=Math.abs(value);}}
   // Outward roundoff/anchor guards affect this PROPOSAL only. The bound is
   // subsequently recomputed from unchanged original rows and all boxes.
   const guard=128*Number.EPSILON*(1+magnitude+Math.abs(at)),lower=row.lower<=-LP_INFINITY?-LP_INFINITY:Math.min(row.lower-hi-guard,at-1e-10),upper=row.upper>=LP_INFINITY?LP_INFINITY:Math.max(row.upper-lo+guard,at+1e-10);
   if(!terms.length&&lower<=0&&upper>=0)continue;
   projection.rows.push({terms:packLinearTerms(terms),lower,upper});originalRows.push(k);
  }
  solver=new VerifiedConicSolver(projection,optimizer.conic!);// Completion budget only: a larger shared model may need more than2s
  // to certify the unchanged projected problem. Never relax an original bound.
  solver.setTimeLimit(30);ready=true;stats.preparationSeconds+=(performance.now()-start)/1000;
 };
 return{stats,upper(q:number,terms:LinearTerm[]):number|undefined{
  if(unavailable||!sensitive.has(q))return;
  // Bounds here are only requested for nonnegative base/complement forms.
  if(terms.some(([i,w])=>!Number.isSafeInteger(i)||i<0||i>=p.upper.length||!Number.isFinite(w)||w<0||p.upper[i]>0&&!keep.has(i)))return;
  const key=JSON.stringify(terms),cached=cache.get(key);if(cached!==undefined)return cached;
  if(!terms.some(([i,w])=>w&&p.upper[i]))return 0;
  const start=performance.now();stats.attempts++;
  try{
   prepare();if(!solver)return;
   const mapped=terms.filter(([i])=>p.upper[i]>0).map(([i,w])=>[index.get(i)!,w]as LinearTerm);
   solver.solve(mapped,true);const dual=solver.rowDuals(),planes=solver.supportPlanes();
   const linear:MarginalProgram={upper:p.upper,cost:p.cost,rows:p.rows.concat(planes.map(plane=>({terms:plane.terms.map(([i,w])=>[columns[i],w]as LinearTerm),lower:-LP_INFINITY,upper:plane.upper})))},multipliers=new Float64Array(linear.rows.length),cost=new Float64Array(p.upper.length);
   for(let k=0;k<dual.length;k++)multipliers[originalRows[k]]+=dual[k];
   planes.forEach((_,k)=>multipliers[p.rows.length+k]=-1);
   let error=0;
   for(const[i,w]of terms){const v=-w,old=cost[i],sum=old+v,part=sum-old,residual=(old-(sum-part))+(v-part);cost[i]=sum;error+=Math.abs(residual)*p.upper[i];}
   if(!Number.isFinite(error))return;
   if(error)error+=64*Number.EPSILON*(1+error);
   const lower=linearDualOuterLower(linear,cost,multipliers);if(lower===undefined)return;
   const at=terms.reduce((s,[i,w])=>s+w*anchor[i],0),bound=Math.max(at,-lower+error+64*Number.EPSILON*(1+Math.abs(lower)+error));
   if(!Number.isFinite(bound))return;
   stats.certifiedOuter++;cache.set(key,bound);if(cache.size>16)cache.delete(cache.keys().next().value!);return bound;
  }catch(error){stats.failures++;stats.lastFailure=String(error);solver?.dispose();solver=undefined;ready=false;return;}finally{stats.seconds+=(performance.now()-start)/1000;}
 },dispose(){solver?.dispose();solver=undefined;cache.clear();originalRows=[];columns=[];index.clear();keep.clear();sensitive.clear();}};
}
