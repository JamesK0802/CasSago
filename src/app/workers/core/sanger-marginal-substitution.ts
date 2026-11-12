import {fixedGeometryWitness} from './sanger-fixed-geometry-witness';
import type {SangerSubstitutionSummary} from '../../models/sanger-analysis';
import type {MarginalTargetState} from './sanger-marginal-engine';
import type {Prepared} from './sanger-path-common';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import {MarginalLinearSolver,LP_INFINITY,validMarginalSolution,type LinearTerm,type MarginalProgram} from './sanger-marginal-lp';
/** Joint observable-edit bounds retain all original measured constraints and
 * the learned departure budget. The auxiliary union mass lies between the
 * largest position departure and their sum, capped by the family mass. */
export function boundMarginalSubstitutions(target:MarginalTargetState,refs:Prepared[],shares:number[],optimizer:SangerLinearOptimizer,progress:(stage:string)=>void=()=>{}):SangerSubstitutionSummary[]{
 const p:MarginalProgram={...target.program,upper:target.program.upper.slice(),cost:target.program.cost.slice(),rows:target.program.rows.slice()},termsByReference=refs.map(()=>[]as LinearTerm[]),lowerRows=new Set<number>(),upperRows=new Set<number>(),unions:{column:number;mass:LinearTerm[];departures:LinearTerm[][]}[]=[];
 // The original m <= family-mass row supplies the physical bound. A loose
 // finite reference-share box avoids numerically duplicating a derived
 // report interval at the optimum, without changing either endpoint.
 for(const g of target.groups){const m=p.upper.length;p.upper.push(shares[g.reference]);p.cost.push(0);termsByReference[g.reference].push([m,1]);upperRows.add(p.rows.length);p.rows.push({terms:[[m,1],...g.mass.map(([i,w])=>[i,-w]as LinearTerm)],lower:-LP_INFINITY,upper:0});const all:LinearTerm[]=[],departures=[...(g.departures?.values()||[])];
  for(const norm of departures){lowerRows.add(p.rows.length);p.rows.push({terms:[...norm,[m,-1]],lower:-LP_INFINITY,upper:0});all.push(...norm);}
  upperRows.add(p.rows.length);p.rows.push({terms:[[m,1],...all.map(([i,w])=>[i,-w]as LinearTerm)],lower:-LP_INFINITY,upper:0});unions.push({column:m,mass:g.mass,departures});
 }
 const evaluate=(terms:LinearTerm[],x:number[])=>terms.reduce((n,[i,w])=>n+w*x[i],0),lo:number[]=[],hi:number[]=[],pending:{reference:number;maximum:boolean}[]=[];
 const budgetRow=p.dualHint?p.rows[p.dualHint.budgetRow]:undefined,budgetTerms=new Map(budgetRow?.terms||[]),budget=budgetRow&&p.cost.every((c,i)=>(c===0||c===1)&&c===(budgetTerms.get(i)||0))?budgetRow.upper:Infinity;
 const seeds=[target.primaryPoint,target.point].filter((x):x is number[]=>!!x&&x.length===target.program.upper.length),zeroSeeds=seeds.map(seed=>{const x=seed.slice();for(const g of target.groups)for(const flow of g.flow?.values()||[]){flow.changed.forEach(i=>x[i]=0);for(const c of flow.unchanged)x[c.column]=evaluate(c.cap,seed);}return x;});
 const lift=(seed:number[],maximum:boolean)=>{const x=seed.slice();while(x.length<p.upper.length)x.push(0);for(const u of unions){const amounts=u.departures.map(terms=>evaluate(terms,x));x[u.column]=maximum?Math.max(0,Math.min(p.upper[u.column],evaluate(u.mass,x),amounts.reduce((n,v)=>n+v,0))):Math.max(0,...amounts);}return x;};
 const witnessedEndpoint=(h:number,maximum:boolean)=>{const endpoint=maximum?Math.min(shares[h],budget):0;for(const seed of maximum?seeds:[...seeds,...zeroSeeds]){const x=lift(seed,maximum),value=evaluate(termsByReference[h],x);if(Math.abs(value-endpoint)<=1e-7&&validMarginalSolution(p,x))return endpoint;}
  // A feasible point attaining a proved global upper bound certifies the
  // endpoint (or the nonnegative lower bound zero). Search the existing geometry first; if unattainable,
  // the unrestricted original conic solve still runs unchanged.
  if(target.point&&target.program.localBlocks?.length){const free=[...target.program.localBlocks.flat(),...unions.map(u=>u.column)];if(fixedGeometryWitness(p,target.point,free,termsByReference[h],endpoint,optimizer))return endpoint;}
  return undefined;};
 // For a fixed physical base flow, the minimizing union is max(position
 // departures), and the maximizing union is min(family mass, sum departures).
 // Thus upper union rows are redundant for MIN, lower rows for MAX. Restore
 // every auxiliary mass and verify the COMPLETE original program afterwards.
 // This objective-aware elimination changes neither endpoint nor hypotheses.
 for(const maximum of [false,true]){const redundant=maximum?lowerRows:upperRows;let solver:MarginalLinearSolver|undefined;
  try{for(let h=0;h<refs.length;h++){progress('target '+((target.groups.find(g=>g.reference===h)?.target??0)+1)+' · '+refs[h].name+' · '+(maximum?'maximum':'minimum'));const witnessed=witnessedEndpoint(h,maximum);if(witnessed!==undefined){(maximum?hi:lo)[h]=witnessed;continue;}try{if(!solver){solver=new MarginalLinearSolver({...p,localBlocks:undefined,rows:p.rows.filter((_,k)=>!redundant.has(k))},optimizer);solver.setTimeSlice(600);}const answer=solver.solveToCompletion(termsByReference[h],maximum),x=lift(answer.x,maximum),value=evaluate(termsByReference[h],x);if(!validMarginalSolution(p,x)||Math.abs(value-answer.value)>1e-7)throw Error('Substitution endpoint did not lift into the original union model.');(maximum?hi:lo)[h]=value;}catch{pending.push({reference:h,maximum});}}}finally{solver?.dispose();}
 }
 if(pending.length){const original=new MarginalLinearSolver({...p,localBlocks:undefined},optimizer);original.setTimeSlice(600);try{for(const q of pending){progress('target '+((target.groups.find(g=>g.reference===q.reference)?.target??0)+1)+' · '+refs[q.reference].name+' · retry original union '+(q.maximum?'maximum':'minimum'));(q.maximum?hi:lo)[q.reference]=original.solveToCompletion(termsByReference[q.reference],q.maximum).value;}}finally{original.dispose();}}
 return refs.map((ref,h)=>{const share=shares[h],t=target.groups.find(g=>g.reference===h)?.target??0;return{reference:ref.name,target:ref.reverse?ref.cuts.length-t:t+1,scope:'all-length-groups',method:'joint-observable-departure-envelope',fractionRange:[Math.max(0,(lo[h]-1e-7)/share),Math.min(1,(hi[h]+1e-7)/share)],zeroFeasible:lo[h]<=1e-7};});
}
