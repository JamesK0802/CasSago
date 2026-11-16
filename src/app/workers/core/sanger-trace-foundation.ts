import type {SignalPathInput} from './sanger-path-common';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import {buildMarginalModel} from './sanger-marginal-model';
import {fitMarginalPoint} from './sanger-marginal-fit';
import {marginalViolation,type LinearTerm} from './sanger-marginal-lp';
/** Intermediate fit checkpoint, NOT a second inference algorithm or UI result.
 * The existing report and this checkpoint use the very same construction/fit.
 * No worst-case ranges, GT checks, sequence enumeration or display selection.
 * Narrow confidence is deliberately NOT inferred from a stable point estimate.
 */
export function fitTraceFoundation(input:SignalPathInput,optimizer:SangerLinearOptimizer,progress:(stage:string)=>void=()=>{}){
 const started=performance.now(),stages:{stage:string;seconds:number}[]=[];
 const tick=(stage:string)=>{stages.push({stage,seconds:(performance.now()-started)/1000});progress(stage);};
 const model=buildMarginalModel(input,tick),fit=fitMarginalPoint(model,optimizer,tick,progress);
 try{
  const evaluate=(terms:LinearTerm[])=>terms.reduce((n,[i,w])=>n+w*fit.point[i],0),violation=marginalViolation(model.program,fit.point);
  if(violation>1e-7)throw Error('Shared fit checkpoint failed the original physical model certificate.');
  const groups=model.drafts.flatMap(target=>target.groups.map(g=>{const ref=model.refs[g.reference],mass=Math.max(0,evaluate(g.mass));return{id:g.id,reference:ref.name,target:ref.reverse?ref.cuts.length-g.target:g.target+1,size:g.size,fraction:mass/model.shares[g.reference],bases:[...g.bases].map(([q,forms])=>({position:ref.reverse?ref.input.sequence.length-ref.offset-q-1:ref.offset+q,estimate:(ref.reverse?forms.slice().reverse():forms).map(f=>mass?Math.max(0,evaluate(f))/mass:0)}))};}));
  return{scope:'Intermediate shared physical-model fit, no range/identity/overall PASS claim.',model,point:fit.point,groups,representativeSelection:fit.representativeSelection,originalViolation:violation,seconds:(performance.now()-started)/1000,stages};
 }finally{fit.solver.dispose();}
}
