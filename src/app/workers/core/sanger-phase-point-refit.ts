import {tryMarginalWitness} from './sanger-witness-proposal';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import {evaluateQuadratic} from './sanger-quadratic-envelope';
import {MarginalLinearSolver,LP_INFINITY,validMarginalSolution,type MarginalProgram,type LinearTerm} from './sanger-marginal-lp';
/** Select a representative closest to the independently measured flank length
 * estimates, in the ORIGINAL joint feasible set. An exact feasible seed is a
 * certified zero-distance optimum. Otherwise solve the same L1 projection.
 * No estimate, origin or uncertainty constraint is fixed in the parent model. */
export function refineMeasuredPhases(program:MarginalProgram,point:number[],seed:{column:number;value:number}[],optimizer:SangerLinearOptimizer){
 const started=performance.now(),distance=(x:number[])=>seed.reduce((s,v)=>s+Math.abs(x[v.column]-v.value),0),beforeDistance=distance(point),clone=():MarginalProgram=>({...program,upper:program.upper.slice(),cost:program.cost.slice(),rows:program.rows.slice()});let reason:string|undefined;
 const finish=(x:number[],optimal:boolean,method:'zero-distance'|'joint-projection')=>({point:x,applied:true,valid:true,optimal,method,seconds:(performance.now()-started)/1000,beforeDistance,afterDistance:distance(x),reason});
 if(beforeDistance<=1e-12)return finish(point,true,'zero-distance');
 const trial=point.slice();for(const s of seed)trial[s.column]=s.value;
 if((program.quadratics||[]).every(q=>evaluateQuadratic(q,trial)<=q.upper+1e-7)){
  const p=clone();p.quadratics=undefined;for(const s of seed)p.rows.push({terms:[[s.column,1]],lower:s.value,upper:s.value});const answer=tryMarginalWitness(p,optimizer,10);if(answer&&validMarginalSolution(program,answer.x)&&distance(answer.x)<=1e-7)return finish(answer.x,true,'zero-distance');reason='Seed witness not completed in the original model; continue the unrestricted joint projection.';
 }
 const p=clone();p.cost=p.cost.map(()=>0);for(const s of seed){const d=p.upper.length;p.upper.push(program.upper[s.column]+Math.abs(s.value));p.cost.push(1);p.rows.push({terms:[[s.column,1],[d,-1]],lower:-LP_INFINITY,upper:s.value},{terms:[[s.column,-1],[d,-1]],lower:-LP_INFINITY,upper:-s.value});}
 const solver=new MarginalLinearSolver(p,optimizer);try{const answer=solver.solve(),x=answer.x.slice(0,program.upper.length);if(!validMarginalSolution(program,x))throw Error('Measured-length projection failed the original model check.');return finish(x,true,'joint-projection');}catch(error){return{point,applied:false,valid:validMarginalSolution(program,point),optimal:false,method:'joint-projection' as const,seconds:(performance.now()-started)/1000,beforeDistance,afterDistance:beforeDistance,reason:String(error)};}finally{solver.dispose();}
}
