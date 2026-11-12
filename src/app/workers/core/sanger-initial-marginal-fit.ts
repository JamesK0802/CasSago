import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import {zeroCostProjection} from './sanger-zero-cost-projection';
import {MarginalLinearSolver,LP_INFINITY,validMarginalSolution,type MarginalProgram} from './sanger-marginal-lp';
/** A verified zero-cost face is a GLOBAL optimum of the same nonnegative
 * departure objective. This presolver never uses sample type or allele count.
 * Failure merely resumes the full original problem with every variable. */
export function verifiedZeroDeparturePoint(program:MarginalProgram,seed:{column:number;value:number}[],optimizer:SangerLinearOptimizer){
 const started=performance.now(),projection=zeroCostProjection(program);if(!projection)return;
 const p=projection.program;for(const s of seed){const column=projection.indices.get(s.column);if(column===undefined)continue;const d=p.upper.length;p.upper.push(p.upper[column]+Math.abs(s.value));p.cost.push(1);p.rows.push({terms:[[column,1],[d,-1]],lower:-LP_INFINITY,upper:s.value},{terms:[[column,-1],[d,-1]],lower:-LP_INFINITY,upper:-s.value});}
 const solver=new MarginalLinearSolver(p,optimizer);try{solver.setTimeSlice(15);const answer=solver.solve(),x=projection.lift(answer.x),value=program.cost.reduce((s,c,i)=>s+c*x[i],0);if(value!==0||!validMarginalSolution(program,x))return;return{x,value:0,seconds:(performance.now()-started)/1000};}catch{return;}finally{solver.dispose();}
}
