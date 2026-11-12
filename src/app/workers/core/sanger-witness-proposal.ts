import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import {FactorizedMarginalSolver} from './sanger-factorized-lp';
import {MarginalLinearSolver,type MarginalProgram} from './sanger-marginal-lp';
/** Bounded work for an OPTIONAL exact witness, never for the parent model.
 * A failed/incomplete proposal leaves its unrestricted solve unchanged. */
export function tryMarginalWitness(p:MarginalProgram,optimizer:SangerLinearOptimizer,seconds=10):{x:number[];value:number}|undefined{
 if(p.localBlocks?.some(b=>b.length)){
  const solver=new FactorizedMarginalSolver(p,q=>{const flat=new MarginalLinearSolver({...q,localBlocks:undefined},optimizer);flat.setTimeSlice(.25);return flat;}),steps=solver.steps(),start=performance.now();
  try{while((performance.now()-start)/1000<seconds){const next=steps.next();if(next.done)return next.value;}}catch{/* An incomplete proposal is not a parent-model infeasibility decision. */}finally{steps.return(undefined as any);solver.dispose();}
  return;
 }
 let solver:MarginalLinearSolver|undefined;try{solver=new MarginalLinearSolver(p,optimizer);solver.setTimeSlice(seconds);return solver.solve();}catch{return;}finally{solver?.dispose();}
}
