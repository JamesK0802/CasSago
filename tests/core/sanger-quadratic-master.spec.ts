import {it,expect} from 'vitest';
import loadCommon from '../probes/load-common-optimizer';
import {MarginalLinearSolver,validMarginalSolution,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
it('preserves an analytic curved phase bound while factoring local base allocations',async()=>{
 const optimizer=await loadCommon(),p:MarginalProgram={upper:[1,1,1,1],cost:[0,0,0,0],localBlocks:[[2,3]],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1},{terms:[[2,1],[3,1],[0,-1]],lower:0,upper:0}],quadratics:[{name:'phase measurement',terms:[[[0,1]],[[1,1]]],upper:.68}]},snapshot=JSON.stringify(p),solver=new MarginalLinearSolver(p,optimizer);
 try{for(const maximum of [false,true]){const answer=solver.solve([[0,1]],maximum);expect(answer.value).toBeCloseTo(maximum?.8:.2,6);expect(validMarginalSolution(p,answer.x)).toBe(true);}expect(JSON.stringify(p)).toBe(snapshot);}finally{solver.dispose();}
});
