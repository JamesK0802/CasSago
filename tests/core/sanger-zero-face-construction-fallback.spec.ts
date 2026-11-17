import {it,expect} from 'vitest';
import loadCommon from '../probes/load-common-optimizer';
import {verifiedZeroDeparturePoint} from '../../src/app/workers/core/sanger-initial-marginal-fit';
import {MarginalLinearSolver,validMarginalSolution,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
const fixture=():MarginalProgram=>({upper:[1,1,1],cost:[0,0,1],rows:[{terms:[[0,1]],lower:1,upper:1},{terms:[[1,1],[2,1]],lower:1,upper:1},{terms:[[1,1],[0,-1]],lower:-1e30,upper:0},{terms:[[2,1]],lower:.05,upper:1}],zeroCostSubstitutions:[{column:1,terms:[[0,1]]}]});
it('keeps the original positive-cost solve when constructing its zero-cost projection fails',async()=>{
 const p=fixture(),before=JSON.stringify(p),o=await loadCommon(),factory=o.createModel!;let created=0,failedDisposed=0;
 o.createModel=()=>{const m=factory.call(o);if(++created===2){m.passModel=()=>{throw Error('Projected model construction failed');};const dispose=m.dispose.bind(m);m.dispose=()=>{failedDisposed++;dispose();};}return m;};
 const solver=new MarginalLinearSolver(p,o);try{const answer=solver.solveToCompletion();expect(answer.value).toBeCloseTo(.05,8);expect(validMarginalSolution(p,answer.x)).toBe(true);expect(solver.objectiveLowerBound()).toBeCloseTo(.05,8);expect(solver.diagnostics.zeroFaceProjection?.failure).toContain('Projected model construction failed');expect(failedDisposed).toBe(1);expect(JSON.stringify(p)).toBe(before);}finally{solver.dispose();}
});
it('returns to the full primary solve if the optional zero-departure model cannot initialize',async()=>{
 const p=fixture(),before=JSON.stringify(p),o=await loadCommon(),factory=o.createModel!;let disposed=0;
 o.createModel=()=>{const m=factory.call(o);m.passModel=()=>{throw Error('Optional zero model unavailable');};const close=m.dispose.bind(m);m.dispose=()=>{disposed++;close();};return m;};
 expect(()=>verifiedZeroDeparturePoint(p,[{column:0,value:1}],o)).not.toThrow();expect(verifiedZeroDeparturePoint(p,[{column:0,value:1}],o)).toBeUndefined();expect(disposed).toBe(2);expect(JSON.stringify(p)).toBe(before);
});
