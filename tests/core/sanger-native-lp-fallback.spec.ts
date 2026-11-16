import {it,expect} from 'vitest';
import loadCommon from '../probes/load-common-optimizer';
import {MarginalLinearSolver,validMarginalSolution,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
it('recovers the identical small LP after a native runtime failure with an original primal/dual proof',async()=>{
 const optimizer=await loadCommon(),create=optimizer.createModel.bind(optimizer);let failures=0;
 optimizer.createModel=()=>{const native=create();return new Proxy(native,{get(target,key){if(key==='run')return()=>{failures++;throw Error('Injected native runtime failure');};const value=target[key];return typeof value==='function'?value.bind(target):value;}});};
 const p:MarginalProgram={upper:[1,1],cost:[0,0],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1},{terms:[[0,1]],lower:.2,upper:.6}]},before=JSON.stringify(p),solver=new MarginalLinearSolver(p,optimizer);
 try{const lo=solver.solve([[0,1]]),hi=solver.solve([[0,1]],true);expect(lo.value).toBeCloseTo(.2,7);expect(hi.value).toBeCloseTo(.6,7);expect(validMarginalSolution(p,lo.x)).toBe(true);expect(validMarginalSolution(p,hi.x)).toBe(true);expect(solver.objectiveLowerBound()).toBeCloseTo(-.6,7);expect(failures).toBe(1);expect(JSON.stringify(p)).toBe(before);}finally{solver.dispose();}
});
