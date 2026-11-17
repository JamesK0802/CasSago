import {it,expect} from 'vitest';
import {VerifiedConicSolver} from '../../src/app/workers/core/sanger-conic-solver';
it('uses equivalent bridge units only after rejected original certificates and certifies the unscaled objective',()=>{
 const scales:number[]=[],solver=new VerifiedConicSolver({upper:[2],cost:[-1],rows:[],quadraticCost:[[[0,1]]],objectiveUnitScale:1e-4},{solve:(...args)=>{const scale=-args[3][0];scales.push(scale);return{status:'unknown',x:[Math.abs(scale-.01)<1e-15?.5:.4],z:Array(args[9]).fill(0),iterations:1};}});
 try{const answer=solver.solve();expect(scales).toEqual([1,.0001,.01]);expect(answer.x[0]).toBe(.5);expect(answer.value).toBeCloseTo(-.25,12);expect(solver.objectiveLowerBound()).toBeCloseTo(-.25,12);expect(solver.stats.maximumGap).toBeLessThanOrEqual(1e-7);expect(solver.stats.maximumViolation).toBeLessThanOrEqual(1e-7);}finally{solver.dispose();}
});
