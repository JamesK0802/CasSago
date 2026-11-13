import {expect,it} from 'vitest';
import {VerifiedConicSolver} from '../../src/app/workers/core/sanger-conic-solver';

it.each([10,.1])('certifies a complete witness only after undoing equivalent objective scale %s',scale=>{
 const seen:{cost:number;tol:number}[]=[],solver=new VerifiedConicSolver({upper:[1],cost:[1],rows:[{terms:[[0,1]],lower:.25,upper:.25}]},{solve:(...a)=>{
  const cost=a[3][0];seen.push({cost,tol:JSON.parse(a[11]).tol_gap_abs});const z=Array(a[9]).fill(0);if(cost===scale)z[0]=-cost;return{status:'unknown',x:[.25],z,iterations:1};
 }});
 try{expect(solver.solve().value).toBe(.25);expect(solver.objectiveLowerBound()).toBeCloseTo(.25,10);expect(seen.slice(0,2)).toEqual([{cost:1,tol:1e-8},{cost:10,tol:1e-10}]);expect(seen.at(-1)?.cost).toBe(scale);expect(solver.stats.scalingRetries).toBe(scale===10?1:2);expect(solver.stats.maximumGap).toBeLessThanOrEqual(1e-7);}finally{solver.dispose();}
});
it('scales the quadratic Hessian and restores its original global certificate',()=>{
 const hessians:number[]=[],solver=new VerifiedConicSolver({upper:[1],cost:[0],quadraticCost:[[[0,1]]],rows:[{terms:[[0,1]],lower:.25,upper:.25}]},{solve:(...a)=>{
  const hessian=a[2][0];hessians.push(hessian);const z=Array(a[9]).fill(0);if(hessian===20)z[0]=-5;return{status:'unknown',x:[.25],z,iterations:1};
 }});
 try{expect(solver.solve().value).toBe(.0625);expect(hessians).toEqual([2,20]);expect(solver.objectiveLowerBound()).toBeCloseTo(.0625,10);}finally{solver.dispose();}
});
it('rejects every scale when finite witnesses do not certify the original objective',()=>{
 let calls=0;const solver=new VerifiedConicSolver({upper:[1],cost:[1],rows:[{terms:[[0,1]],lower:.25,upper:.25}]},{solve:(...a)=>{calls++;return{status:'unknown',x:[.25],z:Array(a[9]).fill(0),iterations:1};}});
 try{expect(()=>solver.solve()).toThrow('did not certify the original model');expect(calls).toBe(7);expect(solver.objectiveLowerBound()).toBeUndefined();}finally{solver.dispose();}
});
