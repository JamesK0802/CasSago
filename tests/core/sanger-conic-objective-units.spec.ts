import {expect,it} from 'vitest';
import {VerifiedConicSolver} from '../../src/app/workers/core/sanger-conic-solver';
import {OptimizationPending,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
// min x + (x-.3)^2 on [0,1]: x=0, objective=.09. The fixed-unit
// column has equality dual -.18 and x's lower-box dual is .4.
const program=():MarginalProgram=>({objectiveUnitScale:1e-4,upper:[1,1],cost:[1,0],rows:[{terms:[[1,1]],lower:1,upper:1}],quadraticCost:[[[0,1],[1,-.3]]]});
const answer=(args:any[],bad=false)=>{const scale=args[3][0];return{status:'unknown',x:[bad?-.01:0,1],z:[-.18*scale,.4*scale,0,0,0],iterations:1};};

it('restores nonzero duals and certifies the original objective after a scoped unit change',()=>{
 const p=program(),before=JSON.stringify(p),calls:any[][]=[];const s=new VerifiedConicSolver(p,{solve:(...args)=>{calls.push(args);const r=answer(args);return calls.length===1?{...r,z:r.z.map(()=>0)}:r;}});
 try{const r=s.solve();expect(r.value).toBeCloseTo(.09,12);expect(s.objectiveLowerBound()).toBeCloseTo(.09,12);expect(calls).toHaveLength(2);expect(calls[1][3][0]).toBe(1e-4);expect(calls[1][2][0]).toBe(2e-4);expect(JSON.parse(calls[1][11]).tol_gap_abs).toBeCloseTo(1e-14,20);expect(JSON.stringify(p)).toBe(before);}finally{s.dispose();}
});
it('rejects an infeasible proposal and recovers in equivalent units',()=>{
 const calls:any[][]=[],s=new VerifiedConicSolver(program(),{solve:(...args)=>{calls.push(args);return answer(args,calls.length===1);}});
 try{expect(s.solve().value).toBeCloseTo(.09,12);expect(calls.map(a=>a[3][0])).toEqual([1,1e-4]);expect(s.objectiveLowerBound()).toBeCloseTo(.09,12);}finally{s.dispose();}
});
it('does not repeatedly multiply the original tolerance on repeated solves',()=>{
 const tolerances:number[]=[],s=new VerifiedConicSolver(program(),{solve:(...args)=>{const r=answer(args);if(args[3][0]===1e-4){tolerances.push(JSON.parse(args[11]).tol_gap_abs);return r;}return{...r,z:r.z.map(()=>0)};}});
 try{s.solve();s.solve();expect(tolerances).toHaveLength(2);expect(tolerances[0]).toBeCloseTo(1e-14,20);expect(tolerances[1]).toBe(tolerances[0]);}finally{s.dispose();}
});
it('keeps a successful established solve as the only backend call',()=>{
 const costs:number[]=[],s=new VerifiedConicSolver(program(),{solve:(...args)=>{costs.push(args[3][0]);return answer(args);}});
 try{expect(s.solve().value).toBeCloseTo(.09,12);expect(costs).toEqual([1]);}finally{s.dispose();}
});
it('tries equivalent units on resuming a timed-out established solve',()=>{
 const costs:number[]=[],s=new VerifiedConicSolver(program(),{solve:(...args)=>{costs.push(args[3][0]);return costs.length===1?{status:'max_iterations',x:[],z:[],iterations:300}:answer(args);}});
 try{expect(()=>s.solve()).toThrow(OptimizationPending);expect(s.solve().value).toBeCloseTo(.09,12);expect(costs).toEqual([1,1e-4]);}finally{s.dispose();}
});
it('keeps linear inference and bound-query units unchanged',()=>{
 const costs:number[]=[],s=new VerifiedConicSolver({objectiveUnitScale:1e-4,upper:[1],cost:[1],rows:[]},{solve:(...args)=>{costs.push(args[3][0]);return{status:'optimal',x:[0],z:[1,0],iterations:1};}});
 try{expect(s.solve().value).toBe(0);expect(costs).toEqual([1]);}finally{s.dispose();}
});
it('never accepts a corrupted original primal in any objective units',()=>{
 const s=new VerifiedConicSolver(program(),{solve:(...args)=>answer(args,true)});
 try{expect(()=>s.solve()).toThrow('did not certify');expect(s.objectiveLowerBound()).toBeUndefined();}finally{s.dispose();}
});
it.each([0,-1,NaN,Infinity])('rejects invalid objective unit scale %s',scale=>{
 const p=program();p.objectiveUnitScale=scale;const s=new VerifiedConicSolver(p,{solve:()=>{throw Error('Must reject before optimization');}});
 try{expect(()=>s.solve()).toThrow('finite and positive');}finally{s.dispose();}
});
