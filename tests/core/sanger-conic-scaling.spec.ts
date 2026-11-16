import {expect,it} from 'vitest';
import loadCommon from '../probes/load-common-optimizer';
import {VerifiedConicSolver} from '../../src/app/workers/core/sanger-conic-solver';
import {marginalViolation,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
for(const quadratic of [false,true])it('certifies original '+(quadratic?'quadratic':'linear')+' objective after an equivalent scale retry',async()=>{
 const o=await loadCommon(),p:MarginalProgram={upper:[1],cost:[-.2],rows:[],quadraticCost:quadratic?[[[0,1]]]:undefined},before=JSON.stringify(p);let calls=0;
 const runtime={solve:(...args:any[])=>++calls===1?{status:'unknown'}:o.conic.solve(...args)},s=new VerifiedConicSolver(p,runtime);
 try{const a=s.solve();expect(a.value).toBeCloseTo(quadratic?-.01:-.2,7);expect(a.x[0]).toBeCloseTo(quadratic?.1:1,6);expect(marginalViolation(p,a.x)).toBeLessThanOrEqual(1e-7);expect(Math.abs(a.value-s.objectiveLowerBound()!)).toBeLessThanOrEqual(1e-7);expect(s.stats.scalingRetries).toBe(1);expect(JSON.stringify(p)).toBe(before);}finally{s.dispose();}
});
it('still rejects a corrupted primal after equivalent scale retries',async()=>{const o=await loadCommon();let calls=0;const s=new VerifiedConicSolver({upper:[1],cost:[-1],rows:[]},{solve:(...args:any[])=>{if(++calls===1)return{status:'unknown'};return{...o.conic.solve(...args),x:[2]};}});try{expect(()=>s.solve()).toThrow('did not certify the original model');}finally{s.dispose();}});
it('limits unavailable numerical representations to the existing fallback',()=>{let calls=0;const s=new VerifiedConicSolver({upper:[1],cost:[-1],rows:[]},{solve:()=>{calls++;return{status:'unknown'};}});try{expect(()=>s.solve()).toThrow('no dual witness');expect(calls).toBe(3);}finally{s.dispose();}});
