import {it,expect} from 'vitest';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import initializeConic,* as conic from '../../public/vendor/clarabel-faer-0.11.1/clarabel.mjs';
import {readFileSync} from 'node:fs';
async function loadCommon(){await initializeConic({module_or_path:readFileSync('public/vendor/clarabel-faer-0.11.1/clarabel_wasm_bg.wasm')});const o=await loadHighs();o.conic=conic;return o;}

import {VerifiedConicSolver} from '../../src/app/workers/core/sanger-conic-solver';
import {OptimizationPending,MarginalLinearSolver,validMarginalSolution,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
it('uses a complete original primal/global-dual proof even when the backend stops early',async()=>{
 const optimizer=await loadCommon(),p:MarginalProgram={upper:[2],cost:[-1],rows:[],quadratics:[{name:'unit ball',terms:[[[0,1]]],upper:1}]};
 for(const status of ['unknown','max_iterations']){const solver=new VerifiedConicSolver(p,{solve:(...args)=>({...optimizer.conic!.solve(...args),status,iterations:300})});try{expect(solver.solve().value).toBeCloseTo(-1,6);expect(solver.objectiveLowerBound()).toBeCloseTo(-1,6);expect(solver.stats.statusIndependentCertificates).toBe(1);}finally{solver.dispose();}}
});
it('rejects uncertified early termination and resumes timed-out work without accepting a partial point',async()=>{
 const optimizer=await loadCommon(),p:MarginalProgram={upper:[2],cost:[-1],rows:[],quadratics:[{name:'unit ball',terms:[[[0,1]]],upper:1}]};
 for(const status of ['unknown','max_iterations']){const solver=new VerifiedConicSolver(p,{solve:(...args)=>{const r=optimizer.conic!.solve(...args);return{...r,status,iterations:300,z:Array(r.z.length).fill(0)};}});try{if(status==='max_iterations')expect(()=>solver.solve()).toThrow(OptimizationPending);else expect(()=>solver.solve()).toThrow('did not certify the original model');expect(solver.objectiveLowerBound()).toBeUndefined();expect(solver.stats.statusIndependentCertificates).toBe(0);}finally{solver.dispose();}}
});
it('tightens numerical precision after an uncertified complete early-stop witness',async()=>{
 const optimizer=await loadCommon(),settings:any[]=[],solver=new VerifiedConicSolver({upper:[2],cost:[-1],rows:[],quadratics:[{name:'unit ball',terms:[[[0,1]]],upper:1}]},{solve:(...args)=>{settings.push(JSON.parse(args[11]));const result=optimizer.conic!.solve(...args);return settings.length===1?{...result,status:'unknown',z:Array(result.z.length).fill(0)}:result;}});
 try{expect(solver.solve().value).toBeCloseTo(-1,6);expect(settings.map(s=>s.tol_gap_abs)).toEqual([1e-8,1e-10]);expect(solver.stats.certificationRetries).toBe(1);expect(solver.stats.maximumGap).toBeLessThanOrEqual(1e-7);}finally{solver.dispose();}
});
it('never accepts an early-stop witness with a corrupted original primal',async()=>{
 const optimizer=await loadCommon(),solver=new VerifiedConicSolver({upper:[2],cost:[-1],rows:[],quadratics:[{name:'unit ball',terms:[[[0,1]]],upper:1}]},{solve:(...args)=>({...optimizer.conic!.solve(...args),status:'unknown',x:[2]})});
 try{expect(()=>solver.solve()).toThrow('did not certify the original model');expect(solver.objectiveLowerBound()).toBeUndefined();expect(solver.stats.statusIndependentCertificates).toBe(0);}finally{solver.dispose();}
});
it.each([false,true])('uses an independent original-box certificate when redundant native multipliers are inaccurate: negative=%s',async negative=>{
 const optimizer=await loadCommon(),target=negative?2:0,p:MarginalProgram={upper:[2],cost:[negative?-3:1],rows:[{terms:[[0,1]],lower:target,upper:target}]};
 optimizer.conic={solve:(...args:any[])=>{const z=Array(args[9]).fill(0);z[0]=negative?10:-10;return{status:'unknown',x:[target],z,iterations:1};}};
 const s=new MarginalLinearSolver(p,optimizer,'conic');try{expect(s.solve().value).toBe(negative?-6:0);expect(s.objectiveLowerBound()).toBeCloseTo(negative?-6:0,10);expect(s.rowDuals()).toEqual([0]);expect(s.diagnostics.conic!.boxCertificates).toBe(1);expect(validMarginalSolution(p,s.solve().x)).toBe(true);}finally{s.dispose();}
});
it('does not accept a feasible suboptimal point on the basis of an unattained box lower bound',async()=>{
 const optimizer=await loadCommon();optimizer.conic={solve:(...args:any[])=>({status:'unknown',x:[.5],z:Array(args[9]).fill(0),iterations:1})};const s=new MarginalLinearSolver({upper:[1],cost:[1],rows:[]},optimizer,'conic');try{expect(()=>s.solve()).toThrow('did not certify');}finally{s.dispose();}
});
it('polishes an unknown feasible proposal until its original objective has a verified global certificate',async()=>{
 const optimizer=await loadCommon(),tolerances:number[]=[];optimizer.conic={solve:(...args:any[])=>{tolerances.push(JSON.parse(args[11]).tol_gap_abs);const z=Array(args[9]).fill(0);if(tolerances.length>1)z[0]=-args[3][0];return{status:'unknown',x:[.25],z,iterations:1};}};
 const s=new MarginalLinearSolver({upper:[1],cost:[1],rows:[{terms:[[0,1]],lower:.25,upper:.25}]},optimizer,'conic');try{expect(s.solve().value).toBe(.25);expect(s.objectiveLowerBound()).toBeCloseTo(.25,10);expect(tolerances).toEqual([1e-8,1e-10]);expect(s.diagnostics.conic!.certificationRetries).toBe(1);}finally{s.dispose();}
});
