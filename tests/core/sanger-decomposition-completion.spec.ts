import {it,expect,vi} from 'vitest';
import {readFileSync} from 'node:fs';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import initializeConic,* as conic from '../../public/vendor/clarabel-cvxjs-0.1.4/clarabel.mjs';
import {MarginalLinearSolver,marginalViolation,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
import {FactorizedMarginalSolver} from '../../src/app/workers/core/sanger-factorized-lp';
it('releases an unfinished decomposition and certifies the identical unrestricted model',async()=>{
 await initializeConic({module_or_path:readFileSync('public/vendor/clarabel-cvxjs-0.1.4/clarabel_wasm_bg.wasm')});const optimizer=await loadHighs(),nativeSolve=conic.solve.bind(conic);let first=true,inDecomposition=false,released=false;
 optimizer.conic={solve:(...args:any[])=>{if(first){first=false;throw Error('Unfinished initial probe');}return(nativeSolve as any)(...args);}};
 const p:MarginalProgram={upper:[1,1],cost:[0,0],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1}],localBlocks:[[1]]},before=JSON.stringify(p),realNow=performance.now.bind(performance);
 const clock=vi.spyOn(performance,'now').mockImplementation(()=>realNow()+(inDecomposition?20000:0));
 const pending=vi.spyOn(FactorizedMarginalSolver.prototype,'steps').mockImplementation(function*(){inDecomposition=true;try{yield;}finally{inDecomposition=false;released=true;}return{x:[],value:0};});
 const solver=new MarginalLinearSolver(p,optimizer);
 try{const answer=solver.solve([[0,1]],true);expect(released).toBe(true);expect(answer.value).toBeCloseTo(1,7);expect(marginalViolation(p,answer.x)).toBeLessThanOrEqual(1e-7);expect(solver.objectiveLowerBound()).toBeCloseTo(-1,7);expect(solver.diagnostics.factorizationFallbackReason).toContain('numerical selection budget');expect(solver.diagnostics.selected).toBe('flat');expect(JSON.stringify(p)).toBe(before);}
 finally{solver.dispose();pending.mockRestore();clock.mockRestore();}
});
