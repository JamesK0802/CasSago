import {it,expect,vi} from 'vitest';
import loadCommon from '../probes/load-common-optimizer';
import {MarginalLinearSolver,marginalViolation,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
import {FactorizedMarginalSolver} from '../../src/app/workers/core/sanger-factorized-lp';
const program=():MarginalProgram=>({upper:[1,1],cost:[0,0],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1}],localBlocks:[[1]]});
it('completes exact decomposition after more than the former wall-clock selection budget',async()=>{
 const p=program(),before=JSON.stringify(p),realNow=performance.now.bind(performance),realSteps=FactorizedMarginalSolver.prototype.steps;let delayed=false,released=false;
 const clock=vi.spyOn(performance,'now').mockImplementation(()=>realNow()+(delayed?20000:0));
 const pending=vi.spyOn(FactorizedMarginalSolver.prototype,'steps').mockImplementation(function*(...args){try{delayed=true;yield;return yield* realSteps.apply(this,args);}finally{released=true;}});
 const solver=new MarginalLinearSolver(p,await loadCommon());
 try{const answer=solver.solve([[0,1]],true);expect(released).toBe(true);expect(answer.value).toBeCloseTo(1,7);expect(marginalViolation(p,answer.x)).toBeLessThanOrEqual(1e-7);expect(solver.objectiveLowerBound()).toBeCloseTo(-1,7);expect(solver.diagnostics.factorizationFallbackReason).toBeUndefined();expect(solver.diagnostics.selected).toBe('factorized');expect(JSON.stringify(p)).toBe(before);}
 finally{solver.dispose();pending.mockRestore();clock.mockRestore();}
});
it('retains independently certified unrestricted recovery after an actual decomposition failure',async()=>{
 const p=program(),before=JSON.stringify(p),failed=vi.spyOn(FactorizedMarginalSolver.prototype,'steps').mockImplementation(function*(){throw Error('Injected decomposition failure');}),solver=new MarginalLinearSolver(p,await loadCommon());
 try{const a=solver.solve([[0,1]],true);expect(a.value).toBeCloseTo(1,7);expect(marginalViolation(p,a.x)).toBeLessThanOrEqual(1e-7);expect(solver.objectiveLowerBound()).toBeCloseTo(-1,7);expect(solver.diagnostics.factorizationFallbackReason).toContain('Injected decomposition failure');expect(JSON.stringify(p)).toBe(before);}finally{solver.dispose();failed.mockRestore();}
});
