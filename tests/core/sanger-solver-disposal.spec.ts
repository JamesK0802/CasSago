import {it,expect} from 'vitest';
import loadCommonOptimizer from '../probes/load-common-optimizer';
import {MarginalLinearSolver} from '../../src/app/workers/core/sanger-marginal-lp';
it('drops completed decomposition and witness storage and disposes native handles only once',async()=>{
 const optimizer=await loadCommonOptimizer(),create=optimizer.createModel.bind(optimizer);optimizer.conic=undefined;let created=0,disposed=0;
 optimizer.createModel=()=>{const model=create();created++;return new Proxy(model,{get(t,k){if(k==='dispose')return()=>{disposed++;return t.dispose();};const v=t[k];return typeof v==='function'?v.bind(t):v;}});};
 const solver=new MarginalLinearSolver({upper:[1,1],cost:[0,0],rows:[{terms:[[0,1],[1,-1]],lower:0,upper:0},{terms:[[0,1]],lower:.4,upper:.4}],localBlocks:[[1]]},optimizer);
 expect((solver as any).factorized).toBeUndefined();expect(solver.solve([[1,1]]).value).toBeCloseTo(.4,7);
 solver.addRow([[1,1]],.4,.5);expect((solver as any).factorized).toBeUndefined();expect(solver.solve([[1,1]]).value).toBeCloseTo(.4,7);solver.dispose();solver.dispose();
 expect(created).toBeGreaterThan(0);expect(disposed).toBe(created);expect((solver as any).factorized).toBeUndefined();expect((solver as any).fullCompetitor).toBeUndefined();expect((solver as any).point).toBeUndefined();
});
it('releases conic certificate arrays after disposal',async()=>{
 const solver=new MarginalLinearSolver({upper:[1],cost:[1],rows:[],quadratics:[{name:'bound',terms:[[[0,1]]],upper:1}]},await loadCommonOptimizer());solver.solve();expect(solver.supportPlanes().length).toBeGreaterThan(0);solver.dispose();solver.dispose();expect(solver.supportPlanes()).toEqual([]);expect(solver.rowDuals()).toEqual([]);
});
