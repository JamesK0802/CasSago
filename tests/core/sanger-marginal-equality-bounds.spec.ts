import {expect,it} from 'vitest';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import {MarginalLinearSolver} from '../../src/app/workers/core/sanger-marginal-lp';
it('reuses only algebraically fixed objectives with a verified feasible witness',async()=>{const s=new MarginalLinearSolver({upper:[1,1,0],cost:[0,0,0],rows:[{terms:[[1,1],[0,1]],lower:.4,upper:.4}]},await loadHighs());try{s.solve();expect(s.solve([[0,1],[1,.5],[1,.5]]).value).toBe(.4);expect(s.solve([[2,1]],true).value).toBe(0);expect(s.bounds([[0,1]])[1]).toBeCloseTo(.4,6);s.setRowBounds([0],[.2],[.2]);expect(s.bounds([[0,1],[1,1]])[0]).toBeCloseTo(.2,6);s.addRow([[0,1],[1,1]],.7,.7);expect(()=>s.solve([[0,1],[1,1]])).toThrow();}finally{s.dispose();}});
