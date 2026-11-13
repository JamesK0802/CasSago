import {it,expect} from 'vitest';
import {readFileSync} from 'node:fs';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import initializeConic,* as conic from '../../public/vendor/clarabel-faer-0.11.1/clarabel.mjs';
async function loadCommon(){await initializeConic({module_or_path:readFileSync('public/vendor/clarabel-faer-0.11.1/clarabel_wasm_bg.wasm')});const o=await loadHighs();o.conic=conic;return o;}

import {cachedLinearBound} from '../../src/app/workers/core/sanger-cached-linear-bound';
import {MarginalLinearSolver,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
it('reuses exact, scaled and nearby objectives with certified box corrections and invalidates changed constraints',async()=>{
 const p:MarginalProgram={immutableRows:true,upper:[1,1],cost:[0,0],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1}]},solver=new MarginalLinearSolver(p,await loadCommon()),cache=cachedLinearBound(p,solver);
 try{const first=cache.lower([[0,1]]);expect(first).toBeCloseTo(0,7);expect(cache.lower([[0,2]])).toBeCloseTo(0,7);expect(cache.lower([[0,1],[1,1e-10]])).toBeLessThanOrEqual(1e-10);expect(cache.stats.hits).toBe(2);expect(cache.stats.transfers).toBe(2);solver.setRowBounds([0],[.4],[.4]);expect(cache.lower([[0,1]],true)).toBeCloseTo(-.4,7);expect(cache.stats.solves).toBe(2);}finally{solver.dispose();}
});
it('rejects a feasible point without an independent lower certificate',()=>{
 const p:MarginalProgram={immutableRows:true,upper:[1],cost:[0],rows:[]},cache=cachedLinearBound(p,{solve:()=>({x:[1],value:1}),objectiveLowerBound:()=>undefined});
 expect(()=>cache.lower([[0,1]])).toThrow('independently certified');expect(cache.stats.hits).toBe(0);
});
it('exposes exact lower certificates for equality-fixed and zero objectives',async()=>{
 const p:MarginalProgram={immutableRows:true,upper:[1,1],cost:[0,0],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1}]},solver=new MarginalLinearSolver(p,await loadCommon());
 try{solver.solve([[0,1]]);solver.solve([[0,1],[1,1]]);expect(solver.objectiveLowerBound()).toBe(1);solver.solve([[0,1],[1,1]],true);expect(solver.objectiveLowerBound()).toBe(-1);solver.solve([]);expect(solver.objectiveLowerBound()).toBe(0);}finally{solver.dispose();}
});

it('does not transfer an unattained bound across a changed objective',async()=>{const p:MarginalProgram={immutableRows:true,upper:[1,1],cost:[0,0],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1}]},solver=new MarginalLinearSolver(p,await loadCommon()),cache=cachedLinearBound(p,solver);try{for(const [a,b]of [[1,0],[0,1],[-2,3],[.7,.2],[.7,.20000000001]]){const lower=cache.lower([[0,a],[1,b]]);expect(lower).toBeLessThanOrEqual(Math.min(a,b)+1e-10);expect(Math.min(a,b)-lower).toBeLessThanOrEqual(1e-7);}}finally{solver.dispose();}});
it('rejects malformed objective coordinates before any cached proof can hide them',()=>{const p:MarginalProgram={immutableRows:true,upper:[1],cost:[0],rows:[]},cache=cachedLinearBound(p,{solve:()=>{throw Error('must not call the backend');},objectiveLowerBound:()=>0});for(const terms of [[[2,1]],[[-1,1]],[[0,NaN]]]as[number,number][][])expect(()=>cache.lower(terms)).toThrow('Invalid cached objective');});
