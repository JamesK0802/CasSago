import{it,expect}from'vitest';
import{attainedBoxObjective}from'../../src/app/workers/core/sanger-attained-box';
import{MarginalLinearSolver,marginalViolation}from'../../src/app/workers/core/sanger-marginal-lp';
import{readFileSync}from'node:fs';
import loadHighs from'../../public/vendor/highs-1.15.4/highs.mjs';
import initialize,*as conic from'../../public/vendor/clarabel-faer-0.11.1/clarabel.mjs';
async function optimizer(){await initialize({module_or_path:readFileSync('public/vendor/clarabel-faer-0.11.1/clarabel_wasm_bg.wasm')});const engine=await loadHighs();engine.conic=conic;return engine;}
it('certifies an attained nonnegative minimum and a mixed-sign box maximum',()=>{
 const a=attainedBoxObjective([2,3],[0,2],[[0,1]]);expect(a?.value).toBe(0);expect(a?.lower).toBeLessThanOrEqual(0);
 const b=attainedBoxObjective([2,3],[2,0],[[0,2],[1,-1]],true);expect(b?.value).toBe(4);expect(b?.lower).toBeLessThanOrEqual(-4);expect(-4-b!.lower).toBeLessThan(1e-7);
});
it('coalesces repeated columns before constructing a global box bound',()=>{const a=attainedBoxObjective([1],[.4],[[0,1],[0,-1]]);expect(a?.value).toBe(0);expect(a?.lower).toBeLessThanOrEqual(0);});
it('does not certify a feasible but suboptimal point',()=>{expect(attainedBoxObjective([1],[.25],[[0,1]])).toBeUndefined();expect(attainedBoxObjective([1],[.25],[[0,1]],true)).toBeUndefined();});
it('rejects malformed and nonfinite objective data',()=>{expect(attainedBoxObjective([1],[NaN],[[0,1]])).toBeUndefined();expect(attainedBoxObjective([Infinity],[1],[[0,1]])).toBeUndefined();expect(attainedBoxObjective([1],[0],[[2,1]])).toBeUndefined();expect(attainedBoxObjective([1],[0],[[0,Infinity]])).toBeUndefined();});
it('reuses only the current analysis witness and exports a fresh independent box certificate',async()=>{
 const o=await optimizer(),backend=o.conic;let calls=0;o.conic={solve:(...a:any[])=>{calls++;return backend.solve(...a);}};const s=new MarginalLinearSolver({upper:[1,1],cost:[0,0],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1}],quadratics:[{name:'ball',terms:[[[0,1]],[[1,1]]],upper:1}]},o,'conic');
 try{s.solve([[0,1]],true);const before=calls,a=s.solve([[1,1]]);expect(calls).toBe(before);expect(a.value).toBeLessThanOrEqual(1e-7);expect(s.objectiveLowerBound()).toBeLessThanOrEqual(a.value);expect(a.value-s.objectiveLowerBound()!).toBeLessThanOrEqual(1e-7);expect(s.rowDuals()).toEqual([0]);expect(s.supportPlanes()).toEqual([]);expect(s.diagnostics.attainedBoxHits).toBe(1);expect(marginalViolation(s.program,a.x)).toBeLessThanOrEqual(1e-7);}finally{s.dispose();}
});
it('checks the original nonlinear constraints before taking a box shortcut',async()=>{
 const o=await optimizer(),backend=o.conic;let calls=0;o.conic={solve:(...a:any[])=>{calls++;return backend.solve(...a);}};const s=new MarginalLinearSolver({upper:[1,1],cost:[0,0],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1}],quadratics:[{name:'original ball',terms:[[[0,1]],[[1,1]]],upper:.6}]},o,'conic');
 try{s.solve([[0,1]],true);(s as any).point=[1,0];const before=calls,a=s.solve([[1,1]]);expect(calls).toBeGreaterThan(before);expect(a.value).toBeGreaterThan(.2);expect(marginalViolation(s.program,a.x)).toBeLessThanOrEqual(1e-7);expect(s.diagnostics.attainedBoxHits).toBe(0);}finally{s.dispose();}
});
it('rejects a stale witness after a change to an original row',async()=>{
 const o=await optimizer(),backend=o.conic;let calls=0;o.conic={solve:(...a:any[])=>{calls++;return backend.solve(...a);}};const s=new MarginalLinearSolver({upper:[1,1],cost:[0,0],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1}]},o,'conic');
 try{s.solve([[0,1]],true);s.setRowBounds([0],[.5],[.5]);(s as any).point=[1,0];const before=calls,a=s.solve([[1,1]]);expect(calls).toBeGreaterThan(before);expect(marginalViolation(s.program,a.x)).toBeLessThanOrEqual(1e-7);expect(s.diagnostics.attainedBoxHits).toBe(0);}finally{s.dispose();}
});
