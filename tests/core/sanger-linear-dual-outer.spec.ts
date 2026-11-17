import{it,expect}from'vitest';
import loadHighs from'../../public/vendor/highs-1.15.4/highs.mjs';
import{linearDualOuterLower,cachedLinearDualOuter}from'../../src/app/workers/core/sanger-linear-dual-outer';
import{LP_INFINITY,MarginalLinearSolver,type MarginalProgram}from'../../src/app/workers/core/sanger-marginal-lp';
it('retains a valid original box/row bound for arbitrary inaccurate multipliers and signed costs',()=>{
 const p:MarginalProgram={upper:[1,1],cost:[0,0],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1},{terms:[[0,-1],[1,.2]],lower:-LP_INFINITY,upper:.1}]};
 for(let k=0;k<80;k++){const cost=[Math.sin(k),Math.cos(3*k)],dual=[1000*Math.sin(7*k),1000*Math.cos(k)];const lower=linearDualOuterLower(p,cost,dual)!;
  for(let j=0;j<=100;j++){const x=j/100,y=1-x;if(-x+.2*y>.1)continue;expect(lower).toBeLessThanOrEqual(cost[0]*x+cost[1]*y);}
 }
});
it('clamps multipliers pointing towards missing row sides, never towards a fictitious bound',()=>{
 const p:MarginalProgram={upper:[1],cost:[-1],rows:[{terms:[[0,1]],lower:.3,upper:LP_INFINITY}]};
 expect(linearDualOuterLower(p,[-1],[-1000])).toBeLessThanOrEqual(-1);expect(linearDualOuterLower(p,[-1],[1])).toBeLessThanOrEqual(-1);
});
it('rejects nonfinite witnesses, malformed physical coefficients, unbounded boxes and nonlinear input',()=>{
 const p:MarginalProgram={upper:[1],cost:[0],rows:[{terms:[[0,1]],lower:0,upper:1}]};
 expect(linearDualOuterLower(p,[0],[NaN])).toBeUndefined();expect(linearDualOuterLower({...p,upper:[LP_INFINITY]},[-1],[0])).toBeUndefined();expect(linearDualOuterLower({...p,rows:[{...p.rows[0],terms:[[0,NaN]]}]},[0],[0])).toBeUndefined();expect(linearDualOuterLower({...p,quadratics:[{name:'not linear',terms:[[[0,1]]],upper:1}]},[0],[0])).toBeUndefined();
});
it('returns only a safe outer proof under an unfinished native proposal, caches exact objectives and invalidates changed boxes',()=>{
 const p:MarginalProgram={immutableRows:true,upper:[1],cost:[0],rows:[{terms:[[0,1]],lower:.2,upper:1}]};let calls=0;
 const fake:any={createModel:()=>({options:{set:()=>{}},passModel:()=>{},changeColsCost:()=>{},zeroAllClocks:()=>{},run:()=>{calls++;return{modelStatus:13};},getSolution:()=>({colValue:[NaN],rowDual:[123]}),dispose:()=>{}})};
 const cache=cachedLinearDualOuter(p,fake);try{expect(cache.lower([[0,1]])).toBeLessThanOrEqual(.2);expect(cache.lower([[0,.5],[0,.5]])).toBeLessThanOrEqual(.2);expect(calls).toBe(1);expect(cache.stats.partial).toBe(1);p.upper[0]=.5;expect(cache.lower([[0,1]])).toBeLessThanOrEqual(.2);expect(calls).toBe(2);}finally{cache.dispose();}
});
it('falls back to an original box proof when the native runtime or its dual is unusable',()=>{
 const p:MarginalProgram={upper:[2],cost:[0],rows:[]},cache=cachedLinearDualOuter(p,{createModel:()=>{throw Error('unavailable');}}as any);
 try{expect(cache.lower([[0,-3]])).toBeLessThanOrEqual(-6);expect(cache.stats.failures).toBe(1);expect(()=>cache.lower([[1,1]])).toThrow('Invalid outer objective');}finally{cache.dispose();}
});
it('bounds actual common-WASM LP optima without requiring its proposed primal to attain the bound',async()=>{
 const p:MarginalProgram={upper:[1,1,1],cost:[0,0,0],rows:[{terms:[[0,1],[1,1],[2,1]],lower:1,upper:1},{terms:[[0,-1],[1,2]],lower:-.2,upper:.3}]},optimizer=await loadHighs(),full=new MarginalLinearSolver(p,optimizer),outer=cachedLinearDualOuter(p,optimizer,.01);
 try{for(let k=0;k<20;k++){const terms:[number,number][]=Array.from({length:3},(_,i)=>[i,Math.sin(3*k+i)]);for(const maximum of[false,true]){const optimum=full.solve(terms,maximum).value,lower=outer.lower(terms,maximum);expect(lower).toBeLessThanOrEqual((maximum?-1:1)*optimum+1e-12);expect(Number.isFinite(lower)).toBe(true);}}}finally{outer.dispose();full.dispose();}
});

it('retains the box proof even if an unusable native workspace also throws during disposal',()=>{
 const p:MarginalProgram={upper:[2],cost:[0],rows:[]},cache=cachedLinearDualOuter(p,{createModel:()=>({options:{set:()=>{}},passModel:()=>{},zeroAllClocks:()=>{},run:()=>{throw Error('unusable numerical workspace');},dispose:()=>{throw Error('already unusable');}})}as any);
 try{expect(cache.lower([[0,-3]])).toBeLessThanOrEqual(-6);expect(cache.stats.failures).toBe(1);}finally{expect(()=>cache.dispose()).not.toThrow();}
});
it('proves a zero objective directly without solving the otherwise large feasibility problem',()=>{
 const p:MarginalProgram={immutableRows:true,upper:[1,1],cost:[0,0],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1}]},cache=cachedLinearDualOuter(p,{createModel:()=>{throw Error('zero needs no optimizer');}}as any);
 try{expect(cache.lower([[0,1],[0,-1]])).toBe(0);expect(cache.lower([])).toBe(0);expect(cache.stats.solves).toBe(0);expect(cache.stats.failures).toBe(0);}finally{cache.dispose();}
});
it('never certifies a falsely zero objective after catastrophic cancellation of duplicate coefficients',()=>{
 const p:MarginalProgram={immutableRows:true,upper:[1],cost:[0],rows:[]},cache=cachedLinearDualOuter(p,{createModel:()=>{throw Error('rounded zero uses its error-corrected box proof');}}as any);
 try{const terms:[number,number][]=[[0,-1e16],[0,-1],[0,1e16]];expect(cache.lower(terms)).toBeLessThanOrEqual(-1);expect(cache.lower([])).toBe(0);expect(cache.stats.solves).toBe(0);}finally{cache.dispose();}
});

it('prepared immutable certificates are bit-identical to complete validation for changing signed objectives',()=>{
 const p:MarginalProgram={immutableRows:true,upper:[1,2,.8],cost:[0,0,0],rows:[{terms:[[0,1],[1,.5],[2,-1]],lower:-.3,upper:.4},{terms:[[0,-.7],[1,.2],[2,1]],lower:0,upper:LP_INFINITY}]};let dual=[0,0];
 const fake:any={createModel:()=>({options:{set:()=>{}},passModel:()=>{},changeColsCost:()=>{},zeroAllClocks:()=>{},run:()=>({modelStatus:7}),getSolution:()=>({rowDual:dual}),dispose:()=>{}})},cache=cachedLinearDualOuter(p,fake);
 try{for(let k=1;k<=40;k++){dual=[3*Math.sin(k),2*Math.cos(2*k)];const terms:[number,number][]=[[0,Math.sin(3*k)],[1,Math.cos(k)],[2,Math.sin(5*k)]];for(const maximum of[false,true]){const cost=terms.map(([,w])=>(maximum?-1:1)*w),expected=Math.max(linearDualOuterLower(p,cost,[0,0])!,linearDualOuterLower(p,cost,dual)!);expect(cache.lower(terms,maximum)).toBe(expected);}}}finally{cache.dispose();}
});

it('invalidates prepared structural validation when numeric boxes or row bounds change',()=>{
 const p:MarginalProgram={immutableRows:true,upper:[1],cost:[0],rows:[{terms:[[0,1]],lower:0,upper:1}]},fake:any={createModel:()=>({options:{set:()=>{}},passModel:()=>{},changeColsCost:()=>{},zeroAllClocks:()=>{},run:()=>({modelStatus:7}),getSolution:()=>({rowDual:[0]}),dispose:()=>{}})},cache=cachedLinearDualOuter(p,fake);
 try{expect(Number.isFinite(cache.lower([[0,-1]]))).toBe(true);p.upper[0]=NaN;expect(()=>cache.lower([[0,-1]])).toThrow();p.upper[0]=1;p.rows[0].lower=2;expect(()=>cache.lower([[0,-1]])).toThrow();}finally{cache.dispose();}
});
