import {it,expect} from 'vitest';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import {MarginalLinearSolver,validMarginalSolution,type MarginalProgram,type LinearTerm} from '../../src/app/workers/core/sanger-marginal-lp';
it('prices omitted origins back in and matches the full model for changing objectives beyond ten origins',async()=>{
 const n=160,p:MarginalProgram={upper:Array(n).fill(1),cost:Array.from({length:n},(_,i)=>i===159?-1:0),rows:[{terms:Array.from({length:n},(_,i)=>[i,1]),lower:1,upper:1}],originPricing:{columns:Array.from({length:n},(_,i)=>i),seed:[0]}},optimizer=await loadHighs(),restricted=new MarginalLinearSolver(structuredClone(p),optimizer),full=new MarginalLinearSolver({...p,originPricing:undefined},optimizer);
 try{expect(restricted.solve().value).toBe(-1);expect(restricted.pricingStats.activated).toBeGreaterThan(0);for(let k=0;k<12;k++){const form=Array.from({length:n},(_,i)=>[i,Math.sin(17*i+3*k)]as LinearTerm);for(const maximize of [false,true]){const answer=restricted.solve(form,maximize);expect(validMarginalSolution(p,answer.x)).toBe(true);expect(answer.value).toBeCloseTo(full.solve(form,maximize).value,7);}}}finally{restricted.dispose();full.dispose();}
});
it('restores the original columns before treating a restricted infeasibility as biological infeasibility',async()=>{
 const p:MarginalProgram={upper:[1,1],cost:[0,0],rows:[{terms:[[1,1]],lower:1,upper:1}],originPricing:{columns:[0,1],seed:[0]}},solver=new MarginalLinearSolver(p,await loadHighs());try{expect(solver.solve().x[1]).toBe(1);expect(solver.pricingStats.fullRestores).toBe(1);}finally{solver.dispose();}
});
