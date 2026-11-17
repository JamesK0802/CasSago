import {it,expect} from 'vitest';
import {tightenBaseBudgetBounds} from '../../src/app/workers/core/sanger-base-budget-bounds';
it('preserves a stronger original departure-budget bound alongside the joint LP envelope',()=>{
 const estimate=[1,0,0,0],lower=[.70,0,0,0],upper=[1,.30,.30,.30];
 const budget=(terms:[number,number][])=>terms.some(([i])=>i===0)?undefined:.019;
 tightenBaseBudgetBounds([[[0,1]],[[1,1]],[[2,1]],[[3,1]]],.1,estimate,lower,upper,budget);
 expect(lower[0]).toBeGreaterThan(.80999);expect(lower[0]).toBeLessThan(.81001);expect(upper[1]).toBeLessThan(.19001);expect(estimate).toEqual([1,0,0,0]);
});
it('retains prior information and remains safe when the family may vanish or no budget certificate exists',()=>{
 const estimate=[.5,.5,0,0],lower=[.3,.2,0,0],upper=[.8,.7,0,0];
 tightenBaseBudgetBounds([[[0,1]],[[1,1]],[],[]],0,estimate,lower,upper,()=>undefined);
 for(let b=0;b<4;b++){expect(lower[b]).toBeCloseTo([.3,.2,0,0][b],12);expect(upper[b]).toBeCloseTo([.8,.7,0,0][b],12);}
});
