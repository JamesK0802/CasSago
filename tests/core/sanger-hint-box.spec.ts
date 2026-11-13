import {it,expect} from 'vitest';
import {prepareHintBox} from '../../src/app/workers/core/sanger-hint-box';
it('preserves full dense-hint box minimization for sparse costs and both coefficient signs',()=>{
 const upper=Array.from({length:1000},(_,i)=>.1+i%17),image=new Map(upper.map((_,i)=>[i,Math.sin(i*7)])),geometric=new Map(upper.flatMap((_,i)=>i%11===0?[[i,i]]:[])),shift=prepareHintBox(image,upper,geometric);
 for(const alpha of [0,1e-8,.13,1,30])for(let k=0;k<20;k++){const reduced=new Map(upper.flatMap((_,i)=>i%31===k?[[i,Math.cos(i*9+k)]]:[])),actual=shift(reduced,alpha);let expected=0;for(let i=0;i<upper.length;i++){const cost=(reduced.get(i)||0)-alpha*(image.get(i)||0);if(geometric.has(i))expect(actual.coefficients.get(i)||0).toBeCloseTo(cost,10);else expected+=Math.min(0,cost)*upper[i];}expect(Math.abs(actual.constant-expected)).toBeLessThan(1e-9*(1+Math.abs(expected)));expect([...actual.coefficients.keys()].every(i=>geometric.has(i))).toBe(true);}
});

it('rejects invalid multipliers rather than reversing the box-minimization inequality',()=>{const shift=prepareHintBox(new Map([[0,1]]),[1],new Map());for(const alpha of [-1,Infinity,NaN])expect(()=>shift(new Map(),alpha)).toThrow('nonnegative finite');});
