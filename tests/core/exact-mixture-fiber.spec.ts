import {expect,it} from 'vitest';
import {exactHaplotypeRepresentatives,exactWeightedEquality,exactFromNumber} from '../../src/app/workers/core/exact-mixture-fiber';
const column=(s:string)=>s.split('').flatMap(c=>'ACGT'.split('').map(b=>Number(b===c)));
it('enumerates insertion haplotype alternatives without ground truth and certifies their complete signal',()=>{const base=[['AAA',.115],['ACC',.125],['CAC',.095],['CCA',.065]].map(([sequence,weight])=>({sequence:String(sequence),weight:Number(weight)})),r=exactHaplotypeRepresentatives(base,32);expect(r.options.length).toBeGreaterThan(1);for(const option of r.options)expect(exactWeightedEquality(base.map(x=>({weight:exactFromNumber(x.weight),column:column(x.sequence)})),option.map(x=>({weight:x.exactWeight,column:column(x.sequence)})))).toBe(true);});
it('does not certify even a one-bit numerical difference as identical',()=>{const a=[{weight:exactFromNumber(1),column:[1,0]}],b=[{weight:exactFromNumber(1),column:[1,Number.MIN_VALUE]}];expect(exactWeightedEquality(a,b)).toBe(false);});
it('retains exact marginals for long tracts without exhaustive haplotype enumeration',()=>{const base=[{sequence:'AAAAAAAAAA',weight:.3},{sequence:'CCCCCCCCCC',weight:.7}],r=exactHaplotypeRepresentatives(base);expect(r.options.length).toBeGreaterThan(1);for(const option of r.options)expect(exactWeightedEquality(base.map(x=>({weight:exactFromNumber(x.weight),column:column(x.sequence)})),option.map(x=>({weight:x.exactWeight,column:column(x.sequence)})))).toBe(true);});

import {exactBoundedMixture,exactNumber}from'../../src/app/workers/core/exact-mixture-fiber';
const bound=(lo:number,hi:number)=>({lower:exactFromNumber(lo),upper:exactFromNumber(hi)});
it('certifies a complete compatible profile even when it is not a listed representative',()=>{
 const source=[{weight:exactFromNumber(.5),column:column('AA')},{weight:exactFromNumber(.5),column:column('CC')}],cols=['AA','AC','CA','CC'].map(column),weights=exactBoundedMixture(cols,source,[bound(.2,.3),bound(.2,.3),bound(.2,.3),bound(.2,.3)]);
 expect(weights).not.toBeNull();expect(exactWeightedEquality(source,cols.map((c,i)=>({column:c,weight:weights![i]})))).toBe(true);expect(weights!.map(exactNumber).every(x=>x>=.2&&x<=.3)).toBe(true);
});
it('rejects a whole profile with incompatible fractions instead of matching its rows independently',()=>{
 expect(exactBoundedMixture(['AA','AC','CA','CC'].map(column),[{weight:exactFromNumber(.5),column:column('AA')},{weight:exactFromNumber(.5),column:column('CC')}],[bound(.4,.45),bound(.3,.35),bound(.05,.1),bound(.1,.15)])).toBeNull();
});
it('never turns numerical simplex tolerance into a one-bit equality certificate',()=>{
 expect(exactBoundedMixture([[1,Number.MIN_VALUE]],[{weight:exactFromNumber(1),column:[1,0]}],[bound(1,1)])).toBeNull();
});
it('uses rational feasibility when a floating vertex violates a much smaller positive bound',()=>{const cols=[[1,0],[1,1],[1,2]],source=[{weight:exactFromNumber(1),column:[1,1e-20]}],weights=exactBoundedMixture(cols,source,[bound(0,1),bound(0,1e-30),bound(0,1)]);expect(weights).not.toBeNull();expect(exactNumber(weights![1])).toBeLessThanOrEqual(1e-30);expect(exactWeightedEquality(source,cols.map((column,i)=>({column,weight:weights![i]})))).toBe(true);});
