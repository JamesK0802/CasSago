import{it,expect}from'vitest';
import loadCommon from'../probes/load-common-optimizer';
import{boundMarginalSubstitutionMass,boundMarginalSubstitutions}from'../../src/app/workers/core/sanger-marginal-substitution';
import{zeroLengthCompositionBounds}from'../../src/app/workers/core/sanger-zero-length-composition';
const ref:any={name:'R',reverse:false,cuts:[100]};
it('retains a mandatory Sub lower bound and the complementary WT upper bound',async()=>{
 const p={upper:[1,.2],cost:[0,1],rows:[{terms:[[0,1]],lower:1,upper:1},{terms:[[1,1]],lower:.2,upper:.2}],dualHint:{budgetRow:1,shadow:[]}},g={id:'zero',reference:0,target:0,size:0,mass:[[0,1]],departures:new Map([[100,[[1,1]]]])},target:any={program:p,groups:[g],point:[1,.2]},before=JSON.stringify(p),o=await loadCommon(),r=boundMarginalSubstitutionMass(target,[ref],[1],o,()=>{},new Set(['zero'])),bounds=zeroLengthCompositionBounds([1,1],.2,1,[r.lower[0],r.upper[0]]);
 expect(bounds.substitutionRange[0]).toBeCloseTo(.2,5);expect(bounds.substitutionRange[1]).toBeCloseTo(.2,5);expect(bounds.unmodifiedRange[0]).toBeCloseTo(.8,5);expect(bounds.unmodifiedRange[1]).toBeCloseTo(.8,5);expect(JSON.stringify(p)).toBe(before);
});
it('keeps genuine multi-site molecular-linkage ambiguity instead of forcing WT largest',async()=>{
 const p={upper:[1,.8,.8],cost:[0,1,1],rows:[{terms:[[0,1]],lower:1,upper:1},{terms:[[1,1]],lower:.8,upper:.8},{terms:[[2,1]],lower:.8,upper:.8}]},g={id:'zero',reference:0,target:0,size:0,mass:[[0,1]],departures:new Map([[100,[[1,1]]],[101,[[2,1]]]])},target:any={program:p,groups:[g],point:[1,.8,.8]},r=boundMarginalSubstitutionMass(target,[ref],[1],await loadCommon(),()=>{},new Set(['zero'])),b=zeroLengthCompositionBounds([1,1],1.6,1,[r.lower[0],r.upper[0]]);
 // Both changed sites can coincide on80% of molecules or cover all molecules.
 expect(b.substitutionRange[0]).toBeCloseTo(.8,5);expect(b.substitutionRange[1]).toBeCloseTo(1,5);expect(b.unmodifiedRange[0]).toBe(0);expect(b.unmodifiedRange[1]).toBeCloseTo(.2,5);
});
it('does not transfer indel-carried Sub into zero-length composition or replace the global metric',async()=>{
 const p={upper:[.6,.4,.1],cost:[0,0,1],rows:[{terms:[[0,1]],lower:.6,upper:.6},{terms:[[1,1]],lower:.4,upper:.4},{terms:[[2,1]],lower:.1,upper:.1}]},zero={id:'zero',reference:0,target:0,size:0,mass:[[0,1]],departures:new Map()},ins={id:'insert',reference:0,target:0,size:1,mass:[[1,1]],departures:new Map([[100,[[2,1]]]])},target:any={program:p,groups:[zero,ins],point:[.6,.4,.1]},o=await loadCommon(),local=boundMarginalSubstitutionMass(target,[ref],[1],o,()=>{},new Set(['zero'])),global=boundMarginalSubstitutions(target,[ref],[1],o);
 expect(local.lower[0]).toBeCloseTo(0,6);expect(local.upper[0]).toBeCloseTo(0,6);expect(global[0].scope).toBe('all-length-groups');expect(global[0].fractionRange[0]).toBeCloseTo(.1,5);expect(global[0].fractionRange[1]).toBeCloseTo(.1,5);expect(global[0].zeroFeasible).toBe(false);
});
it('keeps within-reference normalization separate from molecular mass',()=>{
 const b=zeroLengthCompositionBounds([.15,.15],.1,.25,[.025,.025]);expect(b.substitutionRange[0]).toBeCloseTo(.1,5);expect(b.unmodifiedRange[1]).toBeCloseTo(.5,5);
});
