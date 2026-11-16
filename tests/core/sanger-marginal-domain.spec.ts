import {it,expect} from 'vitest';
import {marginalReferenceSubstitutionAllowed,marginalTargetReadSupport,marginalOriginReadBaseFree} from '../../src/app/workers/core/sanger-marginal-domain';
it('uses the same genomic base support after reverse-complementing a cut and reference',()=>{
 for(const length of [700,901])for(const cut of [170,450])for(let base=cut-48;base<=cut+48;base++)expect(marginalReferenceSubstitutionAllowed(base,cut)).toBe(marginalReferenceSubstitutionAllowed(length-base-1,length-cut));
 expect(marginalReferenceSubstitutionAllowed(undefined,170)).toBe(false);expect(marginalReferenceSubstitutionAllowed(125,170)).toBe(true);expect(marginalReferenceSubstitutionAllowed(214,170)).toBe(true);expect(marginalReferenceSubstitutionAllowed(215,170)).toBe(false);
});
it('covers every surviving genomic Sub base under all allowed local indels and incoming phases',()=>{
 const cut=170;
 for(const incoming of [-20,0,20])for(let size=-10;size<=10;size++){
  const support=marginalTargetReadSupport(cut,incoming,incoming,10,10),lo=size<0?cut+size-3:cut-3;
  for(let start=lo;start<=cut+3;start++)for(let donor=cut-45;donor<cut+45;donor++){
   const end=start+Math.max(0,-size);if(size<0&&donor>=start&&donor<end)continue;
   const q=donor+incoming+(donor>=end?size:0);expect(q).toBeGreaterThanOrEqual(support.first);expect(q).toBeLessThanOrEqual(support.last);
  }
 }
});
it('does not recalibrate the existing nominal noise/prior support when adding phase-registration coverage',()=>{
 const support=marginalTargetReadSupport(150,0,0,10,10);expect(support.first).toBe(105);expect(support.nominalLast).toBe(195);expect(support.last).toBe(204);expect(support.nominalLast-support.first+1).toBe(91);expect(support.last).toBeGreaterThanOrEqual(200);
});

it('uses the identical donor domain for selected-group alternatives and representatives across local/incoming phases',()=>{
 const cut=150;
 for(const incoming of [-20,0,20])for(let size=-10;size<=10;size++){
  const origin={start:147,incoming,outgoing:incoming+size,size};
  for(let q=90+incoming;q<=220+incoming;q++){
   const insertion=q>=origin.start+incoming&&q<origin.start+incoming+Math.max(0,size),donor=q<origin.start+incoming?q-incoming:q-origin.outgoing;
   expect(marginalOriginReadBaseFree(q,origin,cut)).toBe(insertion||marginalReferenceSubstitutionAllowed(donor,cut));
  }
 }
 expect(marginalOriginReadBaseFree(200,{start:150,incoming:0,outgoing:10,size:10},150)).toBe(true);
 expect(marginalOriginReadBaseFree(205,{start:150,incoming:0,outgoing:10,size:10},150)).toBe(false);
});
