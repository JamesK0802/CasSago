import {describe,it,expect} from 'vitest';
import {guideCut,orientSangerReferences,sangerDesignExample} from '../../src/app/utils/sanger-design.utils';
import {reverseComplement} from '../../src/app/utils/biology.utils';
import {SANGER_DEFAULTS,SangerReference} from '../../src/app/models/sanger-design';
import {designSangerPrimers} from '../../src/app/workers/core/sanger-primer-designer';
import {selectedAmplicons,makeAnalysisTemplate} from '../../src/app/utils/sanger-design-export';

const flip=(r:SangerReference):SangerReference=>({...r,sequence:reverseComplement(r.sequence),cuts:r.cuts.map(c=>r.sequence.length-c)});
function fixture(){
 const refs=sangerDesignExample(2),guide=refs[0].sequence.slice(refs[0].cuts[0]-17,refs[0].cuts[0]+3);
 return refs.map(r=>({...r,guides:[guide],targetNames:['First','Second']}));
}
describe('mixed input reference strands',()=>{
 it('orients shared guide sites without mutating input or changing target identities',()=>{
  const refs=fixture(),mixed=[refs[0],flip(refs[1]),refs[2]],before=JSON.stringify(mixed);
  const out=orientSangerReferences(mixed);
  expect(out.map(r=>r.sequence)).toEqual(refs.map(r=>r.sequence));
  expect(out.map(r=>r.cuts)).toEqual(refs.map(r=>r.cuts));
  expect(out[1].inputReversed).toBe(true);expect(out[1].targetNames).toEqual(['First','Second']);
  expect(guideCut(out[1].sequence,out[1].guides![0])).toBe(out[1].cuts[0]);
  expect(orientSangerReferences(out)).toEqual(out);expect(JSON.stringify(mixed)).toBe(before);
 });
 it('anchors either input direction and leaves cuts-only inputs untouched',()=>{
  const refs=fixture(),out=orientSangerReferences([flip(refs[0]),refs[1],refs[2]]);
  expect(out.map(r=>r.sequence)).toEqual(refs.map(r=>reverseComplement(r.sequence)));
  const cutsOnly=refs.map(r=>({...r,guides:undefined}));
  expect(orientSangerReferences(cutsOnly)).toEqual(cutsOnly);
  expect(orientSangerReferences([{...refs[0],guides:['A'.repeat(20)]},flip(refs[1])])[1]).toEqual(flip(refs[1]));
 });
 it('retains physical primer pairs, safety checks and export coordinates',()=>{
  const refs=fixture(),mixed=[refs[0],flip(refs[1]),refs[2]];
  const out=designSangerPrimers(mixed,{...SANGER_DEFAULTS,primerMin:20,primerMax:20,readEnd:1000});
  const c=out.candidates.find(c=>c.sequence==='CCGGCGGGTACGACCATAAC')!;
  expect(c).toBeDefined();expect(c.minOffset).toBe(50);
  expect(c.sites.map(s=>s.cutOffsets)).toEqual([[150,400],[200,450],[255,505]]);
  expect(c.primerSafety!.primary.minimumMismatches).toBeGreaterThanOrEqual(5);
  const amp=selectedAmplicons(out,c)[1];
  expect(amp.inputStrand).toBe('reverse');
  expect(reverseComplement(mixed[1].sequence.slice(amp.inputStart,amp.inputEnd))).toBe(amp.sequence);
  expect(amp.cuts).toEqual(refs[1].cuts.map(cut=>cut-amp.start));
  expect(makeAnalysisTemplate(out,c).amplicons[1]).toEqual(amp);
 });
});
