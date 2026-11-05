import { describe, expect, it } from 'vitest';
import { designSangerPrimers, nearestNeighbourTm, exactBindingSites, complementarity, validateDesign, nearBindingSiteCount } from '../../src/app/workers/core/sanger-primer-designer';
import { SANGER_DEFAULTS, SangerReference } from '../../src/app/models/sanger-design';
import { guideCut, parseSangerFasta, sangerDesignExample, resolveSangerReference } from '../../src/app/utils/sanger-design.utils';
import { reverseComplement } from '../../src/app/utils/biology.utils';

const settings = () => ({...SANGER_DEFAULTS,specificityStrategy:'mismatch' as const,primerMin:20,primerMax:20});
describe('identical-primer Sanger design',() => {
  it('designs finite single-reference PCR pairs without requiring homeolog offsets',()=>{
    const ref=sangerDesignExample()[0],out=designSangerPrimers([ref],{...settings(),workflow:'standard-pcr'});
    expect(out.candidates.length).toBeGreaterThan(0);
    for(const c of out.candidates){expect(c.partner).not.toBeNull();expect(c.minOffset).toBe(0);expect(Number.isFinite(c.traceComparedBases)).toBe(true);expect(c.indelOffsetOverlap).toBe(false);}
    expect(()=>validateDesign([ref],settings())).toThrow(/2–8/);
    expect(()=>validateDesign(sangerDesignExample(),{...settings(),workflow:'standard-pcr'})).toThrow(/one reference/);
  });
  it('finds one shared oligo with distinct cut distances and a common PCR partner',() => {
    const refs = sangerDesignExample();const out = designSangerPrimers(refs,settings());
    const c = out.candidates.find(c => c.sequence === 'CCGGCGGGTACGACCATAAC');
    expect(c).toBeDefined();expect(c!.direction).toBe('forward');
    expect(c!.sites.map(s => s.cutOffsets)).toEqual([[150],[200],[255]]);
    expect(c!.minOffset).toBe(50);expect(c!.partner).not.toBeNull();
    for (let i=0;i<refs.length;i++) {
      const hit = c!.sites[i], partner = c!.partner!.sites[i];
      expect(refs[i].sequence.slice(hit.start,hit.end)).toBe(c!.sequence);
      expect(reverseComplement(refs[i].sequence.slice(partner.start,partner.end))).toBe(c!.partner!.sequence);
      expect(c!.partner!.products[i]).toBe(partner.end-hit.start);
    }
    expect(c!.minTraceDifferences).toBeGreaterThan(0);
  });
  it('limits expected traces to the PCR product in both orientations',() => {
    const refs = sangerDesignExample();
    for (const direction of ['forward','reverse'] as const) {
      const input = direction === 'forward' ? refs : refs.map(r => ({...r,sequence:reverseComplement(r.sequence),cuts:r.cuts.map(c=>r.sequence.length-c)}));
      const out = designSangerPrimers(input,{...settings(),direction});
      expect(out.candidates.length).toBeGreaterThan(0);
      for (const candidate of out.candidates) candidate.sites.forEach((site,i) => {
        expect(site.extension.length).toBeLessThanOrEqual(candidate.partner!.products[i]-candidate.sequence.length);
        const lastWindowEnd=Math.max(...candidate.sites.flatMap(s=>s.cutOffsets))+Math.ceil(settings().window/2)+settings().indelRange;
        expect(candidate.partner!.products[i]-candidate.sequence.length).toBeGreaterThanOrEqual(lastWindowEnd);
      });
    }
  });
  it('accepts full gene references and does not require an offset on the opposite side',() => {
    const refs=sangerDesignExample().map(r=>({...r,sequence:'ACGT'.repeat(6000)+r.sequence,cuts:r.cuts.map(c=>c+24000)}));
    validateDesign(refs,settings());
    const out=designSangerPrimers(refs,settings());
    const c=out.candidates.find(c=>c.sequence==='CCGGCGGGTACGACCATAAC')!;
    expect(c).toBeDefined();expect(c.minOffset).toBe(50);
    expect(new Set(c.partner!.sites.map(s=>s.cutOffsets[0])).size).toBe(1);
    expect(c.sites.every(s=>s.start>20000)).toBe(true);
  });
  it('designs separate PCR pairs with one specific side while allowing a shared opposite primer',() => {
    const refs=sangerDesignExample();
    const out=designSangerPrimers(refs,{...settings(),workflow:'homeolog-pcr'});
    expect(out.coverage!.every(r=>r.candidates>0)).toBe(true);
    expect(out.candidates.some(c=>c.specificity!.primary&&!c.specificity!.partner)).toBe(true);
    for(const c of out.candidates) {
      expect(c.referenceIndex).toBeDefined();expect(c.sites).toHaveLength(1);
      expect(c.specificity!.primary||c.specificity!.partner).toBe(true);
      const own=c.referenceIndex!;
      for(const [seq,specific] of [[c.sequence,c.specificity!.primary],[c.partner!.sequence,c.specificity!.partner]] as const) {
        expect(nearBindingSiteCount(seq,[refs[own].sequence,reverseComplement(refs[own].sequence)],2)).toBe(1);
        if(specific) refs.forEach((r,i)=>{if(i!==own) expect(nearBindingSiteCount(seq,[r.sequence,reverseComplement(r.sequence)],2)).toBe(0);});
      }
    }
    expect(out.complete).toBe(false); // The ranked search is explicitly bounded.
  });
  it('honours both-specific and reverse sequencing constraints',() => {
    const out=designSangerPrimers(sangerDesignExample(),{...settings(),workflow:'homeolog-pcr',specificBoth:true,direction:'reverse'});
    expect(out.candidates.length).toBeGreaterThan(0);
    for(const c of out.candidates) {
      expect(c.direction).toBe('reverse');expect(c.specificity!.primary&&c.specificity!.partner).toBe(true);
      expect(c.partner!.sites[0].direction).toBe('forward');
      expect(c.partner!.products[0]).toBe(c.sites[0].end-c.partner!.sites[0].start);
    }
  });
  it('never labels identical references as homeolog-specific and rejects unknown bases for specificity',() => {
    const ref=sangerDesignExample()[0];
    const out=designSangerPrimers([ref,{...ref,name:'Identical'}],{...settings(),workflow:'homeolog-pcr',specificityMismatches:0});
    expect(out.candidates).toHaveLength(0);expect(out.coverage!.every(r=>!r.candidates)).toBe(true);
    expect(()=>designSangerPrimers([ref,{...ref,name:'Unknown',sequence:ref.sequence+'N'}],{...settings(),workflow:'homeolog-pcr'})).toThrow(/Resolve N/);
  });
  it('uses target-relative distances, not shifted reference coordinates',() => {
    const original = sangerDesignExample()[0];
    const refs = [original,{...original,name:'Shifted coordinates',sequence:'GATC'.repeat(10)+original.sequence,cuts:original.cuts.map(c=>c+40)}];
    const out = designSangerPrimers(refs,{...settings(),requirePartner:false});
    expect(out.candidates).toHaveLength(0);expect(out.stats.geometryPassed).toBeGreaterThan(0);expect(out.stats.offsetPassed).toBe(0);
  });
  it('recovers reverse primers with correct coordinates and read offsets',() => {
    const refs = sangerDesignExample().map(r => ({...r,sequence:reverseComplement(r.sequence),cuts:r.cuts.map(c=>r.sequence.length-c)}));
    const out = designSangerPrimers(refs,{...settings(),direction:'reverse'});
    const c = out.candidates.find(c => c.sequence==='CCGGCGGGTACGACCATAAC');
    expect(c).toBeDefined();expect(c!.sites.map(s=>s.cutOffsets[0])).toEqual([150,200,255]);
    c!.sites.forEach((s,i)=>expect(reverseComplement(refs[i].sequence.slice(s.start,s.end))).toBe(c!.sequence));
    expect(c!.partner!.sites.every(s=>s.direction==='forward')).toBe(true);
  });
  it('rejects primary oligos with another exact site, including reverse sites',() => {
    const refs = sangerDesignExample();const oligo='CCGGCGGGTACGACCATAAC';
    refs[1].sequence += reverseComplement(oligo);
    expect(exactBindingSites(oligo,refs[1])).toHaveLength(2);
    const out=designSangerPrimers(refs,{...settings(),requirePartner:false});
    expect(out.candidates.some(c=>c.sequence===oligo)).toBe(false);
  });
  it('requires every target window and rejects truncated required windows',() => {
    const refs=sangerDesignExample().map(r=>({...r,cuts:[...r.cuts,3]}));
    const out=designSangerPrimers(refs,{...settings(),requirePartner:false});
    expect(out.candidates).toHaveLength(0);
  });
  it('checks corresponding targets separately instead of sorting by strand',() => {
    const refs=sangerDesignExample().map(r=>({...r,cuts:[r.cuts[0],r.cuts[0]+50]}));
    const out=designSangerPrimers(refs,{...settings(),requirePartner:false});
    const c=out.candidates.find(c=>c.sequence==='CCGGCGGGTACGACCATAAC');
    expect(c!.sites.map(s=>s.cutOffsets)).toEqual([[150,200],[200,250],[255,305]]);expect(c!.minOffset).toBe(50);
  });
  it('can return a sequencing primer without promising a PCR pair',() => {
    const out=designSangerPrimers(sangerDesignExample(),{...settings(),requirePartner:false,productMin:19000,productMax:20000});
    expect(out.candidates.length).toBeGreaterThan(0);expect(out.candidates.every(c=>c.partner===null)).toBe(true);
    const paired=designSangerPrimers(sangerDesignExample(),{...settings(),productMin:19000,productMax:20000});
    expect(paired.candidates).toHaveLength(0);expect(paired.stats.offsetPassed).toBeGreaterThan(0);expect(paired.complete).toBe(true);
  });
  it('marks geometric indel-range collisions without claiming deconvolution',() => {
    const out=designSangerPrimers(sangerDesignExample(),{...settings(),requirePartner:false,indelRange:30});
    expect(out.candidates.find(c=>c.sequence==='CCGGCGGGTACGACCATAAC')!.indelOffsetOverlap).toBe(true);
  });
  it('requires exact shared sequence and distinguishes no-candidate stage',() => {
    const refs: SangerReference[]=[{name:'A',sequence:'A'.repeat(500),cuts:[200]},{name:'C',sequence:'C'.repeat(500),cuts:[200]}];
    const out=designSangerPrimers(refs,settings());expect(out.candidates).toHaveLength(0);expect(out.stats.commonUnique).toBe(0);expect(out.complete).toBe(true);
  });
  it('rejects malformed dimensions and unpaired target lists',() => {
    expect(()=>validateDesign(sangerDesignExample(),{...settings(),primerMin:30,primerMax:18})).toThrow();
    expect(()=>validateDesign(sangerDesignExample(),{...settings(),minOffset:NaN})).toThrow();
    expect(()=>validateDesign(sangerDesignExample().map((r,i)=>({...r,cuts:i===0?[1,2]:r.cuts})),settings())).toThrow(/corresponding/);
    expect(()=>validateDesign(sangerDesignExample(),{...settings(),window:90.5})).toThrow(/whole/);
  });
});
describe('Sanger design sequence helpers',()=>{
  it('reads multiline FASTA without silently deleting invalid bases',()=>{
    expect(parseSangerFasta('>A\nacgt\nACGT\n>B\nGGCC')).toEqual([{name:'A',sequence:'ACGTACGT'},{name:'B',sequence:'GGCC'}]);
    expect(()=>parseSangerFasta('>A\nAC!GT')).toThrow();expect(()=>parseSangerFasta('>A\nACGT\n>A\nGGCC')).toThrow(/unique/);
  });
  it('locates forward and reverse SpCas9 cut boundaries with PAM support',()=>{
    const guide='GACCTAGTGCATGACCTGAC', dna='AAAAA'+guide+'AGG'+'TTTT';
    expect(guideCut(dna,guide)).toBe(22);expect(guideCut(reverseComplement(dna),guide)).toBe(dna.length-22);
    expect(guideCut('AAAAA'+guide+'AAA',guide)).toBe(22);
    expect(()=>guideCut('AAAAA'+guide.slice(0,-1)+'AAAA',guide)).toThrow(/No exact/);
  });
  it('does not accept ambiguous guide sites or non-numeric coordinate tokens',()=>{
    const guide='GACCTAGTGCATGACCTGAC';expect(()=>guideCut((guide+'AGG').repeat(2),guide)).toThrow(/multiple/);
    expect(()=>resolveSangerReference('A','ACGT','2bp','cuts')).toThrow();
  });
  it('counts near-binding sites across all seed partitions, strands and reference positions',()=>{
    const oligo='CCGGCGGGTACGACCATAAC';
    for(const positions of [[],[0],[8],[19],[0,19],[5,12]]) {
      const edited=oligo.split('');for(const p of positions) edited[p]=edited[p]==='A'?'C':'A';
      const dna='AAAA'+edited.join('')+'TTTT';
      expect(nearBindingSiteCount(oligo,[dna,reverseComplement(dna)],positions.length)).toBe(1);
      if(positions.length) expect(nearBindingSiteCount(oligo,[dna,reverseComplement(dna)],positions.length-1)).toBe(0);
    }
    expect(nearBindingSiteCount(oligo,[oligo+'AAAA'+oligo],0)).toBe(2);
  });
  it('computes salt/concentration-sensitive NN melting temperatures',()=>{
    const seq='CCGGCGGGTACGACCATAAC', t=nearestNeighbourTm(seq,50,250);
    expect(t).toBeGreaterThan(45);expect(t).toBeLessThan(75);
    expect(nearestNeighbourTm(reverseComplement(seq),50,250)).toBeCloseTo(t,8);
    expect(nearestNeighbourTm(seq,100,250)).toBeGreaterThan(t);
    expect(nearestNeighbourTm(seq,50,500)).toBeGreaterThan(t);
    expect(nearestNeighbourTm('ACNT',50,250)).toBeNaN();
  });
  it('detects complement runs anchored at either oligo 3-prime end',()=>{
    expect(complementarity('ACGTAC','GTACGT').end).toBe(6);
    expect(complementarity('GGGGG','CCCCC').end).toBe(5);
    expect(complementarity('AAAAA','AAAAA').any).toBe(0);
  });
});


// Bounded integration computation; group/fraction assertions remain unchanged.
it('keeps enough usable upstream calibration sequence for single-reference Sanger analysis',()=>{const out=designSangerPrimers([sangerDesignExample()[0]],{...SANGER_DEFAULTS,workflow:'standard-pcr',direction:'forward'});expect(out.candidates.length).toBeGreaterThan(0);expect(Math.min(...out.candidates.flatMap(c=>c.sites.flatMap(s=>s.cutOffsets)))).toBeGreaterThanOrEqual(SANGER_DEFAULTS.readStart+72);},120000);
