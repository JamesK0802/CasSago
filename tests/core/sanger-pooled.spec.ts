import { describe,it,expect } from 'vitest';
import { analyzePooledSanger,checkPooledSeparation,gramCondition,PooledSeparationError } from '../../src/app/workers/core/sanger-pooled';
import { sangerPooledDemo } from '../../src/app/utils/sanger-analysis-demo';
import { parseAbif } from '../../src/app/utils/abif';
import { SANGER_ANALYSIS_DEFAULTS } from '../../src/app/models/sanger-analysis';
import { reverseComplement } from '../../src/app/utils/biology.utils';
const setup=()=>{const d=sangerPooledDemo();return {...d,ampliconConfirmed:true,control:parseAbif(d.control,'control'),sample:parseAbif(d.sample,'sample'),settings:{...SANGER_ANALYSIS_DEFAULTS,scope:'pooled' as const,maxIndel:5}};};
describe('pooled amplicon separation gate',()=>{
  it('checks every amplicon/indel model on the submitted primer axis',()=>{
    const d=setup(),p=checkPooledSeparation(d.amplicons,d.primer,d.settings);
    expect(p.passed).toBe(true);expect(p.modelCount).toBe(22);expect(p.rank).toBe(22);expect(p.condition!).toBeLessThan(100);expect(p.minimumCutOffset).toBe(60);expect(p.amplicons.map(a=>a.cutOffset)).toEqual([230,290]);
  });
  it('flags identical products but returns computed estimates with a reference ambiguity warning',()=>{
    const d=setup();d.amplicons[1]={...d.amplicons[0],name:'Identical'};
    const p=checkPooledSeparation(d.amplicons,d.primer,d.settings);expect(p.passed).toBe(false);expect(p.minimumDifferences).toBe(0);
    const stages:string[]=[],r=analyzePooledSanger(d,s=>stages.push(s));expect(r.mixture!.passesChecks).toBe(false);expect(r.pooled!.amplicons).toHaveLength(2);expect(r.pooled!.amplicons.every(a=>Number.isFinite(a.share)&&a.indel!==null)).toBe(true);expect(r.warnings.join(' ')).toMatch(/Reference QC warning/);expect(stages.some(s=>s.includes('Registering'))).toBe(true);
  });
  it('does not mistake a different FASTA origin for a different sequencing start',()=>{
    const d=setup(),first=d.amplicons[0];d.amplicons=[first,{name:'Padded origin',sequence:'A'.repeat(200)+first.sequence,cut:first.cut+200}];
    const p=checkPooledSeparation(d.amplicons,d.primer,d.settings);expect(p.passed).toBe(false);expect(p.minimumCutOffset).toBe(0);
  });
  it('uses downstream reference markers even outside the reporting window',()=>{
    const d=setup(),first=d.amplicons[0],bases=first.sequence.split('');
    for(let i=first.cut+105;i<first.cut+180;i+=3)bases[i]=bases[i]==='A'?'C':'A';
    d.amplicons=[first,{name:'SNP homeolog',sequence:bases.join(''),cut:first.cut}];
    const narrow=checkPooledSeparation(d.amplicons,d.primer,{...d.settings,maxIndel:1,window:60});expect(narrow.passed).toBe(true);
    const expanded=checkPooledSeparation(d.amplicons,d.primer,{...d.settings,maxIndel:1,window:400});
    expect(expanded.positions).toEqual(narrow.positions);expect(narrow.positions).toContain(first.cut-d.primer.length+150);expect(expanded.minimumCutOffset).toBe(0);expect(expanded.passed).toBe(true);
  });
  it('tests full model rank beyond pairwise sequence differences',()=>{
    expect(gramCondition([[1,0],[0,1]])).toEqual({rank:2,condition:1});
    expect(gramCondition([[1,0],[1,0]])).toEqual({rank:1,condition:null});
    const d=setup(),first=d.amplicons[0],extension=first.sequence.slice(d.primer.length);
    d.amplicons[1]={name:'Shifted two bases',sequence:d.primer+'AC'+extension.slice(0,-2),cut:first.cut+2};
    const p=checkPooledSeparation(d.amplicons,d.primer,{...d.settings,maxIndel:3});
    expect(p.passed).toBe(false);
  });
  it('handles a reverse-complement submitted amplicon with the same sequencing primer',()=>{
    const d=setup(),r=d.amplicons[1];d.amplicons[1]={...r,sequence:reverseComplement(r.sequence),cut:r.sequence.length-r.cut};
    const p=checkPooledSeparation(d.amplicons,d.primer,d.settings);expect(p.passed).toBe(true);expect(p.amplicons[1].direction).toBe('reverse');expect(p.amplicons[1].cutOffset).toBe(290);
  });
  it('rejects missing/duplicate primer sites, unknown bases and clipped target windows',()=>{
    const d=setup();expect(checkPooledSeparation(d.amplicons,'AAAAAAAAAAAA',d.settings).passed).toBe(false);
    const duplicate=d.amplicons.map((r,i)=>i?{...r,sequence:r.sequence+d.primer}:r);expect(checkPooledSeparation(duplicate,d.primer,d.settings).reason).toMatch(/one exact/);
    const unknown=d.amplicons.map((r,i)=>i?{...r,sequence:r.sequence+'N'}:r);expect(checkPooledSeparation(unknown,d.primer,d.settings).passed).toBe(false);
    expect(checkPooledSeparation(d.amplicons,d.primer,{...d.settings,readEnd:300}).passed).toBe(false);
    expect(()=>analyzePooledSanger({...d,ampliconConfirmed:false})).toThrow(/actual full amplicon/);
  });
});
describe('reference-guided pooled trace fitting (synthetic evidence only)',()=>{
  it('recovers known global shares and within-amplicon edit fractions',()=>{
    const r=analyzePooledSanger(setup()),a=r.pooled!.amplicons;
    expect(r.mixture!.r2).toBeGreaterThan(.99);expect(r.pooled!.sampleTrim).toBe(12);expect(r.pooled!.controlTrim).toBe(12);
    expect(a[0].share).toBeCloseTo(.55,2);expect(a[1].share).toBeCloseTo(.45,2);expect(a[0].indel!).toBeCloseTo(.4,2);expect(a[1].indel!).toBeCloseTo(.6,2);
    expect(r.mixture!.indel).toBeCloseTo(.49,2);expect(r.plots[0].cuts).toHaveLength(2);
  });
  it('does not discard mixed signals based on base-call Phred',()=>{
    const d=setup();d.sample.quality!.fill(0);d.control.quality!.fill(0);
    expect(analyzePooledSanger(d).mixture!.r2).toBeGreaterThan(.99);
  });
  it('requires enough actual downstream support at every target',()=>{
    const d=setup();d.sample={...d.sample,sequence:d.sample.sequence.slice(0,300),peaks:d.sample.peaks.slice(0,300),quality:d.sample.quality!.slice(0,300)};
    expect(()=>analyzePooledSanger(d)).toThrow(/usable trace must cover/);
  });
});

describe('calibrated trace result contract',()=>{
 it('returns annotation models for every reference and retains accurate classic peak-normalization estimates',()=>{const d=setup(),enhanced=analyzePooledSanger(d),classic=analyzePooledSanger({...d,settings:{...d.settings,signalModel:'classic'}});expect(enhanced.annotations!.map(a=>a.name)).toEqual(d.amplicons.map(a=>a.name));expect(enhanced.annotations![0].reference).toBe(d.amplicons[0].sequence);expect(enhanced.diagnostics?.sampleCalibration?.validationBases).toBeGreaterThan(20);expect(classic.diagnostics?.signalModel).toBe('classic');expect(Math.abs(classic.pooled!.amplicons[0].indel!-.4)).toBeLessThan(.01);expect(Math.abs(classic.pooled!.amplicons[1].indel!-.6)).toBeLessThan(.01);expect(classic.diagnostics?.discordantPositions).toBeUndefined();});
 it('does not accept featureless four-channel noise as a calibrated trace',()=>{const d=setup();for(const trace of [d.control,d.sample])for(const b of ['A','C','G','T'] as const)trace.channels[b].fill(1000);const result=analyzePooledSanger(d);expect(result.mixture!.passesChecks).toBe(false);expect(Number.isFinite(result.mixture!.indel)).toBe(true);});
});
