import { describe,it,expect } from 'vitest';
import { readFileSync } from 'node:fs';
import { parseAbif } from '../../src/app/utils/abif';
import { sangerAnalysisDemo,writeAbif } from '../../src/app/utils/sanger-analysis-demo';
import { analyzeSanger,alignTraceSequence,fitTraceMixture } from '../../src/app/workers/core/sanger-analysis';
import { SANGER_ANALYSIS_DEFAULTS,TRACE_BASES,SangerTrace } from '../../src/app/models/sanger-analysis';
import { reverseComplement } from '../../src/app/utils/biology.utils';
const setup=()=>{const d=sangerAnalysisDemo();return {referenceName:'Synthetic',reference:d.reference,cut:d.cut,sample:parseAbif(d.sample,'sample.ab1'),control:parseAbif(d.control,'control.ab1'),settings:{...SANGER_ANALYSIS_DEFAULTS},synthetic:true};};
function traceForSequence(sequence:string):SangerTrace {
  const peaks=Array.from({length:sequence.length},(_,i)=>12+10*i),channels=TRACE_BASES.map(()=>new Uint16Array(sequence.length*10+30).fill(1));
  Array.from(sequence).forEach((b,i)=>{const c=TRACE_BASES.indexOf(b as any);if(c>=0)channels[c][peaks[i]]=1000;});
  return parseAbif(writeAbif(sequence,peaks,channels,'Test','GATC'));
}
describe('ABIF analyzed chromatogram reader',()=>{
  it.each([['3100.ab1',795,10303],['3730.ab1',1165,16302]])('reads public instrument fixture %s',(name,bases,points)=>{
    const bytes=readFileSync(new URL('../fixtures/sanger/'+name,import.meta.url));
    const t=parseAbif(bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength),String(name));
    expect(t.sequence.length).toBe(bases);expect(t.peaks.length).toBe(bases);expect(t.channels.A.length).toBe(points);expect(t.quality!.length).toBe(bases);
    expect(new Set(Object.keys(t.channels))).toEqual(new Set(TRACE_BASES));
  });
  it('uses inline FWO dye order rather than assuming a channel order',()=>{
    const seq='ACGT'.repeat(30),t=traceForSequence(seq);
    for(let i=0;i<seq.length;i++){const peak=t.peaks[i];expect(t.channels[seq[i] as keyof typeof t.channels][peak]).toBe(1000);}
  });
  it('rejects corrupt directories, missing processed traces and invalid peak positions',()=>{
    const demo=sangerAnalysisDemo();
    expect(()=>parseAbif(new ArrayBuffer(8))).toThrow(/Truncated/);
    const bad=demo.control.slice(0);new DataView(bad).setUint32(26,0xfffffff0);expect(()=>parseAbif(bad)).toThrow(/offset/);
    const missing=demo.control.slice(0),v=new DataView(missing),bytes=new Uint8Array(missing);const root=v.getUint32(26);
    bytes.set(new TextEncoder().encode('RAWW'),root+5*28);expect(()=>parseAbif(missing)).toThrow(/DATA9/);
    const wrong=sangerAnalysisDemo().control;const w=new DataView(wrong);const peakEntry=w.getUint32(26)+2*28;const peakOffset=w.getUint32(peakEntry+20);w.setUint16(peakOffset+2,w.getUint16(peakOffset));expect(()=>parseAbif(wrong)).toThrow(/positions/);
  });
});
describe('single-locus Sanger analysis',()=>{
  it('recovers the known synthetic WT/deletion/insertion signal mixture',()=>{
    const r=analyzeSanger(setup());expect(r.mixture!.wt).toBeCloseTo(.6,2);expect(r.mixture!.components.find(c=>c.size===-1)!.fraction).toBeCloseTo(.25,2);expect(r.mixture!.components.find(c=>c.size===1)!.fraction).toBeCloseTo(.15,2);expect(r.mixture!.r2).toBeGreaterThan(.99);expect(r.mixture!.passesChecks).toBe(true);expect(r.synthetic).toBe(true);expect(r.plots).toHaveLength(2);
  });
  it('finds the reverse read orientation without changing the sign of indels',()=>{
    const input=setup();input.reference=reverseComplement(input.reference);input.cut=input.reference.length-input.cut;
    const r=analyzeSanger(input);expect(r.direction).toBe('reverse');expect(r.mixture!.components.find(c=>c.size===-1)!.fraction).toBeCloseTo(.25,2);
  });
  it.each(['3100.ab1','3730.ab1'])('reports WT when a real %s trace is compared to itself',name=>{
    const b=readFileSync(new URL('../fixtures/sanger/'+name,import.meta.url)),t=parseAbif(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),name);
    const r=analyzeSanger({referenceName:'Self-check',reference:t.sequence.replace(/[^ACGTN]/g,'N'),cut:250,sample:t,control:t,settings:{...SANGER_ANALYSIS_DEFAULTS,maxIndel:5,controlPhred:10}});
    expect(r.mixture!.wt).toBeGreaterThan(.9999);expect(r.mixture!.r2).toBeCloseTo(1,5);
  });
  it('suggests a supported wider single-locus range using held-out positions without changing the requested fit',()=>{
    const d=setup(),reference=d.control!.sequence,cut=230,control=traceForSequence(reference),sample=traceForSequence(reference.slice(0,cut)+reference.slice(cut+14)),input={...d,reference,cut,control,sample,settings:{...d.settings,maxIndel:10}};
    const r=analyzeSanger(input);expect(r.diagnostics!.requestedMaxIndel).toBe(10);expect(r.diagnostics!.rangeSuggestion!.maxIndel).toBeGreaterThanOrEqual(14);expect(r.mixture!.passesChecks).toBe(false);
    const wider=analyzeSanger({...input,settings:{...input.settings,maxIndel:15}});expect(wider.mixture!.components.find(c=>c.size===-14)!.fraction).toBeGreaterThan(.99);
    expect(analyzeSanger({...input,settings:{...input.settings,checkWiderRange:false}}).diagnostics!.rangeSuggestion).toBeUndefined();
  });
  it('requires a control for mixture analysis but not consensus alignment',()=>{
    const input=setup();expect(()=>analyzeSanger({...input,control:undefined})).toThrow(/control/);
    const r=analyzeSanger({...input,control:undefined,settings:{...input.settings,method:'consensus'}});expect(r.mixture).toBeUndefined();expect(r.warnings.join(' ')).toMatch(/no editing frequency/);
  });
  it('calls a pure deletion and maps reverse-strand alleles to input reference coordinates',()=>{
    const input=setup(),wt=input.control!.sequence,read=wt.slice(0,230)+wt.slice(233),sample=traceForSequence(read);
    const f=analyzeSanger({...input,sample,settings:{...input.settings,method:'consensus'}});
    const deletion=f.variants.find(v=>v.netSize===-3)!;expect(deletion).toBeDefined();expect(deletion.type).toBe('In-frame indel');expect(input.reference.slice(deletion.position-1,deletion.position+2)).toBe(deletion.reference);
    const r=analyzeSanger({...input,sample,reference:reverseComplement(input.reference),cut:input.reference.length-input.cut,settings:{...input.settings,method:'consensus'}});
    const reverseDeletion=r.variants.find(v=>v.netSize===-3)!;expect(reverseDeletion.reference).toBe(reverseComplement(deletion.reference));expect(reverseComplement(input.reference).slice(reverseDeletion.position-1,reverseDeletion.position+2)).toBe(reverseDeletion.reference);
  });
  it('labels uncertain IUPAC base calls as ambiguous instead of substitutions',()=>{
    const input=setup(),wt=input.control!.sequence,sample=traceForSequence(wt.slice(0,240)+'N'+wt.slice(241));
    const r=analyzeSanger({...input,sample,settings:{...input.settings,method:'consensus'}});expect(r.variants.some(v=>v.type==='Ambiguous call')).toBe(true);
  });
  it('returns estimates for low control quality and rejects insufficient registration or wrong loci',()=>{
    const input=setup(),poor={...input.control!,quality:new Uint8Array(input.control!.sequence.length).fill(2)};
    const low=analyzeSanger({...input,control:poor});expect(Number.isFinite(low.mixture!.indel)).toBe(true);expect(low.mixture!.passesChecks).toBe(false);expect(low.warnings.join(' ')).toMatch(/quality|Phred/);
    expect(()=>analyzeSanger({...input,cut:80})).toThrow(/upstream/);
    expect(()=>analyzeSanger({...input,reference:'C'.repeat(input.reference.length)})).toThrow(/align/);
  });
  it('warns rather than disguising a poor downstream fit',()=>{
    const input=setup();const t=input.sample;
    for(let i=250;i<t.peaks.length;i++){const p=t.peaks[i];for(let k=-3;k<=3;k++)for(const b of TRACE_BASES)t.channels[b][p+k]=100;}
    const r=analyzeSanger(input);expect(r.mixture!.passesChecks).toBe(false);expect(r.warnings.join(' ')).toMatch(/R²|split|unexplained|range/);
  });
  it('does not fabricate statistical significance or read counts',()=>{
    const r=analyzeSanger(setup());expect('pValue' in r.mixture!).toBe(false);expect('readCount' in r).toBe(false);
    const fitted=fitTraceMixture([[1,0,1,0],[0,1,0,1]],[.7,.3,.7,.3]);expect(fitted.weights[0]).toBeCloseTo(.7,6);
    expect(alignTraceSequence('AAAAACCGTGGGG','CCGT').identity).toBe(1);
    expect(()=>fitTraceMixture([[1,0,1,0],[1,0,1,0]],[1,0,1,0])).toThrow(/indistinguishable/);
  });
});
