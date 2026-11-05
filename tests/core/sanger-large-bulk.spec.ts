import {it,expect,describe} from 'vitest';
import {readFileSync} from 'node:fs';
import {parseAbif} from '../../src/app/utils/abif';
import {analyzeSanger} from '../../src/app/workers/core/sanger-analysis';
import {shiftedControlQuality} from '../../src/app/workers/core/sanger-signal';
describe('full-junction large-deletion fits',{timeout:30000},()=>{
const base=new URL('../fixtures/sanger/large-bulk/',import.meta.url);
function input(label:string){
 const c=JSON.parse(readFileSync(new URL(label+'.json',base),'utf8'));
 const load=(f:string)=>{const b=readFileSync(new URL(f,base));return parseAbif(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),f);};
 return {c,data:{referenceName:label,reference:c.reference,cut:c.cut,control:load(c.files.control),sample:load(c.files.sample),settings:c.settings}};
}
it.each(['reverse-31','context-200'])('recovers several simultaneous small and large deletion classes: %s',label=>{
 const {c,data}=input(label),r=analyzeSanger(data),p=new Map(r.mixture!.components.map(x=>[x.size,x.fraction]));
 expect(r.mixture!.passesChecks).toBe(true);
 for(const t of c.truth.net)expect(Math.abs((p.get(t.size)||0)-t.fraction)).toBeLessThan(.05);
 expect(r.diagnostics!.modeledRange).toEqual({maxDeletion:225,maxInsertion:6});
 expect(r.mixture!.components).toHaveLength(232);
 if(c.reverseInputReference)expect(r.direction).toBe('reverse');
});
it('returns a QC warning for a strong out-of-range deletion',()=>{
 const {data}=input('outside-300'),r=analyzeSanger(data);
 expect(r.mixture!.passesChecks).toBe(false);expect(r.annotations![0].estimable).toBe(false);
});
it('flags a deletion at the configured deletion boundary',()=>{
 const {data}=input('boundary-225'),r=analyzeSanger(data);
 expect(r.mixture!.components.find(x=>x.size===-225)!.fraction).toBeCloseTo(.5,1);
 expect(r.mixture!.passesChecks).toBe(false);expect(r.warnings.join(' ')).toMatch(/boundary/);
});
it('uses known reference responses with a warning when downstream control donors are unmeasured',()=>{
 const {data}=input('context-200'),control=data.control;
 data.control={...control,sequence:control.sequence.slice(0,420),peaks:control.peaks.slice(0,420),quality:control.quality!.slice(0,420)};
 const r=analyzeSanger(data);expect(r.mixture!.components.every(c=>Number.isFinite(c.fraction))).toBe(true);expect(r.warnings.join(' ')).toMatch(/beyond the measured control/);
});
it('does not invent a bulk edit when a real waveform is compared with itself under the broad range',()=>{
 const {data}=input('context-200');data.sample=data.control;
 const r=analyzeSanger(data);expect(r.mixture!.indel).toBeLessThan(1e-6);expect(r.mixture!.passesChecks).toBe(true);
});
it('evaluates asymmetric control coverage without also demanding an unnecessarily long insertion flank',()=>{
 const q=new Uint8Array(250).fill(40);
 expect(shiftedControlQuality(q,20,6,15,200)).toBe(1);
 expect(shiftedControlQuality(q,20,200,15)).toBe(0);
 q.fill(2,170,225);expect(shiftedControlQuality(q,20,6,15,200)).toBe(0);
});

it('returns a fitted spectrum with a QC warning for indistinguishable repetitive length models',()=>{const {data}=input('repetitive'),r=analyzeSanger(data);expect(r.mixture!.passesChecks).toBe(false);expect(r.mixture!.components.every(c=>Number.isFinite(c.fraction))).toBe(true);expect(r.mixture!.components.reduce((s,c)=>s+c.fraction,0)).toBeCloseTo(1,6);expect(r.warnings.join(' ')).toMatch(/indistinguishable/);});

});
