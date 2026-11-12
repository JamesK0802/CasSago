import {expect,it} from 'vitest';
import {readFileSync} from 'node:fs';
import {SANGER_DEFAULTS,oppositePrimer} from '../../src/app/models/sanger-design';
import {designSangerPrimers} from '../../src/app/workers/core/sanger-primer-designer';
import {sharedMismatchStage} from '../../src/app/workers/core/sanger-shared-primer-designer';
import {resolveSangerReference} from '../../src/app/utils/sanger-design.utils';
import {makeAnalysisTemplate,selectedAmplicons} from '../../src/app/utils/sanger-design-export';
import {synthReferences,generateSangerSynth} from '../../src/app/utils/sanger-synthetic-benchmark';
import {reverseComplement} from '../../src/app/utils/biology.utils';
it('enforces the exact terminal ten and the inclusive three-base single-mismatch position spread',()=>{
 expect(sharedMismatchStage([[5],[8]],20)).toBe('single-aligned');
 expect(sharedMismatchStage([[5],[9]],20)).toBe('single-unrestricted');
 expect(sharedMismatchStage([[10],[10]],20)).toBe('single-aligned');
 expect(sharedMismatchStage([[11],[11]],20)).toBeNull();
 expect(sharedMismatchStage([[1,10],[8,9]],20)).toBe('double');
 expect(sharedMismatchStage([[1,10],[8]],20)).toBeNull();
 expect(sharedMismatchStage([[1,10],[]],20)).toBeNull();
 expect(sharedMismatchStage([[1,10],[8,9],[2]],20)).toBeNull();
 expect(sharedMismatchStage([[1,2,3],[1,2]],20)).toBeNull();
 expect(sharedMismatchStage([[0],[1]],20)).toBeNull();
});
it('finds symmetric CPC/TRY single mismatches, different exact partners and usable primer-incorporated products',()=>{
 const data=JSON.parse(readFileSync(new URL('../fixtures/primer-cpc-try.json',import.meta.url),'utf8'));
 const refs=Object.entries(data).map(([name,sequence])=>resolveSangerReference(name,String(sequence),'AATATCTCTCTATCTCCTC','guides'));
 const out=designSangerPrimers(refs,{...SANGER_DEFAULTS,primerMin:20,primerMax:20,readEnd:1000,sharedPrimerPolicy:'progressive'});
 expect(out.candidates.length).toBeGreaterThan(0);
 for(const c of out.candidates){expect(c.sharedDesign?.stage).toBe('single-aligned');expect(c.sharedDesign?.terminalExact).toBe(10);expect(c.sharedDesign?.positionSpread).toBeLessThanOrEqual(3);expect(c.sites.map(h=>h.mismatchPositions)).toEqual([[6],[6]]);expect(c.primerSafety!.primary.minimumMismatches).toBeGreaterThanOrEqual(5);expect(c.partner?.perReference).toHaveLength(2);expect(c.partnerSafety!.every(p=>p.minimumMismatches>=5)).toBe(true);}
 const c=out.candidates[0],amplicons=selectedAmplicons(out,c);
 expect(c.sequence).toBe('GTGAGGAGTATCGAATGGGA');
 amplicons.forEach((a,i)=>{const primary=c.sites[i],other=c.partner!.sites[i],q=oppositePrimer(c,i)!;expect(refs[i].sequence.slice(other.start,other.end)).toBe(reverseComplement(q.sequence));expect(a.sequence.startsWith(c.sequence)).toBe(true);expect(a.sequence.endsWith(reverseComplement(q.sequence))).toBe(true);expect(a.sequence.length).toBe(c.partner!.products[i]);expect(a.sequence.slice(c.sequence.length,-q.sequence.length)).toBe(refs[i].sequence.slice(primary.end,other.start));});
 const template=makeAnalysisTemplate(out,c);expect(template.oppositePrimer).toBeUndefined();expect(template.amplicons).toHaveLength(2);
 const synth=synthReferences(out,c);expect(synth.map(r=>r.partnerLength)).toEqual(c.sites.map((_,i)=>oppositePrimer(c,i)!.sequence.length));expect(()=>generateSangerSynth(out,c,synth,452)).not.toThrow();
},120000);

it('uses the wider-position tier only after no aligned single-mismatch pair passes',()=>{
 const data=JSON.parse(readFileSync(new URL('../fixtures/primer-cpc-try.json',import.meta.url),'utf8'));
 const refs=Object.entries(data).map(([name,sequence])=>resolveSangerReference(name,String(sequence),'AATATCTCTCTATCTCCTC','guides'));
 const out=designSangerPrimers(refs,{...SANGER_DEFAULTS,primerMin:21,primerMax:21,readEnd:1000,tmMin:55,tmMax:68,sharedPrimerPolicy:'progressive'});
 expect(out.candidates.length).toBeGreaterThan(0);expect(out.candidates.every(c=>c.sharedDesign?.stage==='single-unrestricted')).toBe(true);
 expect(out.candidates.every(c=>c.sharedDesign!.positionSpread>3&&c.sites.every(h=>h.mismatchPositions!.length===1&&h.mismatchPositions![0]<=c.sequence.length-10))).toBe(true);
},120000);

it('does not admit a four-mismatch unintended copy merely because the fast near-site scan stops at three',()=>{
 const data=JSON.parse(readFileSync(new URL('../fixtures/primer-cpc-try.json',import.meta.url),'utf8'));
 const refs=Object.entries(data).map(([name,sequence])=>resolveSangerReference(name,String(sequence),'AATATCTCTCTATCTCCTC','guides'));
 const oligo='GTGAGGAGTATCGAATGGGA',near=[...oligo];for(const i of [0,4,9,15])near[i]='ACGT'[('ACGT'.indexOf(near[i])+1)%4];
 refs[0].sequence+='ACACAC'+near.join('')+'GTGTGT';
 const out=designSangerPrimers(refs,{...SANGER_DEFAULTS,primerMin:20,primerMax:20,readEnd:1000,sharedPrimerPolicy:'progressive'});
 expect(out.candidates.some(c=>c.sequence===oligo)).toBe(false);
 expect(out.candidates.every(c=>c.primerSafety!.primary.minimumMismatches>=5&&c.sites.every(h=>h.mismatchPositions!.length<=2))).toBe(true);
},120000);

function constructedSharedReferences(commonOpposite:boolean,differences=0){
 let state=14923;const random=(n:number)=>Array.from({length:n},()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return 'ACGT'[state>>>30];}).join('');
 const primary='CCGGCGGGTACGACCATAAC',right='CAGCGACTCCTTGAGTGACG';
 return [0,50].map((gap,i)=>{let p=primary;if(i&&differences){const b=[...p];for(let j=0;j<differences;j++)b[j]='ACGT'[('ACGT'.indexOf(b[j])+1)%4];p=b.join('');}const prefix=random(60),before=random(150+gap),after=random(240),partner=commonOpposite?reverseComplement(right):random(20);return {name:'R'+i,sequence:prefix+p+before+after+partner+random(70),cuts:[prefix.length+p.length+before.length]};});
}
it('never relaxes all-pairs cut spacing when the primer mismatch tier is relaxed',()=>{
 const refs=constructedSharedReferences(true).map(r=>({...r,cuts:[r.cuts[0],r.cuts[0]+50]}));
 const out=designSangerPrimers(refs,{...SANGER_DEFAULTS,primerMin:20,primerMax:20,sharedPrimerPolicy:'progressive',sharedSearchStage:3});
 expect(out.candidates).toHaveLength(0);
},120000);
it('prefers both shared exact primers and supports reverse-oriented product handoff',()=>{
 const refs=constructedSharedReferences(true),options={...SANGER_DEFAULTS,primerMin:20,primerMax:20,sharedPrimerPolicy:'progressive' as const};
 const out=designSangerPrimers(refs,options);expect(out.candidates[0]?.sharedDesign?.stage).toBe('both-exact');expect(out.candidates[0].partner?.perReference).toBeUndefined();
 const reversed=refs.map(r=>({...r,sequence:reverseComplement(r.sequence),cuts:r.cuts.map(c=>r.sequence.length-c)}));
 const reverse=designSangerPrimers(reversed,{...options,direction:'reverse'});expect(reverse.candidates[0]?.direction).toBe('reverse');expect(()=>makeAnalysisTemplate(reverse,reverse.candidates[0])).not.toThrow();
},120000);
it('relaxes only the opposite PCR primer when an exact shared sequencing primer exists',()=>{
 const out=designSangerPrimers(constructedSharedReferences(false),{...SANGER_DEFAULTS,primerMin:20,primerMax:20,sharedPrimerPolicy:'progressive'});
 expect(out.candidates[0]?.sharedDesign?.stage).toBe('one-exact');expect(out.candidates[0].sites.every(h=>h.mismatchPositions!.length===0)).toBe(true);expect(out.candidates[0].partner?.perReference).toHaveLength(2);
},120000);
it('uses the final tier for two mismatches without a position-spread constraint and refuses a five-difference target pair',()=>{
 const settings={...SANGER_DEFAULTS,primerMin:20,primerMax:20,sharedPrimerPolicy:'progressive' as const};
 const out=designSangerPrimers(constructedSharedReferences(false,4),settings);expect(out.candidates.length).toBeGreaterThan(0);expect(out.candidates.every(c=>c.sharedDesign?.stage==='double'&&c.sites.every(h=>h.mismatchPositions!.length===2))).toBe(true);
 const none=designSangerPrimers(constructedSharedReferences(false,5),settings);expect(none.candidates).toHaveLength(0);
},120000);
