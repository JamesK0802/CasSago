import {describe,it,expect} from 'vitest';
import {analyzeJointSanger,checkJointSeparation,applyJointEvents} from '../../src/app/workers/core/sanger-joint';
import {SANGER_ANALYSIS_DEFAULTS} from '../../src/app/models/sanger-analysis';
import {parseAbif} from '../../src/app/utils/abif';
import {reverseComplement} from '../../src/app/utils/biology.utils';
import {simulate,random} from '../fixtures/sanger-joint/forward.mjs';
const primer='ACGTTAGCTACGATCGATGC';
const dna=(seed:number)=>{const r=random(seed);return Array.from({length:950},()=> 'ACGT'[Math.floor(r()*4)]).join('');};
function fixture(two=true,large=false,seed=91){
 const seqs=[dna(seed),dna(seed+7)],cuts=two?[[260,370],[305,415]]:[[300],[345]];
 const amps=seqs.map((s,i)=>({name:'Ref'+i,sequence:primer+s,cut:cuts[i][0]+primer.length,cuts:cuts[i].map(c=>c+primer.length)}));
 // Independent molecule edits, no inference proposal generator used.
 const del=(s:string,at:number,n:number)=>s.slice(0,at)+s.slice(at+n);
 const ins=(s:string,at:number,n:number)=>s.slice(0,at)+'T'.repeat(n)+s.slice(at);
 const edited=two?[ins(del(seqs[0],260,1),369,1),del(seqs[1],305,2)]:[del(seqs[0],300,61),seqs[1]];
 const control=simulate(seqs.map((sequence,i)=>({sequence,weight:i?.4:.6})),{seed:1001,trim:10,length:850,preset:'moderate'});
 const sample=simulate([{sequence:seqs[0],weight:.3},{sequence:edited[0],weight:.3},{sequence:seqs[1],weight:.2},{sequence:edited[1],weight:.2}],{seed:1002,trim:8,length:850,preset:'moderate'});
 return {amplicons:amps,primer,ampliconConfirmed:true,settings:{...SANGER_ANALYSIS_DEFAULTS,scope:'pooled' as const,maxIndel:2,window:220,readEnd:850,largeDeletions:large,maxDeletion:80},control:parseAbif(control.buffer),sample:parseAbif(sample.buffer),synthetic:true};
}
describe('joint targets and bounded deletion proposals',()=>{
 it('applies cumulative edits rather than resetting at target 2',()=>{expect(applyJointEvents('AACCGGTTAACC',[{target:0,start:2,end:4,insert:0},{target:1,start:8,end:8,insert:1}])).toBe('AAGGTTNAACC');});
 it('allows overlapping windows and reports the largest separate window',()=>{const d=fixture(),p=checkJointSeparation(d.amplicons,primer,d.settings);expect(p.passed,p.reason).toBe(true);expect(p.amplicons[0].mergedWindows).toBe(true);expect(p.amplicons[0].maximumIndependentWindow).toBe(110);expect(p.joint!.candidateCount).toBeGreaterThan(50);});
 it('blocks overly close targets, >2 targets and invalid primer registration',()=>{const d=fixture();d.amplicons[0].cuts=[280,285];expect(checkJointSeparation(d.amplicons,primer,d.settings).reason).toMatch(/too close/);d.amplicons[0].cuts=[280,380,480];expect(checkJointSeparation(d.amplicons,primer,d.settings).passed).toBe(false);expect(checkJointSeparation(fixture().amplicons,'TTTTTTTTTTTTTTTTTTTT',d.settings).passed).toBe(false);});
 it('requires observed WT reference discrimination even with a large dictionary',()=>{const d=fixture();d.amplicons[1]={...d.amplicons[0],name:'duplicate'};expect(checkJointSeparation(d.amplicons,primer,d.settings).reason).toMatch(/cannot be distinguished/);});
 it('recovers a two-target zero-net edit without labeling it WT',()=>{const d=fixture(),r=analyzeJointSanger(d);console.log('JOINT2',JSON.stringify({qc:r.mixture!.passesChecks,r2:r.joint!.heldOutR2,groups:r.joint!.groups,targets:r.joint!.targets}));expect(r.joint!.heldOutR2).toBeGreaterThan(.9);expect(r.joint!.targets.find(t=>t.reference==='Ref0'&&t.target===1)!.edited).toBeCloseTo(.5,1);expect(r.pooled!.amplicons[0].wt!).toBeLessThan(.6);expect(r.joint!.groups.some(g=>g.reference==='Ref0'&&g.netSize===0&&g.events.length===2)).toBe(true);},120000);
 it('finds a 61 bp deletion outside the small-indel range in an independent sequence',()=>{const d=fixture(false,true,304),r=analyzeJointSanger(d);console.log('LARGE',JSON.stringify({qc:r.mixture!.passesChecks,r2:r.joint!.heldOutR2,groups:r.joint!.groups}));expect(r.joint!.heldOutR2).toBeGreaterThan(.9);expect(r.joint!.groups.some(g=>g.netSize===-61&&g.reference==='Ref0'&&Math.abs(g.fraction-.3)<.05)).toBe(true);},120000);
 it('supports reverse submitted products without changing target labels',()=>{const d=fixture(),amps=d.amplicons.map(a=>({...a,sequence:reverseComplement(a.sequence),cut:a.sequence.length-a.cut,cuts:a.cuts!.map(c=>a.sequence.length-c).sort((a,b)=>a-b)}));const p=checkJointSeparation(amps,primer,d.settings);expect(p.passed,p.reason).toBe(true);expect(p.amplicons.every(a=>a.direction==='reverse')).toBe(true);});
});

describe('joint unsupported-range safeguards',()=>{it('withholds reference editing fractions when deletion is outside the searched range',()=>{const d=fixture(false,true,304);d.settings.maxDeletion=40;const r=analyzeJointSanger(d);expect(r.mixture!.passesChecks).toBe(false);expect(r.joint!.heldOutR2).toBeLessThan(.9);expect(r.joint!.targets.every(t=>t.edited===null)).toBe(true);},120000);it('retains known reference and size when only excision boundaries are unresolved',()=>{const d=fixture();const a=d.amplicons[0].sequence.slice(primer.length),b=d.amplicons[1].sequence.slice(primer.length);d.sample=parseAbif(simulate([{sequence:a,weight:.3},{sequence:a.slice(0,260)+a.slice(370),weight:.3},{sequence:b,weight:.4}],{seed:1010,trim:8,length:850,preset:'moderate'}).buffer);const r=analyzeJointSanger(d),g=r.joint!.groups.find(g=>g.netSize===-110)!;expect(g.reference).toBe('Ref0');expect(g.alternatives.length).toBeGreaterThan(1);expect(g.events).toEqual([]);expect(r.joint!.targets.find(t=>t.reference==='Ref0')!.edited).toBeCloseTo(.5,1);},120000);});

it('refines an asymmetric large-deletion junction instead of forcing a centered breakpoint',()=>{const d=fixture(false,true,503);d.settings.maxDeletion=120;const a=d.amplicons[0].sequence.slice(primer.length),b=d.amplicons[1].sequence.slice(primer.length);d.sample=parseAbif(simulate([{sequence:a,weight:.3},{sequence:a.slice(0,283)+a.slice(368),weight:.3},{sequence:b,weight:.4}],{seed:1002,trim:8,length:850,preset:'moderate'}).buffer);const r=analyzeJointSanger(d);expect(r.joint!.heldOutR2).toBeGreaterThan(.9);expect(r.joint!.groups.some(g=>g.netSize===-85&&g.reference==='Ref0'&&(g.events.some(e=>e.start===303&&e.end===388)||g.alternatives.some(a=>a.events.some(e=>e.start===303&&e.end===388))))).toBe(true);},120000);
it('supports single-reference two-target decomposition on the primer axis',()=>{const d=fixture();d.amplicons=d.amplicons.slice(0,1);const a=d.amplicons[0].sequence.slice(primer.length),edited=a.slice(0,260)+a.slice(261,370)+'T'+a.slice(370);d.control=parseAbif(simulate([{sequence:a,weight:1}],{seed:1234,trim:10,length:850,preset:'moderate'}).buffer);d.sample=parseAbif(simulate([{sequence:a,weight:.5},{sequence:edited,weight:.5}],{seed:1235,trim:8,length:850,preset:'moderate'}).buffer);const r=analyzeJointSanger(d);expect(r.joint!.targets).toHaveLength(2);expect(r.joint!.targets.every(t=>Math.abs(t.edited!-.5)<.03)).toBe(true);},120000);

it('retains asymmetric deletion alleles across new backbones and sizes after extending the usable fit range',()=>{
 for(const [seed,length] of [[704,31],[811,61],[920,85]]){
  const d=fixture(false,true,seed);d.settings.maxDeletion=120;d.settings.window=60;const a=d.amplicons[0].sequence.slice(primer.length),b=d.amplicons[1].sequence.slice(primer.length),start=300-Math.floor(length*.29),edited=a.slice(0,start)+a.slice(start+length),expected=primer+edited;
  d.sample=parseAbif(simulate([{sequence:a,weight:.3},{sequence:edited,weight:.3},{sequence:b,weight:.4}],{seed:seed+10000,trim:8,length:850,preset:'moderate'}).buffer);
  const result=analyzeJointSanger(d),correctMass=result.joint!.groups.filter(g=>g.reference==='Ref0'&&g.netSize===-length).reduce((sum,g)=>{const events=g.alternatives.length?g.alternatives.map(a=>a.events):[g.events];return sum+(events.some(es=>applyJointEvents(d.amplicons[0].sequence,es.map(e=>({target:e.target,start:e.start,end:e.end,insert:e.inserted})))===expected)?g.fraction:0);},0);
  expect(correctMass).toBeGreaterThan(.27);expect(correctMass).toBeLessThan(.33);expect(result.mixture!.r2).toBeGreaterThan(.99);
 }
},120000);
