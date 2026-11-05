import {it,expect} from 'vitest';
import {readFileSync} from 'node:fs';
import {parseAbif} from '../../src/app/utils/abif';
import {analyzeAutomaticSanger} from '../../src/app/workers/core/sanger-automatic';
const analyzePooledSanger=(d:any)=>analyzeAutomaticSanger({references:d.amplicons,primer:d.primer,sample:d.sample,control:d.control,maxDeletion:d.settings.maxIndel,maxInsertion:d.settings.maxIndel});
const base=new URL('../fixtures/sanger/homeolog-single/',import.meta.url);
function input(label:string){
 const c=JSON.parse(readFileSync(new URL(label+'.json',base),'utf8'));
 const load=(f:string)=>{const b=readFileSync(new URL(f,base));return parseAbif(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),f);};
 return {c,data:{amplicons:c.amplicons,primer:c.primer,ampliconConfirmed:true,control:load(c.files.control),sample:load(c.files.sample),settings:c.settings}};
}
it.each(['similar-99','structural-three'])('recovers each homeolog share and simultaneous indel spectrum: %s',label=>{
 const {c,data}=input(label),r=analyzePooledSanger(data);
 expect(r.mixture!.passesChecks).toBe(true);
 for(const truth of c.truth.references){
  const h=r.pooled!.amplicons.find(x=>x.name===truth.name)!;
  expect(Math.abs(h.share-truth.share)).toBeLessThan(.05);
  const expected=new Map<number,number>(truth.net.map((a:any)=>[a.size,a.fraction]));
  const tv=h.components.reduce((sum,a)=>sum+Math.abs(a.fraction-(expected.get(a.size)||0)),0)/2;
  expect(tv).toBeLessThan(.05);
 }
});
it('returns unstable homeolog estimates with QC warnings instead of hiding the numbers',()=>{
 const {data}=input('unstable-three');
 for(let i=430;i<434;i++){const p=data.sample.peaks[i];for(const b of ['A','C','G','T'] as const)for(let k=p-6;k<=p+6;k++)if(b==='C')data.sample.channels[b][k]+=Math.round(25000*Math.exp(-(((k-p)/2)**2)));}
 const r=analyzePooledSanger(data);
 expect(r.mixture!.r2).toBeGreaterThan(.9);
 expect(r.diagnostics!.discordantPositions!).toBeGreaterThan(0);
 expect(r.mixture!.passesChecks).toBe(false);
 expect(r.annotations!.every(a=>!a.estimable)).toBe(true);
 expect(r.pooled!.amplicons.every(a=>a.indel!==null&&a.wt!==null&&a.inFrame!==null&&a.outOfFrame!==null)).toBe(true);
 expect(r.pooled!.amplicons.every(a=>a.components.some(c=>c.fraction>0))).toBe(true);
});
it('returns flagged estimates when a homeolog length difference aliases an edited class',()=>{
 const {data}=input('phase-alias'),r=analyzePooledSanger(data);expect(r.mixture!.passesChecks).toBe(false);expect(r.warnings.join(' ')).toMatch(/rank-deficient|poorly conditioned/);expect(r.pooled!.amplicons.every(a=>a.components.every(c=>Number.isFinite(c.fraction)))).toBe(true);
});
it('reports no edit in an unedited pool even when control and sample reference shares differ',()=>{
 const {data}=input('unedited'),r=analyzePooledSanger(data);
 expect(r.mixture!.passesChecks).toBe(true);
 expect(r.mixture!.indel).toBeLessThan(.01);
});

// Retain the original 5 pp gate on the public automatic analysis contract.
it('recovers unstable-three under automatic full-trace analysis within the original 5 pp gate',()=>{
 const {c,data}=input('unstable-three'),r=analyzePooledSanger(data);
 for(const truth of c.truth.references){const h=r.pooled!.amplicons.find(x=>x.name===truth.name)!;const expected=new Map<number,number>(truth.net.map((a:any)=>[a.size,a.fraction]));const tv=h.components.reduce((sum,a)=>sum+Math.abs(a.fraction-(expected.get(a.size)||0)),0)/2;expect(tv).toBeLessThan(.05);}
});


it.each(['distributed-difficult','junction-only-context'])('known accuracy limit: %s must meet the unchanged 5 pp conditional TV gate',label=>{
 const {c,data}=input(label),r=analyzePooledSanger(data);
 for(const truth of c.truth.references){const h=r.pooled!.amplicons.find(x=>x.name===truth.name)!;const expected=new Map<number,number>(truth.net.map((a:any)=>[a.size,a.fraction]));const tv=h.components.reduce((sum,a)=>sum+Math.abs(a.fraction-(expected.get(a.size)||0)),0)/2;expect(tv).toBeLessThan(.05);}
});

// Bounded integration computation; group/fraction assertions remain unchanged.
it('retains the true groups in a noisy pool without adding unsupported -4 or +9 lengths',()=>{
 const {c,data}=input('distributed-difficult'),r=analyzePooledSanger(data);
 for(const truth of c.truth.references){
  const h=r.pooled!.amplicons.find(x=>x.name===truth.name)!;
  const expected=new Set<number>(truth.net.map((a:any)=>a.size));
  for(const a of truth.net)if(a.fraction>=.02)expect(h.components.find(x=>x.size===a.size)!.fraction).toBeGreaterThanOrEqual(.01);
  for(const a of h.components)if(!expected.has(a.size))expect(a.fraction).toBeLessThan(.01);
 }
},30000);
