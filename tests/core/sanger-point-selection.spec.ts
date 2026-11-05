import{it,expect}from'vitest';import{readFileSync}from'node:fs';
import{parseAbif}from'../../src/app/utils/abif';
import{parseAutomaticSangerTemplate}from'../../src/app/utils/sanger-automatic-template';
import{analyzeAutomaticSanger}from'../../src/app/workers/core/sanger-automatic';
const root=new URL('../fixtures/sanger/calibrated-design/',import.meta.url);
function load(name:string){const b=readFileSync(new URL(name,root));return parseAbif(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),name);}
it('removes a nonpredictive +4 artifact without discarding the real 3% +6 class or measured cycles',()=>{
 const c=JSON.parse(readFileSync(new URL('prediction-null.json',root),'utf8')),t=parseAutomaticSangerTemplate(c.template),sample=load(c.files.sample),r=analyzeAutomaticSanger({references:t.references.map((a:any)=>({name:a.name,sequence:a.sequence,cut:Number(a.targets)})),primer:t.primer,maxDeletion:t.maxDeletion,maxInsertion:t.maxInsertion,control:load(c.files.control),sample});
 expect(r.mixture!.components.find(c=>c.size===4)!.fraction).toBe(0);expect(r.mixture!.components.find(c=>c.size===6)!.fraction).toBeGreaterThan(.01);expect(r.mixture!.bases).toBe(sample.peaks.length);
 const truth=new Map<number,number>(c.truth.references[0].net.map((x:any)=>[x.size,x.fraction]));for(const a of r.mixture!.components){if((truth.get(a.size)||0)>=.02)expect(a.fraction).toBeGreaterThan(.01);if(!truth.has(a.size))expect(a.fraction).toBeLessThan(.01);}
},60000);
it('retains a weak homeolog allocation when the length is observable and independent prediction contains useful information',()=>{
 const c=JSON.parse(readFileSync(new URL('reference-allocation.json',root),'utf8')),r=analyzeAutomaticSanger({...c,control:load(c.files.control),sample:load(c.files.sample)}),a=r.pooled!.amplicons.find(x=>x.name==='A')!;
 expect(a.components.find(x=>x.size===-7)!.fraction).toBeGreaterThan(.01);expect(a.components.find(x=>x.size===2)!.fraction).toBeGreaterThan(.01);
 for(const truth of c.truth.references){const h=r.pooled!.amplicons.find(x=>x.name===truth.name)!;for(const t of truth.net)if(t.size&&t.fraction>=.02)expect(h.components.find(x=>x.size===t.size)!.fraction).toBeGreaterThan(.01);}
},60000);
