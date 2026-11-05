import{it,expect}from'vitest';
import{readFileSync}from'node:fs';
import{parseAbif}from'../../src/app/utils/abif';
import{parseAutomaticSangerTemplate}from'../../src/app/utils/sanger-automatic-template';
import{sangerAnnotationViews}from'../../src/app/utils/sanger-annotation';
import{analyzeAutomaticSanger}from'../../src/app/workers/core/sanger-automatic';
it('distinguishes weak wrong-homeolog allocations from supported classes when each homeolog has different edits',()=>{
 const root=new URL('../fixtures/sanger/calibrated-design/',import.meta.url),c=JSON.parse(readFileSync(new URL('heterogeneous.json',root),'utf8')),t=parseAutomaticSangerTemplate(c.template),load=(f:string)=>{const b=readFileSync(new URL(f,root));return parseAbif(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),f);};
 const r=analyzeAutomaticSanger({references:t.references.map((a:any)=>({name:a.name,sequence:a.sequence,cut:Number(a.targets)})),primer:t.primer,maxDeletion:t.maxDeletion,maxInsertion:t.maxInsertion,control:load(c.files.control),sample:load(c.files.sample)}),views=sangerAnnotationViews(r.annotations!,0);let wrong=0,real=0;
 for(const h of c.truth.references){const expected=new Map<number,number>(h.net.map((a:any)=>[a.size,a.fraction])),view=views.find(v=>v.name===h.name)!;
  for(const row of view.rows){if(!row.size)continue;if(!expected.has(row.size)&&row.fraction>1e-8){wrong++;expect(row.support?.status).toBe('unresolved');}if((expected.get(row.size)||0)>=.02){real++;expect(row.fraction).toBeGreaterThanOrEqual(.01);}expect(row.positionLabel).not.toMatch(/^Confirmed/);}
 }
 expect(real).toBe(11);
 // Keeping these numeric allocations is deliberate. A display threshold must
 // never upgrade their evidence or turn them into confirmed edit locations.
 expect(r.pooled!.amplicons.every(h=>h.components.every(a=>Number.isFinite(a.fraction)))).toBe(true);
},60000);

it('does not certify a long insertion which fits the training signal but fails relearned block prediction',()=>{
 const root=new URL('../fixtures/sanger/calibrated-design/',import.meta.url),c=JSON.parse(readFileSync(new URL('insertion-overfit.json',root),'utf8')),t=parseAutomaticSangerTemplate(c.template),load=(f:string)=>{const b=readFileSync(new URL(f,root));return parseAbif(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),f);},sample=load(c.files.sample),r=analyzeAutomaticSanger({references:t.references.map((a:any)=>({name:a.name,sequence:a.sequence,cut:Number(a.targets)})),primer:t.primer,maxDeletion:t.maxDeletion,maxInsertion:t.maxInsertion,control:load(c.files.control),sample}),candidate=r.annotations![0].components.find(c=>c.size===23)!;
 expect(r.mixture!.components.find(c=>c.size===23)!.fraction).toBe(0);expect(candidate).toBeUndefined();expect(r.mixture!.bases).toBe(sample.peaks.length);
},60000);
