import {it,expect} from 'vitest';
import {readFileSync} from 'node:fs';
import {designSangerPrimers} from '../../src/app/workers/core/sanger-primer-designer';
import {makeAnalysisTemplate} from '../../src/app/utils/sanger-design-export';
import {parseAutomaticSangerTemplate} from '../../src/app/utils/sanger-automatic-template';
import {parseAbif} from '../../src/app/utils/abif';
import {analyzeAutomaticSanger} from '../../src/app/workers/core/sanger-automatic';
import {SANGER_DEFAULTS} from '../../src/app/models/sanger-design';
const base=new URL('../fixtures/sanger/design-pipeline/',import.meta.url);
// Group recovery and handoff wiring are separate from the unchanged 5pp gate,
// which the benchmark still reports as failing for these component-based traces.
it.each(['forward','reverse'])('uses designed PCR products and read coordinates through the Sanger workflow: %s',label=>{
 const c=JSON.parse(readFileSync(new URL(label+'.json',base),'utf8'));
 const design=designSangerPrimers(c.designSource.references,c.designSource.designSettings),candidate=design.candidates.find(p=>p.id===c.designSource.candidate.id)!;
 expect(candidate).toBeDefined();
 expect(candidate.sequence).toBe(c.primer);expect(candidate.partner).not.toBeNull();
 expect(candidate.minOffset).toBeGreaterThanOrEqual(design.settings.minOffset);
 expect(candidate.indelOffsetOverlap).toBe(false);
 expect(candidate.sites.map(s=>s.cutOffsets[0])).toEqual(c.designSource.readCutOffsets);
 const input=parseAutomaticSangerTemplate(makeAnalysisTemplate(design,candidate));
 const references=input.references.map((r:any)=>({name:r.name,sequence:r.sequence,cut:Number(r.targets)}));
 expect(references).toEqual(c.amplicons);
 const load=(role:string)=>{const b=readFileSync(new URL(c.files[role],base));return parseAbif(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),c.files[role]);};
 const result=analyzeAutomaticSanger({references,primer:input.primer,maxDeletion:input.maxDeletion,maxInsertion:input.maxInsertion,control:load('control'),sample:load('sample')});
 for(const truth of c.truth.references){
  const h=result.pooled!.amplicons.find(a=>a.name===truth.name)!,expected=new Map<number,number>(truth.net.map((a:any)=>[a.size,a.fraction]));
  for(const [size,fraction]of expected)if(fraction>=.02)expect(h.components.find(a=>a.size===size)!.fraction).toBeGreaterThanOrEqual(.01);
  for(const a of h.components)if(!expected.has(a.size))expect(a.fraction).toBeLessThan(.01);
 }
},15000);
it('rejects same-offset SNP-only pools in the shared-Sanger primer designer',()=>{
 const old=JSON.parse(readFileSync(new URL('../homeolog-single/similar-99.json',base),'utf8'));
 const design=designSangerPrimers(old.amplicons.map((a:any)=>({name:a.name,sequence:a.sequence,cuts:[a.cut]})),SANGER_DEFAULTS);
 expect(design.stats.commonUnique).toBeGreaterThan(0);
 expect(design.stats.offsetPassed).toBe(0);
 expect(design.candidates).toHaveLength(0);
},15000);
