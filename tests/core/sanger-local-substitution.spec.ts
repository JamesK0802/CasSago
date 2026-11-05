import '@angular/compiler';
import {expect,it} from 'vitest';
import {profileFixture} from '../fixtures/sanger-profile-fixture';
import {generateSangerSynth,synthReferences} from '../../src/app/utils/sanger-synthetic-benchmark';
import {analyzeAutomaticSanger} from '../../src/app/workers/core/sanger-automatic';
it('recovers a single 8.66% SNP that tied behind long indel candidates under raw correlation ranking',()=>{
 const f=profileFixture(3),h=1,result={...f.result,references:[f.result.references[h]]},candidate={...f.candidate,sites:[f.candidate.sites[h]],partner:{...f.candidate.partner!,sites:[f.candidate.partner!.sites[h]],products:[f.candidate.partner!.products[h]]}},refs=synthReferences(result,candidate),r=refs[0],start=r.cuts[1]+1;expect(r.sequence[start]).toBe('G');r.alleles=[{id:1,percent:8.66,events:[{kind:'substitution',target:1,offset:1,length:1,selectionStart:start,before:'G',sequence:'A'}]}];
 const data=generateSangerSynth(result,candidate,refs),references=data.input.references.map(a=>({...a,cut:a.cuts![1],cuts:[a.cuts![1]]})),analysis=analyzeAutomaticSanger({...data.input,references,groupRecovery:true}),variant=analysis.annotations![0].components.find(c=>c.events?.some(e=>e.start===start&&e.end===start+1&&e.insertedSequence==='A'));
 expect(variant).toBeDefined();expect(variant!.fraction).toBeCloseTo(.0866,4);expect(analysis.mixture!.rmse).toBeLessThan(1e-4);
},30000);
