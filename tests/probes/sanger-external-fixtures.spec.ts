import '@angular/compiler';
import {it,expect} from 'vitest';
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {fiveIndelSequenceIdentifiabilityFixture} from '../fixtures/sanger-identifiability-fixture';
import {profileFixture} from '../fixtures/sanger-profile-fixture';
import {generateSangerSynth,synthReferences} from '../../src/app/utils/sanger-synthetic-benchmark';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
import {writeAbif} from '../../src/app/utils/sanger-analysis-demo';
import {TRACE_BASES} from '../../src/app/models/sanger-analysis';
const run=process.env.CASSAGO_EXTERNAL_FIXTURES==='1'?it:it.skip;
run('creates raw-channel-preserving single-reference single-target crops for external tools',()=>{
 const folder='../benchmarks/external-tools-20261104';mkdirSync(folder,{recursive:true});const manifest:any[]=[];
 for(const key of ['ambiguity-A','ambiguity-B','original-mixed','five-deletions']){
  const f=profileFixture(3),h=key.startsWith('ambiguity')?0:1,target=key==='five-deletions'?1:0,result={...f.result,references:[f.result.references[h]]},candidate={...f.candidate,sites:[f.candidate.sites[h]],partner:{...f.candidate.partner!,sites:[f.candidate.partner!.sites[h]],products:[f.candidate.partner!.products[h]]}};
  let ref:any;if(key.startsWith('ambiguity')){const fixture=fiveIndelSequenceIdentifiabilityFixture();ref=structuredClone((key.endsWith('A')?fixture.refs:fixture.alternate)[h]);}else ref=structuredClone(JSON.parse(readFileSync('../benchmarks/variable-percent-20261104/5-both-recipe.json','utf8')).recipe[h]);
  ref.share=100;ref.profiles.forEach((p:any)=>{if(p.target!==target)p.alleles=[];});
  if(key==='five-deletions')ref.profiles[target].alleles=[1,2,3,5,8].map((length,i)=>({id:i,percent:[6.3,8.1,11.7,14.2,9.7][i],events:[{kind:'deletion',target,offset:[-3,-1,0,1,3][i],length,sequence:''}]}));
  const d=generateSangerSynth(result,candidate,[ref],20261104),cut=ref.cuts[target]-ref.primerLength,start=Math.max(0,cut-180),end=Math.min(d.cycles,cut+230),localCut=cut-start,guide=ref.sequence.slice(ref.cuts[target]-17,ref.cuts[target]+3);
  const crop=(buffer:ArrayBuffer,name:string)=>{const t=parseSangerTrace(buffer),left=Math.max(0,t.peaks[start]-7),right=Math.min(t.channels.A.length,t.peaks[end-1]+8),peaks=Array.from(t.peaks.slice(start,end),p=>p-left),channels=TRACE_BASES.map(b=>Uint16Array.from(t.channels[b].slice(left,right)));const out=writeAbif(t.sequence.slice(start,end),peaks,channels,name,'ACGT',60);const check=parseSangerTrace(out);for(const base of TRACE_BASES)expect(Array.from(check.channels[base])).toEqual(Array.from(t.channels[base].slice(left,right)));return out;};
  const control=crop(d.control,'control'),edit=crop(d.edit,'edit');writeFileSync(`${folder}/${key}-control.ab1`,new Uint8Array(control));writeFileSync(`${folder}/${key}-edit.ab1`,new Uint8Array(edit));writeFileSync(`${folder}/${key}-truth.json`,JSON.stringify({recipe:d.recipe,truth:d.truth},null,2));manifest.push({key,reference:ref.sourceReference,target:target+1,cutBoundaryInCrop:localCut,cropStartCycle:start,cropEndCycle:end,guide,guideSource:'20 reference bases spanning the specified cut (17 before + 3 after); matched target, not PAM validation',controlFile:key+'-control.ab1',editFile:key+'-edit.ab1',editSha256:createHash('sha256').update(new Uint8Array(edit)).digest('hex'),profile:ref.profiles[target]});
 }
 expect(manifest[0].editSha256).toBe(manifest[1].editSha256);writeFileSync(folder+'/manifest.json',JSON.stringify(manifest,null,2));console.log(JSON.stringify(manifest.map(({profile,...r})=>r),null,2));
});
