import '@angular/compiler';
import {it} from 'vitest';
import {writeFileSync,mkdirSync} from 'node:fs';
import {profileFixture} from '../fixtures/sanger-profile-fixture';
import {synthReferences,generateSangerSynth,randomSynthEvent} from '../../src/app/utils/sanger-synthetic-benchmark';
import {initSynthProfiles} from '../../src/app/utils/sanger-synth-profiles';
import {analyzeAlleleGroups} from '../../src/app/workers/core/sanger-allele-groups';
import {automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {benchmarkComparison,parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
const run=process.env.CASSAGO_FRESH_CONTEXT==='1'?it:it.skip;
run('tests the engine on three independently randomized reference contexts with five edits per target',()=>{
 let state=20261141;const random=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state/4294967296;},f=profileFixture(3),folder=process.env.CASSAGO_FRESH_FOLDER||'../benchmarks/fresh-context-20261104';mkdirSync(folder,{recursive:true});
 f.result.references.forEach((ref,h)=>{const left=f.candidate.sites[h],right=f.candidate.partner!.sites[h],start=Math.min(left.start,right.start),end=Math.max(left.end,right.end),bases=ref.sequence.split('');for(let p=start;p<end;p++)if(!(p>=left.start&&p<left.end)&&!(p>=right.start&&p<right.end))bases[p]='ACGT'[Math.floor(random()*4)];ref.sequence=bases.join('');ref.name='Fresh-'+(h+1);});
 const refs=synthReferences(f.result,f.candidate);for(const ref of refs)for(const p of initSynthProfiles(ref))p.alleles=Array.from({length:5},(_,i)=>({id:i,percent:10,events:[randomSynthEvent(ref,random,p.target)]}));
 const d=generateSangerSynth(f.result,f.candidate,refs,20261141),input={amplicons:d.input.references,primer:d.input.primer!,control:d.input.control,sample:d.input.sample,ampliconConfirmed:true,settings:{...automaticSangerSettings(d.input),readEnd:1000,editStartRadius:10,maxIndel:10}},start=Date.now(),result=analyzeAlleleGroups(input),elapsedMs=Date.now()-start;
 const truth=parseBenchmarkTruth({truth:d.truth,recipe:d.recipe}),comparison=benchmarkComparison(truth,d.input.references,result),rows=comparison.filter(r=>r.scope==='group'&&(r.truth??0)>0&&r.group!=='WT');
 writeFileSync(folder+'/summary.json',JSON.stringify({elapsedMs,variants:rows.length,sequenceMatches:rows.filter(r=>(r.estimated??0)>0).length,accurate2pp:rows.filter(r=>(r.estimated??0)>0&&Math.abs(r.delta!)<=2).length,passesChecks:result.mixture?.passesChecks,rmse:result.mixture?.rmse,warnings:result.warnings},null,2));writeFileSync(folder+'/result.json',JSON.stringify(result));writeFileSync(folder+'/comparison.json',JSON.stringify(comparison));writeFileSync(folder+'/recipe.json',JSON.stringify({recipe:d.recipe,truth:d.truth,mixing:d.mixing}));writeFileSync(folder+'/template.json',JSON.stringify(d.template));writeFileSync(folder+'/control.ab1',new Uint8Array(d.control));writeFileSync(folder+'/edit.ab1',new Uint8Array(d.edit));
},600000);
