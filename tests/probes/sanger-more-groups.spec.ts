import '@angular/compiler';
import {it} from 'vitest';
import {mkdirSync,writeFileSync} from 'node:fs';
import {profileFixture} from '../fixtures/sanger-profile-fixture';
import {generateSangerSynth,randomSynthEvent} from '../../src/app/utils/sanger-synthetic-benchmark';
import {analyzeAlleleGroups} from '../../src/app/workers/core/sanger-allele-groups';
import {automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {benchmarkComparison,parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
const run=process.env.CASSAGO_MORE_GROUPS==='1'?it:it.skip;
run('blindly tests more groups per target across three references and two targets',()=>{
 const folder=process.env.CASSAGO_GROUP_FOLDER||'../benchmarks/more-groups-20261103';mkdirSync(folder,{recursive:true});const summary=[];
 for(const count of (process.env.CASSAGO_COUNTS||'3,5,8').split(',').map(Number)){
  const seed=20261130+count;let state=seed;const random=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state/4294967296;},f=profileFixture(3);
  for(const ref of f.refs)for(const profile of ref.profiles!){const seen=new Set<string>();profile.alleles=[];for(let i=0;i<count;i++){let event;do{event=randomSynthEvent(ref,random,profile.target);}while(seen.has(JSON.stringify(event)));seen.add(JSON.stringify(event));profile.alleles.push({id:profile.target*100+i,percent:count===3?(i===2?20:15):50/count,events:[event]});}}
  const d=generateSangerSynth(f.result,f.candidate,f.refs,seed),input={amplicons:d.input.references,primer:d.input.primer!,control:d.input.control,sample:d.input.sample,ampliconConfirmed:true,settings:{...automaticSangerSettings(d.input),readEnd:1000,editStartRadius:10,maxIndel:10}};
  writeFileSync(`${folder}/${count}-recipe.json`,JSON.stringify({recipe:d.recipe,truth:d.truth,mixing:d.mixing}));writeFileSync(`${folder}/${count}-control.ab1`,new Uint8Array(d.control));writeFileSync(`${folder}/${count}-edit.ab1`,new Uint8Array(d.edit));writeFileSync(`${folder}/${count}-template.json`,JSON.stringify(d.template));
  const start=Date.now(),result=analyzeAlleleGroups(input),elapsedMs=Date.now()-start;
  // Truth is used only after inference, to score independently generated inputs.
  const truth=parseBenchmarkTruth({truth:d.truth,recipe:d.recipe}),comparison=benchmarkComparison(truth,d.input.references,result),rows=comparison.filter(r=>r.scope==='group'&&(r.truth??0)>0),variants=rows.filter(r=>!r.group.includes('WT'));
  const row={seed,groupsPerTarget:count,targets:6,molecules:d.truth.length,elapsedMs,profileMatches:rows.filter(r=>(r.estimated??0)>0).length,accurate2pp:rows.filter(r=>(r.estimated??0)>0&&Math.abs(r.delta!)<=2).length,total:rows.length,variantMatches:variants.filter(r=>(r.estimated??0)>0).length,variantAccurate2pp:variants.filter(r=>(r.estimated??0)>0&&Math.abs(r.delta!)<=2).length,variantTotal:variants.length,rmse:result.mixture?.rmse,passesChecks:result.mixture?.passesChecks,iterations:result.joint?.iterations,linkageAliases:result.diagnostics?.substitutionLinkage?.length||0,warnings:result.warnings};summary.push(row);writeFileSync(`${folder}/summary.json`,JSON.stringify(summary,null,2));writeFileSync(`${folder}/${count}-result.json`,JSON.stringify(result));writeFileSync(`${folder}/${count}-comparison.json`,JSON.stringify(comparison));console.log(JSON.stringify(row));
 }
},1200000);
