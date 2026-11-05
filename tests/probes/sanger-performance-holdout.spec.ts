import '@angular/compiler';
import {it} from 'vitest';
import {mkdirSync,writeFileSync} from 'node:fs';
import {profileFixture} from '../fixtures/sanger-profile-fixture';
import {generateSangerSynth,randomSynthEvent} from '../../src/app/utils/sanger-synthetic-benchmark';
import {analyzeAlleleGroups} from '../../src/app/workers/core/sanger-allele-groups';
import {automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {benchmarkComparison,parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
const run=process.env.CASSAGO_HOLDOUT==='1'?it:it.skip;
run('compares unchanged and optimized engines on previously unused target recipes and sequences',async()=>{
 const {analyzeAlleleGroups:baseline}=await import('../../../benchmarks/performance-holdout/baseline-engine');
 const output=[];for(const seed of [20261104,20261105]){
  let state=seed;const random=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state/4294967296;},f=profileFixture(3);
  for(const ref of f.refs)for(const profile of ref.profiles!)profile.alleles=[0,1].map((_,i)=>({id:10*profile.target+i,percent:(profile.target?15:10)+10*i,events:[randomSynthEvent(ref,random,profile.target)]}));
  const d=generateSangerSynth(f.result,f.candidate,f.refs,seed),input={amplicons:d.input.references,primer:d.input.primer!,control:d.input.control,sample:d.input.sample,ampliconConfirmed:true,settings:{...automaticSangerSettings(d.input),readEnd:1000,editStartRadius:10,maxIndel:10}},results=[];
  for(const [label,engine]of [['before',baseline],['after',analyzeAlleleGroups]]as const){const start=Date.now(),result=engine(input),elapsedMs=Date.now()-start,truth=parseBenchmarkTruth({truth:d.truth,recipe:d.recipe}),comparison=benchmarkComparison(truth,d.input.references,result),rows=comparison.filter(r=>r.scope==='group'&&(r.truth??0)>0);results.push({label,elapsedMs,profileMatches:rows.filter(r=>(r.estimated??0)>0).length,accurate2pp:rows.filter(r=>(r.estimated??0)>0&&Math.abs(r.delta!)<=2).length,total:rows.length,rmse:result.mixture?.rmse,passesChecks:result.mixture?.passesChecks});}
  output.push({seed,results});writeFileSync('../benchmarks/performance-holdout/summary.json',JSON.stringify(output,null,2));
 }
},1200000);
