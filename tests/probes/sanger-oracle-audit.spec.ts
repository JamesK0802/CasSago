import '@angular/compiler';
import {it} from 'vitest';
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {profileFixture} from '../fixtures/sanger-profile-fixture';
import {generateSangerSynth} from '../../src/app/utils/sanger-synthetic-benchmark';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
import {analyzeAlleleGroupsOracle as analyzeAlleleGroups} from '../../src/app/workers/core/sanger-allele-groups-oracle-prototype';
import {automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {benchmarkComparison,parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
import {materialize} from '../../src/app/workers/core/sanger-path-common';
const run=process.env.CASSAGO_ORACLE_AUDIT==='1'?it:it.skip;
run('tests actual single-cut inference on identical pooled traces and pure-reference counterparts',()=>{
 const source='../benchmarks/variable-percent-20261104',folder=process.env.CASSAGO_SINGLE_FOLDER||'../benchmarks/single-ranking-fix';mkdirSync(folder,{recursive:true});const summary:any[]=[];
 for(const target of [0,1])for(const h of [-1,0,1,2]){
  const key='T'+(target+1)+(h<0?'-pooled':'-ref'+h),oldKey='5-'+(target===0?'T2-WT':'T1-WT'),original=JSON.parse(readFileSync(`${source}/${oldKey}-recipe.json`,'utf8')),f=profileFixture(3);let data:any;if(process.env.CASSAGO_ORACLE_CASES&&!process.env.CASSAGO_ORACLE_CASES.split(',').includes(key))continue;
  if(h<0){const template=JSON.parse(readFileSync(`${source}/${oldKey}-template.json`,'utf8'));data={recipe:original.recipe,truth:original.truth,input:{references:template.amplicons.map((a:any)=>({name:a.name,sequence:a.sequence,cut:a.cuts[0],cuts:a.cuts})),primer:template.primer,control:parseSangerTrace(Uint8Array.from(readFileSync(`${source}/${oldKey}-control.ab1`)).buffer),sample:parseSangerTrace(Uint8Array.from(readFileSync(`${source}/${oldKey}-edit.ab1`)).buffer)}};
  }else{const result={...f.result,references:[f.result.references[h]]},candidate={...f.candidate,sites:[f.candidate.sites[h]],partner:{...f.candidate.partner!,sites:[f.candidate.partner!.sites[h]],products:[f.candidate.partner!.products[h]]}},ref=structuredClone(original.recipe[h]);ref.share=100;data=generateSangerSynth(result,candidate,[ref],20261135);}
  const references=data.input.references.map((r:any)=>({...r,cut:r.cuts[target],cuts:[r.cuts[target]]})),settings=automaticSangerSettings({...data.input,references,maxDeletion:10,maxInsertion:10});let audit:any;const start=Date.now(),result=analyzeAlleleGroups({amplicons:references,primer:data.input.primer,control:data.input.control,sample:data.input.sample,ampliconConfirmed:true,settings:{...settings,readEnd:1000,editStartRadius:10,maxIndel:10}},()=>{},state=>audit=state),elapsedMs=Date.now()-start;
  // Project truth only after inference. The selected cut is target 1 in this run.
  const truth={recipe:data.recipe.map((r:any)=>({...r,cuts:[r.cuts[target]],targetNames:[r.targetNames[target]],profiles:[{...r.profiles[target],target:0,alleles:r.profiles[target].alleles.map((a:any)=>({...a,events:a.events.map((e:any)=>({...e,target:0}))}))}]})),truth:data.truth.map((r:any)=>({...r,pattern:r.pattern.split('/')[target],events:r.events.filter((e:any)=>e.target===target).map((e:any)=>({...e,target:0})),description:r.events.some((e:any)=>e.target===target)?r.description:'WT'}))};
  const rows=benchmarkComparison(parseBenchmarkTruth(truth),references,result).filter(r=>r.scope==='group'&&(r.truth??0)>0&&r.group!=='WT'),candidates:any[]=[];
  for(let i=0;i<data.recipe.length;i++){const ref=data.recipe[i],prepared=audit.refs[i];for(const allele of ref.profiles[target].alleles){const e=allele.events[0],start=(e.kind==='substitution'?e.selectionStart:ref.cuts[target]+e.offset)-prepared.offset,end=start+(e.kind==='insertion'?0:e.kind==='deletion'?e.length:e.before.length),sequence=materialize(prepared.sequence,[{start,end,inserted:e.kind==='deletion'?'':e.sequence}]),matches=audit.models.flatMap((m:any,id:number)=>m.reference===i&&m.sequence===sequence?[id]:[]),selected=matches.filter((id:number)=>audit.ids.includes(id));candidates.push({reference:ref.name,group:allele.id,kind:e.kind,truth:allele.percent,candidates:matches.length,selected:selected.length,estimatedSamplePercent:selected.reduce((s:number,id:number)=>s+audit.weights[audit.ids.indexOf(id)],0)*100});}}
  const row={key,target:target+1,references:references.length,elapsedMs,variants:rows.length,sequenceMatches:rows.filter(r=>(r.estimated??0)>0).length,accurate2pp:rows.filter(r=>(r.estimated??0)>0&&Math.abs(r.delta!)<=2).length,rmse:result.mixture?.rmse,passesChecks:result.mixture?.passesChecks,candidates};summary.push(row);writeFileSync(folder+'/summary.json',JSON.stringify(summary,null,2));writeFileSync(`${folder}/${key}-result.json`,JSON.stringify(result));writeFileSync(`${folder}/${key}-comparison.json`,JSON.stringify(rows));console.log(JSON.stringify(row));
 }
},1200000);
