import '@angular/compiler';
import {it} from 'vitest';
import {readFileSync,writeFileSync,mkdirSync,appendFileSync} from 'node:fs';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
import {analyzeAutomaticSanger} from '../../src/app/workers/core/sanger-automatic';
import {benchmarkComparison,parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
const run=process.env.CASSAGO_PERF==='1'?it:it.skip;
run('blind speed and target-profile accuracy replay',()=>{
 const root='../benchmarks/six-target-profiles',label=process.env.CASSAGO_PERF_LABEL||'after',out='../benchmarks/performance-'+label;mkdirSync(out,{recursive:true});const template=JSON.parse(readFileSync(root+'/inputs/analysis-template.json','utf8')),input={references:template.amplicons.map((a:any)=>({name:a.name,sequence:a.sequence,cut:a.cuts[0],cuts:a.cuts})),primer:template.primer,control:parseSangerTrace(Uint8Array.from(readFileSync(root+'/inputs/SYNTHETIC_control.ab1')).buffer),sample:parseSangerTrace(Uint8Array.from(readFileSync(root+'/inputs/SYNTHETIC_edit.ab1')).buffer),maxDeletion:10,maxInsertion:10,groupRecovery:true};
 const started=Date.now();writeFileSync(out+'/progress.log','');const result=analyzeAutomaticSanger(input,stage=>appendFileSync(out+'/progress.log',JSON.stringify({ms:Date.now()-started,stage})+'\n')),elapsedMs=Date.now()-started;
 // Benchmark labels are read only after analysis. No recipe enters inference.
 const truth=parseBenchmarkTruth(readFileSync(root+'/inputs/truth.json','utf8')),comparison=benchmarkComparison(truth,input.references,result),positive=comparison.filter(r=>r.scope==='group'&&(r.truth??0)>0),variants=positive.filter(r=>r.group!=='WT');
 const summary={elapsedMs,algorithm:result.diagnostics?.algorithmVersion,rmse:result.mixture?.rmse,passesChecks:result.mixture?.passesChecks,iterations:result.joint?.iterations,candidates:result.joint?.candidateModels,groups:result.joint?.selectedModels,profileGroups:positive.length,profileMatches:positive.filter(r=>(r.estimated??0)>0).length,profileAccurate2pp:positive.filter(r=>(r.estimated??0)>0&&Math.abs(r.delta!)<=2).length,variantGroups:variants.length,variantMatches:variants.filter(r=>(r.estimated??0)>0).length,variantAccurate2pp:variants.filter(r=>(r.estimated??0)>0&&Math.abs(r.delta!)<=2).length,warnings:result.warnings};writeFileSync(out+'/summary.json',JSON.stringify(summary,null,2));writeFileSync(out+'/result.json',JSON.stringify({result,comparison}));
},1200000);
