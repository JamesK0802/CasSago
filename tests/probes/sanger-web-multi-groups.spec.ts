import '@angular/compiler';
import {it} from 'vitest';
import {readFileSync,writeFileSync,mkdirSync,appendFileSync} from 'node:fs';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
import {analyzeAutomaticSanger,automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {benchmarkComparison,parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
import {analyzeAlleleGroups} from '../../src/app/workers/core/sanger-allele-groups';
const run=process.env.CASSAGO_WEB_MULTI==='1'?it:it.skip;
run('blind replay of the web-generated three-reference two-target traces',()=>{
 const root='../benchmarks/web-3ref-2target',template=JSON.parse(readFileSync(root+'/inputs/analysis-template.json','utf8'));
 const input={references:template.amplicons.map((r:any)=>({name:r.name,sequence:r.sequence,cut:r.cuts[0],cuts:r.cuts})),primer:template.primer,control:parseSangerTrace(Uint8Array.from(readFileSync(root+'/inputs/SYNTHETIC_control.ab1')).buffer),sample:parseSangerTrace(Uint8Array.from(readFileSync(root+'/inputs/SYNTHETIC_edit.ab1')).buffer),maxDeletion:10,maxInsertion:10,groupRecovery:true};
 writeFileSync(root+'/progress.log','');const start=Date.now(),progress=(s:string)=>appendFileSync(root+'/progress.log',s+'\n'),result=process.env.CASSAGO_AUDIT==='1'?analyzeAlleleGroups({amplicons:input.references,primer:input.primer,control:input.control,sample:input.sample,ampliconConfirmed:true,settings:{...automaticSangerSettings(input),readEnd:1000,editStartRadius:10,maxIndel:10}},progress,a=>writeFileSync(root+'/candidate-audit.json',JSON.stringify({models:a.models,ids:a.ids,weights:a.weights,refs:a.refs.map((r:any)=>({offset:r.offset,reverse:r.reverse,name:r.name}))}))):analyzeAutomaticSanger(input,progress);
 // Truth is read only after blind analysis has returned.
 const truth=parseBenchmarkTruth(readFileSync(root+'/inputs/truth.json','utf8')),comparison=benchmarkComparison(truth,input.references,result);
 mkdirSync(root+'/sequence-engine',{recursive:true});writeFileSync(root+'/sequence-engine/result.json',JSON.stringify({elapsedMs:Date.now()-start,result,comparison}));
 console.log(JSON.stringify({elapsedMs:Date.now()-start,engine:result.diagnostics?.algorithmVersion,warnings:result.warnings,recovered:comparison.filter(r=>r.truth===5&&(r.estimated??0)>0).length,accurate:comparison.filter(r=>r.truth===5&&(r.estimated??0)>0&&Math.abs(r.delta!)<=2).length}));
},1200000);
