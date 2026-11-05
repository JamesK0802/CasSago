import '@angular/compiler';import{it}from'vitest';import{readFileSync,writeFileSync,mkdirSync}from'node:fs';
import{parseSangerTrace}from'../../src/app/utils/sanger-trace';
import{analyzeAutomaticSanger}from'../../src/app/workers/core/sanger-automatic';import{benchmarkComparison,parseBenchmarkTruth}from'../../src/app/utils/sanger-benchmark';
(process.env.CASSAGO_RESEARCH_AUDIT==='1'?it:it.skip)('replays saved measured traces through blind inference',()=>{
 const base='../benchmarks/',saved=JSON.parse(readFileSync(base+'point-display-20261104/result-final.json','utf8')),folder=base+'joint-refinement-20261104/'+(process.env.CASSAGO_REPLAY_LABEL||'alternating');mkdirSync(folder,{recursive:true});
 for(const [key,file] of [['control','SYNTHETIC_control.ab1'],['sample','SYNTHETIC_edit.ab1']]as const){const bytes=readFileSync(base+'web-demo-five-each-20261104/inputs/'+file);saved.input[key]=parseSangerTrace(bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength),file);}
 const before=Date.now(),result=analyzeAutomaticSanger(saved.input,stage=>writeFileSync(folder+'/progress.json',JSON.stringify({stage,ms:Date.now()-before}))),elapsedMs=Date.now()-before;
 const rows=benchmarkComparison(parseBenchmarkTruth(saved.syntheticTruth),saved.input.references,result),gt=rows.filter(r=>r.scope==='group'&&(r.truth??0)>0&&r.group!=='WT'),summary={elapsedMs,algorithm:result.diagnostics?.algorithmVersion,detected:gt.filter(r=>(r.estimated??0)>0).length,accurate:gt.filter(r=>(r.estimated??0)>0&&Math.abs(r.delta!)<=2).length,total:gt.length,failures:gt.filter(r=>!(r.estimated??0)||Math.abs(r.delta!)>2)};
 writeFileSync(folder+'/result.json',JSON.stringify(result));writeFileSync(folder+'/summary.json',JSON.stringify(summary,null,2));console.log(JSON.stringify(summary));
},600000);
