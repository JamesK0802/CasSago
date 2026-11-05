import '@angular/compiler';
import {it} from 'vitest';
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {analyzeConditionalGroups} from '../../src/app/workers/core/sanger-conditional-groups';
import {automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
import {parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
import {assessSangerBenchmark} from '../../src/app/utils/sanger-benchmark-assessment';
(process.env.CASSAGO_CASE_REPLAY==='1'?it:it.skip)('diagnoses saved cases without counting them as prospective passes',()=>{
 const folders=process.env.CASSAGO_CASE_INPUT!.split(','),out=process.env.CASSAGO_CASE_OUTPUT!;mkdirSync(out,{recursive:true});const summary:any[]=[];
 for(const [i,folder] of folders.entries()){
 const template=JSON.parse(readFileSync(folder+'/template.json','utf8')),truth=parseBenchmarkTruth(readFileSync(folder+'/truth.json','utf8')),read=(name:string)=>{const b=readFileSync(folder+'/'+name);return parseSangerTrace(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),name);},input={references:template.amplicons.map((r:any)=>({...r,cut:r.cuts[0]})),primer:template.primer,control:read('control.ab1'),sample:read('edit.ab1'),maxDeletion:10,maxInsertion:10},start=Date.now();
 try{const result=analyzeConditionalGroups({amplicons:input.references,primer:input.primer,sample:input.sample,control:input.control,ampliconConfirmed:true,settings:{...automaticSangerSettings(input),readEnd:1000,editStartRadius:3,maxIndel:10}},stage=>writeFileSync(out+'/progress.json',JSON.stringify({i,stage,ms:Date.now()-start})));const a=assessSangerBenchmark(truth,input.references,result);summary.push({folder,ms:Date.now()-start,...a,rows:undefined,certificates:undefined});writeFileSync(out+'/'+i+'-result.json',JSON.stringify(result));}catch(e){summary.push({folder,ms:Date.now()-start,error:String(e)});}writeFileSync(out+'/summary.json',JSON.stringify(summary,null,2));
 }
},1800000);
