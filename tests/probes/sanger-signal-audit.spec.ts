import '@angular/compiler';
import {it} from 'vitest';
import {readFileSync,writeFileSync,appendFileSync} from 'node:fs';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
import {parseAutomaticSangerTemplate} from '../../src/app/utils/sanger-automatic-template';
import {automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {analyzeAlleleGroups} from '../../src/app/workers/core/sanger-allele-groups';
import {emissionVector,fitNormalizedEmissionMixture} from '../../src/app/workers/core/sanger-emission';
import {rawSignalScore} from '../../src/app/workers/core/sanger-signal';
import {reverseComplement} from '../../src/app/utils/biology.utils';
const run=process.env.CASSAGO_SIGNAL_AUDIT==='1'?it:it.skip;
run('audit blind fit, then compare extraction and model search against withheld truth',()=>{
 const dir='../benchmarks/fresh-five-groups-20261103/',template=parseAutomaticSangerTemplate(JSON.parse(readFileSync(dir+'template.json','utf8'))),buffer=(f:string)=>{const b=readFileSync(dir+f);return b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength);},input={references:template.references.map(r=>({name:r.name,sequence:r.sequence,cut:Number(r.targets),cuts:[Number(r.targets)]})),primer:template.primer,control:parseSangerTrace(buffer('control.ab1')),sample:parseSangerTrace(buffer('edit.ab1')),maxDeletion:10,maxInsertion:10};let audit:any;
 const result=analyzeAlleleGroups({amplicons:input.references,primer:input.primer,control:input.control,sample:input.sample,ampliconConfirmed:true,settings:{...automaticSangerSettings(input),readEnd:1000,editStartRadius:10}},stage=>appendFileSync('../signal-audit-stages.log',Date.now()+' '+stage+'\n'),state=>audit=state);
 // Ground truth is first read after blind inference has completed.
 const truth=JSON.parse(readFileSync(dir+'truth.json','utf8')).truth,rows=truth.map((t:any)=>{const r=audit.refs.find((r:any)=>r.name===t.reference);return {...t,prepared:(r.reverse?reverseComplement(t.inputSequence):t.inputSequence).slice(r.offset)};}),trueColumns=rows.map((t:any)=>audit.positions.flatMap((p:number)=>emissionVector(t.prepared,p,audit.sampleLag,audit.theta))),trueWeights=rows.map((t:any)=>t.totalPercent/100),ideal=audit.positions.flatMap((p:number)=>'ACGT'.split('').map(b=>rows.reduce((sum:number,t:any)=>sum+(t.prepared[p]===b?t.totalPercent/100:0),0))),fit=fitNormalizedEmissionMixture(trueColumns,audit.observed),score=rawSignalScore(trueColumns,audit.rawObserved,trueWeights,audit.gains),fittedTruthScore=rawSignalScore(trueColumns,audit.rawObserved,fit.weights,audit.gains);
 const perGroup=rows.map((t:any,i:number)=>{const candidates=audit.models.flatMap((m:any,id:number)=>m.reference===audit.refs.findIndex((r:any)=>r.name===t.reference)&&m.sequence===t.prepared?[id]:[]);return {reference:t.reference,group:t.allele,description:t.description,percent:t.percent,candidatePresent:candidates.length>0,selected:candidates.filter((id:number)=>audit.ids.includes(id)),inferredPercent:candidates.reduce((sum:number,id:number)=>sum+(audit.ids.includes(id)?audit.weights[audit.ids.indexOf(id)]:0),0)/(t.share/100)*100,refittedTruthPercent:fit.weights[i]/(t.share/100)*100};});
 const report={result,blind:{r2:result.mixture?.r2,rmse:result.mixture?.rmse,groups:result.joint?.selectedModels},extraction:{maxAbsoluteError:Math.max(...ideal.map((v:number,i:number)=>Math.abs(v-audit.rawObserved[i]))),rmse:Math.sqrt(ideal.reduce((sum:number,v:number,i:number)=>sum+(v-audit.rawObserved[i])**2,0)/ideal.length),gains:audit.gains,lag:audit.sampleLag,theta:audit.theta},truthFixedScore:{r2:score.r2,rmse:score.rmse},truthRefitScore:{r2:fittedTruthScore.r2,rmse:fittedTruthScore.rmse},perGroup};writeFileSync(dir+'signal-audit.json',JSON.stringify(report,null,2));
},600000);
