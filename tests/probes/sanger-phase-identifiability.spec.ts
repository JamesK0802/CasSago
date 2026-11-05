import {it} from 'vitest';
import {readFileSync,writeFileSync} from 'node:fs';
import {simplexFit} from '../../src/app/workers/core/sanger-signal';
import {materializeGroup} from '../../src/app/utils/sanger-truth-groups';
const run=process.env.CASSAGO_PHASE_AUDIT==='1'?it:it.skip;
run('measure whether the original wrong groups can explain the ideal noiseless mixture',()=>{
 const dir='../benchmarks/fresh-five-groups-20261103/',truth=JSON.parse(readFileSync(dir+'truth.json','utf8')),old=JSON.parse(readFileSync(dir+'browser-result.json','utf8')).result,template=JSON.parse(readFileSync(dir+'template.json','utf8')),primer=template.primer.length;
 const models=old.joint.groups.map((g:any)=>({reference:g.reference,sequence:materializeGroup(template.amplicons.find((r:any)=>r.name===g.reference).sequence,g.events).slice(primer),events:g.events})),parts=truth.truth.map((t:any)=>({sequence:t.sequence.slice(primer),weight:t.totalPercent/100})),end=Math.min(...parts.map((p:any)=>p.sequence.length),...models.map((m:any)=>m.sequence.length)),columns=models.map((m:any)=>Array.from({length:end},(_,p)=>'ACGT'.split('').map(b=>Number(m.sequence[p]===b))).flat()),ideal=Array.from({length:end},(_,p)=>'ACGT'.split('').map(b=>parts.reduce((sum:number,t:any)=>sum+(t.sequence[p]===b?t.weight:0),0))).flat(),fit=simplexFit(columns,ideal);
 writeFileSync(dir+'phase-identifiability.json',JSON.stringify({rmse:fit.rmse,r2:fit.r2,converged:fit.converged,optimalityGap:fit.optimalityGap,models:models.map((m:any,i:number)=>({...m,sequence:undefined,fraction:fit.weights[i]}))},null,2));
},30000);
