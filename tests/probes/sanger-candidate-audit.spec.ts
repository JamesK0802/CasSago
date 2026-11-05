import '@angular/compiler';
import {it} from 'vitest';
import {readFileSync,writeFileSync} from 'node:fs';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
import {analyzeAlleleGroups} from '../../src/app/workers/core/sanger-allele-groups';
import {automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {materialize} from '../../src/app/workers/core/sanger-path-common';
const run=process.env.CASSAGO_CANDIDATE_AUDIT==='1'?it:it.skip;
run('locates missed target sequence groups in the blinded candidate dictionary after inference',()=>{
 const root='../benchmarks/more-groups-20261103',n=5,t=JSON.parse(readFileSync(`${root}/${n}-template.json`,'utf8')),references=t.amplicons.map((a:any)=>({name:a.name,sequence:a.sequence,cut:a.cuts[0],cuts:a.cuts})),base={references,primer:t.primer,control:parseSangerTrace(Uint8Array.from(readFileSync(`${root}/${n}-control.ab1`)).buffer),sample:parseSangerTrace(Uint8Array.from(readFileSync(`${root}/${n}-edit.ab1`)).buffer),maxDeletion:10,maxInsertion:10};let audit:any;
 const result=analyzeAlleleGroups({amplicons:references,primer:t.primer,control:base.control,sample:base.sample,ampliconConfirmed:true,settings:{...automaticSangerSettings(base),readEnd:1000,editStartRadius:10,maxIndel:10}},()=>{},state=>audit=state);
 // Recipes are read strictly after the blind engine has returned its dictionary.
 const data=JSON.parse(readFileSync(`${root}/${n}-recipe.json`,'utf8')),rows=[];
 for(let h=0;h<data.recipe.length;h++){const ref=data.recipe[h],r=audit.refs[h];for(const profile of ref.profiles)for(const allele of profile.alleles){const e=allele.events[0],start=(e.kind==='substitution'?e.selectionStart:ref.cuts[e.target]+e.offset)-r.offset,end=start+(e.kind==='insertion'?0:e.kind==='deletion'?e.length:e.before.length),inserted=e.kind==='deletion'?'':e.sequence;
 const wanted=materialize(r.sequence,[{start,end,inserted}]),cuts=r.cuts,near=(event:any)=>cuts.map((c:number)=>Math.min(Math.abs(c-event.start),Math.abs(c-event.end))).reduce((best:number,d:number,i:number,ds:number[])=>d<ds[best]?i:best,0),matches:number[]=[];
 for(let id=0;id<audit.models.length;id++){const m=audit.models[id];if(m.reference!==h)continue;const local=m.events.filter((event:any)=>near(event)===profile.target);if(materialize(r.sequence,local)===wanted)matches.push(id);}
 const selected=matches.filter(id=>audit.ids.includes(id)),fraction=selected.reduce((s,id)=>s+audit.weights[audit.ids.indexOf(id)],0);rows.push({reference:ref.name,target:profile.target+1,group:allele.id,kind:e.kind,expected:allele.percent,candidates:matches.length,selected:selected.length,estimatedSamplePercent:fraction*100});
 }}
 writeFileSync(`${root}/5-candidate-audit.json`,JSON.stringify({candidates:result.joint?.candidateModels,rows},null,2));
},600000);
