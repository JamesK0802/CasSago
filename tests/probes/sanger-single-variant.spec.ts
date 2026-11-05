import '@angular/compiler';
import {it} from 'vitest';
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {profileFixture} from '../fixtures/sanger-profile-fixture';
import {generateSangerSynth} from '../../src/app/utils/sanger-synthetic-benchmark';
import {analyzeAlleleGroups} from '../../src/app/workers/core/sanger-allele-groups';
import {automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {benchmarkComparison,parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
const run=process.env.CASSAGO_SINGLE_VARIANT==='1'?it:it.skip;
run('checks every failed-recipe edit alone on one reference and one target',()=>{
 const source=JSON.parse(readFileSync('../benchmarks/variable-percent-20261104/5-both-recipe.json','utf8')).recipe,folder='../benchmarks/single-variant-20261104';mkdirSync(folder,{recursive:true});const rows=[];
 for(let h=0;h<3;h++)for(let target=0;target<2;target++)for(let group=0;group<5;group++){
  if(process.env.CASSAGO_POINT_ONLY==='1'&&(h!==1||target!==1||group!==3))continue;let audit:any;const f=profileFixture(3),r=structuredClone(source[h]);r.share=100;for(const p of r.profiles)p.alleles=p.target===target?[p.alleles[group]]:[];
  const result={...f.result,references:[f.result.references[h]]},candidate={...f.candidate,sites:[f.candidate.sites[h]],partner:{...f.candidate.partner!,sites:[f.candidate.partner!.sites[h]],products:[f.candidate.partner!.products[h]]}},d=generateSangerSynth(result,candidate,[r],20261135),references=d.input.references.map(a=>({...a,cut:a.cuts![target],cuts:[a.cuts![target]]})),analysis=analyzeAlleleGroups({amplicons:references,primer:d.input.primer!,control:d.input.control,sample:d.input.sample,ampliconConfirmed:true,settings:{...automaticSangerSettings({...d.input,references}),readEnd:1000,editStartRadius:10,maxIndel:10}},()=>{},state=>audit=state);
  const truth=d.truth.map(t=>({...t,pattern:t.pattern.split('/')[target],events:t.events.map(e=>({...e,target:0}))})),comparison=benchmarkComparison(parseBenchmarkTruth({truth}),references,analysis).filter(t=>t.scope==='group'&&t.group!=='WT'&&(t.truth??0)>0),e=r.profiles[target].alleles[0].events[0];rows.push({reference:h,target:target+1,group,kind:e.kind,percent:r.profiles[target].alleles[0].percent,sequenceMatched:comparison.some(c=>(c.estimated??0)>0),accurate2pp:comparison.every(c=>(c.estimated??0)>0&&Math.abs(c.delta!)<=2),comparison,rmse:analysis.mixture?.rmse});writeFileSync(folder+'/summary.json',JSON.stringify(rows,null,2));if(!rows.at(-1)!.sequenceMatched){writeFileSync(folder+'/failed-result.json',JSON.stringify({recipe:d.recipe,truth:d.truth,result:analysis}));writeFileSync(folder+'/failed-audit.json',JSON.stringify({refs:audit.refs,positions:audit.positions,models:audit.models,ids:audit.ids,weights:audit.weights,theta:audit.theta,lag:audit.sampleLag,gains:audit.gains,observed:audit.observed,iterations:audit.iterations}));}
 }
},600000);
