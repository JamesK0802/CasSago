import '@angular/compiler';
import {it} from 'vitest';
import {mkdirSync,readFileSync,writeFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {profileFixture} from '../fixtures/sanger-profile-fixture';
import {generateSangerSynth,SynthReference} from '../../src/app/utils/sanger-synthetic-benchmark';
import {analyzeAlleleGroups} from '../../src/app/workers/core/sanger-allele-groups';
import {automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {benchmarkComparison,parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
import {projectTargetAnnotation} from '../../src/app/utils/sanger-target-profile';
import {sangerAnnotationViews} from '../../src/app/utils/sanger-annotation';
const run=process.env.CASSAGO_VARIABLE_PERCENT==='1'?it:it.skip;
function randomGenerator(seed:number){return()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/4294967296;};}
/** Fixed total editing keeps the paired experiment comparable. Values are
 * distinct across profiles; within each profile the gap is at least 0.5 pp. */
function variablePercentages(count:number,random:()=>number,used:Set<number>){
 for(let attempt=0;attempt<100000;attempt++){const raw=Array.from({length:count},()=>.25+1.5*random()),total=raw.reduce((s,x)=>s+x,0),values=raw.map(v=>Math.round(5000*v/total));values[count-1]+=5000-values.reduce((s,v)=>s+v,0);const sorted=[...values].sort((a,b)=>a-b);if(values.some(v=>v<150||used.has(v))||sorted.some((v,i)=>i>0&&v-sorted[i-1]<50))continue;values.forEach(v=>used.add(v));return values.map(v=>v/100);}
 throw Error('Could not allocate distinct benchmark percentages.');
}
run('compares identical edits with unequal percentages and isolates T1/T2 cross-target effects',()=>{
 const folder=process.env.CASSAGO_PERCENT_FOLDER||'../benchmarks/variable-percent-20261104';mkdirSync(folder,{recursive:true});const summary:any[]=[];
 const cases=[{count:3,mode:'both'},{count:5,mode:'both'},{count:8,mode:'both'},{count:5,mode:'T1-WT'},{count:5,mode:'T2-WT'}];
 for(const test of cases){const {count,mode}=test,key=count+'-'+mode;if(process.env.CASSAGO_PERCENT_CASES&&!process.env.CASSAGO_PERCENT_CASES.split(',').includes(key))continue;const seed=20261200+count,random=randomGenerator(seed),used=new Set<number>(),f=profileFixture(3),refs=JSON.parse(readFileSync(`../benchmarks/more-groups-20261103/${count}-recipe.json`,'utf8')).recipe as SynthReference[];
  for(const ref of refs)for(const profile of ref.profiles!){const values=variablePercentages(count,random,used);profile.alleles.forEach((a,i)=>a.percent=values[i]);}
  const fullRecipe=structuredClone(refs);for(const ref of refs)for(const p of ref.profiles!)if(mode==='T1-WT'&&p.target===0||mode==='T2-WT'&&p.target===1)p.alleles=[];
  // The same finite pairing seed and unchanged edit recipes isolate proportion effects.
  const d=generateSangerSynth(f.result,f.candidate,refs,20261130+count),input={amplicons:d.input.references,primer:d.input.primer!,control:d.input.control,sample:d.input.sample,ampliconConfirmed:true,settings:{...automaticSangerSettings(d.input),readEnd:1000,editStartRadius:10,maxIndel:10}};
  writeFileSync(`${folder}/${key}-recipe.json`,JSON.stringify({recipe:d.recipe,truth:d.truth,mixing:d.mixing,fullRecipe},null,2));writeFileSync(`${folder}/${key}-control.ab1`,new Uint8Array(d.control));writeFileSync(`${folder}/${key}-edit.ab1`,new Uint8Array(d.edit));writeFileSync(`${folder}/${key}-template.json`,JSON.stringify(d.template));
  const start=Date.now(),result=analyzeAlleleGroups(input),elapsedMs=Date.now()-start;
  // Ground truth enters evaluation only after inference has completed.
  const comparison=benchmarkComparison(parseBenchmarkTruth({recipe:d.recipe,truth:d.truth}),d.input.references,result),targets=[0,1].map(target=>{const rows=comparison.filter(r=>r.target===target&&r.scope==='group'&&(r.truth??0)>0&&r.group!=='WT'),extra=comparison.filter(r=>r.target===target&&r.scope==='group'&&r.truth===0&&(r.estimated??0)>=1);return{target:target+1,total:rows.length,sequenceMatches:rows.filter(r=>(r.estimated??0)>0).length,accurate2pp:rows.filter(r=>(r.estimated??0)>0&&Math.abs(r.delta!)<=2).length,falsePositiveGroupsAtLeast1Percent:extra.length,falsePositiveSamplePercent:extra.reduce((s,r)=>s+(r.estimated||0)*(result.annotations?.find(a=>a.name===r.reference)?.share||0),0)};});
  const placementAudit=refs.flatMap((ref,h)=>ref.profiles!.flatMap(profile=>{const a=projectTargetAnnotation(result.annotations![h],profile.target),rows=sangerAnnotationViews([a],0,60)[0].rows;return profile.alleles.map(allele=>{const e=allele.events[0],start=e.kind==='substitution'?e.selectionStart!:ref.cuts[e.target]+e.offset,removed=e.kind==='insertion'?0:e.kind==='deletion'?e.length:e.before!.length,inserted=e.kind==='deletion'?'':e.sequence,candidates=rows.flatMap(row=>{const events=row.representative?.events||[];if(events.length!==1)return[];if(row.fraction<.01)return[];const actual=events[0];if(actual.end-actual.start!==removed||actual.insertedSequence!==inserted)return[];return [{fractionPercent:row.fraction*100,start:actual.start,positionDeltaBp:actual.start-start,equivalentSequence:ref.sequence.slice(0,actual.start)+inserted+ref.sequence.slice(actual.end)===ref.sequence.slice(0,start)+inserted+ref.sequence.slice(start+removed)}];}).sort((a,b)=>Math.abs(a.positionDeltaBp)-Math.abs(b.positionDeltaBp));return {reference:ref.name,target:profile.target+1,group:allele.id,kind:e.kind,expectedPercent:allele.percent,expectedStart:start,removed,inserted,locationCandidates:candidates.slice(0,3)};});}));
  const row={count,mode,seed,elapsedMs,molecules:d.truth.length,targets,rmse:result.mixture?.rmse,passesChecks:result.mixture?.passesChecks,iterations:result.joint?.iterations,editSha256:createHash('sha256').update(new Uint8Array(d.edit)).digest('hex'),positionOnlyCandidates:placementAudit.filter(r=>r.locationCandidates.length&&!r.locationCandidates.some(c=>c.equivalentSequence)&&r.locationCandidates[0].positionDeltaBp!==0).map(r=>({reference:r.reference,target:r.target,group:r.group,kind:r.kind,delta:r.locationCandidates[0].positionDeltaBp}))};summary.push(row);writeFileSync(folder+'/summary.json',JSON.stringify(summary,null,2));writeFileSync(`${folder}/${key}-result.json`,JSON.stringify(result));writeFileSync(`${folder}/${key}-comparison.json`,JSON.stringify(comparison));writeFileSync(`${folder}/${key}-placements.json`,JSON.stringify({note:'Same edit payload/length at another location is a diagnostic candidate, not a ground-truth assignment. Repeats can also produce equivalent breakpoints.',rows:placementAudit},null,2));console.log(JSON.stringify(row));
 }
},1800000);
