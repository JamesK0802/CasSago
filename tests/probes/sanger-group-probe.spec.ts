import '@angular/compiler';
import {it} from 'vitest';
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
import {parseAutomaticSangerTemplate} from '../../src/app/utils/sanger-automatic-template';
import {automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {analyzeAlleleGroups as analyzeSignalPaths} from '../../src/app/workers/core/sanger-allele-groups';
import {compareSangerSynth} from '../../src/app/utils/sanger-synthetic-benchmark';
const run=process.env.CASSAGO_GROUP_PROBE==='1'?it:it.skip;
run.each(['signal-local']as const)('full sequence engine probe %s',kind=>{const dir='../benchmarks/random-five-groups-20261103/',templateRaw=JSON.parse(readFileSync(dir+'analysis-template.json','utf8')),template=parseAutomaticSangerTemplate(templateRaw),buffer=(name:string)=>{const b=readFileSync(dir+name);return b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength);},control=parseSangerTrace(buffer('SYNTHETIC_control.ab1'),'control.ab1'),sample=parseSangerTrace(buffer('SYNTHETIC_edit.ab1'),'edit.ab1'),references=template.references.map(r=>({name:r.name,sequence:r.sequence,cut:Number(r.targets.split(',')[0]),cuts:r.targets.split(',').map(Number)})),input={references,primer:template.primer,control,sample,maxDeletion:10,maxInsertion:10,synthetic:true},settings={...automaticSangerSettings(input),jointEngine:kind,editStartRadius:10,readEnd:1000},started=Date.now();let old='';try{const result=analyzeSignalPaths({amplicons:references,primer:template.primer,control,sample,ampliconConfirmed:true,settings,synthetic:true},stage=>{if(stage!==old){console.log(kind,stage);old=stage;}}),truth=JSON.parse(readFileSync(dir+'truth.json','utf8')),comparison=compareSangerSynth({template:templateRaw,truth:truth.truth},result);writeFileSync(dir+'prototype-'+kind+'.json',JSON.stringify({elapsedMs:Date.now()-started,result,comparison}));console.log('PROBE',kind,'r2',result.mixture?.r2,'passes',result.mixture?.passesChecks,comparison.filter(r=>r.scope==='group'&&(r.truth||0)>0).map(r=>({group:r.group,truth:r.truth,estimated:r.estimated,delta:r.delta,placement:r.placement})));}catch(e){writeFileSync(dir+'prototype-'+kind+'-error.json',JSON.stringify({error:String(e),elapsedMs:Date.now()-started}));throw e;}},600000);

import {sangerDesignExample} from '../../src/app/utils/sanger-design.utils';
import {exactBindingSites} from '../../src/app/workers/core/sanger-primer-designer';
import {SANGER_DEFAULTS,SangerDesignResult,SangerCandidate} from '../../src/app/models/sanger-design';
import {synthReferences,generateSangerSynth,randomSynthEvent} from '../../src/app/utils/sanger-synthetic-benchmark';
import {reverseComplement} from '../../src/app/utils/biology.utils';

function setup(reverse=false){const references=sangerDesignExample(),p='CCGGCGGGTACGACCATAAC',q='CAGCGACTCCTTGAGTGACG',sites=(seq:string)=>references.map(r=>{const h=exactBindingSites(seq,r)[0];return {...h,reference:r.name,cutOffsets:r.cuts.map(c=>h.direction==='forward'?c-h.end:h.start-c),extension:h.direction==='forward'?r.sequence.slice(h.end):''};});const primary=sites(reverse?q:p),other=sites(reverse?p:q),candidate={id:'test',sequence:reverse?q:p,direction:reverse?'reverse':'forward',sites:primary,partner:{sequence:reverse?p:q,sites:other,products:references.map((_,i)=>Math.max(primary[i].end,other[i].end)-Math.min(primary[i].start,other[i].start))}} as SangerCandidate,result={references,candidates:[candidate],settings:{...SANGER_DEFAULTS,readEnd:1000,primerMin:20,primerMax:20}} as SangerDesignResult;return{result,candidate,refs:synthReferences(result,candidate)};}
for(const seed of [20261103,20261104,20261105])run('new random five-group dataset '+seed,()=>{
 let state=seed;const random=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state/4294967296;};
 const {result,candidate,refs}=setup();for(const ref of refs)ref.alleles=Array.from({length:5},(_,i)=>({id:i+1,percent:[5,8,11,14,17][i],events:[randomSynthEvent(ref,random)]}));
 const data=generateSangerSynth(result,candidate,refs),settings={...automaticSangerSettings(data.input),jointEngine:'signal-local' as const,editStartRadius:10,readEnd:1000};
 const started=Date.now(),analysis=analyzeSignalPaths({amplicons:data.input.references,primer:data.input.primer!,control:data.input.control,sample:data.input.sample,ampliconConfirmed:true,settings,synthetic:true});
 const dir='../benchmarks/fresh-five-groups-'+seed;mkdirSync(dir,{recursive:true});writeFileSync(dir+'/control.ab1',new Uint8Array(data.control));writeFileSync(dir+'/edit.ab1',new Uint8Array(data.edit));writeFileSync(dir+'/truth.json',JSON.stringify({recipe:data.recipe,truth:data.truth}));writeFileSync(dir+'/template.json',JSON.stringify(data.template));writeFileSync(dir+'/result.json',JSON.stringify({elapsedMs:Date.now()-started,result:analysis}));
},600000);

for(const seed of [20261106,20261107])run('new random five-group dataset '+seed,()=>{
 let state=seed;const random=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state/4294967296;};
 const {result,candidate,refs}=setup();for(const ref of refs){ref.sequence=ref.sequence.slice(0,ref.primerLength)+Array.from({length:ref.sequence.length-ref.primerLength-ref.partnerLength},()=> 'ACGT'[Math.floor(random()*4)]).join('')+ref.sequence.slice(-ref.partnerLength);const source=result.references.find(r=>r.name===ref.sourceReference)!;const site=candidate.sites.find(r=>r.reference===ref.sourceReference)!;source.sequence=source.sequence.slice(0,site.start)+ref.sequence+source.sequence.slice(site.start+ref.sequence.length);}
 for(const ref of refs)ref.alleles=Array.from({length:5},(_,i)=>({id:i+1,percent:[5,8,11,14,17][i],events:[randomSynthEvent(ref,random)]}));
 const data=generateSangerSynth(result,candidate,refs),settings={...automaticSangerSettings(data.input),jointEngine:'signal-local' as const,editStartRadius:10,readEnd:1000};
 const started=Date.now(),analysis=analyzeSignalPaths({amplicons:data.input.references,primer:data.input.primer!,control:data.input.control,sample:data.input.sample,ampliconConfirmed:true,settings,synthetic:true});
 const dir='../benchmarks/fresh-five-groups-'+seed;mkdirSync(dir,{recursive:true});writeFileSync(dir+'/control.ab1',new Uint8Array(data.control));writeFileSync(dir+'/edit.ab1',new Uint8Array(data.edit));writeFileSync(dir+'/truth.json',JSON.stringify({recipe:data.recipe,truth:data.truth}));writeFileSync(dir+'/template.json',JSON.stringify(data.template));writeFileSync(dir+'/result.json',JSON.stringify({elapsedMs:Date.now()-started,result:analysis}));
},600000);
