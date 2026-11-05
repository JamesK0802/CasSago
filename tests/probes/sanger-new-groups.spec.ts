import '@angular/compiler';
import {it} from 'vitest';
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
import {parseAutomaticSangerTemplate} from '../../src/app/utils/sanger-automatic-template';
import {automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {analyzeAlleleGroups as analyzeSignalPaths} from '../../src/app/workers/core/sanger-allele-groups';
import {compareSangerSynth} from '../../src/app/utils/sanger-synthetic-benchmark';
const run=process.env.CASSAGO_NEW_GROUPS==='1'?it:it.skip;
import {sangerDesignExample} from '../../src/app/utils/sanger-design.utils';
import {exactBindingSites} from '../../src/app/workers/core/sanger-primer-designer';
import {SANGER_DEFAULTS,SangerDesignResult,SangerCandidate} from '../../src/app/models/sanger-design';
import {synthReferences,generateSangerSynth,randomSynthEvent} from '../../src/app/utils/sanger-synthetic-benchmark';
import {reverseComplement} from '../../src/app/utils/biology.utils';

function setup(reverse=false){const references=sangerDesignExample(),p='CCGGCGGGTACGACCATAAC',q='CAGCGACTCCTTGAGTGACG',sites=(seq:string)=>references.map(r=>{const h=exactBindingSites(seq,r)[0];return {...h,reference:r.name,cutOffsets:r.cuts.map(c=>h.direction==='forward'?c-h.end:h.start-c),extension:h.direction==='forward'?r.sequence.slice(h.end):''};});const primary=sites(reverse?q:p),other=sites(reverse?p:q),candidate={id:'test',sequence:reverse?q:p,direction:reverse?'reverse':'forward',sites:primary,partner:{sequence:reverse?p:q,sites:other,products:references.map((_,i)=>Math.max(primary[i].end,other[i].end)-Math.min(primary[i].start,other[i].start))}} as SangerCandidate,result={references,candidates:[candidate],settings:{...SANGER_DEFAULTS,readEnd:1000,primerMin:20,primerMax:20}} as SangerDesignResult;return{result,candidate,refs:synthReferences(result,candidate)};}
for(const seed of [20261109,20261110])run('new random five-group dataset '+seed,()=>{
 let state=seed;const random=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state/4294967296;};
 const {result,candidate,refs}=setup();for(const ref of refs)ref.alleles=Array.from({length:5},(_,i)=>({id:i+1,percent:seed===20261109?10:[5,8,11,14,17][i],events:[randomSynthEvent(ref,random)]}));
 const data=generateSangerSynth(result,candidate,refs),settings={...automaticSangerSettings(data.input),jointEngine:'signal-local' as const,editStartRadius:10,readEnd:1000};
 const started=Date.now(),analysis=analyzeSignalPaths({amplicons:data.input.references,primer:data.input.primer!,control:data.input.control,sample:data.input.sample,ampliconConfirmed:true,settings,synthetic:true});
 const dir='../benchmarks/fresh-five-groups-'+seed;mkdirSync(dir,{recursive:true});writeFileSync(dir+'/control.ab1',new Uint8Array(data.control));writeFileSync(dir+'/edit.ab1',new Uint8Array(data.edit));writeFileSync(dir+'/truth.json',JSON.stringify({recipe:data.recipe,truth:data.truth}));writeFileSync(dir+'/template.json',JSON.stringify(data.template));writeFileSync(dir+'/result.json',JSON.stringify({elapsedMs:Date.now()-started,result:analysis}));
},600000);
