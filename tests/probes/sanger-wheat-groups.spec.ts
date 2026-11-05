import '@angular/compiler';
import {it} from 'vitest';
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {exactBindingSites} from '../../src/app/workers/core/sanger-primer-designer';
import {SANGER_DEFAULTS,SangerCandidate,SangerDesignResult} from '../../src/app/models/sanger-design';
import {synthReferences,generateSangerSynth,randomSynthEvent,compareSangerSynth} from '../../src/app/utils/sanger-synthetic-benchmark';
import {reverseComplement} from '../../src/app/utils/biology.utils';
import {analyzeAutomaticSanger} from '../../src/app/workers/core/sanger-automatic';
const run=process.env.CASSAGO_WHEAT_GROUPS==='1'?it:it.skip;
run('random five groups on the existing public wheat TaPDS reference molecules',()=>{
 const config=JSON.parse(readFileSync('public/examples/sanger-wheat/config.json','utf8')),references=config.amplicons.map((a:any)=>({...a,cuts:[a.cut]})),p=config.primer,q=reverseComplement(references[0].sequence.slice(-20));
 const sites=(seq:string)=>references.map((r:any)=>{const h=exactBindingSites(seq,r);if(h.length!==1)throw Error('Primer binding not unique: '+r.name);return {...h[0],reference:r.name,cutOffsets:r.cuts.map((c:number)=>h[0].direction==='forward'?c-h[0].end:h[0].start-c),extension:r.sequence.slice(h[0].end)};});
 const primary=sites(p),opposite=sites(q),candidate={id:'wheat-benchmark',sequence:p,direction:'forward',sites:primary,partner:{sequence:q,sites:opposite,products:references.map((r:any)=>r.sequence.length)}}as SangerCandidate,result={references,candidates:[candidate],settings:{...SANGER_DEFAULTS,readEnd:1000,primerMin:20,primerMax:21}}as SangerDesignResult,refs=synthReferences(result,candidate);
 let state=20261108;const random=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state/4294967296;};for(const r of refs)r.alleles=Array.from({length:5},(_,i)=>({id:i+1,percent:[5,8,11,14,17][i],events:[randomSynthEvent(r,random)]}));
 const data=generateSangerSynth(result,candidate,refs),dir='../benchmarks/fresh-five-groups-20261108-wheat';mkdirSync(dir,{recursive:true});writeFileSync(dir+'/control.ab1',new Uint8Array(data.control));writeFileSync(dir+'/edit.ab1',new Uint8Array(data.edit));writeFileSync(dir+'/truth.json',JSON.stringify({recipe:data.recipe,truth:data.truth}));writeFileSync(dir+'/template.json',JSON.stringify(data.template));const started=Date.now(),analysis=analyzeAutomaticSanger({...data.input,groupRecovery:true});writeFileSync(dir+'/result.json',JSON.stringify({elapsedMs:Date.now()-started,result:analysis,comparison:compareSangerSynth(data,analysis)}));
},600000);
