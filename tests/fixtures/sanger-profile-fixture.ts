import {sangerDesignExample} from '../../src/app/utils/sanger-design.utils';
import {exactBindingSites} from '../../src/app/workers/core/sanger-primer-designer';
import {SANGER_DEFAULTS,SangerCandidate,SangerDesignResult} from '../../src/app/models/sanger-design';
import {synthReferences} from '../../src/app/utils/sanger-synthetic-benchmark';
import {initSynthProfiles} from '../../src/app/utils/sanger-synth-profiles';
export function profileFixture(referenceCount=1){
 const references=sangerDesignExample(2).slice(0,referenceCount),p='CCGGCGGGTACGACCATAAC',q='CAGCGACTCCTTGAGTGACG',sites=(seq:string)=>references.map(r=>{const h=exactBindingSites(seq,r)[0];return {...h,reference:r.name,cutOffsets:r.cuts.map(c=>h.direction==='forward'?c-h.end:h.start-c),extension:h.direction==='forward'?r.sequence.slice(h.end):''};}),a=sites(p),b=sites(q),candidate={id:'profiles',sequence:p,direction:'forward',sites:a,partner:{sequence:q,sites:b,products:references.map((_,i)=>Math.max(a[i].end,b[i].end)-Math.min(a[i].start,b[i].start))}}as SangerCandidate,result={references,candidates:[candidate],settings:{...SANGER_DEFAULTS,readEnd:1000,primerMin:20,primerMax:20}}as SangerDesignResult,refs=synthReferences(result,candidate),r=refs[0];refs.forEach(initSynthProfiles);r.profiles![0].alleles=[{id:1,percent:25,events:[{kind:'deletion',target:0,offset:0,length:3,sequence:''}]}];r.profiles![1].alleles=[{id:2,percent:40,events:[{kind:'insertion',target:1,offset:0,length:2,sequence:'GC'}]}];return {result,candidate,refs};
}

/** Rebuild the design wrapper around actual submitted products, without
 * changing their DNA, cut positions, or exported product names. */
export function submittedProfileFixture(template:{primer:string;oppositePrimer:string;amplicons:{name:string;sourceReference?:string;sequence:string;cuts:number[];start:number}[]}){
 const references=template.amplicons.map(r=>({name:r.sourceReference||r.name,sequence:'N'.repeat(r.start)+r.sequence,cuts:r.cuts.map(c=>c+r.start)})),p=template.primer,q=template.oppositePrimer;
 const sites=(primer:string)=>references.map(r=>{const hits=exactBindingSites(primer,r);if(hits.length!==1)throw Error(r.name+': expected exactly one primer site in the submitted product');const h=hits[0];return {...h,reference:r.name,cutOffsets:r.cuts.map(c=>h.direction==='forward'?c-h.end:h.start-c),extension:h.direction==='forward'?r.sequence.slice(h.end):''};});
 const a=sites(p),b=sites(q),candidate={id:'submitted-profiles',sequence:p,direction:'forward',sites:a,partner:{sequence:q,sites:b,products:references.map((_,i)=>Math.max(a[i].end,b[i].end)-Math.min(a[i].start,b[i].start))}}as SangerCandidate,result={references,candidates:[candidate],settings:{...SANGER_DEFAULTS,readEnd:1000,primerMin:p.length,primerMax:p.length}}as SangerDesignResult,refs=synthReferences(result,candidate);
 refs.forEach(initSynthProfiles);if(refs.some((r,i)=>r.name!==template.amplicons[i].name||r.sequence!==template.amplicons[i].sequence||r.cuts.some((c,t)=>c!==template.amplicons[i].cuts[t])))throw Error('Submitted reference panel was altered.');
 return {result,candidate,refs};
}
