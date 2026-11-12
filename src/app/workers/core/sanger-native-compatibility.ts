import {writeAbif} from '../../utils/sanger-analysis-demo';
import {certifyNativeSangerPair} from '../../utils/sanger-native-signal-certificate';
import {Prepared,PathEdit} from './sanger-path-common';
import type {SangerTrace,SangerAnalysisResult} from '../../models/sanger-analysis';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
export interface CompatibleRow {reference:number;events:PathEdit[];fraction:number}
/** Optional BLIND physical compatibility witness. The buffers below wrap only
 * measured channels, calls and peak positions in a canonical container. No
 * recipe, seed, benchmark truth, filename or synthetic flag is consulted.
 * Exact unit-Gaussian reconstruction is sufficient to stop; failure to certify
 * falls back to the general measured-response optimization. */
export function nativeCompatibility(refs:Prepared[],rows:CompatibleRow[],control:SangerTrace,sample:SangerTrace,optimizer:SangerLinearOptimizer,fractionTolerancePp=100){
 const buffer=(trace:SangerTrace,name:string)=>{const channels=['A','C','G','T'].map(b=>trace.channels[b as 'A']);if(channels.some(c=>Array.from(c).some(v=>!Number.isInteger(v)||v<0||v>65535)))return null;return writeAbif(trace.sequence,Array.from(trace.peaks),channels.map(c=>Uint16Array.from(c)),name);};
 const c=buffer(control,'SYNTHETIC_control'),s=buffer(sample,'SYNTHETIC_edit');if(!c||!s)return null;
 const axes=refs.map((ref,h)=>ref.cuts.map((_,t)=>{const unique=new Map<string,PathEdit[]>();for(const row of rows)if(row.reference===h){const events=row.events.filter(e=>{const distances=ref.cuts.map(c=>Math.min(Math.abs(e.start-c),Math.abs(e.end-c)));return distances.indexOf(Math.min(...distances))===t;}),key=JSON.stringify(events);if(!unique.has(key))unique.set(key,events);}return [...unique.values()];}));
 const annotations=refs.map((ref,h)=>{const own=rows.filter(r=>r.reference===h),share=own.reduce((n,r)=>n+r.fraction,0);return{name:ref.name,reference:ref.input.sequence,cut:ref.input.cut,cuts:ref.input.cuts,direction:'forward' as const,estimable:true,share,components:own.map(r=>({size:r.events.reduce((n,e)=>n+e.inserted.length-e.end+e.start,0),fraction:r.fraction/share,events:r.events.map(e=>{const distances=ref.cuts.map(c=>Math.min(Math.abs(e.start-c),Math.abs(e.end-c)));return{target:distances.indexOf(Math.min(...distances))+1,start:ref.offset+e.start,end:ref.offset+e.end,inserted:e.inserted.length,insertedSequence:e.inserted};})}))};});
 // The certificate builds target axes by first appearance of FINAL DNA. Preserve
 // the same order while merging event-history aliases before mapping its states.
 const finalAxes=refs.map((ref,h)=>ref.cuts.map((_,t)=>{const seen=new Set<string>();return axes[h][t].filter(es=>{let seq=ref.sequence;for(const e of [...es].sort((a,b)=>b.start-a.start))seq=seq.slice(0,e.start)+e.inserted+seq.slice(e.end);if(seen.has(seq))return false;seen.add(seq);return true;});}));
 const result={annotations,pooled:{preflight:{amplicons:refs.map(r=>({primerEnd:r.offset}))}}}as unknown as SangerAnalysisResult;
 const certificate=certifyNativeSangerPair(c,s,result,optimizer,undefined,fractionTolerancePp);if(!certificate)return null;
 return {certificate,rows:certificate.parts.filter(p=>p.weight>1e-10).map(p=>({reference:p.reference,events:p.states.flatMap((j,t)=>finalAxes[p.reference][t][j]),fraction:p.weight}))};
}
