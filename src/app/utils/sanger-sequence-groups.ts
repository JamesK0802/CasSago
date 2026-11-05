import type {SangerJointEvent} from '../models/sanger-analysis';
import {canonicalSangerEvents} from './sanger-placement';
import {referenceChecksum} from './reference-checksum';
export type SequenceState='base-call'|'mixture'|'unphased'|'unresolved';
export interface SequenceChoice {events:SangerJointEvent[];weight:number;state:SequenceState}
export const SANGER_BASE_CALL_THRESHOLD=.99;
export const SANGER_SEQUENCE_CANDIDATE_LIMIT=64;
export function sangerEditedSequence(reference:string,events:SangerJointEvent[]){let sequence=reference;for(const e of [...events].sort((a,b)=>b.start-a.start||b.end-a.end))sequence=sequence.slice(0,e.start)+(e.insertedSequence||'N'.repeat(e.inserted))+sequence.slice(e.end);return sequence;}
export function sangerSequenceKey(reference:string,events:SangerJointEvent[]){return referenceChecksum(sangerEditedSequence(reference,events));}
/** Expand the fitted independent-base emission model, not sequenced molecules.
 * Multiple mixed positions are unphased candidates. Truncation retains its mass
 * in an explicit unresolved remainder; near-certain calls use a 1% tolerance. */
export function sangerSequenceChoices(reference:string,events:SangerJointEvent[],limit=SANGER_SEQUENCE_CANDIDATE_LIMIT):SequenceChoice[]{
 let rows:{events:SangerJointEvent[];weight:number}[]=[{events:events.map(e=>({...e,insertedSequence:e.insertedSequence||'N'.repeat(e.inserted),insertedBaseProbabilities:undefined})),weight:1}],mixed=0,missing=false;
 for(let e=0;e<events.length;e++)for(let j=0;j<events[e].inserted;j++){
  const raw=events[e].insertedBaseProbabilities?.[j],sum=raw?.reduce((a,b)=>a+b,0)||0,valid=raw?.length===4&&raw.every(v=>Number.isFinite(v)&&v>=0)&&sum>0;
  let options:{base:string;weight:number}[];
  if(valid){const p=raw!.map(v=>v/sum),max=Math.max(...p);options=max>=SANGER_BASE_CALL_THRESHOLD?[{base:'ACGT'[p.indexOf(max)],weight:1}]:p.flatMap((weight,b)=>weight>0?[{base:'ACGT'[b],weight}]:[]);if(options.length>1)mixed++;}
  else {const base=events[e].insertedSequence?.[j]||'N';options=[{base:/^[ACGT]$/.test(base)?base:'N',weight:1}];if(options[0].base==='N')missing=true;}
  rows=rows.flatMap(row=>options.map(o=>({weight:row.weight*o.weight,events:row.events.map((event,i)=>i!==e?event:{...event,insertedSequence:event.insertedSequence!.slice(0,j)+o.base+event.insertedSequence!.slice(j+1)})}))).sort((a,b)=>b.weight-a.weight||JSON.stringify(a.events).localeCompare(JSON.stringify(b.events))).slice(0,Math.max(1,limit));
 }
 const state:SequenceState=missing?'unresolved':mixed>1?'unphased':mixed?'mixture':'base-call';
 const out:SequenceChoice[]=rows.map(row=>({...row,events:canonicalSangerEvents(reference,row.events),state})),remaining=1-rows.reduce((sum,row)=>sum+row.weight,0);
 if(remaining>1e-12)out.push({events:events.map(e=>({...e,insertedSequence:e.inserted?'N'.repeat(e.inserted):'',insertedBaseProbabilities:undefined})),weight:remaining,state:'unresolved'});
 return out;
}
