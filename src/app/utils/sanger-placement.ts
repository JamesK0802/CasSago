import type {SangerJointEvent} from '../models/sanger-analysis';

/** One canonical representation for exactly the same final DNA. This does not
 * infer the physical repair history or merge different edited sequences. */
export function leftEquivalentDeletionStart(reference:string,start:number,length:number):number {
 while(start>0&&/^[ACGT]$/.test(reference[start-1])&&reference[start-1]===reference[start+length-1])start--;
 return start;
}
export function canonicalSangerEvents(reference:string,events:SangerJointEvent[]):SangerJointEvent[]{
 if(events.length!==1)return events.map(e=>({...e}));
 const e={...events[0]};
 if(e.end>e.start&&e.inserted===0){const length=e.end-e.start;e.start=leftEquivalentDeletionStart(reference,e.start,length);e.end=e.start+length;}
 else if(e.end===e.start&&e.inserted>0&&e.insertedSequence?.length===e.inserted&&/^[ACGT]+$/.test(e.insertedSequence)){
  let sequence=e.insertedSequence;
  while(e.start>0&&reference[e.start-1]===sequence.at(-1)){sequence=sequence.at(-1)!+sequence.slice(0,-1);e.start--;e.end--;}
  e.insertedSequence=sequence;
 }
 return [e];
}

/** Push a fitted independent-base insertion distribution through DNA-equivalent
 * left normalization exactly. A shift requires successive inserted suffix bases
 * to match the WT bases on the left. No guessed inserted sequence is used. */
export function canonicalSangerPlacements(reference:string,events:SangerJointEvent[]):{events:SangerJointEvent[];weight:number}[]{
 const e=events[0],pi=e?.insertedBaseProbabilities,bases='ACGT';
 if(events.length!==1||e.start!==e.end||!e.inserted||pi?.length!==e.inserted||pi.some(v=>v.length!==4||v.some(x=>!Number.isFinite(x)||x<0)||Math.abs(v.reduce((s,x)=>s+x,0)-1)>1e-6))return [{events:canonicalSangerEvents(reference,events),weight:1}];
 const out:{events:SangerJointEvent[];weight:number}[]=[];let remaining=1,shift=0;
 while(remaining>1e-15&&shift<Math.min(e.inserted,e.start)){
  const b=bases.indexOf(reference[e.start-shift-1]);if(b<0)break;
  const match=pi[e.inserted-shift-1][b],mass=remaining*(1-match);
  if(mass>1e-15)out.push({events:[{...e,start:e.start-shift,end:e.start-shift,insertedSequence:'N'.repeat(e.inserted),insertedBaseProbabilities:undefined}],weight:mass});
  remaining*=match;shift++;
 }
 if(remaining>1e-15){
  let start=e.start-shift;
  if(shift===e.inserted){const sequence=reference.slice(e.start-e.inserted,e.start),normalized=canonicalSangerEvents(reference,[{...e,start,end:start,insertedSequence:sequence,insertedBaseProbabilities:undefined}])[0];start=normalized.start;}
  out.push({events:[{...e,start,end:start,insertedSequence:'N'.repeat(e.inserted),insertedBaseProbabilities:undefined}],weight:remaining});
 }
 return out;
}
