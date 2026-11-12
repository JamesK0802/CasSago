import {originBaseMass} from './sanger-origin-base-support';
import type {MarginalInternalGroup} from './sanger-marginal-engine';
import type {Prepared} from './sanger-path-common';
import type {SangerJointEvent,SangerMarginalGroup} from '../../models/sanger-analysis';

/** One modal proposal from the fitted allocation, without another LP or GT.
 * Uses the same origin-specific Hall residual capacities as lazy seed proposals.
 * The full length-family mass is not assigned to this proposed haplotype. */
export function fittedMarginalRepresentative(g:MarginalInternalGroup,point:number[],ref:Prepared,share:number,event:(start:number,end:number,inserted:string)=>SangerJointEvent):SangerMarginalGroup['representative']{
 const DNA='ACGT';
 let best:SangerMarginalGroup['representative'];
 for(const e of g.edges){
  let support=Math.max(0,point[e.mass]);if(support<=1e-7)continue;
  const word=new Map<number,string>();
  for(const[q]of g.bases){
   const insertion=q>=e.start+e.incoming&&q<e.start+e.incoming+Math.max(0,e.size),free=insertion||Math.abs(q-e.incoming-ref.cuts[e.target])<=45,at=q<e.start+e.incoming?q-e.incoming:q-e.outgoing;
   const values=Array.from({length:4},(_,b)=>originBaseMass(g,e,q,b,ref.sequence,point)),base=free?values.indexOf(Math.max(...values)):DNA.indexOf(ref.sequence[at]);
   if(base<0){support=0;break;}support=Math.min(support,Math.max(0,values[base]));if(free)word.set(q,DNA[base]);
  }
  if(support<=1e-7||best&&support/share<=best.witnessFraction+1e-10)continue;
  const changes:SangerJointEvent[]=[],inserted=Array.from({length:Math.max(0,e.size)},(_,i)=>word.get(e.start+e.incoming+i)||'').join('');
  if(inserted.length!==Math.max(0,e.size))continue;
  if(e.size)changes.push(event(e.start,e.start+Math.max(0,-e.size),inserted));
  for(const[q,base]of word){if(q>=e.start+e.incoming&&q<e.start+e.incoming+Math.max(0,e.size))continue;const at=q<e.start+e.incoming?q-e.incoming:q-e.outgoing;if(base!==ref.sequence[at])changes.push(event(at,at+1,base));}
  changes.sort((a,b)=>a.start-b.start||a.end-b.end);
  best={events:changes,witnessFraction:support/share,source:'fitted-marginal-proposal'};
 }
 return best;
}
