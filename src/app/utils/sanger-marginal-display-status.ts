import type {SangerAnnotation,SangerMarginalGroup} from '../models/sanger-analysis';
import type {SangerAnnotationRow} from './sanger-annotation';
import {sangerEditedSequence,sangerSequenceKey} from './sanger-sequence-groups';
import {targetEvents} from './sanger-target-profile';

export type MarginalComparisonStatus='correct'|'covered'|'mismatch';
/** A compact fitted proposal is not a uniqueness certificate. Colour only
 * follows the retained envelopes, never a GT match or the selected estimate. */
export function marginalSequenceResolved(g:SangerMarginalGroup){
 if(!g.representative||!g.bases.length||g.status!=='bounded')return false;
 const possibleDifferent=g.bases.reduce((n,p)=>n+Math.max(0,1-Math.max(...p.lower)),0);
 return possibleDifferent<=.002&&(!g.unmodifiedRange||g.unmodifiedRange[1]-g.unmodifiedRange[0]<=.002)&&(!g.substitutionRange||g.substitutionRange[1]-g.substitutionRange[0]<=.002);
}

/** Display-only GT coverage. A family can contain several GT alleles: compare
 * their weighted coordinate proportions, not each allele as 100% at every base.
 * This is not a phased-haplotype certificate or a benchmark/quantification pass. */
export function marginalTruthStatus(truth:SangerAnnotation,g:SangerMarginalGroup,row:SangerAnnotationRow):{status:MarginalComparisonStatus;reason:string}{
 const target=g.target-1,components=truth.components.filter(c=>(c.pathSizes?.[target]??targetEvents(c.events||[],target).reduce((n,e)=>n+e.inserted-(e.end-e.start),0))===g.size),mass=components.reduce((n,c)=>n+c.fraction,0);
 if(!components.length||mass<=0||Math.abs(g.fraction-mass)>.02+1e-10)return{status:'mismatch',reason:'No matching length family within 2 percentage points.'};
 if(!g.bases.length||components.some(c=>c.unresolved))return{status:'mismatch',reason:'Base-range coverage is unavailable for this GT family.'};
 const sequences=components.map(c=>({fraction:c.fraction,sequence:sangerEditedSequence(truth.reference,c.events||[])}));
 for(const p of g.bases){const values=[0,0,0,0];for(const c of sequences){const q=truth.direction==='reverse'?p.position+c.sequence.length-truth.reference.length:p.position,b='ACGT'.indexOf(c.sequence[q]);if(b<0)return{status:'mismatch',reason:'GT bases are unavailable at a reported read position.'};values[b]+=c.fraction/mass;}if(values.some((v,b)=>v<p.lower[b]-.002||v>p.upper[b]+.002))return{status:'mismatch',reason:'GT family base proportions fall outside the reported A/C/G/T bounds.'};}
 const exact=!!g.representative&&!!row.sequenceKey&&row.sequenceKey===sangerSequenceKey(truth.reference,g.representative.events)&&Math.abs(row.fraction-g.fraction)<=.02+1e-10;
 return exact?{status:'correct',reason:'The displayed representative DNA and proportion agree with GT; model uncertainty may remain.'}:{status:'covered',reason:'GT family base proportions are within the reported bounds; individual sequence linkage remains unresolved. This is not a quantification pass.'};
}
