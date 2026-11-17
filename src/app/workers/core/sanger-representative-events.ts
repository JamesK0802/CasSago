import type {SangerJointEvent} from '../../models/sanger-analysis';
/** Collapse touching edit histories only when their complete final DNA stays
 * identical. A shifted insertion plus compensating SNPs is one DNA proposal,
 * not evidence for additional biological substitutions. Targets stay separate.
 * No genotype mass, fitted base marginal or alternative is changed. */
export function representativeEvents(reference:string,input:SangerJointEvent[]):SangerJointEvent[]{
 if(input.some(e=>e.insertedSequence?.length!==e.inserted))return input;
 const sorted=input.map(e=>({...e})).sort((a,b)=>a.start-b.start||a.end-b.end),merged:SangerJointEvent[]=[];
 for(const event of sorted){const old=merged.at(-1);if(old&&old.target===event.target&&old.end===event.start){old.end=event.end;old.insertedSequence=(old.insertedSequence||'')+(event.insertedSequence||'');old.inserted=old.insertedSequence.length;}else merged.push(event);}
 const out=merged.map(e=>{let {start,end}=e,word=e.insertedSequence!;while(start<end&&word.length&&reference[start]===word[0]){start++;word=word.slice(1);}while(start<end&&word.length&&reference[end-1]===word.at(-1)){end--;word=word.slice(0,-1);}return{...e,start,end,inserted:word.length,insertedSequence:word};}).filter(e=>e.start!==e.end||e.inserted);
 const apply=(events:SangerJointEvent[])=>{let dna=reference;for(const e of [...events].sort((a,b)=>b.start-a.start||b.end-a.end))dna=dna.slice(0,e.start)+(e.insertedSequence||'')+dna.slice(e.end);return dna;};
 return apply(input)===apply(out)?out:input;
}
