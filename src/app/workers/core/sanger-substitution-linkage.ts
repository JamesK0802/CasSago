import {exactFromNumber,exactWeightedEquality} from './exact-mixture-fiber';
/** Exact fluorescence aliases, independent of truth labels. These witness ranges
 * are not posterior intervals or a claim that every possible alias was found. */
export interface SubstitutionLinkage {reference:number;target:number;blocks:number;fraction:number;wtRange:[number,number];modelSequence:string}
interface Edit {start:number;end:number;inserted:string}
interface Model {reference:number;sequence:string;events:Edit[];fraction:number}
export function substitutionLinkageAliases(refs:{sequence:string;cuts:number[]}[],chosen:Model[],column:(sequence:string)=>number[]):SubstitutionLinkage[]{
 const out:SubstitutionLinkage[]=[];
 const apply=(sequence:string,events:Edit[])=>{for(const e of [...events].sort((a,b)=>b.start-a.start))sequence=sequence.slice(0,e.start)+e.inserted+sequence.slice(e.end);return sequence;};
 for(const model of chosen){const ref=refs[model.reference],share=chosen.filter(m=>m.reference===model.reference).reduce((s,m)=>s+m.fraction,0);if(!share)continue;
  for(let ei=0;ei<model.events.length;ei++){const event=model.events[ei];if(event.end-event.start!==event.inserted.length||!event.inserted.length)continue;
   const blocks:Edit[]=[];for(let p=event.start;p<event.end;){if(ref.sequence[p]===event.inserted[p-event.start]){p++;continue;}const start=p;while(p<event.end&&ref.sequence[p]!==event.inserted[p-event.start])p++;blocks.push({start,end:p,inserted:event.inserted.slice(start-event.start,p-event.start)});}
   if(blocks.length===1&&blocks[0].end-blocks[0].start>1){const block=blocks.pop()!;for(let p=block.start;p<block.end;p++)blocks.push({start:p,end:p+1,inserted:event.inserted[p-event.start]});}
   if(blocks.length<2)continue;
   const target=ref.cuts.map(c=>Math.abs(c-event.start)).reduce((best,d,i,ds)=>d<ds[best]?i:best,0),scaffold=model.events.filter((_,i)=>i!==ei),base=apply(ref.sequence,scaffold),baseMass=chosen.filter(m=>m.reference===model.reference&&m.sequence===base).reduce((s,m)=>s+m.fraction,0),shift=Math.min(model.fraction,baseMass/(blocks.length-1));if(shift<1e-5)continue;
   const joined=column(model.sequence),empty=column(base),parts=blocks.map(e=>column(apply(ref.sequence,[...scaffold,e])));
   if(!exactWeightedEquality([{weight:exactFromNumber(1),column:joined},{weight:exactFromNumber(blocks.length-1),column:empty}],parts.map(column=>({weight:exactFromNumber(1),column}))))continue;
   const wt=chosen.filter(m=>m.reference===model.reference&&!m.events.some(e=>ref.cuts.map(c=>Math.abs(c-e.start)).indexOf(Math.min(...ref.cuts.map(c=>Math.abs(c-e.start))))===target)).reduce((s,m)=>s+m.fraction,0)/share;
   out.push({reference:model.reference,target:target+1,blocks:blocks.length,fraction:shift/share,wtRange:[Math.max(0,wt-(blocks.length-1)*shift/share),wt],modelSequence:model.sequence});
  }
 }
 return out;
}
