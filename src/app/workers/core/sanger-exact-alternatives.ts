import {normalizeMoleculeEdits} from './sanger-edit-normalize';
import {SangerAnnotation,SangerExactAlternatives,SangerExactOption,SangerJointEvent} from '../../models/sanger-analysis';
import {sangerSequenceKey} from '../../utils/sanger-sequence-groups';
import {exactFromNumber,exactHaplotypeRepresentatives,exactSum,exactWeightedEquality} from './exact-mixture-fiber';
import {materialize,Model,PathEdit,Prepared} from './sanger-path-common';
interface Row {model:Model;fraction:number}
/** Same-site insertion haplotype fibers and identical-sequence deletion placements. Every offered representative is
 * certified under the complete fitted emitter, including context and phase.
 * This enumerates verified representatives, not all possible allele solutions. */
export function discoverExactAlternatives(refs:Prepared[],annotations:SangerAnnotation[],rows:Row[][],toInput:(h:number,events:PathEdit[])=>SangerJointEvent[],emit:(sequence:string,p:number)=>number[],optionLimit=16):SangerExactAlternatives[]{
 const sets:SangerExactAlternatives[]=[];
 for(let h=0;h<refs.length;h++){
  const ref=refs[h],source=rows[h],annotation=annotations[h];if(!annotation.share||source.length!==annotation.components.length)continue;
  const base:SangerExactOption={id:'fit',components:annotation.components.map(c=>({...c})),exactFractions:source.map(c=>exactFromNumber(c.fraction))},options=[base],seen=new Set<string>(),columnCache=new Map<string,number[]>();
  const columnLength=Math.max(ref.sequence.length,...source.map(row=>row.model.sequence.length))+2;
  const column=(sequence:string)=>{let col=columnCache.get(sequence);if(!col){col=Array.from({length:columnLength},(_,p)=>emit(sequence,p)).flat();columnCache.set(sequence,col);}return col;};

  for(let t=0;t<ref.cuts.length&&options.length<optionLimit;t++){
   const targetLimit=1+Math.floor((optionLimit-1)*(t+1)/ref.cuts.length);
   const families=new Map<string,{indices:number[];background:PathEdit[];kind:'insert';start:number;end:number;length:number}>();
   source.forEach((row,index)=>{
    const nearest=(e:PathEdit)=>ref.cuts.map(c=>Math.min(Math.abs(e.start-c),Math.abs(e.end-c))).reduce((best,d,i,ds)=>d<ds[best]?i:best,0),local=row.model.events.filter(e=>nearest(e)===t),background=row.model.events.filter(e=>nearest(e)!==t);if(local.length>1||row.model.sequence.includes('N')||row.fraction<=0||annotation.components[index].unresolved)return;
    const e=local[0];if(!e||e.start!==e.end||!e.inserted.length)return;
    const kind='insert' as const,key=JSON.stringify(background)+':'+e.start+':'+e.inserted.length;let family=families.get(key);if(!family){family={indices:[],background,kind,start:e.start,end:e.end,length:e.inserted.length};families.set(key,family);}family.indices.push(index);

   });
   for(const family of families.values()){
    if(options.length>=targetLimit)break;if(family.indices.length<2||!Number.isFinite(family.start)||family.end<family.start)continue;
    const {start,end,background}=family,length=family.length;if(!length||length>40)continue;
    const localWords=family.indices.map(index=>{const row=source[index],local=row.model.events.filter(e=>!background.some(b=>JSON.stringify(b)===JSON.stringify(e))),sequence=materialize(ref.sequence.slice(start,end),local.map(e=>({...e,start:e.start-start,end:e.end-start})));return {sequence:sequence||'',weight:row.fraction};});
    if(localWords.some(w=>w.sequence.length!==length))continue;const reps=exactHaplotypeRepresentatives(localWords,optionLimit),old=family.indices.map(index=>({weight:exactFromNumber(source[index].fraction),column:column(source[index].model.sequence)}));
    for(const rep of reps.options.slice(1)){
     if(options.length>=targetLimit)break;const replacement=rep.map(word=>{const events=normalizeMoleculeEdits(ref.sequence,[...background,{start,end,inserted:word.sequence}]),sequence=materialize(ref.sequence,events);return {word,events,sequence};});if(replacement.some(r=>!r.sequence))continue;
     if(!exactWeightedEquality(old,replacement.map(r=>({weight:r.word.exactWeight,column:column(r.sequence!)}))))continue;
     const remove=new Set(family.indices),keep=source.flatMap((_,i)=>remove.has(i)?[]:[i]),components=keep.map(i=>base.components[i]),exactFractions=keep.map(i=>base.exactFractions[i]);
     for(const r of replacement){const events=toInput(h,r.events),sizes=(ref.input.cuts??[ref.input.cut]).map((_,target)=>events.filter(e=>e.target===target+1||e.target===0).reduce((sum,e)=>sum+e.inserted-(e.end-e.start),0));components.push({size:sizes.reduce((a,b)=>a+b,0),pathSizes:sizes,fraction:r.word.weight/annotation.share,events});exactFractions.push(r.word.exactWeight);}
     const key=JSON.stringify(components.map((c,i)=>[sangerSequenceKey(annotation.reference,c.events||[]),exactFractions[i]]).sort((a,b)=>String(a[0]).localeCompare(String(b[0]))));if(seen.has(key))continue;seen.add(key);options.push({id:'equivalent-'+options.length,components,exactFractions});
    }
   }
  }
  for(let i=0;i<source.length&&options.length<optionLimit;i++)for(const alias of source[i].model.aliases){
   if(annotation.components[i].unresolved||alias.reference!==h||materialize(ref.sequence,alias.events)!==source[i].model.sequence)continue;
   const events=toInput(h,alias.events);if(JSON.stringify(events)===JSON.stringify(base.components[i].events||[]))continue;
   const components=base.components.map((c,j)=>j===i?{...c,events}:c),key='placement:'+JSON.stringify(components.map(c=>c.events));if(seen.has(key))continue;seen.add(key);options.push({id:'placement-'+options.length,components,exactFractions:base.exactFractions});if(options.length>=optionLimit)break;
  }
  if(options.length<2)continue;const ambiguousSequences:Record<string,string[]>={},ambiguousPlacements:Record<string,string[]>={},targets:number[]=[];
  for(let t=0;t<(annotation.cuts?.length||1);t++){
   const maps=options.map(option=>{const values=new Map<string,typeof option.exactFractions>();option.components.forEach((c,i)=>{const key=sangerSequenceKey(annotation.reference,(c.events||[]).filter(e=>e.target===t+1||e.target===0));values.set(key,[...(values.get(key)||[]),option.exactFractions[i]]);});return new Map([...values].map(([key,v])=>[key,exactSum(v)]));}),keys=[...new Set(maps.flatMap(m=>[...m.keys()]))],changed=keys.filter(key=>new Set(maps.map(m=>{const v=m.get(key)||{n:'0',d:'1'};return v.n+'/'+v.d;})).size>1);
   const placementMaps=options.map(option=>{const m=new Map<string,string[]>();for(const c of option.components){const es=(c.events||[]).filter(e=>e.target===t+1||e.target===0),key=sangerSequenceKey(annotation.reference,es);m.set(key,[...(m.get(key)||[]),JSON.stringify(es)]);}return m;}),placementKeys=keys.filter(key=>placementMaps.every(m=>m.has(key))&&new Set(placementMaps.map(m=>JSON.stringify([...new Set(m.get(key)!)].sort()))).size>1);
   if(changed.length||placementKeys.length){targets.push(t+1);ambiguousSequences[String(t+1)]=changed;ambiguousPlacements[String(t+1)]=placementKeys;}
  }
  if(targets.length)sets.push({reference:annotation.name,certificate:'exact-rational-emission',selected:-1,options,targets,ambiguousSequences,ambiguousPlacements,exhaustive:false});
 }
 return sets;
}
