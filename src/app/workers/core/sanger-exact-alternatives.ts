import{composeExactOptions}from'./sanger-compose-exact-options';
import type{SangerLinearOptimizer}from'./sanger-discrete-phase';
import{reverseComplement}from'../../utils/biology.utils';
import {crossClassAlternatives} from './sanger-cross-class-alternatives';
import {crossTargetAlternatives} from './sanger-cross-target-alternatives';
import {normalizeMoleculeEdits} from './sanger-edit-normalize';
import {SangerAnnotation,SangerExactAlternatives,SangerExactOption,SangerJointEvent} from '../../models/sanger-analysis';
import {sangerSequenceKey} from '../../utils/sanger-sequence-groups';
import {exactFromNumber,exactHaplotypeRepresentatives,exactSum,exactWeightedEquality} from './exact-mixture-fiber';
import {materialize,Model,PathEdit,Prepared} from './sanger-path-common';
interface Row {model:Model;fraction:number}
/** Same-site insertion haplotype fibers and distinct final sequence profiles. Every offered representative is
 * certified under the complete fitted emitter, including context and phase.
 * This enumerates verified representatives, not all possible allele solutions. */
export function discoverExactAlternatives(refs:Prepared[],annotations:SangerAnnotation[],rows:Row[][],toInput:(h:number,events:PathEdit[])=>SangerJointEvent[],emit:(sequence:string,p:number)=>number[],optionLimit=16,verify?:(reference:number,base:SangerExactOption,candidate:SangerExactOption)=>boolean,progress:(stage:string)=>void=()=>{},optimizer?:SangerLinearOptimizer):SangerExactAlternatives[]{
 const sets:SangerExactAlternatives[]=[];
 for(let h=0;h<refs.length;h++){
  const ref=refs[h],source=rows[h],annotation=annotations[h];if(!annotation.share||source.length!==annotation.components.length)continue;
  const nearest=(e:PathEdit)=>{const distances=ref.cuts.map(c=>Math.min(Math.abs(e.start-c),Math.abs(e.end-c)));return distances.indexOf(Math.min(...distances));};
  // Cost follows the disjoint target backgrounds actually searched, rather
  // than squaring the whole cross-target molecule product.
  const workload=ref.cuts.length===2?Math.max(...ref.cuts.map((_,target)=>{const backgrounds=new Map<string,number>();for(const row of source){const key=JSON.stringify(row.model.events.filter(e=>nearest(e)!==target));backgrounds.set(key,(backgrounds.get(key)||0)+1);}return [...backgrounds.values()].reduce((n,count)=>n+count*count,0);})):source.length**2;
  if(workload*ref.sequence.length>20_000_000)continue;
  const base:SangerExactOption={id:'fit',components:annotation.components.map(c=>({...c})),exactFractions:source.map(c=>exactFromNumber(c.fraction))},options=[base],seen=new Set<string>(),columnCache=new Map<string,number[]>();
  const columnLength=Math.max(ref.sequence.length,...source.map(row=>row.model.sequence.length))+2;
  const column=(sequence:string)=>{let col=columnCache.get(sequence);if(!col){col=Array.from({length:columnLength},(_,p)=>emit(sequence,p)).flat();columnCache.set(sequence,col);}return col;};

  for(let t=0;t<ref.cuts.length&&options.length<optionLimit;t++){
   progress('Certifying insertion alternatives · reference '+(h+1)+' · target '+(t+1));
   const targetLimit=1+Math.floor((optionLimit-1)*(t+1)/(2*ref.cuts.length));
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
  // These discovery routines already certify the complete fitted emitter.
  // Connect them to production reporting before placement-only aliases use the
  // option budget. This changes uncertainty reporting, never fitted coefficients.
  for(let target=0;target<ref.cuts.length&&options.length<optionLimit;target++){
   progress('Certifying cross-class alternatives · reference '+(h+1)+' · target '+(target+1));
   const remaining=Math.min(Math.ceil((optionLimit-1)/2),optionLimit-options.length),toLocal=(events:PathEdit[])=>toInput(h,events),alternatives=ref.cuts.length===2?crossTargetAlternatives(ref,annotation,source,target,toLocal,emit,remaining):crossClassAlternatives(ref,annotation,source,toLocal,emit,remaining);
   for(const option of alternatives){if(options.length>=optionLimit)break;const key='sequence:'+JSON.stringify(option.components.map((c,i)=>[sangerSequenceKey(annotation.reference,c.events||[]),option.exactFractions[i]]).sort((a,b)=>String(a[0]).localeCompare(String(b[0]))));if(seen.has(key))continue;seen.add(key);options.push(option);}
  }
  // Coordinate aliases producing the same final DNA are one biological group.
  // They must not consume alternative capacity or trigger signal certificates.
  // Independent certified moves can be combined into a complete target
  // profile. A single pair witness may repair only one of several linked
  // substitutions; retaining only those isolated moves hides valid profiles.
  if(optimizer){
   const substitutionCost=(c:SangerExactOption['components'][number])=>(c.events||[]).reduce((n,e)=>{const word=e.insertedSequence||'';return n+(word.length===e.end-e.start?Array.from(word,(b,k)=>Number(b!==annotation.reference[e.start+k])).reduce((a,b)=>a+b,0):0);},0),toPrepared=(c:SangerExactOption['components'][number])=>(c.events||[]).map(e=>({start:(ref.reverse?ref.input.sequence.length-e.end:e.start)-ref.offset,end:(ref.reverse?ref.input.sequence.length-e.start:e.end)-ref.offset,inserted:ref.reverse?reverseComplement(e.insertedSequence||''):e.insertedSequence||''})).sort((a,b)=>a.start-b.start);
   let current=base,choices=options.slice();
   for(let round=0;round<64;round++){
    const valid=verify?choices.filter(o=>o===current||verify(h,current,o)):choices,combined=composeExactOptions(current,valid,annotation.share,substitutionCost,optimizer);if(!combined||verify&&!verify(h,current,combined))break;
    combined.id='combined-equivalent-'+round;options.push(combined);current=combined;progress('Combining certified target alternatives · reference '+(h+1)+' · '+(round+1));
    const a={...annotation,components:current.components},localRows=current.components.map((c,i)=>{const events=toPrepared(c);return{fraction:Number(current.exactFractions[i].n)/Number(current.exactFractions[i].d),model:{reference:c.unresolved?-1:h,events,sequence:c.unresolved?'':materialize(ref.sequence,events)!,aliases:[]}};});choices=[current];
    for(let target=0;target<ref.cuts.length;target++){const found=ref.cuts.length===2?crossTargetAlternatives(ref,a,localRows,target,events=>toInput(h,events),emit,optionLimit,current.exactFractions):crossClassAlternatives(ref,a,localRows,events=>toInput(h,events),emit,optionLimit,3,10,current.exactFractions);for(const o of found){if(verify&&!verify(h,current,o))continue;choices.push(o);options.push({...o,id:'combined-witness-'+round+'-'+options.length});}}
   }
  }
  if(verify)for(let i=options.length-1;i>0;i--){progress('Verifying whole-emitter alternative · reference '+(h+1)+' · option '+i);if(!verify(h,base,options[i]))options.splice(i,1);}
  // Deduplicate only after emission verification: different coordinate histories
  // can have different donor-control proxies, so keep the first VALID witness.
  // Exact summed masses distinguish genuinely different abundance profiles.
  const dnaCache=new Map<string,string>(),profiles=new Set<string>();
  for(let i=0;i<options.length;i++){
   const option=options[i],mass=new Map<string,typeof option.exactFractions>();
   option.components.forEach((c,j)=>{const eventKey=JSON.stringify(c.events||[]);let key=c.unresolved?'__unresolved__':dnaCache.get(eventKey);if(key===undefined){key=sangerSequenceKey(annotation.reference,c.events||[]);dnaCache.set(eventKey,key);}mass.set(key,[...(mass.get(key)||[]),option.exactFractions[j]]);});
   const key=JSON.stringify([...mass].map(([sequence,weights])=>[sequence,exactSum(weights)]as const).filter(([,weight])=>weight.n!=='0').sort((a,b)=>a[0].localeCompare(b[0])));
   if(profiles.has(key)){options.splice(i--,1);continue;}profiles.add(key);
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
