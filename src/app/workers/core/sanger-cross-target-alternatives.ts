import{SangerAnnotation,SangerExactOption,SangerJointEvent}from'../../models/sanger-analysis';
import{Prepared,Model,PathEdit,materialize}from'./sanger-path-common';
import{crossClassAlternatives}from'./sanger-cross-class-alternatives';
import{exactFromNumber,exactSum,exactSubtract,exactNumber}from'./exact-mixture-fiber';
import{sangerSequenceKey}from'../../utils/sanger-sequence-groups';
interface Row{model:Model;fraction:number}
/** Hold the other target's molecule background fixed. Its upstream length
 * shifts the local coordinate system but never changes the full emitter used
 * for certification. No target independence or truth is used. */
export function crossTargetAlternatives(ref:Prepared,annotation:SangerAnnotation,source:Row[],target:number,toInput:(events:PathEdit[])=>SangerJointEvent[],emit:(sequence:string,p:number)=>number[],limit=4):SangerExactOption[]{
 if(ref.cuts.length!==2||ref.reverse)return[];
 const nearest=(e:PathEdit)=>Math.abs(e.start-ref.cuts[0])<=Math.abs(e.start-ref.cuts[1])?0:1,backgrounds=new Map<string,{indices:number[];events:PathEdit[]}>(),weights=source.map(r=>exactFromNumber(r.fraction));
 source.forEach((row,i)=>{if(annotation.components[i].unresolved)return;const events=row.model.events.filter(e=>nearest(e)!==target),key=JSON.stringify(events),group=backgrounds.get(key)||{indices:[],events};group.indices.push(i);backgrounds.set(key,group);});
 const localKey=(events:SangerJointEvent[])=>sangerSequenceKey(annotation.reference,events.filter(e=>e.target===target+1||e.target===0)),sourceKeys=[...new Set(annotation.components.map(c=>localKey(c.events||[])))],baseMass=sourceKeys.map(key=>exactSum(annotation.components.flatMap((c,i)=>localKey(c.events||[])===key?[weights[i]]:[]))),candidates:{option:SangerExactOption;mask:number}[]=[];
 for(const background of [...backgrounds.values()].sort((a,b)=>b.indices.length-a.indices.length)){
  if(background.indices.length<2)continue;
  const sequence=materialize(ref.sequence,background.events);if(!sequence)continue;const shift=background.events.filter(e=>e.end<=ref.cuts[target]).reduce((sum,e)=>sum+e.inserted.length-e.end+e.start,0),cut=ref.cuts[target]+shift,share=background.indices.reduce((sum,i)=>sum+source[i].fraction,0),rows=background.indices.map(i=>{const events=source[i].model.events.filter(e=>nearest(e)===target).map(e=>({...e,start:e.start+shift,end:e.end+shift}));return{fraction:source[i].fraction,model:{...source[i].model,events,aliases:[]}};}),virtual:Prepared={...ref,sequence,cuts:[cut],input:{name:ref.name,sequence,cut},offset:0,reverse:false},localInput=(es:PathEdit[])=>es.map(e=>({target:1,start:e.start,end:e.end,inserted:e.inserted.length,insertedSequence:e.inserted})),a:SangerAnnotation={name:ref.name,reference:sequence,cut,cuts:[cut],share,estimable:true,direction:'forward',components:rows.map(r=>({fraction:r.fraction/share,size:r.model.events.reduce((sum,e)=>sum+e.inserted.length-e.end+e.start,0),events:localInput(r.model.events)}))};
  for(const local of crossClassAlternatives(virtual,a,rows,localInput,emit,limit)){
   const removed=new Set(background.indices),keep=source.flatMap((_,i)=>removed.has(i)?[]:[i]),components=keep.map(i=>annotation.components[i]),exactFractions=keep.map(i=>weights[i]);
   local.components.forEach((c,i)=>{const events=[...background.events,...(c.events||[]).map(e=>({start:e.start-shift,end:e.end-shift,inserted:e.insertedSequence||''}))],converted=toInput(events),sizes=(annotation.cuts||[annotation.cut]).map((_,t)=>converted.filter(e=>e.target===t+1||e.target===0).reduce((sum,e)=>sum+e.inserted-e.end+e.start,0));components.push({fraction:exactNumber(local.exactFractions[i])/annotation.share,size:sizes.reduce((a,b)=>a+b,0),pathSizes:sizes,events:converted});exactFractions.push(local.exactFractions[i]);});
   const mask=sourceKeys.reduce((mask,key,i)=>exactSubtract(exactSum(components.flatMap((c,k)=>localKey(c.events||[])===key?[exactFractions[k]]:[])),baseMass[i]).n==='0'?mask:mask|2**i,0);if(mask)candidates.push({mask,option:{id:'target-cross-'+target+'-'+candidates.length,components,exactFractions}});
  }
 }
 const chosen:typeof candidates=[],count=(n:number)=>n.toString(2).replace(/0/g,'').length;let covered=0;
 while(candidates.length&&chosen.length<limit){candidates.sort((a,b)=>count(b.mask&~covered)-count(a.mask&~covered)||a.option.components.length-b.option.components.length);const next=candidates.shift()!;chosen.push(next);covered|=next.mask;}
 return chosen.map((c,i)=>({...c.option,id:'target-cross-'+target+'-'+i}));
}
