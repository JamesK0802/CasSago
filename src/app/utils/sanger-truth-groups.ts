import {BenchmarkTruth} from './sanger-benchmark';
import {SangerAnnotation,SangerJointEvent} from '../models/sanger-analysis';
import {sangerAnnotationViews,SangerAnnotationRow} from './sanger-annotation';
import {reverseComplement} from './biology.utils';
export interface TruthGroupView {label:string;description:string;sequence:string;row:SangerAnnotationRow;fraction:number;edited:boolean}
/** Display adapter only: ground truth never enters the analysis worker. */
export function sangerTruthGroups(data:BenchmarkTruth|null,annotation:SangerAnnotation,radius:number|{left:number;right:number},focusCut=annotation.cut):TruthGroupView[]{
 if(!data)return [];
 return data.truth.filter(t=>t.reference===annotation.name).map(t=>{
  const recipe=data.recipe?.find(r=>r.name===annotation.name),reverse=recipe?recipe.sequence!==annotation.reference&&recipe.sequence===reverseComplement(annotation.reference):t.sequence!==t.inputSequence&&reverseComplement(t.sequence)===t.inputSequence;
  const cuts=recipe?.cuts||(annotation.cuts||[annotation.cut]).map(c=>reverse?annotation.reference.length-c:c);
  const events:SangerJointEvent[]=t.events.map(e=>{const start=e.kind==='substitution'?e.selectionStart:cuts[e.target]+e.offset,end=Number(start)+(e.kind==='insertion'?0:e.kind==='deletion'?e.length:(e.before||'').length),inserted=e.kind==='deletion'?'':e.sequence;return {target:e.target+1,start:reverse?annotation.reference.length-end:Number(start),end:reverse?annotation.reference.length-Number(start):end,inserted:inserted.length,insertedSequence:reverse?reverseComplement(inserted):inserted};});
  const wt=t.description==='WT',sequence=t.inputSequence||'',valid=sequence&&materializeGroup(annotation.reference,events)===sequence;
  const row=sangerAnnotationViews([{...annotation,sequenceGrouping:false,cut:focusCut,cuts:[focusCut],estimable:true,components:[{size:t.pattern.split('/').map(Number).reduce((sum,v)=>sum+v,0),fraction:t.percent/100,pathSizes:t.pattern.includes('/')?t.pattern.split('/').map(Number):undefined,events:valid?events:[],unresolved:!valid&&!wt}]}],0,radius)[0].rows[0];
  return {label:t.allele,description:t.description,sequence,row,fraction:t.percent/100,edited:!wt};
 }).filter(t=>!!t.row);
}
export function materializeGroup(reference:string,events:SangerJointEvent[]):string{
 let result=reference;for(const e of [...events].sort((a,b)=>b.start-a.start))result=result.slice(0,e.start)+(e.insertedSequence||'')+result.slice(e.end);return result;
}
