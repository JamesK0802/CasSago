import type {SangerAnnotation,SangerJointEvent} from '../models/sanger-analysis';
import type {BenchmarkTruth} from './sanger-benchmark';
import type {SynthTruth} from './sanger-synthetic-benchmark';
import {sangerEditedSequence,sangerSequenceKey} from './sanger-sequence-groups';
import {reverseComplement} from './biology.utils';
export function targetEvents(events:SangerJointEvent[],target:number){return events.filter(e=>e.target===target+1||e.target===0);}
export function projectTargetAnnotation(a:SangerAnnotation,target:number):SangerAnnotation {
 const cut=a.cuts?.[target]??a.cut;return {...a,cut,cuts:[cut],components:a.components.map(c=>{const events=targetEvents(c.events||[],target),size=c.pathSizes?.[target]??events.reduce((s,e)=>s+e.inserted-(e.end-e.start),0);return {...c,size,pathSizes:[size],events};})};
}
/** Marginalize whole-molecule truth for display/evaluation only. */
export function projectTargetTruth(data:BenchmarkTruth,a:SangerAnnotation,target:number):BenchmarkTruth {
 const recipe=data.recipe?.find(r=>r.name===a.name),reverse=recipe?recipe.sequence!==a.reference&&recipe.sequence===reverseComplement(a.reference):a.direction==='reverse',cuts=recipe?.cuts||a.cuts||[a.cut],profile=recipe?.profiles?.find(p=>p.target===target),groups=new Map<string,SynthTruth>();
 for(const t of data.truth.filter(t=>t.reference===a.name)){
  const events=t.events.filter(e=>e.target===target),converted:SangerJointEvent[]=events.map(e=>{const start=Number(e.kind==='substitution'?e.selectionStart:cuts[e.target]+e.offset),end=start+(e.kind==='insertion'?0:e.kind==='deletion'?e.length:(e.before||'').length),inserted=e.kind==='deletion'?'':e.sequence;return {target:e.target+1,start:reverse?a.reference.length-end:start,end:reverse?a.reference.length-start:end,inserted:inserted.length,insertedSequence:reverse?reverseComplement(inserted):inserted};});
  const inputSequence=t.inputSequence?sangerEditedSequence(a.reference,converted):'',pattern=t.pattern.split('/')[target]||'0',key=pattern+':'+(inputSequence||JSON.stringify(events)),old=groups.get(key),index=profile?.alleles.findIndex(g=>JSON.stringify(g.events)===JSON.stringify(events))??-1;
  if(old){old.percent+=t.percent;old.totalPercent+=t.totalPercent;continue;}
  const description=!events.length?'WT':events.map(e=>e.kind==='substitution'?`${e.before}→${e.sequence}`:`${e.kind} ${e.offset>=0?'+':''}${e.offset}: ${e.kind==='insertion'?e.sequence:e.length+' bp'}`).join('; ');
  groups.set(key,{...t,allele:!events.length?'WT':index>=0?'Group '+(index+1):t.allele,events,pattern,description,substitution:events.some(e=>e.kind==='substitution'),inputSequence,sequence:reverse?reverseComplement(inputSequence):inputSequence});
 }
 return {...data,truth:[...groups.values()]};
}
