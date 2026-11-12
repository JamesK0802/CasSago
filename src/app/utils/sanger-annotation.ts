import {sangerSequenceChoices,sangerSequenceKey,SANGER_SEQUENCE_CANDIDATE_LIMIT,SequenceState} from './sanger-sequence-groups';
import {SangerAnnotation, SangerJointEvent, SangerGroupSupport} from '../models/sanger-analysis';
import {AlignmentToken} from '../models/analysis.model';
import {resolveAnnotationCutSite} from '../components/result-dashboard/annotation-cut-site';

export interface SangerPlacementView {
 fraction:number; tokens:AlignmentToken[]; location:string; events:SangerJointEvent[];
}
export interface SangerAnnotationRow {
 marginalRole?:'wt'|'substitution';marginalProfile?:boolean;marginalMixture?:boolean;marginalFamilyId?:string;substitutionMixture?:boolean; marginal?:import('../models/sanger-analysis').SangerMarginalGroup;fractionRange?:[number,number];id:string;sequenceKey?:string;sequenceState:SequenceState;size:number; edited:boolean; patternLabel?:string; fraction:number; classification:string; tokens:AlignmentToken[];
 placements:SangerPlacementView[]; representative:SangerPlacementView|null; positionLabel:string; support:SangerGroupSupport|null;
}
export interface SangerAnnotationView {
 name:string; direction:'forward'|'reverse'; reference:string; start:number; end:number;
 cutIndex:number; cutIndices:number[]; share:number; estimable:boolean;
 rows:SangerAnnotationRow[]; hiddenFraction:number;
}
type Component=SangerAnnotation['components'][number];

export function tokensFor(a:SangerAnnotation,start:number,end:number,events?:SangerJointEvent[]):AlignmentToken[]{
 const tokens:AlignmentToken[]=[];
 const append=(type:AlignmentToken['type'],val:string)=>{
  const previous=tokens.at(-1);if(previous?.type===type)previous.val+=val;else tokens.push({type,val});
 };
 const substitution=(e:SangerJointEvent)=>e.end>e.start&&e.inserted===e.end-e.start&&e.insertedSequence?.length===e.inserted;
 for(let p=start;p<end;p++){
  for(const e of events||[])if(p===e.start&&e.inserted&&!substitution(e))append('insert',e.insertedSequence||'N'.repeat(e.inserted));
  const edit=events?.find(e=>p>=e.start&&p<e.end);
  if(edit&&substitution(edit)){const base=edit.insertedSequence![p-edit.start];append(base===a.reference[p]?'equal':'substitute',base);}else append(edit?'delete':'equal',a.reference[p]);
 }
 for(const e of events||[])if(e.start===end&&e.inserted&&!substitution(e))append('insert',e.insertedSequence||'N'.repeat(e.inserted));
 return tokens;
}
function locationFor(events:SangerJointEvent[]):string{
 return events.map(e=>e.inserted>0&&e.inserted===e.end-e.start?`Substitution bases ${e.start+1}–${e.end}`:e.end>e.start?e.end===e.start+1?`Deletion base ${e.start+1}`:`Deletion bases ${e.start+1}–${e.end}`:`Insertion after base ${e.start}`).join('; ');
}
/** Sequence-specific rows conserve fitted reference/size mass. Mixed-base
 * combinations remain model-derived candidates, not molecule-level phasing. */
export function sangerAnnotationViews(items:SangerAnnotation[],threshold:number,radius:number|{left:number;right:number}=30,showUnconfirmed=true):SangerAnnotationView[]{
 const cutoff=Math.max(0,Math.min(100,Number(threshold)||0))/100;
 return items.map(a=>{
  const cuts=a.cuts?.length?a.cuts:[a.cut],start=Math.max(0,Math.min(...cuts)-(typeof radius==='number'?radius:radius.left)),end=Math.min(a.reference.length,Math.max(...cuts)+(typeof radius==='number'?radius:radius.right));
  type Group={components:Component[];size:number;fraction:number;events:SangerJointEvent[];sequenceKey?:string;state:SequenceState};
  const grouped=new Map<string,Group>();
  const rank:Record<SequenceState,number>={'base-call':0,mixture:1,unphased:2,unresolved:3};
  for(const c of a.components){
   if(c.fraction<=0)continue;
   const pattern=c.pathSizes?.join('/')??`${c.size}:${c.size===0&&c.events?.length?'edited':'net'}`;
   const choices=c.unresolved||!c.events?.length?[{events:[],weight:1,state:c.size===0&&!c.unresolved?'base-call' as const:'unresolved' as const}]:sangerSequenceChoices(a.reference,c.events);
   for(const choice of choices){if(a.sequenceGrouping&&choice.events.length&&choice.state!=='unresolved')choice.state='unphased';const sequenceKey=choice.state==='unresolved'?undefined:sangerSequenceKey(a.reference,choice.events),key=pattern+':'+(sequenceKey||'unresolved'),fraction=c.fraction*choice.weight,previous=grouped.get(key);
    if(previous){previous.fraction+=fraction;previous.components.push(c);if(rank[choice.state]>rank[previous.state])previous.state=choice.state;}
    else grouped.set(key,{components:[c],size:c.size,fraction,events:choice.events,sequenceKey,state:choice.state});
   }
  }
  const ordered=[...grouped.entries()].sort(([,x],[,y])=>y.fraction-x.fraction||Math.abs(x.size)-Math.abs(y.size)||x.size-y.size);
  // Keep a bounded number of concrete rows without discarding small-group mass.
  let concrete=0;const retained:typeof ordered=[],remainders=new Map<string,Group>();
  for(const [key,g] of ordered){if(g.state==='unresolved'||concrete<SANGER_SEQUENCE_CANDIDATE_LIMIT){retained.push([key,g]);if(g.state!=='unresolved')concrete++;}else{const pattern=g.components[0].pathSizes?.join('/')??String(g.size),old=remainders.get(pattern);if(old)old.fraction+=g.fraction;else remainders.set(pattern,{...g,events:[],sequenceKey:undefined,state:'unresolved'});}}
  for(const [pattern,g] of remainders){const old=retained.find(([,r])=>r.state==='unresolved'&&(r.components[0].pathSizes?.join('/')??String(r.size))===pattern);if(old)old[1].fraction+=g.fraction;else retained.push(['remainder:'+pattern,g]);}
  const displayable=a.estimable||(showUnconfirmed&&a.share>0);
  const rows:SangerAnnotationRow[]=displayable?retained.filter(([,g])=>g.fraction+1e-12>=cutoff).sort(([,x],[,y])=>y.fraction-x.fraction||Math.abs(x.size)-Math.abs(y.size)).map(([id,g])=>{
   const c=g.components[0],isBaseline=g.size===0&&!c.events?.length&&!c.unresolved,events=g.state==='unresolved'?[]:g.events;
   const placements:SangerPlacementView[]=events.length?[{fraction:g.fraction,events,tokens:tokensFor(a,start,end,events),location:locationFor(events)}]:[];
   const classification=g.state==='unresolved'?'Sequence / breakpoint unresolved':isBaseline?'No net indel model':g.state==='unphased'?'Unphased sequence candidate':g.size===0?(events.every(e=>e.inserted===e.end-e.start)?'Substitution':'Net-zero multi-target edit'):g.size%3===0?'In-frame indel':'Out-of-frame indel';
   const representative=placements[0]||null,positionLabel=isBaseline?'No net indel placement':!representative?'Sequence / position not estimated':`Best fit: ${representative.location}`;
   return {id,sequenceKey:g.sequenceKey,sequenceState:g.state,size:g.size,edited:!isBaseline,patternLabel:c.pathSizes?.map((v,t)=>`T${t+1}: ${v>0?'+':''}${v} bp`).join(' · '),fraction:g.fraction,classification,placements,representative,positionLabel,support:c.support||null,tokens:representative?.tokens||tokensFor(a,start,end)};
  }):[];
  const reference=a.reference.slice(start,end);
  return {name:a.name,direction:a.direction,reference,start,end,
   cutIndex:resolveAnnotationCutSite({ref_sequence:reference,cut_site_index:a.cut-start})!,cutIndices:cuts.map(c=>c-start),
   share:a.share,estimable:a.estimable,rows,
   hiddenFraction:displayable?retained.map(([,g])=>g).filter(g=>g.fraction+1e-12<cutoff).reduce((s,g)=>s+g.fraction,0):0};
 });
}
