import {SangerAnnotation,SangerExactAlternatives} from '../models/sanger-analysis';
import {SangerAnnotationRow,sangerAnnotationViews} from './sanger-annotation';
import {projectTargetAnnotation} from './sanger-target-profile';

// A group is a final DNA sequence, not a particular edit-coordinate history.
export const alternativeRowKey=(row:SangerAnnotationRow)=>row.sequenceKey?JSON.stringify(['sequence',row.sequenceKey]):JSON.stringify([row.id,(row.representative?.events||[]).map(e=>[e.start,e.end,e.inserted,e.insertedSequence||''])]);
const shownFraction=(row:SangerAnnotationRow)=>(100*row.fraction).toFixed(2);
const rowShape=(row:SangerAnnotationRow)=>JSON.stringify((row.representative?.events||[]).map(e=>[e.start,e.end,e.inserted]));
const distributionKey=(rows:SangerAnnotationRow[])=>JSON.stringify(rows.map(row=>[alternativeRowKey(row),shownFraction(row)]).sort((a,b)=>a[0].localeCompare(b[0])));
export interface GroupAlternative {option:number;rows:SangerAnnotationRow[];key:string;affectedKeys:string[];selected:boolean}
interface OptionView {option:number;rows:SangerAnnotationRow[];key:string;values:Map<string,string>}
/** Presentation-only deduplication at the UI's 0.01 percentage-point precision.
 * Return existing certified option indices. Never synthesize a partial fit or
 * remove/reorder the engine's options. Target-invariant groups stay out of menus. */
export class SangerGroupAlternatives {
 readonly options:OptionView[]=[];
 readonly ambiguousKeys=new Set<string>();
 readonly entryKeys=new Set<string>();
 private readonly current:OptionView;
 private readonly menus=new Map<string,GroupAlternative[]>();
 constructor(annotation:SangerAnnotation,target:number,readonly set:SangerExactAlternatives,bounds:number|{left:number;right:number}){
  const seen=new Map<string,OptionView>();
  set.options.forEach((option,index)=>{
   const changed={...annotation,components:option.components},projected=(annotation.cuts?.length||1)>1?projectTargetAnnotation(changed,target):changed;
   const rows=sangerAnnotationViews([projected],0,bounds,true)[0].rows,key=distributionKey(rows),view={option:index,rows,key,values:new Map(rows.map(row=>[alternativeRowKey(row),shownFraction(row)]))};
   const previous=seen.get(key);if(!previous){seen.set(key,view);this.options.push(view);}else if(index===set.selected){previous.option=index;}
  });
  this.current=this.options.find(o=>o.option===set.selected)||this.options.find(o=>o.option===0)||this.options[0];
  if(!this.current)return;
  for(const option of this.options)for(const key of this.difference(option))this.ambiguousKeys.add(key);
  // One visible entry point per connected equivalence family. Companion WT or
  // indel mass changes belong inside its menu, not in additional yellow rows.
  const parents=new Map<string,string>(),find=(key:string):string=>{const parent=parents.get(key);if(!parent){parents.set(key,key);return key;}if(parent===key)return key;const root=find(parent);parents.set(key,root);return root;};
  for(const option of this.options){const keys=[...this.difference(option)];for(const key of keys)parents.set(find(key),find(keys[0]));}
  const families=new Map<string,SangerAnnotationRow[]>();for(const row of this.current.rows){const key=alternativeRowKey(row);if(!this.ambiguousKeys.has(key))continue;const root=find(key),rows=families.get(root)||[];rows.push(row);families.set(root,rows);}
  const substitution=(row:SangerAnnotationRow)=>!!row.representative?.events.length&&row.representative.events.every(e=>e.inserted===e.end-e.start),burden=(row:SangerAnnotationRow)=>(row.representative?.events||[]).reduce((n,e)=>n+Math.max(e.inserted,e.end-e.start),0);
  for(const rows of families.values()){rows.sort((a,b)=>Number(!a.edited)-Number(!b.edited)||Number(!substitution(a))-Number(!substitution(b))||burden(a)-burden(b)||a.fraction-b.fraction||alternativeRowKey(a).localeCompare(alternativeRowKey(b)));this.entryKeys.add(alternativeRowKey(rows[0]));}

 }
 private difference(option:OptionView){return new Set([...new Set([...this.current.values.keys(),...option.values.keys()])].filter(key=>this.current.values.get(key)!==option.values.get(key)));}
 isAmbiguous(row:SangerAnnotationRow){if(!this.current)return false;return this.ambiguousKeys.has(alternativeRowKey(row))||!!row.sequenceKey&&this.current.rows.some(other=>other.sequenceKey===row.sequenceKey&&this.ambiguousKeys.has(alternativeRowKey(other)))||row.id.startsWith('exact-insertion:')&&this.options.some(option=>option.rows.some(other=>rowShape(other)===rowShape(row)&&this.ambiguousKeys.has(alternativeRowKey(other))));}
 isEntry(row:SangerAnnotationRow){return this.entryKeys.has(alternativeRowKey(row))||row.id.startsWith('exact-insertion:')&&this.current?.rows.some(other=>rowShape(other)===rowShape(row)&&this.entryKeys.has(alternativeRowKey(other)));}
 forGroup(row:SangerAnnotationRow):GroupAlternative[]{
  const rowKey=alternativeRowKey(row),saved=this.menus.get(rowKey);if(saved)return saved;
  if(!this.current||!this.isAmbiguous(row))return [];
  const seeds=new Set(this.current.rows.filter(other=>alternativeRowKey(other)===rowKey||!!row.sequenceKey&&other.sequenceKey===row.sequenceKey||row.id.startsWith('exact-insertion:')&&rowShape(other)===rowShape(row)).map(alternativeRowKey));
  const changes=this.options.map(option=>({option,keys:this.difference(option)})).filter(change=>[...seeds].some(key=>change.keys.has(key)));
  // Independent edits offered in other options must not appear in this group's
  // dropdown. Minimal observed changes isolate its certified replacement family.
  const minimal=changes.filter(change=>!changes.some(other=>other!==change&&other.keys.size<change.keys.size&&[...other.keys].every(key=>change.keys.has(key))));
  const family=new Set([...seeds,...minimal.flatMap(change=>[...change.keys])]);
  const seen=new Set<string>(),candidates:GroupAlternative[]=[];
  for(const option of this.options){const affected=this.difference(option);if([...affected].some(key=>!family.has(key)))continue;
   const rows=option.rows.filter(other=>family.has(alternativeRowKey(other))),key=distributionKey(rows);if(!rows.length||seen.has(key))continue;seen.add(key);
   candidates.push({option:option.option,rows,key,affectedKeys:[...new Set([...family,...affected])],selected:option.key===this.current.key});
  }
  this.menus.set(rowKey,candidates);return candidates;
 }
}

export function alternativeGroupLabel(row:SangerAnnotationRow){
 const events=row.representative?.events||[];
 if(!events.length)return row.edited?row.classification:'WT';
 return events.map(e=>e.start===e.end?'+'+(e.insertedSequence||'N'.repeat(e.inserted))+' @'+e.start:e.inserted===e.end-e.start?'Substitution '+(e.start+1)+'–'+e.end+' → '+(e.insertedSequence||'N'.repeat(e.inserted)):!e.inserted?'−'+(e.end-e.start)+' bp @'+(e.start+1):'−'+(e.end-e.start)+' / +'+(e.insertedSequence||'N'.repeat(e.inserted))+' @'+(e.start+1)).join(' · ');
}
