import{exactSum,exactSubtract,ExactFraction}from'../workers/core/exact-mixture-fiber';
import{canonicalSangerEvents}from'./sanger-placement';
import{sangerSequenceKey}from'./sanger-sequence-groups';
import {SangerAnnotation,SangerExactAlternatives}from'../models/sanger-analysis';
import {SangerAnnotationView,sangerAnnotationViews,tokensFor}from'./sanger-annotation';
import {projectTargetAnnotation,targetEvents}from'./sanger-target-profile';
/** Only an isolated, exactly certified same-site/length insertion fiber may
 * become an N block. Other ambiguities keep separate fitted groups and point
 * fractions; substitutions, changed origins and mass transfers are never pooled. */
export function exactDisplayView(annotation:SangerAnnotation,target:number,radius:number|{left:number;right:number},set?:SangerExactAlternatives):SangerAnnotationView {
 const projected=(annotation.cuts?.length||1)>1?projectTargetAnnotation(annotation,target):annotation,view=sangerAnnotationViews([projected],0,radius,true)[0];
 if(!set||set.selected>=0||set.certificate!=='exact-rational-emission')return view;
 const base=set.options[0];if(!base||base.exactFractions.length!==base.components.length)return view;
 const local=(c:typeof base.components[number])=>canonicalSangerEvents(annotation.reference,targetEvents(c.events||[],target)),geometry=(events:ReturnType<typeof local>)=>events.length===1&&events[0].start===events[0].end&&events[0].inserted>0&&/^[ACGT]+$/.test(events[0].insertedSequence||'')?events[0].start+':'+events[0].inserted:null;
 const equal=(a:Map<string,ExactFraction>,b:Map<string,ExactFraction>)=>a.size===b.size&&[...a].every(([key,v])=>b.has(key)&&exactSubtract(v,b.get(key)!).n==='0'),distribution=(option:typeof base,geo:string)=>{const inside=new Map<string,ExactFraction[]>(),outside=new Map<string,ExactFraction[]>(),words=new Set<string>();option.components.forEach((c,i)=>{const events=local(c),isInside=geometry(events)===geo,key=sangerSequenceKey(annotation.reference,events),map=isInside?inside:outside;map.set(key,[...(map.get(key)||[]),option.exactFractions[i]]);if(isInside)words.add(events[0].insertedSequence!);});return{inside:new Map([...inside].map(([k,v])=>[k,exactSum(v)])),outside:new Map([...outside].map(([k,v])=>[k,exactSum(v)])),mass:exactSum([...inside.values()].flat()),words};};
 let rows=view.rows;
 for(const geo of new Set(base.components.map(c=>geometry(local(c))).filter((g):g is string=>!!g))){
  const original=distribution(base,geo);if(original.words.size<2)continue;const allowed=set.options.slice(1).filter(o=>o.exactFractions.length===o.components.length).map(o=>distribution(o,geo)).filter(d=>equal(original.outside,d.outside)&&exactSubtract(original.mass,d.mass).n==='0'&&!equal(original.inside,d.inside));if(!allowed.length)continue;
  const members=rows.filter(r=>geometry(r.representative?.events||[])===geo);if(members.length<2)continue;const words=[...new Set([original,...allowed].flatMap(d=>[...d.words]))],length=words[0].length,insertedSequence=Array.from({length},(_,i)=>words.every(w=>w[i]===words[0][i])?words[0][i]:'N').join('');if(!insertedSequence.includes('N'))continue;
  const start=Number(geo.split(':')[0]),events=[{target:target+1,start,end:start,inserted:length,insertedSequence}],fraction=members.reduce((sum,r)=>sum+r.fraction,0),tokens=tokensFor(projected,view.start,view.end,events),memberIds=new Set(members.map(r=>r.id));rows=rows.filter(r=>!memberIds.has(r.id));rows.push({...members[0],id:'exact-insertion:'+geo,sequenceKey:undefined,sequenceState:'unphased',fraction,classification:'Signal-identical insertion sequences at the same site',positionLabel:'Insertion sequence unresolved',tokens,placements:[],representative:{fraction,events,tokens,location:'Insertion sequence unresolved'}});
 }
 return rows===view.rows?view:{...view,rows:rows.sort((a,b)=>b.fraction-a.fraction)};
}

export function exactMetricRange(annotation:SangerAnnotation,target:number,set:SangerExactAlternatives,kind:'wt'|'substitution'):[number,number]{
 const values=set.options.map(o=>o.components.reduce((sum,c)=>{const es=targetEvents(c.events||[],target),matches=kind==='wt'?!es.length:!!es.length&&es.every(e=>e.inserted===e.end-e.start);return sum+(matches?c.fraction:0);},0));return [Math.min(...values),Math.max(...values)];
}

export function exactFractionLabel(row:{fraction:number;fractionRange?:[number,number]}):string{return row.fraction<.0001?'<0.01':(row.fraction*100).toFixed(2);}
