import {SangerAnalysisResult,SangerJointEvent} from '../models/sanger-analysis';
export interface SangerWindowSummary {
 window:number;targets:NonNullable<SangerAnalysisResult['joint']>['targets'];
 references:{name:string;share:number;wt:number|null;indel:number|null;inFrame:number|null;outOfFrame:number|null;unresolvedFraction:number;components:{size:number;fraction:number}[]}[];
}
/** Summarise immutable full-molecule fits within editable target windows.
 * No fitting, renormalisation, discarded raw event or new confidence claim. */
export function summarizeSangerWindow(result:SangerAnalysisResult,window:number):SangerWindowSummary|null {
 if(!result.joint||!Number.isInteger(window)||window<1||window>1000)return null;
 const joint=result.joint,within=(e:SangerJointEvent,cut:number)=>e.end>e.start?e.start<cut+Math.ceil(window/2)&&e.end>cut-Math.floor(window/2):e.start>=cut-Math.floor(window/2)&&e.start<cut+Math.ceil(window/2);
 const refs=(result.pooled?.amplicons??[]).map(ref=>{
  const targets=joint.targets.filter(t=>t.reference===ref.name),groups=joint.groups.filter(g=>g.reference===ref.name),bins=new Map<number,number>();let wt=0,edited=0,inFrame=0,outOfFrame=0,unresolved=0;
  for(const g of groups){
   if(g.netSize===null){unresolved+=g.fraction;continue;}
   const alternatives=g.alternatives.length?g.alternatives:[{reference:ref.name,events:g.events}];
   const states=alternatives.map(a=>{if(a.reference!==ref.name)return null;const events=a.events.filter(e=>targets.some(t=>within(e,t.cut))),size=events.reduce((sum,e)=>sum+e.inserted-e.end+e.start,0);return {size,edited:events.length>0};});
   const first=states[0];if(!first||states.some(s=>!s||s.size!==first.size||s.edited!==first.edited)){unresolved+=g.fraction;continue;}
   bins.set(first.size,(bins.get(first.size)||0)+g.fraction);if(!first.edited)wt+=g.fraction;else{edited+=g.fraction;if(first.size%3===0)inFrame+=g.fraction;else outOfFrame+=g.fraction;}
  }
  const estimable=!!result.mixture?.passesChecks&&ref.wt!==null&&ref.share>=.02&&unresolved<1e-8,ratio=(v:number)=>estimable?v/ref.share:null;
  return {name:ref.name,share:ref.share,wt:ratio(wt),indel:ratio(edited),inFrame:ratio(inFrame),outOfFrame:ratio(outOfFrame),unresolvedFraction:ref.share?unresolved/ref.share:0,components:[...bins].sort((a,b)=>a[0]-b[0]).map(([size,fraction])=>({size,fraction:ref.share?fraction/ref.share:0}))};
 });
 const targets=joint.targets.map(target=>{
  if(target.edited===null)return {...target};
  const ref=refs.find(r=>r.name===target.reference);if(!ref||ref.share<.02)return {...target,edited:null};let edited=0;
  for(const g of joint.groups.filter(g=>g.reference===target.reference)){
   const alternatives=g.alternatives.length?g.alternatives:[{reference:g.reference!,events:g.events}],states=alternatives.map(a=>a.reference===target.reference?a.events.some(e=>within(e,target.cut)):null);
   if(states.some(s=>s===null||s!==states[0]))return {...target,edited:null};if(states[0])edited+=g.fraction;
  }
  return {...target,edited:edited/ref.share};
 });
 return {window,targets,references:refs};
}
