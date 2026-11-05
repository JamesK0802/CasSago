import {SangerTrace,SANGER_ANALYSIS_DEFAULTS} from '../models/sanger-analysis';
import {alignTraceSequence} from '../workers/core/sanger-analysis';
import {registerTrace} from '../workers/core/sanger-pooled';
import {exactBindingSites} from '../workers/core/sanger-primer-designer';
import {reverseComplement} from './biology.utils';
export interface PreviewReference {name:string;sequence:string;cuts:number[];targetNames?:string[]}
export interface PreviewSite {name:string;cut:number}
/** Preview registration only: no edit fitting, ground truth or mixture changes. */
export function tracePreviewSites(refs:PreviewReference[],primer:string,trace:SangerTrace,knownTrim?:number):PreviewSite[]{
 if(!refs.length||!trace.sequence.length)return [];
 const label=(r:PreviewReference,t:number)=>r.name.replace(/_\d+-\d+$/,'')+' · '+(r.targetNames?.[t]||'T'+(t+1));
 if(/^[ACGT]{12,40}$/.test(primer)){
  try{const prepared=refs.map(r=>{const hits=exactBindingSites(primer,{...r,cuts:r.cuts});if(hits.length!==1)throw Error('Ambiguous primer site.');const hit=hits[0],sequence=hit.direction==='forward'?r.sequence:reverseComplement(r.sequence),origin=hit.direction==='forward'?hit.end:r.sequence.length-hit.start,cuts=r.cuts.map(c=>(hit.direction==='forward'?c:r.sequence.length-c)-origin);return {name:r.name,sequence:sequence.slice(origin),cut:Math.min(...cuts),cuts};});
   const trim=knownTrim??registerTrace(trace,prepared,{...SANGER_ANALYSIS_DEFAULTS,readStart:0}).trim;
   return prepared.flatMap((r,i)=>r.cuts.map((cut,t)=>({name:label(refs[i],t),cut:cut-trim}))).filter(s=>s.cut>=0&&s.cut<trace.sequence.length-2);
  }catch{if(refs.length>1)return [];}
 }
 if(refs.length!==1)return [];
 const ref=refs[0],sites:PreviewSite[]=[];
 for(const [t,cut] of ref.cuts.entries()){
  const choices:{score:number;cut:number}[]=[];
  for(const reverse of [false,true]){const seq=reverse?reverseComplement(ref.sequence):ref.sequence,c=reverse?ref.sequence.length-cut:cut,anchor=seq.slice(Math.max(0,c-20),c+20);if(anchor.length<30||seq.indexOf(anchor)!==seq.lastIndexOf(anchor))continue;const radius=Math.min(900,trace.sequence.length+100),lo=Math.max(0,c-radius),hi=Math.min(seq.length,c+radius);
   try{const alignment=alignTraceSequence(seq.slice(lo,hi),trace.sequence),mapped=alignment.columns.filter(x=>x.ref!==null&&x.query!==null),before=mapped.filter(x=>x.ref!+lo<c),after=mapped.filter(x=>x.ref!+lo>=c);if(alignment.identity<.8||mapped.length<60||before.length<15||after.length<15)continue;choices.push({score:alignment.score,cut:after[0].query!});}catch{continue;}
  }
  choices.sort((a,b)=>b.score-a.score);if(choices[0]&&(!choices[1]||choices[0].score>choices[1].score))sites.push({name:label(ref,t),cut:choices[0].cut});
 }
 return sites;
}
