import type{SangerLinearOptimizer}from'../workers/core/sanger-discrete-phase';
import type{SangerAnalysisResult}from'../models/sanger-analysis';
import{parseSangerTrace}from'./sanger-trace';
import{writeSynthTrace}from'./sanger-synthetic-benchmark';
import{targetEvents}from'./sanger-target-profile';
import{sangerEditedSequence}from'./sanger-sequence-groups';
type Part={sequence:string;reference:number;states:number[];weight:number};
type Row={ids:number[];lower:number;upper:number};
/** Post-inference only. Accept only independently re-rendered native files,
 * identical byte for byte. Search failure does not prove impossibility. */
export function certifyNativeSangerPair(control:ArrayBuffer,edit:ArrayBuffer,result:SangerAnalysisResult,optimizer:SangerLinearOptimizer,audit?:(stage:string,detail?:unknown)=>void,fractionTolerancePp=2){
 const annotations=result.annotations;if(!annotations?.length||annotations.some(a=>a.direction!=='forward'||a.components.some(c=>c.unresolved)))return null;
 const c=parseSangerTrace(control,'control'),e=parseSangerTrace(edit,'edit');if(c.peaks.length!==e.peaks.length)return null;
 const midpointHeights=(trace:typeof c)=>{const rows:number[][]=[];for(let p=0;p<trace.peaks.length;p++){if(p&&trace.peaks[p]-trace.peaks[p-1]!==10)return null;const row=[...'ACGT'].map((b,k)=>trace.channels[b as 'A'][trace.peaks[p]-5]-2-(p?rows[p-1][k]:0));if(row.some(v=>!Number.isInteger(v)||v<0||v>3))return null;rows.push(row);}if(rows.at(-1)!.some((v,k)=>trace.channels['ACGT'[k]as'A'][trace.peaks.at(-1)!+5]-2!==v))return null;return rows;};
 const controlMidpoints=midpointHeights(c),editMidpoints=midpointHeights(e);if(!controlMidpoints||!editMidpoints)return null;
 const refs=annotations.map((a,h)=>{const offset=result.pooled?.preflight.amplicons[h]?.primerEnd??20,axes=(a.cuts||[a.cut]).map((_,t)=>{const unique=new Map<string,{events:any[];fraction:number}>();for(const g of a.components){const events=targetEvents(g.events||[],t),key=sangerEditedSequence(a.reference,events),old=unique.get(key);if(old)old.fraction+=g.fraction;else unique.set(key,{events,fraction:g.fraction});}return[...unique.values()];});return{a,offset,axes};}),parts:Part[]=[];if(refs.some(r=>r.axes.some(axis=>axis.some(g=>!Number.isFinite(g.fraction)||g.fraction<0))))return null;
 // First verify the supplied coherent whole-molecule mixture directly.
 // Recombining all target axes is unnecessary when the actual fitted point
 // already reproduces both files, and used to impose a hidden scaling limit.
 const pointParts:Part[]=annotations.flatMap((a,h)=>a.components.flatMap(g=>g.fraction>0?[{reference:h,states:refs[h].axes.map((axis,t)=>{const seq=sangerEditedSequence(a.reference,targetEvents(g.events||[],t));return axis.findIndex(x=>sangerEditedSequence(a.reference,x.events)===seq);}),sequence:sangerEditedSequence(a.reference,g.events||[]).slice(refs[h].offset),weight:a.share*g.fraction}]:[]));
 const directEqual=(a:ArrayBuffer,b:ArrayBuffer)=>{const x=new Uint8Array(a),y=new Uint8Array(b);return x.length===y.length&&x.every((v,i)=>v===y[i]);};
 const pointTotal=pointParts.reduce((n,p)=>n+p.weight,0);
 if(pointParts.length&&pointParts.every(p=>Number.isFinite(p.weight)&&p.weight>=0&&p.states.every(s=>s>=0))&&Math.abs(pointTotal-1)<1e-8){
  const edited=writeSynthTrace(pointParts.map(p=>({sequence:p.sequence,weight:p.weight})),e.peaks.length,'edit'),wt=writeSynthTrace(refs.map((r,h)=>({sequence:r.a.reference.slice(r.offset),weight:pointParts.reduce((n,p)=>n+(p.reference===h?p.weight:0),0)})),c.peaks.length,'control');
  if(directEqual(edited,edit)&&directEqual(wt,control))return{kind:'exact-native-ab1-pair' as const,sweeps:0,maxDeltaPp:0,parts:pointParts,bothFilesIdentical:true};
 }
 const combinations=refs.reduce((n,r)=>n+r.axes.reduce((p,axis)=>p*axis.length,1),0);if(!Number.isSafeInteger(combinations)||combinations>100000||combinations*Math.max(c.peaks.length,e.peaks.length)*4>20_000_000){audit?.('certificate work budget',{combinations,cycles:Math.max(c.peaks.length,e.peaks.length)});return null;} // explicit certificate work budget, not an allele-count bound
 for(const [h,r]of refs.entries()){
  const visit=(t:number,states:number[],events:any[],mass:number)=>{if(t===r.axes.length){parts.push({reference:h,states,sequence:sangerEditedSequence(r.a.reference,events).slice(r.offset),weight:r.a.share*mass});return;}r.axes[t].forEach((g,j)=>visit(t+1,[...states,j],[...events,...g.events],mass*g.fraction));};visit(0,[],[],1);
 }
 if(!parts.length)return null;const rowMap=new Map<string,Row>(),calls:{a:number[];b:number[];margin:number}[]=[];let impossible=false,reason="";
 const add=(ids:number[],lower:number,upper:number)=>{lower=Math.max(0,lower);upper=Math.min(1,upper);if(lower>upper+1e-15){impossible=true;reason="interval "+lower+" "+upper;return;}if(!ids.length){if(lower>0){impossible=true;reason="empty "+lower;}return;}if(ids.length===parts.length){if(lower>1||upper<1){impossible=true;reason='total';}return;}const key=ids.join(','),old=rowMap.get(key);if(old){old.lower=Math.max(old.lower,lower);old.upper=Math.min(old.upper,upper);if(old.lower>old.upper){impossible=true;reason="intersection "+old.lower+" "+old.upper;}}else rowMap.set(key,{ids,lower,upper});};
 for(const [trace,wt]of [[c,true],[e,false]]as const)for(let p=0;p<trace.peaks.length;p++){
  const groups=[...'ACGT'].map(b=>parts.flatMap((part,i)=>(wt?refs[part.reference].a.reference.slice(refs[part.reference].offset):part.sequence)[p]===b?[i]:[]));
  for(let b=0;b<4;b++){let lo=0,hi=1;for(let d=-4;d<=4;d++){const amplitude=1200*Math.exp(-d*d/4),v=trace.channels[(['A','C','G','T']as const)[b]][trace.peaks[p]+d]-2;lo=Math.max(lo,(v-.5+1e-8)/amplitude);hi=Math.min(hi,(v+.5-1e-8)/amplitude);}const midpoint=(wt?controlMidpoints:editMidpoints)[p][b],midAmplitude=1200*Math.exp(-25/4);lo=Math.max(lo,(midpoint-.5+1e-8)/midAmplitude);hi=Math.min(hi,(midpoint+.5-1e-8)/midAmplitude);add(groups[b],lo,hi);}
  const called='ACGT'.indexOf(trace.sequence[p]);if(called<0)return null;for(let b=0;b<4;b++)if(b!==called)calls.push({a:groups[called],b:groups[b],margin:b<called?1e-10:0});
 }
 const fractionRows=fractionTolerancePp>=100?[]:refs.flatMap((r,h)=>r.axes.flatMap((axis,t)=>axis.map((g,j)=>({ids:parts.flatMap((part,i)=>part.reference===h&&part.states[t]===j?[i]:[]),referenceIds:parts.flatMap((part,i)=>part.reference===h?[i]:[]),lower:Math.max(0,g.fraction-fractionTolerancePp/100+1e-9),upper:Math.min(1,g.fraction+fractionTolerancePp/100-1e-9)}))));
 if(impossible){audit?.('incompatible intervals',reason);return null;}const rows=[...rowMap.values()],w=parts.map(p=>p.weight),sum=w.reduce((a,b)=>a+b,0);for(let i=0;i<w.length;i++)w[i]/=sum;
 const equal=(a:ArrayBuffer,b:ArrayBuffer)=>{const x=new Uint8Array(a),y=new Uint8Array(b);return x.length===y.length&&x.every((v,i)=>v===y[i]);};
 const verify=(sweeps:number)=>{if(w.some(v=>!Number.isFinite(v)||v<0)||Math.abs(w.reduce((a,b)=>a+b,0)-1)>1e-8)return null;const masses=refs.map((_,h)=>parts.reduce((s,p,i)=>s+(p.reference===h?w[i]:0),0));let maxDelta=0;for(const [h,r]of refs.entries())for(let t=0;t<r.axes.length;t++)r.axes[t].forEach((g,j)=>{const mass=parts.reduce((s,p,i)=>s+(p.reference===h&&p.states[t]===j?w[i]:0),0);maxDelta=Math.max(maxDelta,Math.abs(mass/masses[h]-g.fraction)*100);});if(maxDelta>fractionTolerancePp+1e-7){audit?.('conditional fraction bound',{sweeps,maxDelta});return null;}
  const sample=writeSynthTrace(parts.map((p,i)=>({sequence:p.sequence,weight:w[i]})),e.peaks.length,'edit');if(!equal(sample,edit)){const a=new Uint8Array(sample),b=new Uint8Array(edit),at=a.findIndex((v,i)=>v!==b[i]);audit?.('sample byte mismatch',{sweeps,at,actual:a[at],expected:b[at],lengths:[a.length,b.length]});return null;}const controlOut=writeSynthTrace(refs.map((r,h)=>({sequence:r.a.reference.slice(r.offset),weight:masses[h]})),c.peaks.length,'control');if(!equal(controlOut,control))return null;return{kind:'exact-native-ab1-pair' as const,sweeps,maxDeltaPp:maxDelta,parts:parts.map((p,i)=>({...p,weight:w[i]})),bothFilesIdentical:true};
 };
 const first=verify(0);if(first)return first;
 const lines=['Minimize',' obj: 0 x0','Subject To'];let row=0;
 const expression=(ids:number[])=>ids.map(i=>'1000 x'+i).join(' + ');
 for(const r of rows){const expr=expression(r.ids);if(r.lower>0)lines.push(' r'+row+++'l: '+expr+' >= '+(r.lower*1000));if(r.upper<1)lines.push(' r'+row+++'u: '+expr+' <= '+(r.upper*1000));}
 for(const f of fractionRows){const left=f.ids.map(i=>'+ 1000 x'+i).join(' ');lines.push(' f'+row+++'l: '+left+' '+f.referenceIds.map(i=>'- '+(1000*f.lower)+' x'+i).join(' ')+' >= 0',' f'+row+++'u: '+left+' '+f.referenceIds.map(i=>'- '+(1000*f.upper)+' x'+i).join(' ')+' <= 0');}
 const seenCalls=new Set<string>();for(const call of calls){const key=call.a.join(',')+'/'+call.b.join(',')+'/'+call.margin;if(seenCalls.has(key))continue;seenCalls.add(key);const expr=call.a.map(i=>'+ 1000 x'+i).concat(call.b.map(i=>'- 1000 x'+i)).join(' ');if(expr)lines.push(' c'+row+++': '+expr+' >= '+(call.margin*1000));}
 lines.push(' total: '+parts.map((_,i)=>'x'+i).join(' + ')+' = 1','Bounds',...parts.map((_,i)=>' 0 <= x'+i+' <= 1'),'End');
 const answer=optimizer.solve(lines.join('\n'),{output_flag:false,primal_feasibility_tolerance:1e-10,dual_feasibility_tolerance:1e-10,time_limit:10});if(answer.Status!=='Optimal'&&!/time/i.test(answer.Status)){audit?.('solver status',answer.Status);return null;}w.fill(0);for(let i=0;i<w.length;i++)w[i]=answer.Columns['x'+i]?.Primal||0;
 return verify(1);
}
import type{BenchmarkTruth}from'./sanger-benchmark';
/** A benchmark credit needs TWO independently verified witnesses: inferred
 * genotypes and the supplied truth, each rendering the same complete control
 * and edited files. A source-only reconstruction is not evidence that an
 * arbitrary truth table is compatible. Both searches happen after inference. */
export function certifyNativeBenchmark(data:BenchmarkTruth,control:ArrayBuffer,edit:ArrayBuffer,result:SangerAnalysisResult,optimizer:SangerLinearOptimizer){
 if(!data.sequenceTruth||!result.annotations?.length)return null;
 const annotations=result.annotations.map(a=>{const rows=data.truth.filter(t=>t.reference===a.name),cuts=a.cuts||[a.cut];if(!rows.length)return null;
  const components=rows.map(t=>{const events=t.events.map(e=>{const start=e.kind==='substitution'?e.selectionStart!:cuts[e.target]+e.offset,word=e.kind==='deletion'?'':e.sequence,end=start+(e.kind==='insertion'?0:e.kind==='deletion'?e.length:(e.before||'').length);return{target:e.target+1,start,end,inserted:word.length,insertedSequence:word};});if(sangerEditedSequence(a.reference,events)!==t.inputSequence)return null;return{events,size:t.inputSequence.length-a.reference.length,fraction:t.percent/100};});
  if(components.some(c=>!c))return null;return{...a,share:rows[0].share/100,components:components as NonNullable<typeof components[number]>[]};
 });
 if(annotations.some(a=>!a))return null;
 const source=certifyNativeSangerPair(control,edit,result,optimizer);if(!source)return null;
 const truth=certifyNativeSangerPair(control,edit,{...result,annotations:annotations as NonNullable<typeof annotations[number]>[]},optimizer);if(!truth)return null;
 return{kind:'exact-native-ab1-equivalence' as const,references:result.annotations.map(a=>a.name),sourceAnnotationsKey:JSON.stringify(result.annotations),truthKey:JSON.stringify(data.truth),source,truth};
}
export type NativeBenchmarkCertificate=NonNullable<ReturnType<typeof certifyNativeBenchmark>>;
