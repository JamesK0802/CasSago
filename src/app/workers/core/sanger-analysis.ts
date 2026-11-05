import {analyzeSingleJunction} from './sanger-single-junction';
import {matchedSingleSignals} from './sanger-single-signals';
import {alignConsensusSequence,consensusVariants} from './sanger-consensus';
import {integratedPeaks} from './sanger-signal';
import { reverseComplement } from '../../utils/biology.utils';
import { SangerAnalysisResult, SangerAnalysisSettings, SangerTrace, SangerVariant, TraceAlignment, TraceAlignmentColumn, TracePlot, TRACE_BASES, TraceBase } from '../../models/sanger-analysis';

/** Smith–Waterman local alignment with a linear gap penalty; no N matches. */
export function alignTraceSequence(reference:string,query:string):{score:number;identity:number;columns:TraceAlignmentColumn[]} {
  if(reference.length*query.length>24000000)throw new Error('Alignment region is too large.');
  const m=query.length,n=reference.length,trace=new Uint8Array((n+1)*(m+1));
  let prev=new Int32Array(m+1),curr=new Int32Array(m+1),best=0,bi=0,bj=0;
  for(let i=1;i<=n;i++) {
    curr[0]=0;
    for(let j=1;j<=m;j++) {
      const equal=reference[i-1]===query[j-1]&&/[ACGT]/.test(query[j-1]);
      const diag=prev[j-1]+(equal?3:-4),up=prev[j]-5,left=curr[j-1]-5;
      const value=Math.max(0,diag,up,left);curr[j]=value;
      trace[i*(m+1)+j]=value===0?0:value===diag?1:value===up?2:3;
      if(value>best){best=value;bi=i;bj=j;}
    }
    [prev,curr]=[curr,prev];
  }
  let i=bi,j=bj;const columns:TraceAlignmentColumn[]=[];
  while(i>0&&j>0) {
    const step=trace[i*(m+1)+j];if(!step)break;
    if(step===1){columns.push({ref:i-1,query:j-1,a:reference[i-1],b:query[j-1]});i--;j--;}
    else if(step===2){columns.push({ref:i-1,query:null,a:reference[i-1],b:'-'});i--;}
    else{columns.push({ref:null,query:j-1,a:'-',b:query[j-1]});j--;}
  }
  columns.reverse();const matches=columns.filter(c=>c.a===c.b&&/[ACGT]/.test(c.a)).length;
  return {score:best,identity:columns.length?matches/columns.length:0,columns};
}
function referenceAlignment(reference:string,cut:number,trace:SangerTrace,direction:SangerAnalysisSettings['direction'],consensus=false):TraceAlignment {
  const regionRadius=Math.min(1700,trace.sequence.length+100);let best:TraceAlignment|null=null;
  for(const dir of ['forward','reverse'] as const) {
    if(direction!=='auto'&&direction!==dir)continue;
    const oriented=dir==='forward'?reference:reverseComplement(reference),orientedCut=dir==='forward'?cut:reference.length-cut;
    const lo=Math.max(0,orientedCut-regionRadius),hi=Math.min(reference.length,orientedCut+regionRadius);
    const aligned=(consensus?alignConsensusSequence:alignTraceSequence)(oriented.slice(lo,hi),trace.sequence);
    const columns=aligned.columns.map(c=>({...c,ref:c.ref===null?null:c.ref+lo}));
    const coordinates=columns.filter(c=>c.ref!==null).map(c=>c.ref!);
    const item={...aligned,columns,direction:dir,cut:orientedCut,referenceStart:Math.min(...coordinates),referenceEnd:Math.max(...coordinates)+1};
    if(!best||item.score>best.score)best=item;
  }
  if(!best||best.identity<.8||best.columns.length<60||best.referenceStart>best.cut-15||best.referenceEnd<best.cut+15)throw new Error('The trace does not align reliably across this target. Check reference, target and read direction.');
  return best;
}
function traceCut(alignment:TraceAlignment):number {
  const column=alignment.columns.find(c=>c.ref!==null&&c.ref>=alignment.cut&&c.query!==null);
  if(!column)throw new Error('The cut site is not covered by the trace.');return column.query!;
}
export function normalizedPeaks(trace:SangerTrace):number[][] {
  return Array.from(trace.peaks,peak=>{
    const values=TRACE_BASES.map(b=>{
      let sum=0,count=0;for(let at=Math.max(0,peak-2);at<=Math.min(trace.channels[b].length-1,peak+2);at++){sum+=trace.channels[b][at];count++;}
      return sum/Math.max(1,count);
    });
    const total=values.reduce((a,b)=>a+b,0);return values.map(v=>total>0?v/total:0);
  });
}
/** Cyclic coordinate descent NNLS; normalized output fractions are re-scored. */
export function fitTraceMixture(columns:number[][],observed:number[]):{weights:number[];r2:number;rmse:number;converged:boolean} {
  const p=columns.length,n=observed.length,gram=Array.from({length:p},()=>new Float64Array(p)),rhs=new Float64Array(p);
  for(let a=0;a<p;a++) {
    rhs[a]=columns[a].reduce((v,x,i)=>v+x*observed[i],0);
    for(let b=0;b<=a;b++){let dot=0;for(let i=0;i<n;i++)dot+=columns[a][i]*columns[b][i];gram[a][b]=dot;gram[b][a]=dot;}
  }
  for(let a=0;a<p;a++)for(let b=0;b<a;b++)if(gram[a][a]+gram[b][b]-2*gram[a][b]<1e-8)throw new Error('Indel models are indistinguishable in this sequence window. Extend the usable trace range or reduce the indel range.');
  const weights=new Float64Array(p);let converged=false;
  for(let step=0;step<10000;step++) {
    let delta=0;
    for(let a=0;a<p;a++) {
      let residual=rhs[a];for(let b=0;b<p;b++)if(a!==b)residual-=gram[a][b]*weights[b];
      const value=gram[a][a]>0?Math.max(0,residual/gram[a][a]):0;
      delta=Math.max(delta,Math.abs(value-weights[a]));weights[a]=value;
    }
    if(delta<1e-9){converged=true;break;}
  }
  const total=weights.reduce((a,b)=>a+b,0);if(total<=0)throw new Error('The trace mixture model has no usable signal.');
  const fractions=Array.from(weights,w=>w/total),mean=observed.reduce((a,b)=>a+b,0)/n;
  let error=0,variance=0;
  for(let i=0;i<n;i++){let pred=0;for(let a=0;a<p;a++)pred+=fractions[a]*columns[a][i];error+=(observed[i]-pred)**2;variance+=(observed[i]-mean)**2;}
  return {weights:fractions,r2:variance>1e-12?1-error/variance:0,rmse:Math.sqrt(error/n),converged};
}
export function plotTrace(trace:SangerTrace,cutIndex:number,maxAfter:number,maxBefore=30):TracePlot {
  const before=Math.min(maxBefore,cutIndex),after=Math.min(maxAfter,trace.sequence.length-cutIndex-2),start=-before,end=after;
  const width=(end-start)*13+32,paths={} as Record<TraceBase,string>;
  let max=1;
  const samples:{x:number;values:number[]}[]=[];
  for(let offset=start;offset<end;offset++) {
    const i=cutIndex+offset;
    for(let k=0;k<6;k++) {
      const at=trace.peaks[i]+(trace.peaks[i+1]-trace.peaks[i])*k/6,lo=Math.floor(at),fraction=at-lo;
      const values=TRACE_BASES.map(b=>{const a=trace.channels[b][lo]||0,c=trace.channels[b][lo+1]||0;return a+(c-a)*fraction;});
      max=Math.max(max,...values);samples.push({x:16+(offset-start+k/6)*13,values});
    }
  }
  TRACE_BASES.forEach((b,i)=>paths[b]=samples.map((p,j)=>`${j?'L':'M'}${p.x.toFixed(2)},${(124-98*p.values[i]/max).toFixed(2)}`).join(' '));
  return {cutIndex,name:trace.name,start,end,width,cutX:16+(before-.5)*13,paths,letters:Array.from({length:Math.max(0,end-start)},(_,j)=>({x:16+j*13,base:trace.sequence[cutIndex+start+j],offset:start+j}))};
}
export interface SangerAnalysisInput {referenceName:string;reference:string;cut:number;sample:SangerTrace;control?:SangerTrace;settings:SangerAnalysisSettings;synthetic?:boolean}
export function analyzeSanger(input:SangerAnalysisInput,progress:(stage:string)=>void=()=>{}):SangerAnalysisResult {
  const {reference,cut,sample,control,settings:s}=input;
  if(s.joinNearby!==undefined&&(!Number.isInteger(s.joinNearby)||s.joinNearby<0||s.joinNearby>10))throw new Error('Join nearby changes must be 0–10 matched bases.');
  if(s.signalModel!==undefined&&!['classic','calibrated'].includes(s.signalModel))throw new Error('Invalid signal processing mode.');
  if(!reference||reference.length>2000000||/[^ACGTN]/.test(reference)||!Number.isInteger(cut)||cut<1||cut>=reference.length)throw new Error('Provide one A/C/G/T/N reference and one cut position inside it.');
  if(!['mixture','consensus'].includes(s.method)||!['auto','forward','reverse'].includes(s.direction))throw new Error('Invalid Sanger analysis method or direction.');
  if(!Number.isInteger(s.window)||s.window<60||s.window>500||!Number.isInteger(s.maxIndel)||s.maxIndel<1||s.maxIndel>30||!Number.isInteger(s.controlPhred)||s.controlPhred<0||s.controlPhred>60)throw new Error('Use a 60–500 bp window, 1–30 bp maximum indel and Phred 0–60.');
  if(s.method==='mixture'&&(!Number.isInteger(s.readStart)||!Number.isInteger(s.readEnd)||s.readStart<0||s.readEnd>(s.analysisPolicy==='automatic'?3160:1000)||s.readEnd<=s.readStart))throw new Error('Use a measured read range within 0–1,000 bp.');
  const extended=(s.analysisPolicy==='automatic'||!!s.largeDeletions)&&s.method==='mixture';
  const maxDeletion=extended?(s.maxDeletion??200):s.maxIndel;
  const maxInsertion=extended?(s.maxInsertion??s.maxIndel):s.maxIndel;
  if(extended&&(s.signalModel==='classic'||!Number.isInteger(maxDeletion)||maxDeletion<1||maxDeletion>500||!Number.isInteger(maxInsertion)||maxInsertion<1||maxInsertion>30))throw new Error('Large-deletion decomposition requires calibrated signal, deletion limit 1–500 and insertion limit 1–30 bp.');
  const originRadius=s.editStartRadius??5;
  if(!Number.isInteger(originRadius)||originRadius<0||originRadius>100)throw new Error('Indel origin radius must be 0–100 bp.');
  const warnings:string[]=[];
  progress('Aligning trace to the target');
  const anchor=s.method==='mixture'?control:sample;
  if(!anchor)throw new Error('Trace decomposition requires an unedited control AB1 from the same locus and sequencing primer.');
  const aligned=referenceAlignment(reference,cut,anchor,s.direction,s.method==='consensus'),cutIndex=traceCut(aligned);
  const base:SangerAnalysisResult={method:s.method,referenceName:input.referenceName,sampleName:sample.name,controlName:s.method==='mixture'?control!.name:undefined,
    direction:aligned.direction,cutPosition:cut,alignmentIdentity:aligned.identity,alignment:aligned.columns,variants:[],plots:[],warnings,synthetic:!!input.synthetic};
  if(s.method==='consensus') {
    base.variants=consensusVariants(aligned,reference.length,reference.length*2,s.joinNearby??3);
    base.plots=[plotTrace(sample,cutIndex,Math.min(s.window,120))];
    warnings.push('Base calls describe a consensus. Mixed peaks are not decomposed and no editing frequency is estimated.');
    if(sample.quality===null)warnings.push('This file has no per-base quality values.');
    return base;
  }
  if(aligned.identity<.95)throw new Error('Control-to-reference identity is below 95%. Check that this is an unedited control from the intended locus.');
  progress('Matching control and sample upstream of the cut');
  let anchorStart=Math.max(10,cutIndex-140),anchorEnd=cutIndex-12;
  if(anchorEnd-anchorStart<60)throw new Error('At least 60 bp of usable sequence upstream of the cut is required.');
  let upstream=alignTraceSequence(control!.sequence.slice(anchorStart,anchorEnd),sample.sequence.slice(0,Math.min(sample.sequence.length,cutIndex+180)));
  if(s.signalModel!=='classic')for(let start=10;start+100<=cutIndex-12;start+=10){
    const rows=Array.from({length:100},(_,j)=>start+j);
    if(rows.filter(p=>/[ACGT]/.test(control!.sequence[p])&&(!control!.quality||control!.quality[p]>=s.controlPhred)).length<80)continue;
    const candidate=alignTraceSequence(control!.sequence.slice(start,start+100),sample.sequence.slice(0,Math.min(sample.sequence.length,cutIndex+180))),shifts=new Map<number,number>();
    for(const c of candidate.columns)if(c.ref!==null&&c.query!==null){const shift=c.query-(start+c.ref);shifts.set(shift,(shifts.get(shift)||0)+1);}
    if(candidate.identity>=.9&&Math.max(0,...shifts.values())>=Math.max(60,candidate.columns.length*.9)){anchorStart=start;anchorEnd=start+100;upstream=candidate;break;}
  }
  if(upstream.identity<.9||upstream.columns.filter(c=>c.ref!==null&&c.query!==null).length<60)throw new Error('Control and sample do not share a reliable upstream region. Use the same sequencing primer and direction.');
  const offsets=new Map<number,number>();
  for(const c of upstream.columns)if(c.ref!==null&&c.query!==null){const offset=c.query-(anchorStart+c.ref);offsets.set(offset,(offsets.get(offset)||0)+1);}
  const [offset,count]=[...offsets.entries()].sort((a,b)=>b[1]-a[1])[0];
  if(count<upstream.columns.length*.9)throw new Error('Upstream alignment contains conflicting shifts. Trace decomposition requires an intact shared upstream region.');
  const enhanced=s.signalModel!=='classic',sampleCut=cutIndex+offset;
  let upstreamRows=Array.from({length:anchorEnd-anchorStart},(_,j)=>anchorStart+j).filter(i=>i+offset>=0&&i+offset<sample.peaks.length&&(enhanced||!control!.quality||control!.quality[i]>=s.controlPhred));
  if(upstreamRows.length<40)throw new Error('Too few measured upstream control bases for registration.');
  let matched=enhanced?matchedSingleSignals(control!,sample,upstreamRows,offset):undefined;
  if(matched){
    // Select an intact WT calibration interval from measured concordance. This
    // choice is independent of the requested deletion/insertion search limits;
    // all these positions still participate in the final junction likelihood.
    const windows:[[number,number],[number,number],[number,number]]=[[10,Math.min(110,cutIndex-12)],[Math.max(10,cutIndex-140),cutIndex-12],[Math.max(10,cutIndex-240),Math.max(10,cutIndex-140)]];
    let quality=matched.calibration.validationR2-matched.upstreamError**2;
    for(const [lo,hi] of windows){const rows=Array.from({length:Math.max(0,hi-lo)},(_,j)=>lo+j).filter(i=>i+offset>=0&&i+offset<sample.peaks.length&&(enhanced||!control!.quality||control!.quality[i]>=s.controlPhred));if(rows.length<60)continue;const candidate=matchedSingleSignals(control!,sample,rows,offset),score=candidate.calibration.validationR2-candidate.upstreamError**2;if(score>quality){matched=candidate;upstreamRows=rows;quality=score;}}
  }
  if(matched){
    const origins=new Map<number,number>();for(const c of aligned.columns)if(c.ref!==null&&c.query!==null)origins.set(c.ref-c.query,(origins.get(c.ref-c.query)||0)+1);
    const referenceOrigin=[...origins].sort((a,b)=>b[1]-a[1])[0][0],oriented=aligned.direction==='forward'?reference:reverseComplement(reference);
    return analyzeSingleJunction({base,reference:oriented,referenceOrigin,cutIndex,control:control!,sample,settings:s,offset,upstreamRows,matched,maxDeletion,maxInsertion,plots:[plotTrace(control!,cutIndex,Math.min(s.window,120)),plotTrace(sample,sampleCut,Math.min(s.window,120))]});
  }
  // Classic changes signal extraction only; it uses the same physical junction
  // hypotheses and complete measured read axis as calibrated processing.
  const wt=normalizedPeaks(control!),raw=normalizedPeaks(sample),controlSignal={...integratedPeaks(control!),vectors:wt},sampleSignal={...integratedPeaks(sample),vectors:raw};
  const upstreamError=upstreamRows.reduce((sum,i)=>sum+wt[i].reduce((v,x,b)=>v+Math.abs(x-raw[i+offset][b]),0)/2,0)/upstreamRows.length;
  const calibration={gains:[1,1,1,1],lag:0,used:false,trainingBases:0,validationBases:0,rawValidationR2:1,validationR2:1,weights:[1]};
  const identity={method:'area' as const,controlSignal,sampleSignal,wt,raw,edited:raw,calibration,upstreamError,templateUsable:wt.map(()=>true)};
  const origins=new Map<number,number>();for(const c of aligned.columns)if(c.ref!==null&&c.query!==null)origins.set(c.ref-c.query,(origins.get(c.ref-c.query)||0)+1);
  const referenceOrigin=[...origins].sort((a,b)=>b[1]-a[1])[0][0],oriented=aligned.direction==='forward'?reference:reverseComplement(reference);
  return analyzeSingleJunction({base,reference:oriented,referenceOrigin,cutIndex,control:control!,sample,settings:s,offset,upstreamRows,matched:identity,maxDeletion,maxInsertion,plots:[plotTrace(control!,cutIndex,Math.min(s.window,120)),plotTrace(sample,sampleCut,Math.min(s.window,120))]});
}
