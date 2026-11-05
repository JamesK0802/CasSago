import {targetEvents} from '../../utils/sanger-target-profile';
import {sangerSequenceKey} from '../../utils/sanger-sequence-groups';
import {SangerAnalysisResult,SangerTrace,TRACE_BASES,TracePlot,SangerJointEvent} from '../../models/sanger-analysis';
import {SangerAnnotationRow} from '../../utils/sanger-annotation';
import {materializeGroup} from '../../utils/sanger-truth-groups';
import {reverseComplement} from '../../utils/biology.utils';
import {emissionVector} from '../../workers/core/sanger-emission';
import {pathColumn} from '../../workers/core/sanger-junction-paths';
import {JunctionModel,JunctionReference} from '../../workers/core/sanger-junction-mixture';
import {plotTrace} from '../../workers/core/sanger-analysis';
import {repairSignalTraceCycles} from '../../workers/core/sanger-cycle';
/** Multiply actual measured samples by smoothly interpolated responsibilities.
 * No reference waveform, Gaussian peak or artificial noise is generated. */
export function attributeRawTrace(trace:SangerTrace,responsibilities:number[][]):SangerTrace{
 const channels=Object.fromEntries(TRACE_BASES.map(b=>[b,new Float32Array(trace.channels[b].length)]))as SangerTrace['channels'];
 let cycle=0;
 for(let sample=0;sample<trace.channels.A.length;sample++){
  while(cycle+1<trace.peaks.length&&sample>trace.peaks[cycle+1])cycle++;
  const next=Math.min(cycle+1,trace.peaks.length-1),distance=trace.peaks[next]-trace.peaks[cycle],t=distance?Math.max(0,Math.min(1,(sample-trace.peaks[cycle])/distance)):0;
  TRACE_BASES.forEach((b,j)=>{const a=responsibilities[cycle]?.[j]||0,c=responsibilities[next]?.[j]||0,share=Math.max(0,Math.min(1,a+(c-a)*t));channels[b][sample]=trace.channels[b][sample]*share;});
 }
 const sequence=Array.from(trace.peaks,p=>{const v=TRACE_BASES.map(b=>channels[b][p]||0),peak=Math.max(...v);return peak>0?TRACE_BASES[v.indexOf(peak)]:'N';}).join('');
 return {...trace,name:trace.name+' · raw group contribution',sequence,quality:null,channels};
}
export function rawSeparatedGroupTrace(result:SangerAnalysisResult,sample:SangerTrace|null,index:number,row:SangerAnnotationRow,span:number,targetIndex=0):{plot:TracePlot|null;note:string}{
 const unavailable=(note:string)=>({plot:null,note});
 if(!sample)return unavailable('The original edited AB1/SCF is required for raw-signal separation.');
 const annotation=result.annotations?.[index];if(!annotation||!row.sequenceKey||row.sequenceState==='unresolved')return unavailable('This group has no resolved sequence for raw-signal attribution.');
 const trace=repairSignalTraceCycles(sample).trace,n=trace.peaks.length,lag=result.diagnostics?.sampleCalibration?.lag||0;
 const geometry=(i:number)=>{const a=result.annotations![i],meta=result.pooled?.preflight.amplicons[i],length=a.reference.length,reverse=a.direction==='reverse';let origin:number,trim=result.pooled?.sampleTrim||0;
  if(meta)origin=reverse?length-meta.primerStart+1:meta.primerEnd;
  else{const rawCut=result.plots[1]?.cutIndex,controlCut=result.plots[0]?.cutIndex;if(rawCut===undefined||controlCut===undefined)throw Error('Registered raw-read coordinates are unavailable.');const referenceOrigin=(reverse?length-a.cut:a.cut)-controlCut,offset=rawCut-controlCut,prefix=Math.min(Math.max(0,referenceOrigin),Math.max(0,offset));origin=referenceOrigin-prefix;trim=prefix-offset;}
  const oriented=reverse?reverseComplement(a.reference):a.reference,sequence=origin>=0?oriented.slice(origin):'N'.repeat(-origin)+oriented;
  return {annotation:a,origin,trim,reverse,length,sequence};};
 try{
  const maps=result.annotations!.map((_,i)=>geometry(i)),positions=Array.from({length:n},(_,p)=>p+maps[0].trim),selected=maps[index],isSequenceFit=result.diagnostics?.algorithmVersion?.startsWith('sanger-sequence-groups-v'),theta=result.diagnostics?.contextEmission?.theta||[],denominator=positions.map(()=>[0,0,0,0]);
  const sequenceFor=(i:number,events:SangerJointEvent[])=>{const g=maps[i],full=materializeGroup(g.annotation.reference,events),oriented=g.reverse?reverseComplement(full):full;return g.origin>=0?oriented.slice(g.origin):'N'.repeat(-g.origin)+oriented;};
  const sequenceColumn=(i:number,events:SangerJointEvent[])=>{const seq=sequenceFor(i,events);return positions.flatMap(p=>emissionVector(seq,p,lag,theta));};
  const add=(column:number[],weight:number)=>{positions.forEach((_,p)=>{for(let b=0;b<4;b++)denominator[p][b]+=weight*column[p*4+b];});};
  let numerator:number[];let numeratorWeighted=false;
  if(isSequenceFit){
   numerator=positions.flatMap(()=>[0,0,0,0]);numeratorWeighted=true;
   for(let i=0;i<maps.length;i++)for(const c of maps[i].annotation.components){if(c.fraction<=0)continue;if(c.unresolved)return unavailable('Unresolved fitted molecules prevent reliable channel attribution in this result.');const column=sequenceColumn(i,c.events||[]),weight=maps[i].annotation.share*c.fraction;add(column,weight);const events=(annotation.cuts?.length||1)>1?targetEvents(c.events||[],targetIndex):c.events||[];if(i===index&&sangerSequenceKey(annotation.reference,events)===row.sequenceKey)for(let k=0;k<column.length;k++)numerator[k]+=weight*column[k];}
  }else{
   const diagnostic=result.diagnostics,fit=diagnostic?.junctionFit;if(!fit||!diagnostic?.controlResponse)return unavailable('Registered per-molecule contributions are unavailable for this result.');
   const refs:JunctionReference[]=maps.map(g=>({sequence:g.sequence,cut:(g.reverse?g.length-g.annotation.cut:g.annotation.cut)-g.origin,contextResponse:diagnostic.controlResponse}));
   const modelPath=(m:JunctionModel)=>m.path||(m.size?[{target:1,size:m.size,origins:m.origins?.length?m.origins:[{start:m.start,weight:1,probabilities:m.probabilities||(m.size>0?Array.from({length:m.size},(_,j)=>Array.from('ACGT',b=>/^[ACGT]$/.test(m.inserted[j]||'')?Number(m.inserted[j]===b):.25)):undefined)}]}]:[]);
   for(const model of fit.models){if(model.fraction>0)add(pathColumn(refs[model.reference],{...model,path:modelPath(model as JunctionModel)},positions,lag,[]),model.fraction);}
   const events=(row.representative?.events||[]).map(e=>({target:e.target,start:(selected.reverse?selected.length-e.end:e.start)-selected.origin,removed:e.end-e.start,inserted:selected.reverse?reverseComplement(e.insertedSequence||''):e.insertedSequence||''})).sort((a,b)=>a.start-b.start),path=events.flatMap(e=>[...(e.inserted?[{target:e.target,size:e.inserted.length,origins:[{start:e.start,weight:1,probabilities:Array.from(e.inserted,base=>Array.from('ACGT',b=>Number(base===b)))}]}]:[]),...(e.removed?[{target:e.target,size:-e.removed,origins:[{start:e.start,weight:1}]}]:[])]);
   numerator=pathColumn(refs[index],{reference:index,size:row.size,start:events[0]?.start||0,inserted:'',path},positions,lag,[]);
  }
  const weight=annotation.share*row.fraction,responsibilities=positions.map((_,p)=>TRACE_BASES.map((_,b)=>denominator[p][b]>1e-12?(numeratorWeighted?numerator[p*4+b]:weight*numerator[p*4+b])/denominator[p][b]:0)),attributed=attributeRawTrace(trace,responsibilities);
  const target=annotation.cuts?.[targetIndex]??annotation.cut,orientedCut=selected.reverse?selected.length-target:target,events=(row.representative?.events||[]).map(e=>({start:selected.reverse?selected.length-e.end:e.start,end:selected.reverse?selected.length-e.start:e.end,inserted:e.inserted}));let shift=0,shifted=orientedCut;
  for(const e of events.sort((a,b)=>a.start-b.start)){if(e.start>=orientedCut)break;const removed=e.end-e.start;if(e.inserted===removed)continue;if(e.end<=orientedCut)shift+=e.inserted-removed;else{shifted=e.start;break;}}
  shifted+=shift;
  const cut=Math.round(shifted-selected.origin-selected.trim);if(cut<0||cut>=n-2)return unavailable('The selected target lies outside measured raw cycles.');
  return {plot:plotTrace(attributed,cut,Math.ceil(span/2),Math.floor(span/2)),note:'Measured raw-channel contribution. Shared peaks are apportioned using the fitted mixture; unmodelled channels remain unassigned. This is not an independently measured molecule.'};
 }catch(error){return unavailable(error instanceof Error?error.message:String(error));}
}
