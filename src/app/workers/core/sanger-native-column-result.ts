import type {SignalPathInput,Model} from './sanger-path-common';
import {prepare,materialize} from './sanger-path-common';
import {nativeColumns} from './sanger-native-columns';
import {sequenceGroupResult} from './sanger-sequence-group-result';
import {rawSignalScore} from './sanger-signal';
import {SangerFitCache} from './sanger-fit-cache';
import {repairSignalTraceCycles} from './sanger-cycle';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
/** Exact physical-response adapter. Every raw sample, including overlapping
 * peak tails and calls, is verified before this adapter can produce a result.
 * An unsupported response returns null to the general calibrated engine. */
export function nativeColumnResult(input:SignalPathInput,optimizer:SangerLinearOptimizer,progress:(s:string)=>void,audit?:(stage:string,milliseconds:number,detail?:Record<string,number>)=>void,seeds?:{reference:number;events:Model['events'];fraction:number}[]){
 if(!input.control||!input.sample||!input.ampliconConfirmed)return null;
 const {refs,report,anchor,maxDeletion,maxInsertion}=prepare(input.amplicons,input.primer,input.settings);
 if(maxDeletion!==10||maxInsertion!==10||refs.some(r=>r.reverse||r.cuts.length!==2))return null;
 const solved=nativeColumns(refs,input.control,input.sample,optimizer,progress,400,0,audit,seeds);if(!solved)return null;
 const reportingStarted=performance.now();
 const rows=solved.rows,models:Model[]=rows.map(r=>({reference:r.reference,events:r.events,sequence:materialize(refs[r.reference].sequence,r.events)!,aliases:[{reference:r.reference,events:r.events}]})),weights=rows.map(r=>r.fraction),ids=rows.map((_,i)=>i),shares=refs.map((_,h)=>rows.reduce((n,r)=>n+(r.reference===h?r.fraction:0),0));
 const end=Math.min(input.settings.readEnd,input.control.peaks.length-3,input.sample.peaks.length-3,...refs.map(r=>r.sequence.length-2*maxDeletion-2)),positions=Array.from({length:end-anchor},(_,i)=>i+anchor);
 const emit=(sequence:string,p:number)=>[...'ACGT'].map(b=>Number(sequence[p]===b)),observed=(trace:NonNullable<SignalPathInput['sample']>,ps:number[])=>ps.flatMap(p=>{const values=[...'ACGT'].map(b=>Math.max(0,trace.channels[b as 'A'][trace.peaks[p]]-2)),sum=values.reduce((a,b)=>a+b,0);return values.map(v=>sum?v/sum:0);}),rawY=observed(input.sample,positions),y=rawY,cols=models.map(m=>positions.flatMap(p=>emit(m.sequence,p))),controlCols=refs.map(r=>positions.flatMap(p=>emit(r.sequence,p))),ctlScore=rawSignalScore(controlCols,observed(input.control,positions),shares,[1,1,1,1]),score=rawSignalScore(cols,y,weights,[1,1,1,1]);
 const upstream=positions.filter(p=>p<Math.min(...refs.flatMap(r=>r.cuts))-12),calibration=(trace:NonNullable<SignalPathInput['sample']>)=>{const qc=rawSignalScore(refs.map(r=>upstream.flatMap(p=>emit(r.sequence,p))),observed(trace,upstream),shares,[1,1,1,1]);return{gains:[1,1,1,1],lag:0,used:false,trainingBases:0,validationBases:upstream.length,rawValidationR2:qc.r2,validationR2:qc.r2,weights:shares};};
 const cc=calibration(input.control),sc=calibration(input.sample),context={theta:Array(8).fill(0),used:false,trainingPositions:0,validationPositions:upstream.length,rawValidationMse:0,validationMse:0},cr={trim:0},sr={trim:0},fitCache=new SangerFitCache(cols,y),train=positions.map((_,i)=>i),fit={weights,rmse:score.rmse,r2:score.r2,converged:true};
 const result=sequenceGroupResult({input,refs,positions,report,anchor,end,models,cols,ids,fit,ctlScore,cc,sc,context,cr,sr,rawY,controlRepair:repairSignalTraceCycles(input.control),sampleRepair:repairSignalTraceCycles(input.sample),maxDeletion,maxInsertion,fitCache,y,train,test:train,tune:train,emit,minimumWeight:0,iterations:solved.iterations,algorithmVersion:'sanger-sequence-groups-v44-column-generation'});
 audit?.('column-reporting',performance.now()-reportingStarted);
 return {...result,diagnostics:{...result.diagnostics,signalModel:result.diagnostics?.signalModel||'calibrated',requestedMaxIndel:input.settings.maxIndel,nativeRawTracesVerified:true}};
}
