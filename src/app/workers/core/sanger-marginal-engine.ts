import {fitMarginalPoint} from './sanger-marginal-fit';
import {buildMarginalModel} from './sanger-marginal-model';
import {summarizeMarginalPoint} from './sanger-marginal-point-report';
import {marginalSessions} from './sanger-marginal-state';
export {marginalModelDiagnostics,marginalSessions,marginalInputEvent} from './sanger-marginal-state';
export type {MarginalEdge,MarginalInternalGroup,MarginalMeasurement,MarginalTargetState,MarginalSession} from './sanger-marginal-state';
import type {SignalPathInput} from './sanger-path-common';
import type {SangerAnalysisResult} from '../../models/sanger-analysis';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import {marginalViolation,type LinearTerm} from './sanger-marginal-lp';
import {plotTrace} from './sanger-analysis';
const evaluate=(terms:LinearTerm[],point:number[])=>terms.reduce((n,[i,w])=>n+w*point[i],0);
/** One common physical fit for every reference/target/count. Default reporting
 * reads its point only. No whole-model extrema or genotype enumeration runs
 * before returning a result. The complete model is retained for selected detail.
 */
export function analyzeMarginalSanger(input:SignalPathInput,progress:(s:string)=>void,optimizer:SangerLinearOptimizer):SangerAnalysisResult{
 const timings:{stage:string;seconds:number}[]=[],started=performance.now();let last=started,stage='Registering measured fluorescence';
 const tick=(next:string)=>{const now=performance.now();if(now>last)timings.push({stage,seconds:(now-last)/1000});last=now;stage=next;progress(next);};
 const model=buildMarginalModel(input,tick);
 const {refs,report,observation,calibration,quality,shares,noise,analysisEnd,regionalResponse,initialization,robustDelta,allMeasurements,firstMeasured,heldOutError,heldOutVariance,heldOutRows}=model;
 const fit=fitMarginalPoint(model,optimizer,tick,progress),{point,representativeSelection}=fit;
 let squaredError=0,signalVariance=0,signalRows=0;
 let groups:ReturnType<typeof summarizeMarginalPoint>;
 try{
  if(marginalViolation(model.program,point)>1e-7)throw Error('Shared fit failed the original physical model certificate.');
  for(const m of allMeasurements){const den=evaluate(m.total,point);for(let b=0;b<4;b++){const prediction=evaluate(m.forms[b],point)/den,y=m.observed[b];squaredError+=(y-prediction)**2;signalVariance+=(y-.25)**2;signalRows++;}}
  tick('Preparing fitted representatives and base proportions');
  groups=summarizeMarginalPoint(model,point,input.settings.indelThreshold??2);
 }finally{fit.solver.dispose();}
 const targets=model.drafts;
 const annotations=refs.map((r,h)=>({name:r.name,reference:r.input.sequence,cut:r.input.cut,cuts:r.input.cuts,direction:r.reverse?'reverse' as const:'forward' as const,share:shares[h],components:[],marginalGroups:groups.filter(g=>g.reference===r.name),estimable:true})),first=groups.filter(g=>g.target===1),components=Array.from(new Set(first.map(g=>g.size))).sort((a,b)=>a-b).map(size=>({size,fraction:first.filter(g=>g.size===size).reduce((n,g)=>n+g.fraction*shares[refs.findIndex(r=>r.name===g.reference)],0)})),netZero=components.find(c=>c.size===0)?.fraction||0,wtLower=first.filter(g=>g.size===0).reduce((n,g)=>n+(g.unmodifiedEstimate??0)*shares[refs.findIndex(r=>r.name===g.reference)],0),inFrame=components.filter(c=>c.size!==0&&c.size%3===0).reduce((n,c)=>n+c.fraction,0),outOfFrame=components.filter(c=>c.size%3!==0).reduce((n,c)=>n+c.fraction,0);
 const rmse=Math.sqrt(squaredError/Math.max(1,signalRows)),cvRmse=Math.sqrt(heldOutError/Math.max(1,heldOutRows)),passesChecks=calibration.validationR2>=.8&&cvRmse<=Math.max(.005,3*noise);annotations.forEach(a=>a.estimable=passesChecks);
 const result:SangerAnalysisResult={method:'pooled',referenceName:refs.length===1?refs[0].name:'Mixed amplicons',sampleName:input.sample.name,controlName:input.control.name,direction:'forward',cutPosition:Math.min(...refs.flatMap(r=>r.cuts)),alignmentIdentity:0,alignment:[],variants:[],synthetic:!!input.synthetic,plots:[plotTrace(input.control,refs[0].cuts[0]-observation.cr.trim,input.settings.window),plotTrace(input.sample,refs[0].cuts[0]-observation.sr.trim,input.settings.window)],annotations,marginalGroups:groups,exactAlternatives:[],warnings:['Rows show fitted representatives and A/C/G/T point estimates, not unique molecular sequences. Length-family abundance is not representative haplotype abundance. Feasible ranges and sequence alternatives are calculated only for a selected group.'],diagnostics:{representativeSelection,optimization:[],referenceSignal:{frame:'control-normalized',robustDelta,controlShares:shares.slice(),fittedShares:shares.map((share,h)=>share*regionalResponse(h,initialization.center)/shares.reduce((sum,share,k)=>sum+share*regionalResponse(k,initialization.center),0)),responseScale:shares.map((_,h)=>regionalResponse(h,initialization.center)),drift:{center:initialization.center,span:initialization.span,parameters:initialization.parameters,trainingOnly:true}},representativeRefit:[],algorithmVersion:'sanger-marginals-v118-background-profile',signalModel:'measured-control',requestedMaxIndel:input.settings.maxIndel,pipeline:{version:'v118-background-profile',observation:'measured-control',thresholdStage:'postprocess'},initialization:{status:'complete',candidates:groups.length},sampleCalibration:calibration,measuredValidation:{scope:'model-supported-spectrum',controlStart:quality.start,controlEnd:quality.end,fitStart:firstMeasured,fitEnd:analysisEnd,qualityAvailable:quality.qualityAvailable,noiseRmse:noise,noiseObservable:'integrated-photon-area',contrastNoiseRmse:Math.sqrt(observation.contrastNoise),predictionLimit:Math.max(.005,3*noise),crossValidatedRmse:cvRmse,crossValidatedR2:heldOutVariance?1-heldOutError/heldOutVariance:1,nullR2:0,meanPredictionGain:0,predictionGainSE:0,unresolvedFraction:0,editingSensitivity:[1-netZero,1-netZero]}},mixture:{components,wt:wtLower,indel:1-netZero,inFrame,outOfFrame,r2:signalVariance?1-squaredError/signalVariance:1,rmse,upstreamError:noise,bases:allMeasurements.length,fitStart:firstMeasured,fitEnd:analysisEnd,converged:true,passesChecks},pooled:{preflight:report,controlR2:calibration.validationR2,sampleTrim:observation.sr.trim,controlTrim:observation.cr.trim,amplicons:refs.map((r,h)=>({name:r.name,share:shares[h],wt:null,indel:null,inFrame:null,outOfFrame:null,components:groups.filter(g=>g.reference===r.name&&g.target===1).map(g=>({size:g.size,fraction:g.fraction}))}))}};
 timings.push({stage,seconds:(performance.now()-last)/1000});result.stageTimings=timings;marginalSessions.set(result,{refs,shares,targets,optimizer,pages:new Map()});return result;
}
