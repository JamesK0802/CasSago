import {solveJunctions,JunctionSolution} from './sanger-junction-solution';
import {isWildType,pathKey} from './sanger-junction-paths';
import {fitControlTransfer} from './sanger-control-transfer';
import {refineJunctionSupport,rejectNonpredictiveJunctions} from './sanger-junction-support';
import {fitControlResponse,controlContextEmission,ControlResponseProfile} from './sanger-control-response';
import {junctionPredictionEvidence} from './sanger-junction-evidence';
import {fitTransitionResponse,transitionEmission,TransitionResponse} from './sanger-transition-response';
import {fitJunctionMixture,junctionBase,junctionAnnotationComponents,JunctionModel,JunctionFit,aliasedDeletionMass} from './sanger-junction-mixture';
import {emissionGroupSupport,sequenceYield,emissionVector,fitEmissionMixture} from './sanger-emission';
import {integratedPeaks,signalWeights,simplexFit,calibrateSignal,correctGain,phaseVector,SignalCalibration,PeakSignal,gramCondition,rawSignalScore} from './sanger-signal';
import { PooledAmplicon, PooledPreflight, SangerAnalysisResult, SangerAnalysisSettings, SangerTrace, TRACE_BASES } from '../../models/sanger-analysis';
import { reverseComplement } from '../../utils/biology.utils';
import { exactBindingSites } from './sanger-primer-designer';
import { fitTraceMixture, normalizedPeaks, plotTrace } from './sanger-analysis';
interface Prepared {contextResponse?:ControlResponseProfile;name:string;sequence:string;cut:number;cuts?:number[];transitionResponse?:number[][];baseYields?:number[]}
interface Model {reference:number;size:number;sizes?:number[];path?:JunctionModel['path'];start?:number;inserted?:string;sequence:string;values:number[]}
interface ModelSet {refs:Prepared[];report:PooledPreflight;models:Model[]}
const failed=(reason:string):PooledPreflight=>({passed:false,reason,positions:[],modelCount:0,rank:0,condition:null,minimumDifferences:0,minimumCutOffset:0,amplicons:[]});
/** Eigenvalues of the normalized symmetric Gram matrix (cyclic Jacobi).
 * Geometry gate only; thresholds are prototype design policies, not validated limits. */
export {gramCondition} from './sanger-signal';

function prepare(amplicons:PooledAmplicon[],primer:string,s:SangerAnalysisSettings,measuredCoverage=false):ModelSet {
  const report=failed('');
  if(s.signalModel!==undefined&&!['classic','calibrated'].includes(s.signalModel))throw new Error('Invalid signal processing mode.');
  if(s.checkWiderRange!==undefined&&typeof s.checkWiderRange!=='boolean')throw new Error('Invalid wider-range check setting.');
  if(amplicons.length<1||amplicons.length>4)throw new Error('Use 2–4 submitted amplicons for mixed-amplicon analysis.');
  if(new Set(amplicons.map(r=>r.name)).size!==amplicons.length||amplicons.some(r=>!r.name.trim()))throw new Error('Amplicon names must be unique and non-empty.');
  if(!/^[ACGT]{12,40}$/.test(primer))throw new Error('Enter the one sequencing primer actually submitted to Sanger (12–40 nt, 5′→3′).');
  if(!Number.isInteger(s.editStartRadius??5)||(s.editStartRadius??5)<0||(s.editStartRadius??5)>100)throw new Error('Indel origin radius must be 0–100 bp.');
  if(!Number.isInteger(s.window)||s.window<60||s.window>500||!Number.isInteger(s.maxIndel)||s.maxIndel<1||s.maxIndel>(s.analysisPolicy==='automatic'||amplicons.length===1?30:15)||!Number.isInteger(s.readStart)||!Number.isInteger(s.readEnd)||s.readStart<0||s.readEnd>(s.analysisPolicy==='automatic'?3160:1000)||s.readEnd<=s.readStart)throw new Error('Use a 60–500 bp target window, indel range ±1–15 bp and a usable read range within 0–1,000 bp.');
  const maxDeletion=s.analysisPolicy==='automatic'?s.maxDeletion!:s.maxIndel,maxInsertion=s.analysisPolicy==='automatic'?s.maxInsertion!:s.maxIndel;
  if(!Number.isInteger(maxDeletion)||maxDeletion<1||maxDeletion>500||!Number.isInteger(maxInsertion)||maxInsertion<1||maxInsertion>30)throw Error('Invalid deletion/insertion limits.');
  const refs=amplicons.map(r=>{
    if(!/^[ACGT]+$/.test(r.sequence)||r.sequence.length>20000||!Number.isInteger(r.cut)||r.cut<1||r.cut>=r.sequence.length)throw new Error(`${r.name}: use the full submitted amplicon (A/C/G/T, ≤20 kb) and a cut inside it.`);
    const hits=exactBindingSites(primer,{name:r.name,sequence:r.sequence,cuts:r.cuts||[r.cut]});
    if(hits.length!==1)throw new Error(`${r.name}: the submitted sequencing primer must have one exact binding site across both strands.`);
    const hit=hits[0],oriented=hit.direction==='forward'?r.sequence:reverseComplement(r.sequence),end=hit.direction==='forward'?hit.end:r.sequence.length-hit.start;
    const cuts=(r.cuts||[r.cut]).map(c=>(hit.direction==='forward'?c:r.sequence.length-c)-end),cut=Math.min(...cuts),sequence=oriented.slice(end),lo=cut-Math.floor(s.window/2),hi=Math.max(...cuts)+Math.ceil(s.window/2);
    if(cuts.length>2||cuts.some(c=>!Number.isInteger(c)||c<0||c>=sequence.length)||new Set(cuts).size!==cuts.length)throw Error(r.name+': enter one or two distinct target positions.');
    if(cut<s.readStart+72||(s.analysisPolicy!=='automatic'&&Math.max(...cuts)+30>Math.min(s.readEnd,sequence.length)))throw new Error(`${r.name}: the usable trace must cover the target, at least 60 upstream bases and 30 downstream bases.`);
    report.amplicons.push({name:r.name,direction:hit.direction,primerStart:hit.start+1,primerEnd:hit.end,cutOffset:cut,windowStart:lo,windowEnd:hi,targets:cuts.map(c=>({cutOffset:c,windowStart:c-Math.floor(s.window/2),windowEnd:c+Math.ceil(s.window/2)}))});
    return {name:r.name,sequence,cut,cuts};
  });
  const start=s.readStart,end=measuredCoverage?s.readEnd:Math.min(s.readEnd,...refs.map(r=>r.sequence.length));
  // Fixed measured read axis. Calibrated fits retain all junction observations;
  // the configured indel range changes hypotheses, never the observed cycles.
  const positions=Array.from({length:end-start},(_,j)=>start+j);
  const models:Model[]=[];
  refs.forEach((ref,h)=>{for(let size=-maxDeletion;size<=maxInsertion;size++){
    if(s.analysisPolicy==='automatic'&&size!==0)continue;
    const sequence=positions.map(p=>ref.sequence[p<ref.cut?p:p-size]||'N').join('');
    models.push({reference:h,size,sequence,values:Array.from(sequence).flatMap(base=>TRACE_BASES.map(b=>Number(base===b)))});
  }});
  report.positions=positions;report.modelCount=s.analysisPolicy==='automatic'?refs.length*(maxDeletion+maxInsertion+1):models.length;
  if(positions.length<(s.analysisPolicy==='automatic'?60:Math.max(80,models.length)))throw new Error(s.analysisPolicy==='automatic'?'Too few reference-covered sample cycles. Review trace/reference coverage.':'Too few usable positions for this number of amplicon × indel models. Increase the window or reduce the indel range.');
  const checked=s.analysisPolicy==='automatic'?models.filter(m=>m.size===0):models;
  let minimumDifferences=Infinity;
  for(let a=0;a<checked.length;a++)for(let b=0;b<a;b++){
    let diff=0;for(let j=0;j<positions.length;j++)if(checked[a].sequence[j]!==checked[b].sequence[j])diff++;
    minimumDifferences=Math.min(minimumDifferences,diff);
  }
  report.minimumDifferences=minimumDifferences;
  report.minimumCutOffset=Math.min(...refs.flatMap((r,i)=>refs.slice(i+1).map(t=>Math.abs(r.cut-t.cut))));
  if(minimumDifferences<6){report.reason=s.analysisPolicy==='automatic'?'Submitted WT references have few distinguishing positions. Reference attribution may be unstable.':'At least two amplicon/indel models are identical or too similar on the sequencing-primer axis. Change the window, reduce the indel range or redesign the primer/start positions.';return {refs,report,models};}
  let geometryColumns=checked.map(m=>m.values),requiredRank=checked.length;
  if(s.analysisPolicy==='automatic'&&checked.length>1){
   const pivot=checked.length-1;
   geometryColumns=checked.slice(0,-1).map((m,h)=>positions.flatMap((_,row)=>{
    const v=Array.from({length:4},(_,b)=>checked.reduce((sum,c)=>sum+c.values[4*row+b]/checked.length,0)),total=v.reduce((sum,x)=>sum+x,0),a=m.values.slice(4*row,4*row+4),z=checked[pivot].values.slice(4*row,4*row+4),at=a.reduce((sum,x)=>sum+x,0),zt=z.reduce((sum,x)=>sum+x,0);
    return a.map((x,b)=>total?(x-z[b]-(v[b]/total)*(at-zt))/total:0);
   }));requiredRank=checked.length-1;
  }
  const geometry=gramCondition(geometryColumns);report.rank=geometry.rank+(s.analysisPolicy==='automatic'&&checked.length>1?1:0);report.condition=geometry.condition;
  if(geometry.rank<requiredRank||geometry.condition===null||geometry.condition>100){report.reason='The combined amplicon/indel model is rank-deficient or poorly conditioned after fluorescence normalization. This trace cannot be separated reliably under the current references/window/indel range.';return {refs,report,models};}
  report.passed=true;report.reason='Reference models are distinguishable under the current settings. Experimental trace fit must also pass.';
  return {refs,report,models};
}
export function checkPooledSeparation(amplicons:PooledAmplicon[],primer:string,s:SangerAnalysisSettings):PooledPreflight {
  try{return prepare(amplicons,primer,s).report;}catch(e){return failed(e instanceof Error?e.message:String(e));}
}
export class PooledSeparationError extends Error {constructor(public report:PooledPreflight){super(report.reason);}}
/** Register the first called base to the shared primer extension axis. */
export function registerTrace(trace:SangerTrace,refs:Prepared[],s:SangerAnalysisSettings,signals?:number[][],allowQCWarnings=false):{trim:number;r2:number;warnings?:string[]} {
  const observed=signals||normalizedPeaks(trace),firstCut=Math.min(...refs.map(r=>r.cut));const options:{trim:number;r2:number}[]=[];
  for(let trim=0;trim<=160;trim++){
    const lo=Math.max(s.readStart,trim+10),hi=Math.min(firstCut-12,trim+trace.sequence.length-2);
    if(hi-lo<60)continue;
    const positions=Array.from({length:hi-lo},(_,j)=>lo+j).filter(p=>observed[p-trim]?.some(v=>v>0));if(positions.length<60)continue;
    // Identical upstream WT templates can be merged for registration only.
    const sequences=[...new Set(refs.map(r=>positions.map(p=>r.sequence[p]).join('')))];
    const columns=sequences.map(seq=>Array.from(seq).flatMap(base=>TRACE_BASES.map(b=>Number(base===b))));
    const fit=simplexFit(columns,positions.flatMap(p=>observed[p-trim]));options.push({trim,r2:fit.r2});
  }
  options.sort((a,b)=>b.r2-a.r2);const best=options[0];
  if(!best)throw new Error(`${trace.name}: insufficient upstream coverage for registration.`);
  const warnings:string[]=[];
  if(best.r2<(signals ? .55 : .8)&&!allowQCWarnings)throw new Error(`${trace.name}: the upstream trace cannot be registered to the submitted primer/amplicons (fit R² <0.80 or insufficient coverage).`);
  if(options[1]&&best.r2-options[1].r2<.02&&!allowQCWarnings)throw new Error(`${trace.name}: multiple sequencing start positions explain the upstream trace. Trim/primer registration is ambiguous.`);
  if(best.r2<(signals ? .55 : .8))warnings.push(`${trace.name}: weak upstream registration; the best fitted start is shown.`);
  if(options[1]&&best.r2-options[1].r2<.02)warnings.push(`${trace.name}: sequencing start is ambiguous; the best fitted start is shown.`);
  return {...best,warnings};
}
export interface PooledAnalysisInput {amplicons:PooledAmplicon[];primer:string;ampliconConfirmed:boolean;sample:SangerTrace;control:SangerTrace;settings:SangerAnalysisSettings;synthetic?:boolean}
export function preparePooledSanger(input:PooledAnalysisInput,progress:(stage:string)=>void=()=>{}) {
  if(!input.ampliconConfirmed)throw new PooledSeparationError(failed('Confirm that every reference is the actual full amplicon submitted to Sanger. Genomic/gene reference sequences are not valid here.'));
  progress('Checking amplicon/indel separability');let built:ModelSet;
  try{built=prepare(input.amplicons,input.primer,input.settings);}catch(e){throw new PooledSeparationError(failed(e instanceof Error?e.message:String(e)));}
  let {refs,report,models}=built;const analysisWarnings:string[]=[];
  if(!input.control||!input.sample)throw new Error('Mixed-amplicon analysis requires one unedited pooled control and one pooled sample, using the same sequencing primer.');
  progress('Registering pooled traces to the sequencing primer');
  const enhanced=input.settings.signalModel!=='classic';
  let controlSignal=enhanced?integratedPeaks(input.control,[0,0,0,0]):null,sampleSignal=enhanced?integratedPeaks(input.sample,[0,0,0,0]):null;
  const controlRegistration=registerTrace(input.control,refs,input.settings,controlSignal?.vectors,true),sampleRegistration=registerTrace(input.sample,refs,input.settings,sampleSignal?.vectors,true);
  analysisWarnings.push(...(controlRegistration.warnings||[]),...(sampleRegistration.warnings||[]));
  const firstCut=Math.min(...refs.map(r=>r.cut));
  if(enhanced){controlSignal=integratedPeaks(input.control,[0,0,0,0]);sampleSignal=integratedPeaks(input.sample,[0,0,0,0]);}
  const controlRaw=controlSignal?.vectors||normalizedPeaks(input.control),sampleRaw=sampleSignal?.vectors||normalizedPeaks(input.sample);
  const margin=input.settings.analysisPolicy==='automatic'?0:10,tail=input.settings.analysisPolicy==='automatic'?0:2;
  const actualStart=input.settings.analysisPolicy==='automatic'?Math.max(input.settings.readStart,sampleRegistration.trim):Math.max(input.settings.readStart,controlRegistration.trim+margin,sampleRegistration.trim+margin),actualEnd=input.settings.analysisPolicy==='automatic'?Math.min(input.settings.readEnd,sampleRegistration.trim+sampleRaw.length):Math.min(input.settings.readEnd,controlRegistration.trim+controlRaw.length-tail,sampleRegistration.trim+sampleRaw.length-tail);
  if(input.settings.analysisPolicy==='automatic'||actualStart>input.settings.readStart||report.positions.at(-1)!>=actualEnd){const covered=prepare(input.amplicons,input.primer,{...input.settings,readStart:actualStart,readEnd:actualEnd},input.settings.analysisPolicy==='automatic');({refs,report,models}=covered);}
  const coverageWarnings=refs.flatMap((r,h)=>(r.cuts||[r.cut]).flatMap((cut,t)=>cut+30>actualEnd?[`${r.name} target ${t+1}: this read has fewer than 30 measured downstream bases; estimates use the available signal.`]:[]));
  if(!report.passed)analysisWarnings.push(`Reference QC warning: ${report.reason} Fitted reference attribution may not be unique.`);

  let response:TransitionResponse|undefined,baseYields:number[]|undefined,rawResponse:ReturnType<typeof fitControlResponse>|undefined;
  const calibrate=(raw:number[][],trim:number,theta?:number[],prior?:number[])=>{
    const observed=Array.from({length:raw.length+trim},(_,p)=>raw[p-trim]||[0,0,0,0]);
    const positions=Array.from({length:firstCut-12},(_,p)=>p).filter(p=>p>=Math.max(input.settings.readStart,trim+10)&&p-trim<raw.length-2);
    const templates=refs.map(r=>Array.from(r.sequence,(base,p)=>{if(rawResponse){const value=controlContextEmission(r.sequence,p,0,rawResponse),total=refs.reduce((sum,ref,h)=>sum+(prior?.[h]??1/refs.length)*controlContextEmission(ref.sequence,p,0,rawResponse!).reduce((s,v)=>s+v,0),0);return value.map(v=>v/Math.max(1e-9,total));}if(response?.used){const value=transitionEmission(r.sequence,p,0,theta||Array(8).fill(0),response.table,baseYields||[1,1,1,1]),total=refs.reduce((sum,ref,h)=>sum+(prior?.[h]??1/refs.length)*transitionEmission(ref.sequence,p,0,theta||Array(8).fill(0),response!.table,baseYields||[1,1,1,1]).reduce((s,v)=>s+v,0),0);return value.map(v=>v/Math.max(1e-9,total));}const yieldValue=theta?sequenceYield(r.sequence,p,theta):1,total=theta?refs.reduce((s,ref,h)=>s+(prior?.[h]||0)*sequenceYield(ref.sequence,p,theta),0):1;return [...TRACE_BASES].map(b=>Number(base===b)*yieldValue/Math.max(1e-9,total));}));
    return calibrateSignal(observed,positions,templates,'fine');
  };
  let controlCalibration:SignalCalibration|undefined,sampleCalibration:SignalCalibration|undefined,transfer:ReturnType<typeof fitControlTransfer>|undefined;
  if(enhanced){progress('Calibrating trace gains and phase using held-out upstream WT signal');controlCalibration=calibrate(controlRaw,controlRegistration.trim);sampleCalibration=calibrate(sampleRaw,sampleRegistration.trim);
    if(controlCalibration.validationR2<.8||sampleCalibration.validationR2<.8)analysisWarnings.push('Upstream held-out calibration R² is below0.80; fitted estimates may be unreliable.');}
  let wt=controlRaw.map(v=>controlCalibration?correctGain(v,controlCalibration.gains):v),sample=sampleRaw.map(v=>sampleCalibration?correctGain(v,sampleCalibration.gains):v);
  const emissionRows=Array.from({length:controlRaw.length+controlRegistration.trim},(_,p)=>p).filter(p=>p>=Math.max(input.settings.readStart,controlRegistration.trim+10));
  const theta=Array(8).fill(0),fitter=enhanced?fitEmissionMixture:simplexFit;
  if(!enhanced){
    const rawControl=Array.from({length:controlRaw.length+controlRegistration.trim},(_,p)=>controlRaw[p-controlRegistration.trim]||[0,0,0,0]);
    response=fitTransitionResponse(rawControl,refs.map(r=>r.sequence),emissionRows,0);
    if(response.used){baseYields=Array(16).fill(1);refs=refs.map(ref=>({...ref,transitionResponse:response!.table,baseYields}));}
  }
  if(enhanced){
   const areas=Array.from({length:controlRaw.length+controlRegistration.trim},(_,p)=>controlSignal!.areas?.[p-controlRegistration.trim]||[0,0,0,0]);
   rawResponse=fitControlResponse(refs.map(r=>r.sequence),emissionRows,areas,controlCalibration!.weights);response=undefined;refs=refs.map(r=>({...r,transitionResponse:undefined,contextResponse:rawResponse!}));
   let cp=rawResponse.shares,sp=sampleCalibration!.weights;for(let iteration=0;iteration<3;iteration++){controlCalibration=calibrate(controlRaw,controlRegistration.trim,theta,cp);sampleCalibration=calibrate(sampleRaw,sampleRegistration.trim,theta,sp);cp=controlCalibration.weights;sp=sampleCalibration.weights;}
   const sampleAreas=Array.from({length:sampleRaw.length+sampleRegistration.trim},(_,p)=>sampleSignal!.areas![p-sampleRegistration.trim]||[0,0,0,0]),training=Array.from({length:firstCut-12},(_,p)=>p).filter(p=>p>=Math.max(input.settings.readStart,sampleRegistration.trim+10)&&p-sampleRegistration.trim<sampleRaw.length-2);
   transfer=fitControlTransfer(refs.map(r=>r.sequence),training,sampleAreas,rawResponse,sampleCalibration!.lag,sampleCalibration!.gains,sp);
   if(transfer.used)sampleCalibration={...sampleCalibration!,gains:transfer.gains,weights:transfer.shares,used:true,validationR2:transfer.colourR2,trainingBases:transfer.trainingBases,validationBases:transfer.validationBases};
   wt=controlRaw.map(v=>correctGain(v,controlCalibration!.gains));sample=sampleRaw.map(v=>correctGain(v,sampleCalibration!.gains));
  }
  const positions=report.positions;
  if(positions.some(p=>p-sampleRegistration.trim<0||p-sampleRegistration.trim>=sample.length||!sample[p-sampleRegistration.trim].some(v=>v>0)||(input.settings.analysisPolicy!=='automatic'&&(p-controlRegistration.trim<0||p-controlRegistration.trim>=wt.length||!wt[p-controlRegistration.trim].some(v=>v>0)))))throw new Error(input.settings.analysisPolicy==='automatic'?'Some reference-axis sample cycles have no measured signal. Review trace/reference coverage.':'The actual traces do not cover every checked target-window position. Shorten the window or resequence.');
  const columnsFor=(items:Model[],rows:number[],lag:number)=>items.map(m=>{const ref=refs[m.reference];return rows.flatMap(p=>{if(rawResponse&&m.start===undefined){const start=ref.cut,sequence=ref.sequence.slice(0,start)+(m.size>0?'N'.repeat(m.size):'')+ref.sequence.slice(start+Math.max(0,-m.size));return controlContextEmission(sequence,p,lag,rawResponse);}if(m.start!==undefined)return m.values.slice(4*positions.indexOf(p),4*positions.indexOf(p)+4);if(response?.used){const start=ref.cut,seq=ref.sequence.slice(0,start)+(m.size>0?'N'.repeat(m.size):'')+ref.sequence.slice(start+Math.max(0,-m.size));return transitionEmission(seq,p,lag,theta,response.table,baseYields!);}if(m.start===undefined){const q=p<ref.cut?p:p-m.size;return emissionVector(ref.sequence,q,lag,theta);}const event={reference:m.reference,size:m.size,start:m.start,inserted:m.inserted||''};const local=Array.from({length:4},(_,i)=>junctionBase(ref.sequence,event,p-2+i)).join('');return emissionVector(local,2,lag,theta);});});
  const rawScore=(columns:number[][],rows:number[],weights:number[],raw:number[][],trim:number,gains:number[])=>{
    return rawSignalScore(columns,rows.flatMap(p=>raw[p-trim]),weights,gains);
  };
  const reliability=enhanced?positions.map(p=>wt[p-controlRegistration.trim]?.some(v=>v>0)?signalWeights([controlSignal!,sampleSignal!],[[p-controlRegistration.trim],[p-sampleRegistration.trim]])[0]:signalWeights([sampleSignal!],[[p-sampleRegistration.trim]])[0]):undefined;
  const controlPositions=positions.filter(p=>wt[p-controlRegistration.trim]?.some(v=>v>0)),controlReliability=controlPositions.map(p=>reliability?.[positions.indexOf(p)]??1);
  progress('Checking unedited pooled-control fit');
  const controlColumns=enhanced?columnsFor(models.filter(m=>m.size===0),controlPositions,controlCalibration!.lag):models.filter(m=>m.size===0).map(m=>m.values);
  let controlFit=enhanced?fitter(controlColumns,controlPositions.flatMap(p=>wt[p-controlRegistration.trim]),controlReliability):simplexFit(controlColumns,positions.flatMap(p=>wt[p-controlRegistration.trim]));
  if(enhanced)controlFit={...controlFit,...rawScore(controlColumns,controlPositions,controlFit.weights,controlRaw,controlRegistration.trim,controlCalibration!.gains)};
  if(controlFit.r2<.85)analysisWarnings.push('The pooled control is poorly explained by the submitted amplicons (R² <0.85). Check primer, products, control and trace quality.');
  progress('Fitting calibrated, nonnegative amplicon/indel signals');
  const observed=positions.flatMap(p=>sample[p-sampleRegistration.trim]);
  const fitInput={references:refs,positions,observed,maxDeletion:input.settings.analysisPolicy==='automatic'?input.settings.maxDeletion!:input.settings.maxIndel,maxInsertion:input.settings.analysisPolicy==='automatic'?input.settings.maxInsertion!:input.settings.maxIndel,radius:input.settings.editStartRadius??5,lag:sampleCalibration?.lag||0,theta,positionWeights:reliability,fit:fitter};
  const templateError=(junction:JunctionFit)=>{
    const controlAt=new Map(controlPositions.map((p,i)=>[p,i])),templateScale=Math.max(1,...refs.map((_,h)=>junction.models.reduce((sum,m,i)=>sum+(m.reference===h?junction.weights[i]:0),0)/Math.max(1e-8,controlFit.weights[h])));
    return enhanced?positions.flatMap(p=>{
      const at=controlAt.get(p);if(at===undefined)return [0,0,0,0];
      const value=Array.from({length:4},(_,b)=>controlColumns.reduce((sum,c,h)=>sum+controlFit.weights[h]*c[4*at+b],0)),total=value.reduce((sum,v)=>sum+v,0);
      return value.map((v,b)=>templateScale*(wt[p-controlRegistration.trim][b]-(total?v/total:0)));
    }):undefined;
  };
  const finish=(solution:JunctionSolution):SangerAnalysisResult=>{
  const {junction,unprunedWeights,support,referenceEvidence,predictive,predictionRemoved}=solution;
  models=junction.models.map((m,i)=>({...m,sequence:'',values:junction.columns[i]}));
  const columns=junction.columns;
  let fit:{weights:number[];r2:number;rmse:number;converged:boolean}=junction;
  if(enhanced)fit={...fit,...rawScore(columns,positions,fit.weights,sampleRaw,sampleRegistration.trim,sampleCalibration!.gains)};
  const residuals=enhanced?rawScore(columns,positions,fit.weights,sampleRaw,sampleRegistration.trim,sampleCalibration!.gains).positionRmse:[],discordantPositions=residuals.filter(v=>v>.20).length,maxPositionRmse=residuals.length?Math.max(...residuals):undefined;
  let referenceStability:{jointTVpp:number;maximumConditionalTVpp:number;maximumShareDifferencePp:number}|undefined;
  let splitStabilityPp:number|undefined,rangeSuggestion:{maxIndel:number;validationR2Gain:number;outsideFraction:number}|undefined;
  if(enhanced){
    const stabilityInput=solution.stabilityInput,stabilityColumns=solution.stabilityColumns;
    const halves=[0,1].map(parity=>{const ids=stabilityInput.positions.flatMap((_,i)=>i%2===parity?[i]:[]);return fitter(stabilityColumns.map(c=>ids.flatMap(i=>c.slice(i*4,i*4+4))),ids.flatMap(i=>stabilityInput.observed.slice(i*4,i*4+4)),ids.map(i=>stabilityInput.positionWeights?.[i]??1)).weights;});
    splitStabilityPp=Math.max(...refs.map((_,h)=>{const edit=halves.map(w=>{const share=models.reduce((sum,m,i)=>sum+(m.reference===h?w[i]:0),0),zero=models.reduce((sum,m,i)=>sum+(m.reference===h&&isWildType(m)?w[i]:0),0);return share>=.02?1-zero/share:null;});return edit.some(v=>v===null)?0:100*Math.abs(edit[0]!-edit[1]!);}));
    const share=(weights:number[],h:number)=>models.reduce((sum,m,i)=>sum+(m.reference===h?weights[i]:0),0);
    const jointTVpp=50*halves[0].reduce((sum,w,i)=>sum+Math.abs(w-halves[1][i]),0);
    const maximumShareDifferencePp=Math.max(...refs.map((_,h)=>100*Math.abs(share(halves[0],h)-share(halves[1],h))));
    const maximumConditionalTVpp=Math.max(...refs.map((_,h)=>{const a=share(halves[0],h),b=share(halves[1],h);if(Math.max(a,b)<.02)return 0;if(Math.min(a,b)<.02)return 100;return 50*models.reduce((sum,m,i)=>sum+(m.reference===h?Math.abs(halves[0][i]/a-halves[1][i]/b):0),0);}));
    referenceStability={jointTVpp,maximumConditionalTVpp,maximumShareDifferencePp};
    splitStabilityPp=Math.max(splitStabilityPp,jointTVpp,maximumConditionalTVpp,maximumShareDifferencePp);
    if(input.settings.analysisPolicy!=='automatic'&&input.settings.checkWiderRange!==false&&input.settings.maxIndel<15){
      const expansions:ModelSet[]=[];
      for(let range=Math.min(15,input.settings.maxIndel+5);range<=15;range=Math.min(15,range+5)){
        try{const expanded=prepare(input.amplicons,input.primer,{...input.settings,maxIndel:range});if(expanded.report.passed)expansions.push(expanded);}catch{}
        if(range===15)break;
      }
      if(expansions.length){const common=positions.filter(p=>expansions.every(e=>e.report.positions.includes(p))),train=common.filter((_,i)=>i%4!==0),valid=common.filter((_,i)=>i%4===0);
        if(train.length>=Math.max(...expansions.map(e=>e.models.length))&&valid.length>=30){
          const y=train.flatMap(p=>sample[p-sampleRegistration.trim]),small=fitter(columnsFor(models,train,sampleCalibration!.lag),y),baseScore=rawScore(columnsFor(models,valid,sampleCalibration!.lag),valid,small.weights,sampleRaw,sampleRegistration.trim,sampleCalibration!.gains),accepted:{maxIndel:number;validationR2Gain:number;outsideFraction:number;rmse:number}[]=[];
          for(const expanded of expansions){const fitted=fitter(columnsFor(expanded.models,train,sampleCalibration!.lag),y),score=rawScore(columnsFor(expanded.models,valid,sampleCalibration!.lag),valid,fitted.weights,sampleRaw,sampleRegistration.trim,sampleCalibration!.gains),outside=expanded.models.reduce((sum,m,i)=>sum+(Math.abs(m.size)>input.settings.maxIndel?fitted.weights[i]:0),0);
            if(score.r2-baseScore.r2>=.01&&score.rmse<baseScore.rmse*.85&&outside>=.03)accepted.push({maxIndel:Math.max(...expanded.models.map(m=>m.size)),validationR2Gain:score.r2-baseScore.r2,outsideFraction:outside,rmse:score.rmse});
          }
          if(accepted.length){const best=Math.min(...accepted.map(a=>a.rmse)),pick=accepted.filter(a=>a.rmse<=best*1.02).sort((a,b)=>a.maxIndel-b.maxIndel)[0];rangeSuggestion={maxIndel:pick.maxIndel,validationR2Gain:pick.validationR2Gain,outsideFraction:pick.outsideFraction};}
        }
      }
    }
  }

  const amplicons=refs.map((ref,h)=>{
    const components=models.map((m,i)=>({...m,fraction:fit.weights[i]})).filter(m=>m.reference===h),share=components.reduce((a,b)=>a+b.fraction,0);
    const zero=components.filter(isWildType).reduce((s,c)=>s+c.fraction,0),inside=components.filter(c=>!isWildType(c)&&c.size%3===0).reduce((a,b)=>a+b.fraction,0),outside=components.filter(c=>c.size%3!==0).reduce((a,b)=>a+b.fraction,0);
    return {name:ref.name,share,wt:share>0?zero/share:null,indel:share>0?1-zero/share:null,inFrame:share>0?inside/share:null,outOfFrame:share>0?outside/share:null,components:[...new Set(components.map(c=>c.size))].map(size=>({size,fraction:share>0?components.filter(c=>c.size===size).reduce((s,c)=>s+c.fraction,0)/share:0}))};
  });
  const sizes=[...new Set(models.map(m=>m.size))].sort((a,b)=>a-b),components=sizes.map(size=>({size,fraction:models.reduce((sum,m,i)=>sum+(m.size===size?fit.weights[i]:0),0)}));
  const zero=models.reduce((sum,m,i)=>sum+(isWildType(m)?fit.weights[i]:0),0),inFrame=models.reduce((sum,m,i)=>sum+(!isWildType(m)&&m.size%3===0?fit.weights[i]:0),0),outOfFrame=components.filter(c=>c.size%3!==0).reduce((a,b)=>a+b.fraction,0);
  const warnings:string[]=[...analysisWarnings,...coverageWarnings];if(discordantPositions)warnings.push(`${discordantPositions} local positions have high unexplained signal. Check for trace artifacts or variants outside the modeled classes.`);if(splitStabilityPp!==undefined&&splitStabilityPp>10)warnings.push('Reference shares or within-reference indel spectra differ by more than 10 pp between alternating-position fits. Review component stability.');if(rangeSuggestion)warnings.push(`A wider ±${rangeSuggestion.maxIndel} bp model improves held-out signal fit. Reanalyse with that range before interpreting the indel spectrum.`);if(fit.r2<.9)warnings.push('Pooled sample R² is below 0.90; component estimates are exploratory.');if(!fit.converged)warnings.push('Pooled trace regression did not converge.');
  let automaticGeometryPass=true;
  if(input.settings.analysisPolicy==='automatic'&&aliasedDeletionMass(refs,junction.models,fit.weights,positions,fitInput.maxDeletion,fitInput.radius)>.01){automaticGeometryPass=false;warnings.push('Some fitted deletions have indistinguishable alternative lengths over the measured DNA span. Estimates are shown; attribution may not be unique.');}
  if(input.settings.analysisPolicy==='automatic'){const active=solution.stabilityColumns.filter((_,i)=>Math.max(fit.weights[i],unprunedWeights[i])>.01),g=gramCondition(active);if(g.rank<active.length||g.condition===null||g.condition>100){automaticGeometryPass=false;warnings.push('Fitted components are poorly conditioned; review the reported attribution.');}const aliased=solution.stabilityColumns.some((a,i)=>Math.max(fit.weights[i],unprunedWeights[i])>.01&&solution.stabilityColumns.some((b,j)=>i!==j&&a.every((v,k)=>Math.abs(v-b[k])<1e-8)));if(aliased){automaticGeometryPass=false;warnings.push('At least one fitted component has an indistinguishable alternative. Attribution may not be unique.');}}
  if(amplicons.some(r=>r.share<.02))warnings.push('At least one amplicon contributes less than 2% of fitted signal; its within-reference percentages may be unstable.');
  if(controlFit.weights.some(w=>w<.02))warnings.push('At least one submitted amplicon has less than 2% fitted signal in the control.');
  const plots=[{trace:input.control,trim:controlRegistration.trim},{trace:input.sample,trim:sampleRegistration.trim}].map(({trace,trim})=>{
    const plot=plotTrace(trace,firstCut-trim,Math.min(120,Math.max(...refs.map(r=>r.cut))-firstCut+input.settings.window/2));
    plot.cuts=refs.map(r=>({name:r.name,x:plot.cutX+(r.cut-firstCut)*13})).filter(c=>c.x<=plot.width);return plot;
  });
  const aliases=solution.aliases.map(a=>({reference:refs[a.reference].name,target:a.target,possibleReallocatedFraction:amplicons[a.reference].share?a.shiftableJointFraction/amplicons[a.reference].share:0}));
  if(aliases.some(a=>a.possibleReallocatedFraction>=.01)){automaticGeometryPass=false;warnings.push('Different physical edit mixtures give the same measured signal for '+aliases.map(a=>a.reference+' T'+a.target).join(', ')+'. The best numerical estimates are shown; target attribution is not unique.');}
  const checksPass=(solution.coveragePass??!coverageWarnings.length)&&automaticGeometryPass&&report.passed&&!analysisWarnings.length&&(sampleCalibration?.validationR2??1)>=.8&&(controlCalibration?.validationR2??1)>=.8&&controlFit.r2>=.85&&fit.r2>=.9&&fit.converged&&controlFit.converged&&(splitStabilityPp??0)<=10&&!rangeSuggestion&&!discordantPositions;
  const multi=refs.some(r=>(r.cuts?.length||1)>1);
  const pathComponents=(h:number)=>junctionAnnotationComponents(junction.models.filter(m=>m.reference===h),fit.weights.filter((_,i)=>junction.models[i].reference===h),(report.amplicons[h].direction==='forward'?input.amplicons[h].cut:input.amplicons[h].sequence.length-input.amplicons[h].cut)-(refs[h].cuts?.[0]??refs[h].cut),input.amplicons[h].sequence.length,report.amplicons[h].direction==='reverse');
  const joint=multi?{attributionAmbiguous:aliases,groups:junction.models.flatMap((m,i)=>{if(fit.weights[i]<=1e-8)return [];const components=pathComponents(m.reference).filter(c=>c.pathSizes?.join('/')===m.sizes?.join('/')),best=components.sort((a,b)=>b.fraction-a.fraction)[0];return [{id:m.reference+':'+pathKey(m),reference:refs[m.reference].name,fraction:fit.weights[i],netSize:m.size,pathSizes:m.sizes,events:best?.events||[],alternatives:[]}];}),targets:refs.flatMap((r,h)=>(input.amplicons[h].cuts||[input.amplicons[h].cut]).map((cut,t)=>({reference:r.name,target:t+1,cut,edited:amplicons[h].share?pathComponents(h).filter(c=>c.events?.some(e=>e.target===t+1||e.start<cut&&e.end>cut)).reduce((sum,c)=>sum+c.fraction,0)/amplicons[h].share:0}))),heldOutR2:predictive?.size?Math.min(...[...predictive.values()].flatMap(e=>e.validationR2||[])):null,selectedModels:fit.weights.filter(v=>v>1e-8).length,candidateModels:junction.candidateCount||junction.models.length,iterations:junction.iterations,maxDeletion:fitInput.maxDeletion,ambiguousFraction:0}:undefined;
  if(multi&&!junction.searchComplete)warnings.push('Joint candidate search reached its computation budget; the best numerical estimates are shown.');
  return {joint,method:'pooled',referenceName:'Mixed amplicons',sampleName:input.sample.name,controlName:input.control.name,direction:'forward',cutPosition:firstCut,alignmentIdentity:0,alignment:[],variants:[],plots,warnings,synthetic:!!input.synthetic,diagnostics:{modeledRange:input.settings.analysisPolicy==='automatic'?{maxDeletion:input.settings.maxDeletion!,maxInsertion:input.settings.maxInsertion!}:undefined,transitionResponse:response,junctionFit:junction?{positions,iterations:junction.iterations,models:junction.models.map((m,i)=>({...m,fraction:fit.weights[i]})),duplicateLengths:junction.duplicateLengths}:undefined,referenceStability,algorithmVersion:enhanced?'sanger-junction-v2':'sanger-classic-v1',signalModel:enhanced?'calibrated':'classic',controlCalibration,sampleCalibration,splitStabilityPp,rangeSuggestion,discordantPositions:enhanced?discordantPositions:undefined,maxPositionRmse,requestedMaxIndel:input.settings.maxIndel,controlResponse:rawResponse,controlTransfer:transfer,supportScreen:enhanced?{predictionRemoved:predictionRemoved.map(i=>({reference:models[i].reference,size:models[i].size,foldGains:predictive!.get(i)!.foldGains})),removed: support!.removed.map(i=>({reference:models[i].reference,size:models[i].size})),threshold:support!.threshold,evidence:(referenceEvidence||[]).map(e=>({...e,reference:models[e.group[0]].reference,size:models[e.group[0]].size})),unpruned:models.map((m,i)=>({reference:m.reference,size:m.size,jointFraction:unprunedWeights[i]}))}:undefined,background:enhanced?{method:'integrated-four-channel-area',control:{levels:controlSignal!.background||[],trainingSamples:controlSignal!.backgroundSamples||[]},sample:{levels:sampleSignal!.background||[],trainingSamples:sampleSignal!.backgroundSamples||[]}}:undefined},annotations:input.amplicons.map((a,i)=>({name:a.name,reference:a.sequence,cut:a.cut,cuts:a.cuts,direction:report.amplicons[i].direction,share:amplicons[i].share,components:junction?junctionAnnotationComponents(junction.models.filter(m=>m.reference===i),fit.weights.filter((_,j)=>junction.models[j].reference===i).map(v=>amplicons[i].share>0?v/amplicons[i].share:0),(report.amplicons[i].direction==='forward'?a.cut:a.sequence.length-a.cut)-(refs[i].cuts?.[0]??refs[i].cut),a.sequence.length,report.amplicons[i].direction==='reverse').map(c=>({...c,support:emissionGroupSupport(referenceEvidence,models.findIndex(m=>m.reference===i&&pathKey(m)===(c.pathSizes?.join('/')??String(c.size))),predictive?.get(models.findIndex(m=>m.reference===i&&pathKey(m)===(c.pathSizes?.join('/')??String(c.size)))))})):amplicons[i].components,estimable:checksPass&&amplicons[i].share>=.02})),
    mixture:{components,wt:zero,indel:1-zero,inFrame,outOfFrame,r2:fit.r2,rmse:fit.rmse,upstreamError:1-sampleRegistration.r2,bases:positions.length,fitStart:positions[0],fitEnd:positions[positions.length-1]+1,converged:fit.converged,passesChecks:checksPass},
    pooled:{preflight:report,controlR2:controlFit.r2,sampleTrim:sampleRegistration.trim,controlTrim:controlRegistration.trim,amplicons}};
  };
  const rawObserved=positions.flatMap(p=>sampleRaw[p-sampleRegistration.trim]),mean=rawObserved.reduce((s,v)=>s+v,0)/rawObserved.length,rawVariance=rawObserved.reduce((s,v)=>s+(v-mean)**2,0);
  return {fitInput,templateError,finish,report,enhanced,rawVariance,products:input.amplicons};
}
export function analyzePooledSanger(input:PooledAnalysisInput,progress:(stage:string)=>void=()=>{}):SangerAnalysisResult {
 const prepared=preparePooledSanger(input,progress);
 return prepared.finish(solveJunctions(prepared.fitInput,prepared.templateError,prepared.enhanced));
}
