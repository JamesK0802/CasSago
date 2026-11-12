import {referenceTransfer} from './sanger-reference-transfer';
import {cycleFluorescence} from './sanger-cycle-fluorescence';
import {amplitudeTransfer} from './sanger-amplitude-transfer';
import {measuredIntensityResponses} from './sanger-intensity-response';
import {alleleComplexity} from './sanger-allele-complexity';
import {marginalProposals} from './sanger-marginal-proposals';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import {priceMolecules} from './sanger-molecule-pricing';
import type {SangerAnalysisSettings,SangerTrace} from '../../models/sanger-analysis';
import type {Prepared,PathEdit} from './sanger-path-common';
import {materialize} from './sanger-path-common';
import {registerTrace} from './sanger-pooled';
import {calibrateSignal,correctGain,integratedPeaks,contrastPeaks,mixPhase,phaseNeighbor,rawSignalScore,simplexFit,referenceBackgroundPeaks} from './sanger-signal';
import {calibrateControlMixture} from './sanger-control-calibration';
import {controlUsableRegion} from './sanger-usable-region';
import {measuredLasso,measuredGramCache} from './sanger-measured-lasso';
const bases='ACGT';
export function prepareMeasuredObservation(refs:Prepared[],control:SangerTrace,sample:SangerTrace,settings:SangerAnalysisSettings,progress:(stage:string)=>void){
 progress('Registering measured fluorescence');
 const registered=refs.map(r=>({name:r.name,sequence:r.sequence,cut:Math.min(...r.cuts)})),cr=registerTrace(control,registered,settings,integratedPeaks(control,[0,0,0,0]).vectors),sr=registerTrace(sample,registered,settings,integratedPeaks(sample,[0,0,0,0]).vectors),first=Math.min(...refs.flatMap(r=>r.cuts)),last=Math.max(...refs.flatMap(r=>r.cuts)),quality=controlUsableRegion(control,cr.trim,first,last,settings.readStart,settings.readEnd);
 // The automatic display window is not a limit on available inference data.
 // Use every control-qualified flank cycle within the configured read range;
 // stopping at lastCut+200 discarded information needed to stabilize pools.
 const align=(v:number[][],trim:number)=>Array.from({length:v.length+trim},(_,p)=>v[p-trim]||[0,0,0,0]),cContrast=align(contrastPeaks(control).vectors,cr.trim),sContrast=align(contrastPeaks(sample).vectors,sr.trim),end=Math.min(quality.end,sContrast.length-2,...refs.map(r=>r.sequence.length-settings.maxIndel*r.cuts.length)),start=Math.max(quality.start,first-15),positions=Array.from({length:Math.max(0,end-start)},(_,i)=>start+i);
 if(positions.length<25)throw Error('Insufficient reliable downstream measurements for the requested targets.');
 progress('Calibrating the shared observation model');
 const upstreamStart=Math.max(quality.start+10,cr.trim+10,sr.trim+10),calibrationGuard=Math.max(45,(settings.maxDeletion??settings.maxIndel)+3),upstream=Array.from({length:Math.max(0,first-calibrationGuard-upstreamStart)},(_,i)=>upstreamStart+i);if(upstream.length<40)throw Error('Insufficient reliable upstream calibration outside the modeled edit/substitution window.');
 const controlEnd=Math.min(quality.end,...refs.map(r=>r.sequence.length)),controlPositions=Array.from({length:Math.max(0,controlEnd-quality.start)},(_,i)=>quality.start+i),sequences=refs.map(r=>r.sequence),cBackground=referenceBackgroundPeaks(control,sequences,cr.trim,upstream).background!,sBackground=referenceBackgroundPeaks(sample,sequences,sr.trim,upstream).background!,controlSignal=cycleFluorescence(control,cBackground),sampleSignal=cycleFluorescence(sample,sBackground),c=align(controlSignal.vectors,cr.trim),s=align(sampleSignal.vectors,sr.trim),training=upstream.filter((_,i)=>i%3!==0),validation=upstream.filter((_,i)=>i%3===0),transfer=amplitudeTransfer(controlSignal,sampleSignal,training,cr.trim,sr.trim),gains=transfer.gains,raw=s.map(v=>correctGain(v,gains)),heights=controlSignal.strength.filter(v=>v>0).sort((a,b)=>a-b),scale=heights[Math.floor(heights.length/2)]||1,amplitudes=Array.from({length:c.length},(_,p)=>(controlSignal.strength[p-cr.trim]||0)/scale);
 const expected=(p:number,lag:number)=>{const q=phaseNeighbor(p,lag),values=[0,1,2,3].map(b=>(1-Math.abs(lag))*amplitudes[p]*(c[p]?.[b]||0)+Math.abs(lag)*(amplitudes[q]||0)*(c[q]?.[b]||0)),total=values.reduce((n,v)=>n+v,0);return values.map(v=>total?v/total:0);},error=(rows:number[],lag:number,corrected=true)=>rows.reduce((n,p)=>{const pred=expected(p,lag),y=corrected?raw[p]:s[p];return n+y.reduce((n,v,b)=>n+(v-pred[b])**2,0);},0);
 let lag=0,minimum=error(training,0);for(let k=-8;k<=8;k++){const l=k*.025,e=error(training,l);if(e<minimum-1e-12){lag=l;minimum=e;}}
 const variance=validation.reduce((n,p)=>n+raw[p].reduce((n,v)=>n+(v-.25)**2,0),0),rawVariance=validation.reduce((n,p)=>n+s[p].reduce((n,v)=>n+(v-.25)**2,0),0),calibration={gains,lag,used:Math.abs(lag)>0||gains.some(v=>Math.abs(v-1)>1e-8),trainingBases:training.length,validationBases:validation.length,rawValidationR2:rawVariance?1-error(validation,0,false)/rawVariance:1,validationR2:variance?1-error(validation,lag)/variance:1,weights:[1]};
 // Keep peak-contrast noise as a diagnostic. The inference measurements are
 // integrated photons, so their predictive noise must be measured on the same
 // observable, using unedited prefix rows withheld from gain/transfer fitting.
 const noiseCalibration=calibrateSignal(sContrast,upstream,[cContrast]),noiseRaw=sContrast.map(v=>correctGain(v,noiseCalibration.gains)),calibrationNoise=upstream.reduce((sum,p)=>sum+noiseRaw[p].reduce((e,v,b)=>e+(v-mixPhase(cContrast[p],cContrast[phaseNeighbor(p,noiseCalibration.lag)]||[0,0,0,0],noiseCalibration.lag)[b])**2,0),0)/(4*upstream.length),counts=[...inputControlPeaks(control,quality,cr.trim),...inputControlPeaks(sample,quality,sr.trim)].filter(v=>v>0).sort((a,b)=>a-b),quantizationNoise=1/(6*Math.max(1,counts[Math.floor(counts.length/2)]||1)**2),contrastNoise=Math.max(calibrationNoise,quantizationNoise);
 const controlCalibration=calibrateControlMixture(sequences,controlPositions,c),{shares,means}=controlCalibration,intensity=measuredIntensityResponses(refs,shares,means,c,amplitudes,controlPositions,1/scale),transferByReference=referenceTransfer(intensity.responses,shares,sampleSignal.areas,training,sr.trim,scale,transfer);
 const referencePrefixError=validation.reduce((sum,p)=>{const predicted=[0,1,2,3].map(b=>refs.reduce((n,_,h)=>n+shares[h]*transferByReference.apply(h,p,intensity.responses[h][p])[b],0)),total=predicted.reduce((a,b)=>a+b,0);return sum+predicted.reduce((e,v,b)=>e+(v/Math.max(1e-30,total)-raw[p][b])**2,0);},0)/(4*validation.length),noise=Math.max(quantizationNoise,referencePrefixError),fitDelta=1.5*Math.sqrt(noise);
 return{cr,sr,quality,positions,c,s,raw,noise,means,calibration,shares,start,end,intensity,fitDelta,transferByReference,controlSignal,sampleSignal,noiseCalibration,contrastNoise,referencePrefixError,controlCalibration,trainingPositions:training,validationPositions:validation};
}
export function fitMeasuredMixture(refs:Prepared[],control:SangerTrace,sample:SangerTrace,settings:SangerAnalysisSettings,progress:(stage:string)=>void,optimizer:SangerLinearOptimizer){
 const {cr,sr,quality,positions,c,s,raw,noise,means,calibration,shares,start,end}=prepareMeasuredObservation(refs,control,sample,settings,progress),observed=positions.flatMap(p=>raw[p]);
 const models:{reference:number;events:PathEdit[];sequence:string;column:number[]}[]=[],seen=new Set<string>();
 const add=(h:number,events:PathEdit[])=>{const ref=refs[h],sequence=materialize(ref.sequence,events);if(!sequence)return;const key=h+':'+sequence;if(seen.has(key))return;seen.add(key);if((models.length+1)*observed.length>20_000_000)throw Error('Measured candidate basis exceeds the browser work budget; analysis was not certified.');const donors:number[]=[];let at=0;for(const e of [...events].sort((a,b)=>a.start-b.start)){for(let p=at;p<e.start;p++)donors.push(p);for(let p=0;p<e.inserted.length;p++)donors.push(-1);at=e.end;}for(let p=at;p<ref.sequence.length;p++)donors.push(p);const one=(p:number)=>{const donor=donors[p],base=sequence[p];if(refs.length===1&&donor>=0&&c[donor])return c[donor];if(base==='N')return[0,1,2,3].map(b=>means.reduce((n,row)=>n+row[b]/4,0));const b=bases.indexOf(base);return b<0?[0,0,0,0]:means[b];};models.push({reference:h,events,sequence,column:positions.flatMap(p=>mixPhase(one(p),one(phaseNeighbor(p,calibration.lag)),calibration.lag))});};
 const maxDeletion=settings.maxDeletion??settings.maxIndel,maxInsertion=settings.maxInsertion??settings.maxIndel,range=Array.from({length:maxDeletion+maxInsertion+1},(_,i)=>i-maxDeletion),event=(cut:number,size:number,start=cut-Math.floor(Math.max(0,-size)/2)):PathEdit[]=>size?[{start,end:start+Math.max(0,-size),inserted:size>0?'N'.repeat(size):''}]:[];
 progress('Constructing full-molecule hypotheses');
 for(let h=0;h<refs.length;h++)add(h,[]);
 // The same measured-data initializer is attempted for every input. Its
 // compact hypotheses select a representative; they are not identity proofs.
 const partial:{reference:number;events:PathEdit[];fraction:number}[]=[],heights=Array.from(inputControlPeaks(control,quality,cr.trim)).sort((a,b)=>a-b),proposalNoise=Math.max(Math.sqrt(noise),1/Math.sqrt(6)/Math.max(1,heights[Math.floor(heights.length/2)]||1));let initialization:{status:'complete'|'partial'|'unavailable';candidates:number;reason?:string};
 const nullFit=simplexFit(models.map(m=>m.column),observed);
 if(nullFit.rmse<1e-8){initialization={status:'complete',candidates:models.length};}
 else try{const coded=raw.map(v=>simplexFit(means,v).weights),seed=marginalProposals(refs,coded,start,end,maxDeletion,maxInsertion,proposalNoise,optimizer,progress,undefined,rows=>partial.push(...rows));for(const candidate of seed)add(candidate.reference,candidate.events);initialization={status:'complete',candidates:seed.length};}
 catch(error){for(const candidate of partial)add(candidate.reference,candidate.events);initialization={status:partial.length?'partial':'unavailable',candidates:partial.length,reason:String(error)};}

 models.push({reference:-1,events:[],sequence:'',column:positions.flatMap(()=>[.25,.25,.25,.25])});
 progress('Fitting measured mixtures and blocked validation');
 const shared=positions.map(p=>refs.some(ref=>ref.cuts.some(cut=>p>=cut-maxDeletion-3-maxDeletion*(cut!==ref.cuts[0]?1:0)&&p<=cut+maxInsertion+3+maxInsertion*(cut!==ref.cuts[0]?1:0)))),folds=positions.flatMap((p,i)=>Array(4).fill(shared[i]?-1:Math.floor(p/5)%3));
 if(folds.filter(f=>f===2).length<24)throw Error('Insufficient independent flanking signal for measured prediction validation.');
 const gramState=measuredGramCache();
 const regress=(selection=false)=>measuredLasso(models.map(m=>m.column),observed,models.map(m=>m.reference<0?0:alleleComplexity(m.events,refs[m.reference].cuts)),noise,folds,positions.flatMap(p=>Array(4).fill(Math.floor(p/5))),selection,gramState);
 let regression=regress();
 // Common residual column generation for measured and simulated traces.
 let pricingConverged=false,pricingStop:'noise-consistent'|'optimality'='optimality';
 for(let round=0;round<200;round++){
  const w=regression.developmentFit.weights,predicted=observed.map((_,k)=>models.reduce((n,m,i)=>n+w[i]*m.column[k],0)),residual=observed.map((v,k)=>v-predicted[k]);
  const training=residual.filter((_,k)=>folds[k]!==2),blocks=new Set(positions.filter((_,i)=>folds[4*i]!==2).map(p=>Math.floor(p/5))).size,noiseLimit=noise*(1+2*Math.sqrt(2/Math.max(1,blocks)));
  // Stop expanding a statistically adequate dictionary. This development-only
  // discrepancy rule applies to every observation, including measured reads;
  // the independent outer fold still decides whether the final result passes.
  if(regression.developmentFit.converged&&training.reduce((n,e)=>n+e*e,0)/training.length<=noiseLimit){pricingConverged=true;pricingStop='noise-consistent';break;}
  const penalty=Math.max(1e-10,noise*.1),count=folds.filter(f=>f!==2).length,cost=(m:typeof models[number])=>m.reference<0?0:alleleComplexity(m.events,refs[m.reference].cuts),lambda=models.reduce((n,m,i)=>n+w[i]*(m.column.reduce((v,x,k)=>v-(folds[k]!==2?x*residual[k]/count:0),0)+penalty*cost(m)),0),proposals=priceMolecules(refs,maxDeletion,maxInsertion,positions,residual,folds,means,c,calibration.lag,penalty),before=models.length;
  for(const proposal of proposals)if(proposal.cost<lambda-Math.max(1e-7,noise*.001))add(proposal.reference,proposal.events);
  if(models.length===before){pricingConverged=true;break;}
  progress('Refining whole-molecule dictionary · '+(round+1)+' · '+models.length+' candidates');regression=regress();
 }
 if(!pricingConverged)throw Error('Common molecule pricing exhausted its work budget; no result was certified.');
 regression=regress(true);
 const columns=models.map(m=>m.column),weights=regression.fit.weights,score=rawSignalScore(columns,positions.flatMap(p=>s[p]),weights,calibration.gains),nullScore=rawSignalScore(columns,positions.flatMap(p=>s[p]),regression.nullFit.weights,calibration.gains);
 const predict=(w:number[],k:number)=>columns.reduce((sum,col,i)=>sum+col[k]*w[i],0),blockGains=new Map<number,number[]>();let cvError=0,variance=0,predictorDifference=0;const mean=observed.reduce((a,b)=>a+b,0)/observed.length;
 let validationCount=0;const validationMean=observed.filter((_,k)=>folds[k]===2).reduce((a,b)=>a+b,0)/folds.filter(f=>f===2).length;
 for(let k=0;k<observed.length;k++){if(folds[k]!==2)continue;validationCount++;const prediction=predict(regression.developmentFit.weights,k),nullPrediction=predict(regression.nullFit.weights,k),e=(observed[k]-prediction)**2;predictorDifference+=(prediction-nullPrediction)**2;cvError+=e;variance+=(observed[k]-validationMean)**2;const id=Math.floor(positions[Math.floor(k/4)]/5),g=blockGains.get(id)||[];g.push((observed[k]-nullPrediction)**2-e);blockGains.set(id,g);}
 const gains=[...blockGains.values()].map(g=>g.reduce((a,b)=>a+b,0)/g.length),gain=gains.reduce((a,b)=>a+b,0)/gains.length,blockGainSE=Math.sqrt(gains.reduce((s,v)=>s+(v-gain)**2,0)/Math.max(1,gains.length-1)/gains.length),gainSe=2*Math.sqrt(3*Math.max(noise,cvError/validationCount)*predictorDifference)/validationCount,cvR2=variance?1-cvError/variance:1,editing=(w:number[])=>w.reduce((n,v,i)=>n+(models[i].events.length?v:0),0),unresolved=weights[models.findIndex(m=>m.reference<0)],edited=editing(weights),foldEdited=regression.foldFits.map(f=>editing(f.weights)),hasEditingEvidence=gain>2*gainSe,wtConsistent=nullScore.rmse<=Math.max(.005,3*Math.sqrt(noise)),cvRmse=Math.sqrt(cvError/validationCount),predictionLimit=Math.max(.005,3*Math.sqrt(noise)),profileSearchComplete=initialization.status==='complete',signalPassesChecks=score.rmse<=predictionLimit&&cvRmse<=predictionLimit&&regression.developmentFit.converged&&regression.fit.converged&&regression.foldFits.every(f=>f.converged)&&calibration.validationR2>=.8&&(hasEditingEvidence||wtConsistent),passesChecks=profileSearchComplete&&signalPassesChecks;
 return{profileSearchComplete,signalPassesChecks,pricingStop,controlVectors:c,means,initialization,regularization:regression.alphaScale,blockGainSE,cvRmse,predictionLimit,foldDiagnostics:regression.foldFits.map(f=>({converged:f.converged,gap:f.optimalityGap,iterations:f.iterations})),refs,models,weights,score,nullScore,noise,calibration,quality,positions,controlTrim:cr.trim,sampleTrim:sr.trim,shares,unresolved,edited,foldEdited,cvR2,gain,gainSe,hasEditingEvidence,passesChecks,converged:regression.fit.converged,iterations:regression.fit.iterations};
}

function inputControlPeaks(control:SangerTrace,quality:{start:number;end:number},trim:number){return Array.from(control.peaks).slice(Math.max(0,quality.start-trim),Math.min(control.peaks.length,quality.end-trim)).map(p=>Object.values(control.channels).reduce((n,channel)=>n+channel[p],0));}
