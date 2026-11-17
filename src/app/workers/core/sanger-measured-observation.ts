import {fractionNoiseVariance} from './sanger-fraction-noise';
import {meanQuietBackground} from './sanger-mean-background';
import {prefixPhotonVariance} from './sanger-prefix-residual';
import {referenceTransfer} from './sanger-reference-transfer';
import {cycleFluorescence} from './sanger-cycle-fluorescence';
import {amplitudeTransfer} from './sanger-amplitude-transfer';
import {measuredIntensityResponses} from './sanger-intensity-response';
import type {SangerAnalysisSettings,SangerTrace} from '../../models/sanger-analysis';
import type {Prepared} from './sanger-path-common';
import {registerTrace} from './sanger-pooled';
import {calibrateSignal,correctGain,integratedPeaks,contrastPeaks,mixPhase,phaseNeighbor} from './sanger-signal';
import {calibrateControlMixture} from './sanger-control-calibration';
import {controlUsableRegion} from './sanger-usable-region';

/** Common measured-control observation model. No genotype dictionary, legacy
 * molecule proposal, synthetic selector or group-count limit enters this step. */
export function prepareMeasuredObservation(refs:Prepared[],control:SangerTrace,sample:SangerTrace,settings:SangerAnalysisSettings,progress:(stage:string)=>void,centerWeight=0){
 progress('Registering measured fluorescence');
 const registered=refs.map(r=>({name:r.name,sequence:r.sequence,cut:Math.min(...r.cuts)})),cr=registerTrace(control,registered,settings,integratedPeaks(control,[0,0,0,0]).vectors),sr=registerTrace(sample,registered,settings,integratedPeaks(sample,[0,0,0,0]).vectors),first=Math.min(...refs.flatMap(r=>r.cuts)),last=Math.max(...refs.flatMap(r=>r.cuts)),quality=controlUsableRegion(control,cr.trim,first,last,settings.readStart,settings.readEnd);
 // The automatic display window is not a limit on available inference data.
 // Use every control-qualified flank cycle within the configured read range;
 // stopping at lastCut+200 discarded information needed to stabilize pools.
 const align=(v:number[][],trim:number)=>Array.from({length:v.length+trim},(_,p)=>v[p-trim]||[0,0,0,0]),cContrast=align(contrastPeaks(control).vectors,cr.trim),sContrast=align(contrastPeaks(sample).vectors,sr.trim),end=Math.min(quality.end,sContrast.length-2,...refs.map(r=>r.sequence.length-settings.maxIndel*r.cuts.length)),start=Math.max(quality.start,first-15),positions=Array.from({length:Math.max(0,end-start)},(_,i)=>start+i);
 if(positions.length<25)throw Error('Insufficient reliable downstream measurements for the requested targets.');
 progress('Calibrating the shared observation model');
 const upstreamStart=Math.max(quality.start+10,cr.trim+10,sr.trim+10),calibrationGuard=Math.max(45,(settings.maxDeletion??settings.maxIndel)+3),upstream=Array.from({length:Math.max(0,first-calibrationGuard-upstreamStart)},(_,i)=>upstreamStart+i);if(upstream.length<40)throw Error('Insufficient reliable upstream calibration outside the modeled edit/substitution window.');
 const controlEnd=Math.min(quality.end,...refs.map(r=>r.sequence.length)),controlPositions=Array.from({length:Math.max(0,controlEnd-quality.start)},(_,i)=>quality.start+i),sequences=refs.map(r=>r.sequence),cBackgroundStats=meanQuietBackground(control,sequences,cr.trim,upstream,centerWeight),sBackgroundStats=meanQuietBackground(sample,sequences,sr.trim,upstream,centerWeight),cBackground=cBackgroundStats.background,sBackground=sBackgroundStats.background,controlSignal=cycleFluorescence(control,cBackground,false,centerWeight),sampleSignal=cycleFluorescence(sample,sBackground,true,centerWeight),sampleRawSignal=cycleFluorescence(sample,[0,0,0,0],false,centerWeight),c=align(controlSignal.vectors,cr.trim),s=align(sampleRawSignal.vectors,sr.trim),training=upstream.filter((_,i)=>i%3!==0),validation=upstream.filter((_,i)=>i%3===0),transfer=amplitudeTransfer(controlSignal,sampleSignal,training,cr.trim,sr.trim),gains=transfer.gains,raw=s.map(v=>correctGain(v,gains)),heights=controlSignal.strength.filter(v=>v>0).sort((a,b)=>a-b),scale=heights[Math.floor(heights.length/2)]||1,amplitudes=Array.from({length:c.length},(_,p)=>(controlSignal.strength[p-cr.trim]||0)/scale);
 const cycleWidth=(p:number)=>{const i=p-sr.trim,peak=sample.peaks[i],left=i?sample.peaks[i-1]:peak-12,right=i+1<sample.peaks.length?sample.peaks[i+1]:peak+12,width=Math.min(sample.channels.A.length-1,(peak+right)/2)-Math.max(0,(peak+left)/2);return width;},backgroundAt=(p:number)=>sBackground.map((v,b)=>v*cycleWidth(p)/scale/gains[b]);
 const expected=(p:number,lag:number)=>{const q=phaseNeighbor(p,lag),at=(r:number)=>Math.exp(transfer.logScale+transfer.slope*(2*(r-transfer.first)/transfer.span-1)),background=backgroundAt(p),values=[0,1,2,3].map(b=>(1-Math.abs(lag))*amplitudes[p]*at(p)*(c[p]?.[b]||0)+Math.abs(lag)*(amplitudes[q]||0)*at(q)*(c[q]?.[b]||0)+background[b]),total=values.reduce((n,v)=>n+v,0);return values.map(v=>total?v/total:0);},error=(rows:number[],lag:number,corrected=true)=>rows.reduce((n,p)=>{const pred=expected(p,lag),y=corrected?raw[p]:s[p];return n+y.reduce((n,v,b)=>n+(v-pred[b])**2,0);},0);
 let lag=0,minimum=error(training,0);for(let k=-8;k<=8;k++){const l=k*.025,e=error(training,l);if(e<minimum-1e-12){lag=l;minimum=e;}}
 const variance=validation.reduce((n,p)=>n+raw[p].reduce((n,v)=>n+(v-.25)**2,0),0),rawVariance=validation.reduce((n,p)=>n+s[p].reduce((n,v)=>n+(v-.25)**2,0),0),calibration={gains,lag,used:Math.abs(lag)>0||gains.some(v=>Math.abs(v-1)>1e-8),trainingBases:training.length,validationBases:validation.length,rawValidationR2:rawVariance?1-error(validation,0,false)/rawVariance:1,validationR2:variance?1-error(validation,lag)/variance:1,weights:[1]};
 // Keep peak-contrast noise as a diagnostic. The inference measurements are
 // integrated photons, so their predictive noise must be measured on the same
 // observable, using unedited prefix rows withheld from gain/transfer fitting.
 const noiseCalibration=calibrateSignal(sContrast,upstream,[cContrast]),noiseRaw=sContrast.map(v=>correctGain(v,noiseCalibration.gains)),calibrationNoise=upstream.reduce((sum,p)=>sum+noiseRaw[p].reduce((e,v,b)=>e+(v-mixPhase(cContrast[p],cContrast[phaseNeighbor(p,noiseCalibration.lag)]||[0,0,0,0],noiseCalibration.lag)[b])**2,0),0)/(4*upstream.length),counts=[...inputControlPeaks(control,quality,cr.trim),...inputControlPeaks(sample,quality,sr.trim)].filter(v=>v>0).sort((a,b)=>a-b),quantizationNoise=1/(6*Math.max(1,counts[Math.floor(counts.length/2)]||1)**2),contrastNoise=Math.max(calibrationNoise,quantizationNoise);
 const controlCalibration=calibrateControlMixture(sequences,controlPositions,c),{shares,means}=controlCalibration,intensity=measuredIntensityResponses(refs,shares,means,c,amplitudes,controlPositions,1/scale),transferByReference=referenceTransfer(intensity.responses,shares,sampleSignal.areas,training,sr.trim,scale,transfer,lag);

 const referencePrefixError=prefixPhotonVariance(validation,raw,intensity.responses,shares,transferByReference,lag,backgroundAt),noise=Math.max(quantizationNoise,referencePrefixError),fitDelta=1.5*Math.sqrt(noise);
 const backgroundVarianceAt=(p:number)=>{const width=cycleWidth(p),v=sBackgroundStats.variance.map((x,b)=>x*width*width/scale**2/gains[b]**2),total=sampleRawSignal.areas[p-sr.trim].reduce((n,x,b)=>n+x/scale/gains[b],0);return fractionNoiseVariance(raw[p],v,total);},prefixBackgroundVariance=validation.reduce((n,p)=>n+backgroundVarianceAt(p),0)/validation.length,stableVariance=Math.max(quantizationNoise,noise-prefixBackgroundVariance),noiseAt=(p:number)=>Math.sqrt(stableVariance+backgroundVarianceAt(p));
 return{noiseAt,backgroundVarianceAt,prefixBackgroundVariance,stableVariance,backgroundAt,backgrounds:{control:cBackground,sample:sBackground},cr,sr,quality,positions,c,s,raw,noise,means,calibration,shares,start,end,intensity,fitDelta,transferByReference,controlSignal,sampleSignal,noiseCalibration,contrastNoise,referencePrefixError,controlCalibration,trainingPositions:training,validationPositions:validation};
}

export function inputControlPeaks(control:SangerTrace,quality:{start:number;end:number},trim:number){return Array.from(control.peaks).slice(Math.max(0,quality.start-trim),Math.min(control.peaks.length,quality.end-trim)).map(p=>Object.values(control.channels).reduce((n,channel)=>n+channel[p],0));}
