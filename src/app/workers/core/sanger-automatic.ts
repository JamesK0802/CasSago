import {analyzeUnifiedSequenceGroups} from './sanger-unified-pipeline';
import type{SangerLinearOptimizer}from'./sanger-discrete-phase';
import {analyzePairedSanger} from './sanger-paired';
import {PooledAmplicon,SANGER_ANALYSIS_DEFAULTS,SangerAnalysisResult,SangerAnalysisSettings,SangerTrace} from '../../models/sanger-analysis';
import {analyzeSanger} from './sanger-analysis';
import {analyzePooledSanger} from './sanger-pooled';
export interface AutomaticSangerInput {references:PooledAmplicon[];primer?:string;sample:SangerTrace;control:SangerTrace;maxDeletion:number;maxInsertion:number;additionalReads?:{primer:string;sample:SangerTrace;control:SangerTrace}[];synthetic?:boolean;groupRecovery?:boolean;indelThreshold?:number}
/** One public analysis contract. Instrument calibration, direction and measured
 * coverage are automatic. Reference count selects the observation adapter;
 * requested sequence groups use whole-molecule candidates when supported;
 * the existing length model remains available for other observation adapters. */
export function automaticSangerSettings(input:AutomaticSangerInput):SangerAnalysisSettings {
 if(!Number.isInteger(input.maxDeletion)||input.maxDeletion<1||input.maxDeletion>500||!Number.isInteger(input.maxInsertion)||input.maxInsertion<1||input.maxInsertion>30)throw Error('Use deletion limit 1–500 bp and insertion limit 1–30 bp.');
 if(input.indelThreshold!==undefined&&(!Number.isFinite(input.indelThreshold)||input.indelThreshold<0||input.indelThreshold>100))throw Error('Use an indel threshold from 0 to 100%.');
 return {...SANGER_ANALYSIS_DEFAULTS,indelThreshold:input.indelThreshold??2,analysisPolicy:'automatic',scope:input.references.length===1?'single':'pooled',method:'mixture',signalModel:'calibrated',direction:'auto',readStart:0,readEnd:3160,maxDeletion:input.maxDeletion,maxInsertion:input.maxInsertion,maxIndel:Math.min(30,Math.max(input.maxDeletion,input.maxInsertion)),largeDeletions:false,checkWiderRange:false,jointEngine:'staged',window:200};
}
export function analyzeAutomaticSanger(input:AutomaticSangerInput,progress:(stage:string)=>void=()=>{},optimizer?:SangerLinearOptimizer):SangerAnalysisResult {
 if(!input.sample||!input.control)throw Error('Upload an edited sample and an unedited control chromatogram. Reference DNA cannot replace measured control fluorescence.');
 if(!input.references.length||input.references.length>4)throw Error('Enter 1–4 reference sequences.');
 const targetCounts=input.references.map(r=>r.cuts?.length||1);if(targetCounts.some(n=>n<1||n>2)||new Set(targetCounts).size!==1)throw Error('Use one or two corresponding targets in every reference.');
 const multi=targetCounts[0]>1;if(multi&&!input.primer)throw Error('A sequencing primer is required to register multiple targets.');
 const settings=automaticSangerSettings(input),only=input.references[0];
 if(input.additionalReads?.length){
  if(input.additionalReads.length!==1)throw Error('Add one opposite-direction read pair.');
  const pairs=[{primer:input.primer||'',sample:input.sample,control:input.control},...input.additionalReads];
  if(pairs.some(p=>!p.sample||!p.control||!/^[ACGT]{12,40}$/.test(p.primer.replace(/\s/g,'').toUpperCase())))throw Error('Each read direction needs its sequencing primer, edited trace and WT control.');
  return analyzePairedSanger(pairs.map(p=>({...p,primer:p.primer.replace(/\s/g,'').toUpperCase(),amplicons:input.references,ampliconConfirmed:true,settings,synthetic:input.synthetic})),progress);
 }
 if(input.groupRecovery&&input.primer)return analyzeUnifiedSequenceGroups({amplicons:input.references,primer:input.primer.replace(/\s/g,'').toUpperCase(),ampliconConfirmed:true,sample:input.sample,control:input.control,settings:{...settings,readEnd:1000,editStartRadius:3,maxIndel:input.maxDeletion},synthetic:input.synthetic},progress,optimizer);
 const result=input.references.length===1&&!multi?analyzeSanger({referenceName:only.name,reference:only.sequence,cut:only.cut,sample:input.sample,control:input.control,settings,synthetic:input.synthetic},progress):analyzePooledSanger({amplicons:input.references,primer:input.primer?.replace(/\s/g,'').toUpperCase()||'',ampliconConfirmed:true,sample:input.sample,control:input.control,settings,synthetic:input.synthetic},progress);
 return result;
}
