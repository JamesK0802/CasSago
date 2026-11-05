import type {ExactFraction} from '../workers/core/exact-mixture-fiber';
export const TRACE_BASES = ['A','C','G','T'] as const;
export type TraceBase = typeof TRACE_BASES[number];
export interface SangerTrace {
  synthetic?:boolean; name:string; sampleName:string; sequence:string;
  peaks:Uint32Array; quality:Uint8Array|null;
  channels:Record<TraceBase,Float32Array>;
}
export interface SangerAnalysisSettings {
  analysisPolicy?:'automatic';
  editStartRadius?:number;largeDeletions?:boolean;maxDeletion?:number;maxInsertion?:number;jointEngine?:'staged'|'signal-local'|'signal-path'|'signal-hybrid';
  scope:'single'|'pooled'; readStart:number; readEnd:number;
  method:'mixture'|'consensus'; maxIndel:number; window:number; controlPhred:number;
  direction:'auto'|'forward'|'reverse';joinNearby?:number;signalModel?:'classic'|'calibrated';checkWiderRange?:boolean;
}
export const SANGER_ANALYSIS_DEFAULTS:SangerAnalysisSettings={scope:'single',editStartRadius:5,readStart:35,readEnd:700,method:'mixture',maxIndel:10,window:200,controlPhred:15,direction:'auto',joinNearby:3,signalModel:'calibrated',checkWiderRange:true};
export interface TraceAlignmentColumn {ref:number|null; query:number|null; a:string;b:string}
export interface TraceAlignment {score:number;identity:number;columns:TraceAlignmentColumn[];direction:'forward'|'reverse';cut:number;referenceStart:number;referenceEnd:number}
export interface SangerVariant {position:number;reference:string;sample:string;netSize:number;type:'Substitution'|'In-frame indel'|'Out-of-frame indel'|'Complex replacement'|'Ambiguous call';endPosition?:number;grouped?:boolean;parts?:SangerVariant[];matchedBridgeBases?:number}
export interface TracePlot {cutIndex?:number;name:string;start:number;end:number;paths:Record<TraceBase,string>;width:number;cutX:number;cuts?:{name:string;x:number}[];letters:{x:number;base:string;offset:number}[]}
export interface SangerExactOption {id:string;components:SangerAnnotation['components'];exactFractions:ExactFraction[]}
export interface SangerExactAlternatives {reference:string;certificate:'exact-rational-emission';selected:number;options:SangerExactOption[];targets:number[];ambiguousSequences:Record<string,string[]>;ambiguousPlacements?:Record<string,string[]>;exhaustive:false}
export interface SangerAnalysisResult {
  exactAlternatives?:SangerExactAlternatives[];
  method:'mixture'|'consensus'|'pooled';referenceName:string;sampleName:string;controlName?:string;direction:'forward'|'reverse';
  cutPosition:number;alignmentIdentity:number;alignment:TraceAlignmentColumn[];variants:SangerVariant[];
  plots:TracePlot[];warnings:string[];synthetic:boolean;
  reads?:{primer:string;sampleName:string;controlName:string;preflight:PooledPreflight;sampleTrim:number;controlTrim:number;r2:number;rmse:number;bases:number;fitStart:number;fitEnd:number}[];
  joint?:{attributionAmbiguous?:{reference:string;target:number;possibleReallocatedFraction:number}[];groups:{id:string;reference:string|null;fraction:number;netSize:number|null;pathSizes?:number[];events:SangerJointEvent[];alternatives:{reference:string;events:SangerJointEvent[]}[]}[];targets:{reference:string;target:number;cut:number;edited:number|null}[];reconstructionR2?:number;heldOutR2:number|null;selectedModels:number;candidateModels:number;iterations:number;maxDeletion?:number;ambiguousFraction:number;alignmentSearches?:number};
  annotations?:SangerAnnotation[];
  diagnostics?:{substitutionLinkage?:{reference:string;target:number;blocks:number;fraction:number;wtRange:[number,number];modelSequence:string}[];readAgreementPp?:number;controlTransfer?:{used:boolean;gains:number[];shares:number[];rawMse:number;colourMse:number;colourR2:number;initialRawMse:number;initialColourMse:number;trainingBases:number;validationBases:number};controlResponse?:{used:boolean;capacity:number;table:number[][];endTable?:number[][];knots?:number[][][];start:number;end:number;shares:number[];rawValidationMse:number;validationMse:number;baseValidationMse:number;trainingBases:number;validationBases:number;refitBases:number};transitionResponse?:{used:boolean;table:number[][];trainingCounts:number[];validationCounts:number[];rawValidationMse:number;validationMse:number};junctionFit?:{positions:number[];iterations:number;models:{reference:number;size:number;start:number;inserted:string;fraction:number;probabilities?:number[][];origins?:{start:number;weight:number;probabilities?:number[][]}[]}[];duplicateLengths:number};referenceStability?:{jointTVpp:number;maximumConditionalTVpp:number;maximumShareDifferencePp:number};modeledRange?:{maxDeletion:number;maxInsertion:number};latentInsertion?:{unknownCycles:number;knownTrainingRows:number;knownFitRows:number;knownFlankR2:number;probabilities:{position:number;distribution:number[]}[]};algorithmVersion?:string;signalModel:string;controlCalibration?:SangerCalibration;sampleCalibration?:SangerCalibration;splitStabilityPp?:number;discordantPositions?:number;maxPositionRmse?:number;rangeSuggestion?:{maxIndel:number;validationR2Gain:number;outsideFraction:number};requestedMaxIndel:number;supportScreen?:{predictionRemoved?:{reference:number;size:number;foldGains:number[]}[];removed:{reference:number;size:number}[];threshold:number;evidence?:{group:number[];threshold:number;mass:number;standardError:number;score:number;reference:number;size:number}[];unpruned?:{reference:number;size:number;jointFraction:number}[]};contextEmission?:{source?:'control-strength'|'normalized-control';theta:number[];used:boolean;trainingPositions:number;validationPositions:number;rawValidationMse:number;validationMse:number};background?:{method:string;control:{levels:number[];trainingSamples:number[]};sample:{levels:number[];trainingSamples:number[]}}};
  pooled?:{preflight:PooledPreflight;controlR2:number;sampleTrim:number;controlTrim:number;amplicons:{name:string;share:number;wt:number|null;indel:number|null;inFrame:number|null;outOfFrame:number|null;components:{size:number;fraction:number}[]}[]};
  mixture?:{components:{size:number;fraction:number}[];wt:number;inFrame:number;outOfFrame:number;indel:number;r2:number;rmse:number;upstreamError:number;bases:number;fitStart:number;fitEnd:number;converged:boolean;passesChecks:boolean};
}

export interface PooledAmplicon {forwardPrimer?:string;reversePrimer?:string;name:string;sequence:string;cut:number;cuts?:number[]}
export interface PooledPreflight {passed:boolean;reason:string;positions:number[];modelCount:number;rank:number;condition:number|null;minimumDifferences:number;minimumCutOffset:number;
  amplicons:{name:string;direction:'forward'|'reverse';primerStart:number;primerEnd:number;cutOffset:number;windowStart:number;windowEnd:number;targets?:{cutOffset:number;windowStart:number;windowEnd:number}[];maximumWindow?:number;maximumIndependentWindow?:number;mergedWindows?:boolean}[];joint?:{mode:string;candidateCount:number;ambiguousClasses:number;coverageExcluded:number;maxDeletion:number};}

export interface SangerCalibration {gains:number[];lag:number;used:boolean;trainingBases:number;validationBases:number;rawValidationR2:number;validationR2:number}
export interface SangerJointEvent {target:number;start:number;end:number;inserted:number;insertedSequence?:string;insertedBaseProbabilities?:number[][]}
export interface SangerGroupSupport {prediction?:{foldGains:number[];passed:boolean;lengthFoldGains?:number[];lengthPassed?:boolean};status:'supported'|'unresolved';score:number|null;requiredScore:number|null;scope:'reference-length'|'pooled-length'}
export interface SangerAnnotation {sequenceGrouping?:boolean;name:string;reference:string;cut:number;cuts?:number[];direction:'forward'|'reverse';share:number;components:{size:number;pathSizes?:number[];fraction:number;events?:SangerJointEvent[];unresolved?:boolean;support?:SangerGroupSupport}[];estimable:boolean}
