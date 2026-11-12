import {ReferenceSource} from './reference-source';
export interface SangerReference {
  /** Design coordinates use sequence; true records reversal of the supplied input. */
  inputReversed?: boolean;
  name: string;
  sequence: string;
  /** Cut boundaries: N means the cut is after base N in the supplied strand. */
  cuts: number[];
  targetNames?:string[];
  guides?:string[];
  source?:ReferenceSource;
}
export interface SangerDesignSettings {
  sharedPrimerPolicy?: 'exact' | 'progressive';
  /** UI requests one relaxation stage at a time; absent retains automatic search. */
  sharedSearchStage?:0|1|2|3;
  workflow: 'shared-sanger' | 'shared-pcr' | 'homeolog-pcr' | 'standard-pcr';
  platform?: 'sanger' | 'ngs';
  specificityMismatches: number; specificBoth: boolean;
  specificityStrategy:'three-prime'|'mismatch'; threePrimeWindow:1|3;
  primerMin: number; primerMax: number;
  tmMin: number; tmMax: number; gcMin: number; gcMax: number;
  sodiumMm: number; primerNm: number;
  window: number; readStart: number; readEnd: number;
  /** Minimum read-axis gap between every pair of reference/target cuts. */
  minOffset: number; indelRange: number;
  direction: 'both' | 'forward' | 'reverse';
  requirePartner: boolean; productMin: number; productMax: number;
  maxTmDifference: number;
}
export const SANGER_DEFAULTS: SangerDesignSettings = {
  workflow: 'shared-sanger', platform:'sanger', specificityMismatches: 2, specificBoth: false, specificityStrategy:'three-prime', threePrimeWindow:3,
  primerMin: 18, primerMax: 25, tmMin: 52, tmMax: 68, gcMin: 30, gcMax: 70,
  sodiumMm: 50, primerNm: 250, window: 90, readStart: 50, readEnd: 700,
  minOffset: 30, indelRange: 10, direction: 'both', requirePartner: true,
  productMin: 150, productMax: 1500, maxTmDifference: 5,
};
export interface PrimerSite {
  /** One-based positions on the oligo, counted from its 5-prime end. */
  mismatchPositions?:number[];
  reference: string;
  /** Zero-based, half-open coordinates on the supplied strand. */
  start: number; end: number;
  direction: 'forward' | 'reverse';
  cutOffsets: number[];
  extension: string;
}
export interface PrimerProperties {
  sequence: string; tm: number; gc: number;
  selfEnd: number; selfAny: number; hairpinStem: number;
}
export interface CommonPartner extends PrimerProperties {
  /** Present when each reference uses its own opposite PCR oligo. */
  perReference?:PrimerProperties[];
  sites: PrimerSite[]; products: number[]; pairEnd: number;
}
export interface OppositePrimerChoice {primer:PrimerProperties;site:PrimerSite;product:number;safety:PrimerBindingScreen}
export interface ThreePrimeVariant {background:string;position:number;backgroundPosition:number;targetBase:string;backgroundBase:string;fromThreePrime:number;alignmentIdentity:number}
export interface PrimerBindingScreen {
  searchLimit:3;
  /** Exact minimum over resolved ungapped non-intended positions when minimumExact is true. */
  minimumMismatches:number; minimumExact:boolean; nearSiteCount:number;
  sites:{reference:string;start:number;end:number;direction:'forward'|'reverse';mismatches:number;terminalMismatches:number}[];
  references:{reference:string;minimumMismatches:number;minimumExact:boolean;nearSiteCount:number;hasUnknownBases:boolean}[];
}
export interface SangerCandidate extends PrimerProperties {
  oppositeChoices?:OppositePrimerChoice[][];
  commonOppositeChoices?:CommonPartner[];
  individualOppositeDefault?:CommonPartner;
  oppositeSetChoices?:{partner:CommonPartner;screens:PrimerBindingScreen[]}[];
  partnerSafety?:PrimerBindingScreen[];
  sharedDesign?:{stage:'both-exact'|'one-exact'|'single-aligned'|'single-unrestricted'|'double';positionSpread:number;worstFivePrimePosition:number;terminalExact:10};
  primerSafety?:{primary:PrimerBindingScreen;partner:PrimerBindingScreen|null};
  wtSeparation?:{minimumDifferences:number;comparedBases:number;condition:number|null;warnings:string[]};
  id: string; direction: 'forward' | 'reverse'; sites: PrimerSite[];
  targetDistance?: number;
  minOffset: number; minTraceDifferences: number; traceComparedBases: number;
  indelOffsetOverlap: boolean;
  partner: CommonPartner | null;
  referenceIndex?: number;
  specificity?: { primary: boolean; partner: boolean; mismatches: number;strategy?:'three-prime'|'mismatch';terminalWindow?:number;primaryVariants?:ThreePrimeVariant[];partnerVariants?:ThreePrimeVariant[] };
}
export function oppositePrimer(candidate:SangerCandidate,index:number):PrimerProperties|null {
  return candidate.partner?.perReference?.[index] || candidate.partner || null;
}
export interface SangerDesignResult {
  searchedStage?:0|1|2|3;
  candidates: SangerCandidate[];
  stats: { enumerated: number; qualityPassed: number; commonUnique: number;
    geometryPassed: number; offsetPassed: number; partnerChecked: number;separabilityRejected?:number };
  excludedCandidates?:{id:string;reason:string}[];
  complete: boolean;
  coverage?: { reference: string; candidates: number }[];
  references: SangerReference[];
  settings: SangerDesignSettings;
}
