import {ReferenceSource} from './reference-source';
export interface GuideReference {name:string;sequence:string;selected:boolean;regionStart:number|null;regionEnd:number|null;source?:ReferenceSource}
export interface GuideSettings {mode:'single'|'shared'|'specific';gcMin:number;gcMax:number;mismatches:number;uniqueTargets:boolean;excludePolyT:boolean;maxCandidates:number}
export const GUIDE_DEFAULTS:GuideSettings={mode:'single',gcMin:35,gcMax:65,mismatches:2,uniqueTargets:true,excludePolyT:false,maxCandidates:1000};
export interface GuideSite {reference:number;start:number;end:number;pamStart:number;pamEnd:number;pam:string;strand:'+'|'-';cut:number;sequence:string}
export interface GuideCandidate {id:string;sequence:string;gc:number;polyT:boolean;sites:GuideSite[];otherHits:number[];nearestOtherMismatch:number|null;exampleOtherSites:(GuideSite&{mismatches:number})[]}
export interface GuideResult {references:GuideReference[];backgrounds:{name:string;sequence:string}[];settings:GuideSettings;candidates:GuideCandidate[];eligible:number;screened:number;complete:boolean;sitesScanned:number}
