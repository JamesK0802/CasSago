import {SangerAnalysisResult} from '../models/sanger-analysis';
import {exactDisplayView} from './sanger-exact-display';
import {sangerGroupKey} from './sanger-curation';
export type SangerCurationScope='file'|'gene'|'target';
export const sangerScopeKey=(kind:SangerCurationScope,reference='',target=0)=>JSON.stringify([kind,reference,kind==='target'?target:null]);
export function sangerScopeExcluded(scopes:readonly string[],reference:string,target?:number){return scopes.includes(sangerScopeKey('file'))||scopes.includes(sangerScopeKey('gene',reference))||target!==undefined&&scopes.includes(sangerScopeKey('target',reference,target));}
/** Expand parent exclusions for presentation only. Child selections are retained
 * independently, so restoring a parent cannot undo earlier group curation. */
export function effectiveSangerExclusions(result:SangerAnalysisResult,groups:readonly string[],scopes:readonly string[]){const keys=new Set(groups);for(const a of result.annotations||[])for(let target=0;target<(a.cuts?.length||1);target++)if(sangerScopeExcluded(scopes,a.name,target))for(const row of exactDisplayView(a,target,45,result.exactAlternatives?.find(s=>s.reference===a.name)).rows)keys.add(sangerGroupKey(a.name,target,row.id));return [...keys];}
