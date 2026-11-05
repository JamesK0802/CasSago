import {guardCompetingTargetPaths} from './sanger-path-agreement';
import {guardEquivalentJunctions} from './sanger-path-identifiability';
import {analyzeSignalPaths as known} from './sanger-path-known';
import {analyzeSignalPaths as latent} from './sanger-path-latent';
import {SignalPathInput} from './sanger-path-common';
import {SangerAnalysisResult} from '../../models/sanger-analysis';
export {discoverSignalPath,assertSignalTraceGrid,checkSignalPathSeparation,eventsRespectIndelOrigin} from './sanger-path-common';
export type {PathEdit,DiscoveredPath,SignalPathInput} from './sanger-path-common';
/** First certify known-sequence paths. If that fails, try a bounded insertion
 * nuisance model fitted on known flanks. Neither attempt can relax its QC.
 * These reconstruction checks are not independent validation. */
export function analyzeSignalPaths(input:SignalPathInput,progress:(s:string)=>void=()=>{},audit?:(state:any)=>void):SangerAnalysisResult{
 const primary=known(input,progress,audit);
 if(primary.mixture?.passesChecks)return guardEquivalentJunctions(primary,input.amplicons,input.settings.editStartRadius);
 progress('Checking unresolved insertion scaffolds against known flanks');
 try{
  const alternative=latent(input,progress,audit);
  if(alternative.mixture?.passesChecks){alternative.warnings.push('Known-flank fallback used after known-sequence inference was not certified.');return guardEquivalentJunctions(guardCompetingTargetPaths(primary,alternative,input.amplicons),input.amplicons,input.settings.editStartRadius);}
  // Preserve the first explanation when neither strategy is certified.
  primary.warnings.push('Known-flank fallback also failed certification; estimates remain unresolved.');
 }catch(error){primary.warnings.push('Known-flank fallback did not complete: '+(error instanceof Error?error.message:String(error)));}
 return guardEquivalentJunctions(primary,input.amplicons,input.settings.editStartRadius);
}
