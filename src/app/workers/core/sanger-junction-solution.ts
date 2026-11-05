import {findPathAliases} from './sanger-path-aliases';
import {fitJunctionMixture,JunctionFit,JunctionInput} from './sanger-junction-mixture';
import {refineJunctionSupport,rejectNonpredictiveJunctions} from './sanger-junction-support';
import {junctionPredictionEvidence} from './sanger-junction-evidence';
import {isWildType} from './sanger-junction-paths';
/** Selection/evidence is identical for a single observation block and joined
 * reads. No input read gets to fit a different molecular distribution. */
export function solveJunctions(input:JunctionInput,templateError:(fit:JunctionFit)=>number[]|undefined,evidence=true){
 const junction=fitJunctionMixture(input),unprunedWeights=[...junction.weights];
 const support=evidence?refineJunctionSupport(input,junction,templateError(junction)):undefined;
 const referenceEvidence=support?.referenceScores;
 const predictive=evidence?junctionPredictionEvidence(input,junction,junction.models.flatMap((m,i)=>!isWildType(m)&&junction.weights[i]>1e-8?[i]:[])):undefined;
 const predictionRemoved=predictive?rejectNonpredictiveJunctions(input,junction,predictive,referenceEvidence,support!.scores):[];
 if(predictionRemoved.length)support!.removed.push(...predictionRemoved.filter(i=>!support!.removed.includes(i)));
 const aliases=findPathAliases(input,junction);
 return {aliases,junction,unprunedWeights,support,referenceEvidence,predictive,predictionRemoved,stabilityInput:input,stabilityColumns:junction.columns,coveragePass:undefined as boolean|undefined};
}
export type JunctionSolution=ReturnType<typeof solveJunctions>;
