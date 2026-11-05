import {pathKey,isWildType} from './sanger-junction-paths';
import {fitJunctionMixture,JunctionFit,JunctionInput,moleculeJunctionStarts} from './sanger-junction-mixture';
import {screenEmissionSupport,EmissionSupportEvidence} from './sanger-emission';
function refitSelected(input:JunctionInput,junction:JunctionFit,keep:number[]) {
 const refit=fitJunctionMixture({...input,classes:keep.map(i=>({reference:junction.models[i].reference,size:junction.models[i].size,sizes:junction.models[i].sizes}))}),weights=junction.models.map(()=>0);
 keep.forEach((id,k)=>{weights[id]=refit.weights[k];junction.models[id]=refit.models[k];junction.columns[id]=refit.columns[k];});
 const originalPredict=junction.predictColumns;junction.predictColumns=rows=>{const out=originalPredict(rows),updated=refit.predictColumns(rows);keep.forEach((id,k)=>out[id]=updated[k]);return out;};
 Object.assign(junction,{weights,r2:refit.r2,rmse:refit.rmse,converged:refit.converged,iterations:junction.iterations+refit.iterations,objective:refit.objective,objectiveHistory:refit.objectiveHistory});
 return {weights,r2:refit.r2,rmse:refit.rmse,converged:refit.converged};
}
/** WT-uncertainty selection followed by the SAME junction/base refit. */
export function refineJunctionSupport(input:JunctionInput,junction:JunctionFit,templateResidual?:number[]) {
 const mandatory=junction.models.flatMap((m,i)=>isWildType(m)?[i]:[]),groups=[...new Set(junction.models.map(pathKey))].map(size=>junction.models.flatMap((m,i)=>pathKey(m)===size?[i]:[]));
 const hypotheses=junction.models.map(m=>m.path?m.path.reduce((n,e)=>n*Math.max(1,moleculeJunctionStarts(input.references[m.reference].sequence,input.references[m.reference].cuts![e.target-1],e.size,input.radius).length),1):m.size?moleculeJunctionStarts(input.references[m.reference].sequence,input.references[m.reference].cut,m.size,input.radius).length:1);
 const support=screenEmissionSupport(junction.columns,input.observed,junction,mandatory,input.positionWeights,groups,hypotheses,{templateResidual});
 if(support.removed.length)support.fit=refitSelected(input,junction,junction.models.flatMap((m,i)=>mandatory.includes(i)||support.fit.weights[i]>1e-10?[i]:[]));
 else junction.weights=support.fit.weights;
 const referenceScores=screenEmissionSupport(junction.columns,input.observed,junction,mandatory,input.positionWeights,junction.models.map((_,i)=>[i]),hypotheses,{refit:false,templateResidual}).scores;
 return {...support,referenceScores};
}
/** A class dominated by its null in BOTH independent prediction folds has no
 * predictive reason to remain in the point model. Keep one-sided uncertainty;
 * never use a user display threshold or remove measured observations. */
export function rejectNonpredictiveJunctions(input:JunctionInput,junction:JunctionFit,evidence:Map<number,{foldGains:number[];lengthFoldGains?:number[]}>,precisionEvidence?:EmissionSupportEvidence[],lengthEvidence?:EmissionSupportEvidence[]):number[] {
 const tolerance=1e-8*input.positions.length;
 const removed=[...evidence].flatMap(([i,e])=>{
  if(isWildType(junction.models[i])||e.foldGains.length!==2)return [];
  const precise=(scores:EmissionSupportEvidence[]|undefined)=>{const v=scores?.find(v=>v.group.includes(i));return !!v&&v.score>=v.threshold;},length=e.lengthFoldGains||e.foldGains;
  const lengthRejected=length.length===2&&length.every(g=>g<=tolerance)&&(length.some(g=>g< -tolerance)||!precise(lengthEvidence));
  const attributionRejected=!precise(precisionEvidence)&&e.foldGains.every(g=>g<=tolerance);
  return lengthRejected||attributionRejected?[i]:[];
 });
 if(removed.length)refitSelected(input,junction,junction.models.flatMap((m,i)=>isWildType(m)||junction.weights[i]>1e-10&&!removed.includes(i)?[i]:[]));
 return removed;
}
