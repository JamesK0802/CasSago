import {SangerAnalysisResult} from '../../models/sanger-analysis';
import {PooledAnalysisInput,preparePooledSanger} from './sanger-pooled';
import {JunctionInput,JunctionRead,JunctionFit} from './sanger-junction-mixture';
import {solveJunctions,JunctionSolution} from './sanger-junction-solution';
import {mapPathModel,pathColumn} from './sanger-junction-paths';
import {fitEmissionMixture} from './sanger-emission';
/** Two independent fluorescence calibrations, one shared physical distribution.
 * No per-read mode, clipping parameter, or second analytical engine. */
export function analyzePairedSanger(inputs:PooledAnalysisInput[],progress:(stage:string)=>void=()=>{}):SangerAnalysisResult {
 if(inputs.length!==2)throw Error('Provide the primary read pair and one opposite-direction read pair.');
 const prepared=inputs.map((input,i)=>preparePooledSanger(input,stage=>progress(`Read ${i+1}/2: ${stage}`))),primary=prepared[0];
 inputs.forEach(input=>{if(input.amplicons.length!==inputs[0].amplicons.length||input.amplicons.some((a,h)=>a.name!==inputs[0].amplicons[h].name||a.sequence!==inputs[0].amplicons[h].sequence||JSON.stringify(a.cuts||[a.cut])!==JSON.stringify(inputs[0].amplicons[h].cuts||[inputs[0].amplicons[h].cut])))throw Error('Both directions must use the same reference products and target positions.');});
 let offset=0;
 const reads:JunctionRead[]=prepared.map((p,i)=>{
  const positions=p.fitInput.positions,rowIds=positions.map((_,j)=>offset+j);offset+=positions.length;
  const mappings=p.report.amplicons.map((site,h)=>{const first=primary.report.amplicons[h],length=inputs[0].amplicons[h].sequence.length,origin=(x:typeof site)=>x.direction==='forward'?x.primerEnd:length-x.primerStart+1,reverse=site.direction!==first.direction;
   if(i===1&&!reverse)throw Error('The additional read pair must use the opposite sequencing primer in every reference.');
   return{reverse,offset:reverse?length-origin(first)-origin(site):origin(first)-origin(site)};
  });
  return{positions,rowIds,references:p.fitInput.references,mappings,lag:p.fitInput.lag,theta:p.fitInput.theta};
 });
 const input:JunctionInput={...primary.fitInput,reads,positions:reads.flatMap(r=>r.rowIds),observed:prepared.flatMap(p=>p.fitInput.observed),positionWeights:prepared.flatMap(p=>p.fitInput.positionWeights||p.fitInput.positions.map(()=>1))};
 progress('Jointly fitting the same physical edits to both read directions');
 const solution=solveJunctions(input,fit=>prepared.flatMap(p=>p.templateError(fit)||p.fitInput.observed.map(()=>0)));
 solution.coveragePass=inputs[0].amplicons.every((a,h)=>(a.cuts||[a.cut]).every((_,t)=>prepared.some(p=>{const cut=p.fitInput.references[h].cuts![t],positions=p.fitInput.positions;return positions[0]<=cut-10&&positions.at(-1)!>=cut+29;})));
 const views=prepared.map((p,i)=>{
  const models=solution.junction.models.map(m=>mapPathModel(m,reads[i].mappings[m.reference])),columns=solution.junction.predictColumns(reads[i].rowIds);
  const junction:JunctionFit={...solution.junction,models,columns,predictColumns:positions=>models.map(m=>pathColumn(p.fitInput.references[m.reference],m,positions,p.fitInput.lag,p.fitInput.theta))};
  const projected:JunctionSolution={...solution,junction};return p.finish(projected);
 });
 const result=views[0],bases=views.reduce((n,v)=>n+v.mixture!.bases,0),error=views.reduce((n,v)=>n+v.mixture!.rmse**2*4*v.mixture!.bases,0),variance=prepared.reduce((n,p)=>n+p.rawVariance,0);
 result.reads=views.map((v,i)=>({primer:inputs[i].primer,sampleName:inputs[i].sample.name,controlName:inputs[i].control.name,preflight:v.pooled!.preflight,sampleTrim:v.pooled!.sampleTrim,controlTrim:v.pooled!.controlTrim,r2:v.mixture!.r2,rmse:v.mixture!.rmse,bases:v.mixture!.bases,fitStart:v.mixture!.fitStart,fitEnd:v.mixture!.fitEnd}));
 // Disagreement is evaluated on each read's observable first-target marginal,
 // not arbitrary one-way allocations of compound molecular paths.
 let agreement=0;
 prepared.forEach((p,i)=>{
  const columns=solution.junction.predictColumns(reads[i].rowIds),independent=fitEmissionMixture(columns,p.fitInput.observed,p.fitInput.positionWeights).weights;
  p.fitInput.references.forEach((ref,h)=>{const t=(ref.cuts||[ref.cut]).indexOf(Math.min(...(ref.cuts||[ref.cut]))),ids=solution.junction.models.flatMap((m,j)=>m.reference===h?[j]:[]),share=(w:number[])=>ids.reduce((s,j)=>s+w[j],0),a=share(independent),b=share(solution.junction.weights);
   agreement=Math.max(agreement,100*Math.abs(a-b));if(Math.min(a,b)<.02)return;
   const sizes=[...new Set(ids.map(j=>solution.junction.models[j].sizes?.[t]||0))];agreement=Math.max(agreement,50*sizes.reduce((sum,size)=>sum+Math.abs(ids.reduce((s,j)=>s+((solution.junction.models[j].sizes?.[t]||0)===size?independent[j]/a-solution.junction.weights[j]/b:0),0)),0));
  });
 });
 result.plots=views.flatMap((v,i)=>v.plots.map(p=>({...p,name:`Read ${i+1} · ${p.name}`})));
 result.warnings=[...new Set(views.flatMap((v,i)=>v.warnings.map(w=>`Read ${i+1}: ${w}`)))];
 if(!solution.coveragePass)result.warnings.push('Some targets lack sufficient flanking signal in both directions; numerical estimates remain displayed.');
 if(agreement>10)result.warnings.push('Opposite-direction reads disagree by more than 10 pp in reference shares or first-target indel distributions. Check that both pairs came from the same edited sample and WT control.');
 result.mixture={...result.mixture!,r2:variance?1-error/variance:0,rmse:Math.sqrt(error/(4*bases)),bases,passesChecks:views.every(v=>v.mixture!.passesChecks)&&agreement<=10};
 if(result.joint)result.joint.reconstructionR2=result.mixture.r2;
 result.pooled={...result.pooled!,controlR2:Math.min(...views.map(v=>v.pooled!.controlR2))};
 if(result.diagnostics)result.diagnostics.readAgreementPp=agreement;
 for(const a of result.annotations||[])a.estimable=result.mixture.passesChecks&&a.share>=.02;
 return result;
}
