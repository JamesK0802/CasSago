import type { SangerAnalysisResult } from '../../models/sanger-analysis';
import type { SangerLinearOptimizer } from './sanger-discrete-phase';
import { verifiedZeroDeparturePoint } from './sanger-initial-marginal-fit';
import { MarginalLinearSolver,type LinearTerm } from './sanger-marginal-lp';
import type { MarginalModel } from './sanger-marginal-model';
import { marginalModelDiagnostics } from './sanger-marginal-state';
import { selectMarginalRepresentative } from './sanger-representative-selection';
/** Shared fit stage, independent of sample/ref/target/allele count.
 * The caller owns solver.dispose(); reporting does not participate in fitting. */
export function fitMarginalPoint(model:MarginalModel,optimizer:SangerLinearOptimizer,tick:(stage:string)=>void,progress:(stage:string)=>void){
 const {program,drafts,phaseSeed,noise,allDepartures,states,shares,allEdges}=model;
 tick('Fitting the shared global sequence-flow model');
 const solver=new MarginalLinearSolver(program,optimizer);let point:number[]=[],budget=0,representativeSelection:NonNullable<SangerAnalysisResult['diagnostics']>['representativeSelection'];
 try{
  progress('Fitting the shared global sequence-flow model · checking zero-departure witness');const zeroPoint=verifiedZeroDeparturePoint(program,phaseSeed,optimizer);progress('Fitting the shared global sequence-flow model · verifying full-model optimum');let fitted:{x:number[];value:number}|undefined=zeroPoint,shadow:number[]=[],planes:{terms:LinearTerm[];upper:number}[]=[],minimumLower:number|undefined=zeroPoint?0:undefined;
  if(!fitted){const direct=new MarginalLinearSolver({...program,localBlocks:undefined},optimizer);try{direct.setTimeSlice(600);fitted=direct.solveToCompletion();shadow=direct.rowDuals();planes=direct.supportPlanes();minimumLower=direct.objectiveLowerBound();}finally{direct.dispose();}}
  marginalModelDiagnostics.primaryFit?.({point:fitted.x,value:fitted.value,lower:minimumLower,shadow,planes});
  // Registering an existing genomic locus at a shifted read cycle must not
  // change the nominal noise allowance or representative penalty. Additional
  // phase coverage reclassifies measured rows; it is not a new prior budget.
  point=fitted.x;for(const draft of drafts)draft.primaryPoint=point;budget=fitted.value+3*noise*Math.sqrt(4*drafts.reduce((n,d)=>n+(d.mutationPriorOpportunities??d.positions.length),0));const budgetRow=program.rows.length;solver.addRow(allDepartures,0,budget);program.dualHint={shadow,budgetRow,planes};
  const context={substitutionOpportunities:drafts.reduce((n,d)=>n+(d.mutationPriorOpportunities??d.positions.length),0)/drafts.length,noise,phaseSeed,states,shares,edges:allEdges},edited=drafts.flatMap(d=>d.groups.filter(g=>g.size!==0).flatMap(g=>g.mass));
  tick('Selecting parsimonious sequence-flow representative');
  const baseline=selectMarginalRepresentative(program,point,edited,optimizer,context);point=baseline.point;
  tick('Refining base proportions');
  const selected=selectMarginalRepresentative(program,point,edited,optimizer,{...context,diagnostics:marginalModelDiagnostics.representativeProgram,measurements:drafts.flatMap(d=>d.measurements||[]),baseFlowAnchor:drafts.flatMap(d=>d.groups.flatMap(g=>[...g.bases.values()].map(forms=>({forms,mass:g.mass,share:shares[g.reference]}))))});
  point=selected.point;representativeSelection=selected.summary;marginalModelDiagnostics.representativeFit?.(selected.summary);
  return {solver,point,budget,representativeSelection};
 }catch(error){solver.dispose();throw error;}
}
