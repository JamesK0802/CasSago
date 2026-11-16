import type {SangerAnalysisResult,SangerMarginalGroup} from '../../models/sanger-analysis';
import {marginalSessions} from './sanger-marginal-state';
import {MarginalLinearSolver,type LinearTerm} from './sanger-marginal-lp';
import {zeroLengthCompositionBounds} from './sanger-zero-length-composition';
import {localBoundCertificates} from './sanger-local-bound-certificate';
import {nonnegativeBudgetEnvelope} from './sanger-nonnegative-budget-envelope';
import {withinLengthDisplay} from './sanger-marginal-display';

interface Query {index:number;group:SangerMarginalGroup;mass:[number,number];familyMass:[number,number];seconds:number}
const queries=new WeakMap<SangerAnalysisResult,Map<string,Query>>();
/** On-demand report query only. It reads the retained original joint model;
 * fitting, representative selection and all other groups stay unchanged.
 * Return no partial envelope as a completed range. A time slice preserves the
 * coordinate cursor; it is not a limit on the model or accepted alternatives.
 */
export function marginalGroupBounds(result:SangerAnalysisResult,groupId:string,progress:(s:string)=>void=()=>{},timeSliceSeconds=20){
 const started=performance.now(),session=marginalSessions.get(result);
 if(!session)throw Error('Analyze the original traces before requesting group ranges.');
 const parent=result.marginalGroups?.find(g=>g.id===groupId||g.displayComponents?.some(c=>c.id===groupId)),child=parent?.displayComponents?.find(c=>c.id===groupId),marker=child?.selection.marker;
 if(!parent||child&&!marker)throw Error('Select a length family or a base-marker group for feasible range calculation.');
 if(child?child.rangeStatus==='bounded':parent.status==='bounded')return{group:parent,complete:true,seconds:0,completedCoordinates:parent.bases.length,totalCoordinates:parent.bases.length};
 const target=session.targets.find(t=>t.groups.some(g=>g.id===parent.id))!,g=target.groups.find(g=>g.id===parent.id)!,ref=session.refs[g.reference],share=session.shares[g.reference],childEstimates=new Map(child?.bases.map(p=>[p.position,p.estimate]));
 let map=queries.get(result);if(!map){map=new Map();queries.set(result,map);}let query=map.get(groupId);
 // The flat solver is only the fallback for objectives without a local dual
 // certificate. Do not run thousands of decomposition rounds per nucleotide.
 const solver=new MarginalLinearSolver({...target.program,localBlocks:undefined},session.optimizer),local=localBoundCertificates(target.program,target.point,target.groups,session.optimizer),budgetUpper=nonnegativeBudgetEnvelope(target.program);
 try{
  if(!query){const familyMass=solver.bounds(g.mass),mass=marker?solver.bounds(marker.bases.flatMap(b=>g.bases.get(marker.position)![b])):familyMass;query={index:0,mass,familyMass,seconds:0,group:{...parent,bases:[],fractionRange:[Math.min(1,mass[0]/share),Math.min(1,mass[1]/share)],displayComponents:parent.displayComponents}};map.set(groupId,query);}
  const positions=[...g.bases],evaluate=(terms:LinearTerm[])=>terms.reduce((n,[i,w])=>n+w*target.point[i],0);
  for(;query.index<positions.length;query.index++){
   if(performance.now()-started>=1000*timeSliceSeconds)break;
   const[q,forms]=positions[query.index];progress('Bounding selected group · '+(query.index+1)+' / '+positions.length);
   const position=ref.reverse?ref.input.sequence.length-ref.offset-q-1:ref.offset+q,reported=childEstimates.get(position),estimate=reported?(ref.reverse?reported.slice().reverse():reported.slice()):forms.map(f=>g.fraction?Math.max(0,evaluate(f))/g.fraction:0),lower=[0,0,0,0],upper=[1,1,1,1],dependent=estimate.indexOf(Math.max(...estimate));
   for(let b=0;b<4;b++){
    if(b===dependent)continue;
    if(marker&&q===marker.position&&!marker.bases.includes(b)){lower[b]=upper[b]=0;continue;}
    if(!forms[b].length){lower[b]=upper[b]=0;continue;}
    const lo=local.outerBound(q,forms[b],false)??solver.solve(forms[b]).value,hi=local.outerBound(q,forms[b],true)??solver.solve(forms[b],true).value;
    // For a marker-defined mixture, another coordinate's linkage is not
    // observed. Fréchet bounds retain every compatible coupling rather than
    // treating copied fitted percentages as phased molecular evidence.
    const numeratorLo=marker&&q!==marker.position?Math.max(0,lo+query.mass[0]-query.familyMass[1]-1e-7):Math.max(0,lo-1e-7),numeratorHi=marker?Math.min(hi+1e-7,query.mass[1]):hi+1e-7;
    lower[b]=Math.max(0,Math.min(1,numeratorLo/Math.max(1e-30,query.mass[1])));upper[b]=query.mass[0]>1e-8?Math.min(1,numeratorHi/query.mass[0]):1;
    lower[b]=Math.min(lower[b],estimate[b]);upper[b]=Math.max(upper[b],estimate[b]);
   }
   lower[dependent]=Math.min(estimate[dependent],Math.max(0,1-upper.reduce((n,v,b)=>n+(b===dependent?0:v),0)));
   upper[dependent]=Math.max(estimate[dependent],Math.min(1,1-lower.reduce((n,v,b)=>n+(b===dependent?0:v),0)));
   if(query.mass[0]>1e-8)for(let b=0;b<4;b++){
    const hi=budgetUpper(forms[b]),other=budgetUpper(forms.flatMap((f,a)=>a===b?[]:Array.from(f)));
    if(hi!==undefined)upper[b]=Math.min(upper[b],Math.max(estimate[b],Math.min(1,(hi+1e-7)/query.mass[0])));
    if(other!==undefined)lower[b]=Math.max(lower[b],Math.min(estimate[b],Math.max(0,1-(other+1e-7)/query.mass[0])));
   }
   query.group.bases.push({position,estimate:ref.reverse?estimate.reverse():estimate,lower:ref.reverse?lower.reverse():lower,upper:ref.reverse?upper.reverse():upper});
  }
  const complete=query.index===positions.length;query.seconds+=(performance.now()-started)/1000;
  if(complete){
   query.group.status='bounded';query.group.boundMethod='certified-outer-envelope';
   if(!g.size&&!child){const budgetRow=target.program.dualHint?.budgetRow,budget=budgetRow===undefined?0:target.program.rows[budgetRow].upper;Object.assign(query.group,zeroLengthCompositionBounds(query.mass,budget,share));}
   if(child){child.rangeStatus='bounded';child.bases=query.group.bases;child.fractionRange=query.group.fractionRange;}
   else{
   // Enrich only this family's report; do not refit or replace its point.
   Object.assign(parent,query.group);g.fractionRange=query.mass;
   // Previously this separation ran on EVERY family before a result. Keep
   // its supported profiles available as detail for this requested family.
   const detailStarted=performance.now();
   try{const parts=withinLengthDisplay(session,target,g,parent);if(parts)parent.displayComponents=parts;}
   catch{result.warnings.push('Within-length separation did not complete for '+groupId+'; the full parent mixture is retained.');}
   query.seconds+=(performance.now()-detailStarted)/1000;
   }
   map.delete(groupId);
  }
  return{group:complete?parent:undefined,complete,seconds:query.seconds,completedCoordinates:query.index,totalCoordinates:positions.length};
 }finally{local.dispose();solver.dispose();}
}
