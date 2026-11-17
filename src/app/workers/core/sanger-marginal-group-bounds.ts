import type {SangerAnalysisResult,SangerMarginalGroup} from '../../models/sanger-analysis';
import {marginalSessions} from './sanger-marginal-state';
import {MarginalLinearSolver,type LinearTerm} from './sanger-marginal-lp';
import {zeroLengthCompositionBounds} from './sanger-zero-length-composition';
import {nonnegativeBudgetEnvelope} from './sanger-nonnegative-budget-envelope';
import {tightenBaseBudgetBounds} from './sanger-base-budget-bounds';
import {jointBoundCertificates} from './sanger-joint-bound-certificate';
import {conditionalRatioEnvelope} from './sanger-conditional-ratio-envelope';
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
 // One warm joint linear proposal retains all coordinates. Only an incomplete
 // proposal falls back to a certified full original-model solve.
 const solver=new MarginalLinearSolver({...target.program,localBlocks:undefined},session.optimizer),joint=jointBoundCertificates(target.program,session.optimizer),budgetUpper=nonnegativeBudgetEnvelope(target.program);
 try{
  if(!query){const familyMass=solver.bounds(g.mass),mass=marker?solver.bounds(marker.bases.flatMap(b=>g.bases.get(marker.position)![b])):familyMass;query={index:0,mass,familyMass,seconds:0,group:{...parent,bases:[],fractionRange:[Math.min(1,mass[0]/share),Math.min(1,mass[1]/share)],displayComponents:parent.displayComponents}};map.set(groupId,query);}
  const positions=[...g.bases],evaluate=(terms:LinearTerm[])=>terms.reduce((n,[i,w])=>n+w*target.point[i],0);
  for(;query.index<positions.length;query.index++){
   if(performance.now()-started>=1000*timeSliceSeconds)break;
   const[q,forms]=positions[query.index];progress('Bounding selected group · '+(query.index+1)+' / '+positions.length);
   const position=ref.reverse?ref.input.sequence.length-ref.offset-q-1:ref.offset+q,reported=childEstimates.get(position),estimate=reported?(ref.reverse?reported.slice().reverse():reported.slice()):forms.map(f=>g.fraction?Math.max(0,evaluate(f))/g.fraction:0),lower=[0,0,0,0],upper=[1,1,1,1];
   const denominator=marker?g.bases.get(marker.position)!.flatMap((f,b)=>marker.bases.includes(b)?Array.from(f):[]):g.mass;
   for(let b=0;b<4;b++){
    if(marker&&q===marker.position&&!marker.bases.includes(b)||!forms[b].length){lower[b]=upper[b]=0;continue;}
    // Marker membership at another cycle is unobserved. Frechet inequalities
    // retain every compatible linkage rather than assigning a phased allele.
    const numeratorLower=marker&&q!==marker.position?[...forms[b],...denominator,...g.mass.map(([i,w])=>[i,-w]as LinearTerm)]:forms[b];
    const lo=joint.bound(numeratorLower,false)??solver.solveToCompletion(numeratorLower).value,hi=joint.bound(forms[b],true)??solver.solveToCompletion(forms[b],true).value;
    lower[b]=Math.min(estimate[b],Math.max(0,(lo-1e-7)/Math.max(1e-30,query.mass[1])));
    upper[b]=Math.max(estimate[b],query.mass[0]>1e-8?Math.min(1,(hi+1e-7)/query.mass[0]):1);
    [lower[b],upper[b]]=conditionalRatioEnvelope(numeratorLower,forms[b],denominator,query.mass,estimate[b],[lower[b],upper[b]],(terms,maximum)=>joint.bound(terms,maximum)??solver.solveToCompletion(terms,maximum).value);
   }
   // Every conditional four-base profile sums to one. Endpoints remain
   // marginal outer bounds; independent extrema are not a joint genotype.
   for(let b=0;b<4;b++){
    lower[b]=Math.max(lower[b],Math.min(estimate[b],Math.max(0,1-upper.reduce((n,v,a)=>n+(a===b?0:v),0))));
    upper[b]=Math.min(upper[b],Math.max(estimate[b],Math.min(1,1-lower.reduce((n,v,a)=>n+(a===b?0:v),0))));
   }
   tightenBaseBudgetBounds(forms,query.mass[0],estimate,lower,upper,budgetUpper);
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
 }finally{joint.dispose();solver.dispose();}
}
