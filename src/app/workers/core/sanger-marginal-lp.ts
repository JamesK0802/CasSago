import {VerifiedConicSolver} from './sanger-conic-solver';
import {zeroCostProjection} from './sanger-zero-cost-projection';
import {evaluateQuadratic,quadraticCut,type QuadraticEnvelope} from './sanger-quadratic-envelope';
import {FactorizedMarginalSolver} from './sanger-factorized-lp';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
export type LinearTerm=[number,number];
export interface MarginalProgram {/** Coefficients are immutable; numeric row/box changes still invalidate caches. */immutableRows?:boolean;quadraticCost?:LinearTerm[][];dualHint?:{shadow:number[];budgetRow:number;planes?:{terms:LinearTerm[];upper:number}[]};zeroCostSubstitutions?:{column:number;terms:LinearTerm[]}[];quadratics?:QuadraticEnvelope[];localBlocks?:number[][];originPricing?:{columns:number[];seed:number[]};upper:number[];cost:number[];rows:{terms:LinearTerm[];lower:number;upper:number}[]}
export const LP_INFINITY=1e30;
export class OptimizationPending extends Error {constructor(){super('Optimization is still in progress; no solution was accepted.');}}
const formKey=(terms:LinearTerm[])=>{const m=new Map<number,number>();for(const[i,w]of terms)m.set(i,(m.get(i)||0)+w);return[...m].filter(([,w])=>w!==0).sort((a,b)=>a[0]-b[0]).map(([i,w])=>i+':'+w).join(',');};
export function nativeMarginalProgram(p:MarginalProgram){const starts=[0],indices:number[]=[],values:number[]=[];for(const r of p.rows){for(const[i,v]of r.terms){if(!v)continue;indices.push(i);values.push(v);}starts.push(indices.length);}return{numCols:p.upper.length,numRows:p.rows.length,colCost:p.cost,colLower:p.upper.map(()=>0),colUpper:p.upper,integrality:p.upper.map(()=>0),rowLower:p.rows.map(r=>r.lower),rowUpper:p.rows.map(r=>r.upper),matrix:{format:'csr',numCols:p.upper.length,numRows:p.rows.length,starts,indices,values}};}
export function marginalViolation(p:MarginalProgram,x:number[]){if(x.length!==p.upper.length)return Infinity;let violation=0;for(let i=0;i<x.length;i++){if(!Number.isFinite(x[i]))return Infinity;violation=Math.max(violation,-x[i],x[i]-p.upper[i]);}for(const r of p.rows){const v=r.terms.reduce((n,[i,w])=>n+w*x[i],0);violation=Math.max(violation,r.lower-v,v-r.upper);}for(const q of p.quadratics||[])violation=Math.max(violation,evaluateQuadratic(q,x)-q.upper);return violation;}
export function validMarginalSolution(p:MarginalProgram,x:number[],tolerance=1e-7){return marginalViolation(p,x)<=tolerance;}
/** Convex tangent + arbitrary row multipliers + complete box minimization.
 * Expanding every original row and box by tolerance makes an infeasibility
 * certificate conservative even under the public primal acceptance tolerance. */
export function quadraticTangentLowerBound(p:MarginalProgram,x:number[],shadow:number[],tolerance=0){
 const gradient=p.cost.slice();let constant=0,magnitude=0;
 for(const form of p.quadraticCost||[]){const v=form.reduce((n,[i,w])=>n+w*x[i],0);constant-=v*v;magnitude+=v*v;for(const[i,w]of form)gradient[i]+=2*w*v;}
 let lower=constant;
 for(let k=0;k<p.rows.length;k++){const r=p.rows[k];let y=shadow[k]||0;if(!Number.isFinite(y)||y>0&&r.lower<=-LP_INFINITY||y<0&&r.upper>=LP_INFINITY)y=0;if(!y)continue;const v=y*(y>0?r.lower-tolerance:r.upper+tolerance);lower+=v;magnitude+=Math.abs(v);for(const[i,w]of r.terms){gradient[i]-=y*w;magnitude+=Math.abs(y*w*x[i]);}}
 for(let i=0;i<gradient.length;i++){const v=gradient[i]>=0?-tolerance*gradient[i]:(p.upper[i]+tolerance)*gradient[i];lower+=v;magnitude+=Math.abs(v);}
 return lower-128*Number.EPSILON*(1+magnitude);
}
/** One continuous LP, solved by certified column activation when origin hints
 * are supplied. Hints never change the public program or discard a column.
 * Reduced costs certify the original optimum; an unsuccessful restricted
 * solve restores every column before making any infeasibility decision. */
export class MarginalLinearSolver {
 private linearLower?:number;
 private lowerObjectiveIsCurrent=false;
 private native:any;
 private conic?:VerifiedConicSolver;
 private factorized?:FactorizedMarginalSolver;
 private fullCompetitor?:MarginalLinearSolver;
 private decomposed=false;
 private sliceSeconds?:number;
 private initialSliceSeconds?:number;
 private factorizationFailed=false;
 private directProbeAttempted=false;
 private directProbe={attempts:0,certified:0,seconds:0};
 private faceAttempted=false;
 private zeroFaceProjection?:{variables:number;rows:number;seconds:number;verified:boolean;failure?:string};
 private zeroFaceCertified=false;
 private winner?:'flat'|'factorized';
 private optimizer:SangerLinearOptimizer;
 private factorizationFallbackReason?:string;
 get diagnostics(){return{directProbe:this.directProbe,spectralFeasibility:this.spectralFeasibility,conicFallbackReason:this.conicFallbackReason,conic:this.conic?.stats,zeroFaceProjection:this.fullCompetitor?.zeroFaceProjection||this.zeroFaceProjection,selected:this.conic?'conic' as const:this.winner,factorization:this.factorized?.stats,factorizationFallbackReason:this.factorizationFallbackReason,pricing:this.fullCompetitor?.pricingStats||this.pricingStats};}
 private conicFallbackReason?:string;
 private parentFallbackBudget?:number;
 private coneRecoveryStarted?:number;
 private spectralFeasibility={attempts:0,feasible:0,infeasible:0,fallbacks:0,seconds:0,minimum:undefined as number|undefined,lower:undefined as number|undefined,originalViolation:undefined as number|undefined};
 private costs:number[];
 private point?:number[];
 private equalities?:Map<string,{value:number;row:number}>;
 private shortcut?:{row?:number;multiplier?:number};
 private inactive=new Set<number>();
 private incidence=new Map<number,LinearTerm[]>();
 readonly pricingStats={rounds:0,activated:0,fullRestores:0};
 constructor(readonly program:MarginalProgram,optimizer:SangerLinearOptimizer,representation:'auto'|'conic'='auto'){
  if(!optimizer.createModel)throw Error('Initialize the common optimizer.');
  this.optimizer=optimizer;this.costs=program.cost.slice();
  if(program.quadraticCost?.length){if(!optimizer.conic)throw Error('A verified convex quadratic optimizer is required.');this.conic=new VerifiedConicSolver(program,optimizer.conic);return;}
  if(representation!=='conic'&&program.localBlocks?.some(v=>v.length)){this.decomposed=true;return;}
  if(optimizer.conic&&(representation==='conic'||program.quadratics?.length)){this.conic=new VerifiedConicSolver(program,optimizer.conic);return;}
  this.initializeNative();
 }
 private initializeNative(){const program=this.program,optimizer=this.optimizer;this.costs=program.cost.slice();this.inactive.clear();this.incidence.clear();
  this.native=optimizer.createModel!();
  if(program.originPricing){const seed=new Set(program.originPricing.seed);for(const i of program.originPricing.columns){this.incidence.set(i,[]);if(program.upper[i]>0&&!seed.has(i))this.inactive.add(i);}program.rows.forEach((r,k)=>{for(const[i,w]of r.terms)this.incidence.get(i)?.push([k,w]);});}
  this.native.options.set({output_flag:false,solver:'simplex',simplex_strategy:4,time_limit:300,primal_feasibility_tolerance:program.quadratics?.length?1e-10:1e-8,dual_feasibility_tolerance:program.quadratics?.length?1e-10:1e-8});
  const model=nativeMarginalProgram(program);model.colUpper=program.upper.map((u,i)=>this.inactive.has(i)?0:u);this.native.passModel(model);
 }
 private makeFactorized(){return new FactorizedMarginalSolver(this.program,p=>{const solver=new MarginalLinearSolver({...p,localBlocks:undefined},this.optimizer);solver.setTimeSlice(p.quadratics?.length?30:.25);solver.parentFallbackBudget=120;return solver;});}
 setTimeSlice(seconds:number){this.initialSliceSeconds=seconds;this.sliceSeconds=seconds;if(this.fullCompetitor)this.fullCompetitor.setTimeSlice(seconds);else if(this.conic)this.conic.setTimeLimit(seconds);else if(this.native)this.native.options.set({time_limit:seconds});}
 private pending():never{if(this.sliceSeconds!==undefined&&this.sliceSeconds>0){this.sliceSeconds*=2;if(this.conic)this.conic.setTimeLimit(this.sliceSeconds);else this.native.options.set({time_limit:this.sliceSeconds});}throw new OptimizationPending();}
 private tryZeroCostFace():number[]|undefined{
  if(this.costs.some(v=>v<0)||!this.costs.some(v=>v>0))return;
  const projection=zeroCostProjection(this.program);if(projection){const started=performance.now(),solver=new MarginalLinearSolver(projection.program,this.optimizer);this.zeroFaceProjection={variables:projection.program.upper.length,rows:projection.program.rows.length,seconds:0,verified:false};try{solver.setTimeSlice(30);solver.native.options.set({solver:'ipm',run_crossover:'on'});const answer=projection.lift(solver.solve().x);if(validMarginalSolution(this.program,answer)){this.zeroFaceProjection.verified=true;return answer;}}catch(error){this.zeroFaceProjection.failure=String(error);}finally{this.zeroFaceProjection.seconds=(performance.now()-started)/1000;solver.dispose();}return;}
  if(this.inactive.size){this.pricingStats.fullRestores++;this.activate([...this.inactive]);}
  const indices=this.costs.flatMap((c,i)=>c>0&&this.program.upper[i]>0?[i]:[]);if(!indices.length)return;
  this.native.changeColsBounds({kind:'set',indices},indices.map(()=>0),indices.map(()=>0));this.native.options.set({solver:'ipm',run_crossover:'on'});const started=performance.now();let answer:number[]|undefined;
  try{for(;;){const remaining=10-(performance.now()-started)/1000;if(remaining<=0)break;this.native.options.set({time_limit:remaining});this.native.zeroAllClocks();const status=this.native.run().modelStatus;if(status!==7)break;const x=Array.from(this.native.getSolution().colValue as ArrayLike<number>);if(marginalViolation({...this.program,quadratics:undefined},x)>1e-7)break;const bad=(this.program.quadratics||[]).filter(q=>evaluateQuadratic(q,x)>q.upper+5e-8);if(!bad.length){answer=x;break;}for(const q of bad){const cut=quadraticCut(q,x);this.addRow(cut.terms,-LP_INFINITY,cut.upper);}}}catch{/* A failed initialization never changes or certifies the full problem. */}
  finally{this.native.changeColsBounds({kind:'set',indices},indices.map(()=>0),indices.map(i=>this.program.upper[i]));this.native.options.set({solver:'simplex',time_limit:this.sliceSeconds??300});}
  return answer;
 }
 private activate(indices:number[]){
  if(!indices.length)return;indices=indices.slice().sort((a,b)=>a-b);
  this.native.changeColsBounds({kind:'set',indices},indices.map(()=>0),indices.map(i=>this.program.upper[i]));
  for(const i of indices)this.inactive.delete(i);this.pricingStats.activated+=indices.length;
 }
 /** Completion mode extends execution budgets only; original primal and
  * global dual tolerances and every scientific constraint remain unchanged. */
 solveToCompletion(terms?:LinearTerm[],maximize=false):{x:number[];value:number}{
  if(this.sliceSeconds===undefined||this.sliceSeconds<=0)this.setTimeSlice(300);
  for(;;){try{return this.solve(terms,maximize);}catch(error){if(!(error instanceof OptimizationPending))throw error;}}
 }
 solve(terms?:LinearTerm[],maximize=false):{x:number[];value:number}{
  this.lowerObjectiveIsCurrent=false;this.linearLower=undefined;this.shortcut=undefined;if(this.conicFallbackReason&&this.coneRecoveryStarted===undefined)this.coneRecoveryStarted=performance.now();
  if(terms&&this.point&&!this.program.quadraticCost?.length){
   if(!this.equalities){this.equalities=new Map();for(const[row,r]of this.program.rows.entries())if(r.lower===r.upper&&Math.abs(r.lower)<LP_INFINITY)this.equalities.set(formKey(r.terms),{value:r.lower,row});}
   const fixed=this.equalities.get(formKey(terms));if(fixed!==undefined&&validMarginalSolution(this.program,this.point)){this.linearLower=(maximize?-1:1)*fixed.value;this.shortcut={row:fixed.row,multiplier:maximize?-1:1};return{x:this.point.slice(),value:fixed.value};}
   if(terms.every(([i,w])=>w===0||this.program.upper[i]===0)&&validMarginalSolution(this.program,this.point)){this.linearLower=0;this.shortcut={};return{x:this.point.slice(),value:0};}

  }
  if(this.conic){try{const answer=this.conic.solve(terms,maximize);this.lowerObjectiveIsCurrent=true;this.point=answer.x;if(this.initialSliceSeconds!==undefined){this.sliceSeconds=this.initialSliceSeconds;this.conic.setTimeLimit(this.sliceSeconds);}return answer;}catch(error){if(error instanceof OptimizationPending)this.pending();if(!this.program.quadraticCost?.length&&String(error).includes('Conic optimization returned no dual witness: unknown')){this.conicFallbackReason=String(error);this.conic.dispose();this.conic=undefined;this.initializeNative();this.winner='flat';if(this.sliceSeconds!==undefined)this.native.options.set({time_limit:this.sliceSeconds});return this.solve(terms,maximize);}throw error;}}
  if(this.decomposed){
   // Sequential numeric probe of the identical full program. No decomposition
   // has been allocated yet; unfinished/uncertified proposals are disposed.
   // The short limit selects a representation, never a scientific result.
   if(!this.directProbeAttempted&&!this.fullCompetitor&&this.optimizer.conic){
    this.directProbeAttempted=true;this.directProbe.attempts++;const started=performance.now();let probe:MarginalLinearSolver|undefined;
    try{probe=new MarginalLinearSolver({...this.program,localBlocks:undefined},this.optimizer,'conic');probe.conic!.setTimeLimit(2);const a=probe.conic!.solve(terms,maximize);probe.point=a.x;probe.lowerObjectiveIsCurrent=true;probe.setTimeSlice(this.initialSliceSeconds||600);this.fullCompetitor=probe;probe=undefined;this.factorizationFailed=true;this.factorizationFallbackReason='Identical full program certified during common sequential probe.';this.directProbe.certified++;this.point=a.x;this.winner='flat';return a;}
    catch{/* Retain the complete exact decomposition and unrestricted fallback. */}
    finally{this.directProbe.seconds+=(performance.now()-started)/1000;probe?.dispose();}
   }
   // Complete the common exact decomposition before constructing a second
   // full-model factorization. Racing both representations duplicates a
   // large model's memory while neither scientific hypothesis set changes.
   if(!this.factorizationFailed&&!this.fullCompetitor){try{
    this.factorized??=this.makeFactorized();const started=performance.now(),steps=this.factorized.steps(terms,maximize);
    try{for(;;){const next=steps.next();if(next.done){this.point=next.value.x;this.winner='factorized';return next.value;}
     // A representation-selection budget never accepts an approximate answer
     // or removes a hypothesis. Finish the identical full model below.
     if((performance.now()-started)/1000>15)throw Error('Exact decomposition exceeded its numerical selection budget; completing the unrestricted original parent.');
    }}finally{steps.return(undefined as never);}
   }catch(error){this.factorizationFallbackReason=String(error);this.factorizationFailed=true;this.factorized?.dispose();this.factorized=undefined;}}
   if(!this.fullCompetitor){this.fullCompetitor=new MarginalLinearSolver({...this.program,localBlocks:undefined},this.optimizer,'conic');this.fullCompetitor.setTimeSlice(this.initialSliceSeconds&&this.initialSliceSeconds>0?this.initialSliceSeconds:600);}
   for(;;){try{const answer=this.fullCompetitor.solve(terms,maximize);this.point=answer.x;this.winner='flat';return answer;}catch(error){if(!(error instanceof OptimizationPending))throw error;}}
  }
  if(terms){const next=this.costs.map(()=>0);for(const[i,v]of terms)next[i]+=(maximize?-1:1)*v;const indices=next.flatMap((v,i)=>v!==this.costs[i]?[i]:[]);if(indices.length)this.native.changeColsCost({kind:'set',indices},indices.map(i=>next[i]));this.costs=next;}
  if(terms)this.zeroFaceCertified=false;
  if(!terms&&!this.faceAttempted){this.faceAttempted=true;const x=this.tryZeroCostFace();if(x){this.point=x;this.zeroFaceCertified=true;this.linearLower=0;this.coneRecoveryStarted=undefined;return{x,value:0};}}
  const run=()=>{try{this.native.zeroAllClocks();const status=this.native.run().modelStatus,x=Array.from(this.native.getSolution().colValue as ArrayLike<number>);return{status,x,violation:marginalViolation({...this.program,quadratics:undefined},x),reason:undefined as string|undefined};}catch(error){return{status:-1,x:[]as number[],violation:Infinity,reason:String(error)};}};
  let solved:ReturnType<typeof run>;
  for(;;){
   if(this.conicFallbackReason&&this.parentFallbackBudget!==undefined&&this.coneRecoveryStarted!==undefined&&(performance.now()-this.coneRecoveryStarted)/1000>this.parentFallbackBudget)throw Error('Local numerical completion deferred to the unrestricted original parent optimizer.');
   solved=run();if(solved.status===13&&this.sliceSeconds!==undefined)this.pending();
   // A native numerical/runtime failure is not evidence that this small LP
   // needs a solve of its much larger parent. Retry the IDENTICAL full LP
   // using the independently checked conic backend, including every box.
   if(solved.status===-1&&this.optimizer.conic&&!this.conicFallbackReason){try{this.native.dispose();}catch{/* A failed native instance may already be unusable. */}this.native=undefined;this.conic=new VerifiedConicSolver(this.program,this.optimizer.conic);if(this.sliceSeconds!==undefined)this.conic.setTimeLimit(this.sliceSeconds);return this.solve(terms,maximize);}
   if(solved.status!==7||solved.violation>1e-7){
    if(this.inactive.size){this.pricingStats.fullRestores++;this.activate([...this.inactive]);this.native.clearSolver();continue;}
    this.native.options.set({primal_feasibility_tolerance:1e-10,dual_feasibility_tolerance:1e-10});
    try{this.native.clearSolver();solved=run();if(solved.status===13&&this.sliceSeconds!==undefined)this.pending();if(solved.status!==7&&solved.status!==8||solved.violation>1e-7&&solved.status===7){this.native.clearSolver();this.native.options.set({solver:'ipm'});solved=run();}}
    finally{this.native.options.set({solver:'simplex',primal_feasibility_tolerance:1e-8,dual_feasibility_tolerance:1e-8});}
   }
   if(solved.status===13&&this.sliceSeconds!==undefined)this.pending();
   // Polish an unverified local LP in the same cone representation before
   // abandoning its spatial decomposition. An optimal status with an
   // excessive original residual is not an acceptable local witness.
   if(solved.status!==8&&(solved.status!==7||solved.violation>1e-7)&&this.optimizer.conic&&!this.conicFallbackReason){this.native.dispose();this.native=undefined;this.conic=new VerifiedConicSolver(this.program,this.optimizer.conic);this.conic.setTimeLimit(Math.max(30,this.sliceSeconds||0));return this.solve(terms,maximize);}
   if(solved.status!==7||solved.violation>1e-7)throw Error('Marginal feasible-range optimization did not complete: '+solved.status+'; maximum violation '+solved.violation+(solved.reason?'; reason: '+solved.reason:''));
   const violated=(this.program.quadratics||[]).filter(q=>evaluateQuadratic(q,solved.x)>q.upper+5e-8);
   if(violated.length){for(const q of violated){const cut=quadraticCut(q,solved.x);this.addRow(cut.terms,-LP_INFINITY,cut.upper);}continue;}
   if(!this.inactive.size)break;
   // An attained original box lower bound is an exact global certificate.
   const value=this.costs.reduce((n,c,i)=>n+c*solved.x[i],0),box=this.costs.reduce((n,c,i)=>n+(c<0?c*this.program.upper[i]:0),0);
   if(value-box<=1e-8)break;
   const dual=Array.from(this.native.getSolution().rowDual as ArrayLike<number>);
   if(dual.length!==this.program.rows.length||dual.some(v=>!Number.isFinite(v))){this.activate([...this.inactive]);this.native.clearSolver();continue;}
   const priced=[...this.inactive].map(i=>{const reduced=this.costs[i]-(this.incidence.get(i)||[]).reduce((n,[k,w])=>n+dual[k]*w,0);return{i,reduced,gain:Math.max(0,-reduced)*this.program.upper[i]};}).filter(v=>v.gain>0).sort((a,b)=>b.gain-a.gain),gap=priced.reduce((n,v)=>n+v.gain,0);
   // All other original columns retain their bounds and the native optimal
   // dual certificate. This is the complete extra duality gap due to origins
   // temporarily held at zero, rather than a heuristic candidate cutoff.
   if(gap<=1e-8)break;
   this.pricingStats.rounds++;this.activate(priced.slice(0,64).map(v=>v.i));
  }
  const certificateDual=Array.from(this.native.getSolution().rowDual as ArrayLike<number>),certificateReduced=this.costs.slice();let certificateLower=0,certificateMagnitude=0;for(let k=0;k<this.program.rows.length;k++){const r=this.program.rows[k];let y=certificateDual[k]||0;if(!Number.isFinite(y)||y>0&&r.lower<=-LP_INFINITY||y<0&&r.upper>=LP_INFINITY)y=0;if(!y)continue;const v=y*(y>0?r.lower:r.upper);certificateLower+=v;certificateMagnitude+=Math.abs(v);for(const[i,w]of r.terms)certificateReduced[i]-=y*w;}for(let i=0;i<certificateReduced.length;i++)if(certificateReduced[i]<0){const v=certificateReduced[i]*this.program.upper[i];certificateLower+=v;certificateMagnitude+=Math.abs(v);}this.linearLower=certificateLower-128*Number.EPSILON*(1+certificateMagnitude);
  if(this.conicFallbackReason&&Math.abs(this.costs.reduce((n,c,i)=>n+c*solved.x[i],0)-this.linearLower!)>1e-7)throw Error('Native cone completion did not certify the original global objective.');
  this.coneRecoveryStarted=undefined;this.point=solved.x;if(this.initialSliceSeconds!==undefined){this.sliceSeconds=this.initialSliceSeconds;this.native.options.set({time_limit:this.sliceSeconds});}
  return{x:solved.x,value:(terms||this.program.cost.map((v,i)=>[i,v]as LinearTerm)).reduce((n,[i,v])=>n+v*solved.x[i],0)};
 }
 addRow(terms:LinearTerm[],lower:number,upper:number){
  this.linearLower=undefined;this.lowerObjectiveIsCurrent=false;
  if(this.conic){this.program.rows.push({terms,lower,upper});this.conic.invalidate();this.point=undefined;this.equalities=undefined;return;}
  if(this.decomposed){if(this.fullCompetitor)this.fullCompetitor.addRow(terms,lower,upper);else this.program.rows.push({terms,lower,upper});this.factorized?.dispose();this.factorized=undefined;this.factorizationFailed=false;this.point=undefined;this.equalities=undefined;return;}
  const k=this.program.rows.length;this.program.rows.push({terms,lower,upper});
  if(this.factorized){this.factorized.dispose();this.factorized=undefined;this.factorized=new FactorizedMarginalSolver(this.program,p=>new MarginalLinearSolver({...p,localBlocks:undefined},this.optimizer));this.point=undefined;this.equalities=undefined;return;}for(const[i,w]of terms)this.incidence.get(i)?.push([k,w]);this.equalities=undefined;
  if(this.point){const v=terms.reduce((n,[i,w])=>n+w*this.point![i],0);if(v<lower-1e-7||v>upper+1e-7)this.point=undefined;}
  this.native.addRows({numRows:1,lower:[lower],upper:[upper],matrix:{format:'csr',numRows:1,numCols:this.program.upper.length,starts:[0,terms.length],indices:terms.map(t=>t[0]),values:terms.map(t=>t[1])}});
 }
 /** Batch lazy hypothesis rows so the common decomposition is rebuilt once. */
 addRows(rows:MarginalProgram['rows']){
  if(!rows.length)return;this.lowerObjectiveIsCurrent=false;
  if(this.decomposed){if(this.fullCompetitor)this.fullCompetitor.addRows(rows);else this.program.rows.push(...rows);this.factorized?.dispose();this.factorized=undefined;this.factorizationFailed=false;this.point=undefined;this.equalities=undefined;return;}
  const offset=this.program.rows.length;this.program.rows.push(...rows);this.equalities=undefined;this.point=undefined;
  if(this.conic){this.conic.invalidate();return;}
  if(this.factorized){this.factorized.dispose();this.factorized=undefined;this.factorized=this.makeFactorized();return;}
  const starts=[0],indices:number[]=[],values:number[]=[];rows.forEach((row,k)=>{for(const[i,w]of row.terms){indices.push(i);values.push(w);this.incidence.get(i)?.push([offset+k,w]);}starts.push(indices.length);});this.native.addRows({numRows:rows.length,lower:rows.map(r=>r.lower),upper:rows.map(r=>r.upper),matrix:{format:'csr',numRows:rows.length,numCols:this.program.upper.length,starts,indices,values}});
 }
 bounds(terms:LinearTerm[]):[number,number]{if(!terms.length)return[0,0];const lo=this.solve(terms).value,hi=this.solve(terms,true).value;return[Math.max(0,lo-1e-7),Math.max(0,hi+1e-7)];}
 setRowBounds(indices:number[],lower:number[],upper:number[]){this.linearLower=undefined;this.lowerObjectiveIsCurrent=false;if(this.conic){indices.forEach((id,i)=>{this.program.rows[id].lower=lower[i];this.program.rows[id].upper=upper[i];});this.conic.invalidate();this.point=undefined;this.equalities=undefined;return;}if(this.decomposed){if(this.fullCompetitor)this.fullCompetitor.setRowBounds(indices,lower,upper);else indices.forEach((id,i)=>{this.program.rows[id].lower=lower[i];this.program.rows[id].upper=upper[i];});this.factorized?.dispose();this.factorized=undefined;this.factorizationFailed=false;this.point=undefined;this.equalities=undefined;return;}this.equalities=undefined;this.point=undefined;this.zeroFaceCertified=false;indices.forEach((id,i)=>{this.program.rows[id].lower=lower[i];this.program.rows[id].upper=upper[i];});if(this.factorized){this.factorized.dispose();this.factorized=undefined;this.factorized=new FactorizedMarginalSolver(this.program,p=>new MarginalLinearSolver({...p,localBlocks:undefined},this.optimizer));return;}this.native.changeRowsBounds({kind:'set',indices},lower,upper);}
 useDualForBoundChanges(){if(this.fullCompetitor)this.fullCompetitor.useDualForBoundChanges();else if(this.native)this.native.options.set({simplex_strategy:1});}
 feasibleMaximum(terms:LinearTerm[]){
  if(!terms.length&&!this.program.quadraticCost?.length&&this.optimizer.conic&&this.program.quadratics?.length&&this.program.quadratics.every(q=>q.upper>0&&Number.isFinite(1/Math.sqrt(q.upper)))){
   this.lowerObjectiveIsCurrent=false;this.linearLower=undefined;
   if(this.point&&validMarginalSolution(this.program,this.point)){this.linearLower=0;return{x:this.point.slice(),value:0};}
   const started=performance.now(),proof=this.spectralFeasibility;proof.attempts++;let solver:MarginalLinearSolver|undefined;
   try{
    // This is a proposal/dual-proof problem, not a replacement physical model.
    // A returned witness is checked against EVERY original cap, row and box.
    const p:MarginalProgram={...this.program,upper:this.program.upper.slice(),cost:this.program.cost.map(()=>0),rows:this.program.rows.slice(),localBlocks:undefined,quadratics:undefined,quadraticCost:this.program.quadratics.flatMap(q=>q.terms.map(f=>f.map(([i,w])=>[i,w/Math.sqrt(q.upper)]as LinearTerm)))};
    solver=new MarginalLinearSolver(p,this.optimizer);solver.setTimeSlice(this.initialSliceSeconds||600);const a=solver.solveToCompletion();proof.minimum=a.value;proof.originalViolation=marginalViolation(this.program,a.x);
    if(proof.originalViolation<=1e-7){proof.feasible++;this.point=a.x;this.linearLower=0;return{x:a.x.slice(),value:0};}
    proof.lower=quadraticTangentLowerBound(p,a.x,solver.rowDuals(),1e-7);const cap=this.program.quadratics.reduce((n,q)=>n+1+1e-7/q.upper,0);
    if(proof.lower>cap+1e-7){proof.infeasible++;return null;}
   }catch{/* An incomplete proposal cannot reject or certify the original model. */}
   finally{proof.seconds+=(performance.now()-started)/1000;solver?.dispose();}
   proof.fallbacks++;
  }
  try{return this.solve(terms,true);}catch(error){if(/did not complete: 8(?:;|$)/.test(String(error)))return null;throw error;}
 }
 objectiveLowerBound():number|undefined{if(this.linearLower!==undefined)return this.linearLower;if(this.winner==='flat'&&this.fullCompetitor)return this.fullCompetitor.objectiveLowerBound();return this.lowerObjectiveIsCurrent?this.conic?.objectiveLowerBound():this.linearLower;}
 supportPlanes():{terms:LinearTerm[];upper:number}[]{if(this.shortcut)return[];return this.conic?.supportPlanes()||(this.winner==='flat'?this.fullCompetitor?.supportPlanes():[])||[];}
 rowDuals():number[]{if(this.shortcut){const shadow=this.program.rows.map(()=>0);if(this.shortcut.row!==undefined)shadow[this.shortcut.row]=this.shortcut.multiplier!;return shadow;}if(this.conic)return this.conic.rowDuals();if(this.zeroFaceCertified)return this.program.rows.map(()=>0);return this.fullCompetitor?(this.winner==='flat'?this.fullCompetitor.rowDuals():[]):this.native?Array.from(this.native.getSolution().rowDual as ArrayLike<number>):[];}
 dispose(){
  this.conic?.dispose();this.factorized?.dispose();this.fullCompetitor?.dispose();this.native?.dispose();
  // Completed phases must not retain duplicated block graphs or witnesses
  // while the next phase constructs its full physical model.
  this.factorized=undefined;this.fullCompetitor=undefined;this.native=undefined;
  this.point=undefined;this.equalities=undefined;this.shortcut=undefined;this.costs=[];this.incidence.clear();this.inactive.clear();
 }
}
