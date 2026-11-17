import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import {incrementPrior} from './sanger-increment-prior';
import {MarginalLinearSolver,validMarginalSolution,type MarginalProgram,type LinearTerm} from './sanger-marginal-lp';
interface RepresentativeContext {baseFlowAnchor?:{forms:LinearTerm[][];mass:LinearTerm[];share:number}[];diagnostics?:(problem:{program:MarginalProgram;anchor:number[];photonForms:LinearTerm[][];originalColumns:number;proximityQuadraticCount:number})=>void;measurements?:{noise?:number;forms:LinearTerm[][];total:LinearTerm[];observed:number[]}[];substitutionOpportunities:number;noise:number;phaseSeed:{column:number;value:number}[];states:{size:number;column:number}[][][];shares:number[];edges:Parameters<typeof incrementPrior>[0]}
/** One representative convention for all inputs. Preserve independently
 * measured phase estimates while penalizing variant mass at the measured
 * noise scale. A soft increment prior resolves otherwise unconstrained
 * histories. All correlated/positional alternatives remain in the unchanged
 * parent confidence model used for bounds, GT retrieval and lazy queries. */
export function representativeBurden(program:MarginalProgram,edited:LinearTerm[],substitutionOpportunities=1){if(!(Number.isFinite(substitutionOpportunities)&&substitutionOpportunities>0))throw Error('A positive substitution opportunity count is required.');const cost=program.cost.map(c=>c/substitutionOpportunities);for(const[i,w]of edited)cost[i]+=w;if(cost.some(v=>!Number.isFinite(v)||v<0))throw Error('Representative burden must be finite and nonnegative.');return cost.flatMap((w,i)=>w?[[i,w]as LinearTerm]:[]);}
export function selectMarginalRepresentative(program:MarginalProgram,point:number[],edited:LinearTerm[],optimizer:SangerLinearOptimizer,context:RepresentativeContext){
 const started=performance.now();if(!(Number.isFinite(context.noise)&&context.noise>0))throw Error('Representative selection requires positive measured noise.');if(!validMarginalSolution(program,point))throw Error('Representative selection requires a valid original fit.');
 const terms=representativeBurden(program,edited,context.substitutionOpportunities),cost=program.cost.map(()=>0);for(const[i,w]of terms)cost[i]=w;const p:MarginalProgram={...program,objectiveUnitScale:1e-4,upper:program.upper.slice(),cost,rows:program.rows.slice(),localBlocks:undefined},constant=p.upper.length,scale=1/Math.sqrt(2*context.noise);p.upper.push(1);p.cost.push(0);p.rows.push({terms:[[constant,1]],lower:1,upper:1});const phaseShares=new Map(context.states.flatMap(hs=>hs.flatMap((states,h)=>states.map(s=>[s.column,context.shares[h]]as const))));
 // Linear variant mass scales with reference share. Normalize squared phase
 // deviations by the same share so dilution cannot amplify shrinkage.
 p.quadraticCost=context.phaseSeed.map(s=>{const share=phaseShares.get(s.column);if(share===undefined||!Number.isFinite(share)||share<0)throw Error('Missing valid reference share for a phase state.');const weight=share>0?scale/Math.sqrt(share):0;return[[s.column,weight],[constant,-s.value*weight]];});
 const prior=incrementPrior(context.edges,context.states,context.phaseSeed,context.shares,p.upper.length,context.noise);p.upper.push(...prior.upper);p.cost.push(...prior.upper.map(()=>0));p.rows.push(...prior.rows);p.quadraticCost.push(...prior.forms.map(form=>form.map(([i,w])=>[i,w*scale]as LinearTerm)));
 // Rank feasible allocations by the observed local fluorescence, rather than
 // only selecting an interior point of its noise envelope. This changes the
 // representative objective, never the original measured constraints or the
 // alternative/GT domain. A fixed positive photon denominator from the valid
 // anchor makes the weighted least-squares selection convex; no GT or allele
 // dictionary participates. Mean squared residual per nominal target window
 // avoids changing the penalty merely by adding more measurement rows. The
 // noise scale is measured in the WT prefix, not selected against GT.
 const photonForms=(context.measurements||[]).flatMap(m=>{
  const denominator=m.total.reduce((n,[i,w])=>n+w*point[i],0);
  if(!(denominator>0))throw Error('Representative photon fit needs positive measured intensity.');
  return m.forms.map((form,b)=>{const terms=new Map<number,number>();
   for(const[i,w]of form)terms.set(i,(terms.get(i)||0)+w);
   for(const[i,w]of m.total)terms.set(i,(terms.get(i)||0)-w*m.observed[b]);
   return[...terms].filter(([,w])=>w!==0).map(([i,w])=>[i,w/(denominator*Math.sqrt(8*context.substitutionOpportunities)*(m.noise??context.noise))]as LinearTerm);
  });
 });
 // Keep the least-squares Hessian sparse: signed residuals are represented
 // by two nonnegative auxiliaries, with bounds derived from ORIGINAL boxes.
 // Their difference is exactly the measured residual. This is an objective
 // reformulation, not an added signal constraint or a hypothesis restriction.
 const photonAux=photonForms.map(form=>{
  const lower=form.reduce((n,[i,w])=>n+Math.min(0,w*program.upper[i]),0),upper=form.reduce((n,[i,w])=>n+Math.max(0,w*program.upper[i]),0),positive=p.upper.length,negative=positive+1;
  p.upper.push(Math.max(0,upper),Math.max(0,-lower));p.cost.push(0,0);
  p.rows.push({terms:[[positive,1],[negative,-1],...form.map(([i,w])=>[i,-w]as LinearTerm)],lower:0,upper:0});
  p.quadraticCost!.push([[positive,1],[negative,-1]]);
  return{form,positive,negative};
 });
 // Both photon evidence and this public-base stabilization are averaged
 // over genomic opportunities. Summing only stabilization made arbitrary
 // anchor allocations stronger merely because the target window was longer.
 // A noise-scale proximal step preserves the existing allocation along weak
 // inverse directions instead of following ADC-scale differences to a distant
 // molecular explanation. Penalize PUBLIC base masses, not internal columns
 // or an enumerated genotype count. Reference dilution has the same linear
 // mass scaling as the phase prior. No original hypothesis is removed.
 const proximityForms=(context.baseFlowAnchor||[]).flatMap(g=>{
  const mass=g.mass.reduce((n,[i,w])=>n+w*point[i],0);if(mass<=1e-7)return[];
  if(!(g.share>0))throw Error('Base-flow stabilization needs a positive reference share.');
  const scale=1/Math.sqrt(2*context.noise*g.share*context.substitutionOpportunities);
  return g.forms.map(form=>{const fraction=Math.max(0,form.reduce((n,[i,w])=>n+w*point[i],0))/mass,terms=new Map<number,number>();for(const[i,w]of form)terms.set(i,(terms.get(i)||0)+w);for(const[i,w]of g.mass)terms.set(i,(terms.get(i)||0)-w*fraction);return[...terms].filter(([,w])=>w!==0).map(([i,w])=>[i,w*scale]as LinearTerm);});
 });
 p.quadraticCost.push(...proximityForms);
 const dot=(form:LinearTerm[],x:number[])=>form.reduce((s,[i,w])=>s+w*x[i],0),burden=(x:number[])=>dot(terms,x),objective=(x:number[])=>burden(x)+p.quadraticCost!.reduce((s,f)=>s+dot(f,x)**2,0),anchor=point.slice();anchor.push(1);for(const row of prior.rows)anchor.push(-row.terms.slice(1).reduce((s,[i,w])=>s+w*anchor[i],0));for(const a of photonAux){const v=dot(a.form,point);anchor[a.positive]=Math.max(0,v);anchor[a.negative]=Math.max(0,-v);}const before=objective(anchor),beforeBurden=burden(point);
 context.diagnostics?.({program:p,anchor,photonForms,originalColumns:program.upper.length,proximityQuadraticCount:proximityForms.length});
 const finish=(x:number[],minimum:number,lower:number,certificate:string,solver?:MarginalLinearSolver['diagnostics'])=>({point:x.slice(0,program.upper.length),summary:{method:'noise-regularized-sequence-flow' as const,before,minimum,lower,beforeBurden,afterBurden:burden(x),noise:context.noise,photonRows:photonForms.length,beforePhotonLoss:photonForms.reduce((s,f)=>s+dot(f,point)**2,0),afterPhotonLoss:photonForms.reduce((s,f)=>s+dot(f,x)**2,0),baseFlowStabilization:proximityForms.length,afterProximityLoss:proximityForms.reduce((s,f)=>s+dot(f,x)**2,0),substitutionOpportunities:context.substitutionOpportunities,phaseSquaredDistance:context.phaseSeed.reduce((s,v)=>s+(x[v.column]-v.value)**2,0),incrementSquaredDeviation:prior.forms.reduce((s,f)=>s+dot(f,x)**2,0),seconds:(performance.now()-started)/1000,certificate,solver}});
 if(before<=1e-7)return finish(anchor,before,0,'nonnegative-box-bound');
 const solver=new MarginalLinearSolver(p,optimizer);solver.setTimeSlice(600);try{const answer=solver.solveToCompletion();
 // Rebuild private residual identities from physical columns before checking
 // the ORIGINAL objective and global lower bound. Never certify a rounded
 // auxiliary objective as the molecular objective. No physical column changes.
 for(const row of prior.rows)answer.x[row.terms[0][0]]=-row.terms.slice(1).reduce((sum,[i,w])=>sum+w*answer.x[i],0);for(const a of photonAux){const residual=dot(a.form,answer.x);answer.x[a.positive]=Math.max(0,residual);answer.x[a.negative]=Math.max(0,-residual);}const minimum=objective(answer.x),lower=solver.objectiveLowerBound(),x=answer.x.slice(0,program.upper.length);if(!validMarginalSolution(program,x)||!validMarginalSolution(p,answer.x)||lower===undefined||Math.abs(minimum-lower)>1e-7)throw Error('Representative selection did not certify its original model and global objective.');return finish(answer.x,minimum,lower,'original-primal-and-global-dual',solver.diagnostics);}finally{solver.dispose();}
}
