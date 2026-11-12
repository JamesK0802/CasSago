import {emissionUncertainty} from './sanger-emission-uncertainty';
import {calibrateRegionalResponse} from './sanger-regional-calibration';
import {rebasePhotonForms} from './sanger-photon-state-rebase';
import {verifiedZeroDeparturePoint} from './sanger-initial-marginal-fit';
import {selectMarginalRepresentative} from './sanger-representative-selection';
import {spectralBoxCertificates} from './sanger-spectral-box-certificate';
import {seedQuadraticTangents} from './sanger-quadratic-envelope';
import {regionalLengthInitialization,type RegionalLengthWindow} from './sanger-regional-length-fit';
import {originPrefixes} from './sanger-origin-prefix';
import {localBoundCertificates} from './sanger-local-bound-certificate';
import {withinLengthDisplay} from './sanger-marginal-display';
import {fittedMarginalRepresentative} from './sanger-marginal-representative';
import {marginalEndpointWitnesses} from './sanger-marginal-witness';
import {measuredControlKernels} from './sanger-control-kernels';
import {simultaneousResidualLimit} from './sanger-marginal-noise';
import {boundMarginalSubstitutions} from './sanger-marginal-substitution';
import type {SangerAnalysisResult,SangerMarginalGroup,SangerJointEvent} from '../../models/sanger-analysis';
import type {SignalPathInput,Prepared} from './sanger-path-common';
import {prepare} from './sanger-path-common';
import {prepareMeasuredObservation} from './sanger-measured-mixture';
import {SangerFitCache} from './sanger-fit-cache';
import {mixPhase,phaseNeighbor} from './sanger-signal';
import {plotTrace} from './sanger-analysis';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import {MarginalLinearSolver,LP_INFINITY,type MarginalProgram,type LinearTerm} from './sanger-marginal-lp';
const DNA='ACGT',SUB_RADIUS=45;
export interface MarginalEdge {prefixColumns?:number[];reference:number;target:number;incoming:number;outgoing:number;size:number;start:number;mass:number;free:Map<number,number[]>}
export interface MarginalInternalGroup {sharedInsertions?:boolean;flow?:Map<number,{unchanged:{base:number;column:number;cap:LinearTerm[]}[];changed:number[]}>;dye?:Map<number,LinearTerm[][]>;departures?:Map<number,LinearTerm[]>;id:string;reference:number;target:number;size:number;mass:LinearTerm[];bases:Map<number,LinearTerm[][]>;edges:MarginalEdge[];fraction:number;fractionRange:[number,number]}
export interface MarginalMeasurement {position:number;forms:LinearTerm[][];total:LinearTerm[];observed:number[]}
export interface MarginalTargetState {primaryPoint?:number[];measurements?:MarginalMeasurement[];program:MarginalProgram;groups:MarginalInternalGroup[];edges:MarginalEdge[];point:number[];positions:number[];tolerance:number}
export interface MarginalSession {refs:Prepared[];shares:number[];targets:MarginalTargetState[];optimizer:SangerLinearOptimizer;pages:Map<string,unknown>}
export const marginalSessions=new WeakMap<SangerAnalysisResult,MarginalSession>();
/** Optional scientific diagnostics; no truth or inference policy is supplied. */
export const marginalModelDiagnostics:{observe?:(model:{program:MarginalProgram;refs:Prepared[];shares:number[];targets:MarginalTargetState[];states:{size:number;column:number}[][][];noise:number;phaseSeed:{column:number;value:number}[];measurements:MarginalMeasurement[];spectralPositions:number[];departures:LinearTerm[]})=>void;primaryFit?:(proof:{point:number[];value:number;lower:number|undefined;shadow:number[];planes:{terms:LinearTerm[];upper:number}[]})=>void}={};
const evaluate=(terms:LinearTerm[],x:number[])=>terms.reduce((n,[i,w])=>n+w*x[i],0);
const sumTerms=(lists:LinearTerm[][])=>{const m=new Map<number,number>();for(const list of lists)for(const[i,w]of list)m.set(i,(m.get(i)||0)+w);return[...m].filter(([,w])=>Math.abs(w)>1e-14)as LinearTerm[];};
const identical=(a:LinearTerm[],b:LinearTerm[])=>{const x=new Map(a),y=new Map(b);return x.size===y.size&&[...x].every(([i,v])=>y.get(i)===v);};
export function marginalInputEvent(ref:Prepared,start:number,end:number,inserted:string,target:number):SangerJointEvent{const a=ref.offset+start,b=ref.offset+end;return{target:ref.reverse?ref.cuts.length-target:target+1,start:ref.reverse?ref.input.sequence.length-b:a,end:ref.reverse?ref.input.sequence.length-a:b,inserted:inserted.length,insertedSequence:ref.reverse?[...inserted].reverse().map(c=>({A:'T',C:'G',G:'C',T:'A'}[c]||'N')).join(''):inserted};}
/** The SAME continuous state-flow model for every sample, reference and target.
 * Flanks identify length marginals; local geometry and base flows retain
 * uncertainty. No input allele count, truth, reporting threshold or synth selector
 * enters fitting. A single fitted display proposal is saved per family; further alternatives are lazy. */
export function analyzeMarginalSanger(input:SignalPathInput,progress:(s:string)=>void,optimizer:SangerLinearOptimizer):SangerAnalysisResult{
 const representativeRefit:NonNullable<NonNullable<SangerAnalysisResult['diagnostics']>['representativeRefit']>=[],optimization:NonNullable<NonNullable<SangerAnalysisResult['diagnostics']>['optimization']>=[],timings:{stage:string;seconds:number}[]=[],started=performance.now();let last=started,stage='Registering measured fluorescence';
 const tick=(next:string)=>{const now=performance.now();if(now>last)timings.push({stage,seconds:(now-last)/1000});last=now;stage=next;progress(next);};
 const {refs,report,maxDeletion,maxInsertion}=prepare(input.amplicons,input.primer,input.settings),observation=prepareMeasuredObservation(refs,input.control,input.sample,input.settings,s=>tick(s)),{raw,means,calibration,quality}=observation,shares=observation.shares.slice(),noise=Math.sqrt(observation.noise),analysisEnd=Math.min(observation.end,observation.c.length-maxDeletion*refs[0].cuts.length-1,quality.end-maxDeletion*refs[0].cuts.length-1),targets:MarginalTargetState[]=[],groups:SangerMarginalGroup[]=[],lengths:{size:number;fraction:number}[][][]=[],count=refs[0].cuts.length;let squaredError=0,signalVariance=0,signalRows=0,heldOutError=0,heldOutVariance=0,heldOutRows=0;
 let regionalResponse=(h:number,p:number)=>1;
 const responseScale=shares.map(()=>1),emissions=observation.intensity.responses,emissionMeans=observation.intensity.novel,dyeAt=(h:number,p:number,dye:number[])=>observation.transferByReference.apply(h,p,dye).map(v=>v*responseScale[h]*regionalResponse(h,p)),emit=(h:number,shift:number,p:number)=>dyeAt(h,p,emissions[h][p-shift]||emissionMeans[DNA.indexOf(refs[h].sequence[p-shift])]||[0,0,0,0]),normalized=(v:number[])=>{const n=v.reduce((a,b)=>a+b,0);return v.map(x=>n?x/n:0);};
 tick('Fitting invariant-flank length proportions');
 const flankPositions:number[][]=[],flowWindows:RegionalLengthWindow[]=[];
 for(let t=0;t<count;t++){
  const left=Math.max(...refs.map(r=>r.cuts[t]))+Math.max(SUB_RADIUS+t*maxInsertion,(t+1)*maxInsertion+3)+5,right=Math.min(analysisEnd,...refs.map(r=>t+1<count?r.cuts[t+1]-Math.max(SUB_RADIUS+(t+1)*maxDeletion,(t+2)*maxDeletion+3)-5:analysisEnd)),ps=Array.from({length:Math.max(0,right-left)},(_,i)=>left+i);
  if(ps.length<Math.max(20,Math.ceil(refs.length*((t+1)*(maxDeletion+maxInsertion)+1)/2)))throw Error('Insufficient invariant flank to distinguish the requested reference/length classes.');
  flankPositions.push(ps);flowWindows.push({classes:refs.flatMap((_,h)=>Array.from({length:(t+1)*(maxDeletion+maxInsertion)+1},(_,i)=>({reference:h,size:i-(t+1)*maxDeletion}))).map(c=>({...c,column:ps.flatMap(p=>mixPhase(emit(c.reference,c.size,p),emit(c.reference,c.size,phaseNeighbor(p,calibration.lag)),calibration.lag))})),observed:ps.flatMap(p=>raw[p])});
 }
 const preliminary=regionalLengthInitialization(flowWindows,flankPositions,shares,observation.fitDelta,observation.transferByReference),initialization=calibrateRegionalResponse(flowWindows,flankPositions,shares,preliminary.delta,{positions:observation.trainingPositions,columns:refs.map((_,h)=>observation.trainingPositions.flatMap(p=>mixPhase(emit(h,0,p),emit(h,0,phaseNeighbor(p,calibration.lag)),calibration.lag))),observed:observation.trainingPositions.flatMap(p=>raw[p])}),robustDelta=initialization.delta;
 // Reference brightness is a nuisance distinct from control molecule share.
 // Its flank estimate uses a paired-prefix predictive prior, whose variance
 // includes the actual ADC floor. The same correction applies to every input.
 regionalResponse=initialization.at;
 lengths.push(...initialization.lengths);
 for(let t=0;t<count;t++)for(let i=0;i<flankPositions[t].length;i++)if(Math.floor(flankPositions[t][i]/5)%3===2)for(let b=0;b<4;b++){const y=raw[flankPositions[t][i]][b],prediction=initialization.trainingPredictions[t][4*i+b];heldOutError+=(y-prediction)**2;heldOutVariance+=(y-.25)**2;heldOutRows++;}
 const program:MarginalProgram={upper:[],cost:[],rows:[]},variable=(upper:number,cost=0)=>{const id=program.upper.length;program.upper.push(upper);program.cost.push(cost);return id;},row=(terms:LinearTerm[],lower:number,upper:number)=>program.rows.push({terms:sumTerms([terms]),lower,upper}),states=flowWindows.map(w=>refs.map((_,h)=>w.classes.flatMap(c=>c.reference===h?[{size:c.size,fraction:shares[h],column:variable(shares[h])}]:[]))),drafts:(MarginalTargetState&{qs:number[];positionForms:Map<number,unknown[]>})[]=[],photons=new Map<number,LinearTerm[][]>(),intensityFloor=new Map<number,number>(),allDepartures:LinearTerm[]=[],allEdges:MarginalEdge[]=[],blockColumns=new Map<number,Set<number>>();
 const uncertainEmission=emissionUncertainty(program,shares,refs.map(r=>r.sequence),emissions,observation.intensity.allocationBounds,blockColumns);
 for(let t=0;t<count;t++)for(let h=0;h<refs.length;h++)row(states[t][h].map(c=>[c.column,1]),shares[h],shares[h]);
 for(let target=0;target<count;target++){
  tick('Constructing joint marginal constraints · target '+(target+1));
  const incoming=refs.map((_,h)=>target?states[target-1][h]:[{size:0,fraction:shares[h],column:undefined as number|undefined}]),outgoing=states[target],lo=Math.max(quality.start,Math.min(...refs.map((r,h)=>r.cuts[target]-Math.max(SUB_RADIUS,maxDeletion+3)+Math.min(...incoming[h].filter(c=>c.fraction>0).map(c=>c.size))))-1),hi=Math.min(analysisEnd,Math.max(...refs.map((r,h)=>r.cuts[target]+Math.max(SUB_RADIUS,maxInsertion+3)+Math.max(...incoming[h].filter(c=>c.fraction>0).map(c=>c.size))))+2),qs=Array.from({length:hi-lo},(_,i)=>lo+i),positions=qs.slice(1,-1),tolerance=simultaneousResidualLimit(noise,4*positions.length*count),edges:MarginalEdge[]=[],measurements:MarginalMeasurement[]=[],minimumIntensity=new Map<number,number>(),internal=new Map<string,MarginalInternalGroup>();
  const nominalPhotons=new Map<string,Map<number,LinearTerm[][]>>();
  const id=(h:number,size:number)=>'marginal:'+h+':'+target+':'+size;
  for(let h=0;h<refs.length;h++)for(const a of incoming[h])if(a.fraction>0)for(const b of outgoing[h])if(b.fraction>0){
   const size=b.size-a.size;if(size< -maxDeletion||size>maxInsertion)continue;const cut=refs[h].cuts[target],startMin=size<0?cut+size-3:size>0?cut-3:cut,startMax=size?cut+3:cut,seen=new Set<string>();
   for(let start=startMin;start<=startMax;start++){
    if(start<0||start+Math.max(0,-size)>refs[h].sequence.length)continue;
    if(size<0){const sequence=refs[h].sequence.slice(0,start)+refs[h].sequence.slice(start-size);if(seen.has(sequence))continue;seen.add(sequence);}
    const mass=variable(Math.min(a.fraction,b.fraction)),edge:MarginalEdge={reference:h,target,incoming:a.size,outgoing:b.size,size,start,mass,free:new Map()};edges.push(edge);const key=id(h,size);let g=internal.get(key);if(!g){g={id:key,reference:h,target,size,mass:[],bases:new Map(),departures:new Map(),dye:new Map(),flow:new Map(),edges:[],fraction:0,fractionRange:[0,0]};internal.set(key,g);}g.mass.push([mass,1]);g.edges.push(edge);

   }
  }
  const geometry=new Map([...internal.values()].map(g=>[g,originPrefixes(g,program,variable,row)])),insertionFlows=new Map<MarginalInternalGroup,Map<number,number[]>>();
  // Exact shared insertion colour transport, expressed through sparse origin
  // prefixes instead of repeating every covering origin in every cycle row.
  for(const g of internal.values()){
   const covering=new Map<number,MarginalEdge[]>(),flows=new Map<number,number[]>();insertionFlows.set(g,flows);
   for(const e of g.edges)for(let p=e.start+e.incoming;p<e.start+e.incoming+Math.max(0,e.size);p++){const es=covering.get(p)||[];es.push(e);covering.set(p,es);}
   for(const[p,es]of covering){const flow=Array.from({length:4},()=>variable(shares[g.reference])),capacity=sumTerms(geometry.get(g)!.flatMap(b=>[b.at(p-b.incoming),b.at(p-b.incoming-Math.max(0,g.size)).map(([i,w])=>[i,-w]as LinearTerm)]));row([...flow.map(i=>[i,1]as LinearTerm),...capacity.map(([i,w])=>[i,-w]as LinearTerm)],0,0);flows.set(p,flow);for(const e of es)e.free.set(p,flow);if(es.length>1)g.sharedInsertions=true;}
  }
  for(let h=0;h<refs.length;h++){
   for(const a of incoming[h])if(a.fraction>0){const terms=[...internal.values()].filter(g=>g.reference===h).flatMap(g=>geometry.get(g)!.filter(b=>b.incoming===a.size).map(b=>[b.total,1]as LinearTerm));if(a.column===undefined)row(terms,a.fraction,a.fraction);else row([...terms,[a.column,-1]],0,0);}
   for(const b of outgoing[h])if(b.fraction>0)row([... [...internal.values()].filter(g=>g.reference===h).flatMap(g=>geometry.get(g)!.filter(a=>a.outgoing===b.size).map(a=>[a.total,1]as LinearTerm)),[b.column,-1]],0,0);
  }
  const departures:LinearTerm[]=[],positionForms=new Map<number,unknown[]>();
  for(const g of internal.values())for(const p of qs){
   // Group identical nucleotide/dye contexts; this is a lossless column merge
   // for every input, not a synthetic-data or genotype-count selector.
   const map=new Map<string,{key:string;base:number;dye:number[];cap:LinearTerm[];donor?:number}>();
   const add=(base:number,dye:number[],cap:LinearTerm[],donor?:number)=>{cap=sumTerms([cap]);if(!cap.length)return;const key=base+':'+dye.join(',')+(donor!==undefined&&observation.intensity.allocationBounds[g.reference]?.[donor]?.some(v=>Math.abs(v)>1e-12)?':'+donor:''),old=map.get(key);if(old)old.cap=sumTerms([old.cap,cap]);else map.set(key,{key,base,dye,cap,donor});};
   for(const bucket of geometry.get(g)!){const before=sumTerms([[[bucket.total,1]],bucket.at(p-bucket.incoming).map(([i,w])=>[i,-w]as LinearTerm)]),after=bucket.at(p-bucket.incoming-Math.max(0,g.size));for(const[at,cap]of [[p-bucket.incoming,before],[p-bucket.outgoing,after]]as const)if(cap.length){const b=DNA.indexOf(refs[g.reference].sequence[at]);if(b<0)throw Error('Resolve ambiguous reference bases in the modeled target window.');add(b,emit(g.reference,p-at,p),cap,at);}}
   const free=insertionFlows.get(g)!.get(p);if(free)for(let b=0;b<4;b++)add(b,dyeAt(g.reference,p,emissionMeans[b]),[[free[b],1]]);
   const contexts=[...map.values()].sort((a,b)=>a.key.localeCompare(b.key));minimumIntensity.set(p,Math.min(minimumIntensity.get(p)??Infinity,...contexts.map(c=>c.dye.reduce((n,v)=>n+v,0)*(c.donor===undefined?1:1+(observation.intensity.allocationBounds[g.reference]?.[c.donor]?.[0]||0))),...emissionMeans.map(row=>dyeAt(g.reference,p,row).reduce((n,v)=>n+v,0))));contexts.forEach(c=>c.cap=sumTerms([c.cap]));const expected=Array.from({length:4},(_,b)=>sumTerms(contexts.filter(c=>c.base===b).map(c=>c.cap))),subAllowed=g.edges.some(e=>Math.abs(p-e.incoming-refs[g.reference].cuts[target])<=SUB_RADIUS),signatures=positionForms.get(p)||[];signatures.push([subAllowed,contexts.map(c=>[c.base,c.dye,c.cap])]);positionForms.set(p,signatures);
   const nominal=Array.from({length:4},(_,b)=>sumTerms(contexts.map(c=>c.cap.map(([i,w])=>[i,w*c.dye[b]]as LinearTerm)))),saved=nominalPhotons.get(g.id)||new Map<number,LinearTerm[][]>();saved.set(p,nominal);nominalPhotons.set(g.id,saved);
   const correction=(cap:LinearTerm[],c:typeof contexts[number])=>c.donor===undefined?[]:uncertainEmission(cap,g.reference,c.donor,p);
   if(!subAllowed){g.bases.set(p,expected);g.dye!.set(p,Array.from({length:4},(_,b)=>sumTerms([nominal[b],...contexts.map(c=>correction(c.cap,c).map(([i,w])=>[i,w*c.dye[b]]as LinearTerm))])));continue;}
   const pure=contexts.length===1&&identical(contexts[0].cap,g.mass),unchanged=contexts.map(()=>variable(shares[g.reference])),changed=Array.from({length:4},(_,b)=>variable(pure&&b===contexts[0].base?0:shares[g.reference],1));if(!pure)contexts.forEach((c,j)=>row([[unchanged[j],1],...c.cap.map(([i,w])=>[i,-w]as LinearTerm)],-LP_INFINITY,0));
   row([...unchanged.map(i=>[i,1]as LinearTerm),...changed.map(i=>[i,1]as LinearTerm),...g.mass.map(([i,w])=>[i,-w]as LinearTerm)],0,0);
   // A substitution cannot keep the same nucleotide. These four Hall bounds
   // are equivalent to existence of a nonnegative off-diagonal 4x4 transport,
   // eliminating twelve transition variables without losing any base flow.
   // The mass balance eliminates the long common-mass expression from
   // Hall's inequalities. A single immutable source context needs only the
   // mass balance and a zero bound on its forbidden diagonal substitution.
   if(!pure)for(let b=0;b<4;b++)row([...expected[b],...contexts.flatMap((c,j)=>c.base===b?[[unchanged[j],-1]as LinearTerm]:[]),...changed.flatMap((i,a)=>a!==b?[[i,-1]as LinearTerm]:[])],-LP_INFINITY,0);
   g.flow!.set(p,{unchanged:contexts.map((c,j)=>({base:c.base,column:unchanged[j],cap:c.cap})),changed});
   const norm=changed.map(i=>[i,1]as LinearTerm);departures.push(...norm);g.departures!.set(p,norm);g.bases.set(p,Array.from({length:4},(_,b)=>[[changed[b],1]as LinearTerm,...contexts.flatMap((c,j)=>c.base===b?[[unchanged[j],1]as LinearTerm]:[])]));g.dye!.set(p,Array.from({length:4},(_,b)=>sumTerms([contexts.map((c,j)=>[unchanged[j],c.dye[b]]as LinearTerm),...contexts.map((c,j)=>correction([[unchanged[j],1]],c).map(([i,w])=>[i,w*c.dye[b]]as LinearTerm)),changed.map((i,a)=>[i,dyeAt(g.reference,p,emissionMeans[a])[b]]as LinearTerm)])));
  }
  const junctions=refs.map((_,h)=>{const own=edges.filter(e=>e.reference===h);return{first:Math.min(...own.map(e=>e.start+e.incoming)),last:Math.max(...own.map(e=>e.start+e.incoming+Math.max(0,e.size)))};});
  for(const q of qs){const byReference=refs.map((_,h)=>{const own=[...internal.values()].filter(g=>g.reference===h),original=[0,1,2,3].map(b=>sumTerms(own.map(g=>g.dye!.get(q)?.[b]||[]))),before=q<junctions[h].first,after=q>=junctions[h].last;if(!before&&!after)return original;
   const baseline:LinearTerm[][]=[[],[],[],[]];for(const g of own)for(let b=0;b<4;b++)baseline[b].push(...nominalPhotons.get(g.id)!.get(q)![b]);
   // Before the first cut these outgoing columns sum to the fixed reference
   // share; every column uses the SAME unshifted emission, not its future phase.
   const referenceStates=before&&target?incoming[h]:outgoing[h],replacement=[0,1,2,3].map(b=>referenceStates.map(state=>[state.column!,emit(h,before?(target?state.size:0):state.size,q)[b]]as LinearTerm));return rebasePhotonForms(original,baseline,replacement);
  }),forms=[0,1,2,3].map(b=>sumTerms(byReference.map(forms=>forms[b])));if(photons.has(q))throw Error('The requested editing neighborhoods overlap; independent target windows are not supported by these measurements.');photons.set(q,forms);intensityFloor.set(q,minimumIntensity.get(q)!);}
  allDepartures.push(...departures);allEdges.push(...edges);
  for(const g of internal.values()){for(const[q,flow]of g.flow||[]){const columns=blockColumns.get(q)||new Set<number>();flow.unchanged.forEach(c=>columns.add(c.column));flow.changed.forEach(i=>columns.add(i));blockColumns.set(q,columns);}for(const e of g.edges)for(const[q,free]of e.free){const columns=blockColumns.get(q)||new Set<number>();free.forEach(i=>columns.add(i));blockColumns.set(q,columns);}}
  drafts.push({program,groups:[...internal.values()],edges,point:[],positions,tolerance,measurements,qs,positionForms});
 }
 // A target's bounding interval can contain unedited gaps between references.
 // Those cycles are invariant-flank evidence, not cut-proximal SNP measurements.
 // Retain their exact origin-flow photons and the same flank RSS confidence rule.
 const firstMeasured=Math.min(...drafts.map(d=>d.positions[0])),allMeasurements:MarginalMeasurement[]=[],localMeasurements:MarginalMeasurement[]=[],localReadCycles=new Set(drafts.flatMap((d,t)=>d.positions.filter(p=>refs.some((r,h)=>{const incoming=t?states[t-1][h].map(s=>s.size):[0];return p>=r.cuts[t]-Math.max(SUB_RADIUS,maxDeletion+3)+Math.min(...incoming)&&p<=r.cuts[t]+Math.max(SUB_RADIUS,maxInsertion+3)+Math.max(...incoming);})))),spectralMeasurements=new Map<number,MarginalMeasurement[]>();
 const invariant=(q:number)=>{let t=-1;for(let j=0;j<drafts.length;j++)if(q>drafts[j].positions.at(-1)!)t=j;if(t<0)throw Error('No calibrated cumulative phase is available for the requested measurement.');const forms:LinearTerm[][]=[[],[],[],[]];let floor=Infinity;for(let h=0;h<refs.length;h++)for(const state of states[t][h]){const dye=emit(h,state.size,q);floor=Math.min(floor,dye.reduce((a,b)=>a+b,0)*(1+(observation.intensity.allocationBounds[h]?.[q-state.size]?.[0]||0)));const correction=uncertainEmission([[state.column,1]],h,q-state.size,q);for(let b=0;b<4;b++)forms[b].push([state.column,dye[b]],...correction.map(([i,w])=>[i,w*dye[b]]as LinearTerm));}photons.set(q,forms);intensityFloor.set(q,floor);return forms;};
 for(let p=firstMeasured;p<analysisEnd;p++){const neighbor=phaseNeighbor(p,calibration.lag),spectral=!localReadCycles.has(p)&&(!calibration.lag||!localReadCycles.has(neighbor)),forms=[0,1,2,3].map(b=>sumTerms([[p,1-Math.abs(calibration.lag)],[neighbor,Math.abs(calibration.lag)]].filter(([,w])=>w>0).map(([q,w])=>(photons.get(q)||invariant(q))[b].map(([i,v])=>[i,v*w]as LinearTerm)))),total=sumTerms(forms),floor=Math.min(intensityFloor.get(p)!,intensityFloor.get(neighbor)??intensityFloor.get(p)!);if(!(floor>0))throw Error('Missing positive measured donor fluorescence.');row(total,floor,LP_INFINITY);const m={position:p,forms,total,observed:raw[p]};allMeasurements.push(m);if(spectral){let t=0;for(let j=0;j<drafts.length;j++)if(p>drafts[j].positions.at(-1)!)t=j;const list=spectralMeasurements.get(t)||[];list.push(m);spectralMeasurements.set(t,list);}else localMeasurements.push(m);}
 const tolerance=simultaneousResidualLimit(noise,4*localMeasurements.length);for(const m of localMeasurements)for(let b=0;b<4;b++){row([...m.forms[b],...m.total.map(([i,w])=>[i,-w*(m.observed[b]-tolerance)]as LinearTerm)],0,LP_INFINITY);row([...m.forms[b],...m.total.map(([i,w])=>[i,-w*(m.observed[b]+tolerance)]as LinearTerm)],-LP_INFINITY,0);}
 // A regression confidence region uses all raw flank residuals. Unlike an
 // unjustified per-peak maximum, it permits the isolated errors present in
 // measured reads without deleting them or changing the held-out QC rule.
 const seed=program.upper.map(()=>0);for(let t=0;t<count;t++)for(let h=0;h<refs.length;h++)for(const state of states[t][h])seed[state.column]=lengths[t][h].find(c=>c.size===state.size)!.fraction;
 uncertainEmission.seed(seed);
 program.quadratics=[...spectralMeasurements].map(([t,list])=>{const terms=list.flatMap(m=>{const denominator=Math.max(1e-30,evaluate(m.total,seed));return m.forms.map((form,b)=>sumTerms([form,m.total.map(([i,w])=>[i,-w*m.observed[b]]as LinearTerm)]).map(([i,w])=>[i,w/denominator]as LinearTerm));}),rss=terms.reduce((n,r)=>n+evaluate(r,seed)**2,0),parameters=Math.max(1,states[t].reduce((n,g)=>n+g.length-1,0)),variance=Math.max(observation.noise,rss/Math.max(1,terms.length-parameters)),tail=Math.log(Math.max(1,spectralMeasurements.size)/.0027),allowance=variance*(parameters+2*Math.sqrt(parameters*tail)+2*tail);return{name:'invariant-spectrum-'+(t+1),terms,upper:rss+allowance};});
 for(const q of program.quadratics){const t=Number(q.name.split('-').at(-1))-1;for(const cut of seedQuadraticTangents(q,seed,states[t].map((g,h)=>({columns:g.map(s=>s.column),mass:shares[h]}))))row(cut.terms,-LP_INFINITY,cut.upper);}
 // These are certified consequences of the same spectral ellipsoids and
 // reference-mass equations; exposing them tightens numerical conditioning
 // without narrowing the original feasible set.
 for(const q of program.quadratics){const t=Number(q.name.split('-').at(-1))-1,known=new Set(states[t].flatMap(g=>g.map(s=>s.column)));if(q.terms.some(row=>row.some(([i])=>!known.has(i))))continue;const certificate=spectralBoxCertificates(q,program.upper,states[t].map((g,h)=>({columns:g.map(s=>s.column),mass:shares[h]})));for(const b of certificate.bounds){program.upper[b.column]=Math.max(0,b.upper);if(b.lower>0)row([[b.column,1]],b.lower,LP_INFINITY);}}
 for(const e of allEdges){const out=states[e.target][e.reference].find(s=>s.size===e.outgoing)!,inside=e.target?states[e.target-1][e.reference].find(s=>s.size===e.incoming):undefined,cap=Math.min(program.upper[out.column],inside?program.upper[inside.column]:shares[e.reference]);program.upper[e.mass]=Math.min(program.upper[e.mass],cap);for(const i of e.prefixColumns||[])program.upper[i]=Math.min(program.upper[i],cap);}
 for(const d of drafts){d.tolerance=tolerance;d.measurements=allMeasurements.filter(m=>d.positions.includes(m.position));}
 program.zeroCostSubstitutions=drafts.flatMap(d=>d.groups.flatMap(g=>[...g.flow!.values()].flatMap(f=>f.unchanged.map(c=>({column:c.column,terms:c.cap})))));
 const spectralColumns=new Set(program.quadratics.flatMap(q=>q.terms.flatMap(row=>row.map(([i])=>i))));
 program.localBlocks=[...blockColumns].sort((a,b)=>a[0]-b[0]).map(([,columns])=>[...columns].filter(i=>!spectralColumns.has(i)));
 program.originPricing={columns:allEdges.map(e=>e.mass),seed:allEdges.filter(e=>(e.target===0||lengths[e.target-1][e.reference].some(c=>c.size===e.incoming&&c.fraction>0))&&lengths[e.target][e.reference].some(c=>c.size===e.outgoing&&c.fraction>0)).map(e=>e.mass)};
 const phaseSeed=states.flatMap(hs=>hs.flatMap(g=>g.map(s=>({column:s.column,value:seed[s.column]}))));
 marginalModelDiagnostics.observe?.({program,refs,shares,targets:drafts,states,noise,phaseSeed,measurements:allMeasurements,spectralPositions:[...spectralMeasurements.values()].flatMap(list=>list.map(m=>m.position)),departures:allDepartures});
 tick('Fitting the shared global sequence-flow model');
 const solver=new MarginalLinearSolver(program,optimizer);let point:number[]=[],budget=0,representativeSelection:NonNullable<SangerAnalysisResult['diagnostics']>['representativeSelection'];
 try{
  progress('Fitting the shared global sequence-flow model · checking zero-departure witness');const zeroPoint=verifiedZeroDeparturePoint(program,phaseSeed,optimizer);progress('Fitting the shared global sequence-flow model · verifying full-model optimum');let fitted:{x:number[];value:number}|undefined=zeroPoint,shadow:number[]=[],planes:{terms:LinearTerm[];upper:number}[]=[],minimumLower:number|undefined=zeroPoint?0:undefined;
  if(!fitted){const direct=new MarginalLinearSolver({...program,localBlocks:undefined},optimizer);try{direct.setTimeSlice(600);fitted=direct.solveToCompletion();shadow=direct.rowDuals();planes=direct.supportPlanes();minimumLower=direct.objectiveLowerBound();}finally{direct.dispose();}}
  marginalModelDiagnostics.primaryFit?.({point:fitted.x,value:fitted.value,lower:minimumLower,shadow,planes});
  point=fitted.x;for(const draft of drafts)draft.primaryPoint=point;budget=fitted.value+3*noise*Math.sqrt(4*drafts.reduce((n,d)=>n+d.positions.length,0));const budgetRow=program.rows.length;solver.addRow(allDepartures,0,budget);program.dualHint={shadow,budgetRow,planes};
  tick('Selecting parsimonious sequence-flow representative');const selected=selectMarginalRepresentative(program,point,drafts.flatMap(d=>d.groups.filter(g=>g.size!==0).flatMap(g=>g.mass)),optimizer,{substitutionOpportunities:drafts.reduce((n,d)=>n+d.positions.length,0)/drafts.length,noise,phaseSeed,states,shares,edges:allEdges});point=selected.point;representativeSelection=selected.summary;
  for(const m of allMeasurements){const den=evaluate(m.total,point);for(let b=0;b<4;b++){const prediction=evaluate(m.forms[b],point)/den,y=m.observed[b];squaredError+=(y-prediction)**2;signalVariance+=(y-.25)**2;signalRows++;}}
  for(let target=0;target<drafts.length;target++){const draft=drafts[target],{qs,positions,positionForms}=draft,internal=new Map(draft.groups.map(g=>[g.id,g]));draft.point=point;tick('Bounding length proportions · target '+(target+1));
   const records:{g:MarginalInternalGroup;share:number;belowReporting:boolean;weak:boolean;bases:SangerMarginalGroup['bases']}[]=[];
   for(const g of internal.values()){progress('Bounding length proportions · target '+(target+1)+' · '+refs[g.reference].name+' · '+g.size+' bp');g.fraction=Math.max(0,evaluate(g.mass,point));const share=shares[g.reference],belowReporting=g.fraction/share<(input.settings.indelThreshold??2)/100;g.fractionRange=belowReporting?[0,share]:solver.bounds(g.mass);if(g.fractionRange[1]>=1e-8)records.push({g,share,belowReporting,weak:g.fractionRange[1]<=3*noise||belowReporting,bases:[]});}
   tick('Bounding base proportions · target '+(target+1));const localBounds=localBoundCertificates(program,point,draft.groups,optimizer);let fallbackQueries=0,fallbackSeconds=0;
   const certifiedBound=(q:number,terms:LinearTerm[],maximum=false)=>{const bound=localBounds.outerBound(q,terms,maximum);if(bound!==undefined)return bound;const started=performance.now();fallbackQueries++;try{return solver.solve(terms,maximum).value;}finally{fallbackSeconds+=(performance.now()-started)/1000;}};
  try{
   const endpointWitness=marginalEndpointWitnesses(program,point),signatureIds=new Map<string,number>(),signatures=new Map(qs.map(p=>{const key=calibration.lag===0?JSON.stringify([positions.includes(p)?raw[p]:null,positionForms.get(p)]):String(p);let id=signatureIds.get(key);if(id===undefined){id=signatureIds.size;signatureIds.set(key,id);}return[p,id];})),boundCache=new Map<string,[number,number]>();
   // Intern by exact equality, not a probabilistic hash, then release the
   // long descriptions. Group/base cache keys retain only the short ID.
   signatureIds.clear();
   // Process all groups at one coordinate before releasing its workspace.
   // This preserves warm solves without retaining a native LP per coordinate.
   for(const p of qs){progress('Bounding base proportions · target '+(target+1)+' · read '+(p+1));for(const record of records){const{g,weak,bases}=record,ref=refs[g.reference],forms=g.bases.get(p)!,extremes=weak?[false,false,false,false]:endpointWitness(g,p),estimate=forms.map(f=>g.fraction?Math.max(0,evaluate(f,point))/g.fraction:0),lower:number[]=Array(4).fill(0),upper:number[]=Array(4).fill(1),dependent=estimate.indexOf(Math.max(...estimate));
     if(!weak){for(let b=0;b<4;b++){if(b===dependent)continue;const form=forms[b];if(!form.length){lower[b]=0;upper[b]=0;}else if(identical(form,g.mass)){lower[b]=1;upper[b]=1;}else{const key=g.id+':'+b+':'+signatures.get(p);let bound=boundCache.get(key);if(!bound){const lo=extremes.some((yes,a)=>yes&&a!==b)?0:Math.max(0,Math.min(1,((certifiedBound(p,form))-1e-7)/Math.max(1e-30,g.fractionRange[1]))),hi=extremes[b]?1:g.fractionRange[0]>1e-8?Math.min(1,((certifiedBound(p,form,true))+1e-7)/g.fractionRange[0]):1;bound=[lo,hi];boundCache.set(key,bound);}lower[b]=bound[0];upper[b]=bound[1];}}lower[dependent]=Math.max(0,1-upper.reduce((n,v,b)=>n+(b===dependent?0:v),0));upper[dependent]=Math.min(1,1-lower.reduce((n,v,b)=>n+(b===dependent?0:v),0));}
     bases.push({position:ref.reverse?ref.input.sequence.length-ref.offset-p-1:ref.offset+p,estimate,lower,upper});
   }}
   for(const{g,share,belowReporting,weak,bases}of records){const ref=refs[g.reference];
    const result:SangerMarginalGroup={id:g.id,reference:ref.name,target:ref.reverse?ref.cuts.length-target:target+1,size:g.size,fraction:g.fraction/share,fractionRange:[Math.min(1,g.fractionRange[0]/share),Math.min(1,g.fractionRange[1]/share)],bases,scope:'read-cycle-by-length',boundMethod:belowReporting?'deferred-report-envelope':weak?'conservative-noise-envelope':'certified-outer-envelope',status:'bounded'};
    result.representative=fittedMarginalRepresentative(g,point,ref,share,(start,end,inserted)=>marginalInputEvent(ref,start,end,inserted,target));
    if(g.size===0){
     // Fréchet bounds on unchanged molecules: position marginals do not phase
     // multiple substitutions onto the same molecules.
     const shifts=[...new Set(g.edges.map(e=>e.incoming))];if(shifts.length===1)result.baselineReadOffset=ref.reverse?shifts[0]:-shifts[0];if(shifts.length===1&&!weak){let total=0,largest=0;for(let i=0;i<qs.length;i++){const b=DNA.indexOf(ref.sequence[qs[i]-shifts[0]]);if(b<0)continue;total+=g.fractionRange[1]*(1-bases[i].lower[b]);largest=Math.max(largest,g.fractionRange[0]*(1-bases[i].upper[b]));}total=Math.min(total,budget);result.unmodifiedRange=[Math.max(0,(g.fractionRange[0]-total)/share),Math.max(0,Math.min(1,(g.fractionRange[1]-largest)/share))];result.substitutionRange=[largest/share,Math.min(g.fractionRange[1],total)/share];}else{result.unmodifiedRange=[0,Math.min(1,g.fractionRange[1]/share)];result.substitutionRange=[0,Math.min(1,g.fractionRange[1]/share)];}
    }
    if(ref.reverse)result.bases=result.bases.map(p=>({...p,estimate:[...p.estimate].reverse(),lower:[...p.lower].reverse(),upper:[...p.upper].reverse()}));
    groups.push(result);
   }
  }finally{optimization.push({target:target+1,variables:program.upper.length,rows:program.rows.length,coefficients:program.rows.reduce((n,r)=>n+r.terms.length,0),bounds:{localQueries:localBounds.stats.queries,localSeconds:localBounds.stats.seconds,fallbackQueries,fallbackSeconds},solver:structuredClone(solver.diagnostics)});localBounds.dispose();}targets.push(draft);
  }
 }finally{solver.dispose();}
 tick('Bounding substitution signal across length groups');
 const substitutionSummaries=targets.flatMap(target=>boundMarginalSubstitutions(target,refs,shares,optimizer,s=>progress('Bounding substitution signal across length groups · '+s)));
 tick('Separating supported within-length display profiles');
 for(const target of targets)for(const g of target.groups){const display=groups.find(d=>d.id===g.id);if(!display||display.fraction<(input.settings.indelThreshold??2)/100)continue;try{display.displayComponents=withinLengthDisplay({refs,shares,targets,optimizer,pages:new Map()},target,g,display);}catch{/* Retain the complete parent mixture if a display query is incomplete. */}}
 tick('Preparing length and base-range results');
 const annotations=refs.map((r,h)=>({name:r.name,reference:r.input.sequence,cut:r.input.cut,cuts:r.input.cuts,direction:r.reverse?'reverse' as const:'forward' as const,share:shares[h],components:[],marginalGroups:groups.filter(g=>g.reference===r.name),estimable:true})),first=groups.filter(g=>g.target===1),components=Array.from(new Set(first.map(g=>g.size))).sort((a,b)=>a-b).map(size=>({size,fraction:first.filter(g=>g.size===size).reduce((n,g)=>n+g.fraction*shares[refs.findIndex(r=>r.name===g.reference)],0)})),netZero=components.find(c=>c.size===0)?.fraction||0,wtLower=first.filter(g=>g.size===0).reduce((n,g)=>n+(g.unmodifiedRange?.[0]||0)*shares[refs.findIndex(r=>r.name===g.reference)],0),inFrame=components.filter(c=>c.size!==0&&c.size%3===0).reduce((n,c)=>n+c.fraction,0),outOfFrame=components.filter(c=>c.size%3!==0).reduce((n,c)=>n+c.fraction,0);
 const rmse=Math.sqrt(squaredError/Math.max(1,signalRows)),cvRmse=Math.sqrt(heldOutError/Math.max(1,heldOutRows)),passesChecks=calibration.validationR2>=.8&&cvRmse<=Math.max(.005,3*noise);annotations.forEach(a=>a.estimable=passesChecks);
 const result:SangerAnalysisResult={substitutionSummaries,method:'pooled',referenceName:refs.length===1?refs[0].name:'Mixed amplicons',sampleName:input.sample.name,controlName:input.control.name,direction:'forward',cutPosition:Math.min(...refs.flatMap(r=>r.cuts)),alignmentIdentity:0,alignment:[],variants:[],synthetic:!!input.synthetic,plots:[plotTrace(input.control,refs[0].cuts[0]-observation.cr.trim,input.settings.window),plotTrace(input.sample,refs[0].cuts[0]-observation.sr.trim,input.settings.window)],annotations,marginalGroups:groups,exactAlternatives:[],warnings:['Base ranges are coordinate-wise feasible envelopes, not independently selectable haplotypes. Zero net length is not a WT call. Sequence alternatives are calculated only on request.'],diagnostics:{representativeSelection,optimization,referenceSignal:{frame:'control-normalized',robustDelta,controlShares:shares.slice(),fittedShares:shares.map((share,h)=>share*regionalResponse(h,initialization.center)/shares.reduce((sum,share,k)=>sum+share*regionalResponse(k,initialization.center),0)),responseScale:shares.map((_,h)=>regionalResponse(h,initialization.center)),drift:{center:initialization.center,span:initialization.span,parameters:initialization.parameters,trainingOnly:true}},representativeRefit,algorithmVersion:'sanger-marginals-v114-calibrated-yield',signalModel:'measured-control',requestedMaxIndel:input.settings.maxIndel,pipeline:{version:'v114-wt-calibration-uncertainty',observation:'measured-control',thresholdStage:'postprocess'},initialization:{status:'complete',candidates:groups.length},sampleCalibration:calibration,measuredValidation:{scope:'model-supported-spectrum',controlStart:quality.start,controlEnd:quality.end,fitStart:firstMeasured,fitEnd:analysisEnd,qualityAvailable:quality.qualityAvailable,noiseRmse:noise,noiseObservable:'integrated-photon-area',contrastNoiseRmse:Math.sqrt(observation.contrastNoise),predictionLimit:Math.max(.005,3*noise),crossValidatedRmse:cvRmse,crossValidatedR2:heldOutVariance?1-heldOutError/heldOutVariance:1,nullR2:0,meanPredictionGain:0,predictionGainSE:0,unresolvedFraction:0,editingSensitivity:[1-netZero,1-netZero]}},mixture:{components,wt:wtLower,indel:1-netZero,inFrame,outOfFrame,r2:signalVariance?1-squaredError/signalVariance:1,rmse,upstreamError:noise,bases:allMeasurements.length,fitStart:firstMeasured,fitEnd:analysisEnd,converged:true,passesChecks},pooled:{preflight:report,controlR2:calibration.validationR2,sampleTrim:observation.sr.trim,controlTrim:observation.cr.trim,amplicons:refs.map((r,h)=>({name:r.name,share:shares[h],wt:null,indel:null,inFrame:null,outOfFrame:null,components:groups.filter(g=>g.reference===r.name&&g.target===1).map(g=>({size:g.size,fraction:g.fraction}))}))}};
 timings.push({stage,seconds:(performance.now()-last)/1000});result.stageTimings=timings;marginalSessions.set(result,{refs,shares,targets,optimizer,pages:new Map()});return result;
}
