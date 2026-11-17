import { emissionUncertainty } from './sanger-emission-uncertainty';
import { MARGINAL_INDEL_ORIGIN_RADIUS as ORIGIN_RADIUS,MARGINAL_SUBSTITUTION_RADIUS as SUB_RADIUS,marginalReferenceSubstitutionAllowed,marginalTargetReadSupport } from './sanger-marginal-domain';
import { LP_INFINITY,type LinearTerm,type MarginalProgram } from './sanger-marginal-lp';
import { simultaneousResidualLimit } from './sanger-marginal-noise';
import { marginalModelDiagnostics,type MarginalEdge,type MarginalInternalGroup,type MarginalMeasurement,type MarginalTargetState } from './sanger-marginal-state';
import { prepareMeasuredObservation } from './sanger-measured-observation';
import { originPrefixes } from './sanger-origin-prefix';
import { packLinearTerms } from './sanger-packed-terms';
import type { SignalPathInput } from './sanger-path-common';
import { prepare } from './sanger-path-common';
import { rebasePhotonForms } from './sanger-photon-state-rebase';
import { seedQuadraticTangents } from './sanger-quadratic-envelope';
import { calibrateRegionalResponse } from './sanger-regional-calibration';
import { regionalLengthInitialization,type RegionalLengthWindow } from './sanger-regional-length-fit';
import { mixPhase,phaseNeighbor } from './sanger-signal';
import { spectralBoxCertificates } from './sanger-spectral-box-certificate';
const DNA='ACGT';
const evaluate=(terms:LinearTerm[],x:number[])=>terms.reduce((n,[i,w])=>n+w*x[i],0);
const sumTerms=(lists:LinearTerm[][])=>{const m=new Map<number,number>();for(const list of lists)for(const[i,w]of list)m.set(i,(m.get(i)||0)+w);return[...m].filter(([,w])=>Math.abs(w)>1e-14)as LinearTerm[];};
const identical=(a:LinearTerm[],b:LinearTerm[])=>{const x=new Map(a),y=new Map(b);return x.size===y.size&&[...x].every(([i,v])=>y.get(i)===v);};

/** Shared physical model construction. No GT, allele count or reporting policy. */
export function buildMarginalModel(input:SignalPathInput,tick:(stage:string)=>void){
 const {refs,report,maxDeletion,maxInsertion}=prepare(input.amplicons,input.primer,input.settings),observation=prepareMeasuredObservation(refs,input.control,input.sample,input.settings,s=>tick(s)),{raw,means,calibration,quality}=observation,shares=observation.shares.slice(),noise=Math.sqrt(observation.noise),analysisEnd=Math.min(observation.end,observation.c.length-maxDeletion*refs[0].cuts.length-1,quality.end-maxDeletion*refs[0].cuts.length-1),lengths:{size:number;fraction:number}[][][]=[],count=refs[0].cuts.length;let heldOutError=0,heldOutVariance=0,heldOutRows=0;
 let regionalResponse=(h:number,p:number)=>1;
 const emissions=observation.intensity.responses,emissionMeans=observation.intensity.novelByReference,cycleEnvelope=observation.intensity.cycleEnvelope,dyeAt=(h:number,p:number,dye:number[])=>observation.transferByReference.apply(h,p,dye).map(v=>v*regionalResponse(h,p)*cycleEnvelope(h,p)),emit=(h:number,shift:number,p:number)=>dyeAt(h,p,emissions[h][p-shift]?.map(v=>v/cycleEnvelope(h,p-shift))||emissionMeans[h][DNA.indexOf(refs[h].sequence[p-shift])]||[0,0,0,0]);
 function regionalFit(
  observation: ReturnType<typeof prepareMeasuredObservation>,
  shares: number[],
  analysisEnd: number,
  emit: (reference: number, shift: number, position: number) => number[],
 ) {
  const { raw, calibration } = observation;
 const flankPositions:number[][]=[],flowWindows:RegionalLengthWindow[]=[];
 for(let t=0;t<count;t++){
  const left=Math.max(...refs.map(r=>r.cuts[t]))+Math.max(SUB_RADIUS-1+(t+1)*maxInsertion,(t+1)*maxInsertion+ORIGIN_RADIUS)+5,right=Math.min(analysisEnd,...refs.map(r=>t+1<count?r.cuts[t+1]-Math.max(SUB_RADIUS+(t+1)*maxDeletion,(t+2)*maxDeletion+ORIGIN_RADIUS)-5:analysisEnd)),ps=Array.from({length:Math.max(0,right-left)},(_,i)=>left+i);
  // Hypothesis count alone cannot establish identifiability under known
  // control spectra, nonnegativity and the complete joint constraints.
  // Require independent flank validation; preserve every requested class and
  // let original-model certificates, predictive QC and bounds decide support.
  // The regional initializer checks actual disjoint training and validation
  // blocks; total flank length is not an independent accuracy criterion.
  flankPositions.push(ps);flowWindows.push({classes:refs.flatMap((_,h)=>Array.from({length:(t+1)*(maxDeletion+maxInsertion)+1},(_,i)=>({reference:h,size:i-(t+1)*maxDeletion}))).map(c=>({...c,column:ps.flatMap(p=>mixPhase(emit(c.reference,c.size,p),emit(c.reference,c.size,phaseNeighbor(p,calibration.lag)),calibration.lag))})),observed:ps.flatMap(p=>raw[p]),background:ps.flatMap(p=>observation.backgroundAt(p)),weights:ps.map(p=>observation.noise/observation.noiseAt(p)**2)});
 }
 const preliminary=regionalLengthInitialization(flowWindows,flankPositions,shares,observation.fitDelta,observation.transferByReference),initialization=calibrateRegionalResponse(flowWindows,flankPositions,shares,preliminary.delta,{positions:observation.trainingPositions,columns:refs.map((_,h)=>observation.trainingPositions.flatMap(p=>mixPhase(emit(h,0,p),emit(h,0,phaseNeighbor(p,calibration.lag)),calibration.lag))),observed:observation.trainingPositions.flatMap(p=>raw[p]),background:observation.trainingPositions.flatMap(p=>observation.backgroundAt(p)),weights:observation.trainingPositions.map(p=>observation.noise/observation.noiseAt(p)**2)});
  return { flankPositions, flowWindows, initialization };
 }
 tick('Fitting invariant-flank length proportions');
 const { flankPositions, flowWindows, initialization } = regionalFit(
  observation, shares, analysisEnd, emit,
 );
 const robustDelta = initialization.delta;
 // Reference brightness is a nuisance distinct from control molecule share.
 // Its flank estimate uses a paired-prefix predictive prior, whose variance
 // includes the actual ADC floor. The same correction applies to every input.
 regionalResponse=initialization.at;
 lengths.push(...initialization.lengths);
 for(let t=0;t<count;t++)for(let i=0;i<flankPositions[t].length;i++)if(Math.floor(flankPositions[t][i]/5)%3===2)for(let b=0;b<4;b++){const y=raw[flankPositions[t][i]][b],prediction=initialization.trainingPredictions[t][4*i+b];heldOutError+=(y-prediction)**2;heldOutVariance+=(y-.25)**2;heldOutRows++;}
 const program:MarginalProgram={upper:[],cost:[],rows:[]},variable=(upper:number,cost=0)=>{const id=program.upper.length;program.upper.push(upper);program.cost.push(cost);return id;},row=(terms:LinearTerm[],lower:number,upper:number)=>program.rows.push({terms:packLinearTerms(sumTerms([terms])),lower,upper}),states=flowWindows.map(w=>refs.map((_,h)=>w.classes.flatMap(c=>c.reference===h?[{size:c.size,fraction:shares[h],column:variable(shares[h])}]:[]))),drafts:(MarginalTargetState&{qs:number[];positionForms:Map<number,unknown[]>})[]=[],photons=new Map<number,LinearTerm[][]>(),intensityFloor=new Map<number,number>(),allDepartures:LinearTerm[]=[],allEdges:MarginalEdge[]=[],blockColumns=new Map<number,Set<number>>();
 const uncertainEmission=emissionUncertainty(program,shares,refs.map(r=>r.sequence),emissions,observation.intensity.allocationBounds,blockColumns);
 for(let t=0;t<count;t++)for(let h=0;h<refs.length;h++)row(states[t][h].map(c=>[c.column,1]),shares[h],shares[h]);
 for(let target=0;target<count;target++){
  tick('Constructing joint marginal constraints · target '+(target+1));
  const incoming=refs.map((_,h)=>target?states[target-1][h]:[{size:0,fraction:shares[h],column:undefined as number|undefined}]),outgoing=states[target],supports=refs.map((r,h)=>marginalTargetReadSupport(r.cuts[target],Math.min(...incoming[h].filter(c=>c.fraction>0).map(c=>c.size)),Math.max(...incoming[h].filter(c=>c.fraction>0).map(c=>c.size)),maxDeletion,maxInsertion)),lo=Math.max(quality.start,Math.min(...supports.map(s=>s.first))-1),hi=Math.min(analysisEnd,Math.max(...supports.map(s=>s.last))+2),nominalHi=Math.min(analysisEnd,Math.max(...supports.map(s=>s.nominalLast))+2),mutationPriorOpportunities=Math.max(0,nominalHi-lo-2),qs=Array.from({length:hi-lo},(_,i)=>lo+i),positions=qs.slice(1,-1),tolerance=simultaneousResidualLimit(noise,4*positions.length*count),edges:MarginalEdge[]=[],measurements:MarginalMeasurement[]=[],minimumIntensity=new Map<number,number>(),internal=new Map<string,MarginalInternalGroup>();
  const nominalPhotons=new Map<string,Map<number,LinearTerm[][]>>();
  const id=(h:number,size:number)=>'marginal:'+h+':'+target+':'+size;
  for(let h=0;h<refs.length;h++)for(const a of incoming[h])if(a.fraction>0)for(const b of outgoing[h])if(b.fraction>0){
   const size=b.size-a.size;if(size< -maxDeletion||size>maxInsertion)continue;const cut=refs[h].cuts[target],startMin=size<0?cut+size-ORIGIN_RADIUS:size>0?cut-ORIGIN_RADIUS:cut,startMax=size?cut+ORIGIN_RADIUS:cut,seen=new Set<string>();
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
   const map=new Map<string,{key:string;base:number;dye:number[];cap:LinearTerm[];donor?:number;mutable:boolean}>();
   const add=(base:number,dye:number[],cap:LinearTerm[],donor?:number)=>{cap=sumTerms([cap]);if(!cap.length)return;const mutable=marginalReferenceSubstitutionAllowed(donor,refs[g.reference].cuts[target]),key=base+':'+dye.join(',')+':mutable-'+Number(mutable)+(donor!==undefined&&observation.intensity.allocationBounds[g.reference]?.[donor]?.some(v=>Math.abs(v)>1e-12)?':'+donor:''),old=map.get(key);if(old)old.cap=sumTerms([old.cap,cap]);else map.set(key,{key,base,dye,cap,donor,mutable});};
   for(const bucket of geometry.get(g)!){const before=sumTerms([[[bucket.total,1]],bucket.at(p-bucket.incoming).map(([i,w])=>[i,-w]as LinearTerm)]),after=bucket.at(p-bucket.incoming-Math.max(0,g.size));for(const[at,cap]of [[p-bucket.incoming,before],[p-bucket.outgoing,after]]as const)if(cap.length){const b=DNA.indexOf(refs[g.reference].sequence[at]);if(b<0)throw Error('Resolve ambiguous reference bases in the modeled target window.');add(b,emit(g.reference,p-at,p),cap,at);}}
   const free=insertionFlows.get(g)!.get(p);if(free)for(let b=0;b<4;b++)add(b,dyeAt(g.reference,p,emissionMeans[g.reference][b]),[[free[b],1]]);
   const contexts=[...map.values()].sort((a,b)=>a.key.localeCompare(b.key));minimumIntensity.set(p,Math.min(minimumIntensity.get(p)??Infinity,...contexts.map(c=>c.dye.reduce((n,v)=>n+v,0)*(c.donor===undefined?1:1+(observation.intensity.allocationBounds[g.reference]?.[c.donor]?.[0]||0))),...emissionMeans[g.reference].map(row=>dyeAt(g.reference,p,row).reduce((n,v)=>n+v,0))));contexts.forEach(c=>c.cap=sumTerms([c.cap]));const expected=Array.from({length:4},(_,b)=>sumTerms(contexts.filter(c=>c.base===b).map(c=>c.cap))),subAllowed=contexts.some(c=>c.mutable),signatures=positionForms.get(p)||[];signatures.push([subAllowed,contexts.map(c=>[c.base,c.dye,c.cap])]);positionForms.set(p,signatures);
   const nominal=Array.from({length:4},(_,b)=>sumTerms(contexts.map(c=>c.cap.map(([i,w])=>[i,w*c.dye[b]]as LinearTerm)))),saved=nominalPhotons.get(g.id)||new Map<number,LinearTerm[][]>();saved.set(p,nominal);nominalPhotons.set(g.id,saved);
   const correction=(cap:LinearTerm[],c:typeof contexts[number])=>c.donor===undefined?[]:uncertainEmission(cap,g.reference,c.donor,p);
   if(!subAllowed){g.bases.set(p,expected);g.dye!.set(p,Array.from({length:4},(_,b)=>sumTerms([nominal[b],...contexts.map(c=>correction(c.cap,c).map(([i,w])=>[i,w*c.dye[b]]as LinearTerm))])));continue;}
   const contextRowsStart=program.rows.length;const pure=contexts.length===1&&identical(contexts[0].cap,g.mass),unchanged=contexts.map(()=>variable(shares[g.reference])),changed=Array.from({length:4},(_,b)=>variable(pure&&b===contexts[0].base?0:shares[g.reference],1));if(!pure)contexts.forEach((c,j)=>row([[unchanged[j],1],...c.cap.map(([i,w])=>[i,-w]as LinearTerm)],c.mutable?-LP_INFINITY:0,0));
   const massRow=program.rows.length;row([...unchanged.map(i=>[i,1]as LinearTerm),...changed.map(i=>[i,1]as LinearTerm),...g.mass.map(([i,w])=>[i,-w]as LinearTerm)],0,0);
   // A substitution cannot keep the same nucleotide. These four Hall bounds
   // are equivalent to existence of a nonnegative off-diagonal 4x4 transport,
   // eliminating twelve transition variables without losing any base flow.
   // The mass balance eliminates the long common-mass expression from
   // Hall's inequalities. A single immutable source context needs only the
   // mass balance and a zero bound on its forbidden diagonal substitution.
   if(!pure)for(let b=0;b<4;b++)row([...expected[b],...contexts.flatMap((c,j)=>c.base===b?[[unchanged[j],-1]as LinearTerm]:[]),...changed.flatMap((i,a)=>a!==b?[[i,-1]as LinearTerm]:[])],-LP_INFINITY,0);
   marginalModelDiagnostics.transport?.({mutableSources:[...new Set(contexts.filter(c=>c.mutable).map(c=>c.base))],changed,massRow,capRows:contexts.map((_,j)=>pure?undefined:contextRowsStart+j),unchanged:contexts.map((c,j)=>({base:c.base,column:unchanged[j],cap:c.cap,mutable:c.mutable})),halls:pure?[]:[0,1,2,3].map(base=>({base,row:massRow+1+base})),insertionFirst:free?.[0]});
   g.flow!.set(p,{unchanged:contexts.map((c,j)=>({base:c.base,column:unchanged[j],cap:c.cap})),changed});
   const norm=changed.map(i=>[i,1]as LinearTerm);departures.push(...norm);g.departures!.set(p,norm);g.bases.set(p,Array.from({length:4},(_,b)=>[[changed[b],1]as LinearTerm,...contexts.flatMap((c,j)=>c.base===b?[[unchanged[j],1]as LinearTerm]:[])]));g.dye!.set(p,Array.from({length:4},(_,b)=>sumTerms([contexts.map((c,j)=>[unchanged[j],c.dye[b]]as LinearTerm),...contexts.map((c,j)=>correction([[unchanged[j],1]],c).map(([i,w])=>[i,w*c.dye[b]]as LinearTerm)),changed.map((i,a)=>[i,dyeAt(g.reference,p,emissionMeans[g.reference][a])[b]]as LinearTerm)])));
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
  drafts.push({program,groups:[...internal.values()],edges,point:[],positions,mutationPriorOpportunities,mutationNoiseVarianceSum:qs.slice(1,Math.max(1,nominalHi-lo-1)).reduce((n,p)=>n+observation.noiseAt(p)**2,0),tolerance,measurements,qs,positionForms});
 }
 // A target's bounding interval can contain unedited gaps between references.
 // Those cycles are invariant-flank evidence, not cut-proximal SNP measurements.
 // Retain their exact origin-flow photons and the same flank RSS confidence rule.
 const firstMeasured=Math.min(...drafts.map(d=>d.positions[0])),allMeasurements:MarginalMeasurement[]=[],localMeasurements:MarginalMeasurement[]=[],localReadCycles=new Set(drafts.flatMap((d,t)=>d.positions.filter(p=>refs.some((r,h)=>{const incoming=t?states[t-1][h].map(s=>s.size):[0];const support=marginalTargetReadSupport(r.cuts[t],Math.min(...incoming),Math.max(...incoming),maxDeletion,maxInsertion);return p>=support.first&&p<=support.last;})))),spectralMeasurements=new Map<number,MarginalMeasurement[]>();
 const invariant=(q:number)=>{let t=-1;for(let j=0;j<drafts.length;j++)if(q>drafts[j].positions.at(-1)!)t=j;if(t<0)throw Error('No calibrated cumulative phase is available for the requested measurement.');const forms:LinearTerm[][]=[[],[],[],[]];let floor=Infinity;for(let h=0;h<refs.length;h++)for(const state of states[t][h]){const dye=emit(h,state.size,q);floor=Math.min(floor,dye.reduce((a,b)=>a+b,0)*(1+(observation.intensity.allocationBounds[h]?.[q-state.size]?.[0]||0)));const correction=uncertainEmission([[state.column,1]],h,q-state.size,q);for(let b=0;b<4;b++)forms[b].push([state.column,dye[b]],...correction.map(([i,w])=>[i,w*dye[b]]as LinearTerm));}photons.set(q,forms);intensityFloor.set(q,floor);return forms;};
 for(let p=firstMeasured;p<analysisEnd;p++){let unitTarget=0;for(let j=0;j<drafts.length;j++)if(p>drafts[j].positions.at(-1)!)unitTarget=j;const unit=states[unitTarget].flatMap(g=>g.map(c=>[c.column,1]as LinearTerm)),background=observation.backgroundAt(p),neighbor=phaseNeighbor(p,calibration.lag),spectral=!localReadCycles.has(p)&&(!calibration.lag||!localReadCycles.has(neighbor)),forms=[0,1,2,3].map(b=>sumTerms([[p,1-Math.abs(calibration.lag)],[neighbor,Math.abs(calibration.lag)]].filter(([,w])=>w>0).map(([q,w])=>(photons.get(q)||invariant(q))[b].map(([i,v])=>[i,v*w]as LinearTerm)).concat([unit.map(([i,v])=>[i,v*background[b]]as LinearTerm)]))),total=sumTerms(forms),floor=Math.min(intensityFloor.get(p)!,intensityFloor.get(neighbor)??intensityFloor.get(p)!)+background.reduce((a,b)=>a+b,0);if(!(floor>0))throw Error('Missing positive measured donor fluorescence.');row(total,floor,LP_INFINITY);const m={position:p,forms,total,observed:raw[p],noise:observation.noiseAt(p)};allMeasurements.push(m);if(spectral){let t=0;for(let j=0;j<drafts.length;j++)if(p>drafts[j].positions.at(-1)!)t=j;const list=spectralMeasurements.get(t)||[];list.push(m);spectralMeasurements.set(t,list);}else localMeasurements.push(m);}
 const tolerance=simultaneousResidualLimit(noise,4*localMeasurements.length);for(const m of localMeasurements)for(let b=0;b<4;b++){const tolerance=simultaneousResidualLimit(m.noise!,4*localMeasurements.length);row([...m.forms[b],...m.total.map(([i,w])=>[i,-w*(m.observed[b]-tolerance)]as LinearTerm)],0,LP_INFINITY);row([...m.forms[b],...m.total.map(([i,w])=>[i,-w*(m.observed[b]+tolerance)]as LinearTerm)],-LP_INFINITY,0);}
 // A regression confidence region uses all raw flank residuals. Unlike an
 // unjustified per-peak maximum, it permits the isolated errors present in
 // measured reads without deleting them or changing the held-out QC rule.
 const seed=program.upper.map(()=>0);for(let t=0;t<count;t++)for(let h=0;h<refs.length;h++)for(const state of states[t][h])seed[state.column]=lengths[t][h].find(c=>c.size===state.size)!.fraction;
 uncertainEmission.seed(seed);
 program.quadratics=[...spectralMeasurements].map(([t,list])=>{const terms=list.flatMap(m=>{const denominator=Math.max(1e-30,evaluate(m.total,seed));return m.forms.map((form,b)=>sumTerms([form,m.total.map(([i,w])=>[i,-w*m.observed[b]]as LinearTerm)]).map(([i,w])=>[i,w/denominator*noise/m.noise!]as LinearTerm));}),rss=terms.reduce((n,r)=>n+evaluate(r,seed)**2,0),parameters=Math.max(1,states[t].reduce((n,g)=>n+g.length-1,0)),variance=Math.max(observation.noise,rss/Math.max(1,terms.length-parameters)),tail=Math.log(Math.max(1,spectralMeasurements.size)/.0027),allowance=variance*(parameters+2*Math.sqrt(parameters*tail)+2*tail);return{name:'invariant-spectrum-'+(t+1),terms,upper:rss+allowance};});
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
 // Center-weighted measurements inform only the representative length prior.
 // Keep the integrated-area calibration and constraints as the physical domain.
 tick('Fitting peak-centered flank proportions for the representative');
 const centered = prepareMeasuredObservation(refs, input.control, input.sample, input.settings, () => {}, 0.5);
 const centeredEnd = Math.min(centered.end, centered.c.length - maxDeletion * count - 1,
  centered.quality.end - maxDeletion * count - 1);
 const centeredEmission = (reference: number, shift: number, position: number) => {
  const { responses, novelByReference, cycleEnvelope } = centered.intensity;
  const donor = position - shift;
  const level = responses[reference][donor]?.map(value => value / cycleEnvelope(reference, donor))
   || novelByReference[reference][DNA.indexOf(refs[reference].sequence[donor])] || [0, 0, 0, 0];
  return centered.transferByReference.apply(reference, position, level)
   .map(value => value * cycleEnvelope(reference, position));
 };
 const centeredLengths = regionalFit(centered, centered.shares, centeredEnd, centeredEmission).initialization.lengths;
 const phaseSeed = states.flatMap((byRef, target) => byRef.flatMap((group, reference) =>
  group.map(state => {
   const share = centered.shares[reference];
   const estimate = centeredLengths[target][reference].find(length => length.size === state.size);
   const value = share > 0 && estimate
    ? estimate.fraction * (shares[reference] / share)
    : seed[state.column];
   return {column: state.column, value};
  })));

 marginalModelDiagnostics.observe?.({program,refs,shares,targets:drafts,states,noise,phaseSeed,measurements:allMeasurements,spectralPositions:[...spectralMeasurements.values()].flatMap(list=>list.map(m=>m.position)),departures:allDepartures});
 return {refs,report,maxDeletion,maxInsertion,observation,raw,calibration,quality,shares,noise,analysisEnd,count,lengths,regionalResponse,initialization,robustDelta,program,drafts,states,phaseSeed,allMeasurements,allDepartures,allEdges,uncertainEmission,firstMeasured,heldOutError,heldOutVariance,heldOutRows};
}
export type MarginalModel=ReturnType<typeof buildMarginalModel>;
