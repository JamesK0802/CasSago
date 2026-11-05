import {linearSequenceOracle} from './sanger-linear-sequence-oracle';
import {discoverExactAlternatives} from './sanger-exact-alternatives';
import {substitutionLinkageAliases} from './sanger-substitution-linkage';
import {compoundBaseProposals} from './sanger-compound-refine';
import {SangerFitCache,rawMixtureMse} from './sanger-fit-cache';
import {CandidateSignalScreen} from './sanger-candidate-screen';
import {polishSignalDictionary} from './sanger-dictionary-polish';
import {substitutionProposals} from './sanger-substitution-search';
import {normalizeMoleculeEdits} from './sanger-edit-normalize';
import {repairSignalTraceCycles} from './sanger-cycle';
import {SangerAnalysisResult,SangerJointEvent,SangerTrace,TRACE_BASES} from '../../models/sanger-analysis';
import {reverseComplement} from '../../utils/biology.utils';
import {registerTrace} from './sanger-pooled';
import {contrastPeaks,integratedPeaks,calibrateSignal,correctGain,gramCondition,rawSignalScore} from './sanger-signal';
import {emissionVector,fitContextEmission,fitNormalizedEmissionMixture as fitEmissionMixtureOriginal} from './sanger-emission';
import {plotTrace} from './sanger-analysis';
const fitEmissionMixture=(columns:number[][],observed:number[])=>fitEmissionMixtureOriginal(columns,observed);
import {SignalPathInput,PathEdit,Model,materialize,prepare,assertSignalTraceGrid,eventsRespectIndelOrigin} from './sanger-path-common';
/** Whole-molecule candidate discovery and joint fluorescence refit. These are
 * best-fit cut-local hypotheses, not a posterior over physical haplotypes.
 * All measured positions enter discovery, so no independent held-out score is claimed. */
export function analyzeAlleleGroupsOracle(input:SignalPathInput,progress:(s:string)=>void=()=>{},audit?:(state:any)=>void):SangerAnalysisResult{
 if(!input.ampliconConfirmed||!input.control||!input.sample)throw Error('Confirm submitted products and provide matched control/sample chromatograms.');
 const controlRepair=repairSignalTraceCycles(input.control),sampleRepair=repairSignalTraceCycles(input.sample);input={...input,control:controlRepair.trace,sample:sampleRepair.trace};
 assertSignalTraceGrid(input.control);assertSignalTraceGrid(input.sample);
 const {refs,report,anchor,maxDeletion,maxInsertion}=prepare(input.amplicons,input.primer,input.settings),s=input.settings;
 progress('Registering traces to the sequencing primer');const registered=refs.map(r=>({name:r.name,sequence:r.sequence,cut:Math.min(...r.cuts)}));
 const cr=registerTrace(input.control,registered,s,integratedPeaks(input.control,[0,0,0,0]).vectors),sr=registerTrace(input.sample,registered,s,integratedPeaks(input.sample,[0,0,0,0]).vectors),ctl=contrastPeaks(input.control).vectors,raw=contrastPeaks(input.sample).vectors;
 const calibrate=(vectors:number[][],trim:number)=>{const obs=Array.from({length:vectors.length+trim},(_,p)=>vectors[p-trim]||[0,0,0,0]),rows=Array.from({length:Math.min(...refs.flatMap(r=>r.cuts))-12},(_,p)=>p).filter(p=>p>=Math.max(s.readStart,trim+10)&&p-trim<vectors.length-2);return calibrateSignal(obs,rows,refs.map(r=>Array.from(r.sequence,(_,p)=>TRACE_BASES.map(b=>Number(r.sequence[p]===b)))));};
 const cc=calibrate(ctl,cr.trim),sc=calibrate(raw,sr.trim);if(cc.validationR2<.8||sc.validationR2<.8)throw Error('Signal-path upstream calibration failed.');
 const end=Math.min(Math.max(...report.positions)+1,cr.trim+ctl.length-3,sr.trim+raw.length-3),positions=report.positions.filter(p=>p<end);if(positions.length<120)throw Error('Insufficient actual trace coverage for path inference.');
 const wt=ctl.map(v=>correctGain(v,cc.gains)),sample=raw.map(v=>correctGain(v,sc.gains)),controlObs=Array.from({length:wt.length+cr.trim},(_,p)=>wt[p-cr.trim]||[0,0,0,0]);
 const context=fitContextEmission(refs.map(r=>r.sequence),positions,controlObs,cc.weights,cc.lag),theta=context.theta;
 const emit=(sequence:string,p:number,lag:number)=>{const v=emissionVector(sequence,p,lag,theta);if(sequence[p]==='N')for(let b=0;b<4;b++)v[b]+=(1-Math.abs(lag))/4;return v;};
 const models:Model[]=refs.map((r,h)=>({reference:h,sequence:r.sequence,events:[],aliases:[{reference:h,events:[]}]})),cols=models.map(m=>positions.flatMap(p=>emit(m.sequence,p,sc.lag))),controlCols=models.map(m=>positions.flatMap(p=>emit(m.sequence,p,cc.lag)));
 const ctlFit=fitEmissionMixture(controlCols,positions.flatMap(p=>wt[p-cr.trim])),ctlScore=rawSignalScore(controlCols,positions.flatMap(p=>ctl[p-cr.trim]),ctlFit.weights,cc.gains);if(ctlScore.r2<.85)throw Error('WT product references do not explain the unedited control.');
 const y=positions.flatMap(p=>sample[p-sr.trim]),train=positions.map((_,i)=>i),tune=train,test=train,subset=(v:number[],rows:number[])=>rows.flatMap(i=>v.slice(4*i,4*i+4));
 const positionIndex=new Int32Array(end+maxInsertion+2).fill(-1);positions.forEach((p,i)=>positionIndex[p]=i);
 const rawY=positions.flatMap(p=>raw[p-sr.trim]),fitCache=new SangerFitCache(cols,y);
 const fitIds=(ids:number[])=>fitCache.fit(ids);
 const score=(ids:number[],weights:number[],rows:number[])=>rawSignalScore(ids.map(id=>rows===train?cols[id]:subset(cols[id],rows)),rows===train?rawY:rows.flatMap(i=>raw[positions[i]-sr.trim]),weights,sc.gains);
 const objective=(ids:number[],weights:number[])=>rawMixtureMse(ids.map(id=>cols[id]),rawY,weights,sc.gains)+.0000001*ids.reduce((sum,id,h)=>sum+weights[h]*models[id].events.reduce((a,e)=>a+1+.05*e.inserted.length,0),0);
 let ids=refs.map((_,i)=>i),fit=fitIds(ids),best=objective(ids,fit.weights),iterations=0,alignments=0,substitutionRounds=0,insertionPolishes=0,discoveryLimited=false;
 const multiTarget=refs.some(r=>r.cuts.length>1),screen=new CandidateSignalScreen(cols);

 // Experimental fully corrective column generation. No truth or recipe input.
 if(context.used)throw Error('Unit-emission oracle prototype excludes context-dependent yields.');
 const signature=(sequence:string)=>sequence.slice(positions[0],positions.at(-1)!+1).padEnd(positions.length,'-');
 const seenEdits=refs.map(()=>new Set<string>()),signatures=new Map(models.map((m,id)=>[signature(m.sequence),id]));
  const addVariant=(h:number,events:PathEdit[])=>{events=normalizeMoleculeEdits(refs[h].sequence,events);const editKey=JSON.stringify(events);if(seenEdits[h].has(editKey))return;seenEdits[h].add(editKey);if(!eventsRespectIndelOrigin(events,refs[h].cuts,s.editStartRadius))return;const sequence=materialize(refs[h].sequence,events);if(!sequence)return;const key=signature(sequence),old=signatures.get(key),alias={reference:h,events};if(old!==undefined){if(!models[old].aliases.some(a=>a.reference===h&&JSON.stringify(a.events)===JSON.stringify(events)))models[old].aliases.push(alias);return;}const id=models.length;models.push({reference:h,sequence,events,aliases:[alias]});cols.push(positions.flatMap(p=>emit(sequence,p,sc.lag)));signatures.set(key,id);};

 let complete=false;
 for(let round=0;round<256;round++){
  progress(`Optimizing sequence mixture · round ${round+1}`);
  const cost=Array.from({length:end+2*maxInsertion+4},()=>[0,0,0,0]);
  for(let k=0;k<positions.length;k++){
   const value=[0,0,0,0];for(let j=0;j<ids.length;j++)for(let b=0;b<4;b++)value[b]+=fit.weights[j]*cols[ids[j]][4*k+b];
   const total=value.reduce((a,b)=>a+b,0);if(!total)continue;const prediction=value.map(v=>v/total),residual=prediction.map((v,b)=>v-y[4*k+b]),center=prediction.reduce((sum,v,b)=>sum+v*residual[b],0),p=positions[k],q=p+(sc.lag<0?1:-1);
   for(let b=0;b<4;b++){const gradient=(residual[b]-center)/total;cost[p][b]+=(1-Math.abs(sc.lag))*gradient;if(q>=0&&q<cost.length)cost[q][b]+=Math.abs(sc.lag)*gradient;}
  }
  const proposed:number[]=[];let maximumGap=0;
  for(let h=0;h<refs.length;h++){
   const oracle=linearSequenceOracle({reference:refs[h].sequence,cuts:refs[h].cuts,cost,radius:s.editStartRadius??10,maxDeletion,maxInsertion,maxSubstitution:10});maximumGap=Math.max(maximumGap,-oracle.cost/y.length);
   if(-oracle.cost/y.length<=1e-10)continue;addVariant(h,oracle.events);const id=signatures.get(signature(oracle.sequence));if(id!==undefined&&!ids.includes(id))proposed.push(id);
  }
  if(!proposed.length){complete=maximumGap<=1e-9;break;}
  ids=[...new Set([...ids,...proposed])];fit=fitIds(ids);ids=ids.filter((id,j)=>id<refs.length||fit.weights[j]>1e-9);fit=fitIds(ids);iterations++;
 }
 discoveryLimited=!complete;
 progress('Validating fitted target profiles');
 const qc=score(ids,fit.weights,test),full=qc,geometry=gramCondition(ids.map(id=>cols[id]));
 const net=(events:PathEdit[])=>events.reduce((sum,e)=>sum+e.inserted.length-(e.end-e.start),0),inputEvents=(h:number,events:PathEdit[]):SangerJointEvent[]=>{const r=refs[h];return events.map(e=>{const distances=r.cuts.map(c=>Math.min(Math.abs(e.start-c),Math.abs(e.end-c))),near=distances.indexOf(Math.min(...distances)),span=r.cuts.length===2&&e.start<=r.cuts[0]&&e.end>=r.cuts[1],target=span?0:Math.min(...distances)>60?-1:(r.reverse?r.cuts.length-near:near+1),a=r.offset+e.start,b=r.offset+e.end;return {target,start:r.reverse?r.input.sequence.length-b:a,end:r.reverse?r.input.sequence.length-a:b,inserted:e.inserted.length,insertedSequence:r.reverse?reverseComplement(e.inserted):e.inserted};}).sort((a,b)=>a.start-b.start);};
 const pathSizes=(h:number,events:PathEdit[])=>{const sizes=Array(refs[h].cuts.length).fill(0);for(const e of inputEvents(h,events)){if(e.target>0)sizes[e.target-1]+=e.inserted-(e.end-e.start);}return sizes;};
 const chosen=ids.flatMap((id,i)=>fit.weights[i]>1e-9?[{model:models[id],fraction:fit.weights[i]}]:[]),aliasKey=(a:{reference:number;events:PathEdit[]})=>a.reference+':'+JSON.stringify(inputEvents(a.reference,a.events));
 const ambiguous=chosen.filter(c=>new Set(c.model.aliases.map(a=>a.reference+':'+net(a.events)+':'+Number(a.events.length===0)+':'+[...new Set(inputEvents(a.reference,a.events).map(e=>e.target))].sort().join('/'))).size>1),refAmbiguous=ambiguous.filter(c=>new Set(c.model.aliases.map(a=>a.reference)).size>1);
 const semantic=(a:{reference:number;events:PathEdit[]})=>a.reference+':'+(a.events.length?'edit':'WT')+':'+net(a.events)+':'+[...new Set(inputEvents(a.reference,a.events).map(e=>e.target))].sort().join('/');
 const attributionAmbiguous=ambiguous.filter(c=>new Set(c.model.aliases.map(semantic)).size>1);
 const groups=chosen.map((c,i)=>({id:`P${i+1}`,reference:refAmbiguous.includes(c)?null:refs[c.model.reference].name,fraction:c.fraction,netSize:new Set(c.model.aliases.map(a=>net(a.events))).size>1?null:net(c.model.events),pathSizes:pathSizes(c.model.reference,c.model.events),events:ambiguous.includes(c)?[]:inputEvents(c.model.reference,c.model.events),alternatives:ambiguous.includes(c)?c.model.aliases.map(a=>({reference:refs[a.reference].name,events:inputEvents(a.reference,a.events)})):[]}));
 const substitutionLinkage=substitutionLinkageAliases(refs,chosen.map(c=>({...c.model,fraction:c.fraction})),sequence=>positions.flatMap(p=>emit(sequence,p,sc.lag)));
 const warnings:string[]=[multiTarget?'Best-fit molecules allow one contiguous edit per target within cut ±10 bp. Target profiles are reported separately; molecule linkage is unresolved.':'Best-fit sequence groups assume one contiguous edit per molecule within cut ±10 bp. A single mixed read does not establish unique allele linkage.'];if(s.editStartRadius!==undefined)warnings.push(`Indel initiation constrained to cut ±${s.editStartRadius} bp; off-target-window edits are outside this model.`);if(controlRepair.added+controlRepair.removed+sampleRepair.added+sampleRepair.removed)warnings.push(`Raw-peak-supported cycle metadata repair: control +${controlRepair.added}/-${controlRepair.removed}, sample +${sampleRepair.added}/-${sampleRepair.removed}.`);if(qc.r2<.995)warnings.push('Signal-guided paths failed internal trace reconstruction QC; editing estimates are withheld.');if(ambiguous.length)warnings.push('Some discovered paths have indistinguishable reference or event explanations. Inspect alternatives.');if(discoveryLimited)warnings.push('Path discovery reached its iteration limit.');if(geometry.condition===null||geometry.condition>100)warnings.push('Selected path contributions are poorly conditioned.');
 if(substitutionLinkage.length)warnings.push('Separate and linked substitution groups produce identical measured fluorescence. WT and substitution-group fractions are not uniquely identifiable; inspect annotation ranges.');
 const boundary=chosen.some(c=>c.model.events.some(e=>e.end-e.start>=maxDeletion||e.inserted.length>=maxInsertion));if(boundary)warnings.push('A supported event reaches the configured path-search boundary.');
 const localResidual=full.positionRmse,unknown=positions.map(p=>chosen.some(c=>c.model.sequence[p]==='N'||c.model.sequence[p-1]==='N')),localDiscordant=localResidual.filter((v,i)=>!unknown[i]&&v>(refs.some(r=>r.cuts.length>1)?Math.max(.003,2*(ctlScore.positionRmse[i]||0)+.002):Math.max(.06,(ctlScore.positionRmse[i]||0)+.04))).length;
 if(localDiscordant)warnings.push(`${localDiscordant} positions have unexplained local signal (possible substitution, artifact or missing edit path). Editing estimates are withheld.`);
 const unknownInsertion=chosen.some(c=>c.model.sequence.includes('N'));if(unknownInsertion)warnings.push('Inserted bases are unresolved; known-sequence inference cannot certify these target fractions.');
 const allowed=!unknownInsertion&&qc.r2>=.995&&!localDiscordant&&fit.converged&&!discoveryLimited&&geometry.condition!==null&&geometry.condition<=100;
 const perRef=refs.map((r,h)=>{const rows=chosen.filter(c=>c.model.reference===h&&!attributionAmbiguous.includes(c)),share=rows.reduce((sum,c)=>sum+c.fraction,0),estimable=allowed&&share>=.02&&!attributionAmbiguous.some(c=>c.model.aliases.some(a=>a.reference===h)),bins=new Map<number,number>();rows.forEach(c=>bins.set(net(c.model.events),(bins.get(net(c.model.events))||0)+c.fraction));const zero=rows.filter(c=>!c.model.events.length).reduce((sum,c)=>sum+c.fraction,0);return {name:r.name,share,wt:estimable?zero/share:null,indel:estimable?1-zero/share:null,inFrame:estimable?rows.filter(c=>c.model.events.length&&net(c.model.events)!==0&&net(c.model.events)%3===0).reduce((sum,c)=>sum+c.fraction,0)/share:null,outOfFrame:estimable?rows.filter(c=>net(c.model.events)%3!==0).reduce((sum,c)=>sum+c.fraction,0)/share:null,components:[...bins].map(([size,fraction])=>({size,fraction:share?fraction/share:0})),estimable,rows};});
 const targets=perRef.flatMap((a,h)=>(refs[h].input.cuts??[refs[h].input.cut]).map((cut,t)=>({reference:a.name,target:t+1,cut,edited:a.estimable?a.rows.filter(c=>inputEvents(h,c.model.events).some(e=>e.target===t+1||e.target===0)).reduce((sum,c)=>sum+c.fraction,0)/a.share:null})));
 const annotations=perRef.map((a,h)=>({sequenceGrouping:true,name:a.name,reference:refs[h].input.sequence,cut:refs[h].input.cut,cuts:refs[h].input.cuts,direction:refs[h].reverse?'reverse' as const:'forward' as const,share:a.share,estimable:a.estimable,components:a.rows.map(c=>({size:net(c.model.events),pathSizes:pathSizes(h,c.model.events),fraction:c.fraction/a.share,events:ambiguous.includes(c)?[]:inputEvents(h,c.model.events),unresolved:ambiguous.includes(c)}))}));
 const exactAlternatives=discoverExactAlternatives(refs,annotations,perRef.map(a=>a.rows),inputEvents,(sequence,p)=>emit(sequence,p,sc.lag));
 const bins=new Map<number,number>();chosen.forEach(c=>bins.set(net(c.model.events),(bins.get(net(c.model.events))||0)+c.fraction));const zero=chosen.filter(c=>!c.model.events.length).reduce((sum,c)=>sum+c.fraction,0);
 const first=Math.min(...refs.flatMap(r=>r.cuts)),plots=[{trace:input.control,trim:cr.trim},{trace:input.sample,trim:sr.trim}].map(({trace,trim})=>plotTrace(trace,first-trim,Math.max(...refs.flatMap(r=>r.cuts))-first+s.window/2));report.positions=positions;report.joint!.candidateCount=models.length;
 if(audit)audit({refs,positions,cols,models,ids,weights:fit.weights,theta,sampleLag:sc.lag,gains:sc.gains,observed:y,fitCache:{hits:fitCache.hits,misses:fitCache.misses,linearFits:fitCache.linearFits,nonlinearFits:fitCache.nonlinearFits,pairHits:fitCache.pairHits,pairMisses:fitCache.pairMisses},rawObserved:rawY,anchor,end,iterations,train,test,tune,controlPositionRmse:ctlScore.positionRmse});
 return {exactAlternatives,method:'pooled',referenceName:refs.length===1?refs[0].name:'Mixed amplicons',sampleName:input.sample.name,controlName:input.control.name,direction:'forward',cutPosition:first,alignmentIdentity:0,alignment:[],variants:[],plots,warnings,synthetic:!!input.synthetic,annotations,diagnostics:{substitutionLinkage:substitutionLinkage.map(a=>({...a,reference:refs[a.reference].name})),algorithmVersion:refs.some(r=>r.cuts.length>1)?'sanger-sequence-groups-oracle-prototype':'sanger-sequence-groups-oracle-prototype',signalModel:'calibrated',requestedMaxIndel:s.maxIndel,discordantPositions:localDiscordant,maxPositionRmse:Math.max(...localResidual),controlCalibration:cc,sampleCalibration:sc,contextEmission:context},mixture:{components:[...bins].map(([size,fraction])=>({size,fraction})),wt:zero,indel:1-zero,inFrame:chosen.filter(c=>c.model.events.length&&net(c.model.events)!==0&&net(c.model.events)%3===0).reduce((sum,c)=>sum+c.fraction,0),outOfFrame:chosen.filter(c=>net(c.model.events)%3!==0).reduce((sum,c)=>sum+c.fraction,0),r2:full.r2,rmse:full.rmse,upstreamError:1-sc.validationR2,bases:positions.length,fitStart:anchor,fitEnd:end,converged:fit.converged,passesChecks:allowed&&!attributionAmbiguous.length&&!substitutionLinkage.length},pooled:{preflight:report,controlR2:ctlScore.r2,sampleTrim:sr.trim,controlTrim:cr.trim,amplicons:perRef.map(({rows,estimable,...a})=>a)},joint:{groups,targets,reconstructionR2:qc.r2,heldOutR2:null,selectedModels:chosen.length,candidateModels:models.length,iterations,maxDeletion,ambiguousFraction:ambiguous.reduce((sum,c)=>sum+c.fraction,0),alignmentSearches:alignments}};
}
