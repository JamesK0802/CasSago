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
export function analyzeAlleleGroups(input:SignalPathInput,progress:(s:string)=>void=()=>{},audit?:(state:any)=>void):SangerAnalysisResult{
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
 const modelFamily=(m:Model)=>{const sizes=refs[m.reference].cuts.map(()=>0);for(const e of m.events){const distances=refs[m.reference].cuts.map(c=>Math.min(Math.abs(e.start-c),Math.abs(e.end-c))),t=distances.indexOf(Math.min(...distances));sizes[t]+=e.inserted.length-(e.end-e.start);}return m.reference+':'+sizes.join('/');};
 const seenEdits=refs.map(()=>new Set<string>()),signature=(sequence:string)=>sequence.slice(positions[0],positions.at(-1)!+1).padEnd(positions.length,'-');
 const signatures=new Map(models.map((m,id)=>[signature(m.sequence),id]));
 const refineSubstitutions=()=>{let changed=false;
 for(let sweep=0;sweep<2;sweep++)for(let h=0;h<refs.length;h++){
  const pure=(id:number)=>models[id].reference===h&&models[id].events.length>0&&models[id].events.every(e=>e.inserted.length===e.end-e.start);
  const spans=(refs[h].cuts.length>1?models.map((_,id)=>id):ids).filter(pure).map(id=>({start:Math.min(...models[id].events.map(e=>e.start)),end:Math.max(...models[id].events.map(e=>e.end))})).sort((a,b)=>a.start-b.start),blocks:{start:number;end:number}[]=[];
  for(const span of spans){const previous=blocks.at(-1);if(previous&&span.start<=previous.end){previous.end=Math.max(previous.end,span.end);}else blocks.push({...span});}
  for(const block of blocks){
   progress(`Checking substitution combinations · reference ${h+1}`);
   const inside=(id:number)=>pure(id)&&models[id].events.every(e=>e.start>=block.start-2&&e.end<=block.end+2),baseIds=ids.filter(id=>!inside(id)),baseFit=fitIds(baseIds),prediction=y.map((_,k)=>{let value=0,total=0;for(let j=0;j<baseIds.length;j++){const c=cols[baseIds[j]],at=k-k%4;value+=baseFit.weights[j]*c[k];total+=baseFit.weights[j]*(c[at]+c[at+1]+c[at+2]+c[at+3]);}return total?value/total:0;}),residual=y.map((v,k)=>v-prediction[k]);
   const allCandidates=models.flatMap((m,id)=>{if(!inside(id))return [];let upper=Infinity;for(let k=0;k<y.length;k++)if(cols[id][k]>.5)upper=Math.min(upper,y[k]/cols[id][k]);return upper>.002?[id]:[];});
   const candidates=refs[h].cuts.length>1?allCandidates.map(id=>{let rhs=0,norm=0;for(let k=0;k<y.length;k++){const d=cols[id][k]-cols[h][k];rhs+=d*residual[k];norm+=d*d;}return {id,gain:rhs>0&&norm?rhs*rhs/norm:0};}).filter(c=>c.gain>0).sort((a,b)=>b.gain-a.gain).slice(0,160).map(c=>c.id):allCandidates;
   if(!candidates.length)continue;
   const proposals=substitutionProposals(candidates.map(id=>cols[id]),cols[h],residual,baseFit.weights[baseIds.indexOf(h)]||0,12);
   const value=(candidateIds:number[],weights:number[])=>score(candidateIds,weights,tune).rmse**2+Math.max(1e-10,ctlScore.rmse**2*Math.log(y.length)/y.length)*candidateIds.filter((id,j)=>models[id].events.length&&weights[j]>.001).length;
   let selected={ids,fit,value:value(ids,fit.weights)};const initialValue=selected.value;
   for(const proposal of proposals){const trial=[...baseIds,...proposal.map(i=>candidates[i])],f=fitIds(trial),v=value(trial,f.weights);if(v<selected.value-1e-12)selected={ids:trial,fit:f,value:v};}
   if(selected.value<initialValue-1e-12)changed=true;ids=selected.ids;fit=selected.fit;
  }
 }

 return changed;};
  const addVariant=(h:number,events:PathEdit[])=>{events=normalizeMoleculeEdits(refs[h].sequence,events);const editKey=JSON.stringify(events);if(seenEdits[h].has(editKey))return;seenEdits[h].add(editKey);if(!eventsRespectIndelOrigin(events,refs[h].cuts,s.editStartRadius))return;const sequence=materialize(refs[h].sequence,events);if(!sequence)return;const key=signature(sequence),old=signatures.get(key),alias={reference:h,events};if(old!==undefined){if(!models[old].aliases.some(a=>a.reference===h&&JSON.stringify(a.events)===JSON.stringify(events)))models[old].aliases.push(alias);return;}const id=models.length;models.push({reference:h,sequence,events,aliases:[alias]});cols.push(positions.flatMap(p=>emit(sequence,p,sc.lag)));signatures.set(key,id);};
 const addResidualVariants=(pred:number[][])=>{
  for(let h=0;h<refs.length;h++)for(const cut of refs[h].cuts){
   for(const alpha of [.025,.05,.1,.2,.4].map(f=>f*cc.weights[h])){
    const call=(p:number,column:number[])=>{const k=(positionIndex[p]??-1);if(k<0)return '';const values=TRACE_BASES.map((_,b)=>sample[p-sr.trim][b]-pred[k][b]+alpha*column[k*4+b]);return TRACE_BASES[values.indexOf(Math.max(...values))];};
    for(let start=Math.max(anchor,cut-(s.editStartRadius??10));start<=cut+(s.editStartRadius??10);start++){
     let inserted='';for(let size=1;size<=maxInsertion;size++){inserted+=call(start+size-1,cols[h]);if(inserted.length===size)addVariant(h,[{start,end:start,inserted}]);}
     let replacement='';for(let size=1;size<=Math.min(10,refs[h].sequence.length-start);size++){replacement+=call(start+size-1,cols[h]);if(replacement.length===size&&replacement!==refs[h].sequence.slice(start,start+size))addVariant(h,[{start,end:start+size,inserted:replacement}]);}
    }
   }
  }
 };
 const discover=(rounds:number,phase=1)=>{let complete=false;for(let round=0;round<rounds;round++){
  progress(`${phase===1?'Fitting sequence groups':'Refining redistributed sequence groups'} · round ${round+1}`);
  const pred=positions.map((_,p)=>{const v=TRACE_BASES.map((_,b)=>ids.reduce((sum,id,h)=>sum+fit.weights[h]*cols[id][4*p+b],0)),total=v.reduce((a,b)=>a+b,0);return v.map(x=>total?x/total:0);});

  // Search physical single-event molecules. Bases come from measured residuals,
  // never the synthetic recipe. Each candidate is refitted as a whole molecule.
  if(round===0)for(let h=0;h<refs.length;h++)for(const cut of refs[h].cuts)for(let start=Math.max(anchor,cut-(s.editStartRadius??10));start<=cut+(s.editStartRadius??10);start++)for(let size=1;size<=maxDeletion;size++)addVariant(h,[{start,end:start+size,inserted:''}]);
  addResidualVariants(pred);
  for(const item of ids.map((id,k)=>({id,weight:fit.weights[k]})).filter(x=>x.weight>.001)){
   const model=models[item.id];if(!model.events.some(e=>e.inserted.length))continue;
   const e=model.events[0];if(model.events.length!==1)continue;
   for(const cut of refs[model.reference].cuts)for(let start=Math.max(anchor,cut-(s.editStartRadius??10));start<=cut+(s.editStartRadius??10);start++){
    let inserted='';for(let j=0;j<e.inserted.length;j++){const p=start+j,k=(positionIndex[p]??-1);if(k<0)break;const values=TRACE_BASES.map((_,b)=>sample[p-sr.trim][b]-pred[k][b]+item.weight*cols[item.id][k*4+b]);inserted+=TRACE_BASES[values.indexOf(Math.max(...values))];}
    if(inserted.length===e.inserted.length)addVariant(model.reference,[{start,end:start+e.end-e.start,inserted}]);
   }
   for(let j=0;j<e.inserted.length;j++)for(const base of TRACE_BASES)if(base!==e.inserted[j])addVariant(model.reference,[{...e,inserted:e.inserted.slice(0,j)+base+e.inserted.slice(j+1)}]);
  }
  if(multiTarget){
   const before=models.length;
   for(let h=0;h<refs.length;h++){
    const singles=ids.map((id,j)=>({id,weight:fit.weights[j],model:models[id]})).filter(c=>c.model.reference===h&&c.model.events.length===1&&c.weight>.0005),forTarget=(target:number)=>singles.filter(c=>{const e=c.model.events[0],d=refs[h].cuts.map(cut=>Math.min(Math.abs(e.start-cut),Math.abs(e.end-cut)));return d.indexOf(Math.min(...d))===target;}).sort((a,b)=>b.weight-a.weight).slice(0,8);
    for(const left of forTarget(0))for(const right of forTarget(1))addVariant(h,[...left.model.events,...right.model.events].sort((a,b)=>a.start-b.start));
   }
   if(models.length>before){const trial=[...ids,...Array.from({length:models.length-before},(_,j)=>before+j)],f=fitIds(trial),value=objective(trial,f.weights);if(value<best-1e-11){ids=trial.filter((id,j)=>id<refs.length||f.weights[j]>=.0005);fit=fitIds(ids);best=objective(ids,fit.weights);iterations++;continue;}}
  }
  // Existing rejected candidates remain eligible; no truth-derived edit list is used.
  screen.update(models.map(modelFamily));const additionScores=multiTarget?screen.correlate(y.map((v,k)=>v-pred[k>>2][k%4])):screen.additionScores(y,pred.flat());
  const available=models.map((_,id)=>id).filter(id=>!ids.includes(id)).sort((a,b)=>additionScores[b]-additionScores[a]).slice(0,20);let next:{ids:number[];fit:ReturnType<typeof fitEmissionMixture>;rmse:number}|undefined;
  for(const id of available){const trial=[...ids,id],f=fitIds(trial);if(f.weights.at(-1)!<(multiTarget ? .0005 : .005))continue;const rmse=objective(trial,f.weights);if(rmse<best-1e-10&&(!next||rmse<next.rmse))next={ids:trial,fit:f,rmse};}
  for(const old of ids){if(old<refs.length)continue;const weight=fit.weights[ids.indexOf(old)],netSize=(m:Model)=>m.events.reduce((sum,e)=>sum+e.inserted.length-(e.end-e.start),0);
   const replacements=models.map((m,id)=>({m,id})).filter(x=>!ids.includes(x.id)&&x.m.reference===models[old].reference&&(multiTarget?screen.sameFamily(x.id,old):netSize(x.m)===netSize(models[old]))).map(({id})=>{if(multiTarget)return {id,gain:screen.replacementGain(id,old,weight)};let gain=0;for(let k=0;k<y.length;k++){const delta=cols[id][k]-cols[old][k];gain+=2*weight*delta*(y[k]-pred[k>>2][k%4])-weight*weight*delta*delta;}return {id,gain};}).sort((a,b)=>b.gain-a.gain).slice(0,5);
   for(const {id}of replacements){const trial=ids.map(x=>x===old?id:x),f=fitIds(trial),value=objective(trial,f.weights);if(value<best-1e-11&&(!next||value<next.rmse))next={ids:trial,fit:f,rmse:value};}
  }
  if(!next){
   if(refs.some(r=>r.cuts.length>1)&&substitutionRounds++<2&&refineSubstitutions()){best=objective(ids,fit.weights);continue;}
   if(refs.some(r=>r.cuts.length>1)&&insertionPolishes++<2){
    const before=models.length;
    // Retain competing inserted-base strings rather than committing to one
    // argmax base at each cycle of a contaminated group residual.
    for(const [j,id]of ids.entries()){
     const model=models[id],e=model.events[0],weight=fit.weights[j];if(model.events.length!==1||!e||e.end!==e.start||!e.inserted.length||weight<.003)continue;
     const proposals:{event:PathEdit;gain:number}[]=[];
     for(const cut of refs[model.reference].cuts)for(let start=Math.max(anchor,cut-(s.editStartRadius??10));start<=cut+(s.editStartRadius??10);start++){
      let beam=[{inserted:'',score:0}];
      for(let b=0;b<e.inserted.length;b++){
       const k=(positionIndex[start+b]??-1);if(k<0){beam=[];break;}
       const values=TRACE_BASES.map((_,c)=>(y[k*4+c]-pred[k][c]+weight*cols[id][k*4+c])/weight),maximum=Math.max(...values);
       beam=beam.flatMap(state=>TRACE_BASES.flatMap((base,c)=>values[c]>=maximum-.7?[{inserted:state.inserted+base,score:state.score+values[c]}]:[])).sort((a,b)=>b.score-a.score).slice(0,32);
      }
      for(const state of beam){const event={start,end:start,inserted:state.inserted},sequence=materialize(refs[model.reference].sequence,[event]);if(!sequence)continue;let gain=0;for(let k=0;k<positions.length;k++){const p=positions[k];if(sequence[p]===model.sequence[p]&&sequence[p-1]===model.sequence[p-1]&&sequence[p+1]===model.sequence[p+1])continue;const column=emit(sequence,p,sc.lag);for(let b=0;b<4;b++){const d=column[b]-cols[id][k*4+b];gain+=2*weight*d*(y[k*4+b]-pred[k][b])-weight*weight*d*d;}}if(gain>1e-10)proposals.push({event,gain});}
     }
     for(const proposal of proposals.sort((a,b)=>b.gain-a.gain).slice(0,16))addVariant(model.reference,[proposal.event]);
    }
    if(models.length>before)continue;
   }
   complete=true;break;
  }ids=next.ids;fit=next.fit;best=next.rmse;iterations++;
 }
 discoveryLimited ||= !complete;};
 discover(80);
 // Replacement moves let a simpler scaffold compete directly with a selected
 // complex path. Adding it as a correlated extra column can otherwise leave
 // the unpenalized refit trapped on the more complex explanation.
 for(let sweep=0;sweep<3;sweep++){let changed=false;for(const id of [...ids]){const model=models[id];if(!model.events.length)continue;for(let event=0;event<model.events.length;event++){
  const events=normalizeMoleculeEdits(refs[model.reference].sequence,model.events.filter((_,k)=>k!==event)),sequence=materialize(refs[model.reference].sequence,events);if(!sequence)continue;const key=signature(sequence);let replacement=signatures.get(key);
  if(replacement===undefined){replacement=models.length;models.push({reference:model.reference,sequence,events,aliases:[{reference:model.reference,events}]});cols.push(positions.flatMap(p=>emit(sequence,p,sc.lag)));signatures.set(key,replacement);}
  const trial=[...new Set(ids.map(x=>x===id?replacement!:x))],f=fitIds(trial),value=objective(trial,f.weights);
  if(value<best-1e-8){ids=trial;fit=f;best=value;changed=true;break;}
 }}if(!changed)break;}
 for(let sweep=0;sweep<3;sweep++){let removed=false;for(const id of [...ids]){if(id<refs.length)continue;const trial=ids.filter(x=>x!==id),f=fitIds(trial),value=objective(trial,f.weights);if(value<=best+1e-11){ids=trial;fit=f;best=value;removed=true;}}if(!removed)break;}
 // Revisit substitution blocks jointly after indel positions and sequences
 // stabilize. Single-column greedy updates can be trapped on several wrong
 // substitution molecules even when the correct joint candidates exist.
 refineSubstitutions();
 const polishDictionary=()=>{
  progress('Jointly refitting competing sequence groups');
  const families=models.map(modelFamily);
  const polished=polishSignalDictionary(cols,y,families,ids,fit.weights),f=fitIds(polished.ids),value=objective(polished.ids,f.weights);if(value<objective(ids,fit.weights)-1e-12){ids=polished.ids;fit=f;best=value;}
 }
 if(multiTarget){polishDictionary();progress('Rechecking sequence candidates after joint redistribution');const pred=positions.map((_,p)=>{const v=TRACE_BASES.map((_,b)=>ids.reduce((sum,id,h)=>sum+fit.weights[h]*cols[id][4*p+b],0)),total=v.reduce((a,b)=>a+b,0);return v.map(x=>total?x/total:0);});addResidualVariants(pred);const residual=y.map((v,k)=>v-pred[k>>2][k%4]);for(const [j,id]of ids.entries())for(const events of compoundBaseProposals(models[id].events,cols[id],residual,fit.weights[j],positionIndex))addVariant(models[id].reference,events);polishDictionary();}

 ids=ids.filter((id,i)=>id<refs.length||fit.weights[i]>=(multiTarget ? .0005 : .005));fit=fitIds(ids);
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
 return {exactAlternatives,method:'pooled',referenceName:refs.length===1?refs[0].name:'Mixed amplicons',sampleName:input.sample.name,controlName:input.control.name,direction:'forward',cutPosition:first,alignmentIdentity:0,alignment:[],variants:[],plots,warnings,synthetic:!!input.synthetic,annotations,diagnostics:{substitutionLinkage:substitutionLinkage.map(a=>({...a,reference:refs[a.reference].name})),algorithmVersion:refs.some(r=>r.cuts.length>1)?'sanger-sequence-groups-v4':'sanger-sequence-groups-v2',signalModel:'calibrated',requestedMaxIndel:s.maxIndel,discordantPositions:localDiscordant,maxPositionRmse:Math.max(...localResidual),controlCalibration:cc,sampleCalibration:sc,contextEmission:context},mixture:{components:[...bins].map(([size,fraction])=>({size,fraction})),wt:zero,indel:1-zero,inFrame:chosen.filter(c=>c.model.events.length&&net(c.model.events)!==0&&net(c.model.events)%3===0).reduce((sum,c)=>sum+c.fraction,0),outOfFrame:chosen.filter(c=>net(c.model.events)%3!==0).reduce((sum,c)=>sum+c.fraction,0),r2:full.r2,rmse:full.rmse,upstreamError:1-sc.validationR2,bases:positions.length,fitStart:anchor,fitEnd:end,converged:fit.converged,passesChecks:allowed&&!attributionAmbiguous.length&&!substitutionLinkage.length},pooled:{preflight:report,controlR2:ctlScore.r2,sampleTrim:sr.trim,controlTrim:cr.trim,amplicons:perRef.map(({rows,estimable,...a})=>a)},joint:{groups,targets,reconstructionR2:qc.r2,heldOutR2:null,selectedModels:chosen.length,candidateModels:models.length,iterations,maxDeletion,ambiguousFraction:ambiguous.reduce((sum,c)=>sum+c.fraction,0),alignmentSearches:alignments}};
}
