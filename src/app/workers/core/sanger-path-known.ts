import {normalizeMoleculeEdits} from './sanger-edit-normalize';
import {validateSangerProductBounds} from '../../utils/sanger-product-bounds';
import {repairSignalTraceCycles} from './sanger-cycle';
import {PooledAmplicon,PooledPreflight,SangerAnalysisSettings,SangerAnalysisResult,SangerJointEvent,SangerTrace,TRACE_BASES} from '../../models/sanger-analysis';
import {exactBindingSites} from './sanger-primer-designer';
import {reverseComplement} from '../../utils/biology.utils';
import {registerTrace} from './sanger-pooled';
import {contrastPeaks,integratedPeaks,calibrateSignal,correctGain,gramCondition,rawSignalScore} from './sanger-signal';
import {emissionVector,fitContextEmission,fitNormalizedEmissionMixture as fitEmissionMixtureOriginal} from './sanger-emission';
import {plotTrace} from './sanger-analysis';
const DNA='ACGT';
const fitEmissionMixture=(columns:number[][],observed:number[])=>fitEmissionMixtureOriginal(columns,observed);
import {SignalPathInput,PathEdit,Model,materialize,prepare,assertSignalTraceGrid,eventsRespectIndelOrigin,discoverSignalPath} from './sanger-path-common';
export function analyzeSignalPaths(input:SignalPathInput,progress:(s:string)=>void=()=>{},audit?:(state:any)=>void):SangerAnalysisResult{
 if(!input.ampliconConfirmed||!input.control||!input.sample)throw Error('Confirm submitted products and provide matched control/sample chromatograms.');
 const controlRepair=repairSignalTraceCycles(input.control),sampleRepair=repairSignalTraceCycles(input.sample);input={...input,control:controlRepair.trace,sample:sampleRepair.trace};
 assertSignalTraceGrid(input.control);assertSignalTraceGrid(input.sample);
 const {refs,report,anchor,maxDeletion,maxInsertion}=prepare(input.amplicons,input.primer,input.settings),s=input.settings;
 progress('Registering traces to the sequencing primer');const registered=refs.map(r=>({name:r.name,sequence:r.sequence,cut:anchor}));
 const cr=registerTrace(input.control,registered,s,integratedPeaks(input.control,[0,0,0,0]).vectors),sr=registerTrace(input.sample,registered,s,integratedPeaks(input.sample,[0,0,0,0]).vectors),ctl=contrastPeaks(input.control).vectors,raw=contrastPeaks(input.sample).vectors;
 const calibrate=(vectors:number[][],trim:number)=>{const obs=Array.from({length:vectors.length+trim},(_,p)=>vectors[p-trim]||[0,0,0,0]),rows=Array.from({length:anchor-12},(_,p)=>p).filter(p=>p>=Math.max(s.readStart,trim+10)&&p-trim<vectors.length-2);return calibrateSignal(obs,rows,refs.map(r=>Array.from(r.sequence,(_,p)=>TRACE_BASES.map(b=>Number(r.sequence[p]===b)))));};
 const cc=calibrate(ctl,cr.trim),sc=calibrate(raw,sr.trim);if(cc.validationR2<.8||sc.validationR2<.8)throw Error('Signal-path upstream calibration failed.');
 const end=Math.min(Math.max(...report.positions)+1,cr.trim+ctl.length-3,sr.trim+raw.length-3),positions=report.positions.filter(p=>p<end);if(positions.length<120)throw Error('Insufficient actual trace coverage for path inference.');
 const wt=ctl.map(v=>correctGain(v,cc.gains)),sample=raw.map(v=>correctGain(v,sc.gains)),controlObs=Array.from({length:wt.length+cr.trim},(_,p)=>wt[p-cr.trim]||[0,0,0,0]);
 const context=fitContextEmission(refs.map(r=>r.sequence),positions,controlObs,cc.weights,cc.lag),theta=context.theta;
 const emit=(sequence:string,p:number,lag:number)=>{const v=emissionVector(sequence,p,lag,theta);if(sequence[p]==='N')for(let b=0;b<4;b++)v[b]+=(1-Math.abs(lag))/4;return v;};
 const models:Model[]=refs.map((r,h)=>({reference:h,sequence:r.sequence,events:[],aliases:[{reference:h,events:[]}]})),cols=models.map(m=>positions.flatMap(p=>emit(m.sequence,p,sc.lag))),controlCols=models.map(m=>positions.flatMap(p=>emit(m.sequence,p,cc.lag)));
 const ctlFit=fitEmissionMixture(controlCols,positions.flatMap(p=>wt[p-cr.trim])),ctlScore=rawSignalScore(controlCols,positions.flatMap(p=>ctl[p-cr.trim]),ctlFit.weights,cc.gains);if(ctlScore.r2<.85)throw Error('WT product references do not explain the unedited control.');
 if(s.jointEngine==='signal-hybrid')refs.forEach((r,h)=>r.cuts.forEach(cut=>{for(let size=-Math.min(s.maxIndel,maxDeletion);size<=Math.min(s.maxIndel,maxInsertion);size++){if(!size)continue;const events=[{start:cut,end:cut+Math.max(0,-size),inserted:'N'.repeat(Math.max(0,size))}],sequence=materialize(r.sequence,events);if(sequence){models.push({reference:h,sequence,events,aliases:[{reference:h,events}]});cols.push(positions.flatMap(p=>emit(sequence,p,sc.lag)));}}}));
 if(s.jointEngine==='signal-hybrid')refs.forEach((r,h)=>r.cuts.forEach(cut=>{for(const shift of [-1,0,1])for(const base of TRACE_BASES){const at=cut+shift,b=TRACE_BASES.indexOf(base);if((sample[at-sr.trim]?.[b]||0)<.02||!eventsRespectIndelOrigin([{start:at,end:at,inserted:base}],r.cuts,s.editStartRadius))continue;const events=[{start:at,end:at,inserted:base}],sequence=materialize(r.sequence,events);if(sequence){models.push({reference:h,sequence,events,aliases:[{reference:h,events}]});cols.push(positions.flatMap(p=>emit(sequence,p,sc.lag)));}}}));
 const y=positions.flatMap(p=>sample[p-sr.trim]),train=positions.flatMap((_,i)=>i%5<3?[i]:[]),tune=positions.flatMap((_,i)=>i%5===3?[i]:[]),test=positions.flatMap((_,i)=>i%5===4?[i]:[]),subset=(v:number[],rows:number[])=>rows.flatMap(i=>v.slice(4*i,4*i+4)),trainY=subset(y,train),tuneY=subset(y,tune);
 const trained=Array(end+1).fill(false);train.forEach(i=>trained[positions[i]]=true);
 const fitIds=(ids:number[])=>fitEmissionMixture(ids.map(id=>subset(cols[id],train)),trainY);
 const score=(ids:number[],weights:number[],rows:number[])=>rawSignalScore(ids.map(id=>subset(cols[id],rows)),rows.flatMap(i=>raw[positions[i]-sr.trim]),weights,sc.gains);
 const objective=(ids:number[],weights:number[])=>score(ids,weights,tune).rmse**2+.0001*ids.reduce((sum,id,h)=>sum+weights[h]*models[id].events.reduce((a,e)=>a+1+.05*e.inserted.length,0),0);
 let ids=refs.map((_,i)=>i),fit=fitIds(ids),best=objective(ids,fit.weights),iterations=0,alignments=0;
 const signatures=new Map(models.map((m,id)=>[positions.map(p=>m.sequence[p]||'-').join(''),id]));
 for(let round=0;round<20;round++){
  progress(`Discovering ${s.jointEngine==='signal-local'?'local':'alignment'} paths · round ${round+1}`);
  const pred=positions.map((_,p)=>{const v=TRACE_BASES.map((_,b)=>ids.reduce((sum,id,h)=>sum+fit.weights[h]*cols[id][4*p+b],0)),total=v.reduce((a,b)=>a+b,0);return v.map(x=>total?x/total:0);});
  const candidates:number[]=[];
  for(let h=0;h<refs.length;h++)for(const seed of [
    ...[.01,.06,.18].map(alpha=>({alpha,column:cols[h],normalize:false})),
    ...ids.map((id,k)=>({id,weight:fit.weights[k]})).filter(x=>models[x.id].reference===h&&x.weight>.02).sort((a,b)=>b.weight-a.weight).slice(0,3).map(x=>({alpha:x.weight,column:cols[x.id],normalize:true}))
   ]){const alpha=seed.alpha;
   const profile=Array.from({length:end+1},()=>[0,0,0,0]);positions.forEach((p,k)=>{if(trained[p]){for(let b=0;b<4;b++)profile[p][b]=sample[p-sr.trim][b]-pred[k][b]+alpha*seed.column[k*4+b];if(seed.normalize){const total=profile[p].reduce((sum,v)=>sum+Math.max(0,v),0);if(total>.001)profile[p]=profile[p].map(v=>Math.max(0,v)/total);}}});
   const path=discoverSignalPath(refs[h].sequence,profile,anchor,end,refs[h].cuts,{maxDeletion,maxInsertion,alpha,local:s.jointEngine==='signal-local',training:trained,startRadius:s.editStartRadius});alignments++;
   if(!path||!path.events.length)continue;
   const variants:PathEdit[][]=[path.events,...path.events.map(e=>[e]),...path.events.map((_,i)=>path.events.filter((_,j)=>j!==i))];
   // Reduce unsupported local changes and refine insertion positions/bases using
   // training-signal nucleotides. No exponential enumeration of inserted strings.
   path.events.forEach((e,k)=>{if(e.inserted.length!==1||e.end!==e.start)return;for(let offset=-3;offset<=3;offset++){
    const start=e.start+offset,p=start+path.events.filter(x=>x.end<=start&&x!==e).reduce((sum,x)=>sum+x.inserted.length-(x.end-x.start),0);
    if(start<anchor||p<anchor||p>=end)continue;const bases=trained[p]?DNA.split('').filter((b,i)=>profile[p][i]>.005):['N'];
    for(const base of bases)variants.push(path.events.map((x,i)=>i===k?{start,end:start,inserted:base}:x));
   }});
   for(const rawEvents of variants){const normalized=normalizeMoleculeEdits(refs[h].sequence,rawEvents);if(!eventsRespectIndelOrigin(normalized,refs[h].cuts,s.editStartRadius))continue;let shift=0;const events=normalized.filter(e=>{const observedStart=e.start+shift;shift+=e.inserted.length-(e.end-e.start);return observedStart<end;});const sequence=materialize(refs[h].sequence,events);if(!sequence)continue;const signature=positions.map(p=>sequence[p]||'-').join(''),alias={reference:h,events},old=signatures.get(signature);
    if(old!==undefined){if(!models[old].aliases.some(a=>a.reference===h&&JSON.stringify(a.events)===JSON.stringify(events)))models[old].aliases.push(alias);continue;}
    const id=models.length;models.push({reference:h,sequence,events,aliases:[alias]});cols.push(positions.flatMap(p=>emit(sequence,p,sc.lag)));signatures.set(signature,id);candidates.push(id);
   }
  }
  // Insertion bases have no reference at held-out cycles. Refine the discovered
  // gap scaffold using all observed channels, while keeping coefficient splits.
  // The resulting reconstruction score is an internal check, not independent CV.
  const addVariant=(h:number,events:PathEdit[])=>{events=normalizeMoleculeEdits(refs[h].sequence,events);if(!eventsRespectIndelOrigin(events,refs[h].cuts,s.editStartRadius))return;const sequence=materialize(refs[h].sequence,events);if(!sequence)return;const signature=positions.map(p=>sequence[p]||'-').join(''),old=signatures.get(signature),alias={reference:h,events};if(old!==undefined){if(!models[old].aliases.some(a=>a.reference===h&&JSON.stringify(a.events)===JSON.stringify(events)))models[old].aliases.push(alias);return;}const id=models.length;models.push({reference:h,sequence,events,aliases:[alias]});cols.push(positions.flatMap(p=>emit(sequence,p,sc.lag)));signatures.set(signature,id);};
  for(const item of ids.map((id,k)=>({id,weight:fit.weights[k]})).filter(x=>x.weight>.01)){
   const model=models[item.id];if(!model.events.some(e=>e.inserted.length))continue;
   let shift=0;const events=model.events.map(e=>{let inserted='';for(let j=0;j<e.inserted.length;j++){const p=e.start+shift+j,k=positions.indexOf(p);if(k<0){inserted+='N';continue;}const values=TRACE_BASES.map((_,b)=>sample[p-sr.trim][b]-pred[k][b]+item.weight*cols[item.id][k*4+b]),best=values.indexOf(Math.max(...values));inserted+=values[best]>.005?TRACE_BASES[best]:'N';}shift+=e.inserted.length-(e.end-e.start);return {...e,inserted};});addVariant(model.reference,events);
  }
  // Screen bounded insertion lengths by downstream residual support; construct
  // only the top supported lengths, with inserted bases from observed signal.
  if(round<3)for(let h=0;h<refs.length;h++)for(const cut of refs[h].cuts){
   const ranked=Array.from({length:maxInsertion},(_,i)=>i+1).map(size=>{let support=0,count=0;for(const k of train){const p=positions[k];if(p<cut+size+2||p>cut+size+100)continue;const v=emit(refs[h].sequence,p-size,sc.lag);support+=v.reduce((sum,x,b)=>sum+x*(y[k*4+b]-pred[k][b]),0);count++;}return {size,support:count?support/count:-Infinity};}).sort((a,b)=>b.support-a.support).slice(0,3);
   for(const {size,support}of ranked){if(support<.005)continue;let inserted='';for(let j=0;j<size;j++){const p=cut+j,k=positions.indexOf(p);if(k<0){inserted+='N';continue;}const values=TRACE_BASES.map((_,b)=>sample[p-sr.trim][b]-pred[k][b]+.08*cols[h][k*4+b]),best=values.indexOf(Math.max(...values));inserted+=values[best]>.005?TRACE_BASES[best]:'N';}addVariant(h,[{start:cut,end:cut,inserted}]);}
  }
  // Existing rejected candidates remain eligible; no truth-derived edit list is used.
  const available=models.map((_,id)=>id).filter(id=>!ids.includes(id)).map(id=>({id,score:train.reduce((sum,i)=>sum+TRACE_BASES.reduce((v,_,b)=>v+cols[id][i*4+b]*(y[i*4+b]-pred[i][b]),0),0)})).sort((a,b)=>b.score-a.score).slice(0,12).map(x=>x.id);let next:{ids:number[];fit:ReturnType<typeof fitEmissionMixture>;rmse:number}|undefined;
  for(const id of available){const trial=[...ids,id],f=fitIds(trial);if(f.weights.at(-1)!<.005)continue;const rmse=objective(trial,f.weights);if(rmse<best-Math.max(.0000001,best*.015)&&(!next||rmse<next.rmse))next={ids:trial,fit:f,rmse};}
  if(!next)break;ids=next.ids;fit=next.fit;best=next.rmse;iterations++;
 }
 // Replacement moves let a simpler scaffold compete directly with a selected
 // complex path. Adding it as a correlated extra column can otherwise leave
 // the unpenalized refit trapped on the more complex explanation.
 for(let sweep=0;sweep<3;sweep++){let changed=false;for(const id of [...ids]){const model=models[id];if(!model.events.length)continue;for(let event=0;event<model.events.length;event++){
  const events=normalizeMoleculeEdits(refs[model.reference].sequence,model.events.filter((_,k)=>k!==event)),sequence=materialize(refs[model.reference].sequence,events);if(!sequence)continue;const signature=positions.map(p=>sequence[p]||'-').join('');let replacement=signatures.get(signature);
  if(replacement===undefined){replacement=models.length;models.push({reference:model.reference,sequence,events,aliases:[{reference:model.reference,events}]});cols.push(positions.flatMap(p=>emit(sequence,p,sc.lag)));signatures.set(signature,replacement);}
  const trial=[...new Set(ids.map(x=>x===id?replacement!:x))],f=fitIds(trial),value=objective(trial,f.weights);
  if(value<best-1e-8){ids=trial;fit=f;best=value;changed=true;break;}
 }}if(!changed)break;}
 for(let sweep=0;sweep<3;sweep++){let removed=false;for(const id of [...ids]){if(id<refs.length)continue;const trial=ids.filter(x=>x!==id),f=fitIds(trial),value=objective(trial,f.weights);if(value<=best+1e-7){ids=trial;fit=f;best=value;removed=true;}}if(!removed)break;}
 ids=ids.filter((id,i)=>id<refs.length||fit.weights[i]>=.005);fit=fitIds(ids);
 const qc=score(ids,fit.weights,test),full=score(ids,fit.weights,positions.map((_,i)=>i)),geometry=gramCondition(ids.map(id=>cols[id]));
 const net=(events:PathEdit[])=>events.reduce((sum,e)=>sum+e.inserted.length-(e.end-e.start),0),inputEvents=(h:number,events:PathEdit[]):SangerJointEvent[]=>{const r=refs[h];return events.map(e=>{const distances=r.cuts.map(c=>Math.min(Math.abs(e.start-c),Math.abs(e.end-c))),near=distances.indexOf(Math.min(...distances)),span=r.cuts.length===2&&e.start<=r.cuts[0]&&e.end>=r.cuts[1],target=span?0:Math.min(...distances)>60?-1:(r.reverse?r.cuts.length-near:near+1),a=r.offset+e.start,b=r.offset+e.end;return {target,start:r.reverse?r.input.sequence.length-b:a,end:r.reverse?r.input.sequence.length-a:b,inserted:e.inserted.length,insertedSequence:r.reverse?reverseComplement(e.inserted):e.inserted};}).sort((a,b)=>a.start-b.start);};
 const chosen=ids.flatMap((id,i)=>fit.weights[i]>1e-9?[{model:models[id],fraction:fit.weights[i]}]:[]),aliasKey=(a:{reference:number;events:PathEdit[]})=>a.reference+':'+JSON.stringify(inputEvents(a.reference,a.events));
 const ambiguous=chosen.filter(c=>new Set(c.model.aliases.map(aliasKey)).size>1),refAmbiguous=ambiguous.filter(c=>new Set(c.model.aliases.map(a=>a.reference)).size>1);
 const semantic=(a:{reference:number;events:PathEdit[]})=>a.reference+':'+(a.events.length?'edit':'WT')+':'+net(a.events)+':'+[...new Set(inputEvents(a.reference,a.events).map(e=>e.target))].sort().join('/');
 const attributionAmbiguous=ambiguous.filter(c=>new Set(c.model.aliases.map(semantic)).size>1);
 const groups=chosen.map((c,i)=>({id:`P${i+1}`,reference:refAmbiguous.includes(c)?null:refs[c.model.reference].name,fraction:c.fraction,netSize:new Set(c.model.aliases.map(a=>net(a.events))).size>1?null:net(c.model.events),events:ambiguous.includes(c)?[]:inputEvents(c.model.reference,c.model.events),alternatives:ambiguous.includes(c)?c.model.aliases.map(a=>({reference:refs[a.reference].name,events:inputEvents(a.reference,a.events)})):[]}));
 const warnings:string[]=['Experimental signal-guided inference: a good trace fit does not establish unique allele identities.'];if(s.editStartRadius!==undefined)warnings.push(`Indel initiation constrained to cut ±${s.editStartRadius} bp; off-target-window edits are outside this model.`);if(controlRepair.added+controlRepair.removed+sampleRepair.added+sampleRepair.removed)warnings.push(`Raw-peak-supported cycle metadata repair: control +${controlRepair.added}/-${controlRepair.removed}, sample +${sampleRepair.added}/-${sampleRepair.removed}.`);if(qc.r2<.995)warnings.push('Signal-guided paths failed internal trace reconstruction QC; editing estimates are withheld.');if(ambiguous.length)warnings.push('Some discovered paths have indistinguishable reference or event explanations. Inspect alternatives.');if(iterations===20)warnings.push('Path discovery reached its iteration limit.');if(geometry.condition===null||geometry.condition>100)warnings.push('Selected path contributions are poorly conditioned.');
 const boundary=chosen.some(c=>c.model.events.some(e=>e.end-e.start>=maxDeletion||e.inserted.length>=maxInsertion));if(boundary)warnings.push('A supported event reaches the configured path-search boundary.');
 const localResidual=full.positionRmse,unknown=positions.map(p=>chosen.some(c=>c.model.sequence[p]==='N'||c.model.sequence[p-1]==='N')),localDiscordant=localResidual.filter((v,i)=>!unknown[i]&&v>Math.max(.06,(ctlScore.positionRmse[i]||0)+.04)).length;
 if(localDiscordant)warnings.push(`${localDiscordant} positions have unexplained local signal (possible substitution, artifact or missing edit path). Editing estimates are withheld.`);
 const unknownInsertion=chosen.some(c=>c.model.sequence.includes('N'));if(unknownInsertion)warnings.push('Inserted bases are unresolved; known-sequence inference cannot certify these target fractions.');
 const allowed=!unknownInsertion&&qc.r2>=.995&&!localDiscordant&&fit.converged&&iterations<20&&!boundary&&geometry.condition!==null&&geometry.condition<=100;
 const perRef=refs.map((r,h)=>{const rows=chosen.filter(c=>c.model.reference===h&&!attributionAmbiguous.includes(c)),share=rows.reduce((sum,c)=>sum+c.fraction,0),estimable=allowed&&share>=.02&&!attributionAmbiguous.some(c=>c.model.aliases.some(a=>a.reference===h)),bins=new Map<number,number>();rows.forEach(c=>bins.set(net(c.model.events),(bins.get(net(c.model.events))||0)+c.fraction));const zero=rows.filter(c=>!c.model.events.length).reduce((sum,c)=>sum+c.fraction,0);return {name:r.name,share,wt:estimable?zero/share:null,indel:estimable?1-zero/share:null,inFrame:estimable?rows.filter(c=>c.model.events.length&&net(c.model.events)%3===0).reduce((sum,c)=>sum+c.fraction,0)/share:null,outOfFrame:estimable?rows.filter(c=>net(c.model.events)%3!==0).reduce((sum,c)=>sum+c.fraction,0)/share:null,components:[...bins].map(([size,fraction])=>({size,fraction:share?fraction/share:0})),estimable,rows};});
 const targets=perRef.flatMap((a,h)=>(refs[h].input.cuts??[refs[h].input.cut]).map((cut,t)=>({reference:a.name,target:t+1,cut,edited:a.estimable?a.rows.filter(c=>inputEvents(h,c.model.events).some(e=>e.target===t+1||e.target===0)).reduce((sum,c)=>sum+c.fraction,0)/a.share:null})));
 const annotations=perRef.map((a,h)=>({name:a.name,reference:refs[h].input.sequence,cut:refs[h].input.cut,cuts:refs[h].input.cuts,direction:refs[h].reverse?'reverse' as const:'forward' as const,share:a.share,estimable:a.estimable,components:a.rows.map(c=>({size:net(c.model.events),fraction:c.fraction/a.share,events:ambiguous.includes(c)?[]:inputEvents(h,c.model.events),unresolved:ambiguous.includes(c)}))}));
 const bins=new Map<number,number>();chosen.forEach(c=>bins.set(net(c.model.events),(bins.get(net(c.model.events))||0)+c.fraction));const zero=chosen.filter(c=>!c.model.events.length).reduce((sum,c)=>sum+c.fraction,0);
 const first=Math.min(...refs.flatMap(r=>r.cuts)),plots=[{trace:input.control,trim:cr.trim},{trace:input.sample,trim:sr.trim}].map(({trace,trim})=>plotTrace(trace,first-trim,Math.max(...refs.flatMap(r=>r.cuts))-first+s.window/2));report.positions=positions;report.joint!.candidateCount=models.length;
 if(audit)audit({refs,positions,cols,models,ids,weights:fit.weights,theta,sampleLag:sc.lag,gains:sc.gains,observed:y,rawObserved:positions.flatMap(p=>raw[p-sr.trim]),anchor,end,iterations,train,test,tune,controlPositionRmse:ctlScore.positionRmse});
 return {method:'pooled',referenceName:refs.length===1?refs[0].name:'Mixed amplicons',sampleName:input.sample.name,controlName:input.control.name,direction:'forward',cutPosition:first,alignmentIdentity:0,alignment:[],variants:[],plots,warnings,synthetic:!!input.synthetic,annotations,diagnostics:{algorithmVersion:'sanger-'+(s.jointEngine||'signal-path')+'-v4-known',signalModel:'calibrated',requestedMaxIndel:s.maxIndel,discordantPositions:localDiscordant,maxPositionRmse:Math.max(...localResidual),controlCalibration:cc,sampleCalibration:sc,contextEmission:context},mixture:{components:[...bins].map(([size,fraction])=>({size,fraction})),wt:zero,indel:1-zero,inFrame:chosen.filter(c=>c.model.events.length&&net(c.model.events)%3===0).reduce((sum,c)=>sum+c.fraction,0),outOfFrame:chosen.filter(c=>net(c.model.events)%3!==0).reduce((sum,c)=>sum+c.fraction,0),r2:full.r2,rmse:full.rmse,upstreamError:1-sc.validationR2,bases:positions.length,fitStart:anchor,fitEnd:end,converged:fit.converged,passesChecks:allowed&&!attributionAmbiguous.length},pooled:{preflight:report,controlR2:ctlScore.r2,sampleTrim:sr.trim,controlTrim:cr.trim,amplicons:perRef.map(({rows,estimable,...a})=>a)},joint:{groups,targets,reconstructionR2:qc.r2,heldOutR2:qc.r2,selectedModels:chosen.length,candidateModels:models.length,iterations,maxDeletion,ambiguousFraction:ambiguous.reduce((sum,c)=>sum+c.fraction,0),alignmentSearches:alignments}};
}
