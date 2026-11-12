import {nativeColumnResult} from './sanger-native-column-result';
import {nativeCompatibility,CompatibleRow} from './sanger-native-compatibility';
import{rawIntervalScorer}from './sanger-raw-interval-fit';
import{SangerSequenceFitCache,sequenceFitLoss,sequenceFitTolerance}from'./sanger-sequence-fit-cache';
import{refineDiscretePhase,SangerLinearOptimizer}from'./sanger-discrete-phase';
import {frequencyProfiles} from './sanger-stratified-paths';
import{initializeOffsetProfiles}from'./sanger-offset-groups';
import{SignalPathInput,Prepared,PathEdit,Model,prepare,materialize}from'./sanger-path-common';
import{repairSignalTraceCycles}from'./sanger-cycle';
import{contrastPeaks,pairedReferencePeaks,correctGain,rawSignalScore}from'./sanger-signal';
import{fitNormalizedEmissionMixture}from'./sanger-emission';
import{SangerFitCache}from'./sanger-fit-cache';
import{nonnegativeLeastSquares}from'./sanger-nonnegative-ls';
import{decodeFrequencyPaths}from'./sanger-frequency-path-decoder';
import{sequenceGroupResult}from'./sanger-sequence-group-result';
const bases='ACGT';
interface Family {size:number;fraction:number}
/** Prefix measurements determine upstream length marginals. A convolution
 * inverse proposes downstream profiles; it is ONLY an initializer. Final
 * whole-molecule pairing coefficients are free and scored on measured signal.
 * Unsupported observations are rejected, never silently treated as independent. */
export function analyzeConditionalGroups(input:SignalPathInput,progress:(s:string)=>void=()=>{},optimizer?:SangerLinearOptimizer,audit?:(stage:string,milliseconds:number,detail?:Record<string,number>)=>void){
 if(!optimizer)throw Error('The local discrete optimizer must be initialized before conditional inference.');
 // A quickly certified local factorization preserves its compact representative.
 // The time budget switches search strategy; it never limits the allele count.
 const started=Date.now(),localStarted=performance.now();const candidates:CompatibleRow[]=[];let preliminaryError:unknown,preliminaryResult:ReturnType<typeof conditionalProfiles>|undefined;
 try{const result=conditionalProfiles(input,stage=>{if(Date.now()-started>30000)throw Error('Local initialization work budget exhausted');progress(stage);},optimizer,rows=>candidates.push(...rows));preliminaryResult=result;if(result.diagnostics?.nativeRawTracesVerified){audit?.('local-total',performance.now()-localStarted);return result;}}
 catch(error){preliminaryError=error;}
 audit?.('local-total',performance.now()-localStarted);
 progress('Switching to full-molecule column generation');
 const exact=nativeColumnResult(input,optimizer,progress,audit,candidates);if(exact)return exact;
 if(preliminaryError)throw preliminaryError;
 return preliminaryResult!;
}
function conditionalProfiles(input:SignalPathInput,progress:(s:string)=>void,optimizer:SangerLinearOptimizer,onCandidates?:(rows:CompatibleRow[])=>void){
 if(!optimizer)throw Error('The local discrete optimizer must be initialized before conditional inference.');
 const controlRepair=repairSignalTraceCycles(input.control!),sampleRepair=repairSignalTraceCycles(input.sample!);input={...input,control:controlRepair.trace,sample:sampleRepair.trace};
 const s=input.settings,{refs,report,anchor,maxDeletion,maxInsertion}=prepare(input.amplicons,input.primer,s);
 if(maxDeletion>10||maxInsertion>10||refs.some(r=>r.reverse||r.cuts.length!==2||r.input.cuts?.[0]!==r.input.cut||r.input.cuts[0]>=r.input.cuts[1]))throw Error('Conditional profiles currently require two ordered forward targets and ≤10 bp length limits.');
 const prefixEnd=Math.min(...refs.map(r=>r.cuts[1]))-maxDeletion-5;
 progress('Fitting first-target measured prefix');
 const initialized=initializeOffsetProfiles({...input,amplicons:input.amplicons.map(a=>({...a,cuts:[a.cuts![0]]})),settings:{...s,readEnd:prefixEnd,editStartRadius:3}},stage=>progress('Conditional prefix · '+stage),undefined,optimizer,onCandidates),prefix=initialized.result,firstOptions=initialized.profiles,firstSelected=[...initialized.selected],cc=prefix.diagnostics!.controlCalibration!,sc=prefix.diagnostics!.sampleCalibration!,context=prefix.diagnostics!.contextEmission!;
 if(!prefix.mixture?.passesChecks||prefix.annotations!.some(a=>!a.estimable)||cc.lag!==0||sc.lag!==0||context.used||context.theta.some(v=>v!==0))throw Error('Conditional profiles require validated unit-emission prefix measurements.');
 const cr={trim:prefix.pooled!.controlTrim},sr={trim:prefix.pooled!.sampleTrim},peaks=pairedReferencePeaks(input.control!,input.sample!,refs.map(r=>r.sequence),cr.trim,sr.trim,Math.min(...refs.flatMap(r=>r.cuts)),optimizer),ctl=peaks.control.vectors,raw=peaks.sample.vectors,end=Math.min(s.readEnd,cr.trim+ctl.length-3,sr.trim+raw.length-3,...refs.map(r=>r.sequence.length-2*maxDeletion-2)),positions=Array.from({length:end-anchor},(_,i)=>i+anchor),rawY=positions.flatMap(p=>raw[p-sr.trim]),y=positions.flatMap(p=>correctGain(raw[p-sr.trim],sc.gains));
 const emit=(sequence:string,p:number,_lag=0)=>[...bases].map(b=>Number(sequence[p]===b)),controlCols=refs.map(r=>positions.flatMap(p=>emit(r.sequence,p))),controlFit=fitNormalizedEmissionMixture(controlCols,positions.flatMap(p=>correctGain(ctl[p-cr.trim],cc.gains))),ctlScore=rawSignalScore(controlCols,positions.flatMap(p=>ctl[p-cr.trim]),controlFit.weights,cc.gains);
 if(ctlScore.r2<.995)throw Error('The full control does not validate the conditional unit-emission model.');
 const first=prefix.annotations!.map(a=>{const bins=new Map<number,number>();for(const c of a.components)bins.set(c.size,(bins.get(c.size)||0)+c.fraction*a.share);return[...bins].map(([size,fraction])=>({size,fraction}));}),shares=first.map(fs=>fs.reduce((sum,f)=>sum+f.fraction,0));
 progress('Fitting downstream total-length marginals');
 const tail=positions.filter(p=>p>Math.max(...refs.map(r=>r.cuts[1]))+35);if(tail.length<80)throw Error('Insufficient downstream invariant signal for conditional profiles.');
 const totalClasses=refs.flatMap((_,h)=>Array.from({length:2*(maxDeletion+maxInsertion)+1},(_,i)=>({reference:h,size:i-2*maxDeletion}))),tailColumns=totalClasses.map(f=>tail.flatMap(p=>[...bases].map(b=>Number(refs[f.reference].sequence[p-f.size]===b)))),tailY=tail.flatMap(p=>correctGain(raw[p-sr.trim],sc.gains)),tailFit=new SangerFitCache(tailColumns,tailY).fit(totalClasses.map((_,i)=>i)),totals=refs.map((_,h)=>{const own=totalClasses.flatMap((f,i)=>f.reference===h?[{size:f.size,fraction:tailFit.weights[i]}]:[]),mass=own.reduce((sum,f)=>sum+f.fraction,0);return own.map(f=>({...f,fraction:f.fraction*shares[h]/mass}));});
 if(!tailFit.converged)throw Error('Downstream length fit did not converge.');
 const options:{events:PathEdit[];fraction:number}[][][]=[],proposalDetails:any[]=[];
 for(let h=0;h<refs.length;h++){
  progress(`Proposing second-target profiles · reference ${h+1} / ${refs.length}`);
  const r=refs[h],share=shares[h],u=first[h].map(f=>({...f,fraction:f.fraction/share})),sizes=Array.from({length:maxDeletion+maxInsertion+1},(_,i)=>i-maxDeletion),C=sizes.map(size=>totals[h].map(c=>u.find(a=>a.size===c.size-size)?.fraction||0)),lengthFit=nonnegativeLeastSquares(C,totals[h].map(c=>c.fraction));if(!lengthFit.converged)throw Error('Conditional length proposal did not converge.');
  const mass=lengthFit.weights.reduce((a,b)=>a+b,0),v=sizes.map((size,i)=>({size,fraction:lengthFit.weights[i]*share/mass})),cut=r.cuts[1],q=Array.from({length:17},(_,i)=>cut-3+i),from=q[0]+Math.min(...u.map(a=>a.size))-2,to=q.at(-1)!+Math.max(...u.map(a=>a.size))+3,local=Array.from({length:to-from},(_,i)=>from+i),K=q.map(j=>local.map(p=>u.find(a=>a.size===p-j)?.fraction||0));
  const observed=local.map(p=>{const row=correctGain(raw[p-sr.trim],sc.gains);for(let k=0;k<refs.length;k++)if(k!==h){if(p>=refs[k].cuts[1]-maxDeletion-3&&p<=refs[k].cuts[1]+maxInsertion+13)throw Error('Conditional target windows overlap.');const fs=p<refs[k].cuts[1]-13?first[k]:totals[k];for(const f of fs){const b=bases.indexOf(refs[k].sequence[p-f.size]);if(b>=0)row[b]-=f.fraction;}}return row;});
  const known=local.map(p=>{const row=[0,0,0,0];for(const a of u){const j=p-a.size;if(j>=q[0]&&j<=q.at(-1)!)continue;if(j<q[0]){const b=bases.indexOf(r.sequence[j]);if(b>=0)row[b]+=a.fraction*share;}else for(const f of v){const b=bases.indexOf(r.sequence[j-f.size]);if(b>=0)row[b]+=a.fraction*f.fraction;}}return row;}),baseFits=[...bases].map((_,b)=>nonnegativeLeastSquares(K,observed.map((row,i)=>row[b]-known[i][b])));if(baseFits.some(f=>!f.converged))throw Error('Conditional base proposal did not converge.');
  const F=q.map((_,i)=>{const row=baseFits.map(f=>f.weights[i]),total=row.reduce((a,b)=>a+b,0);return row.map(x=>total?x*share/total:0);}),localPositions=Array.from({length:end-30},(_,i)=>i+30),localObserved=localPositions.flatMap(p=>p>=q[0]&&p<=q.at(-1)!?F[p-q[0]]:p<q[0]?[...bases].map(b=>share*Number(r.sequence[p]===b)):[...bases].map(b=>v.reduce((sum,f)=>sum+f.fraction*Number(r.sequence[p-f.size]===b),0))),keep=v.filter(f=>f.size!==0&&f.fraction>=share*.01),families=[{reference:0,size:0,fraction:share-keep.reduce((sum,f)=>sum+f.fraction,0)},...keep.map(f=>({...f,reference:0}))],single:Prepared={...r,cuts:[cut]},decoded=frequencyProfiles(single,families,localPositions,localObserved,.001,200,stage=>progress('Reference '+(h+1)+' / '+refs.length+' · '+stage),0,optimizer);
  if(!decoded.options.length)throw Error('No supported second-target profile proposal · reference '+r.name+' '+JSON.stringify(decoded.details));options.push(decoded.options);onCandidates?.(decoded.options.flatMap(option=>option.map(g=>({reference:h,events:g.events,fraction:g.fraction}))));proposalDetails.push(decoded.details);
 }
 progress('Jointly fitting free target pairings');
 const models:Model[]=[],sequences:string[]=[],keys=new Map<string,number>(),profileIds=new Map<string,number[]>(),add=(h:number,events:PathEdit[])=>{const sequence=materialize(refs[h].sequence,events);if(!sequence)throw Error('Invalid conditional molecule.');const key=h+':'+sequence,old=keys.get(key),alias={reference:h,events};if(old!==undefined){if(!models[old].aliases.some(a=>JSON.stringify(a.events)===JSON.stringify(events)))models[old].aliases.push(alias);return old;}const id=models.length;keys.set(key,id);models.push({reference:h,sequence,events,aliases:[alias]});sequences.push(sequence);return id;},build=(h:number,index:number,firstIndex=firstSelected[h])=>{const key=h+':'+firstIndex+':'+index,hit=profileIds.get(key);if(hit)return hit;const ids:number[]=[];for(const firstEvents of firstOptions[h][firstIndex])for(const b of options[h][index])ids.push(add(h,[...firstEvents,...b.events]));const unique=[...new Set(ids)];profileIds.set(key,unique);return unique;};
 const selected=refs.map(()=>0);selected.forEach((index,h)=>build(h,index));const fitCache=new SangerSequenceFitCache(sequences,positions,y);let previousFit:{ids:number[];weights:number[]}|undefined;const fitProfiles=(indices:number[],firstIndices=firstSelected,cache=fitCache)=>{const ids=indices.flatMap((index,h)=>build(h,index,firstIndices[h]));const old=new Map(previousFit?.ids.map((id,i)=>[id,previousFit!.weights[i]])||[]),fit=cache.fit(ids,ids.map(id=>old.get(id)||0));previousFit={ids,weights:fit.weights};return{ids,fit};};let best=fitProfiles(selected),iterations=0;
 const optimizeProfiles=(cache=fitCache)=>{if(best.fit.rawResponse&&sequenceFitLoss(best.fit)<1e-7)return;for(let pass=0;pass<4;pass++){iterations++;
  let changed=false;
  for(let h=0;h<refs.length;h++){
   progress(`Comparing joint target hypotheses · reference ${h+1} / ${refs.length} · pass ${pass+1}`);
   let choice=selected[h],firstChoice=firstSelected[h];
   for(let a=0;a<firstOptions[h].length;a++)for(let b=0;b<options[h].length;b++){
    if((a*options[h].length+b)%32===0)progress(`Comparing joint target hypotheses · reference ${h+1} / ${refs.length} · ${a*options[h].length+b+1} / ${firstOptions[h].length*options[h].length}`);const firstTrial=[...firstSelected],secondTrial=[...selected];firstTrial[h]=a;secondTrial[h]=b;const candidate=fitProfiles(secondTrial,firstTrial,cache);
    if(candidate.fit.converged&&sequenceFitLoss(candidate.fit)<sequenceFitLoss(best.fit)-sequenceFitTolerance(best.fit)){choice=b;firstChoice=a;best=candidate;if(best.fit.rawResponse&&sequenceFitLoss(best.fit)<1e-7){selected[h]=choice;firstSelected[h]=firstChoice;return;}}
   }
   if(choice!==selected[h]||firstChoice!==firstSelected[h]){selected[h]=choice;firstSelected[h]=firstChoice;changed=true;}
  }
  if(!changed)return;
 }
 throw Error('Conditional profile search did not stabilize within its budget.');};
 const rawScorer=peaks.method==='quantized-peak-shape'?rawIntervalScorer(input.control!,input.sample!,refs.map(r=>r.sequence),positions,cr.trim,sr.trim,peaks.control.background!,peaks.sample.background!,peaks.calibration,optimizer,y):null;
 const finalCache=rawScorer?new SangerSequenceFitCache(sequences,positions,y,512,rawScorer):fitCache;
 let certifiedRows:CompatibleRow[]|null=null;
 const compatible=(rows:CompatibleRow[])=>{if(!rawScorer)return false;const candidates=[...rows];for(let h=0;h<refs.length;h++){candidates.push({reference:h,events:[],fraction:0});for(const option of firstOptions[h])for(const events of option)candidates.push({reference:h,events,fraction:0});for(const option of options[h])for(const g of option)candidates.push({reference:h,events:g.events,fraction:0});}const proof=nativeCompatibility(refs,candidates,input.control!,input.sample!,optimizer,.25);if(!proof)return false;certifiedRows=proof.rows;return true;};
 const certifiedResult=()=>{const masses=new Map<number,number>();for(const row of certifiedRows!){const id=add(row.reference,row.events);masses.set(id,(masses.get(id)||0)+row.fraction);}const ids=[...masses.keys()],weights=ids.map(id=>masses.get(id)!),chosen=new Set(ids),cols=models.map((_,id)=>chosen.has(id)?fitCache.column(id):[]);let error=0;const mean=y.reduce((a,b)=>a+b,0)/y.length,variance=y.reduce((n,v)=>n+(v-mean)**2,0);for(let k=0;k<y.length;k++)error+=(weights.reduce((n,w,i)=>n+w*cols[ids[i]][k],0)-y[k])**2;const fit={weights,rmse:Math.sqrt(error/y.length),r2:variance?1-error/variance:1,converged:true},train=positions.map((_,i)=>i);const result=sequenceGroupResult({input,refs,positions,report,anchor,end,models,cols,ids,fit,ctlScore,cc,sc,context,cr,sr,rawY,controlRepair,sampleRepair,maxDeletion,maxInsertion,fitCache:finalCache,y,train,test:train,tune:train,emit,iterations:iterations+1,algorithmVersion:'sanger-sequence-groups-v44-threshold-profiles'});return{...result,diagnostics:{...result.diagnostics,signalModel:result.diagnostics?.signalModel||'calibrated',requestedMaxIndel:result.diagnostics?.requestedMaxIndel||s.maxIndel,nativeRawTracesVerified:true}};};
 optimizeProfiles();
 best=fitProfiles(selected);
 if(compatible(best.ids.flatMap((id,i)=>best.fit.weights[i]>1e-9?[{reference:models[id].reference,events:models[id].events,fraction:best.fit.weights[i]}]:[]))){progress('Certified complete raw-trace compatibility');return certifiedResult();}
 progress('Resolving joint sequence phase');
 if(rawScorer){best=fitProfiles(selected,firstSelected,finalCache);optimizeProfiles(finalCache);}
 if(compatible(best.ids.flatMap((id,i)=>best.fit.weights[i]>1e-9?[{reference:models[id].reference,events:models[id].events,fraction:best.fit.weights[i]}]:[]))){progress('Certified complete raw-trace compatibility');return certifiedResult();}

 progress(rawScorer?'Raw ADC interval phase · centroid bound '+rawScorer.maximumCentroidMse:'Using measured centroid phase');
 for(let round=0;round<(rawScorer?2:1);round++){
  let changed=false;
  let targets=rawScorer?[0,1]:[1];
  if(rawScorer){const measured=rawScorer(best.ids.map(id=>sequences[id]),best.fit.weights);if(measured&&!measured.lowerBoundOnly){if((measured.loss??Infinity)<1e-7)break;const error=(target:number)=>(measured.residuals||[]).reduce((n,r)=>n+(refs.some((ref,h)=>{const shifts=target===0?[0]:firstOptions[h][firstSelected[h]].map(es=>es.reduce((s,e)=>s+e.inserted.length-e.end+e.start,0));return r.position>=ref.cuts[target]-3+Math.min(...shifts)&&r.position<=ref.cuts[target]+13+Math.max(...shifts);})?r.error:0),0);targets=targets.sort((a,b)=>error(b)-error(a)||a-b);}}
  for(const target of targets){
   progress(`Resolving full-molecule phase · round ${round+1} · target ${target+1}`);
   const firstEvents=firstSelected.map((index,h)=>firstOptions[h][index]),secondEvents=selected.map((index,h)=>options[h][index].map(g=>g.events)),phased=refineDiscretePhase(refs,target===1?firstEvents:secondEvents,target===1?secondEvents:firstEvents,positions,y,optimizer,progress,target,true,rawScorer||undefined,compatible);
   if(certifiedRows){progress("Certified complete raw-trace compatibility");return certifiedResult();}
   if(phased.changes){phased.profiles.forEach((profile,h)=>{if(target===1){options[h].push(profile.map(events=>({events,fraction:0})));selected[h]=options[h].length-1;}else{firstOptions[h].push(profile);firstSelected[h]=firstOptions[h].length-1;}});best=fitProfiles(selected,firstSelected,finalCache);if(!rawScorer)optimizeProfiles();changed=true;}
  }
  if(!changed)break;
 }
 if(!best.fit.converged)throw Error('Final whole-signal coefficient fit did not converge.');
 progress('Validating fitted conditional profiles');
 const chosenIds=new Set(best.ids),cols=models.map((_,id)=>chosenIds.has(id)?fitCache.column(id):[]),train=positions.map((_,i)=>i),result=sequenceGroupResult({input,refs,positions,report,anchor,end,models,cols,ids:best.ids,fit:best.fit,ctlScore,cc,sc,context,cr,sr,rawY,controlRepair,sampleRepair,maxDeletion,maxInsertion,fitCache:finalCache,y,train,test:train,tune:train,emit,iterations:iterations+1,algorithmVersion:'sanger-sequence-groups-v44-threshold-profiles'});
 return result;
}
