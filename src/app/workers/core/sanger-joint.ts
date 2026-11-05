import {PooledAmplicon,PooledPreflight,SangerAnalysisSettings,SangerAnalysisResult,SangerTrace,TRACE_BASES,SangerJointEvent} from '../../models/sanger-analysis';
import {exactBindingSites} from './sanger-primer-designer';
import {reverseComplement} from '../../utils/biology.utils';
import {registerTrace} from './sanger-pooled';
import {contrastPeaks,integratedPeaks,calibrateSignal,correctGain,simplexFit,rawSignalScore,gramCondition} from './sanger-signal';
import {emissionVector,fitContextEmission,sequenceYield,fitEmissionMixture} from './sanger-emission';
import {plotTrace} from './sanger-analysis';

interface Ref {name:string;sequence:string;cuts:number[];input:PooledAmplicon;reverse:boolean;primerEnd:number}
interface Event {target:number;start:number;end:number;insert:number}
interface Proposal {reference:number;events:Event[];sequence:string;size:number;aliases:{reference:number;events:Event[]}[]}
interface Design {refs:Ref[];models:Proposal[];report:PooledPreflight;positions:number[]}
const blank=(reason=''):PooledPreflight=>({passed:false,reason,positions:[],modelCount:0,rank:0,condition:null,minimumDifferences:0,minimumCutOffset:0,amplicons:[]});
export const usesJointSanger=(amplicons:PooledAmplicon[],s:SangerAnalysisSettings)=>!!s.largeDeletions||!!(s.jointEngine&&s.jointEngine!=='staged')||amplicons.some(a=>(a.cuts?.length||1)>1);
/** Edit complete primer-oriented molecules. Downstream coordinates include ALL upstream edits. */
export function applyJointEvents(sequence:string,events:Event[]):string|null {
 let at=0,result='';for(const e of [...events].sort((a,b)=>a.start-b.start||a.end-b.end)){
  if(e.start<at||e.start<0||e.end<e.start||e.end>sequence.length)return null;
  result+=sequence.slice(at,e.start)+'N'.repeat(e.insert);at=e.end;
 }return result+sequence.slice(at);
}
function prepare(amplicons:PooledAmplicon[],primer:string,s:SangerAnalysisSettings):Design{
 const report=blank();
 if(amplicons.length<1||amplicons.length>4||new Set(amplicons.map(a=>a.name)).size!==amplicons.length)throw Error('Use 1–4 uniquely named submitted amplicons.');
 if(!/^[ACGT]{12,40}$/.test(primer))throw Error('Enter the actual sequencing primer (12–40 nt, 5′→3′).');
 if(s.method!=='mixture'||s.signalModel==='classic')throw Error('Joint targets / large deletions require calibrated trace decomposition.');
 const range=s.maxIndel,large=s.largeDeletions?(s.maxDeletion??200):range;
 if(!Number.isInteger(range)||range<1||range>15||!Number.isInteger(large)||large<range||large>500)throw Error('Use small indels ±1–15 bp and a deletion search limit from that range to 500 bp.');
 if(!Number.isInteger(s.window)||s.window<60||s.window>500||!Number.isInteger(s.readStart)||!Number.isInteger(s.readEnd)||s.readStart<0||s.readEnd>1000||s.readEnd<=s.readStart)throw Error('Use a 60–500 bp window and a usable primer extension range within 0–1,000 bp.');
 const refs=amplicons.map(a=>{
  const cuts=a.cuts??[a.cut];if(!a.name.trim()||!/^[ACGT]+$/.test(a.sequence)||a.sequence.length>20000||cuts.length<1||cuts.length>2||new Set(cuts).size!==cuts.length||cuts.some(c=>!Number.isInteger(c)||c<=0||c>=a.sequence.length))throw Error(a.name+': enter one or two distinct cuts inside the actual submitted amplicon.');
  const hits=exactBindingSites(primer,{name:a.name,sequence:a.sequence,cuts});if(hits.length!==1)throw Error(a.name+': sequencing primer must bind exactly once across both strands.');
  const h=hits[0],reverse=h.direction==='reverse',oriented=reverse?reverseComplement(a.sequence):a.sequence,end=reverse?a.sequence.length-h.start:h.end;
  const offsets=cuts.map(c=>(reverse?a.sequence.length-c:c)-end).sort((a,b)=>a-b),sequence=oriented.slice(end);
  const windows=offsets.map(c=>({start:c-Math.floor(s.window/2),end:c+Math.ceil(s.window/2)}));
  const maxWindow=Math.min(500,Math.max(0,Math.floor(Math.min(...offsets.map(c=>2*Math.min(c-s.readStart,Math.min(s.readEnd,sequence.length)-c))))));
  if(offsets[0]<s.readStart+72||offsets.at(-1)!+30>Math.min(s.readEnd,sequence.length))throw Error(`${a.name}: target windows exceed the primer-relative read range; maximum centered window is ${maxWindow} bp, with ≥72 upstream bases required.`);
  if(offsets.length===2&&offsets[1]-offsets[0]<=2*range+20)throw Error(`${a.name}: targets are too close for independent small-indel attribution at ±${range} bp. Reduce the indel range or use one target. Maximum non-overlapping centered window: ${offsets[1]-offsets[0]} bp.`);
  const lo=windows[0].start,hi=windows.at(-1)!.end;
  report.amplicons.push({name:a.name,direction:h.direction,primerStart:h.start+1,primerEnd:h.end,cutOffset:offsets[0],windowStart:lo,windowEnd:hi,targets:offsets.map((cut,i)=>({cutOffset:cut,windowStart:windows[i].start,windowEnd:windows[i].end})),maximumWindow:maxWindow,maximumIndependentWindow:offsets.length===2?offsets[1]-offsets[0]:undefined,mergedWindows:windows.length===2&&windows[0].end>windows[1].start});
  return {name:a.name,sequence,cuts:offsets,input:a,reverse,primerEnd:end};
 });
 const start=s.readStart,end=Math.min(s.readEnd,Math.max(...refs.map(r=>r.sequence.length)));
 
 // Do not give insertion identity an invented base. Mask every possible insertion
 // position, including displacement of target 2 caused by target 1.
 const positions=Array.from({length:end-start},(_,i)=>start+i).filter(p=>!refs.some(r=>r.cuts.some((c,i)=>p>=c-(i?range:0)-5&&p<c+(i?2*range:range)+5)));
 if(positions.length<100)throw Error('Too few joint observation positions. Increase the window or reduce the indel range.');
 const models:Proposal[]=[],seen=new Map<string,Proposal>();let skipped=0;
 const add=(reference:number,events:Event[])=>{
  const seq=applyJointEvents(refs[reference].sequence,events);if(events.some(e=>e.start<s.readStart+84)||seq===null||seq.length<s.readStart+84||positions.some(p=>p<seq.length&&!/[ACGT]/.test(seq[p]))){skipped++;return;}
  const key=positions.map(p=>seq[p]).join(''),alias={reference,events};const previous=seen.get(key);
  if(previous){previous.aliases.push(alias);return;}
  const m={reference,events,sequence:seq,size:events.reduce((sum,e)=>sum+e.insert-(e.end-e.start),0),aliases:[alias]};models.push(m);seen.set(key,m);
 };
 refs.forEach((ref,h)=>{
  const choices=ref.cuts.map((cut,target)=>Array.from({length:2*range+1},(_,j)=>j-range).map(size=>size===0?[]:[{target,start:cut,end:cut+Math.max(0,-size),insert:Math.max(0,size)}]));
  for(const a of choices[0]){if(choices.length===1)add(h,a);else for(const b of choices[1])add(h,[...a,...b]);}
  if(s.largeDeletions){ref.cuts.forEach((cut,target)=>{for(let length=range+1;length<=large;length++)for(const left of [0,Math.floor(length/2),length])add(h,[{target,start:cut-left,end:cut-left+length,insert:0}]);});}
  if(ref.cuts.length===2){const [left,right]=ref.cuts;if(right-left<=large||!s.largeDeletions)for(let dl=-range;dl<=range;dl++)for(let dr=-range;dr<=range;dr++){
   const e={target:-1,start:left+dl,end:right+dr,insert:0};if(e.end>e.start&&(!s.largeDeletions||e.end-e.start<=large))for(let insert=0;insert<=range;insert++)add(h,[{...e,insert}]);
  }}
 });
 if(models.length*positions.length*4>12_000_000)throw Error('Joint proposal dictionary exceeds the browser memory budget. Reduce the small-indel range, deletion limit or target window.');
 const wt=refs.map((_,h)=>models.find(m=>m.aliases.some(a=>a.reference===h&&!a.events.length))!);
 let minimum=Infinity;for(let i=0;i<wt.length;i++)for(let j=0;j<i;j++)minimum=Math.min(minimum,positions.filter(p=>wt[i].sequence[p]!==wt[j].sequence[p]).length);
 report.positions=positions;report.modelCount=models.length;report.minimumDifferences=Number.isFinite(minimum)?minimum:positions.length;report.minimumCutOffset=refs.length>1?Math.min(...refs.flatMap((r,i)=>refs.slice(i+1).map(t=>Math.abs(r.cuts[0]-t.cuts[0])))):0;
 const geometry=gramCondition(wt.map(m=>positions.flatMap(p=>TRACE_BASES.map(b=>Number(m.sequence[p]===b)))));report.rank=geometry.rank;report.condition=geometry.condition;
 if(report.minimumDifferences<3||geometry.rank<refs.length||geometry.condition===null||geometry.condition>100)throw Error('WT references cannot be distinguished on the submitted primer axis. Change the primer, window or references.');
 report.passed=true;report.reason='Joint proposal geometry checked. Selected components require held-out trace validation.';
 report.joint={mode:'staged',candidateCount:models.length,ambiguousClasses:models.filter(m=>m.aliases.length>1).length,coverageExcluded:skipped,maxDeletion:large};
 return {refs,models,report,positions};
}
export function checkJointSeparation(amplicons:PooledAmplicon[],primer:string,s:SangerAnalysisSettings):PooledPreflight{try{return prepare(amplicons,primer,s).report;}catch(e){return blank(e instanceof Error?e.message:String(e));}}
export interface JointInput {amplicons:PooledAmplicon[];primer:string;ampliconConfirmed:boolean;sample:SangerTrace;control:SangerTrace;settings:SangerAnalysisSettings;synthetic?:boolean}
export function analyzeJointSanger(input:JointInput,progress:(s:string)=>void=()=>{}):SangerAnalysisResult{
 if(!input.ampliconConfirmed)throw Error('Confirm the actual full amplicons submitted to Sanger.');
 const s=input.settings;let {refs,models,report,positions}=prepare(input.amplicons,input.primer,s);
 if(!input.sample||!input.control)throw Error('Joint analysis requires an unedited control and a sample using the same sequencing primer.');
 progress('Registering joint traces on the sequencing-primer axis');
 const anchor=Math.min(...models.flatMap(m=>m.events.map(e=>e.start)));
 const prepared=refs.map(r=>({name:r.name,sequence:r.sequence,cut:anchor}));
 const ctlReg=registerTrace(input.control,prepared,s,integratedPeaks(input.control,[0,0,0,0]).vectors),sampleReg=registerTrace(input.sample,prepared,s,integratedPeaks(input.sample,[0,0,0,0]).vectors);
 const ctl=contrastPeaks(input.control).vectors,raw=contrastPeaks(input.sample).vectors,first=Math.min(...refs.map(r=>r.cuts[0]));
 const actualStart=Math.max(s.readStart,ctlReg.trim+10,sampleReg.trim+10),actualEnd=Math.min(s.readEnd,ctlReg.trim+ctl.length-2,sampleReg.trim+raw.length-2);
 if(actualStart>s.readStart||positions.at(-1)!>=actualEnd)({refs,models,report,positions}=prepare(input.amplicons,input.primer,{...s,readStart:actualStart,readEnd:actualEnd}));

 const calibrate=(vectors:number[][],trim:number)=>{const observed=Array.from({length:vectors.length+trim},(_,p)=>vectors[p-trim]||[0,0,0,0]),rows=Array.from({length:anchor-12},(_,p)=>p).filter(p=>p>=Math.max(s.readStart,trim+10)&&p-trim<vectors.length-2);return calibrateSignal(observed,rows,refs.map(r=>Array.from(r.sequence,(_,p)=>TRACE_BASES.map(b=>Number(r.sequence[p]===b)))));};
 const cc=calibrate(ctl,ctlReg.trim),sc=calibrate(raw,sampleReg.trim);if(cc.validationR2<.8||sc.validationR2<.8)throw Error('Joint upstream calibration failed held-out validation.');
 if(positions.some(p=>p-ctlReg.trim<0||p-sampleReg.trim<0||p-ctlReg.trim>=ctl.length||p-sampleReg.trim>=raw.length))throw Error('Actual chromatograms do not cover the checked joint windows.');
 const corrected=raw.map(v=>correctGain(v,sc.gains)),wtCorrected=ctl.map(v=>correctGain(v,cc.gains));
 const observedControl=Array.from({length:ctl.length+ctlReg.trim},(_,p)=>wtCorrected[p-ctlReg.trim]||[0,0,0,0]);
 const context=fitContextEmission(refs.map(r=>r.sequence),positions,observedControl,cc.weights,cc.lag),theta=context.theta;
 const cols=models.map(m=>positions.flatMap(p=>emissionVector(m.sequence,p,sc.lag,theta)));
 const wtIds=refs.map((_,h)=>models.findIndex(m=>m.aliases.some(a=>a.reference===h&&!a.events.length)));
 const controlCols=wtIds.map(i=>positions.flatMap(p=>emissionVector(models[i].sequence,p,cc.lag,theta))),controlFit=fitEmissionMixture(controlCols,positions.flatMap(p=>wtCorrected[p-ctlReg.trim]));
 const controlScore=rawSignalScore(controlCols,positions.flatMap(p=>ctl[p-ctlReg.trim]),controlFit.weights,cc.gains);if(controlScore.r2<.85)throw Error('Submitted WT amplicons do not explain the pooled control.');
 const y=positions.flatMap(p=>corrected[p-sampleReg.trim]);
 // Three disjoint folds: selection training, proposal tuning, untouched final QC.
 const train=positions.flatMap((_,i)=>i%5<3?[i]:[]),tune=positions.flatMap((_,i)=>i%5===3?[i]:[]),test=positions.flatMap((_,i)=>i%5===4?[i]:[]);
 if(tune.length<20||test.length<20)throw Error('Not enough independent positions for joint candidate selection and held-out QC.');
 const subset=(c:number[],rows:number[])=>rows.flatMap(i=>c.slice(i*4,i*4+4)),trainCols=cols.map(c=>subset(c,train)),trainY=subset(y,train),tuneY=subset(y,tune);
 const fitIds=(ids:number[])=>fitEmissionMixture(ids.map(i=>trainCols[i]),trainY);
 const mse=(ids:number[],weights:number[],rows:number[],obs:number[])=>{let error=0;for(let j=0;j<rows.length;j++){const pred=[0,0,0,0];ids.forEach((id,h)=>{for(let b=0;b<4;b++)pred[b]+=weights[h]*cols[id][4*rows[j]+b];});const total=pred.reduce((a,b)=>a+b,0);for(let b=0;b<4;b++)error+=(obs[j*4+b]-(total?pred[b]/total:0))**2;}return error/obs.length;};
 let ids=[...new Set(wtIds)],fit=fitIds(ids),bestTune=mse(ids,fit.weights,tune,tuneY),iterations=0;
 progress(`Screening ${models.length} joint reference/edit proposals`);
 let refinementSkipped=0;
 for(let stage=0;stage<2;stage++){
 for(let step=0;step<24;step++,iterations++){
  const prediction=Array.from({length:trainY.length},(_,j)=>ids.reduce((sum,id,h)=>sum+fit.weights[h]*trainCols[id][j],0));
  const residual=trainY.map((v,j)=>v-prediction[j]);
  const scores=models.flatMap((m,id)=>ids.includes(id)?[]:[{id,score:trainCols[id].reduce((sum,v,j)=>sum+v*residual[j],0)}]).sort((a,b)=>b.score-a.score).slice(0,8);
  let next:{ids:number[];fit:ReturnType<typeof fitEmissionMixture>;error:number}|undefined;
  for(const candidate of scores){const trial=[...ids,candidate.id],f=fitIds(trial);if(f.weights.at(-1)!<.003)continue;const error=mse(trial,f.weights,tune,tuneY);if(error<bestTune-Math.max(.000005,bestTune*.01)&&(!next||error<next.error))next={ids:trial,fit:f,error};}
  if(!next)break;ids=next.ids;fit=next.fit;bestTune=next.error;
 }
 if(stage===0&&s.largeDeletions){
  progress('Refining supported large-deletion classes at cut-touching junctions');
  const seen=new Map(models.map(m=>[positions.map(p=>m.sequence[p]).join(''),m]));
  const classes=new Map<string,{reference:number;target:number;length:number}>();
  ids.forEach((id,i)=>{const m=models[id];if(fit.weights[i]<.01||m.events.length!==1)return;const e=m.events[0],length=e.end-e.start;if(e.target>=0&&!e.insert&&length>s.maxIndel)classes.set(`${m.reference}/${e.target}/${length}`,{reference:m.reference,target:e.target,length});});
  for(const c of classes.values())for(let start=refs[c.reference].cuts[c.target]-c.length;start<=refs[c.reference].cuts[c.target];start++){
   if(start<s.readStart+84)continue;
   const events=[{target:c.target,start,end:start+c.length,insert:0}],sequence=applyJointEvents(refs[c.reference].sequence,events);
   if(!sequence||sequence.length<=positions.at(-1)!)continue;
   const key=positions.map(p=>sequence[p]).join(''),old=seen.get(key),alias={reference:c.reference,events};
   if(old){if(!old.aliases.some(a=>a.reference===alias.reference&&JSON.stringify(a.events)===JSON.stringify(events)))old.aliases.push(alias);continue;}
   if((models.length+1)*positions.length*4>12_000_000){refinementSkipped++;continue;}
   const m={reference:c.reference,events,sequence,size:-c.length,aliases:[alias]},column=positions.flatMap(p=>emissionVector(sequence,p,sc.lag,theta));models.push(m);cols.push(column);trainCols.push(subset(column,train));seen.set(key,m);
  }
 }
 }
 // A forward-only search can retain several wrong junctions of one net-size
 // class after refinement. Refit each supported class against a single complete
 // candidate while retaining the other classes. Accept only tuning improvement;
 // genuine multi-junction mixtures are retained when the collapsed class fits worse.
 const selection=[...train,...tune].sort((a,b)=>a-b),selectionY=subset(y,selection);
 const classes=new Map<string,number[]>();for(const id of ids)if(models[id].events.length){const key=models[id].reference+'/'+models[id].size;classes.set(key,[...(classes.get(key)||[]),id]);}
 for(const members of classes.values()){
  const seed=models[members[0]],rest=ids.filter(id=>!members.includes(id)),mass=ids.reduce((sum,id,k)=>sum+(members.includes(id)?fit.weights[k]:0),0),weights=rest.map(id=>fit.weights[ids.indexOf(id)]),candidates=models.flatMap((m,id)=>m.reference===seed.reference&&m.size===seed.size&&m.events.length?[id]:[]);
  const classObjective=mse(ids,fit.weights,selection,selectionY);
  const ranked=candidates.map(id=>({id,error:mse([...rest,id],[...weights,mass],selection,selectionY)})).sort((a,b)=>a.error-b.error).slice(0,8);
  let replacement:{ids:number[];fit:ReturnType<typeof fitEmissionMixture>;error:number}|undefined;
  for(const candidate of ranked){const trial=[...rest,candidate.id],f=fitIds(trial),error=mse(trial,f.weights,selection,selectionY);if(error<classObjective-1e-7&&(!replacement||error<replacement.error))replacement={ids:trial,fit:f,error};}
  if(replacement){ids=replacement.ids;fit=replacement.fit;bestTune=mse(ids,fit.weights,tune,tuneY);}
 }
 // Drop numerical dust and refit, with no use of final QC fold for selection.
 ids=ids.filter((_,i)=>fit.weights[i]>=.003||wtIds.includes(ids[i]));fit=fitIds(ids);
 const finalTest=rawSignalScore(ids.map(id=>subset(cols[id],test)),test.flatMap(i=>raw[positions[i]-sampleReg.trim]),fit.weights,sc.gains);
 const full=rawSignalScore(ids.map(id=>cols[id]),positions.flatMap(p=>raw[p-sampleReg.trim]),fit.weights,sc.gains);
 const chosen=ids.flatMap((id,i)=>fit.weights[i]>1e-10?[{model:models[id],weight:fit.weights[i]}]:[]);
 const aliasKey=(a:{reference:number;events:Event[]})=>a.reference+':'+a.events.map(e=>`${e.target}/${e.start}/${e.end}/${e.insert}`).join(';');
 const ambiguous=chosen.filter(c=>new Set(c.model.aliases.map(aliasKey)).size>1);
 const semanticKey=(a:{reference:number;events:Event[]})=>a.reference+':'+(a.events.length?'edit':'WT')+':'+refs[a.reference].cuts.map((_,t)=>a.events.filter(e=>e.target===t||e.target===-1).reduce((sum,e)=>sum+e.insert-(e.end-e.start),0)).join('/');
 const attributionAmbiguous=ambiguous.filter(c=>new Set(c.model.aliases.map(semanticKey)).size>1);
 const net=(events:Event[])=>events.reduce((sum,e)=>sum+e.insert-(e.end-e.start),0);
 const netAmbiguous=(m:Proposal)=>new Set(m.aliases.map(a=>net(a.events))).size>1;
 const selectedGeometry=gramCondition(ids.map(id=>cols[id]));
 const warnings:string[]=[];if(refinementSkipped)warnings.push('Deletion junction refinement reached the browser memory budget; exact breakpoint search is incomplete.');if(ambiguous.length)warnings.push('Some supported models have indistinguishable explanations. Only reference/target/size assignments shared by every explanation are reported; exact breakpoints remain unresolved.');
 if(finalTest.r2<.9)warnings.push('Joint model failed untouched held-out trace QC. Results are exploratory.');
 if(report.joint!.coverageExcluded)warnings.push('Some deletion proposals exceed product/trace coverage and were excluded; the search is not exhaustive.');
 const boundary=chosen.some(c=>c.model.events.some(e=>e.end-e.start===(s.maxDeletion??200)&&s.largeDeletions));if(boundary)warnings.push('Supported deletion reaches the search boundary. Increase the deletion limit and reanalyse.');
 if(selectedGeometry.condition===null||selectedGeometry.condition>100)warnings.push('Selected components are poorly conditioned. Individual fractions are unstable.');
 const eventsInInput=(r:Ref,events:Event[]):SangerJointEvent[]=>events.map(e=>{const start=r.primerEnd+e.start,end=r.primerEnd+e.end;return {target:e.target<0?0:(r.reverse?r.cuts.length-e.target:e.target+1),start:r.reverse?r.input.sequence.length-end:start,end:r.reverse?r.input.sequence.length-start:end,inserted:e.insert};}).sort((a,b)=>a.start-b.start);
 const groups=chosen.map((c,i)=>({id:`J${i+1}`,reference:new Set(c.model.aliases.map(a=>a.reference)).size>1?null:refs[c.model.reference].name,fraction:c.weight,netSize:netAmbiguous(c.model)?null:c.model.size,events:ambiguous.includes(c)?[]:eventsInInput(refs[c.model.reference],c.model.events),alternatives:ambiguous.includes(c)?c.model.aliases.map(a=>({reference:refs[a.reference].name,events:eventsInInput(refs[a.reference],a.events)})):[]}));
 const perRef=refs.map((r,h)=>{
  const resolved=chosen.filter(c=>!attributionAmbiguous.includes(c)&&c.model.reference===h),crossAmbiguous=attributionAmbiguous.some(c=>c.model.aliases.some(a=>a.reference===h)),share=resolved.reduce((sum,c)=>sum+c.weight,0);
  const components=new Map<number,number>();resolved.forEach(c=>components.set(c.model.size,(components.get(c.model.size)||0)+c.weight));
  const noedit=resolved.filter(c=>!c.model.events.length).reduce((sum,c)=>sum+c.weight,0),estimable=!crossAmbiguous&&share>=.02&&finalTest.r2>=.9&&fit.converged;
  return {name:r.name,share,wt:estimable?noedit/share:null,indel:estimable?1-noedit/share:null,inFrame:estimable?resolved.filter(c=>c.model.events.length&&c.model.size%3===0).reduce((sum,c)=>sum+c.weight,0)/share:null,outOfFrame:estimable?resolved.filter(c=>c.model.size%3!==0).reduce((sum,c)=>sum+c.weight,0)/share:null,components:[...components].map(([size,fraction])=>({size,fraction:estimable?fraction/share:0})),estimable,groups:resolved};
 });
 const annotations=perRef.map((a,h)=>({name:a.name,reference:refs[h].input.sequence,cut:refs[h].input.cut,cuts:refs[h].input.cuts,direction:refs[h].reverse?'reverse' as const:'forward' as const,share:a.share,estimable:a.estimable,components:a.groups.map(c=>({size:c.model.size,fraction:a.share?c.weight/a.share:0,events:ambiguous.includes(c)?[]:eventsInInput(refs[h],c.model.events),unresolved:ambiguous.includes(c)}))}));
 const totalNoEdit=chosen.filter(c=>!c.model.events.length).reduce((sum,c)=>sum+c.weight,0),sizes=new Map<number,number>();chosen.forEach(c=>sizes.set(c.model.size,(sizes.get(c.model.size)||0)+c.weight));
 const targetResults=perRef.flatMap((a,h)=>refs[h].input.cuts?.map((cut,t)=>({reference:a.name,target:t+1,cut,edited:a.estimable?a.groups.filter(c=>eventsInInput(refs[h],c.model.events).some(e=>e.target===t+1||e.target===0)).reduce((sum,c)=>sum+c.weight,0)/a.share:null}))||[{reference:a.name,target:1,cut:refs[h].input.cut,edited:a.indel}]);
 const plots=[{trace:input.control,trim:ctlReg.trim},{trace:input.sample,trim:sampleReg.trim}].map(({trace,trim})=>{const p=plotTrace(trace,first-trim,Math.min(400,Math.max(...refs.flatMap(r=>r.cuts))-first+s.window/2));p.cuts=refs.flatMap(r=>r.cuts.map((cut,t)=>({name:`${r.name} T${r.reverse?r.cuts.length-t:t+1}`,x:p.cutX+(cut-first)*13}))).filter(c=>c.x<=p.width);return p;});
 return {method:'pooled',referenceName:refs.length===1?refs[0].name:'Mixed amplicons',sampleName:input.sample.name,controlName:input.control.name,direction:'forward',cutPosition:first,alignmentIdentity:0,alignment:[],variants:[],plots,warnings,synthetic:!!input.synthetic,annotations,
  diagnostics:{algorithmVersion:'sanger-joint-v1',signalModel:'calibrated',controlCalibration:cc,sampleCalibration:sc,requestedMaxIndel:s.maxIndel,contextEmission:context},
  mixture:{components:[...sizes].map(([size,fraction])=>({size,fraction})),wt:totalNoEdit,indel:1-totalNoEdit,inFrame:chosen.filter(c=>c.model.events.length&&c.model.size%3===0).reduce((sum,c)=>sum+c.weight,0),outOfFrame:chosen.filter(c=>c.model.size%3!==0).reduce((sum,c)=>sum+c.weight,0),r2:full.r2,rmse:full.rmse,upstreamError:1-sc.validationR2,bases:positions.length,fitStart:positions[0],fitEnd:positions.at(-1)!+1,converged:fit.converged,passesChecks:fit.converged&&finalTest.r2>=.9&&!attributionAmbiguous.length&&!boundary&&selectedGeometry.condition!==null&&selectedGeometry.condition<=100},
  pooled:{preflight:report,controlR2:controlScore.r2,sampleTrim:sampleReg.trim,controlTrim:ctlReg.trim,amplicons:perRef.map(({estimable,groups,...a})=>a)},
  joint:{groups,targets:targetResults,heldOutR2:finalTest.r2,selectedModels:chosen.length,candidateModels:models.length,iterations,maxDeletion:s.largeDeletions?(s.maxDeletion??200):undefined,ambiguousFraction:ambiguous.reduce((sum,c)=>sum+c.weight,0)}};
}
