import{refineDiscretePhase,SangerLinearOptimizer}from'./sanger-discrete-phase';
import {proposeStratifiedPaths} from './sanger-stratified-paths';
import{sequenceGroupResult}from'./sanger-sequence-group-result';
import{fitLengthMarginals}from'./sanger-length-marginals';import{decodeFrequencyPaths}from'./sanger-frequency-path-decoder';
import {SangerFitCache} from './sanger-fit-cache';
import {normalizeMoleculeEdits} from './sanger-edit-normalize';
import {repairSignalTraceCycles} from './sanger-cycle';
import {SangerAnalysisResult,TRACE_BASES} from '../../models/sanger-analysis';
import {registerTrace} from './sanger-pooled';
import {contrastPeaks,integratedPeaks,pairedReferencePeaks,calibrateSignal,correctGain,gramCondition,rawSignalScore} from './sanger-signal';
import {emissionVector,fitContextEmission,fitNormalizedEmissionMixture as fitEmissionMixture} from './sanger-emission';
import {SignalPathInput,PathEdit,Model,materialize,prepare,assertSignalTraceGrid,eventsRespectIndelOrigin} from './sanger-path-common';
/** Invariant-flank length marginals, local frequency/path decoding and one
 * whole-molecule fluorescence refit. These are
 * best-fit cut-local hypotheses, not a posterior over physical haplotypes.
 * All measured positions enter discovery, so no independent held-out score is claimed. */
export function analyzeOffsetGroups(input:SignalPathInput,progress:(s:string)=>void=()=>{},audit?:(state:any)=>void):SangerAnalysisResult{return initializeOffsetProfiles(input,progress,audit).result;}
export function initializeOffsetProfiles(input:SignalPathInput,progress:(s:string)=>void=()=>{},audit?:(state:any)=>void,optimizer?:SangerLinearOptimizer){
 if(!input.ampliconConfirmed||!input.control||!input.sample)throw Error('Confirm submitted products and provide matched control/sample chromatograms.');
 const controlRepair=repairSignalTraceCycles(input.control),sampleRepair=repairSignalTraceCycles(input.sample);input={...input,control:controlRepair.trace,sample:sampleRepair.trace};
 assertSignalTraceGrid(input.control);assertSignalTraceGrid(input.sample);
 const {refs,report,anchor,maxDeletion,maxInsertion}=prepare(input.amplicons,input.primer,input.settings),s=input.settings;
 progress('Registering offset traces to the sequencing primer');const registered=refs.map(r=>({name:r.name,sequence:r.sequence,cut:Math.min(...r.cuts)}));
 const cr=registerTrace(input.control,registered,s,integratedPeaks(input.control,[0,0,0,0]).vectors),sr=registerTrace(input.sample,registered,s,integratedPeaks(input.sample,[0,0,0,0]).vectors),peaks=pairedReferencePeaks(input.control,input.sample,refs.map(r=>r.sequence),cr.trim,sr.trim,Math.min(...refs.flatMap(r=>r.cuts)),optimizer),ctl=peaks.control.vectors,raw=peaks.sample.vectors;
 const calibrate=(vectors:number[][],trim:number)=>{const obs=Array.from({length:vectors.length+trim},(_,p)=>vectors[p-trim]||[0,0,0,0]),rows=Array.from({length:Math.min(...refs.flatMap(r=>r.cuts))-12},(_,p)=>p).filter(p=>p>=Math.max(s.readStart,trim+10)&&p-trim<vectors.length-2);return calibrateSignal(obs,rows,refs.map(r=>Array.from(r.sequence,(_,p)=>TRACE_BASES.map(b=>Number(r.sequence[p]===b)))));};
 const cc=calibrate(ctl,cr.trim),sc=calibrate(raw,sr.trim);if(cc.validationR2<.8||sc.validationR2<.8)throw Error('Signal-path upstream calibration failed.');
 const end=Math.min(Math.max(...report.positions)+1,cr.trim+ctl.length-3,sr.trim+raw.length-3),positions=report.positions.filter(p=>p<end);if(positions.length<120)throw Error('Insufficient actual trace coverage for path inference.');
 const wt=ctl.map(v=>correctGain(v,cc.gains)),sample=raw.map(v=>correctGain(v,sc.gains)),controlObs=Array.from({length:wt.length+cr.trim},(_,p)=>wt[p-cr.trim]||[0,0,0,0]);
 const context=fitContextEmission(refs.map(r=>r.sequence),positions,controlObs,cc.weights,cc.lag),theta=context.theta;
 const emit=(sequence:string,p:number,lag:number)=>{const v=emissionVector(sequence,p,lag,theta);if(sequence[p]==='N')for(let b=0;b<4;b++)v[b]+=(1-Math.abs(lag))/4;return v;};
 const models:Model[]=refs.map((r,h)=>({reference:h,sequence:r.sequence,events:[],aliases:[{reference:h,events:[]}]})),cols=models.map(m=>positions.flatMap(p=>emit(m.sequence,p,sc.lag))),controlCols=models.map(m=>positions.flatMap(p=>emit(m.sequence,p,cc.lag)));
 const ctlFit=fitEmissionMixture(controlCols,positions.flatMap(p=>wt[p-cr.trim])),ctlScore=rawSignalScore(controlCols,positions.flatMap(p=>ctl[p-cr.trim]),ctlFit.weights,cc.gains);if(ctlScore.r2<.85)throw Error('WT product references do not explain the unedited control.');
 const y=positions.flatMap(p=>sample[p-sr.trim]),train=positions.map((_,i)=>i),tune=train,test=train;
 const rawY=positions.flatMap(p=>raw[p-sr.trim]),fitCache=new SangerFitCache(cols,y),fitIds=(ids:number[])=>fitCache.fit(ids);
 let ids=refs.map((_,i)=>i),fit=fitIds(ids),iterations=0;const alignments=0,discoveryLimited=false;

 if(refs.some(r=>r.cuts.length!==1)||context.used||Math.abs(sc.lag)>1e-12)throw Error('Offset decoder requires one target and unit cycle emissions.');
 progress('Estimating indel length fractions from invariant flanks');
 const lengths=fitLengthMarginals(refs,positions,y,(seq,p)=>emit(seq,p,sc.lag),maxDeletion,maxInsertion,s.editStartRadius??10),families=refs.flatMap((r,h)=>{const all=lengths.families.filter(f=>f.reference===h),share=all.reduce((sum,f)=>sum+f.fraction,0),keep=all.filter(f=>f.size!==0&&f.fraction>=share*.01);return[{reference:h,size:0,fraction:share-keep.reduce((sum,f)=>sum+f.fraction,0)},...keep];});
 progress('Decoding offset local edits');
 const profileOptions:PathEdit[][][][]=[],selected:number[]=[];
 const decoded:{groups:{reference:number;events:PathEdit[];fraction:number}[];details:any[]}={groups:[],details:[]};
 for(let h=0;h<refs.length;h++){
  const rows=positions.flatMap((p,i)=>refs.every((r,k)=>k===h||p<r.cuts[0]-3||p>=r.cuts[0]+13)?[i]:[]),ownPositions=rows.map(i=>positions[i]),ownObserved=rows.flatMap(i=>{
   const p=positions[i],v=y.slice(4*i,4*i+4);
   for(const f of families)if(f.reference!==h){const r=refs[f.reference],b='ACGT'.indexOf(r.sequence[p<r.cuts[0]-3?p:p-f.size]);if(b>=0)v[b]-=f.fraction;}
   return v;
  }),own=families.filter(f=>f.reference===h).map(f=>({...f,reference:0})),proposals=proposeStratifiedPaths(refs[h],own,ownPositions,ownObserved,lengths.rmse,5,200,4,stage=>progress('First target reference '+(h+1)+' / '+refs.length+' · '+stage));
  if(optimizer){
   const share=own.reduce((sum,f)=>sum+f.fraction,0),target=ownObserved.map(v=>v/share);
   for(let i=0;i<proposals.options.length;i++){
    progress(`Refining first-target allocations · reference ${h+1} / ${refs.length} · ${i+1} / ${proposals.options.length}`);
    const phased=refineDiscretePhase([refs[h]],[[[]]],[proposals.options[i].map(g=>g.events)],ownPositions,target,optimizer,()=>{},0,false);
    proposals.options[i]=phased.profiles[0].map(events=>({events,fraction:0}));
    const columns=proposals.options[i].map(g=>{const seq=materialize(refs[h].sequence,g.events)!;return ownPositions.flatMap(p=>'ACGT'.split('').map(b=>Number(seq[p]===b)));}),fit=new SangerFitCache(columns,target).fit(columns.map((_,j)=>j));
    proposals.options[i].forEach((g,j)=>g.fraction=share*fit.weights[j]);
    let maximum=0;for(const [row,p]of ownPositions.entries())if(p>=refs[h].cuts[0]-3&&p<refs[h].cuts[0]+13){let error=0;for(let b=0;b<4;b++){const predicted=columns.reduce((sum,col,j)=>sum+col[4*row+b]*fit.weights[j],0);error+=(predicted-target[4*row+b])**2;}maximum=Math.max(maximum,share*Math.sqrt(error/4));}
    proposals.scores[i]={rmse:share*fit.rmse,maximum};
   }
  }
  let best=-1;for(let i=0;i<proposals.options.length;i++)if(proposals.scores[i].maximum<=Math.max(1e-5,4*lengths.rmse)&&(best<0||proposals.scores[i].rmse<proposals.scores[best].rmse))best=i;
  if(best<0)throw Error('No validated first-target profile for '+refs[h].name+' '+JSON.stringify(proposals.details));
  profileOptions.push(proposals.options.map(option=>option.map(g=>g.events)));selected.push(best);
  decoded.groups.push(...proposals.options[best].map(g=>({...g,reference:h})));decoded.details.push(...proposals.details);
 }
 const signature=(sequence:string)=>sequence.slice(positions[0],positions.at(-1)!+1).padEnd(positions.length,'-');
 const seenEdits=refs.map(()=>new Set<string>()),signatures=new Map(models.map((m,id)=>[signature(m.sequence),id]));
  const addVariant=(h:number,events:PathEdit[])=>{events=normalizeMoleculeEdits(refs[h].sequence,events);const editKey=JSON.stringify(events);if(seenEdits[h].has(editKey))return;seenEdits[h].add(editKey);if(!eventsRespectIndelOrigin(events,refs[h].cuts,s.editStartRadius))return;const sequence=materialize(refs[h].sequence,events);if(!sequence)return;const key=signature(sequence),old=signatures.get(key),alias={reference:h,events};if(old!==undefined){if(!models[old].aliases.some(a=>a.reference===h&&JSON.stringify(a.events)===JSON.stringify(events)))models[old].aliases.push(alias);return;}const id=models.length;models.push({reference:h,sequence,events,aliases:[alias]});cols.push(positions.flatMap(p=>emit(sequence,p,sc.lag)));signatures.set(key,id);};

 const proposed:number[]=[];for(const group of decoded.groups){if(group.fraction<=1e-9)continue;addVariant(group.reference,group.events);const id=signatures.get(signature(materialize(refs[group.reference].sequence,group.events)!));if(id!==undefined)proposed.push(id);}
 ids=[...new Set([...ids,...proposed])];fit=fitIds(ids);
 // Equivalent deletion placements annotate the same molecule, not extra mass.
 for(const group of decoded.groups)for(const e of group.events)if(!e.inserted.length)for(let start=refs[group.reference].cuts[0]-3;start<=refs[group.reference].cuts[0]+3;start++)addVariant(group.reference,[{start,end:start+e.end-e.start,inserted:''}]);
 iterations=1;
 progress('Validating offset profiles');
 return {result:sequenceGroupResult({input,refs,positions,report,anchor,end,models,cols,ids,fit,ctlScore,cc,sc,context,cr,sr,rawY,controlRepair,sampleRepair,maxDeletion,maxInsertion,fitCache,y,audit,train,test,tune,emit,iterations}),profiles:profileOptions,selected};
}
