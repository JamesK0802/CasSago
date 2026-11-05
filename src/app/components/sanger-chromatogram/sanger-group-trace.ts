import type{SangerAnalysisResult,TracePlot}from'../../models/sanger-analysis';
import type{SangerAnnotationRow}from'../../utils/sanger-annotation';
import{reverseComplement}from'../../utils/biology.utils';
import{pathColumn}from'../../workers/core/sanger-junction-paths';
import type{JunctionModel,JunctionReference}from'../../workers/core/sanger-junction-mixture';
/** Read-only rendering of already fitted models. No optimization or new calls. */
export function fittedGroupTrace(result:SangerAnalysisResult,index:number,row:SangerAnnotationRow,span:number):{plot:TracePlot|null;note:string}{
 const diagnostic=result.diagnostics,fit=diagnostic?.junctionFit,annotation=result.annotations?.[index],response=diagnostic?.controlResponse;
 if(!fit||!annotation||!response)return{plot:null,note:'A calibrated per-group signal is not available in this result.'};
 const meta=result.pooled?.preflight.amplicons[index],length=annotation.reference.length,oriented=annotation.direction==='reverse'?reverseComplement(annotation.reference):annotation.reference;
 let origin:number,cut:number;
 if(meta){origin=meta.direction==='forward'?meta.primerEnd:length-meta.primerStart+1;cut=(meta.direction==='reverse'?length-annotation.cut:annotation.cut)-origin;}
 else{const controlCut=result.plots[0]?.cutIndex,sampleCut=result.plots[1]?.cutIndex;if(controlCut===undefined||sampleCut===undefined)return{plot:null,note:'Registered read coordinates are unavailable.'};const genomicCut=annotation.direction==='reverse'?length-annotation.cut:annotation.cut,referenceOrigin=genomicCut-controlCut,prefix=Math.min(Math.max(0,referenceOrigin),Math.max(0,sampleCut-controlCut));origin=referenceOrigin-prefix;cut=controlCut+prefix;}
 const sequence=origin>=0?oriented.slice(origin):'N'.repeat(-origin)+oriented,ref:JunctionReference={sequence,cut,contextResponse:response},pattern=(m:JunctionModel)=>m.sizes?.map((v,t)=>`T${t+1}: ${v>0?'+':''}${v} bp`).join(' · ');
 const models=fit.models.filter(m=>m.reference===index&&m.size===row.size&&(!row.patternLabel||pattern(m as JunctionModel)===row.patternLabel)&&m.fraction>0),mass=models.reduce((s,m)=>s+m.fraction,0);
 if(!mass)return{plot:null,note:'No retained model column matches this group.'};
 const lo=Math.max(fit.positions[0]||0,Math.floor(cut-span/2)),hi=Math.min(fit.positions.at(-1)??cut+span/2,Math.ceil(cut+span/2)),positions=fit.positions.filter(p=>p>=lo&&p<=hi),values=positions.map(()=>[0,0,0,0]);
 if(!positions.length)return{plot:null,note:'No measured model cycles are available in this window.'};
 if(row.sequenceState==='unresolved')return{plot:null,note:'No unique sequence is available for this remainder.'};
 if(row.representative){
  const events=row.representative.events.map(e=>({target:e.target,start:(annotation.direction==='reverse'?length-e.end:e.start)-origin,removed:e.end-e.start,size:e.inserted-(e.end-e.start),inserted:annotation.direction==='reverse'?reverseComplement(e.insertedSequence||''):e.insertedSequence||''})).sort((a,b)=>a.start-b.start);
  if(events.some(e=>e.start<0))return{plot:null,note:'The candidate junction is upstream of the sequenced span.'};
  const path=events.flatMap(e=>[...(e.inserted?[{target:e.target,size:e.inserted.length,origins:[{start:e.start,weight:1,probabilities:Array.from(e.inserted,base=>Array.from('ACGT',b=>Number(b===base)))}]}]:[]),...(e.removed?[{target:e.target,size:-e.removed,origins:[{start:e.start,weight:1,probabilities:undefined}]}]:[])]);
  const column=pathColumn(ref,{reference:index,size:row.size,start:events[0]?.start||0,inserted:'',path},positions,diagnostic?.sampleCalibration?.lag||0,[]);positions.forEach((_,p)=>{for(let b=0;b<4;b++)values[p][b]=column[4*p+b];});
 }else{
 for(const model of models){const m=model as JunctionModel&{fraction:number};let path=m.path;
  if(!path)path=m.size?[{target:1,size:m.size,origins:(m.origins?.length?m.origins:[{start:m.start,weight:1,probabilities:m.probabilities}]).map(o=>({...o,probabilities:o.probabilities||(m.size>0?Array.from({length:m.size},(_,j)=>Array.from('ACGT',b=>m.inserted[j]&&m.inserted[j]!=='N'?Number(b===m.inserted[j]):.25)):undefined)}))}]:[];
  const column=pathColumn(ref,{...m,path},positions,diagnostic?.sampleCalibration?.lag||0,[]);positions.forEach((_,p)=>{for(let b=0;b<4;b++)values[p][b]+=column[4*p+b]*m.fraction/mass;});
 }
 }
 const peak=Math.max(1e-12,...values.flat()),width=Math.max(200,(hi-lo+1)*10+32),paths={}as TracePlot['paths'];for(const[b,base]of[...'ACGT'].entries())paths[base as keyof typeof paths]=positions.map((p,j)=>`${j?'L':'M'}${16+(p-lo)*10},${(124-98*values[j][b]/peak).toFixed(2)}`).join(' ');
 return{plot:{name:'Fitted group signal',start:lo-cut,end:hi-cut,width,cutX:16+(cut-lo)*10,paths,letters:[]},note:'Conditional sequence-candidate model · cycle-integrated signal. Shape is normalized for viewing; this is not a separately measured chromatogram.'};
}
