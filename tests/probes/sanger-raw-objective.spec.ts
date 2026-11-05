import{SangerSequenceFitCache}from'../../src/app/workers/core/sanger-sequence-fit-cache';
import{rawIntervalScorer}from'../../src/app/workers/core/sanger-raw-interval-fit';
import '@angular/compiler';import loadHighs from'../../public/vendor/highs-1.15.4/highs.mjs';import{it}from'vitest';import{readFileSync,writeFileSync}from'node:fs';import{parseSangerTrace}from'../../src/app/utils/sanger-trace';import{pairedReferencePeaks}from'../../src/app/workers/core/sanger-signal';import{materialize}from'../../src/app/workers/core/sanger-path-common';
(process.env.CASSAGO_RAW_OBJECTIVE==='1'?it:it.skip)('audits raw ADC interval fitting after inference only',async()=>{
 const path=process.env.CASSAGO_RAW_INPUT!,template=JSON.parse(readFileSync(path+'/template.json','utf8')),truth=JSON.parse(readFileSync(path+'/truth.json','utf8')),result=JSON.parse(readFileSync(process.env.CASSAGO_RAW_RESULT||path+'/result.json','utf8')),trace=(name:string)=>{const b=readFileSync(path+'/'+name);return parseSangerTrace(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),name);},control=trace('control.ab1'),sample=trace('edit.ab1'),offset=template.primer.length,refs=template.amplicons.map((a:any)=>a.sequence.slice(offset)),ct=result.pooled.controlTrim,st=result.pooled.sampleTrim,peaks=pairedReferencePeaks(control,sample,refs,ct,st,Math.min(...template.amplicons.map((a:any)=>a.cuts[0]-offset))),optimizer=await loadHighs();
 const inferred=result.annotations.flatMap((a:any)=>a.components.map((c:any)=>materialize(a.reference,c.events.map((e:any)=>({start:e.start,end:e.end,inserted:e.insertedSequence})))!.slice(offset))),known=truth.truth.map((t:any)=>t.inputSequence.slice(offset));
 const fit=(seqs:string[],sigma:number)=>{
  const lines=['Minimize',' obj: v','Subject To',' flux: '+seqs.map((_:any,i:number)=>'x'+i).join(' + ')+' - '+refs.map((_:any,i:number)=>'c'+i).join(' - ')+' = 0'];let row=0;const compact=new Map<string,{terms:string;d:number;lo:number;hi:number}>();
  for(const [trace,dict,trim,bg,prefix]of [[control,refs,ct,peaks.control.background!,'c'],[sample,seqs,st,peaks.sample.background!,'x']]as const){
   const end=Math.min(trace.peaks.length+trim-3,...dict.map((s:string)=>s.length-22));
   for(let p=Math.max(30,trim+10);p<end;p++)for(let b=0;b<4;b++){
    const terms=dict.flatMap((s:string,i:number)=>s[p]==='ACGT'[b]?[prefix+i]:[]).join(' + ');
    for(let d=0;d<=4;d++){const value=trace.channels['ACGT'[b]as'A'][trace.peaks[p-trim]+d]-bg[b],key=prefix+'/'+terms+'/'+d,old=compact.get(key);if(old){old.lo=Math.max(old.lo,value-.5);old.hi=Math.min(old.hi,value+.5);}else compact.set(key,{terms,d,lo:value-.5,hi:value+.5});}
   }
  }
  for(const {terms,d,lo,hi}of compact.values()){const shape=Math.exp(-d*d/sigma),sum=terms?terms.split(' + ').map(t=>shape+' '+t).join(' + '):'';lines.push(' r'+row+++': '+(sum?sum+' + ':'')+'v >= '+lo,' r'+row+++': '+(sum?sum+' - v':'- v')+' <= '+hi);}
  lines.push('End');const r=optimizer.solve(lines.join('\n'),{output_flag:false,primal_feasibility_tolerance:1e-9,dual_feasibility_tolerance:1e-9});return{sigma,status:r.Status,violation:r.Columns.v?.Primal,amplitude:refs.reduce((s:number,_:any,i:number)=>s+(r.Columns['c'+i]?.Primal||0),0)};
 };
 const output:any={note:'Post-inference oracle diagnostic only; neither truth nor its dictionary is an inference input.',inferred:[],truth:[]};
 for(const sigma of (process.env.CASSAGO_RAW_GRID==='1'?[3.995,3.998,3.999,4,4.001,4.002,4.005]:[])){output.inferred.push(fit(inferred,sigma));output.truth.push(fit(known,sigma));}
 const positions=Array.from({length:Math.min(control.peaks.length+ct-3,sample.peaks.length+st-3,...refs.map((s:string)=>s.length-22))-30},(_,i)=>i+30),score=rawIntervalScorer(control,sample,refs,positions,ct,st,peaks.control.background!,peaks.sample.background!,(peaks as any).calibration,optimizer)!;
 const signal=positions.flatMap(p=>peaks.sample.vectors[p-st]),bounded=rawIntervalScorer(control,sample,refs,positions,ct,st,peaks.control.background!,peaks.sample.background!,(peaks as any).calibration,optimizer,signal)!;
 output.centroidBound=bounded?.maximumCentroidMse;const centroidTruth=new SangerSequenceFitCache(known,positions,signal).fit(known.map((_:any,i:number)=>i)),centroidInferred=new SangerSequenceFitCache(inferred,positions,signal).fit(inferred.map((_:any,i:number)=>i));output.centroidTruthMse=centroidTruth.rmse**2;output.centroidInferredMse=centroidInferred.rmse**2;
 output.learnedInferred=score(inferred,inferred.map(()=>1/inferred.length));output.learnedTruth=score(known,known.map(()=>1/known.length));
 writeFileSync(process.env.CASSAGO_RAW_OUTPUT!,JSON.stringify(output,null,2));
},120000);
