import '@angular/compiler';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import{it}from'vitest';import{readFileSync,writeFileSync}from'node:fs';
import{parseSangerTrace}from'../../src/app/utils/sanger-trace';
import{automaticSangerSettings}from'../../src/app/workers/core/sanger-automatic';
import{prepare,materialize}from'../../src/app/workers/core/sanger-path-common';
import{pairedReferencePeaks,correctGain}from'../../src/app/workers/core/sanger-signal';
import{SangerSequenceFitCache}from'../../src/app/workers/core/sanger-sequence-fit-cache';
(process.env.CASSAGO_OBJECTIVE_DIAGNOSTIC==='1'?it:it.skip)('compares post-inference truth and inferred dictionaries on the same objective; never used for inference',async()=>{
 const path=process.env.CASSAGO_OBJECTIVE_INPUT!,template=JSON.parse(readFileSync(path+'/template.json','utf8')),truth=JSON.parse(readFileSync(path+'/truth.json','utf8')),result=JSON.parse(readFileSync(path+'/result.json','utf8'));
 const trace=(f:string)=>{const b=readFileSync(path+'/'+f);return parseSangerTrace(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),f);},control=trace('control.ab1'),sample=trace('edit.ab1'),input:any={references:template.amplicons.map((a:any)=>({...a,cut:a.cuts[0]})),primer:template.primer,control,sample,maxDeletion:10,maxInsertion:10},settings={...automaticSangerSettings(input),readEnd:1000,maxIndel:10},prepared=prepare(input.references,input.primer,settings),ct=result.pooled.controlTrim,st=result.pooled.sampleTrim,peaks=pairedReferencePeaks(control,sample,prepared.refs.map(r=>r.sequence),ct,st,Math.min(...prepared.refs.flatMap(r=>r.cuts))),end=Math.min(settings.readEnd,ct+peaks.control.vectors.length-3,st+peaks.sample.vectors.length-3,...prepared.refs.map(r=>r.sequence.length-22)),positions=Array.from({length:end-prepared.anchor},(_,i)=>i+prepared.anchor),observed=positions.flatMap(p=>correctGain(peaks.sample.vectors[p-st],result.diagnostics.sampleCalibration.gains));
 const inferred=result.annotations.flatMap((a:any)=>a.components.map((c:any)=>materialize(a.reference,c.events.map((e:any)=>({start:e.start,end:e.end,inserted:e.insertedSequence})))!.slice(template.primer.length))),known=truth.truth.map((t:any)=>t.inputSequence.slice(template.primer.length));
 const fit=(sequences:string[],signal=observed)=>{const cache=new SangerSequenceFitCache(sequences,positions,signal),r=cache.fit(sequences.map((_,i)=>i));return{rmse:r.rmse,r2:r.r2,converged:r.converged,gap:r.optimalityGap};};
 const calibration=(peaks as any).calibration,response=Array.from({length:9},(_,i)=>calibration.amplitude*Math.exp(-((i-4)**2)/calibration.sigma)),norm=response.reduce((s:number,x:number)=>s+x*x,0),rawObserved=positions.flatMap(p=>[...'ACGT'].map((b,j)=>response.reduce((s:number,x:number,i:number)=>s+x*(sample.channels[b as 'A'][sample.peaks[p-st]+i-4]-peaks.sample.background![j]),0)/norm));
 const optimizer=await loadHighs();
 const intervalFit=(sequences:string[])=>{
  const lines=['Minimize',' obj: v','Subject To',' total: '+sequences.map((_,i)=>'x'+i).join(' + ')+' = 1'];let row=0;
  for(const p of positions)for(let b=0;b<4;b++){
   const sum=sequences.flatMap((seq,i)=>seq[p]==='ACGT'[b]?['x'+i]:[]).join(' + ');if(!sum)continue;
   let lo=0,hi=1;for(let k=0;k<9;k++){const value=sample.channels['ACGT'[b] as 'A'][sample.peaks[p-st]+k-4]-peaks.sample.background![b];lo=Math.max(lo,(value-.5)/response[k]);hi=Math.min(hi,(value+.5)/response[k]);}
   lines.push(' r'+row+++': '+sum+' + v >= '+lo,' r'+row+++': '+sum+' - v <= '+hi);
  }
  lines.push('End');const r=optimizer.solve(lines.join('\n'),{output_flag:false});return{status:r.Status,maxViolation:r.Columns.v?.Primal};
 };
 writeFileSync(process.env.CASSAGO_OBJECTIVE_OUTPUT!,JSON.stringify({inferred:fit(inferred),truth:fit(known),calibration,intervalInferred:intervalFit(inferred),intervalTruth:intervalFit(known),rawInferred:fit(inferred,rawObserved),rawTruth:fit(known,rawObserved),note:'Oracle diagnostic only: truth read after blind result; no candidates are returned to inference.'},null,2));
});
