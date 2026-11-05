import {SangerTrace,TRACE_BASES} from '../../models/sanger-analysis';
/** Repair isolated PLOC metadata errors only when an actual raw fluorescence
 * peak supports the missing cycle, or a spurious call lies in a raw trough.
 * No reference, guide, mutation hypothesis or fitted fraction is consulted. */
export function repairSignalTraceCycles(trace:SangerTrace):{trace:SangerTrace;added:number;removed:number}{
 const peaks=Array.from(trace.peaks);let added=0,removed=0;
 const strength=(p:number)=>TRACE_BASES.reduce((sum,b)=>sum+(trace.channels[b][p]||0),0);
 for(let i=9;i<peaks.length-9;i++){
  const spacings=peaks.slice(i-8,i+9).slice(1).map((p,k)=>p-peaks[i-8+k]).sort((a,b)=>a-b),median=spacings[8],gap=peaks[i]-peaks[i-1];if(!median)continue;
  if(gap>1.8*median&&gap<2.2*median){
   const center=(peaks[i]+peaks[i-1])/2;let p=Math.round(center);
   for(let j=Math.ceil(center-.18*median);j<=Math.floor(center+.18*median);j++)if(strength(j)>strength(p))p=j;
   const neighbours=Math.min(strength(peaks[i-1]),strength(peaks[i])),trough=Math.max(strength(Math.round((peaks[i-1]+p)/2)),strength(Math.round((p+peaks[i])/2)));
   if(strength(p)>.25*neighbours&&trough<.65*strength(p)){peaks.splice(i,0,p);added++;i++;}
  }else if(gap<.6*median&&peaks[i+1]-peaks[i]<.6*median&&Math.abs(peaks[i+1]-peaks[i-1]-median)<.22*median&&strength(peaks[i])<.45*Math.min(strength(peaks[i-1]),strength(peaks[i+1]))){peaks.splice(i,1);removed++;i--;}
  if(added+removed>5)throw Error(trace.name+': more than five cycle repairs required; re-basecall the trace.');
 }
 if(!added&&!removed)return {trace,added,removed};
 const sequence=peaks.map(p=>{const values=TRACE_BASES.map(b=>trace.channels[b][p]||0);return TRACE_BASES[values.indexOf(Math.max(...values))];}).join('');
 return {trace:{...trace,peaks:Uint32Array.from(peaks),sequence,quality:null},added,removed};
}
