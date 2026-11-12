/** Integrate the entire base cycle between adjacent PLOC midpoints.
 * The same photon-area observation is used for measured and constructed reads. */
import type {SangerTrace} from '../../models/sanger-analysis';
export function cycleFluorescence(trace:SangerTrace,background:number[]){
 const areas:number[][]=[],strength:number[]=[],vectors:number[][]=[];
 for(let i=0;i<trace.peaks.length;i++){const p=trace.peaks[i],left=i?trace.peaks[i-1]:p-12,right=i+1<trace.peaks.length?trace.peaks[i+1]:p+12,a=Math.max(0,(p+left)/2),b=Math.min(trace.channels.A.length-1,(p+right)/2);
  const values=['A','C','G','T'].map((base,k)=>{const y=trace.channels[base as keyof SangerTrace['channels']];let total=0;for(let j=Math.floor(a);j<b;j++){const lo=Math.max(a,j),hi=Math.min(b,j+1),yl=y[j]+(y[j+1]-y[j])*(lo-j)-background[k],yh=y[j]+(y[j+1]-y[j])*(hi-j)-background[k];total+=(Math.max(0,yl)+Math.max(0,yh))*(hi-lo)/2;}return total;}),sum=values.reduce((n,v)=>n+v,0);areas.push(values);strength.push(sum);vectors.push(values.map(v=>sum?v/sum:0));
 }
 return{areas,strength,vectors};
}
