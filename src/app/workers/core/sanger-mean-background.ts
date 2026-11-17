import type {SangerTrace} from '../../models/sanger-analysis';
/** First moment of fluorescence where every submitted WT is quiet. Includes
 * rectified electronic/background noise: a cycle minimum estimates its lower
 * envelope, not its expected additive contribution. Unedited TRAIN rows only. */
export function meanQuietBackground(trace:SangerTrace,sequences:string[],trim:number,upstream:number[],centerWeight=0){
 const rows=Array.from({length:4},()=>[] as number[]);
 for(const p of upstream.filter((_,i)=>i%3!==0)){const i=p-trim;if(i<1||i>=trace.peaks.length-1)continue;
  const a=Math.max(0,(trace.peaks[i-1]+trace.peaks[i])/2),z=Math.min(trace.channels.A.length-1,(trace.peaks[i]+trace.peaks[i+1])/2);if(z<=a)continue;
  [...'ACGT'].forEach((base,b)=>{if(!sequences.every(s=>[-1,0,1].every(d=>/[ACGT]/.test(s[p+d]||'')&&s[p+d]!==base)))return;const y=trace.channels[base as keyof SangerTrace['channels']];let area=0;for(let j=Math.floor(a);j<z;j++){const lo=Math.max(a,j),hi=Math.min(z,j+1);area+=(y[j]+(y[j+1]-y[j])*(lo-j)+y[j]+(y[j+1]-y[j])*(hi-j))*(hi-lo)/2;}const mean = area / (z - a); rows[b].push(centerWeight === 0 ? mean : (1 - centerWeight) * mean + centerWeight * y[trace.peaks[i]]);});
 }
 const means=rows.map(v=>v.length>=8?v.reduce((a,b)=>a+b,0)/v.length:0);
 return{background:means,variance:rows.map((v,b)=>v.length>=8?v.reduce((n,x)=>n+(x-means[b])**2,0)/Math.max(1,v.length-1):0),samples:rows.map(v=>v.length)};
}
