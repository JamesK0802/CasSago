import type {SangerTrace} from '../../models/sanger-analysis';
/** Define measurement coverage from the WT control only. Low edited-read
 * Phred is not a reason to discard a genuine mixed signal. A rolling Q30
 * control interval retains high-quality downstream measurements; independent
 * measured-signal validation still decides whether the inferred fit passes. */
export function controlUsableRegion(control:SangerTrace,trim:number,firstCut:number,lastCut:number,readStart:number,readEnd:number){
 const lo=Math.max(readStart,trim),hi=Math.min(readEnd,trim+control.peaks.length),quality=control.quality;
 if(!quality?.some(q=>q>0))return{start:lo,end:hi,qualityAvailable:false};
 const windows:{start:number;end:number}[]=[];let start=-1;
 for(let p=lo;p<=hi-10;p++){const mean=Array.from(quality.slice(p-trim,p-trim+10)).reduce((a,b)=>a+b,0)/10,good=mean>=30;if(good&&start<0)start=p;if((!good||p===hi-10)&&start>=0){windows.push({start,end:good?p+10:p+9});start=-1;}}
 const valid=windows.filter(w=>w.start<=firstCut-30&&w.end>=lastCut+20).sort((a,b)=>(b.end-b.start)-(a.end-a.start));
 if(!valid.length)throw Error('The WT control has no continuous reliable interval spanning the targets and flanks.');
 return{...valid[0],qualityAvailable:true};
}
