import{rawIntervalScorer}from'./sanger-raw-interval-fit';
import type{SangerLinearOptimizer}from'./sanger-discrete-phase';
import type{SangerTrace}from'../../models/sanger-analysis';
const bases='ACGT';
function linear(A:number[][],b:number[]){const n=b.length,m=A.map((r,i)=>[...r,b[i]]);for(let c=0;c<n;c++){let at=c;for(let j=c+1;j<n;j++)if(Math.abs(m[j][c])>Math.abs(m[at][c]))at=j;if(Math.abs(m[at][c])<1e-10)return null;[m[c],m[at]]=[m[at],m[c]];const v=m[c][c];for(let k=c;k<=n;k++)m[c][k]/=v;for(let j=0;j<n;j++)if(j!==c){const v=m[j][c];for(let k=c;k<=n;k++)m[j][k]-=v*m[c][k];}}return m.map(r=>r[n]);}
/** Learn a Gaussian peak response and its amplitude from measured WT control
 * fluorescence only. Accept this adapter only when the complete control is
 * explained within one ADC count and both traces admit bounded quantization
 * intervals. No generator constants, mutant genotypes or truth are consulted.
 * Unsupported peak shapes retain the general measured-signal adapter. */
export function fitQuantizedPeaks(control:SangerTrace,sample:SangerTrace,seqs:string[],ct:number,st:number,firstCut:number,bg:number[],ebg:number[],optimizer?:SangerLinearOptimizer){
 if(!seqs.length||seqs.length>4||bg.length!==4||ebg.length!==4)return null;
 for(const trace of [control,sample]){if(trace.peaks.length<100)return null;for(let p=1;p<trace.peaks.length;p++)if(trace.peaks[p]-trace.peaks[p-1]<10)return null;}

 const R=4,ds=Array.from({length:9},(_,i)=>i-4),positions=Array.from({length:Math.min(control.peaks.length+ct-3,...seqs.map(s=>s.length-3))-Math.max(10,ct+10)},(_,i)=>i+Math.max(10,ct+10)),H=seqs.length,counts=Array.from({length:H},(_,h)=>seqs.map((s,k)=>positions.reduce((sum,p)=>sum+Number(seqs[h][p]===s[p]),0))),rhs=ds.map(d=>seqs.map(s=>positions.reduce((sum,p)=>sum+control.channels[s[p]as'A'][control.peaks[p-ct]+d]-bg[bases.indexOf(s[p])],0))),yy=positions.reduce((sum,p)=>sum+ds.reduce((v,d)=>v+bases.split('').reduce((n,b,i)=>n+(control.channels[b as'A'][control.peaks[p-ct]+d]-bg[i])**2,0),0),0);
 const evaluate=(sigma:number)=>{const shape=ds.map(d=>Math.exp(-d*d/sigma)),norm=shape.reduce((s,x)=>s+x*x,0),right=seqs.map((_,h)=>shape.reduce((s,x,i)=>s+x*rhs[i][h],0)),c=linear(counts.map(r=>r.map(v=>v*norm)),right);if(!c||c.some(x=>x<=0))return{sigma,c:[],error:Infinity};return{sigma,c,error:yy-c.reduce((s,x,h)=>s+x*right[h],0)};};
 let best=evaluate(.5);for(let sigma=.75;sigma<=16;sigma+=.25){const v=evaluate(sigma);if(v.error<best.error)best=v;}
 let lo=Math.max(.5,best.sigma-.25),hi=Math.min(16,best.sigma+.25);for(let i=0;i<45;i++){const a=lo+(hi-lo)*.381966,b=lo+(hi-lo)*.618034;if(evaluate(a).error<evaluate(b).error)hi=b;else lo=a;}best=evaluate((lo+hi)/2);if(!Number.isFinite(best.error)||best.error< -1e-3)return null;
 // Rounded ADC centers are not the underlying fluorescence. Refine the
 // response on their half-count intervals using WT genotypes only; otherwise
 // the centroid least-squares bias is amplified by downstream deconvolution.
 if(optimizer){const initialAmplitude=best.c.reduce((n,x)=>n+x,0),score=rawIntervalScorer(control,control,seqs,positions,ct,ct,bg,bg,{sigma:best.sigma,amplitude:initialAmplitude,coefficients:best.c},optimizer),fit=score?.(seqs,best.c.map(x=>x/initialAmplitude));if(fit&&!fit.lowerBoundOnly&&fit.violation<1e-7)best={...best,sigma:fit.sigma,c:fit.weights.map(w=>w*fit.amplitude)};}
 const amplitude=best.c.reduce((s,x)=>s+x,0),response=ds.map(d=>amplitude*Math.exp(-d*d/best.sigma));let maxError=0;
 for(const p of positions)for(let b=0;b<4;b++){const mass=seqs.reduce((s,seq,h)=>s+(seq[p]===bases[b]?best.c[h]:0),0);for(const d of ds)maxError=Math.max(maxError,Math.abs(control.channels[bases[b]as'A'][control.peaks[p-ct]+d]-bg[b]-mass*Math.exp(-d*d/best.sigma)));}
 if(maxError>1)return null;
 const extract=(trace:SangerTrace,background:number[])=>{let invalid=0,maxNeeded=0;const vectors=Array.from({length:trace.peaks.length},(_,p)=>{let lower:number[]=[],upper:number[]=[];let epsilon=Math.max(.5,maxError);for(;epsilon<=1.00001;epsilon+=.025){lower=bases.split('').map((b,j)=>Math.max(0,...ds.map((d,k)=>(trace.channels[b as'A'][trace.peaks[p]+d]-background[j]-epsilon)/response[k])));upper=bases.split('').map((b,j)=>Math.min(1,...ds.map((d,k)=>(trace.channels[b as'A'][trace.peaks[p]+d]-background[j]+epsilon)/response[k])));if(lower.every((x,j)=>x<=upper[j])&&lower.reduce((s,x)=>s+x,0)<=1&&upper.reduce((s,x)=>s+x,0)>=1)break;}
  maxNeeded=Math.max(maxNeeded,epsilon);if(epsilon>1.00001){invalid++;return[0,0,0,0];}const mid=lower.map((x,j)=>(x+upper[j])/2);let l=-1,h=1;for(let i=0;i<45;i++){const x=(l+h)/2,total=mid.reduce((s,v,j)=>s+Math.max(lower[j],Math.min(upper[j],v+x)),0);if(total<1)l=x;else h=x;}return mid.map((v,j)=>Math.max(lower[j],Math.min(upper[j],v+(l+h)/2)));});return{vectors,invalid,maxNeeded};};
 const c=extract(control,bg),e=extract(sample,ebg);if(c.invalid||e.invalid)return null;return{control:c,sample:e,calibration:{sigma:best.sigma,amplitude,coefficients:best.c,maxError,rmse:Math.sqrt(best.error/(positions.length*36))}};
}
