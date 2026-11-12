import type{SangerLinearOptimizer}from'./sanger-discrete-phase';
import{fitQuantizedPeaks}from'./sanger-quantized-peaks';
import {SangerTrace,TRACE_BASES} from '../../models/sanger-analysis';
export interface PeakSignal {vectors:number[][];strength:number[];noise:number[];areas?:number[][];background?:number[];backgroundSamples?:number[]}
const median=(a:number[])=>{if(!a.length)return 0;const s=[...a].sort((a,b)=>a-b);return s[Math.floor(s.length/2)];};
/** Integrate peaks relative to cycle spacing and subtract a local inter-peak baseline.
 * No sequence identity, allele fraction or fitted variant is used. */
export function integratedPeaks(trace:SangerTrace,fixedBackground?:number[]):PeakSignal{
 const n=trace.peaks.length,mid=TRACE_BASES.map(b=>Array.from({length:Math.max(0,n-1)},(_,i)=>trace.channels[b][Math.round((trace.peaks[i]+trace.peaks[i+1])/2)]));
 const strength:number[]=[],noise:number[]=[],vectors:number[][]=[],areas:number[][]=[];
 for(let i=0;i<n;i++){
  const spacing=Math.min(i?trace.peaks[i]-trace.peaks[i-1]:12,i+1<n?trace.peaks[i+1]-trace.peaks[i]:12),radius=Math.max(2,Math.min(5,Math.floor(spacing*.34)));
  const local=mid.map(values=>values.slice(Math.max(0,i-8),Math.min(values.length,i+8))),baseline=fixedBackground||local.map(median),spread=local.map((v,b)=>1.4826*median(v.map(x=>Math.abs(x-baseline[b]))));
  const values=TRACE_BASES.map((b,j)=>{let v=0,count=0;for(let p=Math.max(0,trace.peaks[i]-radius);p<=Math.min(trace.channels[b].length-1,trace.peaks[i]+radius);p++){v+=Math.max(0,trace.channels[b][p]-baseline[j]);count++;}return v/Math.max(1,count);}),total=values.reduce((a,b)=>a+b,0);
  areas.push(values);vectors.push(values.map(v=>total?v/total:0));strength.push(total);noise.push(spread.reduce((a,b)=>a+b,0));
 }return {vectors,strength,noise,areas};
}
/** Learn background only where every supplied WT sequence has no emitting base
 * in this channel at the current or adjacent cycle. Uses upstream TRAIN rows only,
 * never edited fractions/residuals. Sample the quiet-cycle lower envelope,
 * not peak-centre cross-colour photons. A midpoint tail alone is not a
 * background oracle; the entire cycle must be WT-quiet in this channel. */
export function referenceBackgroundPeaks(trace:SangerTrace,sequences:string[],trim:number,positions:number[],minimumSamples=8):PeakSignal{
 const training=positions.filter((_,i)=>i%3!==0),samples=TRACE_BASES.map(()=>[] as number[]);
 for(const p of training){const at=p-trim;if(at<1||at>=trace.peaks.length-1)continue;
  TRACE_BASES.forEach((base,b)=>{if(sequences.every(s=>[-1,0,1].every(d=>/[ACGT]/.test(s[p+d]||'')&&s[p+d]!==base))){const channel=trace.channels[base],first=Math.max(0,Math.ceil((trace.peaks[at-1]+trace.peaks[at])/2)),last=Math.min(channel.length-1,Math.floor((trace.peaks[at]+trace.peaks[at+1])/2));let floor=Infinity;for(let j=first;j<=last;j++)floor=Math.min(floor,channel[j]);if(Number.isFinite(floor))samples[b].push(Math.max(0,floor));}});
 }
 const background=samples.map(v=>v.length>=minimumSamples?median(v):0),corrected=integratedPeaks(trace,background);
 return {...corrected,background,backgroundSamples:samples.map(v=>v.length)};
}
/** Estimate each trace's baseline only in channels where the submitted WT
 * sequences cannot emit. All control positions are known WT; sample positions
 * are restricted to the unedited prefix. No fitted genotype enters this step. */
export function pairedReferencePeaks(control:SangerTrace,sample:SangerTrace,sequences:string[],controlTrim:number,sampleTrim:number,firstCut:number,optimizer?:SangerLinearOptimizer){
 const end=Math.min(...sequences.map(s=>s.length-2),control.peaks.length+controlTrim-2),known=Array.from({length:Math.max(0,end-Math.max(10,controlTrim+10))},(_,i)=>i+Math.max(10,controlTrim+10)),start=Math.max(10,sampleTrim+10),prefix=Array.from({length:Math.max(0,Math.min(firstCut-13,sample.peaks.length+sampleTrim-2)-start)},(_,i)=>i+start),c=referenceBackgroundPeaks(control,sequences,controlTrim,known),e=referenceBackgroundPeaks(sample,sequences,sampleTrim,prefix,3);
 if(c.backgroundSamples!.some(n=>n<8)||e.backgroundSamples!.some(n=>n<3))return{control:contrastPeaks(control),sample:contrastPeaks(sample),method:'contrast' as const};
 const quantized=fitQuantizedPeaks(control,sample,sequences,controlTrim,sampleTrim,firstCut,c.background!,e.background!,optimizer);
 if(quantized)return{control:{...c,vectors:quantized.control.vectors},sample:{...e,vectors:quantized.sample.vectors},method:'quantized-peak-shape' as const,calibration:quantized.calibration};
 return{control:c,sample:e,method:'reference-background' as const};
}
/** Per-cycle fluorescence contrast. Regress each dye on the summed peak shape,
 * with an intercept, linear background trend and peak derivative as nuisance.
 * Slopes recover channel amplitudes without subtracting inter-peak tails. */
export function contrastPeaks(trace:SangerTrace):PeakSignal{
 const vectors:number[][]=[],strength:number[]=[],noise:number[]=[];
 for(let i=0;i<trace.peaks.length;i++){
  const spacing=Math.min(i?trace.peaks[i]-trace.peaks[i-1]:12,i+1<trace.peaks.length?trace.peaks[i+1]-trace.peaks[i]:12),radius=Math.max(3,Math.min(6,Math.floor(spacing*.4))),center=trace.peaks[i],at=Array.from({length:2*radius+1},(_,k)=>center-radius+k).filter(p=>p>0&&p+1<trace.channels.A.length),n=at.length,coords=at.map(p=>p-center),ys=TRACE_BASES.map(b=>at.map(p=>trace.channels[b][p])),total=at.map((_,k)=>ys.reduce((s,c)=>s+c[k],0)),derivative=at.map(p=>TRACE_BASES.reduce((s,b)=>s+(trace.channels[b][p+1]-trace.channels[b][p-1])/2,0));
  const detrend=(v:number[])=>{const mt=coords.reduce((a,b)=>a+b,0)/n,m=v.reduce((a,b)=>a+b,0)/n,den=coords.reduce((s,t)=>s+(t-mt)**2,0),slope=den?coords.reduce((s,t,k)=>s+(t-mt)*(v[k]-m),0)/den:0;return v.map((x,k)=>x-m-slope*(coords[k]-mt));};
  const x=detrend(total),d=detrend(derivative),xx=x.reduce((s,v)=>s+v*v,0),dd=d.reduce((s,v)=>s+v*v,0),xd=x.reduce((s,v,k)=>s+v*d[k],0),det=xx*dd-xd*xd;
  const slopes=ys.map(v=>{const y=detrend(v),xy=x.reduce((s,v,k)=>s+v*y[k],0),dy=d.reduce((s,v,k)=>s+v*y[k],0);return Math.max(0,det>1e-8?(xy*dd-dy*xd)/det:xx?xy/xx:0);}),sum=slopes.reduce((a,b)=>a+b,0),v=slopes.map(x=>sum?x/sum:0);vectors.push(v);const amp=Math.max(...total)-Math.min(...total);strength.push(amp);
  noise.push(Math.sqrt(ys.reduce((s,y,b)=>{const yy=detrend(y),xy=x.reduce((s,v,k)=>s+v*yy[k],0),dy=d.reduce((s,v,k)=>s+v*yy[k],0),a=det>1e-8?(xy*dd-dy*xd)/det:xx?xy/xx:0,beta=det>1e-8?(dy*xx-xy*xd)/det:0;return s+yy.reduce((s,y,k)=>s+(y-a*x[k]-beta*d[k])**2,0);},0)/Math.max(1,4*(n-4))));
 }return {vectors,strength,noise};
}
/** Position confidence is based on peak signal/background only, not edited residuals. */
export function signalWeights(signals:PeakSignal[],indices:number[][]):number[]{
 return indices[0].map((_,j)=>Math.max(.25,Math.min(...signals.map((s,i)=>{const at=indices[i][j],snr=s.strength[at]/Math.max(1,s.noise[at]);return Math.min(1,snr/30);})))) ;
}
export interface SimplexFit {weights:number[];r2:number;rmse:number;converged:boolean;optimalityGap?:number}
/** Nonnegative, sum-to-one least squares. Scored on ALL unweighted signal rows. */
function denseSimplexFit(columns:number[][],observed:number[],positionWeights?:number[],initialWeights?:number[]):SimplexFit{
 const p=columns.length,n=observed.length;if(!p||!n||columns.some(c=>c.length!==n))throw new Error('Invalid trace regression dimensions.');
 // Subtracting a common per-row vector leaves A*w-y unchanged on the
 // simplex. It removes shared WT signal from the Lipschitz bound rather than
 // letting that uninformative common mode slow homeolog attribution.
 const center=observed.map((_,k)=>columns.reduce((sum,c)=>sum+c[k],0)/p);
 const x=columns.map(c=>c.map((v,k)=>v-center[k])),target=observed.map((v,k)=>v-center[k]);
 const gram=Array.from({length:p},()=>new Float64Array(p)),rhs=new Float64Array(p),weight=(k:number)=>positionWeights?.[Math.floor(k/4)]??1;
 for(let a=0;a<p;a++){for(let k=0;k<n;k++)rhs[a]+=x[a][k]*target[k]*weight(k);for(let b=0;b<=a;b++){let v=0;for(let k=0;k<n;k++)v+=x[a][k]*x[b][k]*weight(k);gram[a][b]=gram[b][a]=v;}}
 const L=Math.max(...gram.map(row=>row.reduce((sum,v)=>sum+Math.abs(v),0)));
 if(!L){const w=columns.map(()=>1/p);return {...scoreMixture(columns,observed,w),weights:w,converged:true,optimalityGap:0};}
 const gap=(w:number[])=>{const g=Array.from({length:p},(_,a)=>gram[a].reduce((sum,v,b)=>sum+v*w[b],0)-rhs[a]);return Math.max(0,w.reduce((sum,v,a)=>sum+v*g[a],0)-Math.min(...g))/n;};
 const initial=initialWeights?.length===p&&initialWeights.every(v=>Number.isFinite(v)&&v>=0)&&initialWeights.some(v=>v>0)?initialWeights:undefined,total=initial?.reduce((s,v)=>s+v,0)||1;
 let w=initial?initial.map(v=>v/total):Array(p).fill(1/p),converged=false;
 if(gap(w)<=1e-10)return {...scoreMixture(columns,observed,w),weights:w,converged:true,optimalityGap:gap(w)};
 // Pairwise exact line minimization on the simplex. Incremental gradients
 // avoid a dense matrix-vector product at every step in nearly aliased fits.
 let gradient=Array.from({length:p},(_,a)=>gram[a].reduce((sum,v,b)=>sum+v*w[b],0)-rhs[a]);
 for(let step=0;step<200000;step++){
  if(step%100===0)gradient=Array.from({length:p},(_,a)=>gram[a].reduce((sum,v,b)=>sum+v*w[b],0)-rhs[a]);
  let enter=0,leave=-1,lambda=0;for(let a=0;a<p;a++){lambda+=w[a]*gradient[a];if(gradient[a]<gradient[enter])enter=a;if(w[a]>0&&(leave<0||gradient[a]>gradient[leave]))leave=a;}
  if((lambda-gradient[enter])/n<=1e-10){converged=true;break;}
  if(leave<0||leave===enter)break;
  const curvature=gram[leave][leave]+gram[enter][enter]-2*gram[leave][enter],amount=Math.min(w[leave],curvature>0?(gradient[leave]-gradient[enter])/curvature:w[leave]);
  if(!(amount>0))break;w[leave]-=amount;w[enter]+=amount;
  for(let a=0;a<p;a++)gradient[a]+=amount*(gram[a][enter]-gram[a][leave]);
 }
 // The simplex dual gap bounds remaining quadratic objective error; it is
 // separate from signal-fit quality and from reference identifiability.
 const optimalityGap=gap(w);converged=converged||optimalityGap<=1e-10;
 return {...scoreMixture(columns,observed,w),weights:w,converged,optimalityGap};
}
/** Numerical identity on every observed emission coordinate. Fingerprints only
 * accelerate lookup; a full coordinate comparison verifies every equivalence. */
export function emissionFamilies(columns:number[][],tolerance=1e-12):number[][] {
 const buckets=new Map<string,number[][]>(),groups:number[][]=[];
 columns.forEach((column,i)=>{
  const step=Math.max(1,Math.floor(column.length/24)),key=Array.from({length:Math.ceil(column.length/step)},(_,j)=>Math.round(column[j*step]*Math.min(1e6,1/(tolerance*1000)))).join(','),candidates=buckets.get(key)||[];
  const match=candidates.find(group=>column.every((v,k)=>Math.abs(v-columns[group[0]][k])<=tolerance));
  if(match)match.push(i);else{const group=[i];candidates.push(group);groups.push(group);buckets.set(key,candidates);}
 });return groups;
}
/** Working-set simplex solver for large candidate families. Every candidate is
 * tested against the full residual at each step; the global simplex dual gap
 * certifies completion. This changes computation, not the fitted objective. */
export function simplexFit(columns:number[][],observed:number[],positionWeights?:number[],initialWeights?:number[]):SimplexFit {
 if(columns.length<128)return denseSimplexFit(columns,observed,positionWeights,initialWeights);
 const n=observed.length,p=columns.length;if(!n||columns.some(c=>c.length!==n))throw Error('Invalid trace regression dimensions.');
 const weight=(k:number)=>positionWeights?.[Math.floor(k/4)]??1;
 const costs=columns.map(c=>c.reduce((s,v,k)=>s+weight(k)*(v*v-2*v*observed[k]),0));
 let active=initialWeights?.length===p?initialWeights.flatMap((v,i)=>v>1e-10?[i]:[]):[];if(!active.length)active=[costs.indexOf(Math.min(...costs))];
 let previous=initialWeights;
 for(let step=0;step<p;step++){
  const fit=denseSimplexFit(active.map(i=>columns[i]),observed,positionWeights,previous&&active.map(i=>previous![i])),residual=observed.map(v=>-v);
  active.forEach((id,j)=>{if(fit.weights[j])for(let k=0;k<n;k++)residual[k]+=fit.weights[j]*columns[id][k];});
  const gradient=columns.map(c=>c.reduce((s,v,k)=>s+weight(k)*v*residual[k],0)),lambda=active.reduce((s,id,j)=>s+fit.weights[j]*gradient[id],0),optimalityGap=Math.max(0,lambda-Math.min(...gradient))/n;
  const weights=Array(p).fill(0);active.forEach((id,j)=>weights[id]=fit.weights[j]);
  if(optimalityGap<=1e-10)return {...scoreMixture(columns,observed,weights),weights,converged:true,optimalityGap};
  const included=new Set(active),add=gradient.map((value,id)=>({value,id})).filter(x=>!included.has(x.id)&&x.value<lambda-n*1e-10).sort((a,b)=>a.value-b.value).slice(0,4);
  // An active-set convergence stall cannot be resolved by allocating the
  // entire quadratic dictionary matrix. Preserve the actual measured gap.
  if(!add.length)return {...scoreMixture(columns,observed,weights),weights,converged:fit.converged&&optimalityGap<=1e-10,optimalityGap};
  previous=weights;active=active.filter((_,j)=>fit.weights[j]>1e-12).concat(add.map(x=>x.id));
 }
 // Degenerate or unusually dense solutions retain the original full solver.
 return denseSimplexFit(columns,observed,positionWeights,previous);
}
export function scoreMixture(columns:number[][],observed:number[],weights:number[]):{r2:number;rmse:number}{const mean=observed.reduce((a,b)=>a+b,0)/observed.length;let error=0,variance=0;for(let k=0;k<observed.length;k++){let value=0;for(let a=0;a<columns.length;a++)value+=columns[a][k]*weights[a];error+=(observed[k]-value)**2;variance+=(observed[k]-mean)**2;}return {r2:variance>1e-12?1-error/variance:0,rmse:Math.sqrt(error/observed.length)};}
const normalize=(v:number[])=>{const sum=v.reduce((a,b)=>a+b,0);return v.map(x=>sum?x/sum:0);};
export function correctGain(v:number[],gains:number[]){return normalize(v.map((x,b)=>x/gains[b]));}
/** Signed relative phase: positive = previous cycle, negative = next cycle.
 * Both are convex fluorescence mixtures, never negative-intensity subtraction. */
export function phaseNeighbor(p:number,lag:number):number{return p+(lag<0?1:-1);}
export function mixPhase(current:number[],neighbor:number[],lag:number):number[]{const a=Math.abs(lag);return current.map((v,b)=>(1-a)*v+a*(neighbor[b]||0));}
export function phaseVector(sequence:string,p:number,lag:number):number[]{return mixPhase([...TRACE_BASES].map(b=>Number(sequence[p]===b)),[...TRACE_BASES].map(b=>Number(sequence[phaseNeighbor(p,lag)]===b)),lag);}
export interface SignalCalibration {gains:number[];lag:number;used:boolean;trainingBases:number;validationBases:number;rawValidationR2:number;validationR2:number;weights:number[]}
/** Learn only from PRE-CUT WT positions. Every third position is held out.
 * Four bounded gains + one lag parameter; no edited positions or truth labels. */
export function calibrateSignal(observed:number[][],positions:number[],templates:number[][][],precision:'standard'|'fine'='standard'):SignalCalibration{
 const train=positions.filter((_,i)=>i%3!==0),valid=positions.filter((_,i)=>i%3===0),make=(rows:number[],lag:number)=>templates.map(t=>rows.flatMap(p=>mixPhase(t[p],t[phaseNeighbor(p,lag)]||[0,0,0,0],lag)));
 const evaluate=(logs:number[],lag:number)=>{const avg=logs.reduce((a,b)=>a+b,0)/4,gains=logs.map(x=>Math.exp(x-avg)),columns=make(train,lag),y=train.flatMap(p=>correctGain(observed[p],gains)),fit=simplexFit(columns,y);
  const rawPred=(rows:number[])=>rows.flatMap(p=>normalize(Array.from({length:4},(_,b)=>gains[b]*templates.reduce((s,t,h)=>s+fit.weights[h]*((1-Math.abs(lag))*t[p][b]+Math.abs(lag)*(t[phaseNeighbor(p,lag)]?.[b]||0)),0))));
  const score=(rows:number[])=>{const pred=rawPred(rows),truth=rows.flatMap(p=>observed[p]),mean=truth.reduce((a,b)=>a+b,0)/truth.length;let error=0,variance=0;for(let k=0;k<truth.length;k++){error+=(truth[k]-pred[k])**2;variance+=(truth[k]-mean)**2;}return {mse:error/truth.length,r2:variance?1-error/variance:0};};
  return {logs:logs.map(x=>x-avg),lag,gains,weights:fit.weights,training:score(train),validation:score(valid)};
 };
 const identity=evaluate([0,0,0,0],0);if(train.length<40||valid.length<20)return {gains:[1,1,1,1],lag:0,used:false,trainingBases:train.length,validationBases:valid.length,rawValidationR2:identity.validation.r2,validationR2:identity.validation.r2,weights:identity.weights};
 let best=identity;const objective=(r:ReturnType<typeof evaluate>)=>r.training.mse+.0001*r.logs.reduce((s,v)=>s+v*v,0);
 for(const lag of [-.2,-.175,-.15,-.125,-.1,-.075,-.05,-.025,0,.025,.05,.075,.1,.125,.15,.175,.2]){const trial=evaluate(best.logs,lag);if(objective(trial)<objective(best))best=trial;}
 for(const step of (precision==='fine'?[.2,.1,.05,.025,.0125,.00625,.003125,.0015625]:[.2,.1,.05,.025]))for(let round=0;round<2;round++){
  for(let b=0;b<4;b++)for(const direction of [-1,1]){const logs=[...best.logs];logs[b]+=step*direction;const avg=logs.reduce((a,b)=>a+b,0)/4;if(Math.max(...logs.map(x=>x-avg))>Math.log(1.8)||Math.min(...logs.map(x=>x-avg))<Math.log(.55))continue;const trial=evaluate(logs,best.lag);if(objective(trial)<objective(best))best=trial;}
  for(const direction of [-1,1]){const lag=best.lag+direction*(precision==='fine'?Math.min(.0125,step/2):.0125);if(lag<-.2||lag>.2)continue;const trial=evaluate(best.logs,lag);if(objective(trial)<objective(best))best=trial;}
 }
 const used=best.validation.mse<identity.validation.mse*.95&&best.validation.mse+1e-5<identity.validation.mse;const chosen=used?best:identity;
 return {gains:chosen.gains,lag:chosen.lag,used,trainingBases:train.length,validationBases:valid.length,rawValidationR2:identity.validation.r2,validationR2:chosen.validation.r2,weights:chosen.weights};
}

export function gramCondition(columns:number[][]):{rank:number;condition:number|null} {
  const n=columns.length;if(!n)return {rank:0,condition:null};
  const a=Array.from({length:n},()=>new Float64Array(n));
  const norms=columns.map(c=>Math.sqrt(c.reduce((s,v)=>s+v*v,0)));
  for(let i=0;i<n;i++)for(let j=0;j<=i;j++){let sum=0;for(let k=0;k<columns[i].length;k++)sum+=columns[i][k]*columns[j][k];a[i][j]=a[j][i]=norms[i]&&norms[j]?sum/(norms[i]*norms[j]):0;}
  for(let sweep=0;sweep<35;sweep++) {
    let largest=0;
    for(let p=0;p<n;p++)for(let q=p+1;q<n;q++) {
      const off=a[p][q];largest=Math.max(largest,Math.abs(off));if(Math.abs(off)<1e-11)continue;
      const tau=(a[q][q]-a[p][p])/(2*off),t=(tau>=0?1:-1)/(Math.abs(tau)+Math.sqrt(1+tau*tau)),c=1/Math.sqrt(1+t*t),s=t*c;
      const pp=a[p][p],qq=a[q][q];a[p][p]=pp-t*off;a[q][q]=qq+t*off;a[p][q]=a[q][p]=0;
      for(let k=0;k<n;k++)if(k!==p&&k!==q){const kp=a[k][p],kq=a[k][q];a[k][p]=a[p][k]=c*kp-s*kq;a[k][q]=a[q][k]=s*kp+c*kq;}
    }
    if(largest<1e-9)break;
  }
  const values=a.map((row,i)=>Math.max(0,row[i])),max=Math.max(...values),min=Math.min(...values),rank=values.filter(v=>v>max*1e-8).length;
  return {rank,condition:min>max*1e-8?Math.sqrt(max/min):null};
}

/** Mask isolated poor control calls without erasing every neighbouring model row.
 * Requires a quality-passing central call and ≥80% passing shifted donors.
 * Not applicable to a pooled-control consensus: mixed peaks can lower Phred. */
export function shiftedControlQuality(quality:Uint8Array|null,center:number,radius:number,minimum:number,forwardRadius=radius):number{
 if(!quality)return 1;
 if(center-radius<0||center+forwardRadius>=quality.length||quality[center]<minimum)return 0;
 let passed=0;for(let p=center-radius;p<=center+forwardRadius;p++)if(quality[p]>=minimum)passed++;
 const fraction=passed/(radius+forwardRadius+1);return fraction>=.8?fraction:0;
}
/** Score normalized emitted signal in the original dye space, without quality
 * weighting. Shared by single and pooled QC; fractions are never altered here. */
export function rawSignalScore(columns:number[][],observed:number[],weights:number[],gains:number[]){
 const mean=observed.reduce((a,b)=>a+b,0)/observed.length;let error=0,variance=0;const positionRmse:number[]=[];
 for(let j=0;j<observed.length/4;j++){let localError=0;const predicted=Array.from({length:4},(_,b)=>gains[b]*columns.reduce((sum,c,h)=>sum+weights[h]*c[j*4+b],0)),total=predicted.reduce((a,b)=>a+b,0);for(let b=0;b<4;b++){const value=observed[j*4+b],residual=(value-(total?predicted[b]/total:0))**2;error+=residual;localError+=residual;variance+=(value-mean)**2;}positionRmse.push(Math.sqrt(localError/4));}
 return {r2:variance>1e-12?1-error/variance:0,rmse:Math.sqrt(error/observed.length),positionRmse};
}
