import type{SangerTrace}from'../../models/sanger-analysis';
import type{SangerLinearOptimizer}from'./sanger-discrete-phase';
export interface RawIntervalFit{weights:number[];violation:number;amplitude:number;sigma:number;activePositions?:number[];residuals?:{position:number;error:number}[];loss?:number;lowerBoundOnly?:boolean}
export type RawSequenceScorer=((sequences:string[],weights:number[])=>RawIntervalFit|null)&{maximumCentroidMse?:number;boundsFor?:(positions:number[],response:{amplitude:number;sigma:number})=>{lower:number[][];upper:number[][]}};
interface Row{ids:number[];d:number;lo:number;hi:number}
// Math.round bins are [k-.5,k+.5), not closed on both sides. A
// sub-count numerical interior prevents distinct ADC bins sharing one endpoint.
const ADC_UPPER_INTERIOR=1e-5;
/** Only enabled for the separately validated, integer-ADC Gaussian response.
 * Fit fluorescence amplitudes against rounding intervals in BOTH original
 * traces. Width is learned from genotype-independent peak ratios. A sequential
 * linearization optimizes the narrow width interval; every returned objective
 * is recomputed using the nonlinear response, not the LP surrogate. */
export function rawIntervalScorer(control:SangerTrace,sample:SangerTrace,refs:string[],positions:number[],ct:number,st:number,bg:number[],ebg:number[],calibration:{sigma:number;amplitude:number;coefficients:number[]},optimizer:SangerLinearOptimizer,observed?:number[]):RawSequenceScorer|null{
 let lower=1/16,upper=2;
 for(const [trace,background]of [[control,bg],[sample,ebg]]as const)for(const p of trace.peaks)for(let b=0;b<4;b++){
  const center=trace.channels['ACGT'[b]as'A'][p]-background[b];if(center<=.5)continue;
  for(let d=1;d<=4;d++){const side=trace.channels['ACGT'[b]as'A'][p+d]-background[b];if(side<0||!Number.isInteger(trace.channels['ACGT'[b]as'A'][p+d]))return null;lower=Math.max(lower,Math.log((center-.5)/(side+.5))/(d*d));if(side>.5)upper=Math.min(upper,Math.log((center+.5)/(side-.5))/(d*d));}
 }
 if(lower>upper||!Number.isFinite(lower)||!Number.isFinite(upper))return null;
 const rows=(trace:SangerTrace,dict:string[],trim:number,background:number[]):Row[]=>{
  const map=new Map<string,Row>();
  for(const p of positions)for(let b=0;b<4;b++){
   const ids=dict.flatMap((s,i)=>s[p]==='ACGT'[b]?[i]:[]),pattern=ids.join(',');
   for(let d=-4;d<=4;d++){const value=trace.channels['ACGT'[b]as'A'][trace.peaks[p-trim]+d]-background[b],key=pattern+'/'+Math.abs(d),old=map.get(key);if(old){old.lo=Math.max(old.lo,value-.5);old.hi=Math.min(old.hi,value+.5-ADC_UPPER_INTERIOR);}else map.set(key,{ids,d:Math.abs(d),lo:value-.5,hi:value+.5-ADC_UPPER_INTERIOR});}
  }
  return[...map.values()];
 };
 // A raw-compatible unit emission must lie within these conservative bounds,
 // allowing every amplitude/width consistent with the WT and measured shapes.
 // Its distance to the extracted centroid therefore cannot exceed this bound.
 let ampLo=0,ampHi=Infinity;
 for(const p of positions){const same=refs.every(s=>s[p]===refs[0][p]);let lo=0,hi=0;if(same){const value=control.channels[refs[0][p]as'A'][control.peaks[p-ct]]-bg['ACGT'.indexOf(refs[0][p])];lo=value-.5;hi=value+.5-ADC_UPPER_INTERIOR;}else for(let b=0;b<4;b++){const value=control.channels['ACGT'[b]as'A'][control.peaks[p-ct]]-bg[b];lo+=Math.max(0,value-.5);hi+=Math.max(0,value+.5-ADC_UPPER_INTERIOR);}ampLo=Math.max(ampLo,lo);ampHi=Math.min(ampHi,hi);}
 if(!(ampLo>0)||ampLo>ampHi)return null;
 let centroidBound=0;const featureBounds:{value:number;low:number;high:number}[][]=[];
 if(observed){for(let k=0;k<positions.length;k++)for(let b=0;b<4;b++){let lo=0,hi=1;const feature:{value:number;low:number;high:number}[]=[];const p=sample.peaks[positions[k]-st];for(let d=-4;d<=4;d++){const value=sample.channels['ACGT'[b]as'A'][p+d]-ebg[b],responseLo=ampLo*Math.exp(-d*d*upper),responseHi=ampHi*Math.exp(-d*d*lower);feature.push({value,low:responseLo,high:responseHi});lo=Math.max(lo,(value-.5)/responseHi);hi=Math.min(hi,(value+.5-ADC_UPPER_INTERIOR)/responseLo);}if(lo>hi)return null;centroidBound+=Math.max((observed[4*k+b]-lo)**2,(observed[4*k+b]-hi)**2);featureBounds.push(feature);}centroidBound/=observed.length;}
 const controlRows=rows(control,refs,ct,bg);let incumbent=Infinity;
 const score:RawSequenceScorer=(sequences,weights)=>{
  const n=sequences.length,all=[...rows(sample,sequences,st,ebg),...controlRows.map(r=>({...r,ids:r.ids.map(i=>i+n)}))];
  // If an identical emission pattern is observed at incompatible ADC heights,
  // no choice of its weights or Gaussian width can reconcile those measurements.
  // This is a lower bound, not a fabricated fitted solution.
  const bound=all.reduce((v,row)=>Math.max(v,(row.lo-row.hi)/2),0);
  if(bound>incumbent+1e-8)return{weights,violation:bound,amplitude:calibration.amplitude,sigma:calibration.sigma,lowerBoundOnly:true};
  const scale=calibration.amplitude,radius=(upper-lower)/2||1;let inverse=Math.max(lower,Math.min(upper,1/calibration.sigma)),coefficients=[...weights,...calibration.coefficients.map(v=>v/scale)],best:RawIntervalFit|null=null;
  const persistent=optimizer.createModel?.();let basis:unknown;
  const options={output_flag:false,solver:'simplex',primal_feasibility_tolerance:1e-10,dual_feasibility_tolerance:1e-10,time_limit:5};persistent?.options.set(options);
  try{
  for(let step=0;step<5;step++){
   const lines:string[]=[];
   if(!persistent){
    lines.push('Minimize',' obj: error','Subject To',' flux: '+sequences.map((_,i)=>'x'+i).join(' + ')+' - '+refs.map((_,i)=>'x'+(n+i)).join(' - ')+' = 0');let k=0;
   for(const row of all){const shape=Math.exp(-row.d*row.d*inverse),mass=row.ids.reduce((s,i)=>s+coefficients[i],0),derivative=-row.d*row.d*shape*mass*radius,terms=row.ids.map(i=>shape+' x'+i).join(' + '),expression=(terms?terms+' ':'')+(derivative<0?'- ':'+ ')+Math.abs(derivative)+' delta';const error='error';lines.push(' r'+k+++': '+expression+' + '+error+' >= '+(row.lo/scale),' r'+k+++': '+expression+' - '+error+' <= '+(row.hi/scale));}
   lines.push('Bounds',` ${(lower-inverse)/radius} <= delta <= ${(upper-inverse)/radius}`,'End');
   }
   let values:ArrayLike<number>,delta:number;
   if(persistent){
    const count=n+refs.length,dc=count,ec=count+1,numCols=ec+1,starts=[0],indices:number[]=[],matrixValues:number[]=[],rowLower=[0],rowUpper=[0];
    const addRow=(ids:number[],coefficients:number[],lo:number,hi:number)=>{indices.push(...ids);matrixValues.push(...coefficients);starts.push(indices.length);rowLower.push(lo);rowUpper.push(hi);};
    indices.push(...Array.from({length:count},(_,i)=>i));matrixValues.push(...Array.from({length:count},(_,i)=>i<n?1:-1));starts.push(indices.length);
    all.forEach((row,r)=>{const shape=Math.exp(-row.d*row.d*inverse),mass=row.ids.reduce((s,i)=>s+coefficients[i],0),derivative=-row.d*row.d*shape*mass*radius,ids=[...row.ids,dc,ec],v=[...row.ids.map(()=>shape),derivative,1];addRow(ids,v,row.lo/scale,1e30);addRow(ids,[...v.slice(0,-1),-1],-1e30,row.hi/scale);});
    const colLower=Array(numCols).fill(0),colUpper=Array(numCols).fill(1e30);colLower[dc]=(lower-inverse)/radius;colUpper[dc]=(upper-inverse)/radius;
    persistent.passModel({numCols,numRows:rowLower.length,colCost:Array.from({length:numCols},(_,i)=>i>=ec?1:0),colLower,colUpper,rowLower,rowUpper,matrix:{format:'csr',numCols,numRows:rowLower.length,starts,indices,values:matrixValues}});
    if(basis)persistent.setBasis(basis);const solved=persistent.run();if(solved.modelStatus!==7)return best;values=persistent.getSolution().colValue;delta=values[dc];basis=persistent.getBasis();
   }else{const result=optimizer.solve(lines.join('\n'),options);if(result.Status!=='Optimal')return best;values=Array.from({length:n+refs.length},(_,i)=>result.Columns['x'+i]?.Primal||0);delta=result.Columns['delta']?.Primal||0;}
   const next=Array.from({length:n+refs.length},(_,i)=>Math.max(0,values[i])),nextInverse=inverse+delta*radius,amplitude=next.slice(0,n).reduce((s,v)=>s+v,0)*scale;if(!(amplitude>0))return best;
   let violation=0,loss=0;for(const row of all){const prediction=scale*Math.exp(-row.d*row.d*nextInverse)*row.ids.reduce((s,i)=>s+next[i],0),error=Math.max(0,row.lo-prediction,prediction-row.hi);violation=Math.max(violation,error);loss+=error;}
   const fit={weights:next.slice(0,n).map(v=>v*scale/amplitude),violation,loss:violation,amplitude,sigma:1/nextInverse};if(!best||violation<best.loss!)best=fit;
   const widthChange=Math.abs(nextInverse-inverse);inverse=nextInverse;coefficients=next;if(violation<1e-8||widthChange<1e-12)break;
  }
  }finally{persistent?.dispose();}
  if(best){incumbent=Math.min(incumbent,best.loss!);best.activePositions=[];best.residuals=[];
   for(const p of positions){const mass=[0,0,0,0];sequences.forEach((seq,i)=>mass['ACGT'.indexOf(seq[p])]+=best!.weights[i]);let error=0;for(let b=0;b<4;b++)for(let d=-4;d<=4;d++){const measured=sample.channels['ACGT'[b]as'A'][sample.peaks[p-st]+d]-ebg[b],prediction=best.amplitude*mass[b]*Math.exp(-d*d/best.sigma);error=Math.max(error,Math.max(measured-.5-prediction,prediction-measured-.5+ADC_UPPER_INTERIOR));}if(error>Math.max(1e-7,scale*1e-9)){best.activePositions.push(p);best.residuals.push({position:p,error});}}
  }return best;
 };
 score.boundsFor=(at,response)=>{const low:number[][]=[],high:number[][]=[];for(const p of at){const l=[0,0,0,0],u=[1,1,1,1];for(let b=0;b<4;b++)for(let d=-4;d<=4;d++){const shape=response.amplitude*Math.exp(-d*d/response.sigma),value=sample.channels['ACGT'[b]as'A'][sample.peaks[p-st]+d]-ebg[b];l[b]=Math.max(l[b],(value-.5)/shape);u[b]=Math.min(u[b],(value+.5-ADC_UPPER_INTERIOR)/shape);}low.push(l);high.push(u);}return{lower:low,upper:high};};
 // A zero-error box only rejects exact compatibility, not improvements to a
 // still-imperfect incumbent. Expand every ADC interval by the incumbent's
 // maximum error, and permit the corresponding control-amplitude uncertainty.
 // This bound can reject only dictionaries unable to improve that SAME loss.
 let lastError=Infinity,lastBound=Infinity;
 if(observed)Object.defineProperty(score,'maximumCentroidMse',{get:()=>{
  if(!Number.isFinite(incumbent))return Infinity;if(lastError===incumbent)return lastBound;
  lastError=incumbent;const lowAmplitude=ampLo-4*incumbent,highAmplitude=ampHi+4*incumbent;if(lowAmplitude<=0)return lastBound=Infinity;
  let total=0;for(let i=0;i<featureBounds.length;i++){let lo=0,hi=1;for(const f of featureBounds[i]){lo=Math.max(lo,(f.value-.5-incumbent)/(f.high*highAmplitude/ampHi));hi=Math.min(hi,(f.value+.5-ADC_UPPER_INTERIOR+incumbent)/(f.low*lowAmplitude/ampLo));}total+=Math.max((observed[i]-lo)**2,(observed[i]-hi)**2);}return lastBound=total/observed.length;
 }});return score;
}
