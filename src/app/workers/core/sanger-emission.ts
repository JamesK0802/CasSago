import {SangerGroupSupport} from '../../models/sanger-analysis';
import {simplexFit,SimplexFit,mixPhase,phaseNeighbor} from './sanger-signal';
export interface ContextEmission {theta:number[];used:boolean;trainingPositions:number;validationPositions:number;rawValidationMse:number;validationMse:number}
const bases='ACGT';
export function sequenceYield(sequence:string,p:number,theta:number[]):number{const l=bases.indexOf(sequence[p-1]),b=bases.indexOf(sequence[p]),r=bases.indexOf(sequence[p+1]);return Math.exp((l>=0?theta[l]:0)+(r>=0?theta[4+r]:0)+(l>=0&&b>=0&&r>=0?(theta[8+16*l+4*b+r]||0):0));}
export function emissionVector(sequence:string,p:number,lag:number,theta:number[]):number[]{const q=phaseNeighbor(p,lag),current=sequenceYield(sequence,p,theta),neighbor=sequenceYield(sequence,q,theta);return mixPhase([...bases].map(b=>Number(sequence[p]===b)*current),[...bases].map(b=>Number(sequence[q]===b)*neighbor),lag);}
/** Bounded categorical neighbour effects, learned from unedited control only.
 * The exact same model is scored on every-third held-out control position. */
export function fitContextEmission(sequences:string[],positions:number[],observed:number[][],shares:number[],lag:number):ContextEmission{
 const train=positions.filter((_,i)=>i%3!==0),valid=positions.filter((_,i)=>i%3===0),score=(theta:number[],rows:number[])=>{let error=0;for(const p of rows){const predicted=[0,0,0,0];sequences.forEach((s,h)=>emissionVector(s,p,lag,theta).forEach((v,b)=>predicted[b]+=shares[h]*v));const sum=predicted.reduce((a,b)=>a+b,0);for(let b=0;b<4;b++)error+=(observed[p][b]-(sum?predicted[b]/sum:0))**2;}return error/Math.max(1,4*rows.length);};
 const zero=Array(8).fill(0);let theta=[...zero];const objective=(v:number[])=>score(v,train)+.00002*v.reduce((s,x)=>s+x*x,0);
 if(train.length>=80&&valid.length>=30)for(const step of [.08,.04,.02,.01,.005])for(let sweep=0;sweep<2;sweep++)for(let k=0;k<8;k++)for(const sign of [-1,1]){const trial=[...theta];trial[k]+=sign*step;for(const offset of [0,4]){const mean=trial.slice(offset,offset+4).reduce((a,b)=>a+b,0)/4;for(let j=offset;j<offset+4;j++)trial[j]-=mean;}if(trial.some(x=>Math.abs(x)>.25))continue;if(objective(trial)<objective(theta))theta=trial;}
 const raw=score(zero,valid),fitted=score(theta,valid),used=fitted<raw*.95&&fitted+1e-6<raw;
 return {theta:used?theta:zero,used,trainingPositions:train.length,validationPositions:valid.length,rawValidationMse:raw,validationMse:used?fitted:raw};
}
/** Normalized fluorescence from unnormalized template emissions. Iteratively
 * update cycle totals; no truth fractions or fixed allele sparsity threshold. */
export function fitEmissionMixture(columns:number[][],observed:number[],positionWeights?:number[],maxIterations=24,initialWeights?:number[]):SimplexFit{
 const n=observed.length,totals=columns.map(c=>Array.from({length:n/4},(_,p)=>c.slice(4*p,4*p+4).reduce((s,v)=>s+v,0)));
 const evaluate=(w:number[])=>{const predicted=Array(n).fill(0),sum=Array(n/4).fill(0);for(let i=0;i<columns.length;i++)if(w[i])for(let p=0;p<n/4;p++){sum[p]+=w[i]*totals[i][p];for(let b=0;b<4;b++)predicted[4*p+b]+=w[i]*columns[i][4*p+b];}let loss=0,raw=0;for(let k=0;k<n;k++){predicted[k]=sum[k>>2]?predicted[k]/sum[k>>2]:0;const residual=predicted[k]-observed[k];raw+=residual*residual;loss+=(positionWeights?.[k>>2]??1)*residual*residual;}return {predicted,sum,loss,raw};};
 let weights=initialWeights?[...initialWeights]:simplexFit(columns.map((c,i)=>c.map((v,k)=>totals[i][Math.floor(k/4)]?v/totals[i][Math.floor(k/4)]:0)),observed,positionWeights).weights,current=evaluate(weights),converged=false,optimalityGap=Infinity;
 for(let iteration=0;iteration<maxIterations;iteration++){
  // The simplex-constrained Gauss–Newton problem uses the derivative of the
  // NORMALIZED four-channel emission, including each molecule's cycle total.
  const linear=columns.map((c,i)=>c.map((v,k)=>current.predicted[k]+(current.sum[k>>2]?(v-current.predicted[k]*totals[i][k>>2])/current.sum[k>>2]:0)));
  const gradient=linear.map(c=>c.reduce((sum,v,k)=>sum+(positionWeights?.[k>>2]??1)*(v-current.predicted[k])*(current.predicted[k]-observed[k]),0));
  optimalityGap=(weights.reduce((s,w,i)=>s+w*gradient[i],0)-Math.min(...gradient))/Math.max(1,n);
  if(optimalityGap<1e-10){converged=true;break;}
  const proposal=simplexFit(linear,observed,positionWeights,weights);let step=1,accepted=false;
  while(step>=1/4096){const next=weights.map((w,i)=>w+step*(proposal.weights[i]-w)),value=evaluate(next);if(value.loss<current.loss-1e-14){weights=next;current=value;accepted=true;break;}step*=.5;}
  if(!accepted){converged=optimalityGap<1e-8;break;}
 }
 const mean=observed.reduce((s,v)=>s+v,0)/n,variance=observed.reduce((s,v)=>s+(v-mean)**2,0);
 return {weights,r2:variance?1-current.raw/variance:0,rmse:Math.sqrt(current.raw/n),converged,optimalityGap};
}

export interface EmissionSupportOptions {refit?:boolean;templateResidual?:number[]}
export interface EmissionSupportEvidence {group:number[];threshold:number;mass:number;standardError:number;score:number}
/** Support on the simplex tangent space, with cycle-block residual covariance.
 * A newly proposed length pays for its searched junction origins. Reference
 * allocation remains an estimate: ambiguity between homeologs is not evidence
 * that an edit is absent. Refit retained lengths without a fraction penalty.
 * This engineering screen is not a confidence interval. */
export function screenEmissionSupport(columns:number[][],observed:number[],fit:SimplexFit,mandatory:number[],positionWeights?:number[],groups?:number[][],originHypotheses?:number[],options:EmissionSupportOptions={}):{fit:SimplexFit;removed:number[];threshold:number;scores?:EmissionSupportEvidence[]}{
 const {templateResidual}=options;
 const n=observed.length,active=fit.weights.flatMap((w,i)=>w>1e-9||mandatory.includes(i)?[i]:[]),pivot=active.reduce((a,b)=>fit.weights[a]>fit.weights[b]?a:b),indices=active.filter(i=>i!==pivot),p=indices.length,threshold=Math.sqrt(Math.log(Math.max(3,n)));
 if(p<1)return {fit,removed:[],threshold};
 const jac=indices.map(()=>new Float64Array(n)),residual=new Float64Array(n),rowWeight=(k:number)=>positionWeights?.[Math.floor(k/4)]??1;let sse=0;
 for(let j=0;j<n/4;j++){const vector=Array.from({length:4},(_,b)=>columns.reduce((s,c,h)=>s+c[j*4+b]*fit.weights[h],0)),total=vector.reduce((a,b)=>a+b,0),totals=columns.map(c=>c.slice(j*4,j*4+4).reduce((a,b)=>a+b,0));if(!total)continue;for(let b=0;b<4;b++){const pred=vector[b]/total;residual[j*4+b]=observed[j*4+b]-pred;sse+=rowWeight(j*4+b)*residual[j*4+b]**2;indices.forEach((i,k)=>{jac[k][j*4+b]=(columns[i][j*4+b]-pred*totals[i]-columns[pivot][j*4+b]+pred*totals[pivot])/total;});}}
 const matrix=Array.from({length:p},(_,i)=>Float64Array.from({length:2*p},(_,j)=>j<p?jac[i].reduce((s,v,k)=>s+rowWeight(k)*v*jac[j][k],0):Number(j-p===i)));
 for(let k=0;k<p;k++){let pivotRow=k;for(let i=k+1;i<p;i++)if(Math.abs(matrix[i][k])>Math.abs(matrix[pivotRow][k]))pivotRow=i;[matrix[k],matrix[pivotRow]]=[matrix[pivotRow],matrix[k]];const value=matrix[k][k];if(Math.abs(value)<1e-8)return {fit,removed:[],threshold};for(let j=0;j<2*p;j++)matrix[k][j]/=value;for(let i=0;i<p;i++)if(i!==k){const scale=matrix[i][k];for(let j=0;j<2*p;j++)matrix[i][j]-=scale*matrix[k][j];}}
 const variance=sse/Math.max(1,n-active.length);
 // Include the dependent pivot when testing group mass. Cycle-block scores
 // retain correlations between dyes and nearby bases, unlike IID residuals.
 const scores:EmissionSupportEvidence[]=[];const removed=groups?groups.filter(group=>!group.some(i=>mandatory.includes(i))).flatMap(group=>{
  const coefficients=indices.map(i=>Number(group.includes(i))-Number(group.includes(pivot)));
  const v=coefficients.reduce((sum,a,k)=>sum+a*coefficients.reduce((s,b,l)=>s+b*matrix[k][p+l],0),0);
  const mass=group.reduce((sum,i)=>sum+fit.weights[i],0);
  const direction=indices.map((_,k)=>coefficients.reduce((s,a,l)=>s+matrix[k][p+l]*a,0));
  const blockLength=Math.max(2,Math.ceil(Math.cbrt(n/4))),blocks=Array.from({length:Math.ceil(n/4/blockLength)},()=>0);
  for(let k=0;k<n;k++)blocks[Math.floor(k/4/blockLength)]+=rowWeight(k)*residual[k]*direction.reduce((sum,a,l)=>sum+a*jac[l][k],0);
  const mean=blocks.reduce((s,x)=>s+x,0)/blocks.length;
  const robustVariance=blocks.reduce((s,x)=>s+(x-mean)**2,0)*blocks.length/Math.max(1,blocks.length-1);
  // WT-control model error is independent of the edited-sample fit. A false
  // class can absorb systematic template error and leave very little sample
  // residual, so residual-only precision can be spuriously high.
  let templateVariance=0;
  if(templateResidual){
   const controlBlocks=blocks.map(()=>0);
   for(let k=0;k<n;k++)controlBlocks[Math.floor(k/4/blockLength)]+=rowWeight(k)*templateResidual[k]*direction.reduce((sum,a,l)=>sum+a*jac[l][k],0);
   const bias=controlBlocks.reduce((s,x)=>s+x,0),average=bias/controlBlocks.length;
   templateVariance=Math.max(bias*bias,controlBlocks.reduce((s,x)=>s+(x-average)**2,0)*controlBlocks.length/Math.max(1,controlBlocks.length-1));
  }
  const groupThreshold=Math.sqrt(threshold**2+2*Math.log(Math.max(1,...group.map(i=>originHypotheses?.[i]??1))));scores.push({group,threshold:groupThreshold,mass,standardError:Math.sqrt(Math.max(0,Math.max(variance*v,robustVariance,templateVariance))),score:mass/Math.max(1e-12,Math.sqrt(Math.max(0,Math.max(variance*v,robustVariance,templateVariance))))});return mass<groupThreshold*Math.sqrt(Math.max(0,Math.max(variance*v,robustVariance,templateVariance)))?group.filter(i=>active.includes(i)):[];
 }):indices.filter((i,k)=>!mandatory.includes(i)&&fit.weights[i]<threshold*Math.sqrt(Math.max(0,variance*matrix[k][p+k])));
 const uniqueRemoved=[...new Set(removed)];removed.splice(0,removed.length,...uniqueRemoved);const keep=active.filter(i=>!removed.includes(i));
 if(options.refit===false||!removed.length)return {fit,removed,threshold,scores};
 const refit=fitEmissionMixture(keep.map(i=>columns[i]),observed,positionWeights),weights=Array(columns.length).fill(0);keep.forEach((i,k)=>weights[i]=refit.weights[k]);return {fit:{...refit,weights},removed,threshold,scores};
}

/** Polish the actual normalized-fluorescence objective when fixed-point scaling
 * stalls, notably when short edited products terminate within the fit region.
 * Uses projected gradients with backtracking; no truth or extra edit models. */
export function fitNormalizedEmissionMixture(columns:number[][],observed:number[]):SimplexFit{
 const initial=fitEmissionMixture(columns,observed,undefined,64);if(initial.converged)return initial;
 const project=(values:number[])=>{const sorted=[...values].sort((a,b)=>b-a);let sum=0,threshold=0;for(let k=0;k<sorted.length;k++){sum+=sorted[k];threshold=(sum-1)/(k+1);if(k===sorted.length-1||sorted[k+1]<=threshold)break;}return values.map(v=>Math.max(0,v-threshold));};
 const totals=columns.map(c=>Array.from({length:observed.length/4},(_,p)=>c.slice(4*p,4*p+4).reduce((s,v)=>s+v,0)));
 const evaluate=(weights:number[],gradient=false)=>{let error=0;const g=weights.map(()=>0);for(let p=0;p<observed.length/4;p++){const v=[0,0,0,0];for(let h=0;h<columns.length;h++)for(let b=0;b<4;b++)v[b]+=weights[h]*columns[h][4*p+b];const total=v.reduce((s,x)=>s+x,0);for(let b=0;b<4;b++){const predicted=total?v[b]/total:0,residual=predicted-observed[4*p+b];error+=residual*residual;if(gradient&&total>1e-9)for(let h=0;h<columns.length;h++)g[h]+=2*residual*(columns[h][4*p+b]-predicted*totals[h][p])/total;}}return {loss:error/observed.length,gradient:g.map(v=>v/observed.length)};};
 let weights=[...initial.weights],value=evaluate(weights).loss,step=1,converged=false;
 for(let iteration=0;iteration<600;iteration++){
  const current=evaluate(weights,true);let accepted=false,next=weights,nextValue=value;
  for(let trial=0;trial<24;trial++){next=project(weights.map((w,k)=>w-step*current.gradient[k]));nextValue=evaluate(next).loss;const descent=next.reduce((sum,w,k)=>sum+current.gradient[k]*(w-weights[k]),0);if(nextValue<=value+1e-4*descent+1e-15){accepted=true;break;}step*=.5;}
  if(!accepted)break;const delta=Math.max(...next.map((w,k)=>Math.abs(w-weights[k])));weights=next;value=nextValue;
  const mapping=project(weights.map((w,k)=>w-current.gradient[k]));if(Math.max(...mapping.map((w,k)=>Math.abs(w-weights[k])))<1e-7&&delta<1e-6){converged=true;break;}step=Math.min(100,step*1.4);
 }
 const mean=observed.reduce((s,v)=>s+v,0)/observed.length,variance=observed.reduce((s,v)=>s+(v-mean)**2,0)/observed.length;
 return {weights,rmse:Math.sqrt(value),r2:variance?1-value/variance:0,converged};
}

/** Residual-based engineering evidence; not a probability of a real edit.
 * Missing/singular evidence stays unresolved and never inherits a QC pass. */
export function emissionGroupSupport(scores:EmissionSupportEvidence[]|undefined,index:number,prediction?:{foldGains:number[];passed:boolean}):SangerGroupSupport {
 const evidence=scores?.find(s=>s.group.includes(index));
 return {status:evidence&&evidence.mass>0&&evidence.score>=evidence.threshold&&(!prediction||prediction.passed)?'supported':'unresolved',score:evidence?.score??null,requiredScore:evidence?.threshold??null,scope:'reference-length',...(prediction?{prediction}:{})};
}
