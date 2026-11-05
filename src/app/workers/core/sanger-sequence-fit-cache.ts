import type{RawSequenceScorer}from './sanger-raw-interval-fit';
import{SangerFitCache}from'./sanger-fit-cache';
import{activeSimplex}from'./simplex-active-set';import type{SimplexFit}from'./sanger-signal';
/** Exact specialization of simplex least squares for unit A/C/G/T emissions.
 * A column is a base per cycle. Its Gram product is the number of agreements
 * minus its two observed-base sums, plus the common observed squared norm.
 * No waveform approximation or biological independence assumption is added. */
export interface SequenceFit extends SimplexFit{selectionLoss?:number;selectionTolerance?:number;rawResponse?:{amplitude:number;sigma:number};rawActivePositions?:number[];rawResiduals?:{position:number;error:number}[]}
export const sequenceFitTolerance=(fit:SequenceFit)=>fit.selectionTolerance??Math.max(1e-15,2*(fit.optimalityGap||0));
export const sequenceFitLoss=(fit:SequenceFit)=>fit.selectionLoss??fit.rmse**2;
export class SangerSequenceFitCache{
 private encoded=new Map<number,Uint8Array>();private observedSum=new Map<number,number>();private gram=new Map<string,number>();private fits=new Map<string,SequenceFit>();private yy:number;private variance:number;
 hits=0;misses=0;linearFits=0;nonlinearFits=0;pairHits=0;pairMisses=0;
 constructor(private sequences:string[],private positions:number[],private observed:number[],private limit=512,private rawScorer?:RawSequenceScorer){if(observed.length!==positions.length*4)throw Error('Invalid unit-sequence regression dimensions.');this.yy=observed.reduce((s,v)=>s+v*v,0);const mean=observed.reduce((s,v)=>s+v,0)/observed.length;this.variance=observed.reduce((s,v)=>s+(v-mean)**2,0);}
 private bases(id:number){let code=this.encoded.get(id);if(code)return code;code=Uint8Array.from(this.positions,p=>{const b='ACGT'.indexOf(this.sequences[id][p]);if(b<0)throw Error('Sequence regression requires complete unit A/C/G/T emissions.');return b;});this.encoded.set(id,code);this.observedSum.set(id,code.reduce((s,b,p)=>s+this.observed[4*p+b],0));return code;}
 private product(a:number,b:number){const key=Math.min(a,b)+'/'+Math.max(a,b),old=this.gram.get(key);if(old!==undefined){this.pairHits++;return old;}this.pairMisses++;const x=this.bases(a),y=this.bases(b);let same=0;for(let p=0;p<x.length;p++)same+=Number(x[p]===y[p]);const v=same-this.observedSum.get(a)!-this.observedSum.get(b)!+this.yy;this.gram.set(key,v);if(this.gram.size>200000)this.gram.delete(this.gram.keys().next().value!);return v;}
 column(id:number){return Array.from(this.bases(id),b=>[0,1,2,3].map(k=>Number(b===k))).flat();}
 predict(ids:number[],weights:number[]){const y=Array(this.observed.length).fill(0);ids.forEach((id,i)=>{const code=this.bases(id);for(let p=0;p<code.length;p++)y[4*p+code[p]]+=weights[i];});return y;}
 fit(ids:number[],initial?:number[]):SequenceFit{const key=ids.join('/'),hit=this.fits.get(key);if(hit){this.hits++;return hit;}this.misses++;this.linearFits++;const G=ids.map((a,i)=>Float64Array.from({length:ids.length},(_,j)=>j<=i?this.product(a,ids[j]):0));for(let i=0;i<ids.length;i++)for(let j=0;j<i;j++)G[j][i]=G[i][j];const solved=activeSimplex(G,this.observed.length,1e-13,initial),codes=ids.map(id=>this.bases(id));let error=0;
 for(let p=0;p<this.positions.length;p++){const v=[0,0,0,0];for(let j=0;j<ids.length;j++)v[codes[j][p]]+=solved.weights[j];for(let b=0;b<4;b++)error+=(v[b]-this.observed[4*p+b])**2;}
 const result:SequenceFit=solved.converged?{weights:solved.weights,rmse:Math.sqrt(error/this.observed.length),r2:this.variance?1-error/this.variance:0,converged:solved.converged,optimalityGap:solved.optimalityGap}:new SangerFitCache(ids.map(id=>this.column(id)),this.observed).fit(ids.map((_,i)=>i));if(this.rawScorer){
 const raw=this.rawScorer.maximumCentroidMse!==undefined&&result.rmse**2>this.rawScorer.maximumCentroidMse+1e-15?null:this.rawScorer(ids.map(id=>this.sequences[id]),result.weights);
 if(raw){let measuredError=0;const prediction=this.predict(ids,raw.weights);for(let i=0;i<prediction.length;i++)measuredError+=(prediction[i]-this.observed[i])**2;result.weights=raw.weights;result.rmse=Math.sqrt(measuredError/this.observed.length);result.r2=this.variance?1-measuredError/this.variance:0;result.selectionLoss=raw.loss??raw.violation;result.selectionTolerance=1e-7;if(!raw.lowerBoundOnly)result.rawResponse={amplitude:raw.amplitude,sigma:raw.sigma};result.rawActivePositions=raw.activePositions;result.rawResiduals=raw.residuals;result.converged=!raw.lowerBoundOnly;result.optimalityGap=undefined;}else{result.converged=false;result.selectionLoss=Infinity;}
 }this.fits.set(key,result);if(this.fits.size>this.limit)this.fits.delete(this.fits.keys().next().value!);return result;}
}
