import{activeSimplex}from'./simplex-active-set';
import {SimplexFit} from './sanger-signal';
import {fitNormalizedEmissionMixture} from './sanger-emission';
/** Immutable waveforms/observations allow reusable fits and Gram products.
 * When all candidates have the SAME cycle total, normalization is independent
 * of the mixture weights and the full problem is exactly simplex least squares.
 * Unequal dye yields or product termination retain the general nonlinear fit. */
export class SangerFitCache {
 private fits=new Map<string,SimplexFit>();private gram=new Map<number,number>();private compatible:boolean[]=[];private inverse:number[];private variance:number;
 hits=0;misses=0;linearFits=0;nonlinearFits=0;pairHits=0;pairMisses=0;
 constructor(private columns:number[][],private observed:number[],private limit=512){this.inverse=Array.from({length:observed.length/4},(_,p)=>{const c=columns[0],k=4*p,total=c[k]+c[k+1]+c[k+2]+c[k+3];return total?1/total:0;});const mean=observed.reduce((s,v)=>s+v,0)/observed.length;this.variance=observed.reduce((s,v)=>s+(v-mean)**2,0);}
 private sharedTotal(id:number){if(this.compatible[id]!==undefined)return this.compatible[id];const c=this.columns[id],base=this.columns[0];for(let k=0;k<c.length;k+=4)if(Math.abs(c[k]+c[k+1]+c[k+2]+c[k+3]-(base[k]+base[k+1]+base[k+2]+base[k+3]))>1e-12)return this.compatible[id]=false;return this.compatible[id]=true;}
 private product(a:number,b:number){const key=Math.max(a,b)*1048576+Math.min(a,b),old=this.gram.get(key);if(old!==undefined){this.pairHits++;this.gram.delete(key);this.gram.set(key,old);return old;}this.pairMisses++;const x=this.columns[a],y=this.columns[b];let value=0;for(let k=0;k<this.observed.length;k++)value+=(x[k]*this.inverse[k>>2]-this.observed[k])*(y[k]*this.inverse[k>>2]-this.observed[k]);this.gram.set(key,value);if(this.gram.size>200000)this.gram.delete(this.gram.keys().next().value!);return value;}
 private linear(ids:number[],initial?:number[]):SimplexFit{
  const p=ids.length,n=this.observed.length,G=ids.map((a,i)=>Float64Array.from({length:p},(_,j)=>j<=i?this.product(a,ids[j]):0));for(let i=0;i<p;i++)for(let j=0;j<i;j++)G[j][i]=G[i][j];
  const solved=activeSimplex(G,n,1e-13,initial),weights=solved.weights,converged=solved.converged,gap=solved.optimalityGap;
  let error=0;for(let k=0;k<n;k++){let predicted=0;for(let j=0;j<p;j++)predicted+=weights[j]*this.columns[ids[j]][k];error+=(this.observed[k]-predicted*this.inverse[k>>2])**2;}
  return {weights,rmse:Math.sqrt(error/n),r2:this.variance?1-error/this.variance:0,converged,optimalityGap:gap};
 }
 fit(ids:number[],initial?:number[]):SimplexFit{const key=ids.join('/'),hit=this.fits.get(key);if(hit){this.hits++;this.fits.delete(key);this.fits.set(key,hit);return hit;}this.misses++;let result:SimplexFit;if(ids.every(id=>id<1048576&&this.sharedTotal(id))){this.linearFits++;result=this.linear(ids,initial);if(!result.converged){this.nonlinearFits++;result=fitNormalizedEmissionMixture(ids.map(id=>this.columns[id]),this.observed);}}else{this.nonlinearFits++;result=fitNormalizedEmissionMixture(ids.map(id=>this.columns[id]),this.observed);}this.fits.set(key,result);if(this.fits.size>this.limit)this.fits.delete(this.fits.keys().next().value!);return result;}
}
/** Same raw-dye loss as rawSignalScore, without allocating QC statistics for
 * every trial. Detailed per-cycle QC is still computed for the final fit. */
export function rawMixtureMse(columns:number[][],observed:number[],weights:number[],gains:number[]){let error=0;for(let p=0;p<observed.length;p+=4){let a=0,c=0,g=0,t=0;for(let j=0;j<columns.length;j++){const col=columns[j],w=weights[j];a+=w*col[p];c+=w*col[p+1];g+=w*col[p+2];t+=w*col[p+3];}a*=gains[0];c*=gains[1];g*=gains[2];t*=gains[3];const total=a+c+g+t;if(total){a/=total;c/=total;g/=total;t/=total;}error+=(observed[p]-a)**2+(observed[p+1]-c)**2+(observed[p+2]-g)**2+(observed[p+3]-t)**2;}return error/observed.length;}
