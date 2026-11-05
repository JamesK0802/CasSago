import {SimplexFit,simplexFit,scoreMixture} from './sanger-signal';
const median=(values:number[])=>{const s=[...values].sort((a,b)=>a-b);return s.length?s[Math.floor(s.length/2)]:0;};
export function upstreamNoiseScale(residuals:number[]):number {const center=median(residuals);return Math.max(1e-5,1.4826*median(residuals.map(v=>Math.abs(v-center))));}
/** Positive sparse regression on contrasts to the WT template. The column penalty
 * is standardized by measured contrast norm, so a SNP is not penalized merely
 * because it changes fewer cycles than a long indel. Noise comes from unchanged
 * upstream observations, never the expected edit or a comparator's answer. */
export function fitSparseSingle(columns:number[][],observed:number[],wtIndex:number,sigma:number,positionWeights?:number[]):SimplexFit {
 const p=columns.length,n=observed.length;if(!p||wtIndex<0||wtIndex>=p||columns.some(c=>c.length!==n)||!Number.isFinite(sigma)||sigma<0)throw Error('Invalid single-trace sparse fit.');
 if(p===1)return {weights:[1],...scoreMixture(columns,observed,[1]),converged:true};
 const ids=columns.flatMap((_,i)=>i===wtIndex?[]:[i]),base=columns[wtIndex],d=ids.map(i=>columns[i].map((x,k)=>x-base[k])),y=observed.map((x,k)=>x-base[k]),q=d.length,gram=Array.from({length:q},()=>new Float64Array(q)),rhs=new Float64Array(q),weight=(k:number)=>positionWeights?.[Math.floor(k/4)]??1;
 for(let a=0;a<q;a++){for(let k=0;k<n;k++)rhs[a]+=weight(k)*d[a][k]*y[k];for(let b=0;b<=a;b++){let value=0;for(let k=0;k<n;k++)value+=weight(k)*d[a][k]*d[b][k];gram[a][b]=gram[b][a]=value;}}
 const penalty=Array.from({length:q},(_,i)=>sigma*Math.sqrt(2*Math.log(Math.max(3,p)))*Math.sqrt(Math.max(0,gram[i][i]))),L=Math.max(...gram.map(r=>r.reduce((s,x)=>s+Math.abs(x),0)));
 if(!L)return {weights:columns.map((_,i)=>Number(i===wtIndex)),...scoreMixture(columns,observed,columns.map((_,i)=>Number(i===wtIndex))),converged:true};
 const project=(v:number[])=>{const nonnegative=v.map(x=>Math.max(0,x));if(nonnegative.reduce((s,x)=>s+x,0)<=1)return nonnegative;const sorted=[...v].sort((a,b)=>b-a);let sum=0,threshold=0;for(let k=0;k<sorted.length;k++){sum+=sorted[k];threshold=(sum-1)/(k+1);if(k===sorted.length-1||sorted[k+1]<=threshold)break;}return v.map(x=>Math.max(0,x-threshold));};
 let w=Array(q).fill(0),z=[...w],t=1,converged=false;
 for(let step=0;step<5000;step++){const next=project(z.map((v,a)=>{let g=penalty[a]-rhs[a];for(let b=0;b<q;b++)g+=gram[a][b]*z[b];return v-g/L;})),delta=Math.max(...next.map((x,i)=>Math.abs(x-w[i]))),nt=(1+Math.sqrt(1+4*t*t))/2;z=next.map((x,i)=>x+(t-1)/nt*(x-w[i]));w=next;t=nt;if(step>30&&delta<1e-9){converged=true;break;}}
 const keep=[wtIndex,...ids.filter((_,i)=>w[i]>1e-8)];const refit=simplexFit(keep.map(i=>columns[i]),observed,positionWeights),weights=Array(p).fill(0);keep.forEach((id,i)=>weights[id]=refit.weights[i]);return {weights,r2:refit.r2,rmse:refit.rmse,converged:converged&&refit.converged};
}
