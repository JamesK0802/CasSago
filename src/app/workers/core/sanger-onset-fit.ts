/** Ordered deletion onsets can be parameterised by a cumulative distribution.
 * Adjacent onset templates differ at only a few junction cycles. Optimising
 * their CDF with isotonic projection avoids the ill-conditioned dense Gram
 * matrix of almost-identical whole-read templates. This solves the same LS
 * problem, retaining every observation and both deletion directions. */
export function fitOrderedOnsets(columns:number[][],observed:number[],positionWeights:number[]|undefined,initial:number[]):number[]|null {
 const count=columns.length,n=observed.length;if(count<2)return [1];
 const normalized=columns.map(c=>[...c]);
 for(let k=0;k<n;k+=4){const total=columns[0].slice(k,k+4).reduce((s,v)=>s+v,0);if(total<=0)return null;for(let j=0;j<count;j++){const sum=columns[j].slice(k,k+4).reduce((s,v)=>s+v,0);if(Math.abs(sum-total)>1e-8*Math.max(1,total))return null;for(let b=0;b<4;b++)normalized[j][k+b]/=total;}}
 const differences=normalized.slice(0,-1).map((c,j)=>{const indices:number[]=[],values:number[]=[];for(let k=0;k<n;k++){const v=c[k]-normalized[j+1][k];if(Math.abs(v)>1e-14){indices.push(k);values.push(v);}}return {indices,values};});
 const rowSum=new Float64Array(n);for(const d of differences)d.indices.forEach((k,i)=>rowSum[k]+=Math.abs(d.values[i]));
 const weight=(k:number)=>positionWeights?.[Math.floor(k/4)]??1;
 const L=Math.max(...differences.map(d=>d.indices.reduce((s,k,i)=>s+Math.abs(d.values[i])*rowSum[k]*weight(k),0)));if(!L)return [...initial];
 const project=(v:number[])=>{const blocks:{lo:number;hi:number;sum:number;count:number}[]=[];for(let i=0;i<v.length;i++){blocks.push({lo:i,hi:i,sum:v[i],count:1});while(blocks.length>1){const b=blocks[blocks.length-1],a=blocks[blocks.length-2];if(a.sum/a.count<=b.sum/b.count)break;blocks.pop();blocks.pop();blocks.push({lo:a.lo,hi:b.hi,sum:a.sum+b.sum,count:a.count+b.count});}}const out=Array(v.length).fill(0);for(const b of blocks){const value=Math.max(0,Math.min(1,b.sum/b.count));for(let i=b.lo;i<=b.hi;i++)out[i]=value;}return out;};
 let sum=0,cdf=initial.slice(0,-1).map(v=>(sum+=v)),y=[...cdf],t=1;
 for(let iteration=0;iteration<2000;iteration++){
  const residual=normalized[count-1].map((v,k)=>v-observed[k]);differences.forEach((d,j)=>d.indices.forEach((k,i)=>residual[k]+=y[j]*d.values[i]));
  const gradient=differences.map(d=>d.indices.reduce((s,k,i)=>s+d.values[i]*residual[k]*weight(k),0)),next=project(y.map((v,j)=>v-gradient[j]/L)),delta=Math.max(...next.map((v,j)=>Math.abs(v-cdf[j])));
  if(iteration>20&&delta<1e-9){cdf=next;break;}
  const restart=next.reduce((s,v,j)=>s+(v-cdf[j])*(y[j]-v),0)>0,nt=restart?1:(1+Math.sqrt(1+4*t*t))/2;y=restart?[...next]:next.map((v,j)=>v+(t-1)/nt*(v-cdf[j]));cdf=next;t=nt;
 }
 const weights=[cdf[0],...cdf.slice(1).map((v,j)=>v-cdf[j]),1-cdf[cdf.length-1]].map(v=>Math.max(0,v)),mass=weights.reduce((s,v)=>s+v,0);return weights.map(v=>v/mass);
}
