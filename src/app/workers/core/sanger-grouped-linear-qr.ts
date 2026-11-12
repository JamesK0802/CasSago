import type {SimplexFit} from './sanger-signal';
export const groupedLinearDiagnostics:{observe?:(event:Record<string,unknown>)=>void}={};
/** Nonnegative least squares with fixed group sums. QR spans the active
 * affine face, retaining every candidate for global KKT pricing. Dependent
 * directions are only a numerical basis choice; no allele is pruned. */
export function groupedLinearQr(columns:number[][],observed:number[],labels:number[],masses:number[],initial?:number[]):SimplexFit{
 const n=columns.length,m=observed.length,groups=masses.map((_,g)=>labels.flatMap((v,i)=>v===g?[i]:[])),w=Array(n).fill(0);
 for(let g=0;g<groups.length;g++){const ids=groups[g],total=ids.reduce((s,i)=>s+Math.max(0,initial?.[i]||0),0);if(total){for(const i of ids)w[i]=Math.max(0,initial![i])*masses[g]/total;}else if(ids.length){let best=ids[0],loss=Infinity;for(const i of ids){let error=0;for(let k=0;k<m;k++)error+=(masses[g]*columns[i][k]-observed[k])**2;if(error<loss){loss=error;best=i;}}w[best]=masses[g];}}
 const active=new Set(w.flatMap((v,i)=>v>0?[i]:[]));let gap=Infinity;
 const residual=()=>{const r=Float64Array.from(observed);for(let i=0;i<n;i++)if(w[i])for(let k=0;k<m;k++)r[k]-=w[i]*columns[i][k];return r;};
 for(let iteration=0;iteration<300000;iteration++){
  if(iteration%500===0)groupedLinearDiagnostics.observe?.({fit:'linear-qr',iteration,columns:n,rows:m,active:active.size,gap});
  const r=residual(),anchors=groups.map(ids=>ids.filter(i=>active.has(i)).reduce((a,i)=>a<0||w[i]>w[a]?i:a,-1)),free=[...active].filter(i=>i!==anchors[labels[i]]),vectors=free.map(i=>Float64Array.from(columns[i],(v,k)=>v-columns[anchors[labels[i]]][k])),remaining=new Set(free.map((_,i)=>i)),order:number[]=[],triangular:Float64Array[]=[],rhs:number[]=[],norm=(v:Float64Array)=>v.reduce((s,x)=>s+x*x,0),largest=vectors.reduce((s,v)=>Math.max(s,norm(v)),0);
  while(remaining.size){let chosen=-1,squared=0;for(const j of remaining){const value=norm(vectors[j]);if(value>squared){squared=value;chosen=j;}}if(chosen<0||squared<=largest*1e-26)break;remaining.delete(chosen);const q=vectors[chosen],length=Math.sqrt(squared);for(let k=0;k<m;k++)q[k]/=length;const row=new Float64Array(free.length);row[chosen]=length;order.push(chosen);rhs.push(q.reduce((s,v,k)=>s+v*r[k],0));for(const j of remaining){let coefficient=0;for(let pass=0;pass<2;pass++){const projection=q.reduce((s,v,k)=>s+v*vectors[j][k],0);coefficient+=projection;for(let k=0;k<m;k++)vectors[j][k]-=projection*q[k];}row[j]=coefficient;}triangular.push(row);}
  const delta=Array(n).fill(0);for(let k=order.length-1;k>=0;k--){let value=rhs[k];for(let j=k+1;j<order.length;j++)value-=triangular[k][order[j]]*delta[free[order[j]]];delta[free[order[k]]]=value/triangular[k][order[k]];}
  for(const i of free)delta[anchors[labels[i]]]-=delta[i];let step=1;for(const i of active)if(delta[i]<0)step=Math.min(step,w[i]/-delta[i]);
  const blockers=new Set([...active].filter(i=>delta[i]<0&&w[i]/-delta[i]===step));
  let moved=0;for(const i of active){const change=step*delta[i];w[i]=Math.max(0,w[i]+change);moved=Math.max(moved,Math.abs(change));}
  if(step<1){for(const i of [...active])if(blockers.has(i)||w[i]===0){w[i]=0;active.delete(i);}if(moved>1e-15)continue;}
  const error=residual(),gradient=columns.map(c=>-c.reduce((s,v,k)=>s+v*error[k],0));gap=0;let enter=-1,leave=-1,bestGain=0;
  for(let g=0;g<groups.length;g++){const ids=groups[g];if(!ids.length||!masses[g])continue;const a=ids.reduce((a,i)=>gradient[i]<gradient[a]?i:a,ids[0]),positive=ids.filter(i=>w[i]>0);if(!positive.length)continue;const b=positive.reduce((a,i)=>gradient[i]>gradient[a]?i:a,positive[0]);gap+=ids.reduce((s,i)=>s+w[i]*gradient[i],0)-masses[g]*gradient[a];if(gradient[b]-gradient[a]>bestGain){bestGain=gradient[b]-gradient[a];enter=a;leave=b;}}
  gap=Math.max(0,gap)/Math.max(1,m);if(gap<=1e-14)break;
  if(enter<0||leave<0)break;if(!active.has(enter)){active.add(enter);continue;}
  // Correct residual roundoff on a degenerate active face without introducing
  // a regularizer or removing dependent candidate columns.
  let curve=0;for(let k=0;k<m;k++)curve+=(columns[enter][k]-columns[leave][k])**2;const amount=Math.min(w[leave],curve?bestGain/curve:w[leave]);if(!(amount>0))break;w[enter]+=amount;w[leave]=Math.max(0,w[leave]-amount);if(!w[leave])active.delete(leave);
 }
 const error=residual(),mean=observed.reduce((a,b)=>a+b,0)/Math.max(1,m),squared=error.reduce((s,v)=>s+v*v,0),variance=observed.reduce((s,v)=>s+(v-mean)**2,0);return{weights:w,rmse:Math.sqrt(squared/Math.max(1,m)),r2:variance?1-squared/variance:1,converged:gap<=1e-10,optimalityGap:gap};
}
