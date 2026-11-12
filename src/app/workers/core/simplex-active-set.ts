/** Equality-constrained active-set QP on the probability simplex. Singular
 * tangent systems retain free coordinates at their current values. A dual gap
 * independently checks the returned solution; callers can retain a fallback. */
export function activeSimplex(G:ArrayLike<number>[],observations:number,tolerance=1e-13,initial?:number[],maxSteps?:number,ridge=0){
 const n=G.length,w=initial?.length===n&&initial.every(x=>Number.isFinite(x)&&x>=0)?[...initial]:Array(n).fill(0);if(!w.some(x=>x>0)){let v=0;for(let i=1;i<n;i++)if(G[i][i]<G[v][v])v=i;w[v]=1;}else{const sum=w.reduce((a,b)=>a+b,0);for(let i=0;i<n;i++)w[i]/=sum;}
 let active=w.flatMap((x,i)=>x>0?[i]:[]),gap=Infinity,steps=0;
 const gradient=()=>{const out=Array(n).fill(0);for(let i=0;i<n;i++){let value=0;for(let j=0;j<n;j++)value+=G[i][j]*w[j];out[i]=value;}return out;};
 const solve=(matrix:number[][],rhs:number[])=>{
  const m=rhs.length,rows=matrix.map((r,i)=>[...r,rhs[i]]),pivots:number[]=[];let scale=1;for(const row of matrix)for(const value of row)scale=Math.max(scale,Math.abs(value));const epsilon=1e-12*scale;let at=0;
  for(let c=0;c<m&&at<m;c++){let pivot=at;for(let r=at+1;r<m;r++)if(Math.abs(rows[r][c])>Math.abs(rows[pivot][c]))pivot=r;if(Math.abs(rows[pivot][c])<=epsilon)continue;[rows[at],rows[pivot]]=[rows[pivot],rows[at]];for(let r=at+1;r<m;r++){const f=rows[r][c]/rows[at][c];for(let k=c;k<=m;k++)rows[r][k]-=f*rows[at][k];}pivots.push(c);at++;}
  const x=Array(m).fill(0);for(let r=at-1;r>=0;r--){const c=pivots[r];let value=rows[r][m];for(let j=c+1;j<m;j++)value-=rows[r][j]*x[j];x[c]=value/rows[r][c];}return x;
 };
 for(;steps<(maxSteps??Math.max(100,4*n*n));steps++){
  const g=gradient(),lambda=w.reduce((s,x,i)=>s+x*g[i],0);let enter=0;for(let i=1;i<n;i++)if(g[i]<g[enter])enter=i;gap=Math.max(0,lambda-g[enter])/observations;if(gap<=tolerance)break;if(!active.includes(enter))active.push(enter);
  let moved=false;
  for(let inner=0;inner<=n;inner++){
   if(active.length<2)break;const pivot=active[0],rest=active.slice(1),H=rest.map(i=>rest.map(j=>G[i][j]-G[i][pivot]-G[pivot][j]+G[pivot][pivot]+(i===j?ridge:0))),rhs=rest.map(i=>g[pivot]-g[i]),d=solve(H,rhs),direction=Array(n).fill(0);rest.forEach((id,i)=>direction[id]=d[i]);direction[pivot]=-d.reduce((a,b)=>a+b,0);
   if(Math.max(...direction.map(Math.abs))<1e-15)break;let alpha=1;for(const id of active)if(direction[id]<0)alpha=Math.min(alpha,-w[id]/direction[id]);for(const id of active)w[id]=Math.max(0,w[id]+alpha*direction[id]);moved||=alpha>0;
   active=active.filter(id=>w[id]>1e-14);if(alpha>=1-1e-12)break;
   // Recompute the gradient before another constrained step.
   const updated=gradient();for(let i=0;i<n;i++)g[i]=updated[i];
  }
  if(!moved)break;
 }
 const g=gradient(),lambda=w.reduce((s,x,i)=>s+x*g[i],0);gap=Math.max(0,lambda-Math.min(...g))/observations;
 return{weights:w,converged:gap<=tolerance,optimalityGap:gap,steps};
}
