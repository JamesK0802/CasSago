/** Nonnegative least squares for the small local convolution systems. Cyclic
 * exact coordinate minimization; the projected gradient is the stop criterion. */
export function nonnegativeLeastSquares(columns:number[][],observed:number[],maxIterations=10000){
 const n=columns.length;if(!n||columns.some(c=>c.length!==observed.length))throw Error('Invalid local least-squares dimensions.');
 const gram=columns.map(a=>columns.map(b=>a.reduce((s,v,k)=>s+v*b[k],0))),rhs=columns.map(c=>c.reduce((s,v,k)=>s+v*observed[k],0)),weights=Array(n).fill(0),gradient=rhs.map(x=>-x);let converged=false;
 for(let iteration=0;iteration<maxIterations;iteration++){for(let i=0;i<n;i++){if(gram[i][i]<=0)continue;const next=Math.max(0,weights[i]-gradient[i]/gram[i][i]),delta=next-weights[i];if(delta){weights[i]=next;for(let j=0;j<n;j++)gradient[j]+=delta*gram[j][i];}}const violation=Math.max(...gradient.map((g,i)=>weights[i]>0?Math.abs(g):Math.max(0,-g)));if(violation<1e-10){converged=true;break;}}
 const residual=observed.map((v,k)=>v-columns.reduce((s,c,i)=>s+c[k]*weights[i],0));return{weights,rmse:Math.sqrt(residual.reduce((s,v)=>s+v*v,0)/observed.length),converged};
}
