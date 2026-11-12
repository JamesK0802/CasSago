import {activeSimplex} from './simplex-active-set';
/** Accelerated projected gradient for a convex simplex quadratic. A simplex
 * dual gap, evaluated on the original unregularized objective, is mandatory.
 * Rank-deficient dictionaries need no artificial diagonal or inverse. */
export function mixtureQuadratic(gram:number[][],linear:number[],initial:number[],tolerance=1e-8){
 const n=linear.length,sum=initial.reduce((a,b)=>a+b,0);let x=sum?initial.map(v=>v/sum):Array.from({length:n},(_,i)=>Number(i===0)),y=[...x],t=1,gap=Infinity;
 const lipschitz=Math.max(1e-12,...gram.map(row=>row.reduce((s,v)=>s+Math.abs(v),0))),gradient=(v:number[])=>{const g=Array(n);for(let i=0;i<n;i++){let a=linear[i];for(let j=0;j<n;j++)a+=gram[i][j]*v[j];g[i]=a;}return g;};
 const project=(v:number[])=>{const sorted=[...v].sort((a,b)=>b-a);let total=0,theta=0;for(let i=0;i<n;i++){total+=sorted[i];const candidate=(total-1)/(i+1);if(i===n-1||sorted[i+1]<=candidate){theta=candidate;break;}}return v.map(a=>Math.max(0,a-theta));};
 let iterations=0;for(;iterations<20000;iterations++){
  if(iterations===100){const polished=activeSimplex(gram.map((row,i)=>row.map((v,j)=>v+linear[i]+linear[j])),1,tolerance,x,500,1e-10);if(polished.converged)return{...polished,iterations};x=polished.weights;y=[...x];t=1;}
  const g=gradient(y),next=project(y.map((v,i)=>v-g[i]/lipschitz));
  if(iterations%20===0){const q=gradient(next);gap=next.reduce((s,v,i)=>s+v*q[i],0)-Math.min(...q);if(gap<=tolerance){x=next;break;}}
  const u=(1+Math.sqrt(1+4*t*t))/2;y=next.map((v,i)=>v+(t-1)/u*(v-x[i]));x=next;t=u;
 }
 const g=gradient(x);gap=x.reduce((s,v,i)=>s+v*g[i],0)-Math.min(...g);return{weights:x,optimalityGap:gap,converged:gap<=tolerance,iterations};
}
