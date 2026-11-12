import {fitGroupedIntensities} from './sanger-grouped-mixture';
import type {RegionalLengthWindow} from './sanger-regional-length-fit';
/** Estimate residual detector drift jointly with invariant length spectra.
 * The same smooth response must explain an unedited prefix and every flank.
 * Only training blocks fit it; held-out flank photons remain validation data. */
export function calibrateRegionalResponse(windows:RegionalLengthWindow[],positions:number[][],shares:number[],delta:number,prefix:{positions:number[];columns:number[][];observed:number[]}){
 const H=shares.length,center=prefix.positions.reduce((a,b)=>a+b,0)/Math.max(1,prefix.positions.length),span=Math.max(1,...positions.flat())-center,theta=Array(2*Math.max(0,H-1)).fill(0),x=(p:number)=>(p-center)/span,at=(h:number,p:number,params=theta)=>h===H-1?1:Math.exp(params[2*h]+params[2*h+1]*x(p)),train=windows.map((w,t)=>{const ids=positions[t].flatMap((p,i)=>Math.floor(p/5)%3!==2?[i]:[]);return{classes:w.classes.map(c=>({...c,column:ids.flatMap(i=>c.column.slice(4*i,4*i+4))})),observed:ids.flatMap(i=>w.observed.slice(4*i,4*i+4)),positions:ids.map(i=>positions[t][i])};});
 let weights:number[][]=[];
 const phase=(w:RegionalLengthWindow,ps:number[],params:number[])=>{const fit=fitGroupedIntensities(w.classes.map(c=>c.column.map((v,k)=>v*at(c.reference,ps[k>>2],params))),w.observed,w.classes.map(c=>c.reference),shares,delta);if(!fit.converged)throw Error('Joint regional phase calibration did not converge.');return fit;};
 const signals=()=>{const rows:{p:number;y:number[];v:number[][]}[]=[];for(let i=0;i<prefix.positions.length;i++)rows.push({p:prefix.positions[i],y:prefix.observed.slice(4*i,4*i+4),v:shares.map((s,h)=>prefix.columns[h].slice(4*i,4*i+4).map(v=>v*s))});for(let t=0;t<train.length;t++){const w=train[t];for(let i=0;i<w.positions.length;i++)rows.push({p:w.positions[i],y:w.observed.slice(4*i,4*i+4),v:shares.map((_,h)=>[0,1,2,3].map(b=>w.classes.reduce((s,c,j)=>s+(c.reference===h?weights[t][j]*c.column[4*i+b]:0),0)))});}return rows;};
 const predict=(r:ReturnType<typeof signals>[number],params:number[])=>{const v=r.v.map((a,h)=>a.map(v=>v*at(h,r.p,params))),total=v.reduce((s,a)=>s+a.reduce((a,b)=>a+b,0),0),y=[0,1,2,3].map(b=>v.reduce((s,a)=>s+a[b],0)/Math.max(1e-30,total));return{v,total,y};};
 const loss=(rows:ReturnType<typeof signals>,params:number[])=>rows.reduce((s,r)=>s+predict(r,params).y.reduce((s,v,b)=>{const e=Math.abs(v-r.y[b]);return s+(e<=delta?e*e:2*delta*e-delta*delta);},0),0);
 let converged=theta.length===0,iterations=0;
 for(;iterations<80;iterations++){
  weights=train.map(w=>phase(w,w.positions,theta).weights);if(!theta.length)break;const rows=signals(),before=loss(rows,theta),n=theta.length,a=Array.from({length:n},(_,i)=>Array.from({length:n+1},(_,j)=>i===j?1e-8:0));
  for(const r of rows){const p=predict(r,theta);for(let b=0;b<4;b++){const e=p.y[b]-r.y[b],weight=Math.min(1,delta/Math.max(1e-30,Math.abs(e))),jac=theta.map((_,j)=>(p.v[j>>1][b]-p.y[b]*p.v[j>>1].reduce((a,b)=>a+b,0))/Math.max(1e-30,p.total)*(j%2?x(r.p):1));for(let i=0;i<n;i++){a[i][n]-=weight*jac[i]*e;for(let j=0;j<n;j++)a[i][j]+=weight*jac[i]*jac[j];}}}
  for(let i=0;i<n;i++){let pivot=i;for(let j=i+1;j<n;j++)if(Math.abs(a[j][i])>Math.abs(a[pivot][i]))pivot=j;[a[i],a[pivot]]=[a[pivot],a[i]];const d=a[i][i];for(let j=i;j<=n;j++)a[i][j]/=d;for(let k=0;k<n;k++)if(k!==i){const d=a[k][i];for(let j=i;j<=n;j++)a[k][j]-=d*a[i][j];}}
  const step=a.map(r=>r[n]);let alpha=1,accepted=false;while(alpha>=1/65536){const next=theta.map((v,j)=>v+alpha*step[j]);if(loss(rows,next)<before){theta.splice(0,theta.length,...next);accepted=true;break;}alpha/=2;}
  if(!accepted||Math.max(...step.map(v=>Math.abs(alpha*v)))<1e-6){converged=true;break;}
 }
 if(!converged)throw Error('Joint detector drift calibration did not converge.');
 const fits=windows.map((w,t)=>phase(w,positions[t],theta)),trainingFits=train.map(w=>phase(w,w.positions,theta)),predictions=windows.map((w,t)=>{const out:number[]=[];for(let k=0;k<w.observed.length;k+=4){const values=[0,1,2,3].map(b=>w.classes.reduce((s,c,j)=>s+trainingFits[t].weights[j]*c.column[k+b]*at(c.reference,positions[t][k>>2]),0)),sum=values.reduce((a,b)=>a+b,0);out.push(...values.map(v=>v/Math.max(1e-30,sum)));}return out;});
 return{at,parameters:theta.slice(),center,span,iterations,delta,lengths:windows.map((w,t)=>shares.map((_,h)=>w.classes.flatMap((c,j)=>c.reference===h?[{size:c.size,fraction:fits[t].weights[j]}]:[]))),trainingPredictions:predictions};
}
