import {mixtureQuadratic} from './sanger-mixture-qp';
/** Working-set nonnegative regression. Geometry is independent of sample data;
 * junction observations are shared training rows. Independent flanks are split
 * into inner folds0/1 and an outer validation fold2. The reporting cutoff is
 * never used. Lazy Gram columns avoid a full dictionary-sized square matrix. */
export interface MeasuredGramCache {size:number;counts:number[];rhs:number[][];gram:Map<number,{size:number;rows:Float64Array[]}>;weights:Map<string,number[]>}
export const measuredGramCache=():MeasuredGramCache=>({size:0,counts:[0,0,0,0],rhs:Array.from({length:4},()=>[]),gram:new Map(),weights:new Map()});
export function measuredLasso(columns:number[][],observed:number[],costs:number[],noise:number,folds:number[],blocks:number[],selection=true,state=measuredGramCache()){
 const tolerance=Math.max(1e-7,noise*.001),n=columns.length,m=observed.length,counts=state.counts,rhs=state.rhs;
 if(!state.size)for(let k=0;k<m;k++)counts[folds[k]+1]++;
 for(let i=state.size;i<n;i++){for(let c=0;c<4;c++)rhs[c][i]=0;for(let k=0;k<m;k++)rhs[folds[k]+1][i]+=columns[i][k]*observed[k];}state.size=n;
 const ensure=(i:number)=>{const hit=state.gram.get(i);if(hit?.size===n)return hit.rows;const rows=Array.from({length:4},(_,c)=>{const a=new Float64Array(n);if(hit)a.set(hit.rows[c]);return a;});for(let j=hit?.size||0;j<n;j++){const mirror=state.gram.get(j);if(mirror&&mirror.size>i){for(let c=0;c<4;c++)rows[c][j]=mirror.rows[c][i];}else for(let k=0;k<m;k++)rows[folds[k]+1][j]+=columns[i][k]*columns[j][k];}state.gram.set(i,{size:n,rows});return rows;};
 const solve=(categories:number[],alphaScale:number,allowed?:number[])=>{
  const count=categories.reduce((n,c)=>n+counts[c],0),right=Float64Array.from({length:n},(_,i)=>categories.reduce((v,c)=>v+rhs[c][i],0)/count),alpha=Math.max(1e-10,noise*alphaScale),candidates=allowed||Array.from({length:n},(_,i)=>i),key=categories.join(',')+':'+alpha+':'+!!allowed,previous=state.weights.get(key)||[...state.weights.entries()].reverse().find(([k])=>k.endsWith(':'+!!allowed)&&k.startsWith(categories.join(',')+':'))?.[1]||[...state.weights.values()].at(-1)||[],weights=Array.from({length:n},(_,i)=>candidates.includes(i)?previous[i]||0:0),working=candidates.filter(i=>costs[i]===0||weights[i]>0),gram=new Map<number,Float64Array>();
  const column=(i:number)=>{let g=gram.get(i);if(!g){const raw=ensure(i);g=Float64Array.from({length:n},(_,j)=>categories.reduce((v,c)=>v+raw[c][j],0)/count);gram.set(i,g);}return g;};
  let optimalityGap=Infinity,iterations=0;
  for(;iterations<Math.max(100,n);iterations++){
   if(Math.abs(weights.reduce((a,b)=>a+b,0)-1)<1e-10){const g=Float64Array.from(right,(v,i)=>-v+alpha*costs[i]);for(const i of working)if(weights[i]>0){const c=column(i);for(const j of candidates)g[j]+=weights[i]*c[j];}const gap=working.reduce((v,i)=>v+weights[i]*g[i],0)-Math.min(...candidates.map(i=>g[i]));if(gap<=tolerance){optimalityGap=gap;break;}}
   const center=costs.lastIndexOf(0),g0=column(center),offset=g0[center],restricted=mixtureQuadratic(working.map(i=>working.map(j=>column(i)[j]-g0[i]-g0[j]+offset)),working.map(i=>-right[i]+alpha*costs[i]+g0[i]),working.map(i=>weights[i]),tolerance/10);weights.fill(0);working.forEach((id,j)=>weights[id]=restricted.weights[j]);
   const retained=working.filter(i=>weights[i]>1e-12);working.splice(0,working.length,...retained);
   const gradient=Float64Array.from(right,(v,i)=>-v+alpha*costs[i]);for(const i of working)if(weights[i]>0){const g=column(i);for(const j of candidates)gradient[j]+=weights[i]*g[j];}
   const lambda=working.reduce((v,i)=>v+weights[i]*gradient[i],0),ordered=[...candidates].sort((a,b)=>gradient[a]-gradient[b]);optimalityGap=lambda-gradient[ordered[0]];if(optimalityGap<=tolerance)break;
   const entering=ordered.filter(i=>!working.includes(i)&&gradient[i]<lambda-tolerance/2).slice(0,4);if(!entering.length)break;working.push(...entering);
  }
  state.weights.set(key,[...weights]);
  return{weights,converged:optimalityGap<=tolerance,optimalityGap,iterations};
 };
 if(!selection){const fit=solve([0,1,2],.1);return{fit,developmentFit:fit,foldFits:[fit,fit],nullFit:fit,alphaScale:.1,selection:[]};}
 const predict=(w:number[],k:number)=>columns.reduce((n,col,i)=>n+w[i]*col[k],0),options=[.03,.1,.3,1,3].map(scale=>{const a=solve([0,1],scale),b=solve([0,2],scale),losses=new Map<number,number[]>();for(let k=0;k<m;k++)if(folds[k]===0||folds[k]===1){const error=(observed[k]-predict(folds[k]===0?b.weights:a.weights,k))**2,g=losses.get(blocks[k])||[];g.push(error);losses.set(blocks[k],g);}const values=[...losses.values()].map(g=>g.reduce((a,b)=>a+b,0)/g.length),mean=values.reduce((a,b)=>a+b,0)/Math.max(1,values.length),se=Math.sqrt(values.reduce((n,v)=>n+(v-mean)**2,0)/Math.max(1,values.length-1)/Math.max(1,values.length));return{scale,a,b,mean,se};}),best=options.filter(x=>x.a.converged&&x.b.converged).sort((a,b)=>a.mean-b.mean)[0];
 if(!best)throw Error('Measured regression did not converge in inner validation.');
 // Choose regularization on development observations only. The one-SE
 // preference for a compact explanation cannot override the independent
 // control's noise bound by systematically underfitting the training signal.
 const eligible=options.filter(x=>x.a.converged&&x.b.converged&&x.mean<=best.mean+best.se).sort((a,b)=>b.scale-a.scale),limit=Math.max(.005**2,9*noise);
 let selected=best,developmentFit=solve([0,1,2],best.scale),developmentMse=Infinity;
 for(const candidate of eligible){const dev=solve([0,1,2],candidate.scale);let error=0,count=0;for(let k=0;k<m;k++)if(folds[k]!==2){error+=(observed[k]-predict(dev.weights,k))**2;count++;}if(dev.converged&&error/count<=limit){selected=candidate;developmentFit=dev;developmentMse=error/count;break;}}
 const fit=solve([0,1,2,3],selected.scale),nullFit=solve([0,1,2],selected.scale,costs.flatMap((c,i)=>c===0?[i]:[]));
 return{fit,developmentFit,developmentMse,foldFits:[selected.a,selected.b],nullFit,alphaScale:selected.scale,selection:options.map(x=>({scale:x.scale,mean:x.mean,se:x.se}))};
}
