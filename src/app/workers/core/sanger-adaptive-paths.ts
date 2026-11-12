import {Prepared,PathEdit,materialize} from './sanger-path-common';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import {SangerFitCache} from './sanger-fit-cache';
/** Column generation over single contiguous cut-local edits. Group count is
 * inferred by convex fitting, never supplied by the caller. Pricing minimizes
 * the residual independently at each inserted base, avoiding 4^length words
 * and 7^groups origin combinations. A computational budget is NOT convergence. */
export function adaptivePaths(ref:Prepared,families:{size:number;fraction:number}[],positions:number[],observed:number[],noise:number,budget:number,progress:(stage:string)=>void,thresholdPercent=0,optimizer?:SangerLinearOptimizer){
 const share=families.reduce((n,f)=>n+f.fraction,0);if(!(share>0))throw Error('No positive reference share.');
 const y=observed.map(v=>v/share),models:{events:PathEdit[];sequence:string}[]=[],columns:number[][]=[],seen=new Set<string>(),bases='ACGT';
 const add=(events:PathEdit[])=>{const sequence=materialize(ref.sequence,events);if(!sequence||seen.has(sequence))return false;seen.add(sequence);models.push({events,sequence});columns.push(positions.flatMap(p=>[...bases].map(b=>Number(sequence[p]===b))));return true;};
 add([]);
 const geometries:{start:number;removed:number;inserted:number}[]=[];
 for(const f of families)if(f.fraction>0)for(let start=ref.cuts[0]-3;start<=ref.cuts[0]+3;start++){
  if(f.size<0){if(start>=0&&start-f.size<=ref.sequence.length)add([{start,end:start-f.size,inserted:''}]);}
  else geometries.push({start,removed:f.size===0?10:0,inserted:f.size===0?10:f.size});
 }
 const cache=new SangerFitCache(columns,y),rowAt=new Map(positions.map((p,i)=>[p,i]));let iterations=0,gap=Infinity;
 let fit=cache.fit(columns.map((_,i)=>i));
 for(;iterations<budget;iterations++){
  if(!fit.converged)throw Error('Adaptive restricted coefficient fit did not converge.');
  const residual=y.map((v,k)=>columns.reduce((n,c,i)=>n+c[k]*fit.weights[i],0)-v);
  const costs=columns.map(c=>c.reduce((n,v,k)=>n+v*residual[k],0)),lambda=costs.reduce((n,c,i)=>n+c*fit.weights[i],0);
  let minimum=Math.min(...costs);const proposals:{events:PathEdit[];cost:number}[]=[];
  for(const g of geometries){
   const word=Array.from({length:g.inserted},(_,j)=>{const row=rowAt.get(g.start+j);if(row===undefined)return ref.sequence[g.start+j]||'A';let best=0;for(let b=1;b<4;b++)if(residual[4*row+b]<residual[4*row+best])best=b;return bases[best];}).join('');
   const events=[{start:g.start,end:g.start+g.removed,inserted:word}],sequence=materialize(ref.sequence,events);if(!sequence)continue;
   let cost=0;for(let i=0;i<positions.length;i++){const b=bases.indexOf(sequence[positions[i]]);if(b>=0)cost+=residual[4*i+b];}
   minimum=Math.min(minimum,cost);if(cost<lambda-1e-10*y.length&&!seen.has(sequence))proposals.push({events,cost});
  }
  gap=Math.max(0,lambda-minimum)/y.length;
  if(gap<=1e-10)break;
  proposals.sort((a,b)=>a.cost-b.cost);let added=0;for(const p of proposals)if(add(p.events))added++;
  if(!added)throw Error('Adaptive pricing stalled before convergence; unresolved signal.');
  progress(`Adaptive sequence discovery · iteration ${iterations+1} · ${models.length} candidate sequences`);
  fit=cache.fit(columns.map((_,i)=>i));
 }
 const converged=iterations<budget&&gap<=1e-10;
 if(!converged)throw Error(`Adaptive sequence discovery budget exhausted after ${iterations} iterations; gap ${gap}; candidates ${models.length}`);
 // Expand words from the fitted length-class marginal, rather than assuming
 // residual-pricing vertices include every threshold-feasible physical allele.
 // This uses measured/fitted emissions only; no benchmark truth is available.
 if(thresholdPercent>0&&optimizer){
  const weights=[...fit.weights],sizes=models.map(m=>m.events.reduce((n,e)=>n+e.inserted.length-e.end+e.start,0)),minimum=thresholdPercent/100;
  let generated=0;
  for(const g of geometries){
   const size=g.inserted-g.removed,own=models.flatMap((_,i)=>sizes[i]===size?[i]:[]);if(!own.length)continue;
   const choices=Array.from({length:g.inserted},(_,j)=>{const at=rowAt.get(g.start+j);if(at===undefined)return[ref.sequence[g.start+j]||'A'];const v=[0,0,0,0];for(const i of own)for(let b=0;b<4;b++)v[b]+=columns[i][4*at+b]*weights[i];const keep=[...bases].filter((_,b)=>v[b]>=minimum-Math.max(1e-5,noise/share));return keep.length?keep:[bases[v.indexOf(Math.max(...v))]];});
   const words=choices.reduce((n,c)=>n*c.length,1);if(!Number.isFinite(words)||generated+words>20000)throw Error('Measured word expansion budget exhausted; signal remains unresolved.');
   const visit=(at:number,word:string)=>{if(at===choices.length){generated++;add([{start:g.start,end:g.start+g.removed,inserted:word}]);return;}for(const b of choices[at])visit(at+1,word+b);};visit(0,'');
  }
  progress('Expanded measured length-class words · '+models.length+' candidates');
 }
 let weights=models.map((_,i)=>fit.weights[i]||0);
 if(thresholdPercent>0&&optimizer){
  // Select and refit support jointly. Deleting low fitted weights without a
  // refit destroys measured fluorescence and can discard true alleles whose
  // contribution was split across several equivalent proposals.
  const lower=thresholdPercent/100,lines=['Minimize',' obj: '+models.map((_,i)=>'z'+i).join(' + '),'Subject To'];let row=0;
  const expression=(values:number[])=>values.flatMap((v,i)=>v?[String(v)+' w'+i]:[]).join(' + ')||'0 w0';
  lines.push(' sum: '+models.map((_,i)=>'w'+i).join(' + ')+' = 1');
  let maximumResidual=0;for(let k=0;k<y.length;k++)maximumResidual=Math.max(maximumResidual,Math.abs(columns.reduce((n,c,i)=>n+c[k]*weights[i],0)-y[k]));
  const epsilon=maximumResidual+Math.max(1e-5,noise/share);
  const seenRows=new Set<string>();
  for(let k=0;k<y.length;k++){const values=columns.map(c=>c[k]),key=values.join('')+':'+y[k];if(seenRows.has(key))continue;seenRows.add(key);const expr=expression(values);lines.push(' lo'+row+': '+expr+' >= '+(y[k]-epsilon),' hi'+row+++': '+expr+' <= '+(y[k]+epsilon));}
  models.forEach((m,i)=>{lines.push(' up'+i+': w'+i+' - z'+i+' <= 0');if(m.events.length)lines.push(' min'+i+': w'+i+' - '+lower+' z'+i+' >= 0');});
  lines.push('Bounds',...models.map((_,i)=>' 0 <= w'+i+' <= 1'),'Binary',...models.map((_,i)=>' z'+i),'End');
  progress('Joint threshold-constrained candidate selection · '+models.length+' candidates');
  const answer=optimizer.solve(lines.join('\n'),{output_flag:false,time_limit:60,mip_rel_gap:0.01,primal_feasibility_tolerance:1e-9,mip_feasibility_tolerance:1e-9});
  const proposal=models.map((_,i)=>answer.Columns['w'+i]?.Primal??NaN);
  let error=0;for(let k=0;k<y.length;k++)error=Math.max(error,Math.abs(columns.reduce((n,c,i)=>n+c[k]*proposal[i],0)-y[k]));
  if(proposal.some((w,i)=>!Number.isFinite(w)||w< -1e-8||models[i].events.length&&w>1e-8&&w+1e-8<lower)||Math.abs(proposal.reduce((a,b)=>a+b,0)-1)>1e-7||!Number.isFinite(error)||error>epsilon+1e-7)throw Error('No signal-valid threshold-constrained support; candidate vocabulary requires expansion. '+answer.Status);
  weights=proposal.map(w=>Math.max(0,w));
 }
 const groups=models.flatMap((m,i)=>weights[i]>1e-9?[{events:m.events,fraction:weights[i]*share}]:[]);
 let maximum=0;for(let i=0;i<positions.length;i++){let e=0;for(let b=0;b<4;b++){const k=4*i+b,pred=columns.reduce((n,c,j)=>n+c[k]*weights[j],0);e+=(pred-y[k])**2;}maximum=Math.max(maximum,share*Math.sqrt(e/4));}
 return {options:[groups],scores:[{rmse:share*fit.rmse,maximum}],details:[{method:'adaptive-column-generation',iterations,candidates:models.length,groups:groups.length,optimalityGap:gap,converged,exhaustive:false}]};
}
