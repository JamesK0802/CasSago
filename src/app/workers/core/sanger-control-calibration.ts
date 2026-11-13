import {simplexFit} from './sanger-signal';

/** Fit reference abundance and dye mixing against the unedited control only.
 * Both blocks minimize the same squared photon-fraction residual on simplices;
 * estimating abundance against ideal one-hot dyes first biases mixed controls.
 * No sample signal, edit hypotheses or truth labels enter this calibration. */
export function calibrateControlMixture(sequences:string[],positions:number[],control:number[][]){
 const bases='ACGT',observed=positions.flatMap(p=>control[p]),codes=sequences.map(s=>positions.map(p=>bases.indexOf(s[p])));
 // Monochromatic WT cycles identify a dye spectrum independently of
 // abundance or reference-specific base yield. Mixed rows must not distort
 // these directly observed spectra to compensate for an incorrect yield.
 const pure=Array.from({length:4},(_,b)=>positions.filter((p,i)=>codes.every(code=>code[i]===b)&&control[p]?.every(v=>Number.isFinite(v)&&v>=0)&&control[p].reduce((a,b)=>a+b,0)>0)),direct=pure.map(ps=>ps.length?[0,1,2,3].map(d=>ps.reduce((n,p)=>n+control[p][d]/control[p].reduce((a,b)=>a+b,0),0)/ps.length):undefined);
 let means=Array.from({length:4},(_,b)=>Array.from({length:4},(_,k)=>(direct[b]?.[k]??Number(b===k)))),shares=Array(sequences.length).fill(1/sequences.length),previous=Infinity,converged=false,iterations=0;
 for(;iterations<500;iterations++){
  shares=simplexFit(codes.map(code=>code.flatMap(b=>means[b]||[0,0,0,0])),observed,undefined,shares).weights;
  const design=Array.from({length:4},(_,b)=>positions.map((_,i)=>codes.reduce((n,code,h)=>n+shares[h]*Number(code[i]===b),0)));
  // Exact row-simplex coordinate minimization preserves nonnegativity and a
  // unit total response, without a post-fit row normalization changing loss.
  for(let b=0;b<4;b++){
   if(direct[b]||!design[b].some(v=>v>0))continue;
   const residual=observed.map((v,k)=>v-design.reduce((n,col,a)=>n+(a===b?0:col[Math.floor(k/4)]*means[a][k%4]),0));
   means[b]=simplexFit(Array.from({length:4},(_,d)=>positions.flatMap((_,i)=>Array.from({length:4},(_,k)=>k===d?design[b][i]:0))),residual,undefined,means[b]).weights;
  }
  const loss=observed.reduce((n,v,k)=>n+(v-codes.reduce((s,code,h)=>s+shares[h]*(means[code[Math.floor(k/4)]]?.[k%4]||0),0))**2,0)/Math.max(1,observed.length);
  if(Math.abs(previous-loss)<=1e-12*Math.max(1,previous===Infinity?loss:previous)){converged=true;iterations++;break;}
  previous=loss;
 }
 return{shares,means,converged,iterations};
}
