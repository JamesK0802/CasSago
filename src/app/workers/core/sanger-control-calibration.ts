import {simplexFit} from './sanger-signal';

/** Fit reference abundance and dye mixing against the unedited control only.
 * Both blocks minimize the same squared photon-fraction residual on simplices;
 * estimating abundance against ideal one-hot dyes first biases mixed controls.
 * No sample signal, edit hypotheses or truth labels enter this calibration. */
export function calibrateControlMixture(sequences:string[],positions:number[],control:number[][]){
 const bases='ACGT',observed=positions.flatMap(p=>control[p]),codes=sequences.map(s=>positions.map(p=>bases.indexOf(s[p])));
 let means=Array.from({length:4},(_,b)=>Array.from({length:4},(_,k)=>Number(b===k))),shares=Array(sequences.length).fill(1/sequences.length),previous=Infinity,converged=false,iterations=0;
 for(;iterations<500;iterations++){
  shares=simplexFit(codes.map(code=>code.flatMap(b=>means[b]||[0,0,0,0])),observed,undefined,shares).weights;
  const design=Array.from({length:4},(_,b)=>positions.map((_,i)=>codes.reduce((n,code,h)=>n+shares[h]*Number(code[i]===b),0)));
  // Exact row-simplex coordinate minimization preserves nonnegativity and a
  // unit total response, without a post-fit row normalization changing loss.
  for(let b=0;b<4;b++){
   if(!design[b].some(v=>v>0))continue;
   const residual=observed.map((v,k)=>v-design.reduce((n,col,a)=>n+(a===b?0:col[Math.floor(k/4)]*means[a][k%4]),0));
   means[b]=simplexFit(Array.from({length:4},(_,d)=>positions.flatMap((_,i)=>Array.from({length:4},(_,k)=>k===d?design[b][i]:0))),residual,undefined,means[b]).weights;
  }
  const loss=observed.reduce((n,v,k)=>n+(v-codes.reduce((s,code,h)=>s+shares[h]*(means[code[Math.floor(k/4)]]?.[k%4]||0),0))**2,0)/Math.max(1,observed.length);
  if(Math.abs(previous-loss)<=1e-12*Math.max(1,previous===Infinity?loss:previous)){converged=true;iterations++;break;}
  previous=loss;
 }
 return{shares,means,converged,iterations};
}
