import type {Prepared} from './sanger-path-common';
import {nonnegativeLeastSquares} from './sanger-nonnegative-ls';
import {measuredControlKernels} from './sanger-control-kernels';
import {controlReferenceYields} from './sanger-control-yield';
/** Control-only fluorescence decomposition. At each cycle, known WT bases
 * partition the measured photons; coincident WT colours are allocated using
 * a separately cross-validated WT context-yield model. A single reference uses
 * identical construction. It never reads edited fractions or genotype truth. */
export function measuredIntensityResponses(refs:Prepared[],shares:number[],means:number[][],control:number[][],amplitudes:number[],positions:number[],quantum:number){
 const bare=control.map((_,p)=>Array.from({length:4},(_,b)=>refs.reduce((n,r,h)=>n+(r.sequence[p]==='ACGT'[b]?shares[h]:0),0))),weights=control.map((v,p)=>{const active=bare[p].flatMap((v,b)=>v>0?[b]:[]);if(!active.length)return[0,0,0,0];const fit=nonnegativeLeastSquares(active.map(b=>means[b]),v),total=fit.weights.reduce((n,v)=>n+v,0),row=[0,0,0,0];active.forEach((b,i)=>row[b]=total?fit.weights[i]/total:bare[p][b]);const count=amplitudes[p]/Math.max(1e-30,quantum);for(let b=0;b<4;b++)row[b]=(row[b]*count+bare[p][b])/(count+1);return row;}),kernels=measuredControlKernels(refs,shares,means,control,weights),yieldModel=controlReferenceYields(refs.map(r=>r.sequence),shares,amplitudes,weights,positions),responses=refs.map((ref,h)=>control.map((_,p)=>{const b='ACGT'.indexOf(ref.sequence[p]);if(b<0||!bare[p][b])return[0,0,0,0];const yields=refs.map((r,j)=>r.sequence[p]===ref.sequence[p]?shares[j]*yieldModel.at(j,p):0),den=yields.reduce((a,b)=>a+b,0),a=amplitudes[p]*weights[p][b]*yieldModel.at(h,p)/Math.max(1e-30,den);return kernels[p][b].map(v=>a*v);})),baseStrength=[0,1,2,3].map(b=>{const den=positions.reduce((n,p)=>n+bare[p][b],0);return den?positions.reduce((n,p)=>n+amplitudes[p]*weights[p][b],0)/den:1;}),novel=means.map((row,b)=>row.map(v=>v*baseStrength[b]));
 return{responses,novel,novelByReference:referenceNovelResponses(refs.map(r=>r.sequence),responses,positions,novel,yieldModel.cycleEnvelope),cycleEnvelope:yieldModel.cycleEnvelope,kernels,baseStrength,yieldDiagnostic:yieldModel.diagnostic,allocationBounds:refs.map((_,h)=>control.map((_,p)=>yieldModel.bounds(h,p)))};
}

/** New bases must use the same per-reference photon units as WT donors.
 * Remove the WT-learned read-cycle envelope before averaging; reapply it at
 * the observed cycle, so an indel does not transport detector decay.
 * Averaging the reference's control-derived responses preserves the physical
 * product share*response under an otherwise arbitrary abundance/yield scale.
 * Only WT controls are read; a missing base retains the existing fallback. */
export function referenceNovelResponses(sequences:string[],responses:number[][][],positions:number[],fallback:number[][],cycleEnvelope:(h:number,p:number)=>number=()=>1){
 return sequences.map((sequence,h)=>fallback.map((row,b)=>{
  const ps=positions.filter(p=>sequence[p]==='ACGT'[b]&&responses[h][p]?.some(v=>v>0));
  return ps.length?[0,1,2,3].map(d=>ps.reduce((n,p)=>n+responses[h][p][d]/cycleEnvelope(h,p),0)/ps.length):row.slice();
 }));
}
