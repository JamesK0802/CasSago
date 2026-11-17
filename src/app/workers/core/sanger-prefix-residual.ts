import {mixPhase,phaseNeighbor} from './sanger-signal';
/** Held-out paired-prefix error on the exact phased photon observable used by
 * every main-model measurement. Rows are never used to fit the transfer. */
export function prefixPhotonVariance(positions:number[],raw:number[][],responses:number[][][],shares:number[],transfer:{apply(h:number,p:number,dye:number[]):number[]},lag:number,backgroundAt:(p:number)=>number[]=()=>[0,0,0,0]){
 const foreground=(p:number)=>[0,1,2,3].map(b=>responses.reduce((n,r,h)=>n+shares[h]*transfer.apply(h,p,r[p])[b],0));
 return positions.reduce((sum,p)=>{const prediction=mixPhase(foreground(p),foreground(phaseNeighbor(p,lag)),lag).map((v,b)=>v+backgroundAt(p)[b]),total=prediction.reduce((a,b)=>a+b,0);return sum+prediction.reduce((e,v,b)=>e+(v/Math.max(1e-30,total)-raw[p][b])**2,0);},0)/(4*positions.length);
}
