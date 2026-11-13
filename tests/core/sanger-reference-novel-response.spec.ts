import {it,expect} from 'vitest';
import {referenceNovelResponses,measuredIntensityResponses} from '../../src/app/workers/core/sanger-intensity-response';
it('keeps inserted-base photon predictions invariant under a reference abundance/yield rescaling',()=>{
 const sequences=['ACGTACGTACGT','TGCATGCATGCA'],positions=Array.from({length:12},(_,i)=>i),shares=[.35,.65],next=[.6,.4],fallback=Array.from({length:4},(_,b)=>[0,1,2,3].map(d=>Number(d===b))),responses=sequences.map((s,h)=>[...s].map((base,p)=>[0,1,2,3].map(d=>Number(d==='ACGT'.indexOf(base))*(1+.5*h+.02*p))));
 const a=referenceNovelResponses(sequences,responses,positions,fallback),scaled=responses.map((rows,h)=>rows.map(row=>row.map(v=>v*shares[h]/next[h]))),b=referenceNovelResponses(sequences,scaled,positions,fallback);
 for(let h=0;h<2;h++)for(let base=0;base<4;base++)for(let d=0;d<4;d++)expect(next[h]*b[h][base][d]).toBeCloseTo(shares[h]*a[h][base][d],12);
});
it('retains each reference base yield instead of applying a pooled average to new bases',()=>{
 let seed=78124;const dna=()=>Array.from({length:800},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return 'ACGT'[seed>>>30];}).join(''),sequences=[dna(),dna()],refs=sequences.map(sequence=>({sequence})),shares=[.5,.5],yieldByReference=[[.8,1.1,2,.9],[1.4,.7,.5,1.2]],positions=Array.from({length:796},(_,i)=>i+2),means=Array.from({length:4},(_,b)=>[0,1,2,3].map(d=>Number(d===b))),photons=positions;
 const areas=Array.from({length:800},(_,p)=>[0,1,2,3].map(b=>sequences.reduce((n,s,h)=>n+shares[h]*yieldByReference[h][b]*Number(s[p]==='ACGT'[b]),0))),amplitudes=areas.map(row=>row.reduce((a,b)=>a+b,0)),vectors=areas.map((row,p)=>row.map(v=>v/amplitudes[p]));
 const response=measuredIntensityResponses(refs as any,shares,means,vectors,amplitudes,positions,1e-12);
 for(let h=0;h<2;h++)for(let b=0;b<4;b++){const strength=response.novelByReference[h][b].reduce((a,b)=>a+b,0);expect(strength).toBeCloseTo(yieldByReference[h][b],4);}
 expect(Math.abs(response.novel[2].reduce((a,b)=>a+b,0)-yieldByReference[0][2])).toBeGreaterThan(.1);
});
