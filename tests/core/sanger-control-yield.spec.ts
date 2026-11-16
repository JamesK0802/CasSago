import {it,expect} from 'vitest';
import {measuredIntensityResponses} from '../../src/app/workers/core/sanger-intensity-response';
it('learns unresolved WT contributions from identifiable WT colours without losing measured photons',()=>{
 let seed=8231;const dna=()=>Array.from({length:500},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return 'ACGT'[seed>>>30];}).join(''),refs=[{sequence:dna()},{sequence:dna()}],shares=[.5,.5],yieldAt=(h:number,p:number)=>Math.exp((h?.3:-.2)+.2*('ACGT'.indexOf(refs[h].sequence[p-1]||'A')-1.5)),control=Array.from({length:500},(_,p)=>[0,1,2,3].map(b=>refs.reduce((s,r,h)=>s+shares[h]*yieldAt(h,p)*Number(r.sequence[p]==='ACGT'[b]),0))),amplitudes=control.map(row=>row.reduce((a,b)=>a+b,0)),vectors=control.map((row,p)=>row.map(v=>v/amplitudes[p])),response=measuredIntensityResponses(refs as any,shares,[0,1,2,3].map(b=>[0,1,2,3].map(k=>Number(b===k))),vectors,amplitudes,Array.from({length:496},(_,p)=>p+2),1e-8);
 expect(response.yieldDiagnostic.status).toBe('WT-cross-validated');let error=0,count=0;for(let p=2;p<498;p++){for(let b=0;b<4;b++)expect(response.responses.reduce((s,v,h)=>s+shares[h]*v[p][b],0)).toBeCloseTo(control[p][b],8);if(refs[0].sequence[p]===refs[1].sequence[p]){const predicted=response.responses[0][p].reduce((a,b)=>a+b,0);error+=(predicted-yieldAt(0,p))**2;count++;}}
 expect(Math.sqrt(error/count)).toBeLessThan(.01);
});
