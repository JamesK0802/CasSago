import {it,expect} from 'vitest';
import {measuredIntensityResponses} from '../../src/app/workers/core/sanger-intensity-response';
it('reconstructs paired control photons for single and mixed references with unequal brightness',()=>{
 const means=[[.97,.01,.01,.01],[.01,.97,.01,.01],[.01,.01,.97,.01],[.01,.01,.01,.97]];
 for(const refs of [[{sequence:'A'}],[{sequence:'A'},{sequence:'G'}]]){
  const shares=refs.length===1?[1]:[.6,.4],control=refs.length===1?[[.98,.005,.01,.005]]:[[.79,.01,.18,.02]],response=measuredIntensityResponses(refs as any,shares,means,control,[2],[0],1e-4);
  for(let b=0;b<4;b++)expect(response.responses.reduce((n,rows,h)=>n+shares[h]*rows[0][b],0)).toBeCloseTo(2*control[0][b],10);
  for(const rows of response.responses)expect(rows[0].every(v=>v>=0&&Number.isFinite(v))).toBe(true);
 }
});
