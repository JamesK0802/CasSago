import {it,expect} from 'vitest';
import {groupedLinearQr} from '../../src/app/workers/core/sanger-grouped-linear-qr';
import {legacyGroupedLinear} from '../probes/experiments/sanger-linear-coordinate-oracle';
it('matches an independent coordinate-descent objective on noisy fixed-mass mixtures',()=>{
 let state=71;const random=()=>((state=Math.imul(state,1664525)+1013904223>>>0)/4294967296);
 for(let k=0;k<16;k++){const columns=Array.from({length:8},()=>Array.from({length:30},random)),y=Array.from({length:30},random),labels=[0,0,0,0,1,1,1,1],masses=[.6,.4],fit=groupedLinearQr(columns,y,labels,masses),old=legacyGroupedLinear(columns,y,labels,masses);expect(fit.converged).toBe(true);expect(fit.rmse).toBeCloseTo(old.rmse,8);expect(fit.weights.slice(0,4).reduce((a,b)=>a+b,0)).toBeCloseTo(.6,12);}
});
it('keeps exact small reference masses and all duplicate candidates available',()=>{
 const columns=[[1,0],[0,1],[.5,.5],[.5,.5],[.3,.7]],fit=groupedLinearQr(columns,[.4,.6],[0,0,0,0,1],[1,1e-18]);expect(fit.converged).toBe(true);expect(fit.weights[4]).toBe(1e-18);expect(fit.weights).toHaveLength(5);expect(fit.rmse).toBeLessThan(1e-10);
});
