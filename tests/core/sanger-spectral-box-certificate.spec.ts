import {it,expect} from 'vitest';
import {spectralBoxCertificates} from '../../src/app/workers/core/sanger-spectral-box-certificate';
it('retrieves exact full-rank ellipsoid coordinate bounds without excluding any feasible point',()=>{
 const q={name:'paired mass',terms:[[[0,1],[1,-1]]as[number,number][]],upper:.04},fit=spectralBoxCertificates(q,[1,1],[{columns:[0,1],mass:1}]);expect(fit.rank).toBe(1);for(const b of fit.bounds){expect(b.lower).toBeCloseTo(.4,9);expect(b.upper).toBeCloseTo(.6,9);}for(let j=0;j<=20;j++){const x=[.4+.01*j,.6-.01*j];for(const b of fit.bounds){expect(x[b.column]).toBeGreaterThanOrEqual(b.lower-1e-10);expect(x[b.column]).toBeLessThanOrEqual(b.upper+1e-10);}}
});
it('retains indistinguishable states when the spectrum is rank deficient',()=>{
 const q={name:'unresolved',terms:[[[0,1],[1,1],[2,-1]]as[number,number][]],upper:.04},fit=spectralBoxCertificates(q,[1,1,1],[{columns:[0,1,2],mass:1}]);expect(fit.rank).toBe(1);for(let j=0;j<=20;j++)for(let k=0;k<=20-j;k++){const x=[j/20,k/20,1-(j+k)/20];if((x[0]+x[1]-x[2])**2<=q.upper+1e-12)for(const b of fit.bounds){expect(x[b.column]).toBeGreaterThanOrEqual(b.lower-1e-9);expect(x[b.column]).toBeLessThanOrEqual(b.upper+1e-9);}}
});
