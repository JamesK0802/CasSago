import {describe,it,expect} from 'vitest';
import {fitGroupedIntensities,fitGroupedLinear} from '../../src/app/workers/core/sanger-grouped-mixture';
describe('reference constrained physical mixtures',()=>{
 it('retrieves molecule fractions when brighter molecules contribute more photons',()=>{
  const fit=fitGroupedIntensities([[4,0,0,0],[0,1,0,0]],[1.2/1.9,.7/1.9,0,0],[0,0],[1]);
  expect(fit.converged).toBe(true);expect(fit.weights[0]).toBeCloseTo(.3,5);expect(fit.rmse).toBeLessThan(1e-6);
 });
 it('honors reference shares during fitting even with overlapping emissions',()=>{
  const columns=[[1,0,0,0],[0,1,0,0],[1,0,0,0],[0,0,1,0]],labels=[0,0,1,1];
  for(const fit of [fitGroupedLinear(columns,[.4,.3,.3,0],labels,[.6,.4]),fitGroupedIntensities(columns,[.4,.3,.3,0],labels,[.6,.4])]){
   expect(fit.converged).toBe(true);expect(fit.weights[0]+fit.weights[1]).toBeCloseTo(.6,10);expect(fit.weights[2]+fit.weights[3]).toBeCloseTo(.4,10);expect(fit.weights[1]).toBeCloseTo(.3,5);expect(fit.weights[3]).toBeCloseTo(.3,5);
  }
 });
});
it('limits the influence of a corrupted cycle while reporting its full unweighted signal error',()=>{
 const columns=[[1,0,0,0],[0,1,0,0]].map(v=>Array.from({length:40},()=>v).flat()),observed=Array.from({length:40},(_,p)=>p?[.3,.7,0,0]:[1,0,0,0]).flat(),snapshot=JSON.stringify(observed),ordinary=fitGroupedIntensities(columns,observed,[0,0],[1]),robust=fitGroupedIntensities(columns,observed,[0,0],[1],.01);
 expect(robust.converged).toBe(true);expect(Math.abs(robust.weights[0]-.3)).toBeLessThan(Math.abs(ordinary.weights[0]-.3));expect(robust.rmse).toBeGreaterThan(.05);expect(JSON.stringify(observed)).toBe(snapshot);
});
it('uses measured control contribution information only when the spectrum cannot distinguish references',()=>{
 const repeated=Array.from({length:20},()=>[1,0,0,0]).flat(),ambiguous=fitGroupedIntensities([repeated,repeated],repeated,[0,0],[1],Infinity,undefined,{groups:[0,1],expected:[.4,.6],weight:10});expect(ambiguous.converged).toBe(true);expect(ambiguous.weights[0]).toBeCloseTo(.4,6);
 const a=Array.from({length:100},()=>[1,0,0,0]).flat(),b=Array.from({length:100},()=>[0,1,0,0]).flat(),y=Array.from({length:100},()=>[.65,.35,0,0]).flat(),informed=fitGroupedIntensities([a,b],y,[0,0],[1],Infinity,undefined,{groups:[0,1],expected:[.5,.5],weight:.01});expect(informed.converged).toBe(true);expect(informed.weights[0]).toBeCloseTo(.65,3);
});
