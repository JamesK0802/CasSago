import {it,expect} from 'vitest';
import {conditionalRatioEnvelope} from '../../src/app/workers/core/sanger-conditional-ratio-envelope';
import type {LinearTerm} from '../../src/app/workers/core/sanger-marginal-lp';

it('never excludes any attainable ratio over positive-mass convex polytopes',()=>{
 let seed=1701;const random=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/4294967296;};
 for(let example=0;example<80;example++){
  const vertices=Array.from({length:6},()=>{const mass=.01+.98*random();return[mass*random(),mass];}),massRange:[number,number]=[Math.min(...vertices.map(v=>v[1])),Math.max(...vertices.map(v=>v[1]))],estimate=vertices.reduce((n,v)=>n+v[0],0)/vertices.reduce((n,v)=>n+v[1],0),oracle=(terms:LinearTerm[],maximum:boolean)=>{const values=vertices.map(v=>terms.reduce((n,[i,w])=>n+w*v[i],0));return maximum?Math.max(...values):Math.min(...values);};
  const[lo,hi]=conditionalRatioEnvelope([[0,1]],[[0,1]],[[1,1]],massRange,estimate,[0,1],oracle);
  // Linear-fractional extrema on a positive-denominator convex hull occur
  // at vertices. The oracle is independent of the implementation's algebra.
  for(const[n,m]of vertices){expect(lo).toBeLessThanOrEqual(n/m+1e-12);expect(hi).toBeGreaterThanOrEqual(n/m-1e-12);}
  expect(lo).toBeLessThanOrEqual(estimate);expect(hi).toBeGreaterThanOrEqual(estimate);
 }
});

it('retains every unobserved marker/base coupling allowed by Frechet bounds',()=>{
 const vertices=[[.65,.3,1],[.25,.45,.7],[.4,.2,.8]],massRange:[number,number]=[.2,.45],oracle=(terms:LinearTerm[],maximum:boolean)=>{const values=vertices.map(v=>terms.reduce((n,[i,w])=>n+w*v[i],0));return maximum?Math.max(...values):Math.min(...values);},[lo,hi]=conditionalRatioEnvelope([[0,1],[1,1],[2,-1]],[[0,1]],[[1,1]],massRange,.5,[0,1],oracle);
 for(const[f,m,g]of vertices){expect(lo).toBeLessThanOrEqual(Math.max(0,f+m-g)/m+1e-12);expect(hi).toBeGreaterThanOrEqual(Math.min(f,m)/m-1e-12);}
});
