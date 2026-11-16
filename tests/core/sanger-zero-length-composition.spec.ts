import {describe,it,expect} from 'vitest';
import {zeroLengthCompositionBounds} from '../../src/app/workers/core/sanger-zero-length-composition';

describe('zero-length composition across incoming phases',()=>{
 it('retains WT support when upstream indels change read-cycle bases without local substitutions',()=>{
  // Two different upstream phases, same locally unmodified target sequence.
  const r=zeroLengthCompositionBounds([.3,.4],0,1);
  expect(r.unmodifiedRange[0]).toBeCloseTo(.3-1e-7,12);
  expect(r.substitutionRange[1]).toBeCloseTo(1e-7,12);
 });
 it('uses absolute molecular budget before within-reference normalization',()=>{
  const r=zeroLengthCompositionBounds([.075,.1],.005,.25);
  expect(r.unmodifiedRange[0]).toBeCloseTo((.075-.005-1e-7)/.25,12);
  expect(r.substitutionRange[1]).toBeCloseTo((.005+1e-7)/.25,12);
 });
 it('contains coherent mixtures with multiple phases and multiple substitutions per molecule',()=>{
  const molecules=[{mass:.08,phase:-3,changes:0},{mass:.04,phase:2,changes:0},{mass:.02,phase:0,changes:1},{mass:.01,phase:4,changes:3}];
  const share=.25,mass=molecules.reduce((s,m)=>s+m.mass,0),budget=molecules.reduce((s,m)=>s+m.mass*m.changes,0),r=zeroLengthCompositionBounds([mass,mass],budget,share);
  const wt=molecules.filter(m=>!m.changes).reduce((s,m)=>s+m.mass,0)/share,sub=molecules.filter(m=>m.changes).reduce((s,m)=>s+m.mass,0)/share;
  expect(wt).toBeGreaterThanOrEqual(r.unmodifiedRange[0]);expect(wt).toBeLessThanOrEqual(r.unmodifiedRange[1]);
  expect(sub).toBeGreaterThanOrEqual(r.substitutionRange[0]);expect(sub).toBeLessThanOrEqual(r.substitutionRange[1]);
 });
 it('keeps the full composition range when the existing budget permits it',()=>{
  expect(zeroLengthCompositionBounds([0,.2],.5,1)).toEqual({unmodifiedRange:[0,.2],substitutionRange:[0,.2]});
 });
});
