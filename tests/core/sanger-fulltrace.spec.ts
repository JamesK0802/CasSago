import {describe,it,expect} from 'vitest';import {prepare} from '../../src/app/workers/core/sanger-path-common';import {materialize} from '../../src/app/workers/core/sanger-path-common';import {SangerAnalysisSettings} from '../../src/app/models/sanger-analysis';
describe('retain flanking evidence for pooled Sanger attribution',()=>{
 it('uses discriminating downstream bases while retaining the requested editing window',()=>{
  let seed=903517;const random=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/4294967296;},primer='ACTGACGATCGTAGCTACGA',a=Array.from({length:1050},()=> 'ACGT'[Math.floor(random()*4)]).join(''),b=Array.from(a,(base,p)=>[180,190,200,460].includes(p)?'ACGT'[('ACGT'.indexOf(base)+1)%4]:base).join(''),edit=[{start:330,end:330,inserted:'A'}],ea=materialize(a,edit)!,eb=materialize(b,edit)!,difference=(end:number)=>Array.from({length:end-119},(_,i)=>i+119).reduce((sum,p)=>sum+'ACGT'.split('').reduce((s,base)=>s+Math.abs(Number(a[p]===base)+Number(eb[p]===base)-Number(ea[p]===base)-Number(b[p]===base)),0),0);
  // WT A + edited B and edited A + WT B have identical marginal signals in
  // the short region, even though the WT references themselves are distinct.
  expect(difference(375)).toBe(0);expect(difference(700)).toBeGreaterThan(0);
  const settings={method:'mixture',signalModel:'calibrated',readStart:35,readEnd:700,window:90,maxIndel:5,largeDeletions:true,maxDeletion:300,maxInsertion:30,editStartRadius:5} as SangerAnalysisSettings;
  const prepared=prepare([{name:'A',sequence:primer+a,cut:primer.length+330},{name:'B',sequence:primer+b,cut:primer.length+330}],primer,settings);
  expect(prepared.report.positions).toContain(460);expect(prepared.report.amplicons[0].windowEnd).toBe(375);expect(prepared.report.positions.at(-1)).toBe(699);
 });
});
