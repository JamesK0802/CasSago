import {it,expect} from 'vitest';
import loadCommon from '../probes/load-common-optimizer';
import {analyzeAutomaticSanger} from '../../src/app/workers/core/sanger-automatic';
import {marginalModelDiagnostics} from '../../src/app/workers/core/sanger-marginal-engine';
import {writeSynthTrace} from '../../src/app/utils/sanger-synthetic-benchmark';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
it('keeps pre-cut emission phase-invariant and treats separated reference gaps as flank evidence',async()=>{
 let seed=9301;const dna=(n:number)=>Array.from({length:n},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return 'ACGT'[seed>>>30];}).join(''),primer='ACGTCGTAGCTAGGATCCGA',refs=[200,360].map((cut,i)=>({name:'ref'+i,sequence:primer+dna(650),cut:cut+primer.length,cuts:[cut+primer.length]})),trace=parseSangerTrace(writeSynthTrace(refs.map(r=>({sequence:r.sequence.slice(primer.length),weight:.5})),650,'control'));
 const stop=new Error('Original physical model inspected before optimization.');let checked=false,maximumPhaseDifference=Infinity,gapIsSpectral=false;
 marginalModelDiagnostics.observe=model=>{const q=280,m=model.measurements.find(m=>m.position===q)!;expect(m).toBeDefined();const columns=model.states[0][1].map(s=>s.column);maximumPhaseDifference=Math.max(...m.forms.map(form=>{const coefficients=new Map(form),values=columns.map(i=>coefficients.get(i)||0);return Math.max(...values)-Math.min(...values);}));gapIsSpectral=(model as any).spectralPositions?.includes(q)===true;checked=true;throw stop;};
 const optimizer=await loadCommon();
 try{expect(()=>analyzeAutomaticSanger({references:refs,primer,control:trace,sample:trace,maxDeletion:2,maxInsertion:2,indelThreshold:2,groupRecovery:true},()=>{},optimizer)).toThrow(stop);expect(checked).toBe(true);expect(maximumPhaseDifference).toBeLessThan(1e-12);expect(gapIsSpectral).toBe(true);}finally{marginalModelDiagnostics.observe=undefined;}
},120000);
