import '@angular/compiler';
import {it,expect} from 'vitest';
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import initializeConic,* as conic from '../../public/vendor/clarabel-cvxjs-0.1.4/clarabel.mjs';
async function loadOptimizer(){await initializeConic({module_or_path:readFileSync('public/vendor/clarabel-cvxjs-0.1.4/clarabel_wasm_bg.wasm')});const o=await loadHighs();o.conic=conic;return o;}
import {writeSynthTrace} from '../../src/app/utils/sanger-synthetic-benchmark';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
import {analyzeAutomaticSanger} from '../../src/app/workers/core/sanger-automatic';
import {marginalModelDiagnostics} from '../../src/app/workers/core/sanger-marginal-engine';
import {marginalViolation,MarginalLinearSolver} from '../../src/app/workers/core/sanger-marginal-lp';
const root=process.env.CASSAGO_SHORT_FLANK_AUDIT;
function fixture(){const primer='ACGTCGTAGCTAGGATCCGA',cuts=[143,184,225];let seed=83172;const dna=()=>Array.from({length:324},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return 'ACGT'[seed>>>30];}).join('');return{primer,amplicons:cuts.map((cut,h)=>({name:'Reference '+(h+1),sequence:primer+dna(),cuts:[primer.length+cut]}))};}
async function run(scramble:boolean){
 if(root)mkdirSync(root,{recursive:true});const t=fixture(),parts=t.amplicons.map((r:any)=>({sequence:r.sequence.slice(t.primer.length),weight:1/3})),n=324,control=parseSangerTrace(writeSynthTrace(parts,n,'qualified-control'),'control.ab1'),sample=scramble?parseSangerTrace(writeSynthTrace(parts.map((p:any)=>({...p,sequence:p.sequence.slice(0,275)+[...p.sequence.slice(275)].map((b:string)=>({A:'C',C:'G',G:'T',T:'A'}[b])).join('')})),n,'corrupted-tail'),'corrupted.ab1'):control,optimizer=await loadOptimizer(),stop=new Error('Verified blind representative observed'),stages:string[]=[],started=performance.now();let model:any,proof:any,captured:any;

 marginalModelDiagnostics.observe=m=>{model=m;};marginalModelDiagnostics.primaryFit=p=>{proof={value:p.value,lower:p.lower};};let error:any;
 try{analyzeAutomaticSanger({references:t.amplicons.map((r:any)=>({...r,cut:r.cuts[0]})),primer:t.primer,control,sample,maxDeletion:20,maxInsertion:20,indelThreshold:2,groupRecovery:true},stage=>{stages.push(stage);if(stage!=='Bounding length proportions · target 1')return;const x=model.targets[0].point;captured={originalViolation:marginalViolation(model.program,x),phaseClassCounts:model.states[0].map((s:any[])=>s.length),references:model.refs.map((r:any,h:number)=>({name:r.name,wt:model.targets[0].groups.filter((g:any)=>g.reference===h&&g.size===0).reduce((a:number,g:any)=>a+g.mass.reduce((n:number,[i,w]:number[])=>n+w*x[i],0)/model.shares[h],0)}))};throw stop;},optimizer);}catch(e){error=e;}finally{marginalModelDiagnostics.observe=undefined;marginalModelDiagnostics.primaryFit=undefined;}
 const row={scramble,seconds:(performance.now()-started)/1000,stages,proof,captured,error:String(error)};if(root)writeFileSync(root+'/'+(scramble?'corrupted':'valid-wt')+'.json',JSON.stringify(row,null,2));return{row,error,stop};
}
it('accepts a valid 3-reference WT counterexample despite28<62 dimensional heuristic, retaining every length class',async()=>{
 const {row,error,stop}=await run(false);expect(error).toBe(stop);expect(row.captured.originalViolation).toBeLessThanOrEqual(1e-7);expect(row.captured.phaseClassCounts).toEqual([41,41,41]);for(const ref of row.captured.references)expect(ref.wt).toBeGreaterThan(.999);
},600000);
it('still rejects corrupted measured flanks through original model checks, not the old hypothesis-count guard',async()=>{
 const {row,error,stop}=await run(true);expect(error).toBeDefined();expect(error).not.toBe(stop);expect(String(error)).not.toContain('Insufficient invariant flank');expect(row.stages).toContain('Fitting invariant-flank length proportions');expect(row.captured).toBeUndefined();
},600000);
