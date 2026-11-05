import {it,expect} from 'vitest';
import {fitJunctionMixture,junctionBase,junctionStarts,junctionAnnotationComponents} from '../../src/app/workers/core/sanger-junction-mixture';
import {fitOrderedOnsets} from '../../src/app/workers/core/sanger-onset-fit';
import {simplexFit} from '../../src/app/workers/core/sanger-signal';
import {sangerAnalysisDemo,sangerPooledDemo} from '../../src/app/utils/sanger-analysis-demo';
import {parseAbif} from '../../src/app/utils/abif';
import {analyzeSanger} from '../../src/app/workers/core/sanger-analysis';
import {analyzePooledSanger} from '../../src/app/workers/core/sanger-pooled';
import {SANGER_ANALYSIS_DEFAULTS} from '../../src/app/models/sanger-analysis';

it('models left, spanning and right deletions using the actual removed interval',()=>{
 const sequence='ACGTGCTAGTCA',cut=6,size=-3;
 const starts=junctionStarts(cut,size,2,sequence.length);
 for(const start of [3,5,7]){
  expect(starts).toContain(start);
  const model={reference:0,size,start,inserted:''};
  expect(Array.from({length:sequence.length+size},(_,p)=>junctionBase(sequence,model,p)).join('')).toBe(sequence.slice(0,start)+sequence.slice(start-size));
 }
});
it('solves ordered onsets with the same least-squares objective as a full simplex fit',()=>{
 const columns=Array.from({length:12},(_,onset)=>Array.from({length:24},(_,p)=>[Number(p<onset),Number(p>=onset),0,0]).flat());
 const truth=columns.map((_,i)=>i===2?.2:i===7?.5:i===10?.3:0),observed=columns[0].map((_,k)=>columns.reduce((s,c,i)=>s+c[k]*truth[i],0));
 const fast=fitOrderedOnsets(columns,observed,undefined,columns.map(()=>1/columns.length))!,generic=simplexFit(columns,observed).weights;
 const error=(w:number[])=>observed.reduce((s,y,k)=>s+(y-columns.reduce((a,c,i)=>a+c[k]*w[i],0))**2,0);
 expect(fast.every(w=>w>=0)).toBe(true);expect(fast.reduce((s,w)=>s+w,0)).toBeCloseTo(1,10);
 expect(error(fast)).toBeLessThan(1e-10);expect(Math.abs(error(fast)-error(generic))).toBeLessThan(1e-9);
});
it('preserves both same-length junction alternatives and maps reverse coordinates',()=>{
 const models=[{reference:0,size:-4,start:5,inserted:'',origins:[{start:5,weight:.3},{start:12,weight:.7}]}];
 const forward=junctionAnnotationComponents(models,[.6],10,100),reverse=junctionAnnotationComponents(models,[.6],10,100,true);
 expect(forward.map(c=>c.events[0].start)).toEqual([15,22]);expect(forward.map(c=>c.fraction)).toEqual([.18,.42]);
 expect(reverse.map(c=>[c.events[0].start,c.events[0].end])).toEqual([[81,85],[74,78]]);
});
it('retains the same measured cycles when single or pooled indel bounds increase',()=>{
 const d=sangerAnalysisDemo(),single={referenceName:'range check',reference:d.reference,cut:d.cut,control:parseAbif(d.control),sample:parseAbif(d.sample)},p=sangerPooledDemo(),pooled={...p,ampliconConfirmed:true,control:parseAbif(p.control),sample:parseAbif(p.sample)};
 for(const analyze of [(maxIndel:number)=>analyzeSanger({...single,settings:{...SANGER_ANALYSIS_DEFAULTS,maxIndel,checkWiderRange:false}}),(maxIndel:number)=>analyzePooledSanger({...pooled,settings:{...SANGER_ANALYSIS_DEFAULTS,scope:'pooled',maxIndel,checkWiderRange:false}})]){
  const a=analyze(5),b=analyze(10),rows=a.diagnostics!.junctionFit!.positions;
  expect(b.diagnostics!.junctionFit!.positions).toEqual(rows);
  expect(rows.length).toBeGreaterThan(300);expect(a.mixture!.components.every(c=>Number.isFinite(c.fraction))).toBe(true);
 }
});


it('uses the junction radius to search origins without changing observations',()=>{
 let seed=40197;const sequence=Array.from({length:160},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return 'ACGT'[seed>>>30];}).join(''),cut=65,start=cut+8,edited=sequence.slice(0,start)+'GC'+sequence.slice(start),positions=Array.from({length:130},(_,i)=>i+10),observed=positions.flatMap(p=>Array.from('ACGT',b=>.6*Number(sequence[p]===b)+.4*Number(edited[p]===b)));
 const search=(radius:number)=>fitJunctionMixture({references:[{sequence,cut}],positions,observed,maxDeletion:2,maxInsertion:2,radius,lag:0,theta:Array(8).fill(0),fit:simplexFit});
 const narrow=search(1),wide=search(10);
 expect(wide.objective).toBeLessThan(narrow.objective*.1);
 expect(wide.weights[wide.models.findIndex(m=>m.size===2)]).toBeCloseTo(.4,2);
 expect(wide.columns[0]).toHaveLength(observed.length);expect(narrow.columns[0]).toHaveLength(observed.length);
});
it('rejects invalid junction radii in ordinary single and homeolog modes',()=>{
 const d=sangerAnalysisDemo(),p=sangerPooledDemo();
 expect(()=>analyzeSanger({referenceName:'bad radius',reference:d.reference,cut:d.cut,control:parseAbif(d.control),sample:parseAbif(d.sample),settings:{...SANGER_ANALYSIS_DEFAULTS,editStartRadius:-1}})).toThrow(/radius/);
 expect(()=>analyzePooledSanger({...p,ampliconConfirmed:true,control:parseAbif(p.control),sample:parseAbif(p.sample),settings:{...SANGER_ANALYSIS_DEFAULTS,scope:'pooled',editStartRadius:101}})).toThrow(/radius/);
});

it('calibrates leading-cycle fluorescence without erasing real weak +1 and -1 edits after the cut',async()=>{
 const {calibrateSignal,phaseVector}=await import('../../src/app/workers/core/sanger-signal');let seed=771219;const sequence=Array.from({length:440},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return 'ACGT'[seed>>>30];}).join(''),cut=180,deleted=sequence.slice(0,cut-1)+sequence.slice(cut),inserted=sequence.slice(0,cut)+'C'+sequence.slice(cut),lag=-.075;
 const sample=Array.from({length:400},(_,p)=>phaseVector(sequence,p,lag).map((v,b)=>.975*v+.02*phaseVector(deleted,p,lag)[b]+.005*phaseVector(inserted,p,lag)[b])),template=Array.from(sequence,base=>Array.from('ACGT',b=>Number(b===base))),calibration=calibrateSignal(sample,Array.from({length:135},(_,p)=>p+10),[template],'fine'),positions=Array.from({length:370},(_,p)=>p+10),fit=fitJunctionMixture({references:[{sequence,cut}],positions,observed:positions.flatMap(p=>sample[p]),maxDeletion:2,maxInsertion:2,radius:5,lag:calibration.lag,theta:Array(8).fill(0),fit:simplexFit});
 expect(calibration.lag).toBeCloseTo(lag,3);expect(fit.weights[fit.models.findIndex(m=>m.size===-1)]).toBeGreaterThan(.019);expect(fit.weights[fit.models.findIndex(m=>m.size===1)]).toBeGreaterThan(.004);expect(fit.columns[0]).toHaveLength(positions.length*4);
});

it('retains predictive evidence for a real0.5% insertion without reusing full-fit insertion bases',async()=>{
 const {junctionPredictionEvidence}=await import('../../src/app/workers/core/sanger-junction-evidence');let seed=449117;const sequence=Array.from({length:440},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return 'ACGT'[seed>>>30];}).join(''),cut=180,edited=sequence.slice(0,cut)+'C'+sequence.slice(cut),positions=Array.from({length:370},(_,p)=>p+10),input={references:[{sequence,cut}],positions,observed:positions.flatMap(p=>Array.from('ACGT',b=>.995*Number(sequence[p]===b)+.005*Number(edited[p]===b))),maxDeletion:2,maxInsertion:2,radius:5,lag:0,theta:Array(8).fill(0)},full=fitJunctionMixture({...input,fit:simplexFit}),before=JSON.stringify(full.weights),index=full.models.findIndex(m=>m.size===1),evidence=junctionPredictionEvidence(input,full,[index]).get(index)!;
 expect(evidence.passed).toBe(true);expect(evidence.foldGains.every(g=>g>0)).toBe(true);expect(JSON.stringify(full.weights)).toBe(before);expect(full.weights[index]).toBeGreaterThan(.004);
});
