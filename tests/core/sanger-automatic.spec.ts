import {it,expect,vi,beforeEach} from 'vitest';
import {simplexFit} from '../../src/app/workers/core/sanger-signal';
import {analyzeAutomaticSanger,automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
const route=vi.hoisted(()=>({run:vi.fn(()=>({method:'pooled',marginalGroups:[]}))}));
vi.mock('../../src/app/workers/core/sanger-unified-pipeline',()=>({analyzeUnifiedSequenceGroups:route.run}));
const input=(references=1,targets=1)=>({references:Array.from({length:references},(_,h)=>({name:'reference-'+h,sequence:'ACGT'.repeat(200),cut:170,cuts:Array.from({length:targets},(_,t)=>170+t*250+h*50)})),primer:'ACGTCGTAGCTAGGATCCGA',sample:{}as any,control:{}as any,maxDeletion:10,maxInsertion:10});
beforeEach(()=>route.run.mockClear());
it.each([[1,1],[1,2],[3,1],[3,2]])('routes %s references and %s targets through the same model regardless of saved legacy flags', (refs,targets)=>{
 const optimizer={}as any;
 for(const flag of [true,false,undefined]){
  const d={...input(refs,targets),groupRecovery:flag};const r=analyzeAutomaticSanger(d,()=>{},optimizer);expect(r.marginalGroups).toEqual([]);expect(route.run).toHaveBeenCalledTimes([true,false,undefined].indexOf(flag)+1);const [actual,,used]=route.run.mock.calls.at(-1)!as any;expect(actual.amplicons).toBe(d.references);expect(actual.primer).toBe(d.primer);expect(actual.settings.maxDeletion).toBe(10);expect(actual.settings.maxInsertion).toBe(10);expect(used).toBe(optimizer);
 }
});
it('keeps asymmetric indel size limits without selecting a different engine',()=>{
 const d={...input(),maxDeletion:225,maxInsertion:6},s=automaticSangerSettings(d);analyzeAutomaticSanger(d,()=>{},{}as any);expect(s.maxDeletion).toBe(225);expect(s.maxInsertion).toBe(6);expect(s.analysisPolicy).toBe('automatic');expect((route.run.mock.calls[0]as any)[0].settings.maxIndel).toBe(225);
});
it('requires the actual sequencing primer and control rather than silently selecting a legacy model',()=>{
 expect(()=>analyzeAutomaticSanger({...input(),primer:undefined})).toThrow(/sequencing primer/);expect(()=>analyzeAutomaticSanger({...input(),control:undefined as any})).toThrow(/control/);expect(route.run).not.toHaveBeenCalled();
});
it('rejects unsupported target layouts before inference',()=>{
 expect(()=>analyzeAutomaticSanger(input(1,3))).toThrow(/one or two/);const d=input(2,2);d.references[1].cuts=[220];expect(()=>analyzeAutomaticSanger(d)).toThrow(/one or two/);expect(route.run).not.toHaveBeenCalled();
});
it('validates optional read-pair completeness before any secondary adapter',()=>{
 const d=input(),pair={primer:d.primer,sample:d.sample,control:d.control};expect(()=>analyzeAutomaticSanger({...d,additionalReads:[pair,pair]})).toThrow(/one opposite/);expect(()=>analyzeAutomaticSanger({...d,additionalReads:[{...pair,control:undefined as any}]})).toThrow(/Each read direction/);
});

it('solves a large redundant candidate family with the same objective as its convex-hull basis',()=>{let seed=89107;const random=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/4294967296;};const base=Array.from({length:8},()=>Array.from({length:120},()=>random())),columns=[...base,...Array.from({length:140},()=>{const w=base.map(()=>random()),sum=w.reduce((a,b)=>a+b,0);return base[0].map((_,k)=>base.reduce((a,c,i)=>a+c[k]*w[i]/sum,0));})],y=base[0].map((_,k)=>.3*base[1][k]+.7*base[5][k]+.03*(random()-.5)),a=simplexFit(base,y),b=simplexFit(columns,y);expect(Math.abs(a.rmse-b.rmse)).toBeLessThan(1e-7);expect(b.optimalityGap!).toBeLessThanOrEqual(1e-10);expect(b.weights.every(w=>w>=0)).toBe(true);expect(b.weights.reduce((a,b)=>a+b,0)).toBeCloseTo(1,8);});
