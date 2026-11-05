import {it,expect} from 'vitest';
import {readFileSync} from 'node:fs';
import {parseAbif} from '../../src/app/utils/abif';
import {analyzeSanger} from '../../src/app/workers/core/sanger-analysis';
import {matchedSingleSignals} from '../../src/app/workers/core/sanger-single-signals';
const base=new URL('../fixtures/sanger/bulk-overlap/',import.meta.url);
const c=JSON.parse(readFileSync(new URL('case.json',base),'utf8'));
function load(name:string){const b=readFileSync(new URL(name+'.ab1',base));return parseAbif(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),name+'.ab1');}
it('decomposes WT plus three concurrent indels despite unequal overlapping dye shapes',()=>{
 const r=analyzeSanger({referenceName:'Overlap',reference:c.primer+c.reference,cut:c.primer.length+c.cut,control:load('control'),sample:load('sample'),settings:c.settings});
 const pred=new Map(r.mixture!.components.map(x=>[x.size,x.fraction]));
 expect(r.mixture!.passesChecks).toBe(true);
 for(const a of c.truth.net)expect(Math.abs((pred.get(a.size)||0)-a.fraction)).toBeLessThan(.05);
 expect(r.diagnostics!.background!.method).toBe('local-peak-area');
});
it('selects the representation from intact upstream data, independently of the edited fit',()=>{
 const control=load('control'),sample=load('sample');
 const cut=c.cut-10,positions=Array.from({length:100},(_,i)=>cut-130+i);
 const a=matchedSingleSignals(control,sample,positions,2);
 const changed={...sample,channels:Object.fromEntries(Object.entries(sample.channels).map(([b,v])=>[b,new Float32Array(v)]))} as typeof sample;
 // Replace the downstream raw signal; upstream windows are identical.
 for(const values of Object.values(changed.channels))values.fill(0,sample.peaks[cut+30]);
 const b=matchedSingleSignals(control,changed,positions,2);
 expect(b.method).toBe(a.method);expect(b.calibration).toEqual(a.calibration);
 expect(b.upstreamError).toBe(a.upstreamError);
});
