import {it,expect} from 'vitest';
import {polishSignalDictionary} from '../../src/app/workers/core/sanger-dictionary-polish';
const col=(s:string)=>[...s].flatMap(b=>[...'ACGT'].map(c=>Number(c===b)));
it('jointly reassigns competing fluorescence columns while preserving measured family mass',()=>{
 const columns=['AAAAAA','AACAAA','AAAGAA','CCCCCC','CCCCGC'].map(col),truth=[.4,.1,.1,.2,.2],observed=columns[0].map((_,k)=>columns.reduce((s,c,i)=>s+c[k]*truth[i],0));
 const fit=polishSignalDictionary(columns,observed,['A','A','A','B','B'],[0,3],[.6,.4]);
 const w=columns.map((_,id)=>fit.weights[fit.ids.indexOf(id)]||0);expect(w.slice(0,3).reduce((a,b)=>a+b)).toBeCloseTo(.6,8);expect(w.slice(3).reduce((a,b)=>a+b)).toBeCloseTo(.4,8);for(let i=0;i<w.length;i++)expect(w[i]).toBeCloseTo(truth[i],4);
});
it('does not allocate a new reference family absent from the signal-derived initialization',()=>{const columns=['AAAA','ACAA','CCCC'].map(col),fit=polishSignalDictionary(columns,col('AAAA'),['A','A','B'],[0],[1]);expect(fit.ids).toEqual([0]);expect(fit.weights).toEqual([1]);});
