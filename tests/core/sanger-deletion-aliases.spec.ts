import {it,expect} from 'vitest';
import {aliasedDeletionMass,fitJunctionMixture} from '../../src/app/workers/core/sanger-junction-mixture';
import {fitEmissionMixture} from '../../src/app/workers/core/sanger-emission';
it('finds alternative deletion lengths even when the inactive model uses another origin mixture',()=>{
 const sequence='CGGCAC'.repeat(20)+'AT'.repeat(200),cut=118,models=[{reference:0,size:-10,start:118,inserted:'',origins:[{start:118,weight:.6},{start:117,weight:.4}]}];
 expect(aliasedDeletionMass([{sequence,cut}],models,[.3],Array.from({length:250},(_,i)=>i+10),20,0)).toBeCloseTo(.3);
 // Observing the actual molecule ends breaks this alias.
 expect(aliasedDeletionMass([{sequence,cut}],models,[.3],Array.from({length:520},(_,i)=>i),20,0)).toBe(0);
});
it('keeps the observations and finite output when the deletion limit exceeds the entire reference',()=>{
 const sequence='ACGT'.repeat(30),positions=Array.from({length:100},(_,i)=>i),observed=positions.flatMap(p=>Array.from('ACGT',b=>Number(sequence[p]===b))),fit=fitJunctionMixture({references:[{sequence,cut:55}],positions,observed,maxDeletion:150,maxInsertion:1,radius:5,lag:0,theta:Array(8).fill(0),fit:fitEmissionMixture});
 expect(fit.weights.every(Number.isFinite)).toBe(true);expect(fit.weights.reduce((s,v)=>s+v,0)).toBeCloseTo(1);
 fit.models.forEach((m,i)=>{expect(fit.columns[i]).toHaveLength(observed.length);if(m.size< -sequence.length)expect(fit.columns[i].every(v=>v===0)).toBe(true);});
});
