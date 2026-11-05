import {it,expect} from 'vitest';
import {compoundBaseProposals} from '../../src/app/workers/core/sanger-compound-refine';
import {materialize,PathEdit} from '../../src/app/workers/core/sanger-path-common';
it.each([-3,4])('proposes the downstream inserted bases at the actual shifted read cycle after T1 size %s',size=>{
 const reference='ACGT'.repeat(130),first:PathEdit=size<0?{start:150,end:150-size,inserted:''}:{start:150,end:150,inserted:'TTAA'},wrong=[first,{start:350,end:350,inserted:'GA'}],actual=[first,{start:350,end:350,inserted:'GC'}],emit=(seq:string)=>Array.from({length:480},(_,p)=>[...'ACGT'].map(b=>Number(seq[p]===b))).flat(),wt=emit(reference),column=emit(materialize(reference,wrong)!),truth=emit(materialize(reference,actual)!),weight=.2,pred=wt.map((v,k)=>.8*v+weight*column[k]),observed=wt.map((v,k)=>.8*v+weight*truth[k]),residual=observed.map((v,k)=>v-pred[k]),index=Int32Array.from({length:480},(_,i)=>i);
 const proposals=compoundBaseProposals(wrong,column,residual,weight,index);expect(proposals.some(events=>events[1].inserted==='GC'&&JSON.stringify(events[0])===JSON.stringify(first))).toBe(true);expect(wrong[1].inserted).toBe('GA');
});
