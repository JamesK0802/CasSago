import loadHighs from'../../public/vendor/highs-1.15.4/highs.mjs';
import{it,expect}from'vitest';import{refineDiscretePhase,assignSequenceBases}from'../../src/app/workers/core/sanger-discrete-phase';import{materialize,PathEdit}from'../../src/app/workers/core/sanger-path-common';
it('recovers coupled insertion/substitution phase from measured signal with a joint discrete step',async()=>{
 let seed=271828;let sequence=Array.from({length:320},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return'ACGT'[seed>>>30];}).join('');const cut=180,region='TACGGACTAATCCCGGCC';sequence=sequence.slice(0,cut-3)+region+sequence.slice(cut-3+region.length);
 const prefix:PathEdit[][]=[[],[{start:80,end:80,inserted:'TGCAAATC'}]],initial:PathEdit[][]=[[],[{start:cut+3,end:cut+12,inserted:''}],[{start:cut+1,end:cut+2,inserted:'T'}],[{start:cut+1,end:cut+7,inserted:'TAGACC'}],[{start:cut,end:cut,inserted:'AGCAG'}],[{start:cut-2,end:cut-2,inserted:'AGTTACTGCA'}]],truth:PathEdit[][]=[[],initial[1],[{start:cut+1,end:cut+3,inserted:'TC'}],[{start:cut+3,end:cut+7,inserted:'GACC'}],[{start:cut,end:cut,inserted:'ATAAG'}],initial[5]],fractions=[.4044,.1474,.1368,.1115,.1338,.0661],coupling=[0,0,.006,0,-.006,0],positions=Array.from({length:250},(_,i)=>i+30),observed=positions.flatMap(p=>[...'ACGT'].map(base=>prefix.reduce((sum,first,i)=>sum+truth.reduce((s,second,j)=>s+Number(materialize(sequence,[...first,...second])![p]===base)*(fractions[j]*(i===0?.7:.3)+(i===0?1:-1)*coupling[j]),0),0)));
 const result=refineDiscretePhase([{sequence,cuts:[80,cut]}],[prefix],[initial],positions,observed,await loadHighs()),found=result.profiles[0].map(es=>materialize(sequence,es));
 for(const j of [2,3,4])expect(found).toContain(materialize(sequence,truth[j]));expect(result.rmse).toBeLessThan(result.history[0]);expect(result.changes).toBeGreaterThan(0);expect(JSON.stringify(result.profiles)).not.toContain('N');
},60000);

it.each([0,1])('refits coupled deletion origins at target %s without truth candidates or fixed pairings',async target=>{
 let seed=91473;const sequence=Array.from({length:360},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return'ACGT'[seed>>>30];}).join(''),cut=210;
 const otherCut=target===0?290:75,prefix:PathEdit[][]=[[],[{start:otherCut,end:otherCut,inserted:'ATGC'}]],initial:PathEdit[][]=[[],[{start:cut+2,end:cut+9,inserted:''}],[{start:cut-3,end:cut+3,inserted:''}]],truth:PathEdit[][]=[[],[{start:cut-2,end:cut+5,inserted:''}],[{start:cut,end:cut+6,inserted:''}]],weights=[[.39,.07,.04],[.31,.09,.10]],positions=Array.from({length:300},(_,i)=>i+20);
 const observed=positions.flatMap(p=>[...'ACGT'].map(base=>prefix.reduce((sum,first,i)=>sum+truth.reduce((s,second,j)=>s+weights[i][j]*Number(materialize(sequence,[...first,...second])![p]===base),0),0)));
 const result=refineDiscretePhase([{sequence,cuts:target===0?[cut,otherCut]:[otherCut,cut]}],[prefix],[initial],positions,observed,await loadHighs(),()=>{},target);
 expect(result.rmse).toBeLessThan(1e-7);
 for(const events of truth)expect(result.profiles[0].map(es=>materialize(sequence,es))).toContain(materialize(sequence,events));
},60000);

it('uses measured emission intervals when centroid conversion favors the opposite insertion/substitution phase',async()=>{
 let seed=99881;const sequence=Array.from({length:300},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return'ACGT'[seed>>>30];}).join(''),cut=160,truth:PathEdit[][]=[[],[{start:cut,end:cut+3,inserted:'AGC'}],[{start:cut,end:cut,inserted:'TTA'}]],wrong:PathEdit[][]=[[],[{start:cut,end:cut+3,inserted:'TTA'}],[{start:cut,end:cut,inserted:'AGC'}]],weights=[.7,.151,.149],positions=Array.from({length:16},(_,i)=>cut-3+i),signal=(events:PathEdit[][])=>positions.map(p=>[...'ACGT'].map(b=>events.reduce((s,es,j)=>s+weights[j]*Number(materialize(sequence,es)![p]===b),0))),measured=signal(truth),biased=signal(wrong),groups=wrong.map((events,j)=>({events,size:events.reduce((s,e)=>s+e.inserted.length-e.end+e.start,0),sequence:materialize(sequence,events)!,weights:new Map<number,number>([[0,weights[j]]])})),answer=assignSequenceBases(sequence,cut,groups,positions,biased,await loadHighs(),{lower:measured.map(row=>row.map(v=>v-1e-5)),upper:measured.map(row=>row.map(v=>v+1e-5))});
 expect(answer).not.toBeNull();const molecules=answer!.map(es=>materialize(sequence,es));for(const events of truth)expect(molecules).toContain(materialize(sequence,events));

},60000);
