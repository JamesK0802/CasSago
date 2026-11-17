import {it,expect} from 'vitest';
import {fittedMarginalRepresentative} from '../../src/app/workers/core/sanger-marginal-representative';
it('does not turn upstream phase mixtures into substitutions in an unchanged second target',()=>{
 const letters=Array(80).fill('A');letters[8]=letters[9]='C';
 const edges=[0,1,2].map((incoming,mass)=>({reference:0,target:1,incoming,outgoing:incoming,size:0,start:10,mass,free:new Map()}));
 const point=[.4,.35,.25,.4,.6,0,0,0,0,.75,.25,0,0,0,0];
 const group:any={reference:0,target:1,size:0,edges,mass:[[0,1],[1,1],[2,1]],bases:new Map([[10,[[[3,1]],[[4,1]],[],[]]],[11,[[[9,1]],[[10,1]],[],[]]]]),flow:new Map([[10,{unchanged:[{base:0,column:3,cap:[[0,1]]},{base:1,column:4,cap:[[1,1],[2,1]]}],changed:[5,6,7,8]}],[11,{unchanged:[{base:0,column:9,cap:[[0,1],[1,1]]},{base:1,column:10,cap:[[2,1]]}],changed:[11,12,13,14]}]])};
 const ref:any={sequence:letters.join(''),cuts:[0,10],input:{sequence:letters.join('')}};
 const answer=fittedMarginalRepresentative(group,point,ref,1,(start,end,word)=>({target:2,start,end,inserted:word.length,insertedSequence:word}));
 expect(answer?.events).toEqual([]);expect(answer?.witnessFraction).toBeCloseTo(.4,10);
});
it('retains a genuine modal substitution using the same per-origin support',()=>{
 const edge:any={reference:0,target:0,incoming:0,outgoing:0,size:0,start:10,mass:0,free:new Map()};
 const point=[1,.3,0,.7,0,0],group:any={reference:0,target:0,size:0,edges:[edge],mass:[[0,1]],bases:new Map([[10,[[[1,1]],[[3,1]],[],[]]]]),flow:new Map([[10,{unchanged:[{base:0,column:1,cap:[[0,1]]}],changed:[2,3,4,5]}]])},ref:any={sequence:'A'.repeat(80),cuts:[10],input:{sequence:'A'.repeat(80)}};
 const answer=fittedMarginalRepresentative(group,point,ref,1,(start,end,word)=>({target:1,start,end,inserted:word.length,insertedSequence:word}));
 expect(answer?.events).toEqual([{target:1,start:10,end:11,inserted:1,insertedSequence:'C'}]);expect(answer?.witnessFraction).toBeCloseTo(.7,10);
});
