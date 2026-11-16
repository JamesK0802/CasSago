import {it,expect} from 'vitest';
import {fittedMarginalRepresentative} from '../../src/app/workers/core/sanger-marginal-representative';
import type {MarginalInternalGroup,MarginalEdge} from '../../src/app/workers/core/sanger-marginal-engine';
import type {Prepared} from '../../src/app/workers/core/sanger-path-common';
import type {SangerJointEvent} from '../../src/app/models/sanger-analysis';
function fixture(){
 const point=[.25],edge:MarginalEdge={reference:0,target:0,start:150,incoming:0,outgoing:10,size:10,mass:0,free:new Map()},bases:MarginalInternalGroup['bases']=new Map();
 for(let q=150;q<160;q++){const ids=Array.from({length:4},(_,b)=>{const i=point.length;point.push(b===0?.25:0);return i;});edge.free.set(q,ids);bases.set(q,ids.map(i=>[[i,1]]));}
 const unchanged=point.length;point.push(0);const changed=Array.from({length:4},(_,b)=>{const i=point.length;point.push(b===1?.25:0);return i;});bases.set(200,changed.map(i=>[[i,1]]));
 const group:MarginalInternalGroup={id:'test',reference:0,target:0,size:10,mass:[[0,1]],bases,edges:[edge],fraction:.25,fractionRange:[.25,.25],flow:new Map([[200,{unchanged:[{base:0,column:unchanged,cap:[[0,1]]}],changed}]])};
 const ref={sequence:'A'.repeat(400),cuts:[150]}as Prepared;
 return{point,group,ref};
}
it('retains a supported genomic cut+40 substitution shifted past read cut+45 by an insertion',()=>{
 const {point,group,ref}=fixture(),r=fittedMarginalRepresentative(group,point,ref,1,(start,end,inserted)=>({target:1,start,end,inserted:inserted.length,insertedSequence:inserted}as SangerJointEvent));
 expect(r?.witnessFraction).toBe(.25);expect(r?.events).toEqual([{target:1,start:150,end:150,inserted:10,insertedSequence:'AAAAAAAAAA'},{target:1,start:190,end:191,inserted:1,insertedSequence:'C'}]);
});
