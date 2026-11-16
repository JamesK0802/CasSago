import {it,expect} from 'vitest';
import {summarizeMarginalPoint} from '../../src/app/workers/core/sanger-marginal-point-report';

function fixture(minor:number){
 const point=[1,.8,.2,1-minor,minor,0,0,0],edge={reference:0,target:0,incoming:0,outgoing:0,size:0,start:1,mass:0,free:new Map()},changed=(g:number)=>[5,6,g,7],group={id:'zero',reference:0,target:0,size:0,mass:[[0,1]],bases:new Map([[2,[[],[[1,1]],[[2,1]],[]]],[3,[[],[[3,1]],[[4,1]],[]]]]),edges:[edge],fraction:1,fractionRange:[0,1],flow:new Map([[2,{unchanged:[{base:1,column:1,cap:[[0,1]]}],changed:changed(2)}],[3,{unchanged:[{base:1,column:3,cap:[[0,1]]}],changed:changed(4)}]])},ref={name:'R',sequence:'AACCAA',offset:0,reverse:false,cuts:[1],input:{name:'R',sequence:'AACCAA',cut:1,cuts:[1]}};
 return{point,model:{program:{upper:[1,1,1,1,1,0,0,0]},refs:[ref],shares:[1],drafts:[{groups:[group],point:[]}] }as any};
}
it('keeps a reportable SNP marker split despite minor departures while conserving every base marginal',()=>{
 const{model,point}=fixture(.001),g=summarizeMarginalPoint(model,point,2)[0],parts=g.displayComponents!;expect(parts).toHaveLength(2);expect(parts.map(p=>p.fraction)).toEqual([.8,.2]);expect(parts.every(p=>p.mixture&&p.rangeStatus==='estimated')).toBe(true);
 for(const p of g.bases)for(let b=0;b<4;b++)expect(parts.reduce((n,c)=>n+c.fraction*c.bases.find(q=>q.position===p.position)!.estimate[b],0)).toBeCloseTo(g.fraction*p.estimate[b],12);
 expect(g.unmodifiedEstimate).toBeCloseTo(.799,12);expect(parts[0].events).toEqual([]);expect(parts[1].events).toEqual([{target:1,start:2,end:3,inserted:1,insertedSequence:'G'}]);
});
it('retains the full mixture when two reportable SNP sites have unresolved linkage',()=>{
 const{model,point}=fixture(.1),g=summarizeMarginalPoint(model,point,2)[0];expect(g.displayComponents).toBeUndefined();expect(g.bases[0].estimate).toEqual([0,.8,.2,0]);expect(g.bases[1].estimate).toEqual([0,.9,.1,0]);expect(g.unmodifiedEstimate).toBeCloseTo(.7,12);
});
it('keeps other modal substitutions in the representative of a low-abundance marker partition',()=>{
 const{model}=fixture(.001),point=[.03,.006,.024,.012,.018,0,0,0],g=summarizeMarginalPoint(model,point,2)[0],parts=g.displayComponents!;expect(parts).toHaveLength(2);expect(parts.every(p=>p.mixture)).toBe(true);expect(parts.every(p=>p.events.some(e=>e.start===3&&e.insertedSequence==='G'))).toBe(true);expect(g.unmodifiedEstimate).toBe(0);
});
