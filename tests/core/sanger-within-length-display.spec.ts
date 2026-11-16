import {expect,it} from 'vitest';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import {withinLengthDisplay} from '../../src/app/workers/core/sanger-marginal-display';
import type {MarginalInternalGroup,MarginalEdge} from '../../src/app/workers/core/sanger-marginal-engine';
import type {SangerMarginalGroup} from '../../src/app/models/sanger-analysis';
function fixture(sequence:string,size:number,edges:MarginalEdge[],point:number[],bases:Map<number,[number,number][]>[]){
 const forms=new Map(bases.map((b,q)=>[q,Array.from({length:4},(_,i)=>b.get(i)||[])])),g:MarginalInternalGroup={id:'g',reference:0,target:0,size,mass:edges.map(e=>[e.mass,1]),edges,bases:forms,fraction:.4,fractionRange:[.4,.4]},ref={sequence,input:{sequence,name:'R'},offset:0,reverse:false,cuts:[1]},display:SangerMarginalGroup={id:'g',reference:'R',target:1,size,fraction:.4,fractionRange:[.4,.4],scope:'read-cycle-by-length',status:'bounded',bases:[...forms].map(([position,forms])=>{const estimate=forms.map(f=>f.reduce((n,[i,w])=>n+w*point[i],0)/.4);return{position,estimate,lower:estimate.slice(),upper:estimate.slice()};}),representative:{events:[{target:1,start:edges[0].start,end:edges[0].start+Math.max(0,-size),inserted:Math.max(0,size),insertedSequence:size>0?'A'.repeat(size):''}],witnessFraction:.25,source:'fitted-marginal-proposal'}};
 const target={groups:[g],edges,point,positions:[0,1],tolerance:0,program:{upper:point.map(()=>1),cost:point.map(()=>0),rows:point.map((v,i)=>({terms:[[i,1]as[number,number]],lower:v,upper:v}))}};
 return{g,ref,display,target};
}
const edge=(start:number,mass:number,size:number,free=new Map<number,number[]>()):MarginalEdge=>({reference:0,target:0,incoming:0,outgoing:size,size,start,mass,free});
it('separates two same-length inserted bases using their fitted mass, with no GT and no changed signal',async()=>{
 const f=fixture('CCGT',1,[edge(1,0,1,new Map([[1,[1,2,3,4]]]))],[.4,.25,0,.15,0],[new Map([[1,[[0,1]]]]),new Map([[0,[[1,1]]],[2,[[3,1]]]])]),snapshot=JSON.stringify(f),optimizer=await loadHighs(),rows=withinLengthDisplay({refs:[f.ref],shares:[1],optimizer}as any,f.target,f.g,f.display)!;
 expect(rows.filter(r=>!r.mixture)).toHaveLength(2);expect(rows.filter(r=>r.fraction>1e-7).map(r=>r.fraction)).toEqual([.25,.15]);expect(rows.map(r=>r.events[0].insertedSequence)).toEqual(['A','G']);expect(rows.reduce((n,r)=>n+r.fraction,0)).toBeCloseTo(.4);expect(JSON.stringify(f)).toBe(snapshot);
});
it('separates same-size deletions at distinct supported positions',async()=>{
 const f=fixture('ACGTAC',-1,[edge(1,0,-1),edge(3,1,-1)],[.25,.15],[new Map([[0,[[0,1],[1,1]]]]),new Map([[2,[[0,1]]],[1,[[1,1]]]])]),rows=withinLengthDisplay({refs:[f.ref],shares:[1],optimizer:await loadHighs()}as any,f.target,f.g,f.display)!;expect(rows).toHaveLength(2);expect(rows.map(r=>r.events[0].start)).toEqual([1,3]);expect(rows.map(r=>r.fraction)).toEqual([.25,.15]);
});
it('merges identical deletion DNA rather than emitting placement aliases',async()=>{
 const f=fixture('AAAAAA',-1,[edge(1,0,-1),edge(2,1,-1)],[.25,.15],[new Map([[0,[[0,1],[1,1]]]])]);expect(withinLengthDisplay({refs:[f.ref],shares:[1],optimizer:await loadHighs()}as any,f.target,f.g,f.display)).toBeUndefined();
});
it('does not manufacture linkage between two independently mixed insertion sites',async()=>{
 const f=fixture('CCGT',2,[edge(1,0,2,new Map([[1,[1,2,3,4]],[2,[5,6,7,8]]]))],[.4,.25,0,.15,0,.25,0,.15,0],[new Map([[1,[[0,1]]]]),new Map([[0,[[1,1]]],[2,[[3,1]]]]),new Map([[0,[[5,1]]],[2,[[7,1]]]])]);expect(withinLengthDisplay({refs:[f.ref],shares:[1],optimizer:await loadHighs()}as any,f.target,f.g,f.display)).toBeUndefined();
});
it('retains an arbitrary non-identifiable allocation as a mixture',async()=>{
 const f=fixture('ACGTAC',-1,[edge(1,0,-1),edge(3,1,-1)],[.25,.15],[new Map([[0,[[0,1],[1,1]]]])]);f.target.program.rows=[{terms:[[0,1],[1,1]],lower:.4,upper:.4}];expect(withinLengthDisplay({refs:[f.ref],shares:[1],optimizer:await loadHighs()}as any,f.target,f.g,f.display)).toBeUndefined();
});
it('retains bases absent from the fitted point in a compact residual when they are still feasible',async()=>{
 const f=fixture('CCGT',1,[edge(1,0,1,new Map([[1,[1,2,3,4]]]))],[.4,.25,0,.15,0],[new Map([[1,[[0,1]]]]),new Map([[0,[[1,1]]],[1,[[2,1]]],[2,[[3,1]]],[3,[[4,1]]]])]);
 f.target.program.rows=[{terms:[[0,1]],lower:.4,upper:.4},{terms:[[1,1]],lower:.23,upper:.25},{terms:[[3,1]],lower:.13,upper:.15},{terms:[[1,1],[2,1],[3,1],[4,1]],lower:.4,upper:.4}];const rows=withinLengthDisplay({refs:[f.ref],shares:[1],optimizer:await loadHighs()}as any,f.target,f.g,f.display)!,rest=rows.find(r=>r.mixture)!;expect(rest.fraction).toBe(0);expect(rest.fractionRange[1]).toBeGreaterThan(.0399);expect(rest.selection.marker?.bases).toEqual([1,3]);expect(rows.reduce((n,r)=>n+r.fraction,0)).toBeCloseTo(.4);
});
it('proves ambiguous partitions without constructing a full parent optimization problem',async()=>{
 const f=fixture('ACGTAC',-1,[edge(1,0,-1),edge(3,1,-1)],[.25,.15],[new Map([[0,[[0,1],[1,1]]]])]);f.target.point.push(...Array(198).fill(0));f.target.program.upper.push(...Array(198).fill(1));f.target.program.cost.push(...Array(198).fill(0));f.target.program.rows=[{terms:[[0,1],[1,1]],lower:.4,upper:.4}];const optimizer=await loadHighs(),create=optimizer.createModel.bind(optimizer);let largest=0;optimizer.createModel=()=>{const native=create();return new Proxy(native,{get(t,k){if(k==='passModel')return(model:any)=>{largest=Math.max(largest,model.numCols);expect(model.numCols).toBeLessThan(10);return t.passModel(model);};const v=t[k];return typeof v==='function'?v.bind(t):v;}});};const before=JSON.stringify(f);expect(withinLengthDisplay({refs:[f.ref],shares:[1],optimizer}as any,f.target,f.g,f.display)).toBeUndefined();expect(largest).toBeLessThan(10);expect(JSON.stringify(f)).toBe(before);
});
