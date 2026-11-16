import {it,expect} from 'vitest';
import loadCommon from '../probes/load-common-optimizer';
import {summarizeMarginalPoint} from '../../src/app/workers/core/sanger-marginal-point-report';
import {marginalSessions} from '../../src/app/workers/core/sanger-marginal-state';
import {marginalGroupBounds} from '../../src/app/workers/core/sanger-marginal-group-bounds';
import type {SangerAnalysisResult} from '../../src/app/models/sanger-analysis';

it('bounds only a selected SNP-marker mixture and retains unresolved cross-position linkage',async()=>{
 const point=[1,.8,.2,.999,.001,0,0,0],program={upper:[1,1,1,1,1,0,0,0],cost:Array(8).fill(0),rows:[{terms:[[0,1]],lower:1,upper:1},{terms:[[1,1],[2,1],[0,-1]],lower:0,upper:0},{terms:[[3,1],[4,1],[0,-1]],lower:0,upper:0},{terms:[[2,1]],lower:.19,upper:.21},{terms:[[4,1]],lower:.0005,upper:.002}]},edge={reference:0,target:0,incoming:0,outgoing:0,size:0,start:1,mass:0,free:new Map()},changed=(g:number)=>[5,6,g,7],raw={id:'zero',reference:0,target:0,size:0,mass:[[0,1]],bases:new Map([[2,[[],[[1,1]],[[2,1]],[]]],[3,[[],[[3,1]],[[4,1]],[]]]]),edges:[edge],fraction:1,fractionRange:[0,1],flow:new Map([[2,{unchanged:[{base:1,column:1,cap:[[0,1]]}],changed:changed(2)}],[3,{unchanged:[{base:1,column:3,cap:[[0,1]]}],changed:changed(4)}]])},ref={name:'R',sequence:'AACCAA',offset:0,reverse:false,cuts:[1],input:{name:'R',sequence:'AACCAA',cut:1,cuts:[1]}},target={program,groups:[raw],point,positions:[2,3]},model={program,refs:[ref],shares:[1],drafts:[target]}as any;
 const groups=summarizeMarginalPoint(model,point),parent=groups[0],child=parent.displayComponents!.find(c=>c.events.some(e=>e.start===2))!,other=parent.displayComponents!.find(c=>c!==child)!,beforeOther=JSON.stringify(other),beforeParentBases=JSON.stringify(parent.bases),beforeProgram=JSON.stringify(program),result={marginalGroups:groups,warnings:[]}as unknown as SangerAnalysisResult;
 marginalSessions.set(result,{refs:[ref],shares:[1],targets:[target],optimizer:await loadCommon(),pages:new Map()}as any);
 expect(marginalGroupBounds(result,child.id,()=>{},0).complete).toBe(false);
 let answer=marginalGroupBounds(result,child.id);while(!answer.complete)answer=marginalGroupBounds(result,child.id);
 expect(answer.group).toBe(parent);expect(child.rangeStatus).toBe('bounded');expect(parent.status).toBe('estimated');expect(child.fractionRange[0]).toBeCloseTo(.19,5);expect(child.fractionRange[1]).toBeCloseTo(.21,5);expect(child.bases.find(p=>p.position===2)!.lower[2]).toBe(1);expect(child.bases.find(p=>p.position===2)!.upper[2]).toBe(1);
 const linked=child.bases.find(p=>p.position===3)!;expect(linked.lower[2]).toBe(0);expect(linked.upper[2]).toBeGreaterThanOrEqual(.002/.19);expect(linked.upper[2]).toBeLessThan(.011);expect(linked.lower[1]).toBeGreaterThan(.989);expect(linked.upper[1]).toBe(1);
 expect(JSON.stringify(other)).toBe(beforeOther);expect(JSON.stringify(parent.bases)).toBe(beforeParentBases);expect(JSON.stringify(program)).toBe(beforeProgram);expect(marginalGroupBounds(result,child.id).seconds).toBe(0);marginalSessions.delete(result);
});
