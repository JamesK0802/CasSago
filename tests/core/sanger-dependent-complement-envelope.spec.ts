import {it,expect} from 'vitest';
import loadCommon from '../probes/load-common-optimizer';
import {marginalGroupBounds} from '../../src/app/workers/core/sanger-marginal-group-bounds';
import {marginalSessions} from '../../src/app/workers/core/sanger-marginal-state';
import type {SangerAnalysisResult,SangerMarginalGroup} from '../../src/app/models/sanger-analysis';
it('uses the JOINT complementary-base budget instead of summing unattainable separate maxima',async()=>{
 // The measured companion change equals A+C+G, so the global0.1 event budget
 // permits at most0.05 non-T mass in this0.5 family. True T lower bound:90%.
 const program={
  upper:[.5,.5,.5,.5,.5,.5,0],cost:[0,1,1,1,0,1,0],rows:[{terms:[[0,1]],lower:.5,upper:.5},{terms:[[1,1],[2,1],[3,1],[4,1],[0,-1]],lower:0,upper:0},{terms:[[5,1],[1,-1],[2,-1],[3,-1]],lower:0,upper:0},{terms:[[1,1],[2,1],[3,1],[5,1]],lower:0,upper:.1}],dualHint:{budgetRow:3,shadow:[]}
 },point=[.5,.01,.01,.005,.475,.025,0],g={id:'parent',reference:0,target:0,size:-1,mass:[[0,1]],bases:new Map([[2,[[[1,1]],[[2,1]],[[3,1]],[[4,1]]]]]),edges:[{reference:0,target:0,incoming:0,outgoing:-1,size:-1,start:1,mass:0,free:new Map()}],fraction:.5,fractionRange:[0,1],flow:new Map([[2,{unchanged:[{base:3,column:4,cap:[[0,1]]}],changed:[1,2,3,5,6]}]])
 },ref={name:'R',sequence:'TTTTTT',offset:0,reverse:false,cuts:[2],input:{name:'R',sequence:'TTTTTT',cut:2,cuts:[2]}
 },target={program,groups:[g],point,positions:[2]},group:SangerMarginalGroup={id:'parent',reference:'R',target:1,size:-1,fraction:.5,fractionRange:[0,1],bases:[{position:2,estimate:[.02,.02,.01,.95],lower:[0,0,0,0],upper:[1,1,1,1]}],scope:'read-cycle-by-length',status:'estimated'},result={marginalGroups:[group],warnings:[]}as unknown as SangerAnalysisResult,before=JSON.stringify(program),estimate=group.bases[0].estimate.slice();

 marginalSessions.set(result,{refs:[ref],shares:[1],targets:[target],optimizer:await loadCommon(),pages:new Map()}as any);
 const detail=marginalGroupBounds(result,'parent',()=>{},60);
 expect(detail.complete).toBe(true);
 expect(group.bases[0].lower[3]).toBeGreaterThan(.89999);
 expect(group.bases[0].lower[3]).toBeLessThan(.90001);
 expect(group.bases[0].estimate).toEqual(estimate);
 expect(target.point).toEqual(point);
 expect(JSON.stringify(program)).toBe(before);
});

it('retains a fixed conditional fraction when numerator and group mass vary together',async()=>{
 const program={upper:[.9,.45,0,0,.45],cost:[0,0,0,0,0],rows:[
  {terms:[[0,1]],lower:.1,upper:.9},
  {terms:[[1,1],[0,-.5]],lower:0,upper:0},
  {terms:[[4,1],[0,-.5]],lower:0,upper:0},
  {terms:[[1,1],[4,1],[0,-1]],lower:0,upper:0},
 ]},point=[.5,.25,0,0,.25],g={id:'coupled',reference:0,target:0,size:-1,mass:[[0,1]],bases:new Map([[2,[[[1,1]],[],[],[[4,1]]]]]),edges:[{reference:0,target:0,incoming:0,outgoing:-1,size:-1,start:1,mass:0,free:new Map()}],fraction:.5,fractionRange:[0,1],flow:new Map([[2,{unchanged:[{base:3,column:4,cap:[[0,1]]}],changed:[1,2,3]}]])},ref={name:'R',sequence:'TTTTTT',offset:0,reverse:false,cuts:[2],input:{name:'R',sequence:'TTTTTT',cut:2,cuts:[2]}},target={program,groups:[g],point,positions:[2]},group:SangerMarginalGroup={id:'coupled',reference:'R',target:1,size:-1,fraction:.5,fractionRange:[0,1],bases:[{position:2,estimate:[.5,0,0,.5],lower:[0,0,0,0],upper:[1,1,1,1]}],scope:'read-cycle-by-length',status:'estimated'},result={marginalGroups:[group],warnings:[]}as unknown as SangerAnalysisResult,before=JSON.stringify(program);
 marginalSessions.set(result,{refs:[ref],shares:[1],targets:[target],optimizer:await loadCommon(),pages:new Map()}as any);
 expect(marginalGroupBounds(result,'coupled',()=>{},60).complete).toBe(true);
 expect(group.bases[0].lower[0]).toBeGreaterThan(.49999);
 expect(group.bases[0].upper[0]).toBeLessThan(.50001);
 expect(group.bases[0].estimate).toEqual([.5,0,0,.5]);
 expect(JSON.stringify(program)).toBe(before);
});

it('retains the mandatory coupled signal at another coordinate in a selected range query',async()=>{
 // The measured companion change equals A+C+G, so the global0.1 event budget
 // permits at most0.05 non-T mass in this0.5 family. True T lower bound:90%.
 const program={
  upper:[.5,.5,.5,.5,.5,.5,0],cost:[0,1,1,1,0,1,0],rows:[{terms:[[0,1]],lower:.5,upper:.5},{terms:[[1,1],[2,1],[3,1],[4,1],[0,-1]],lower:0,upper:0},{terms:[[5,1],[1,-1],[2,-1],[3,-1]],lower:0,upper:0},{terms:[[1,1],[2,1],[3,1],[5,1]],lower:0,upper:.1}],dualHint:{budgetRow:3,shadow:[]}
 },point=[.5,.01,.01,.005,.475,.025,0],g={id:'parent',reference:0,target:0,size:-1,mass:[[0,1]],bases:new Map([[2,[[[1,1]],[[2,1]],[[3,1]],[[4,1]]]]]),edges:[{reference:0,target:0,incoming:0,outgoing:-1,size:-1,start:1,mass:0,free:new Map()}],fraction:.5,fractionRange:[0,1],flow:new Map([[2,{unchanged:[{base:3,column:4,cap:[[0,1]]}],changed:[1,2,3,6]}]])
 },ref={name:'R',sequence:'TTTTTT',offset:0,reverse:false,cuts:[2],input:{name:'R',sequence:'TTTTTT',cut:2,cuts:[2]}
 },target={program,groups:[g],point,positions:[2]},group:SangerMarginalGroup={id:'parent',reference:'R',target:1,size:-1,fraction:.5,fractionRange:[0,1],bases:[{position:2,estimate:[.02,.02,.01,.95],lower:[0,0,0,0],upper:[1,1,1,1]}],scope:'read-cycle-by-length',status:'estimated'},result={marginalGroups:[group],warnings:[]}as unknown as SangerAnalysisResult,before=JSON.stringify(program),estimate=group.bases[0].estimate.slice();

 marginalSessions.set(result,{refs:[ref],shares:[1],targets:[target],optimizer:await loadCommon(),pages:new Map()}as any);
 const detail=marginalGroupBounds(result,'parent',()=>{},60);
 expect(detail.complete).toBe(true);
 expect(group.bases[0].lower[3]).toBeGreaterThan(.89999);
 expect(group.bases[0].lower[3]).toBeLessThan(.90001);
 expect(group.bases[0].estimate).toEqual(estimate);
 expect(target.point).toEqual(point);
 expect(JSON.stringify(program)).toBe(before);
});
