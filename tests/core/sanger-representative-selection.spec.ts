import {it,expect} from 'vitest';
import loadCommon from '../probes/load-common-optimizer';
import {selectMarginalRepresentative,representativeBurden} from '../../src/app/workers/core/sanger-representative-selection';
import {validMarginalSolution,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
const context=()=>({substitutionOpportunities:1,noise:.02,phaseSeed:[{column:0,value:.4},{column:1,value:.6}],states:[[[{size:0,column:0},{size:1,column:1}]]],shares:[1],edges:[]});
it('preserves measured phase proportions while regularizing the representative, with every parent alternative retained',async()=>{const p:MarginalProgram={upper:[1,1,1],cost:[0,0,1],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1}]},before=JSON.stringify(p),a=selectMarginalRepresentative(p,[0,1,0],[[1,1]],await loadCommon(),context());expect(a.point[1]).toBeCloseTo(.59,6);expect(a.summary.minimum).toBeCloseTo(.595,6);expect(validMarginalSolution(p,a.point)).toBe(true);expect(validMarginalSolution(p,[0,1,0])).toBe(true);expect(JSON.stringify(p)).toBe(before);});
it.each([5,6,12,100])('does not penalize splitting the same edited mass into %s classes',n=>{const p:MarginalProgram={upper:Array(n).fill(1),cost:Array(n).fill(0),rows:[]},edits=Array.from({length:n},(_,i)=>[i,1]as[number,number]),terms=representativeBurden(p,edits);expect(terms.reduce((s,[,w])=>s+w*.6/n,0)).toBeCloseTo(.6,14);});
it('retains required Sub support without changing the original uncertainty model',async()=>{const p:MarginalProgram={upper:[1,1,1],cost:[0,0,1],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1},{terms:[[2,1]],lower:.05,upper:.08}]},before=JSON.stringify(p),a=selectMarginalRepresentative(p,[.4,.6,.07],[[1,1]],await loadCommon(),context());expect(a.point[2]).toBeCloseTo(.05,7);expect(a.point[1]).toBeCloseTo(.59,6);expect(JSON.stringify(p)).toBe(before);});
it('certifies an exact zero-burden/zero-distance fit without numerical optimization',()=>{const p:MarginalProgram={upper:[1],cost:[0],rows:[{terms:[[0,1]],lower:1,upper:1}]},a=selectMarginalRepresentative(p,[1],[],{createModel:()=>{throw Error('No optimization is necessary.');}}as any,{substitutionOpportunities:1,noise:.02,phaseSeed:[{column:0,value:1}],states:[[[{size:0,column:0}]]],shares:[1],edges:[]});expect(a.summary.minimum).toBe(0);expect(a.summary.certificate).toBe('nonnegative-box-bound');expect(a.point).toHaveLength(1);});

it.each([1,5,100])('keeps phase inference invariant when the same nuisance signal occupies %s bases',async n=>{const p:MarginalProgram={upper:Array(n+2).fill(1),cost:[0,0,...Array(n).fill(1)],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1},...Array.from({length:n},(_,j)=>({terms:[[j+2,1],[0,-.2]]as[number,number][],lower:0,upper:0}))]},point=[.4,.6,...Array(n).fill(.08)],a=selectMarginalRepresentative(p,point,[[1,1]],await loadCommon(),{...context(),substitutionOpportunities:n});expect(a.point[1]).toBeCloseTo(.592,6);});

it.each([1,.5,.1])('does not increase variant shrinkage when an independent reference is diluted to share %s',async share=>{
 const p:MarginalProgram={upper:[share,share,1-share],cost:[0,0,0],rows:[{terms:[[0,1],[1,1]],lower:share,upper:share},{terms:[[2,1]],lower:1-share,upper:1-share}]};
 const a=selectMarginalRepresentative(p,[.4*share,.6*share,1-share],[[1,1]],await loadCommon(),{substitutionOpportunities:1,noise:.02,phaseSeed:[{column:0,value:.4*share},{column:1,value:.6*share},{column:2,value:1-share}],states:[[[{size:0,column:0},{size:1,column:1}],[{size:0,column:2}]]],shares:[share,1-share],edges:[]});
 expect(a.point[1]/share).toBeCloseTo(.59,6);expect(validMarginalSolution(p,a.point)).toBe(true);
});

it('uses measured fluorescence to resolve an observable within-length allocation without constraining alternatives',async()=>{
 const p:MarginalProgram={upper:[1,1,1],cost:[0,0,0],rows:[{terms:[[0,1],[1,1],[2,-1]],lower:0,upper:0},{terms:[[2,1]],lower:1,upper:1}]},before=JSON.stringify(p);
 const context={substitutionOpportunities:1,noise:.01,phaseSeed:[{column:2,value:1}],states:[[[{size:1,column:2}]]],shares:[1],edges:[],measurements:[{forms:[[[0,1]],[[1,1]],[],[]]as[number,number][][],total:[[0,1],[1,1]]as[number,number][],observed:[.8,.2,0,0]}]};
 const a=selectMarginalRepresentative(p,[.5,.5,1],[[2,1]],await loadCommon(),context);
 expect(a.point[0]).toBeCloseTo(.8,6);expect(a.point[1]).toBeCloseTo(.2,6);expect(a.summary.afterPhotonLoss).toBeLessThan(1e-7);expect(a.summary.beforePhotonLoss).toBeGreaterThan(1);expect(validMarginalSolution(p,a.point)).toBe(true);expect(validMarginalSolution(p,[.5,.5,1])).toBe(true);expect(JSON.stringify(p)).toBe(before);
});

it('does not claim two identical fluorescence columns identify their molecular allocation',async()=>{
 const p:MarginalProgram={upper:[1,1,1],cost:[0,0,0],rows:[{terms:[[0,1],[1,1],[2,-1]],lower:0,upper:0},{terms:[[2,1]],lower:1,upper:1}]};
 const context={substitutionOpportunities:1,noise:.01,phaseSeed:[{column:2,value:1}],states:[[[{size:1,column:2}]]],shares:[1],edges:[],measurements:[{forms:[[[0,1],[1,1]],[],[],[]]as[number,number][][],total:[[0,1],[1,1]]as[number,number][],observed:[1,0,0,0]}]};
 const a=selectMarginalRepresentative(p,[.5,.5,1],[[2,1]],await loadCommon(),context);
 expect(a.summary.afterPhotonLoss).toBeLessThan(1e-7);expect(validMarginalSolution(p,[0,1,1])).toBe(true);expect(validMarginalSolution(p,[1,0,1])).toBe(true);expect(a.point[0]+a.point[1]).toBeCloseTo(1,7);
});

it('stabilizes a sub-noise allocation while preserving every original alternative',async()=>{
 const p:MarginalProgram={upper:[1,1,1],cost:[0,0,0],rows:[{terms:[[0,1],[1,1],[2,-1]],lower:0,upper:0},{terms:[[2,1]],lower:1,upper:1}]},eps=1e-6,before=JSON.stringify(p);
 const context={substitutionOpportunities:1,noise:.01,phaseSeed:[{column:2,value:1}],states:[[[{size:1,column:2}]]],shares:[1],edges:[],measurements:[{forms:[[[0,1],[1,1-eps]],[[1,eps]],[],[]]as[number,number][][],total:[[0,1],[1,1]]as[number,number][],observed:[1-.8*eps,.8*eps,0,0]}],baseFlowAnchor:[{forms:[[[0,1]],[[1,1]],[],[]]as[number,number][][],mass:[[2,1]]as[number,number][],share:1}]};
 const a=selectMarginalRepresentative(p,[.5,.5,1],[[2,1]],await loadCommon(),context);
 expect(a.point[0]).toBeCloseTo(.5,5);expect(a.point[1]).toBeCloseTo(.5,5);expect(validMarginalSolution(p,[0,1,1])).toBe(true);expect(validMarginalSolution(p,[1,0,1])).toBe(true);expect(JSON.stringify(p)).toBe(before);
});

it('still recovers a resolved signal in the noise-scale proximal refinement',async()=>{
 const p:MarginalProgram={upper:[1,1,1],cost:[0,0,0],rows:[{terms:[[0,1],[1,1],[2,-1]],lower:0,upper:0},{terms:[[2,1]],lower:1,upper:1}]};
 const context={substitutionOpportunities:1,noise:.01,phaseSeed:[{column:2,value:1}],states:[[[{size:1,column:2}]]],shares:[1],edges:[],measurements:[{forms:[[[0,1]],[[1,1]],[],[]]as[number,number][][],total:[[0,1],[1,1]]as[number,number][],observed:[.8,.2,0,0]}],baseFlowAnchor:[{forms:[[[0,1]],[[1,1]],[],[]]as[number,number][][],mass:[[2,1]]as[number,number][],share:1}]};
 const a=selectMarginalRepresentative(p,[.5,.5,1],[[2,1]],await loadCommon(),context);
 expect(a.point[0]).toBeCloseTo(.78846153846,6);expect(a.point[1]).toBeCloseTo(.21153846154,6);expect(a.summary.afterPhotonLoss).toBeLessThan(a.summary.beforePhotonLoss);expect(validMarginalSolution(p,a.point)).toBe(true);
});
