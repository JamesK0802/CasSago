import {it,expect} from 'vitest';
import loadCommon from '../probes/load-common-optimizer';
import {normTreeRuntime} from '../probes/conic-norm-tree';
import {VerifiedConicSolver} from '../../src/app/workers/core/sanger-conic-solver';
import {validMarginalSolution,type MarginalProgram,type LinearTerm} from '../../src/app/workers/core/sanger-marginal-lp';
it('maps a decomposed cone back to a certified original primal and dual at known analytic extrema',async()=>{
 const engine=await loadCommon(),runtime=normTreeRuntime(engine.conic,8),terms:LinearTerm[][]=Array.from({length:256},(_,j)=>[[j%8,1/Math.sqrt(32)],[8,-.3/Math.sqrt(32)]]),p:MarginalProgram={upper:Array(9).fill(1),cost:Array(9).fill(0),rows:[{terms:[[8,1]],lower:1,upper:1}],quadratics:[{name:'sphere',terms,upper:.01}]},objective:LinearTerm[]=Array.from({length:8},(_,i)=>[i,i+1]),center=.3*36,radius=.1*Math.sqrt(204),solver=new VerifiedConicSolver(p,runtime);
 try{for(const maximum of [false,true]){const answer=solver.solve(objective,maximum);expect(validMarginalSolution(p,answer.x)).toBe(true);expect(answer.value).toBeCloseTo(center+(maximum?radius:-radius),6);expect(Math.abs((maximum?-answer.value:answer.value)-solver.objectiveLowerBound()!)).toBeLessThan(1e-7);for(const plane of solver.supportPlanes()){const witness=[...Array(8).fill(.3),1];expect(plane.terms.reduce((s,[i,w])=>s+w*witness[i],0)).toBeLessThanOrEqual(plane.upper+1e-10);}}expect(runtime.transformed).toBeGreaterThan(0);}finally{solver.dispose();}
},120000);
