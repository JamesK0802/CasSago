import {it,expect} from 'vitest';import loadCommon from '../probes/load-common-optimizer';
import {verifiedZeroDeparturePoint} from '../../src/app/workers/core/sanger-initial-marginal-fit';
import {validMarginalSolution,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
const fixture=():MarginalProgram=>({upper:[1,1,1],cost:[0,0,1],rows:[{terms:[[0,1]],lower:1,upper:1},{terms:[[1,1],[2,1]],lower:1,upper:1},{terms:[[1,1],[0,-1]],lower:-1e30,upper:0}],quadratics:[{name:'geometry',terms:[[[0,1]]],upper:1}],zeroCostSubstitutions:[{column:1,terms:[[0,1]]}]});
it('certifies the original global zero-departure optimum after exact elimination',async()=>{const p=fixture(),snapshot=JSON.stringify(p),answer=verifiedZeroDeparturePoint(p,[{column:0,value:1}],await loadCommon());expect(answer).toBeDefined();expect(answer!.value).toBe(0);expect(validMarginalSolution(p,answer!.x)).toBe(true);expect(JSON.stringify(p)).toBe(snapshot);});
it('does not accept a zero-departure shortcut when a measured substitution is required',async()=>{const p=fixture();p.rows.push({terms:[[2,1]],lower:.05,upper:1});expect(verifiedZeroDeparturePoint(p,[{column:0,value:1}],await loadCommon())).toBeUndefined();});

it('preserves scaled duplicate rows, original cones and the complete lifted zero face',async()=>{
 const {zeroCostProjection}=await import('../../src/app/workers/core/sanger-zero-cost-projection');
 const p:MarginalProgram={upper:[1,2,1],cost:[0,0,1],rows:[{terms:[[0,1],[1,-.5],[2,-1]],lower:0,upper:0},{terms:[[0,1]],lower:.1,upper:.9},{terms:[[0,-2]],lower:-1.2,upper:-.4}],quadratics:[{name:'original cone',terms:[[[0,1]]],upper:.25}],zeroCostSubstitutions:[{column:1,terms:[[0,2]]}]},before=JSON.stringify(p),q=zeroCostProjection(p)!;
 for(let i=0;i<=20;i++){const x=[i/20],lifted=q.lift(x);expect(validMarginalSolution(q.program,x)).toBe(validMarginalSolution(p,lifted));expect(lifted[1]).toBe(2*x[0]);expect(lifted[2]).toBe(0);}expect(JSON.stringify(p)).toBe(before);
});
