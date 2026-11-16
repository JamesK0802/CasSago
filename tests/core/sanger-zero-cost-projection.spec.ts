import {it,expect} from 'vitest';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import {zeroCostProjection} from '../../src/app/workers/core/sanger-zero-cost-projection';
import {MarginalLinearSolver,validMarginalSolution,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
it('preserves the complete zero-cost face while eliminating implied base-flow variables',async()=>{
 const p:MarginalProgram={upper:[1,1,1,1,1],cost:[0,0,0,1,0],rows:[{terms:[[0,1],[1,1]],lower:1,upper:1},{terms:[[2,1],[0,-1]],lower:-1e30,upper:0},{terms:[[4,1],[1,-1]],lower:-1e30,upper:0},{terms:[[2,1],[4,1],[3,1]],lower:1,upper:1}],zeroCostSubstitutions:[{column:2,terms:[[0,1]]},{column:4,terms:[[1,1]]}]},projected=zeroCostProjection(p)!,original={...p,upper:p.upper.map((u,i)=>p.cost[i]>0?0:u)},optimizer=await loadHighs(),full=new MarginalLinearSolver(original,optimizer),small=new MarginalLinearSolver(projected.program,optimizer);try{expect(projected.program.upper.length).toBe(2);for(let j=0;j<5;j++){const terms:[[number,number],[number,number]]=[[0,j-2],[1,2-j]];for(const maximize of [false,true]){const a=full.solve(terms,maximize),b=small.solve(terms,maximize);expect(b.value).toBeCloseTo(a.value,9);expect(validMarginalSolution(p,projected.lift(b.x))).toBe(true);}}}finally{full.dispose();small.dispose();}
});
