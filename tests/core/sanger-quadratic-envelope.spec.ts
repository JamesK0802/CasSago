import {it,expect} from 'vitest';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import {MarginalLinearSolver,validMarginalSolution,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
it('enforces the original convex residual envelope when optimizing linear summaries',async()=>{
 const p:MarginalProgram={upper:[1,1],cost:[0,0],rows:[],quadratics:[{name:'unit residual ball',terms:[[[0,1]],[[1,1]]],upper:1}]},solver=new MarginalLinearSolver(p,await loadHighs());try{const answer=solver.solve([[0,1],[1,1]],true);expect(answer.value).toBeCloseTo(Math.sqrt(2),6);expect(validMarginalSolution(p,answer.x)).toBe(true);expect(p.rows.length).toBeGreaterThan(0);}finally{solver.dispose();}
});
it('keeps every feasible spectral point when seeding tangent planes',async()=>{
 const {seedQuadraticTangents}=await import('../../src/app/workers/core/sanger-quadratic-envelope'),q={name:'spectral',terms:[[[0,1],[1,-1]]as[number,number][]],upper:.04},cuts=seedQuadraticTangents(q,[.5,.5],[{columns:[0,1],mass:1}]);expect(cuts.length).toBe(2);for(let j=0;j<=20;j++){const x=[.4+j*.01,.6-j*.01];for(const c of cuts)expect(c.terms.reduce((n,[i,w])=>n+w*x[i],0)).toBeLessThanOrEqual(c.upper+1e-9);}
});

it('cuts at the residual boundary even when a trial point is far outside it',async()=>{
 const {quadraticCut}=await import('../../src/app/workers/core/sanger-quadratic-envelope'),q={name:'correlated residual',terms:[[[0,2],[1,-1]]as[number,number][],[[1,3]]as[number,number][]],upper:.04},cut=quadraticCut(q,[10,20]);expect(cut.terms.reduce((s,[i,w])=>s+w*[10,20][i],0)).toBeGreaterThan(cut.upper);for(let i=-20;i<=20;i++)for(let j=-20;j<=20;j++){const x=[i/100,j/100];if((2*x[0]-x[1])**2+9*x[1]**2<=q.upper)expect(cut.terms.reduce((s,[i,w])=>s+w*x[i],0)).toBeLessThanOrEqual(cut.upper+1e-12);}
});
