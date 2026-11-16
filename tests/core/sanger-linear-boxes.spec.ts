import {it,expect} from 'vitest';
import {tightenAuxiliaryBoxes} from '../../src/app/workers/core/sanger-linear-boxes';
import {FactorizedMarginalSolver} from '../../src/app/workers/core/sanger-factorized-lp';
import {MarginalLinearSolver,LP_INFINITY,validMarginalSolution,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
import loadCommon from '../probes/load-common-optimizer';
it('uses a shared budget to tighten huge auxiliary boxes without changing physical columns or rows',()=>{
 const p:MarginalProgram={upper:[1,1000,500],cost:[0,0,0],rows:[{terms:[[1,1],[2,1]],lower:0,upper:.25}]},rows=JSON.stringify(p.rows);tightenAuxiliaryBoxes(p,1);expect(p.upper[0]).toBe(1);expect(p.upper[1]).toBeCloseTo(.25,9);expect(p.upper[2]).toBeCloseTo(.25,9);expect(JSON.stringify(p.rows)).toBe(rows);for(let i=0;i<=25;i++)for(let j=0;j<=25-i;j++)expect(validMarginalSolution(p,[.7,i/100,j/100])).toBe(true);
});
it('preserves feasible points for signed lower/upper rows and leaves infinite sides alone',()=>{
 const p:MarginalProgram={upper:[1,5,5],cost:[0,0,0],rows:[{terms:[[1,2],[2,-1]],lower:-1,upper:2},{terms:[[0,-1],[1,1],[2,1]],lower:-LP_INFINITY,upper:1}]},original=structuredClone(p),points:number[][]=[];for(let i=0;i<=10;i++)for(let j=0;j<=20;j++)for(let k=0;k<=20;k++){const x=[i/10,j/10,k/10];if(validMarginalSolution(original,x))points.push(x);}tightenAuxiliaryBoxes(p,1);expect(points.length).toBeGreaterThan(50);for(const x of points)expect(validMarginalSolution(p,x)).toBe(true);
});
it('retains all physical alternatives when master cost links have a much larger box sum than their shared budget',async()=>{
 const n=100,p:MarginalProgram={upper:Array(n+1).fill(1),cost:Array(n+1).fill(1),rows:[{terms:Array.from({length:n+1},(_,i)=>[i,1]as[number,number]),lower:0,upper:.2}],localBlocks:[Array.from({length:n},(_,i)=>i),[n]]},before=JSON.stringify(p),optimizer=await loadCommon(),factor=new FactorizedMarginalSolver(p,q=>new MarginalLinearSolver({...q,localBlocks:undefined},optimizer));try{const b=(factor as any).base;expect(Math.max(...b.upper)).toBeLessThan(.200001);for(const i of [0,50,100]){const a=factor.solve([[i,1]],true);expect(a.value).toBeCloseTo(.2,6);expect(validMarginalSolution(p,a.x)).toBe(true);}expect(JSON.stringify(p)).toBe(before);}finally{factor.dispose();}
});
