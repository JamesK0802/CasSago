import{it,expect}from'vitest';
import loadHighs from'../../public/vendor/highs-1.15.4/highs.mjs';
import{MarginalLinearSolver,OptimizationPending,marginalViolation,type MarginalProgram}from'../../src/app/workers/core/sanger-marginal-lp';
import{FactorizedMarginalSolver}from'../../src/app/workers/core/sanger-factorized-lp';

it('resumes a pending elastic feasibility model without rebuilding its infeasible parent or losing its workspace',async()=>{
 const p:MarginalProgram={upper:[1,.5],cost:[0,0],localBlocks:[[1]],rows:[{terms:[[1,1],[0,-1]],lower:.2,upper:1}]},optimizer=await loadHighs();
 let interruptions=0,secondaryModels=0;const secondaryCalls:number[]=[];
 const factor=new FactorizedMarginalSolver(p,q=>{
  const solver=new MarginalLinearSolver({...q,localBlocks:undefined},optimizer);
  // This block's hard physical row is infeasible at the first master point.
  // Its all-elastic phase I adds a slack for each finite row side.
  if(q.upper.length>1&&q.cost.some(c=>c===1)){
   const model=++secondaryModels,solve=solver.solve.bind(solver);solver.solve=(...args)=>{
    secondaryCalls.push(model);if(interruptions++<2)throw new OptimizationPending();return solve(...args);
   };
  }
  return solver;
 });
 try{
  const steps=factor.steps([[0,1]],true);let answer:any,yields=0;
  for(;;){const next=steps.next();if(next.done){answer=next.value;break;}yields++;expect(yields).toBeLessThan(20);}
  expect(answer.value).toBeCloseTo(.3,7);expect(marginalViolation(p,answer.x)).toBeLessThanOrEqual(1e-7);
  expect(secondaryCalls.slice(0,3)).toEqual([1,1,1]);
 }finally{factor.dispose();}
});
