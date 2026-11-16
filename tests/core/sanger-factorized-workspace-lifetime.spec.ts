import {it,expect} from 'vitest';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import {FactorizedMarginalSolver} from '../../src/app/workers/core/sanger-factorized-lp';
import {MarginalLinearSolver,OptimizationPending,validMarginalSolution,LP_INFINITY,type MarginalProgram,type LinearTerm} from '../../src/app/workers/core/sanger-marginal-lp';
it('keeps local native workspaces bounded and resumes pending solves without dropping blocks',async()=>{
 const optimizer=await loadHighs(),create=optimizer.createModel.bind(optimizer);let live=0,peak=0,pauses=0;optimizer.createModel=()=>{const model=create();live++;peak=Math.max(peak,live);let closed=false;return new Proxy(model,{get(t,k){if(k==='dispose')return()=>{if(!closed){closed=true;live--;}return t.dispose();};const v=t[k];return typeof v==='function'?v.bind(t):v;}});};
 const count=30,p:MarginalProgram={upper:Array(1+2*count).fill(1),cost:Array(1+2*count).fill(0),rows:[{terms:[[0,1]],lower:.4,upper:.4}],localBlocks:[]},terms:LinearTerm[]=[];for(let q=0;q<count;q++){const i=1+2*q;p.localBlocks!.push([i,i+1]);p.rows.push({terms:[[i,1],[i+1,1],[0,-1]],lower:0,upper:0},{terms:[[i,1]],lower:.1,upper:LP_INFINITY});terms.push([i,1]);}
 const factor=new FactorizedMarginalSolver(p,program=>{const solver=new MarginalLinearSolver({...program,localBlocks:undefined},optimizer);let first=true;return new Proxy(solver,{get(t,k){if(k==='solve')return(...args:any[])=>{if(first){first=false;pauses++;throw new OptimizationPending();}return t.solve(...args);};const v=(t as any)[k];return typeof v==='function'?v.bind(t):v;}});});
 expect((factor as any).elastic.every((b:any)=>b.program===undefined)).toBe(true);
 try{const steps=factor.steps(terms);let answer:any;for(let n=0;n<1000;n++){const next=steps.next();if(next.done){answer=next.value;break;}}expect(answer).toBeDefined();expect(answer.value).toBeCloseTo(.1*count,7);expect(validMarginalSolution(p,answer.x)).toBe(true);expect(pauses).toBeGreaterThan(count);expect(peak).toBeLessThanOrEqual(4);expect((factor as any).elastic.every((b:any)=>b.program===undefined)).toBe(true);}finally{factor.dispose();}expect(live).toBe(0);
});
