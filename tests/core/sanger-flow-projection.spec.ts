import {it,expect} from 'vitest';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import {MarginalLinearSolver,validMarginalSolution,LP_INFINITY,type MarginalProgram,type LinearTerm} from '../../src/app/workers/core/sanger-marginal-lp';
import {expandInsertionTransports} from '../../src/app/workers/core/sanger-insertion-transport';
it('preserves every linear optimum when replacing Hall inequalities with their mass-balanced form',async()=>{
 const optimizer=await loadHighs(),row=(terms:LinearTerm[],lower:number,upper=lower)=>({terms,lower,upper}),base=[row([[0,1],[1,1],[2,1],[3,1],[4,1],[5,1]],1),row([[0,1]],0,.4),row([[1,1]],0,.6),...[[[0,1],[2,1]],[[1,1],[3,1]],[[4,1]],[[5,1]]].map((terms,i)=>row(terms as LinearTerm[],[.25,.35,.1,0][i],[.4,.6,.25,.15][i]))],old=base.slice(),reduced=base.slice();
 for(let b=0;b<4;b++){const expected=[.4,.6,0,0][b];old.push(row([[b+2,1],...(b!==0?[[0,1]as LinearTerm]:[]),...(b!==1?[[1,1]as LinearTerm]:[])],-LP_INFINITY,1-expected));reduced.push(row([...(b===0?[[0,-1]as LinearTerm]:[]),...(b===1?[[1,-1]as LinearTerm]:[]),...[0,1,2,3].flatMap(a=>a!==b?[[a+2,-1]as LinearTerm]:[])],-LP_INFINITY,-expected));}
 const a=new MarginalLinearSolver({upper:Array(6).fill(1),cost:Array(6).fill(0),rows:old},optimizer),b=new MarginalLinearSolver({upper:Array(6).fill(1),cost:Array(6).fill(0),rows:reduced},optimizer);
 try{for(let k=0;k<16;k++){const terms=Array.from({length:6},(_,i)=>[i,Math.sin(7*k+3*i)]as LinearTerm);for(const maximize of [false,true])expect(a.solve(terms,maximize).value).toBeCloseTo(b.solve(terms,maximize).value,8);}}finally{a.dispose();b.dispose();}
});
it('lifts shared insertion colours back to origin-specific flows without losing mass or imposing phase',()=>{
 const p:MarginalProgram={upper:Array(6).fill(1),cost:Array(6).fill(0),rows:[{terms:[[0,1],[1,1]],lower:1,upper:1},{terms:[[2,1],[3,1],[4,1],[5,1],[0,-1],[1,-1]],lower:0,upper:0}]},point=[.25,.75,.45,0,.55,0],edges=[0,1].map((mass,i)=>({reference:0,target:0,incoming:0,outgoing:2,size:2,start:i,mass,free:new Map([[1,[2,3,4,5]]])})),g:any={id:'g',sharedInsertions:true,edges,mass:[[0,1],[1,1]]},target:any={program:p,point,groups:[g],edges,positions:[1]},before=JSON.stringify(target),expanded=expandInsertionTransports(target,g);expect(expanded.target.program.rows[0]).not.toBe(p.rows[0]);expect(expanded.target.program.rows[0].terms).toBe(p.rows[0].terms);
 expect(validMarginalSolution(expanded.target.program,expanded.target.point)).toBe(true);expect(expanded.group.edges[0].free.get(1)).not.toEqual(expanded.group.edges[1].free.get(1));for(let b=0;b<4;b++)expect(expanded.group.edges.reduce((n,e)=>n+expanded.target.point[e.free.get(1)![b]],0)).toBeCloseTo(point[b+2],10);expect(JSON.stringify(target)).toBe(before);
});
