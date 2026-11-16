import {it,expect} from 'vitest';
import loadCommonOptimizer from '../probes/load-common-optimizer';
import {fixedGeometryWitness} from '../../src/app/workers/core/sanger-fixed-geometry-witness';
import {validMarginalSolution,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
it('finds full-model witnesses while leaving the parent geometry and constraints unchanged',async()=>{
 const p:MarginalProgram={upper:[1,1,1],cost:[0,0,0],rows:[{terms:[[1,1],[2,1],[0,-1]],lower:0,upper:0}],quadratics:[{name:'geometry',terms:[[[0,1]]],upper:.25}]},before=JSON.stringify(p),optimizer=await loadCommonOptimizer();for(const value of [0,.5]){const x=fixedGeometryWitness(p,[.5,.25,.25],[1,2],[[1,1]],value,optimizer);expect(x).toBeDefined();expect(validMarginalSolution(p,x!)).toBe(true);expect(x![0]).toBe(.5);expect(x![1]).toBeCloseTo(value,7);}expect(JSON.stringify(p)).toBe(before);
});
it('never accepts an unreachable endpoint or a witness violating the original cone',async()=>{
 const p:MarginalProgram={upper:[1,1],cost:[0,0],rows:[{terms:[[1,1],[0,-1]],lower:0,upper:0}],quadratics:[{name:'geometry',terms:[[[0,1]]],upper:.25}]},optimizer=await loadCommonOptimizer();expect(fixedGeometryWitness(p,[.5,.5],[1],[[1,1]],0,optimizer)).toBeUndefined();expect(fixedGeometryWitness(p,[1,1],[1],[[1,1]],1,optimizer)).toBeUndefined();expect(fixedGeometryWitness(p,[.5,.5],[0,1],[[1,1]],0,optimizer)).toBeUndefined();
});
