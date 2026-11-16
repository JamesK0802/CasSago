import {it,expect} from 'vitest';
import {emissionUncertainty} from '../../src/app/workers/core/sanger-emission-uncertainty';
import {validMarginalSolution,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
it('retains physical mass-response products and enforces exact WT photon conservation',()=>{
 const p:MarginalProgram={upper:[.5,.5],cost:[0,0],rows:[]},f=emissionUncertainty(p,[.5,.5],['A','A'],[[[1,0,0,0]],[[1,0,0,0]]],[[[-.4,.4]],[[-.4,.4]]],new Map()),a=f([[0,1]],0,0,0),b=f([[1,1]],1,0,0),value=(form:typeof a,x:number[])=>form.reduce((s,[i,w])=>s+w*x[i],0);
 for(const ma of [0,.1,.3,.5])for(const mb of [0,.2,.5])for(const delta of [-.4,-.1,0,.2,.4]){const x=[ma,mb,delta+.4,-delta+.4,ma*(delta+.4),mb*(-delta+.4)];expect(validMarginalSolution(p,x)).toBe(true);expect(value(a,x)).toBeCloseTo(ma*delta,12);expect(value(b,x)).toBeCloseTo(-mb*delta,12);}
 const wt=[.5,.5,.6,.2,.3,.1];expect(value(a,wt)+value(b,wt)).toBeCloseTo(0,12);expect(validMarginalSolution(p,[.5,.5,.6,.6,.3,.3])).toBe(false);
 const seed=p.upper.map(()=>0);seed[0]=.3;seed[1]=.2;f.seed(seed);expect(validMarginalSolution(p,seed)).toBe(true);expect(value(a,seed)).toBeCloseTo(0,12);expect(value(b,seed)).toBeCloseTo(0,12);
});
