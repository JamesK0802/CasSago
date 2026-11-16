import {it,expect} from 'vitest';
import {referenceTransfer} from '../../src/app/workers/core/sanger-reference-transfer';
it('learns per-reference dye response and read trend only from the unedited training prefix',()=>{
 const shares=[.4,.6],responses=shares.map((_,h)=>Array.from({length:120},(_,p)=>[0,1,2,3].map(b=>b===(h?(3*p+1)%4:p%4)?1:0))),alpha=[[.08,-.03,.02,.04],[-.04,.07,-.02,.03]],beta=[.03,-.02],at=(p:number)=>2*p/79-1,sample=Array.from({length:120},(_,p)=>[0,1,2,3].map(b=>shares.reduce((n,w,h)=>n+w*responses[h][p][b]*Math.exp(alpha[h][b]+beta[h]*at(p)),0))),positions=Array.from({length:80},(_,p)=>p),initial={gains:[1,1,1,1],slope:0,logScale:0,first:0,span:79},fit=referenceTransfer(responses,shares,sample,positions,0,1,initial);
 for(let p=0;p<120;p++)for(let b=0;b<4;b++)expect(shares.reduce((n,w,h)=>n+w*fit.apply(h,p,responses[h][p])[b],0)).toBeCloseTo(sample[p][b],6);
 const altered=sample.map((v,p)=>p>=80?v.map(x=>100*x):v),again=referenceTransfer(responses,shares,altered,positions,0,1,initial);expect(again.parameters).toEqual(fit.parameters);
});
it('keeps finite ADC uncertainty even for a perfectly matching paired prefix',()=>{
 const responses=[Array.from({length:80},(_,p)=>[0,1,2,3].map(b=>Number(b===p%4)))],areas=responses[0].map(r=>r.map(v=>10000*v)),fit=referenceTransfer(responses,[1],areas,Array.from({length:80},(_,p)=>p),0,10000,{gains:[1,1,1,1],slope:0,logScale:0,first:0,span:79});expect(fit.logVariance).toBeGreaterThan(0);expect(fit.quantizationLogVariance).toBeCloseTo(2e-8,12);expect(fit.logVariance).toBeGreaterThanOrEqual(fit.quantizationLogVariance);
});
