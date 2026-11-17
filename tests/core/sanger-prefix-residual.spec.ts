import {it,expect} from 'vitest';
import {prefixPhotonVariance} from '../../src/app/workers/core/sanger-prefix-residual';
import {mixPhase,phaseNeighbor} from '../../src/app/workers/core/sanger-signal';
it.each([0,.15,-.15])('uses the same phased absolute foreground before normalization for lag %s',lag=>{
 const responses=[[[2,0,0,0],[0,4,0,0],[0,0,1,0],[0,0,0,3]],[[0,1,0,0],[0,0,3,0],[0,0,0,2],[5,0,0,0]]],shares=[.3,.7],transfer={apply:(h:number,p:number,d:number[])=>d.map((v,b)=>v*(1+h+.1*p+.05*b))},foreground=(p:number)=>[0,1,2,3].map(b=>responses.reduce((s,r,h)=>s+shares[h]*transfer.apply(h,p,r[p])[b],0)),raw=responses[0].map((_,p)=>{const q=Math.max(0,Math.min(3,phaseNeighbor(p,lag))),v=mixPhase(foreground(p),foreground(q),lag),total=v.reduce((a,b)=>a+b,0);return v.map(x=>x/total);});
 expect(prefixPhotonVariance([1,2],raw,responses,shares,transfer,lag)).toBeLessThan(1e-28);
 if(lag)expect(prefixPhotonVariance([1,2],raw,responses,shares,transfer,0)).toBeGreaterThan(1e-4);
 raw[1][0]+=.02;raw[1][1]-=.02;
 expect(prefixPhotonVariance([1,2],raw,responses,shares,transfer,lag)).toBeCloseTo(.0001,12);
});
