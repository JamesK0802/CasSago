import {it,expect} from 'vitest';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import {localBoundCertificates} from '../../src/app/workers/core/sanger-local-bound-certificate';
import {LP_INFINITY,type MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
it('releases coordinate workspaces without losing bounds when a coordinate is revisited',async()=>{
 const optimizer=await loadHighs(),create=optimizer.createModel.bind(optimizer);let live=0,peak=0;
 optimizer.createModel=()=>{const native=create();live++;peak=Math.max(peak,live);return new Proxy(native,{get(t,k){if(k==='dispose')return()=>{live--;return t.dispose();};const v=t[k];return typeof v==='function'?v.bind(t):v;}});};
 const count=80,p:MarginalProgram={upper:Array(2*count).fill(1),cost:Array(2*count).fill(0),rows:[]},point:number[]=[],flow=new Map();for(let q=0;q<count;q++){const i=2*q;point.push(0,1);flow.set(q,{unchanged:[{column:i+1}],changed:[i]});p.rows.push({terms:[[i,1],[i+1,1]],lower:1,upper:1},{terms:[[i,1]],lower:-LP_INFINITY,upper:.2});}
 const bounds=localBoundCertificates(p,point,[{flow}]as any,optimizer);try{for(const q of [...Array(count).keys(),0])expect(bounds.bound(q,[[2*q,1]],true)).toBeCloseTo(.2,7);expect(peak).toBeLessThanOrEqual(2);}finally{bounds.dispose();}expect(live).toBe(0);
});
