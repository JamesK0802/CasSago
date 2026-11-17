import {it,expect} from 'vitest';
import {jointBoundCertificates} from '../../src/app/workers/core/sanger-joint-bound-certificate';
import type {MarginalProgram} from '../../src/app/workers/core/sanger-marginal-lp';
it('never presents an unfinished native LP proposal as completed selected-group detail',()=>{
 const p:MarginalProgram={upper:[1],cost:[0],rows:[]};let calls=0;
 const optimizer:any={createModel:()=>({options:{set:()=>{}},passModel:()=>{},zeroAllClocks:()=>{},run:()=>{calls++;return{modelStatus:13};},getSolution:()=>({rowDual:[]}),dispose:()=>{}})},bound=jointBoundCertificates(p,optimizer);
 try{expect(bound.bound([[0,1]])).toBeUndefined();expect(bound.bound([[0,-1]])).toBeUndefined();expect(calls).toBe(2);}finally{bound.dispose();}
});
it('requests the original-model fallback when a native proposal workspace fails',()=>{
 const p:MarginalProgram={upper:[1],cost:[0],rows:[]},bound=jointBoundCertificates(p,{createModel:()=>{throw Error('Unavailable native workspace');}}as any);
 try{expect(bound.bound([[0,-1]])).toBeUndefined();expect(bound.stats.failures).toBe(1);}finally{bound.dispose();}
});
it('adds only implied proposal support planes and leaves the original cone model intact',()=>{
 const p:MarginalProgram={upper:[1],cost:[1],rows:[{terms:[[0,1]],lower:0,upper:1}],quadratics:[{name:'spectrum',terms:[[[0,1]]],upper:1}]},before=JSON.stringify(p);let rows=0;
 const optimizer:any={createModel:()=>({options:{set:()=>{}},passModel:(m:any)=>{rows=m.numRows;},zeroAllClocks:()=>{},run:()=>({modelStatus:7}),getSolution:()=>({rowDual:Array(rows).fill(0)}),dispose:()=>{}})},bound=jointBoundCertificates(p,optimizer);
 try{expect(Number.isFinite(bound.bound([[0,-1]]))).toBe(true);expect(rows).toBeGreaterThan(p.rows.length);expect(JSON.stringify(p)).toBe(before);}finally{bound.dispose();}
});

it('retries the identical LP to completion without accepting the stopped first proposal',()=>{
 const p:MarginalProgram={upper:[1],cost:[0],rows:[]};let calls=0;const limits:number[]=[];
 const optimizer:any={createModel:()=>({options:{set:(v:any)=>{if(v.time_limit!==undefined)limits.push(v.time_limit);}},passModel:()=>{},zeroAllClocks:()=>{},clearSolver:()=>{},run:()=>({modelStatus:++calls===1?13:7}),getSolution:()=>({rowDual:[]}),dispose:()=>{}})},bound=jointBoundCertificates(p,optimizer);
 try{expect(bound.bound([[0,-1]])).toBeLessThanOrEqual(-1);expect(calls).toBe(2);expect(bound.stats.retries).toBe(1);expect(bound.stats.partial).toBe(0);expect(limits).toContain(30);}finally{bound.dispose();}
});

it('falls back on an unusable optimal-status dual instead of calling its box-only interval finished detail',()=>{
 const p:MarginalProgram={upper:[1],cost:[0],rows:[{terms:[[0,1]],lower:0,upper:1}]},optimizer:any={createModel:()=>({options:{set:()=>{}},passModel:()=>{},zeroAllClocks:()=>{},run:()=>({modelStatus:7}),getSolution:()=>({rowDual:[NaN]}),dispose:()=>{}})},bound=jointBoundCertificates(p,optimizer);
 try{expect(bound.bound([[0,-1]])).toBeUndefined();expect(bound.stats.invalidDuals).toBe(1);}finally{bound.dispose();}
});
