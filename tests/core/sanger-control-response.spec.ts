import{it,expect}from'vitest';
import{controlKnownContextEmission,controlContextProbability,ControlResponseProfile}from'../../src/app/workers/core/sanger-control-response';
import{fitEmissionMixture}from'../../src/app/workers/core/sanger-emission';
const bases='ACGT',table=()=>Array.from({length:64},(_,c)=>Array.from(bases,(_,b)=>.01*(c%5+1)+(Math.floor(c/4)%4===b?1:.02)));
it('reuses the same context/read-cycle response and preserves the absence of fluorescence after a molecule terminates',()=>{
 const profile:ControlResponseProfile={table:table(),endTable:table().map(v=>v.map(x=>x*.7)),start:10,end:180,sites:{80:{6:[.1,.2,.3,.4]}}},pi=(s:string)=>Array.from(bases,b=>s==='N'?.25:Number(s===b));
 for(const p of [5,80,190])for(const l of ['',...bases])for(const b of ['',...bases])for(const r of ['',...bases]){
  const expected=controlContextProbability(pi(l),pi(b),pi(r),profile,p),actual=controlKnownContextEmission(l,b,r,profile,p);
  actual.forEach((v,d)=>expect(v).toBeCloseTo(expected[d],12));if(!b)expect(actual).toEqual([0,0,0,0]);
 }
});
it('recovers molecule fractions with unequal context yields and early product termination on all observed cycles',()=>{
 const rows=120,columns=[0,1].map(h=>Array.from({length:rows},(_,p)=>Array.from(bases,(_,b)=>h&&p>=105?0:Number(b===(p+(h&&p>=35?1:0))%4)*(h?.65+.1*Math.cos(p*.2):1.4+.3*Math.sin(p*.13)))).flat()),truth=[.62,.38],observed=Array.from({length:rows},(_,p)=>{const a=Array.from(bases,(_,b)=>columns.reduce((s,c,h)=>s+truth[h]*c[4*p+b],0)),total=a.reduce((s,x)=>s+x,0);return a.map(x=>x/total);}).flat(),fit=fitEmissionMixture(columns,observed,undefined,32,[.99,.01]);
 expect(fit.weights[0]).toBeCloseTo(truth[0],6);expect(fit.weights[1]).toBeCloseTo(truth[1],6);expect(fit.r2).toBeGreaterThan(.99999999);expect(fit.converged).toBe(true);expect(columns[0]).toHaveLength(observed.length);
});

it('solves weighted, nonnegative molecule shares without bias from nearly redundant candidates',async()=>{
 const {simplexFit}=await import('../../src/app/workers/core/sanger-signal');let seed=710349;const rng=()=>((seed=(Math.imul(seed,1664525)+1013904223)>>>0)+.5)/4294967296;
 const columns=Array.from({length:5},()=>Array.from({length:70},()=>{const v=Array.from(bases,()=>.1+rng()),s=v.reduce((a,b)=>a+b,0);return v.map(x=>x/s);}).flat()),truth=[.05,.1,.15,.2,.5],observed=columns[0].map((_,k)=>columns.reduce((sum,c,i)=>sum+truth[i]*c[k],0)),precision=Array.from({length:70},(_,p)=>p%5===0?.25:1),fit=simplexFit(columns,observed,precision);
 expect(fit.rmse).toBeLessThan(1e-5);expect(fit.optimalityGap!).toBeLessThanOrEqual(1e-10);truth.forEach((w,i)=>expect(fit.weights[i]).toBeCloseTo(w,5));
});

it('learns relative dye strength from absolute WT peaks even when normalized WT colours contain no gain information',async()=>{
 const {fitControlTransfer}=await import('../../src/app/workers/core/sanger-control-transfer');let seed=84219;const rng=()=>((seed=(Math.imul(seed,1664525)+1013904223)>>>0)+.5)/4294967296,sequence=Array.from({length:240},()=>bases[Math.floor(rng()*4)]).join(''),profile:ControlResponseProfile={table:Array.from({length:64},(_,id)=>Array.from(bases,(_,b)=>Number(b===Math.floor(id/4)%4))),start:0,end:240},raw=[1.25,.9,.8,1.1],mean=raw.reduce((s,x)=>s+Math.log(x),0)/4,truth=raw.map(x=>Math.exp(Math.log(x)-mean)),areas=Array.from(sequence,b=>Array.from(bases,(d,i)=>Number(b===d)*truth[i]*9000)),positions=Array.from({length:190},(_,p)=>p+10),fit=fitControlTransfer([sequence],positions,areas,profile,0,[1,1,1,1],[1]);
 expect(fit.used).toBe(true);truth.forEach((g,i)=>expect(Math.abs(fit.gains[i]-g)).toBeLessThan(.01));expect(fit.rawMse).toBeLessThan(fit.initialRawMse*.1);expect(fit.colourMse).toBeLessThan(1e-12);
});
