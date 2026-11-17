import {it,expect} from 'vitest';
import {meanQuietBackground} from '../../src/app/workers/core/sanger-mean-background';
import {cycleFluorescence} from '../../src/app/workers/core/sanger-cycle-fluorescence';
import {fitGroupedIntensities} from '../../src/app/workers/core/sanger-grouped-mixture';
it('subtracts the expected rectified noise without reclipping its negative residuals',()=>{
 const n=130,peaks=Array.from({length:n},(_,i)=>14+10*i),sequence='ACGT'.repeat(40),channels=Object.fromEntries([...'ACGT'].map(b=>[b,Float64Array.from({length:n*10+30},(_,j)=>j%2?4:0)])) as any;
 for(let p=0;p<n;p++)for(let d=-4;d<=4;d++)channels[sequence[p]][peaks[p]+d]+=1200*Math.exp(-d*d/4);
 const trace={name:'noise fixture',peaks,sequence:sequence.slice(0,n),channels,quality:null} as any,upstream=Array.from({length:100},(_,i)=>i+10),background=meanQuietBackground(trace,[sequence],0,upstream);expect(background.background).toEqual([2,2,2,2]);
 const signed=cycleFluorescence(trace,background.background,true),clipped=cycleFluorescence(trace,background.background);for(let p=12;p<100;p++){const quiet=[...'ACGT'].find(b=>[-1,0,1].every(d=>sequence[p+d]!==b))!,b='ACGT'.indexOf(quiet);expect(signed.areas[p][b]).toBeCloseTo(0,12);expect(clipped.areas[p][b]).toBeGreaterThan(5);}
 const altered={...trace,channels:Object.fromEntries([...'ACGT'].map(b=>[b,Float64Array.from(channels[b],(v,j)=>j>peaks[111]?v+500:v)]))};expect(meanQuietBackground(altered as any,[sequence],0,upstream).background).toEqual(background.background);
});
it.each([Infinity,.06])('weighted ratio fitting matches repeating identical measurements at robust delta %s',delta=>{
 const columns=[[2,.1,.1,.1,.1,1,.1,.1,.2,.1,2,.1,.1,.1,.1,2],[.1,2,.1,.1,.1,.1,1,.1,2,.1,.1,.1,.1,2,.1,.1]],weights=[1,3,2,1],observed=[.68,.23,.05,.04,.05,.56,.35,.04,.2,.05,.68,.07,.08,.25,.05,.62],repeat=(v:number[])=>weights.flatMap((w,p)=>Array.from({length:w},()=>v.slice(4*p,4*p+4)).flat()),weighted=fitGroupedIntensities(columns,observed,[0,0],[1],delta,undefined,undefined,weights),expanded=fitGroupedIntensities(columns.map(repeat),repeat(observed),[0,0],[1],delta);expect(weighted.converged).toBe(true);expect(expanded.converged).toBe(true);weighted.weights.forEach((w,i)=>expect(w).toBeCloseTo(expanded.weights[i],7));
});
