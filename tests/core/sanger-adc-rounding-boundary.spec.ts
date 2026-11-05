import{it,expect}from'vitest';
import loadHighs from'../../public/vendor/highs-1.15.4/highs.mjs';
import{rawIntervalScorer}from'../../src/app/workers/core/sanger-raw-interval-fit';
import type{SangerTrace}from'../../src/app/models/sanger-analysis';
it('does not collapse two adjacent ADC bins to their shared excluded upper boundary',async()=>{
 let seed=78947;let reference=Array.from({length:320},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return'ACGT'[seed>>>30];}).join('');reference=reference.slice(0,100)+'A'+reference.slice(101,110)+'A'+reference.slice(111);
 const mutate=(positions:number[])=>[...reference].map((b,p)=>positions.includes(p)?'C':b).join(''),truth=[reference,mutate([100]),mutate([110])],trace=(sequences:string[],weights:number[]):SangerTrace=>{const peaks=Uint32Array.from({length:300},(_,i)=>14+i*10),channels=Object.fromEntries([...'ACGT'].map(b=>[b,new Float32Array(3030).fill(2)]))as SangerTrace['channels'];for(let p=0;p<300;p++)for(const b of'ACGT'){const mass=sequences.reduce((n,s,i)=>n+weights[i]*Number(s[p]===b),0);for(let d=-4;d<=4;d++)channels[b as'A'][peaks[p]+d]+=Math.round(1000*mass*Math.exp(-d*d/4));}return{name:'measured',sampleName:'measured',sequence:'A'.repeat(300),peaks,quality:null,channels};},control=trace([reference],[1]),sample=trace(truth,[.799,.10049,.10051]),positions=Array.from({length:260},(_,i)=>i+20),optimizer=await loadHighs(),score=rawIntervalScorer(control,sample,[reference],positions,0,0,[2,2,2,2],[2,2,2,2],{amplitude:1000,sigma:4,coefficients:[1000]},optimizer)!;
 expect(sample.channels.C[sample.peaks[100]]).toBe(102);expect(sample.channels.C[sample.peaks[110]]).toBe(103);
 const collapsed=score([reference,mutate([100,110])],[.8995,.1005])!;
 expect(collapsed.violation).toBeGreaterThan(1e-6);
 const recovered=score(truth,[.799,.10049,.10051])!;
 expect(recovered.violation).toBeLessThan(1e-7);
});
