import{it,expect}from'vitest';import loadHighs from'../../public/vendor/highs-1.15.4/highs.mjs';import{SangerSequenceFitCache}from'../../src/app/workers/core/sanger-sequence-fit-cache';import{rawIntervalScorer}from'../../src/app/workers/core/sanger-raw-interval-fit';import type{SangerTrace}from'../../src/app/models/sanger-analysis';
it.each([true,false])('learns a nonstandard response and rejects inconsistent ADC dictionaries (persistent %s)',async persistent=>{
 const highs=await loadHighs(),optimizer=persistent?highs:{solve:highs.solve.bind(highs)};
 let seed=61237;const dna=()=>Array.from({length:330},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return'ACGT'[seed>>>30];}).join(''),refs=[dna(),dna(),dna()],edited=refs.map(s=>s.slice(0,160)+(s[160]==='A'?'C':'A')+s.slice(161)),sequences=refs.flatMap((s,i)=>[s,edited[i]]),trace=(seqs:string[],weights:number[]):SangerTrace=>{const peaks=Uint32Array.from({length:300},(_,i)=>14+i*10),channels=Object.fromEntries([...'ACGT'].map(b=>[b,new Float32Array(3030).fill(9)]))as SangerTrace['channels'];for(let p=0;p<300;p++)for(const b of 'ACGT'){const mass=seqs.reduce((sum,s,i)=>sum+weights[i]*Number(s[p]===b),0);for(let d=-4;d<=4;d++)channels[b as'A'][peaks[p]+d]+=Math.round(947*mass*Math.exp(-d*d/3.47));}return{name:'measured',sampleName:'measured',sequence:'A'.repeat(300),peaks,quality:null,channels};},control=trace(refs,[.2,.3,.5]),sample=trace(sequences,[.11,.09,.17,.13,.31,.19]),positions=Array.from({length:260},(_,i)=>i+20),score=rawIntervalScorer(control,sample,refs,positions,0,0,[9,9,9,9],[9,9,9,9],{sigma:3.46,amplitude:948,coefficients:[190,284,474]},optimizer);
 expect(score).not.toBeNull();const wrong=score!(refs,[.2,.3,.5])!;expect(wrong.violation).toBeGreaterThan(.1);expect(wrong.loss).toBeGreaterThan(.1);expect(wrong.lowerBoundOnly).not.toBe(true);expect(wrong.activePositions).toContain(160);expect(wrong.loss).toBeCloseTo(wrong.violation,10);
 const recovered=score!(sequences,sequences.map(()=>1/6))!;expect(recovered.violation).toBeLessThan(1e-6);expect(recovered.sigma).toBeCloseTo(3.47,2);expect(recovered.amplitude).toBeCloseTo(947,0);expect(recovered.weights.reduce((s,v)=>s+v,0)).toBeCloseTo(1,10);
 const point=positions.flatMap(p=>[...'ACGT'].map(b=>sequences.reduce((sum,seq,i)=>sum+[.11,.09,.17,.13,.31,.19][i]*Number(seq[p]===b),0)+.0003)),bounded=rawIntervalScorer(control,sample,refs,positions,0,0,[9,9,9,9],[9,9,9,9],{sigma:3.46,amplitude:948,coefficients:[190,284,474]},optimizer,point)!;
 // An imperfect starting model must remain evaluable: an exact-compatibility
 // rejection is not a valid pruning rule while the incumbent has positive loss.
 expect(bounded.maximumCentroidMse).toBe(Infinity);
 const evolving=new SangerSequenceFitCache([...refs,...edited],positions,point,512,bounded),start=evolving.fit([0,1,2]);
 expect(start.converged).toBe(true);expect(start.selectionLoss).toBeGreaterThan(.1);
 expect(bounded.maximumCentroidMse).toBeGreaterThanOrEqual(start.rmse**2-1e-15);
 const improved=evolving.fit([0,3,1,4,2,5]);expect(improved.converged).toBe(true);expect(improved.selectionLoss).toBeLessThan(1e-6);
 expect(Number.isFinite(bounded.maximumCentroidMse)).toBe(true);

 sample.channels.A[sample.peaks[150]+2]+=30;expect(rawIntervalScorer(control,sample,refs,positions,0,0,[9,9,9,9],[9,9,9,9],{sigma:3.46,amplitude:948,coefficients:[190,284,474]},optimizer)).toBeNull();
},60000);
