import{it,expect}from'vitest';
import{SangerTrace,TRACE_BASES}from'../../src/app/models/sanger-analysis';
import{sangerPooledDemo}from'../../src/app/utils/sanger-analysis-demo';
import{analyzeAutomaticSanger}from'../../src/app/workers/core/sanger-automatic';
function trace(molecules:{sequence:string;weight:number}[],name:string):SangerTrace {
 const trim=10,n=Math.min(700,Math.max(...molecules.map(m=>m.sequence.length)))-trim-3,peaks=Uint32Array.from({length:n},(_,i)=>20+12*i),channels=Object.fromEntries(TRACE_BASES.map(b=>[b,new Float32Array(12*n+50)]))as SangerTrace['channels'];
 for(let i=0;i<n;i++)for(const m of molecules){const b=m.sequence[i+trim] as typeof TRACE_BASES[number];if(!b)continue;for(let d=-5;d<=5;d++)channels[b][peaks[i]+d]+=9000*m.weight*Math.exp(-d*d/4);}
 const sequence=Array.from(peaks,p=>TRACE_BASES.reduce((a,b)=>channels[a][p]>=channels[b][p]?a:b)).join('');return {name,sampleName:name,sequence,peaks,quality:null,channels};
}
it('uses the full measured mixed trace after a shorter amplicon and its control template terminate',()=>{
 const d=sangerPooledDemo(),refs=d.amplicons.map((a,h)=>({...a,sequence:h?a.sequence:a.sequence.slice(0,480)})),seq=refs.map(a=>a.sequence.slice(d.primer.length)),cuts=refs.map(a=>a.cut-d.primer.length),control=trace(seq.map(sequence=>({sequence,weight:.5})),'control'),sample=trace([{sequence:seq[0],weight:.24},{sequence:seq[0].slice(0,cuts[0]-1)+seq[0].slice(cuts[0]),weight:.16},{sequence:seq[1],weight:.3},{sequence:seq[1].slice(0,cuts[1])+'GC'+seq[1].slice(cuts[1]),weight:.3}],'sample'),r=analyzeAutomaticSanger({references:refs,primer:d.primer,control,sample,maxDeletion:5,maxInsertion:5});
 expect(r.mixture!.bases).toBe(sample.peaks.length);expect(r.diagnostics!.junctionFit!.positions.at(-1)).toBe(r.pooled!.sampleTrim+sample.peaks.length-1);expect(r.diagnostics!.controlResponse!.refitBases).toBeGreaterThan(seq[0].length);
 expect(r.pooled!.amplicons[0].components.find(c=>c.size===-1)!.fraction).toBeCloseTo(.4,2);expect(r.pooled!.amplicons[1].components.find(c=>c.size===2)!.fraction).toBeCloseTo(.5,2);
},60000);
it('does not certify reference abundance from molecule endpoints when normalized WT colours cannot distinguish the references',()=>{
 const primer='GCTAGTCAGACTCAGTACGA',seq=['ACGT'.repeat(120),'ACGT'.repeat(190)],refs=seq.map((s,h)=>({name:String(h),sequence:primer+s,cut:primer.length+180+h*120})),control=trace(seq.map(sequence=>({sequence,weight:.5})),'control'),sample=trace(seq.map(sequence=>({sequence,weight:.5})),'sample'),r=analyzeAutomaticSanger({references:refs,primer,control,sample,maxDeletion:3,maxInsertion:2});
 expect(r.mixture!.bases).toBe(sample.peaks.length);expect(r.pooled!.preflight.passed).toBe(false);expect(r.mixture!.passesChecks).toBe(false);expect(r.pooled!.amplicons.every(a=>a.components.every(c=>Number.isFinite(c.fraction)))).toBe(true);
},60000);
