import {describe,it,expect} from 'vitest';import {referenceBackgroundPeaks,integratedPeaks,simplexFit,calibrateSignal,signalWeights} from '../../src/app/workers/core/sanger-signal';import {sangerAnalysisDemo} from '../../src/app/utils/sanger-analysis-demo';import {parseAbif} from '../../src/app/utils/abif';
const rng=(seed:number)=>()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/4294967296;};
describe('trace signal optimisation',()=>{
 it('resolves rare homeolog markers without a large shared signal slowing convergence',()=>{
  const columns=Array.from({length:3},(_,h)=>Array.from({length:2400},(_,k)=>k<2388?Number(k%4===0):Number(k-2388===h*4)));
  const truth=[.2,.3,.5],observed=columns[0].map((_,k)=>columns.reduce((s,c,h)=>s+truth[h]*c[k],0));
  const fit=simplexFit(columns,observed);
  expect(fit.converged).toBe(true);expect(fit.optimalityGap!).toBeLessThanOrEqual(1e-10);
  fit.weights.forEach((v,h)=>expect(v).toBeCloseTo(truth[h],5));
 });
 it('preserves the simplex solution when the same background vector is added to all models and observations',()=>{
  const columns=[[1,0,0,0,0,1,0,0],[0,1,0,0,0,0,1,0],[0,0,1,0,0,0,0,1]],truth=[.1,.6,.3];
  const observed=columns[0].map((_,k)=>columns.reduce((s,c,h)=>s+truth[h]*c[k],0)),common=observed.map((_,k)=>50+k*.7);
  const before=simplexFit(columns,observed,[1,.3]),after=simplexFit(columns.map(c=>c.map((v,k)=>v+common[k])),observed.map((v,k)=>v+common[k]),[1,.3]);
  after.weights.forEach((v,h)=>expect(v).toBeCloseTo(before.weights[h],6));
  expect(after.converged).toBe(true);
 });
 it('fits normalized nonnegative fractions directly on the simplex',()=>{const c=[[1,0,0,0,0,1,0,0],[0,1,0,0,0,0,1,0]],y=c[0].map((v,i)=>.7*v+.3*c[1][i]),fit=simplexFit(c,y);expect(fit.weights[0]).toBeCloseTo(.7,5);expect(fit.weights.reduce((a,b)=>a+b,0)).toBeCloseTo(1,10);expect(fit.r2).toBeGreaterThan(.9999);});
 it('subtracts only inter-peak background and keeps signal-based weights bounded',()=>{const d=sangerAnalysisDemo(),trace=parseAbif(d.control),signal=integratedPeaks(trace);expect(signal.vectors).toHaveLength(trace.sequence.length);expect(signal.vectors.every(v=>v.every(x=>Number.isFinite(x)&&x>=0))).toBe(true);expect(signal.vectors.every(v=>Math.abs(v.reduce((a,b)=>a+b,0)-1)<1e-6)).toBe(true);const w=signalWeights([signal,signal],[[1,2,3],[1,2,3]]);expect(w.every(x=>x>=.25&&x<=1)).toBe(true);});
 it('learns dye gain/lag from upstream only and does not inspect edited positions',()=>{
  const random=rng(211),templates=Array.from({length:3},()=>Array.from({length:260},()=>{const base=Math.floor(random()*4);return Array.from({length:4},(_,b)=>Number(b===base));})),shares=[.5,.3,.2],gains=[1.3,.8,1.1,.9],lag=.1;
  const observed=Array.from({length:260},(_,p)=>{const v=Array.from({length:4},(_,b)=>gains[b]*templates.reduce((sum,t,h)=>sum+shares[h]*((1-lag)*t[p][b]+lag*(t[p-1]?.[b]||0)),0));return v.map(x=>x/v.reduce((a,b)=>a+b,0));}),positions=Array.from({length:170},(_,i)=>i+1),fit=calibrateSignal(observed,positions,templates);
  expect(fit.used).toBe(true);expect(fit.validationR2).toBeGreaterThan(.99);expect(fit.lag).toBeCloseTo(.1,1);expect(fit.gains.every(g=>g>=.55&&g<=1.8)).toBe(true);
  const changed=observed.map((v,p)=>p>180?[1,0,0,0]:v);expect(calibrateSignal(changed,positions,templates)).toEqual(fit);
 });
 it('keeps identity calibration for a perfect upstream model',()=>{const random=rng(531),template=Array.from({length:150},()=>{const base=Math.floor(random()*4);return Array.from({length:4},(_,b)=>Number(b===base));});const fit=calibrateSignal(template,Array.from({length:130},(_,i)=>i+1),[template]);expect(fit.used).toBe(false);expect(fit.gains).toEqual([1,1,1,1]);expect(fit.lag).toBe(0);});
});

describe('reference-certified background',()=>{
 it('does not treat a weak edited peak as background or learn from edited positions',()=>{const d=sangerAnalysisDemo(),trace=parseAbif(d.control),positions=Array.from({length:180},(_,i)=>i+10),first=referenceBackgroundPeaks(trace,[trace.sequence],0,positions),altered=structuredClone(trace);for(const b of ['A','C','G','T'] as const)for(let p=altered.peaks[230];p<altered.channels[b].length;p++)altered.channels[b][p]+=125;const second=referenceBackgroundPeaks(altered,[trace.sequence],0,positions);expect(first.background).toEqual(second.background);expect(first.backgroundSamples!.every(n=>n>=8)).toBe(true);expect(first.vectors.every(v=>Math.abs(v.reduce((a,b)=>a+b,0)-1)<1e-9)).toBe(true);});
 it('falls back to zero subtraction when a channel has no certified empty observations',()=>{const d=sangerAnalysisDemo(),trace=parseAbif(d.control),refs=['A','C','G','T'].map(b=>b.repeat(trace.sequence.length)),r=referenceBackgroundPeaks(trace,refs,0,Array.from({length:160},(_,i)=>i+10));expect(r.background).toEqual([0,0,0,0]);expect(r.backgroundSamples).toEqual([0,0,0,0]);});
});

it('learns preceding or following cycle contamination from WT only with nonnegative emission',async()=>{
 const {phaseVector}=await import('../../src/app/workers/core/sanger-signal');let seed=90131;const sequence=Array.from({length:360},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return 'ACGT'[seed>>>30];}).join(''),templates=[Array.from(sequence,base=>Array.from('ACGT',b=>Number(b===base)))],positions=Array.from({length:220},(_,i)=>i+2);
 for(const lag of [-.075,.075]){const observed=Array.from({length:sequence.length},(_,p)=>phaseVector(sequence,p,lag)),fit=calibrateSignal(observed,positions,templates,'fine');expect(fit.used).toBe(true);expect(fit.lag).toBeCloseTo(lag,3);expect(fit.validationR2).toBeGreaterThan(.9999);expect(observed.every(v=>v.every(x=>x>=0))).toBe(true);const changed=observed.map((v,p)=>p>250?[0,1,0,0]:v);expect(calibrateSignal(changed,positions,templates,'fine')).toEqual(fit);}
});
