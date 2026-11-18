import type { SangerTrace } from '../../models/sanger-analysis';

const dyes = ['A', 'C', 'G', 'T'] as const;
type Dye = typeof dyes[number];
type Interval = [number, number];

function traceAt(trace: SangerTrace, dye: Dye, index: number, offset: number) {
 return trace.channels[dye][trace.peaks[index] + offset];
}

/** Fit a common peak shape only where its measured control samples support it.
 * The ADC intervals propagate rounding on both the apex and neighboring samples.
 * A displaced or overlapping peak falls back to the integrated-area measurement. */
export function calibratedPeakAmplitudes(
 control: SangerTrace, sample: SangerTrace, sequences: string[], controlTrim: number, sampleTrim: number,
 training: number[], validation: number[], controlBackground: number[], sampleBackground: number[],
) {
 const offsets = [-4,-3,-2,-1,0,1,2,3,4];
 const inside = (trace: SangerTrace, i: number) => i > 0 && i + 1 < trace.peaks.length &&
  trace.peaks[i] - trace.peaks[i-1] > 8 && trace.peaks[i+1] - trace.peaks[i] > 8;
 const shapes = dyes.map((dye,b) => {
  const background = controlBackground[b];
  const fit = (point: number, bounds: Interval[]) => {
   const index = point - controlTrim;
   const center = traceAt(control,dye,index,0) - background;
   for (const [k,d] of offsets.entries()) {
    if (!d) continue;
    const value = traceAt(control,dye,index,d) - background;
    bounds[k][0] = Math.max(bounds[k][0], (value-.5)/(center+.5));
    bounds[k][1] = Math.min(bounds[k][1], (value+.5)/(center-.5));
   }
  };
  const clean = (p: number) => inside(control,p-controlTrim) && sequences.every(s=>s[p]===dye) &&
   traceAt(control,dye,p-controlTrim,0) > background+2;
  const train = training.filter(clean), test = validation.filter(clean);
  if (!train.length || !test.length) return undefined;
  const bounds = offsets.map(d=>d ? [0, 1] as Interval : [1, 1] as Interval);
  train.forEach(p=>fit(p,bounds));
  const checked = bounds.map(([lo,hi])=>[lo,hi] as Interval);
  test.forEach(p=>fit(p,checked));
  if (checked.some(([lo,hi])=>lo>hi)) return undefined;
  return bounds;
 });
 if (shapes.some(shape=>shape===undefined)) return undefined;
 return (position: number) => {
  const index = position - sampleTrim;
  if (!inside(sample,index)) return undefined;
  const intervals = dyes.map((dye,b)=>{
   const shape = shapes[b]!;
   let lower = 0, upper = Infinity;
   for (const [k,d] of offsets.entries()) {
    const value = traceAt(sample,dye,index,d) - sampleBackground[b];
    const [a,z] = shape[k];
    if (z>0) lower=Math.max(lower,(value-.5)/z);
    if (a>0) upper=Math.min(upper,(value+.5)/a);
   }
   return {lower,upper};
  });
  if (intervals.some(({lower,upper})=>upper<lower)) return undefined;
  // A pure peak has no colour ratio that can be resolved by the shape.
  // Keep its area-based observation in the regular representative fit.
  if (intervals.filter(({lower})=>lower>0).length<=1) return undefined;
  return intervals.map(({lower,upper})=>(lower+upper)/2);
 };
}
