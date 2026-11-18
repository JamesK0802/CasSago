import { describe, expect, it } from 'vitest';
import type { SangerTrace } from '../../src/app/models/sanger-analysis';
import { calibratedPeakAmplitudes } from '../../src/app/workers/core/sanger-peak-shape';

const bases = ['A', 'C', 'G', 'T'] as const;
const sequence = 'TACGTACGTACGTA';
const background = [2, 2, 2, 2];
const training = [1, 2, 3, 4, 5, 6, 7, 8];
const validation = [9, 10, 11, 12];
const peaks = Uint32Array.from(sequence, (_, i) => 10 + 10*i);
const shape = [.02, .1, .37, .78, 1, .78, .37, .1, .02];

function trace() {
 const channels = { A:new Float32Array(150), C:new Float32Array(150), G:new Float32Array(150), T:new Float32Array(150) };
 for (const dye of bases) channels[dye].fill(2);
 const data = { name: 'calibration trace', sampleName:'control', channels, peaks, sequence, quality: null } satisfies SangerTrace;
 for (let i=0; i<sequence.length; i++) {
  shape.forEach((factor,k) => { channels[sequence[i] as keyof typeof channels][peaks[i]+k-4]=Math.round(1100*factor)+2; });
 }
 return data;
}

function estimate(control: SangerTrace, sample: SangerTrace) {
 return calibratedPeakAmplitudes(control, sample, [sequence], 0, 0, training, validation, background, background);
}

describe('control-supported peak shapes', () => {
 it('estimates amplitudes for several dyes in a sample cycle from one calibrated shape', () => {
  const control = trace();
  const sample = trace();
  const index = 4;
  const expected = [120.25, 270.75, 126.5, 451];
  for (const [b,dye] of bases.entries()) {
   shape.forEach((factor,k) => { sample.channels[dye][peaks[index]+k-4]=Math.round(expected[b]*factor)+2; });
  }
  const measured = estimate(control,sample)!(index)!;
  measured.forEach((v,b)=>expect(Math.abs(v-expected[b])).toBeLessThan(.51));
  sample.channels.T[peaks[index]+1] += 100;
  expect(estimate(control,sample)!(index)).toBeUndefined();
 });

 it('keeps the area observation for a single-dye peak', () => {
  const control = trace();
  expect(estimate(control,trace())!(1)).toBeUndefined();
 });

 it('declines the shape when independently held-out control cycles disagree', () => {
  const control = trace();
  control.channels.A[peaks[9]+1] -= 100;
  expect(estimate(control,trace())).toBeUndefined();
 });
});
