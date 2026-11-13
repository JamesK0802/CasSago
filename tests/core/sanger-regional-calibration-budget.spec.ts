import {it,expect} from 'vitest';
import {readFileSync} from 'node:fs';
import {calibrateRegionalResponse} from '../../src/app/workers/core/sanger-regional-calibration';
it('finishes a valid detector calibration that converges just beyond the old80-step budget',()=>{
 const input=JSON.parse(readFileSync(new URL('../fixtures/sanger-regional-calibration-budget.json',import.meta.url),'utf8'));
 const r=calibrateRegionalResponse(input.windows,input.positions,input.shares,input.delta,input.prefix);
 expect(r.iterations).toBeGreaterThanOrEqual(80);expect(r.iterations).toBeLessThan(256);
 const independentlyExtended=[-.018160536886704846,.026232725643114515,-.022235301757648332,.2569448831912775];
 r.parameters.forEach((v,i)=>expect(v).toBeCloseTo(independentlyExtended[i],5));
 for(let h=0;h<input.shares.length;h++){const total=r.lengths[0][h].reduce((n,g)=>n+g.fraction,0);expect(total).toBeCloseTo(input.shares[h],7);}
 expect(r.trainingPredictions[0].every(v=>Number.isFinite(v)&&v>=0&&v<=1)).toBe(true);
},30000);
