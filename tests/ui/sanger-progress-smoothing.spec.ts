import {expect,it} from 'vitest';
import {sangerProgress,smoothSangerProgress} from '../../src/app/utils/sanger-progress';
it('does not drift above a reported stage even after ten minutes of timer ticks',()=>{const stage=sangerProgress('Comparing joint target hypotheses · reference 1 / 3 · 193 / 1050');let shown=0;for(let tick=0;tick<2400;tick++)shown=smoothSangerProgress(shown,stage.percent);expect(shown).toBeCloseTo(stage.percent,10);expect(shown).toBeLessThan(67);});
it('stays monotonic for late stage messages and reserves 100 for a result',()=>{let shown=65;for(let i=0;i<1000;i++){const next=smoothSangerProgress(shown,100);expect(next).toBeGreaterThanOrEqual(shown);expect(next).toBeLessThan(100);shown=next;}expect(smoothSangerProgress(shown,40)).toBe(shown);});
