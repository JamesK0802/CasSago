import {it,expect} from 'vitest';
import {cycleFluorescence} from '../../src/app/workers/core/sanger-cycle-fluorescence';
it('preserves integrated photon mixtures before normalization, including uneven cycle widths',()=>{
 const make=(phase:number)=>({peaks:[8,19,34],channels:Object.fromEntries(['A','C','G','T'].map((base,b)=>[base,Array.from({length:44},(_,j)=>2+(b+1)*Math.exp(-((j-12-phase)**2)/18))]))} as any),a=make(0),b=make(4),mix={...a,channels:Object.fromEntries(['A','C','G','T'].map(base=>[base,a.channels[base].map((v:number,j:number)=>.3*v+.7*b.channels[base][j])]))},x=cycleFluorescence(a,[0,0,0,0]),y=cycleFluorescence(b,[0,0,0,0]),z=cycleFluorescence(mix,[0,0,0,0]);
 for(let p=0;p<3;p++)for(let dye=0;dye<4;dye++)expect(z.areas[p][dye]).toBeCloseTo(.3*x.areas[p][dye]+.7*y.areas[p][dye],10);
});
