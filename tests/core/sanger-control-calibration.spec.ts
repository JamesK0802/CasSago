import {expect,it} from 'vitest';
import {calibrateControlMixture} from '../../src/app/workers/core/sanger-control-calibration';
const bases='ACGT',means=[[.85,.08,.04,.03],[.02,.92,.04,.02],[.08,.02,.86,.04],[.01,.04,.10,.85]];
const sequences=Array.from({length:3},(_,h)=>Array.from({length:256},(_,i)=>bases[(i>>(2*h))%4]).join('')),positions=Array.from({length:256},(_,i)=>i);
for(const shares of [[1],[.3,.7],[.15,.35,.5]])it('recovers control-only reference shares with asymmetric dye mixing '+shares.join('/'),()=>{
 const refs=sequences.slice(0,shares.length),control=positions.map(p=>bases.split('').map((_,d)=>refs.reduce((n,s,h)=>n+shares[h]*means[bases.indexOf(s[p])][d],0))),fit=calibrateControlMixture(refs,positions,control);
 expect(fit.converged).toBe(true);fit.shares.forEach((v,h)=>expect(v).toBeCloseTo(shares[h],5));fit.means.forEach((row,b)=>row.forEach((v,d)=>expect(v).toBeCloseTo(means[b][d],5)));
});
it('retains an ideal pure-dye control without creating cross-colour signal',()=>{const control=positions.map(p=>bases.split('').map(b=>Number(sequences[0][p]===b))),fit=calibrateControlMixture([sequences[0]],positions,control);expect(fit.shares).toEqual([1]);expect(fit.means).toEqual(bases.split('').map((_,b)=>bases.split('').map((_,d)=>Number(b===d))));});
