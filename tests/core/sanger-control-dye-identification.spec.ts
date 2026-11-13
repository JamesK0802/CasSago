import {expect,it} from 'vitest';
import {calibrateControlMixture} from '../../src/app/workers/core/sanger-control-calibration';
const bases='ACGT',means=[[.85,.08,.04,.03],[.02,.92,.04,.02],[.08,.02,.86,.04],[.01,.04,.10,.85]];
const sequences=Array.from({length:3},(_,h)=>Array.from({length:256},(_,i)=>bases[(i>>(2*h))%4]).join('')),positions=Array.from({length:256},(_,i)=>i);
for(const shares of [[1],[.3,.7],[.15,.35,.5]])it('recovers control-only reference shares with asymmetric dye mixing '+shares.join('/'),()=>{
 const refs=sequences.slice(0,shares.length),control=positions.map(p=>bases.split('').map((_,d)=>refs.reduce((n,s,h)=>n+shares[h]*means[bases.indexOf(s[p])][d],0))),fit=calibrateControlMixture(refs,positions,control);
 expect(fit.converged).toBe(true);fit.shares.forEach((v,h)=>expect(v).toBeCloseTo(shares[h],5));fit.means.forEach((row,b)=>row.forEach((v,d)=>expect(v).toBeCloseTo(means[b][d],5)));
});
it('retains an ideal pure-dye control without creating cross-colour signal',()=>{const control=positions.map(p=>bases.split('').map(b=>Number(sequences[0][p]===b))),fit=calibrateControlMixture([sequences[0]],positions,control);expect(fit.shares).toEqual([1]);expect(fit.means).toEqual(bases.split('').map((_,b)=>bases.split('').map((_,d)=>Number(b===d))));});

it('uses WT monochromatic cycles to identify dye spectra independently of reference-specific base yields',()=>{
 const refs=sequences.slice(0,2),abundance=[.5,.5],yieldByReference=[[.8,1.1,2,.9],[1.4,.7,.5,1.2]],control=positions.map(p=>{const row=[0,1,2,3].map(b=>refs.reduce((n,s,h)=>n+abundance[h]*yieldByReference[h][b]*Number(s[p]===bases[b]),0)),total=row.reduce((a,b)=>a+b,0);return row.map(v=>v/total);}),fit=calibrateControlMixture(refs,positions,control);
 expect(fit.converged).toBe(true);for(let b=0;b<4;b++)for(let d=0;d<4;d++)expect(fit.means[b][d]).toBeCloseTo(Number(b===d),10);
});

it('does not treat a missing pure-control measurement as a zero dye response',()=>{const fit=calibrateControlMixture(['AACCGGTT'],Array.from({length:8},(_,i)=>i),[[0,0,0,0],[1,0,0,0],[0,1,0,0],[0,1,0,0],[0,0,1,0],[0,0,1,0],[0,0,0,1],[0,0,0,1]]);expect(fit.means).toEqual([0,1,2,3].map(b=>[0,1,2,3].map(d=>Number(b===d))));});
