import {it,expect} from 'vitest';
import {conditionalProfileRows} from '../probes/conditional-profile-rows';
const check=(mass:number,ratio:number)=>{const rows=conditionalProfileRows([[[1,1]]],[[0,1]],[.6]),x=[mass,mass*ratio];return rows.every(r=>{const v=r.terms.reduce((n,[i,w])=>n+w*x[i],0);return v>=r.lower-1e-12&&v<=r.upper+1e-12;});};
it.each([.01,.1,.5,1])('enforces exactly two conditional percentage points at candidate mass %s',mass=>{expect(check(mass,.619)).toBe(true);expect(check(mass,.621)).toBe(false);expect(check(mass,.581)).toBe(true);expect(check(mass,.579)).toBe(false);});
it('rejects an error a larger reporting-point mass would wrongly permit',()=>{const candidateMass=.1,pointMass=.5,ratio=.66;expect(Math.abs(candidateMass*(ratio-.6))).toBeLessThan(.02*pointMass);expect(check(candidateMass,ratio)).toBe(false);});
it('accepts a valid error a smaller reporting-point mass would wrongly reject',()=>{const candidateMass=.5,pointMass=.1,ratio=.615;expect(Math.abs(candidateMass*(ratio-.6))).toBeGreaterThan(.02*pointMass);expect(check(candidateMass,ratio)).toBe(true);});
