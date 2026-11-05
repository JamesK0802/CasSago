import {it,expect} from 'vitest';
import {CandidateSignalScreen} from '../../src/app/workers/core/sanger-candidate-screen';
it('preserves exact raw-fluorescence screening scores and replacement gains as candidates are appended',()=>{
 const columns=[[1,0,0,0,0,1,0,0],[.8,.2,0,0,0,1,0,0],[0,0,1,0,0,0,0,1]],families=['A','A','B'],residual=[.1,-.1,.03,.01,-.04,.08,.01,-.1],screen=new CandidateSignalScreen(columns);screen.update(families);
 columns.push([.7,.1,.2,0,0,.9,.1,0]);families.push('A');screen.update(families);const correlation=screen.correlate(residual);for(let id=0;id<columns.length;id++)expect(correlation[id]).toBeCloseTo(columns[id].reduce((s,v,k)=>s+v*residual[k],0),12);
 const weight=.15,difference=columns[3].map((v,k)=>v-columns[1][k]),expected=difference.reduce((s,v,k)=>s+2*weight*v*residual[k]-weight*weight*v*v,0);expect(screen.replacementGain(3,1,weight)).toBeCloseTo(expected,12);expect(screen.sameFamily(0,2)).toBe(false);expect(()=>screen.replacementGain(0,2,weight)).toThrow();
});

it('ranks a localized substitution above long paths with the same residual correlation',()=>{const wt=Array.from({length:30},()=>[0,0,1,0]).flat(),snp=[...wt];snp.splice(40,4,1,0,0,0);const bad=Array.from({length:25},(_,i)=>{const col=[...snp];col.splice(4*i,4,0,1,0,0);return col;}),columns=[wt,...bad,snp],y=wt.map((v,k)=>.9134*v+.0866*snp[k]),screen=new CandidateSignalScreen(columns);screen.update(columns.map((_,i)=>'f'+i));const scores=screen.additionScores(y,wt);expect(scores.at(-1)).toBeGreaterThan(Math.max(...scores.slice(1,-1)));});
