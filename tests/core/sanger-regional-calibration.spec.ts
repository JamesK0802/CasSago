import {it,expect} from 'vitest';
import {calibrateRegionalResponse} from '../../src/app/workers/core/sanger-regional-calibration';
it('fits one continuous reference drift to prefix and flanks while preserving within-reference phases',()=>{
 const positions=Array.from({length:120},(_,i)=>220+i),prefixPositions=Array.from({length:60},(_,i)=>30+i),shares=[.5,.5],weights=[.35,.15,.1,.4],emit=(h:number,size:number,p:number)=>{const row=Array.from({length:4},(_,b)=>.1+((p*17+h*31+size*13+b*29+p*h*7)%61)/61),total=row.reduce((a,b)=>a+b,0);return row.map(v=>v/total);},gain=(h:number,p:number)=>h?1:Math.exp(.08+.4*(p-59.5)/(339-59.5)),classes=[0,1,2,3].map(j=>({reference:j>>1,size:j%2,column:positions.flatMap(p=>emit(j>>1,j%2,p))})),observed=positions.flatMap(p=>{const v=[0,1,2,3].map(b=>weights.reduce((s,w,j)=>s+w*emit(j>>1,j%2,p)[b]*gain(j>>1,p),0)),total=v.reduce((a,b)=>a+b,0);return v.map(v=>v/total);}),prefix={positions:prefixPositions,columns:[0,1].map(h=>prefixPositions.flatMap(p=>emit(h,0,p))),observed:prefixPositions.flatMap(p=>{const v=[0,1,2,3].map(b=>shares.reduce((s,w,h)=>s+w*emit(h,0,p)[b]*gain(h,p),0)),total=v.reduce((a,b)=>a+b,0);return v.map(v=>v/total);})};
 const fit=calibrateRegionalResponse([{classes,observed}],[positions],shares,.02,prefix);
 expect(fit.at(0,60)).toBeCloseTo(gain(0,60),3);expect(fit.at(0,300)).toBeCloseTo(gain(0,300),3);
 expect(fit.lengths[0][0][0].fraction).toBeCloseTo(.35,3);expect(fit.lengths[0][1][1].fraction).toBeCloseTo(.4,3);
 const changed=observed.map((v,k)=>Math.floor(positions[k>>2]/5)%3===2?1-v:v),again=calibrateRegionalResponse([{classes,observed:changed}],[positions],shares,.02,prefix);
 expect(again.parameters).toEqual(fit.parameters);
},120000);
