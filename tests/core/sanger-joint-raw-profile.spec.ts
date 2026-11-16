import {it,expect,beforeAll} from 'vitest';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import {jointRawProfile} from '../probes/experiments/sanger-joint-raw-profile';
import {materialize,type PathEdit} from '../../src/app/workers/core/sanger-path-common';
let optimizer:any;beforeAll(async()=>{optimizer=await loadHighs();});
it.each([5,6])('jointly recovers %s classes from an asymmetric measured pulse',count=>{
 let seed=90219;const sequence=Array.from({length:160},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return 'ACGT'[seed>>>30];}).join(''),cut=42,event=(word:string,offset=0):PathEdit[]=>[{start:cut+offset,end:cut+offset,inserted:word}],words=['C','GT','AGC','TGCA','GGACT','CTAGGC'],truth=[[],...words.slice(0,count).map((w,i)=>event(w,i%3-1))],initial=[[],...words.slice(0,count).map(w=>event('A'.repeat(w.length)))],fractions=[1-Array.from({length:count},(_,i)=>(i+3)/100).reduce((a,b)=>a+b,0),...Array.from({length:count},(_,i)=>(i+3)/100)],sequences=truth.map(e=>materialize(sequence,e)!),positions=Array.from({length:120},(_,i)=>i+10),kernel=[.03,.08,.21,.62,1,.71,.36,.12,.04],background=[7,13,5,11],peaks=Array.from({length:160},(_,i)=>14+i*10),channels=Object.fromEntries([...'ACGT'].map((b,j)=>[b,new Uint16Array(1640).fill(background[j])])) as any;
 for(let p=0;p<150;p++)for(const[b,base]of[...'ACGT'].entries())for(let k=0;k<9;k++)channels[base][peaks[p]+k-4]+=Math.round(1000*kernel[k]*sequences.reduce((n,s,j)=>n+fractions[j]*Number(s[p]===base),0));
 const observed=positions.map(p=>[...'ACGT'].map(b=>sequences.reduce((n,s,j)=>n+fractions[j]*Number(s[p]===b),0))),result=jointRawProfile(sequence,cut,0,[[]],initial,positions,observed,[fractions],optimizer,60,.001,()=>{},{outside:[],sample:{channels,peaks},trim:0,background,scale:1000,kernel});
 expect(result?.valid).toBe(true);expect(result?.objective).toBeLessThan(1e-5);const recovered=result!.events!.map(e=>materialize(sequence,e)!);for(let i=0;i<sequences.length;i++){const found=recovered.indexOf(sequences[i]);expect(found).toBeGreaterThanOrEqual(0);expect(result!.weights![0][found]).toBeCloseTo(fractions[i],2);}
},90000);
