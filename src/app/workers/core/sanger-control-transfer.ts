import {ControlResponseProfile,controlContextEmission} from './sanger-control-response';
import {fitEmissionMixture} from './sanger-emission';
import {correctGain} from './sanger-signal';
/** Relative dye gains require absolute WT peak areas: normalizing a pure WT
 * peak first removes nearly all information about its emitting dye's gain.
 * A smooth scalar lane envelope is nuisance. Validation is WT-only and blocked. */
export function fitControlTransfer(sequences:string[],positions:number[],areas:number[][],profile:ControlResponseProfile,lag:number,initialGains:number[],initialShares:number[]) {
 const rows=positions.filter(p=>areas[p]?.some(x=>x>0)),block=Math.max(2,Math.ceil(Math.cbrt(rows.length))),train=rows.filter(p=>Math.floor(p/block)%5!==0),valid=rows.filter(p=>Math.floor(p/block)%5===0);
 const totals=rows.map(p=>areas[p].reduce((s,x)=>s+x,0)).sort((a,b)=>a-b),scale=totals[Math.floor(totals.length/2)]||1,Y=new Map(rows.map(p=>[p,areas[p].map(v=>v/scale)])),E=new Map(rows.map(p=>[p,sequences.map(s=>controlContextEmission(s,p,lag,profile))])),bandwidth=Math.max(30,(rows.at(-1)!-rows[0])/6);
 const fit=(learn:boolean)=>{
  let gains=[...initialGains],shares=[...initialShares],envelope=new Map(rows.map(p=>[p,1]));
  const predict=(p:number)=>Array.from({length:4},(_,b)=>E.get(p)!.reduce((sum,v,h)=>sum+shares[h]*v[b]*gains[b],0));
  for(let iteration=0;iteration<12;iteration++){
   const ratio=train.map(p=>Math.log(Math.max(1e-10,Y.get(p)!.reduce((s,x)=>s+x,0))/Math.max(1e-10,predict(p).reduce((s,x)=>s+x,0))));
   for(const p of rows){let mass=0,sum=0;for(let j=0;j<train.length;j++){const w=Math.exp(-.5*((p-train[j])/bandwidth)**2);mass+=w;sum+=w*ratio[j];}envelope.set(p,Math.exp(sum/Math.max(1e-10,mass)));}
   if(learn){
    const next=Array.from({length:4},(_,b)=>{let xy=0,xx=0;for(const p of train){const x=envelope.get(p)!*E.get(p)!.reduce((s,v,h)=>s+shares[h]*v[b],0);xy+=x*Y.get(p)![b];xx+=x*x;}return Math.max(1e-10,xy/Math.max(1e-10,xx));}),average=next.reduce((s,x)=>s+Math.log(x),0)/4;
    gains=next.map(x=>Math.exp(Math.log(x)-average));
   }
   shares=fitEmissionMixture(sequences.map((_,h)=>train.flatMap(p=>E.get(p)![h])),train.flatMap(p=>correctGain(Y.get(p)!,gains))).weights;
  }
  let raw=0,colour=0,variance=0;const normalized=valid.flatMap(p=>{const v=Y.get(p)!,s=v.reduce((a,b)=>a+b,0);return v.map(x=>x/s);}),mean=normalized.reduce((s,x)=>s+x,0)/Math.max(1,normalized.length);
  for(const p of valid){const v=predict(p),sum=v.reduce((s,x)=>s+x,0),y=Y.get(p)!,total=y.reduce((s,x)=>s+x,0);for(let b=0;b<4;b++){raw+=(y[b]-envelope.get(p)!*v[b])**2;colour+=(y[b]/total-(sum?v[b]/sum:0))**2;variance+=(y[b]/total-mean)**2;}}
  return {gains,shares,rawMse:raw/Math.max(1,4*valid.length),colourMse:colour/Math.max(1,4*valid.length),colourR2:variance?1-colour/variance:0};
 };
 const baseline=fit(false),candidate=train.length>=40&&valid.length>=20?fit(true):baseline,used=candidate!==baseline&&candidate.gains.every(x=>x>=.55&&x<=1.8)&&candidate.rawMse<baseline.rawMse*.95&&candidate.colourMse<=baseline.colourMse+1e-10;
 return {...(used?candidate:baseline),used,initialRawMse:baseline.rawMse,initialColourMse:baseline.colourMse,trainingBases:train.length,validationBases:valid.length};
}
