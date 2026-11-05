import {fitEmissionMixture} from './sanger-emission';
import {mixPhase,phaseNeighbor} from './sanger-signal';
export interface ControlResponseProfile {
 table:number[][];endTable?:number[][];knots?:number[][][];start:number;end:number;
 sites?:Record<number,Record<number,number[]>>;
}
const bases='ACGT';
function context(sequence:string,p:number):number {
 const l=bases.indexOf(sequence[p-1]),b=bases.indexOf(sequence[p]),r=bases.indexOf(sequence[p+1]);
 return l<0||b<0||r<0?-1:16*l+4*b+r;
}
function weightsAt(p:number,start:number,end:number,count:number):[number,number,number] {
 const x=Math.max(0,Math.min(count-1,(p-start)/Math.max(1,end-start)*(count-1))),a=Math.floor(x),b=Math.min(count-1,a+1);
 return [a,b,x-a];
}
function responseAt(profile:ControlResponseProfile,id:number,p:number):number[] {
 const tables=profile.knots||(profile.endTable?[profile.table,profile.endTable]:[profile.table]);
 const [a,b,t]=weightsAt(p,profile.start,profile.end,tables.length);
 return tables[a][id].map((v,d)=>(1-t)*v+t*tables[b][id][d]);
}
const responseCache=new WeakMap<ControlResponseProfile,Map<number,number[][]>>();
export function controlKnownContextEmission(left:string,current:string,right:string,profile:ControlResponseProfile,p:number):number[] {
 if(!current)return [0,0,0,0];
 const a=left.length===1?bases.indexOf(left):-1,b=bases.indexOf(current),c=right.length===1?bases.indexOf(right):-1;
 if(a<0||b<0||c<0){const pi=(base:string)=>Array.from(bases,x=>base==='N'?.25:Number(base===x));return controlContextProbability(pi(left),pi(current),pi(right),profile,p);}
 const id=16*a+4*b+c,local=profile.sites?.[p]?.[id];if(local)return local;
 let cache=responseCache.get(profile);if(!cache){cache=new Map();responseCache.set(profile,cache);}let row=cache.get(p);if(!row){row=Array.from({length:64},(_,context)=>responseAt(profile,context,p));cache.set(p,row);}return row[id];
}
export function controlContextProbability(left:number[],current:number[],right:number[],profile:ControlResponseProfile,p:number):number[] {
 const out=[0,0,0,0],l=left.some(x=>x)?left:[.25,.25,.25,.25],r=right.some(x=>x)?right:[.25,.25,.25,.25];
 for(let a=0;a<4;a++)if(l[a])for(let b=0;b<4;b++)if(current[b])for(let c=0;c<4;c++)if(r[c]){
  const id=16*a+4*b+c,amount=l[a]*current[b]*r[c],v=profile.sites?.[p]?.[id]||responseAt(profile,id,p);
  for(let d=0;d<4;d++)out[d]+=amount*v[d];
 }
 return out;
}
export function controlContextEmission(sequence:string,p:number,lag:number,profile:ControlResponseProfile):number[] {
 const probabilities=(at:number)=>Array.from(bases,b=>sequence[at]==='N'?.25:Number(sequence[at]===b));
 const one=(at:number)=>context(sequence,at)<0&&profile.sites?.[at]?.[-1]?profile.sites[at][-1]:controlContextProbability(probabilities(at-1),probabilities(at),probabilities(at+1),profile,at);
 return lag===0?one(p):mixPhase(one(p),one(phaseNeighbor(p,lag)),lag);
}
/** Nonnegative ridge regression of absolute four-channel fluorescence. */
function solveResponse(X:number[][],targets:number[][],prior:number[][],penalty:number):number[][] {
 const n=prior.length,gram=Array.from({length:n},()=>new Float64Array(n)),rhs=Array.from({length:n},()=>[0,0,0,0]);
 X.forEach((x,j)=>{const active=x.flatMap((v,i)=>v?[[i,v]]:[]);for(const [a,v]of active){for(let d=0;d<4;d++)rhs[a][d]+=v*targets[j][d];for(const [b,w]of active)gram[a][b]+=v*w;}});
 for(let a=0;a<n;a++){gram[a][a]+=penalty;for(let d=0;d<4;d++)rhs[a][d]+=penalty*prior[a][d];}
 const neighbors=gram.map((row,a)=>Array.from(row,(v,b)=>({v,b})).filter(x=>x.b!==a&&x.v!==0)),B=prior.map(v=>[...v]);
 for(let iteration=0;iteration<500;iteration++){
  let delta=0;for(let a=0;a<n;a++)for(let d=0;d<4;d++){
   let value=rhs[a][d];for(const {v,b}of neighbors[a])value-=v*B[b][d];
   const next=Math.max(0,value/Math.max(1e-10,gram[a][a]));delta=Math.max(delta,Math.abs(next-B[a][d]));B[a][d]=next;
  }
  if(delta<1e-7)break;
 }
 return B;
}
/** One WT-only calibration for single and pooled traces. Context responses are
 * evaluated at the actual read cycle, never at the pre-edit donor coordinate.
 * A smooth scalar envelope separates read brightness from context yield. */
export function fitControlResponse(sequences:string[],positions:number[],observed:number[][],priorShares:number[]) {
 const rows=positions.filter(p=>observed[p]?.some(v=>v>0)&&sequences.some(s=>context(s,p)>=0)&&sequences.every(s=>!s[p]||context(s,p)>=0));
 if(rows.length<80)throw Error('Too few WT control cycles for response calibration.');
 const totals=rows.map(p=>observed[p].reduce((s,v)=>s+v,0)).sort((a,b)=>a-b),scale=totals[Math.floor(totals.length/2)];
 const Y=new Map(rows.map(p=>[p,observed[p].map(v=>v/scale)])),block=Math.max(2,Math.ceil(Math.cbrt(rows.length)));
 const train=rows.filter(p=>Math.floor(p/block)%5!==0),valid=rows.filter(p=>Math.floor(p/block)%5===0),start=rows[0],end=rows.at(-1)!,bandwidth=Math.max(30,(end-start)/6);
 const ids=new Map(rows.map(p=>[p,sequences.map(s=>context(s,p))]));
 const basePrior=Array.from({length:4},(_,b)=>Array.from(bases,(_,d)=>Number(b===d)));
 const profile=(table:number[][],capacity:number):ControlResponseProfile=>{
  const knots=capacity===4?[Array.from({length:64},(_,c)=>table[Math.floor(c/4)%4])]:Array.from({length:capacity/64},(_,i)=>table.slice(i*64,(i+1)*64));
  return {table:knots[0],endTable:knots.length>1?knots.at(-1):undefined,knots,start,end};
 };
 const fit=(capacity:number,penalty:number,priorTable:number[][],training=train)=>{
  let shares=[...priorShares],table=priorTable.map(v=>[...v]),envelope=new Map(rows.map(p=>[p,1]));
  const features=(p:number)=>{
   const out=Array(capacity).fill(0),[a,b,t]=weightsAt(p,start,end,capacity===4?1:capacity/64);
   ids.get(p)!.forEach((c,h)=>{if(c<0)return;if(capacity===4)out[Math.floor(c/4)%4]+=shares[h];else{out[64*a+c]+=shares[h]*(1-t);out[64*b+c]+=shares[h]*t;}});
   return out;
  };
  const predict=(p:number)=>{const x=features(p),v=[0,0,0,0];for(let c=0;c<capacity;c++)if(x[c])for(let d=0;d<4;d++)v[d]+=x[c]*table[c][d];return v;};
  for(let iteration=0;iteration<8;iteration++){
   table=solveResponse(training.map(features),training.map(p=>Y.get(p)!.map(v=>v/envelope.get(p)!)),priorTable,penalty);
   const ratio=training.map(p=>Math.log(Math.max(1e-8,Y.get(p)!.reduce((s,v)=>s+v,0))/Math.max(1e-8,predict(p).reduce((s,v)=>s+v,0))));
   for(const p of rows){let sum=0,mass=0;for(let j=0;j<training.length;j++){const w=Math.exp(-.5*((p-training[j])/bandwidth)**2);sum+=w*ratio[j];mass+=w;}envelope.set(p,Math.exp(sum/Math.max(1e-10,mass)));}
   const expanded=profile(table,capacity),cols=sequences.map(s=>training.flatMap(p=>controlContextEmission(s,p,0,expanded))),y=training.flatMap(p=>{const a=Y.get(p)!,sum=a.reduce((s,v)=>s+v,0);return a.map(v=>v/sum);});
   shares=fitEmissionMixture(cols,y).weights;
  }
  let raw=0,colour=0,marker=0,markers=0;
  for(const p of valid){
   const v=predict(p),sum=v.reduce((s,x)=>s+x,0),y=Y.get(p)!,total=y.reduce((s,x)=>s+x,0),informative=new Set(ids.get(p)!.filter(id=>id>=0)).size>1;
   for(let d=0;d<4;d++){raw+=(y[d]-v[d]*envelope.get(p)!)**2;const error=(y[d]/total-(sum?v[d]/sum:0))**2;colour+=error;if(informative)marker+=error;}
   if(informative)markers++;
  }
  return {...profile(table,capacity),shares,markerValidationMse:markers>=6?marker/(4*markers):undefined,markerValidationBases:markers,rawValidationMse:raw/(4*valid.length),validationMse:colour/(4*valid.length),penalty,capacity,envelope,priorTable};
 };
 const base=fit(4,.1,basePrior);let best=base;
 // More context/read detail is accepted only if it predicts withheld WT cycles.
 for(const capacity of [64,128,256]){
  const count=capacity/64,prior=Array.from({length:count},(_,i)=>Array.from({length:64},(_,c)=>responseAt(best,c,start+(end-start)*(count>1?i/(count-1):.5)))).flat();
  const candidates=[.1,1,10].map(lambda=>fit(capacity,lambda,prior)).filter(c=>c.validationMse<best.validationMse*.95&&c.rawValidationMse<best.rawValidationMse&&(best.markerValidationMse===undefined||c.markerValidationMse!==undefined&&c.markerValidationMse<=best.markerValidationMse));
  const score=(c:typeof base)=>c.validationMse/Math.max(1e-12,base.validationMse)+(c.markerValidationMse??c.validationMse)/Math.max(1e-12,base.markerValidationMse??base.validationMse);
  best=candidates.reduce((a,b)=>score(a)<score(b)?a:b,best);
 }
 const final=fit(best.capacity,best.penalty,best.priorTable,rows),sites:Record<number,Record<number,number[]>>={};
 // Reuse a measured WT context at this read cycle. For a mixed control the
 // minimum nonnegative correction respects the independently fitted WT shares.
 for(const p of rows){
  const contexts=new Map<number,number>();ids.get(p)!.forEach((id,h)=>{if(id>=0)contexts.set(id,(contexts.get(id)||0)+final.shares[h]);});
  const keys=[...contexts.keys()],weights=keys.map(k=>contexts.get(k)!),prior=keys.map(k=>responseAt(final,k,p)),target=Y.get(p)!.map(v=>v/final.envelope.get(p)!),local=prior.map(v=>[...v]);
  for(let b=0;b<4;b++){
   let active=keys.map((_,i)=>i);
   while(active.length){
    const mass=active.reduce((s,i)=>s+weights[i]**2,0),prediction=active.reduce((s,i)=>s+weights[i]*prior[i][b],0),delta=(target[b]-prediction)/Math.max(1e-12,mass),negative=active.filter(i=>prior[i][b]+delta*weights[i]<0);
    if(!negative.length){for(const i of active)local[i][b]=prior[i][b]+delta*weights[i];break;}
    for(const i of negative)local[i][b]=0;active=active.filter(i=>!negative.includes(i));
   }
  }
  sites[p]=Object.fromEntries(keys.map((k,i)=>[k,local[i]]));
 }
 // Ambiguous WT bases still have measured fluorescence; keep it on its own
 // unchanged WT cycle rather than replacing it with a fabricated base call.
 if(sequences.length===1)for(let p=0;p<observed.length;p++)if(context(sequences[0],p)<0&&observed[p]?.some(v=>v>0)){
  const nearest=Math.max(start,Math.min(end,p)),gain=final.envelope.get(nearest)||1;sites[p]={[-1]:observed[p].map(v=>v/(scale*gain))};
 }
 return {...best,envelope:undefined,priorTable:undefined,table:final.table,endTable:final.endTable,knots:final.knots,shares:final.shares,sites,used:best!==base,baseValidationMse:base.validationMse,trainingBases:train.length,validationBases:valid.length,refitBases:rows.length};
}
