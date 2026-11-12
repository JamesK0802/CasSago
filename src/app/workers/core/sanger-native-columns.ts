import type {Prepared,PathEdit} from './sanger-path-common';
import {normalizeMoleculeEdits} from './sanger-edit-normalize';
import {materialize} from './sanger-path-common';
import type {SangerTrace} from '../../models/sanger-analysis';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
interface Feature{control:boolean;p:number;b:number;other?:number;lo:number;hi:number}
interface Geometry{start:number;removed:number;inserted:number;delta:number;word?:string}
/** Full-molecule column generation. No fixed allele count; pricing enumerates
 * legal event geometry, and chooses each inserted base from the LP dual.
 * Rows stay fixed and new columns extend the existing simplex basis. */
export function nativeColumns(refs:Prepared[],control:SangerTrace,sample:SangerTrace,optimizer:SangerLinearOptimizer,progress:(s:string)=>void,budget=400,_displayThresholdPercent=2,audit?:(stage:string,milliseconds:number,detail?:Record<string,number>)=>void,seeds:{reference:number;events:PathEdit[];fraction:number}[]=[]){
 const setupStarted=performance.now();
 const N=sample.peaks.length;if(!optimizer.createModel||N!==control.peaks.length||refs.some(r=>r.reverse||r.cuts.length!==2))return null;
 const alphabet='ACGT',features:Feature[]=[];
 for(const [trace,ct]of [[control,true],[sample,false]]as const){let previous=[0,0,0,0];
  if(trace.peaks[0]!==14||Object.values(trace.channels).some(c=>c.length!==N*10+30))return null;
  for(let p=0;p<N;p++){
   if(trace.peaks[p]!==14+10*p)return null;
   const heights=[...alphabet].map((b,k)=>trace.channels[b as'A'][trace.peaks[p]-5]-2-previous[k]);if(heights.some(v=>!Number.isInteger(v)||v<0||v>3))return null;previous=heights;
   for(let b=0;b<4;b++){let lo=0,hi=1;for(let d=-4;d<=4;d++){const v=trace.channels[alphabet[b]as'A'][trace.peaks[p]+d]-2,a=1200*Math.exp(-d*d/4);lo=Math.max(lo,(v-.5+1e-4)/a);hi=Math.min(hi,(v+.5-1e-4)/a);}const a=1200*Math.exp(-25/4);lo=Math.max(lo,(heights[b]-.5+1e-4)/a);hi=Math.min(hi,(heights[b]+.5-1e-4)/a);if(lo>hi+1e-10)return null;features.push({control:ct,p,b,lo,hi});}
   const localBounds=features.slice(-4),called=alphabet.indexOf(trace.sequence[p]);if(called<0)return null;
   // At a zero-slack feasible solution these call rows are already implied by
   // the four measured ADC intervals. Keep only genuinely ambiguous calls.
   for(let b=0;b<4;b++)if(b!==called&&localBounds[called].lo<localBounds[b].hi+(b<called?1e-8:0))features.push({control:ct,p,b:called,other:b,lo:b<called?1e-8:0,hi:1e30});
  }
  if(previous.some((v,b)=>trace.channels[alphabet[b]as'A'][trace.peaks[N-1]+5]-2!==v))return null;
 }
 const codes=(seq:string)=>Uint8Array.from({length:N},(_,p)=>alphabet.indexOf(seq[p]));
 const referenceCodes=refs.map(r=>codes(r.sequence));
 // The display threshold never removes inference hypotheses or measured mass.
 let models:{h:number;events:PathEdit[];sequence:string;codes:Uint8Array}[]=[];const seen=new Set<string>();
 const add=(h:number,events:PathEdit[])=>{events=normalizeMoleculeEdits(refs[h].sequence,events);const sequence=materialize(refs[h].sequence,events);if(!sequence)return false;const key=h+':'+sequence;if(seen.has(key))return false;seen.add(key);models.push({h,events,sequence,codes:codes(sequence)});return true;};refs.forEach((_,h)=>add(h,[]));for(const row of seeds)if(row.reference>=0&&row.reference<refs.length)add(row.reference,row.events);
 // Preserve locally supported words when the global search changes basis.
 // Geometry and novel words remain available; this is initialization, not a cap.
 const geometries=refs.map((ref,h)=>ref.cuts.map((cut,t)=>{const out:Geometry[]=[{start:cut,removed:0,inserted:0,delta:0}],seenGeometry=new Set<string>();
  for(const row of seeds)if(row.reference===h){const local=row.events.filter(e=>ref.cuts.map(c=>Math.min(Math.abs(e.start-c),Math.abs(e.end-c))).reduce((best,d,i,ds)=>d<ds[best]?i:best,0)===t);if(local.length===1){const e=local[0],g={start:e.start,removed:e.end-e.start,inserted:e.inserted.length,delta:e.inserted.length-e.end+e.start,word:e.inserted},key=JSON.stringify(g);if(!seenGeometry.has(key)){seenGeometry.add(key);out.push(g);}}}
  for(let start=cut-3;start<=cut+3;start++)for(let delta=-10;delta<=10;delta++){const removed=delta<0?-delta:delta===0?10:0,inserted=delta>0?delta:delta===0?10:0;if(start>=0&&start+removed<=ref.sequence.length)out.push({start,removed,inserted,delta});}return out;}));
 const model=optimizer.createModel();if(!model.addCols||!model.deleteCols){model.dispose();return null;}model.options.set({output_flag:false,solver:'simplex',simplex_strategy:4,primal_feasibility_tolerance:1e-10,dual_feasibility_tolerance:1e-10,time_limit:300});
 const coefficient=(m:typeof models[number],f:Feature)=>{const b=(f.control?referenceCodes[m.h]:m.codes)[f.p];return Number(b===f.b)-(f.other===undefined?0:Number(b===f.other));};
 const lower:number[]=[1],upper:number[]=[1],dualFeatures:(Feature|null)[]=[null],slackSigns:number[]=[0];
 for(const f of features){lower.push(f.lo);upper.push(1e30);dualFeatures.push(f);slackSigns.push(1);if(f.hi<1e29){lower.push(-1e30);upper.push(f.hi);dualFeatures.push(f);slackSigns.push(-1);}}
 const numRows=lower.length;
 model.passModel({numCols:1,numRows,colCost:[1],colLower:[0],colUpper:[2],rowLower:lower,rowUpper:upper,matrix:{format:'csc',numCols:1,numRows,starts:[0,numRows-1],indices:Array.from({length:numRows-1},(_,i)=>i+1),values:slackSigns.slice(1)}});
 let loaded=0;
 audit?.('column-setup',performance.now()-setupStarted,{features:features.length});
 try{for(let iteration=0;iteration<budget;iteration++){
  const matrixStarted=performance.now(),n=models.length,starts=[0],indices:number[]=[],values:number[]=[];
  for(let j=loaded;j<n;j++){indices.push(0);values.push(1);for(let r=1;r<numRows;r++){const v=coefficient(models[j],dualFeatures[r]!);if(v){indices.push(r);values.push(v);}}starts.push(indices.length);}
  if(n>loaded)model.addCols({cost:Array(n-loaded).fill(0),lower:Array(n-loaded).fill(0),upper:Array(n-loaded).fill(1),matrix:{format:'csc',numCols:n-loaded,numRows,starts,indices,values}});loaded=n;
  audit?.('column-matrix',performance.now()-matrixStarted,{iteration,columns:n,rows:numRows});
  const solveStarted=performance.now();const status=model.run();audit?.('column-lp',performance.now()-solveStarted,{iteration,columns:n});if(status.modelStatus!==7){progress('LP status '+status.modelStatus+' after '+iteration+' iterations');return null;}const solution=model.getSolution()as any,weights=Array.from({length:n},(_,i)=>Math.max(0,solution.colValue[i+1])),violation=solution.colValue[0];
  progress('Full-signal dictionary · '+iteration+' iterations · '+n+' molecules · interval error '+violation);
  if(violation<1e-7){
   const validationStarted=performance.now();
   const valid=(w:number[])=>{for(const [trace,ct]of [[control,true],[sample,false]]as const){const channels=Array.from({length:4},()=>new Uint16Array(N*10+30).fill(2));for(let p=0;p<N;p++){const mass=[0,0,0,0];models.forEach((m,j)=>{const b=(ct?referenceCodes[m.h]:m.codes)[p];if(b<4)mass[b]+=w[j];});const called=mass.indexOf(Math.max(...mass));if(alphabet[called]!==trace.sequence[p]){progress('Call mismatch '+ct+' '+p+' '+mass.join(','));return false;}for(let b=0;b<4;b++)for(let d=-5;d<=5;d++)channels[b][14+p*10+d]+=Math.round(1200*mass[b]*Math.exp(-d*d/4));}for(let b=0;b<4;b++){const mismatch=channels[b].findIndex((v,i)=>v!==trace.channels[alphabet[b]as'A'][i]);if(mismatch>=0){progress('ADC mismatch '+ct+' '+b+' '+mismatch+' '+channels[b][mismatch]+' vs '+trace.channels[alphabet[b]as'A'][mismatch]);return false;}}}return true;};
   const mass=weights.reduce((a,b)=>a+b,0);let normalized=weights.map(w=>w/mass);if(valid(normalized)){audit?.('column-validation',performance.now()-validationStarted);return{rows:models.flatMap((m,i)=>normalized[i]>0?[{reference:m.h,events:m.events,fraction:normalized[i]}]:[]),iterations:iteration,candidates:n};}
   for(const digits of [14,12,10,8]){const scale=10**digits,rounded=normalized.map(w=>Math.round(w*scale)/scale),total=rounded.reduce((a,b)=>a+b,0),w=rounded.map(v=>v/total);if(valid(w)){audit?.('column-validation',performance.now()-validationStarted);return{rows:models.flatMap((m,i)=>w[i]>0?[{reference:m.h,events:m.events,fraction:w[i]}]:[]),iterations:iteration,candidates:n};}}
   progress('Raw reconstruction failed independent validation');return null;
  }
  const pricingStarted=performance.now();
  const sampleCosts=new Float64Array(N*4),controlCosts=new Float64Array(N*4);for(let r=1;r<dualFeatures.length;r++){const f=dualFeatures[r]!;const cost=f.control?controlCosts:sampleCosts,value=-solution.rowDual[r];cost[4*f.p+f.b]+=value;if(f.other!==undefined)cost[4*f.p+f.other]-=value;}
  const lambda=solution.rowDual[0],minimumPrefix=new Float64Array(N+1),bestBases=new Uint8Array(N);for(let p=0;p<N;p++){let b=0;for(let k=1;k<4;k++)if(sampleCosts[4*p+k]<sampleCosts[4*p+b])b=k;bestBases[p]=b;minimumPrefix[p+1]=minimumPrefix[p]+sampleCosts[4*p+b];}
  const basis=model.getBasis()as {colStatus:ArrayLike<number>},removed:number[]=[];models=models.filter((m,j)=>{if(!m.events.length||weights[j]>1e-10||basis.colStatus[j+1]===1)return true;removed.push(j+1);seen.delete(m.h+':'+m.sequence);return false;});if(removed.length)model.deleteCols({kind:'set',indices:removed});loaded=models.length;
  let added=0,firstGeometries=0,admittedFirstGeometries=0,geometryPairs=0;
  for(let h=0;h<refs.length;h++){
   const ref=refs[h],prefix=new Map<number,Float64Array>();for(let shift=-20;shift<=20;shift++){const a=new Float64Array(N+1);for(let p=0;p<N;p++){const b=alphabet.indexOf(ref.sequence[p-shift]);a[p+1]=a[p]+(b>=0?sampleCosts[4*p+b]:0);}prefix.set(shift,a);}
   let controlCost=0;for(let p=0;p<N;p++)controlCost+=controlCosts[4*p+referenceCodes[h][p]];
   firstGeometries+=geometries[h][0].length;admittedFirstGeometries+=geometries[h][0].length;
   const wordCost=(g:Geometry,shift:number)=>g.word===undefined?minimumPrefix[Math.max(0,Math.min(N,g.start+shift+g.inserted))]-minimumPrefix[Math.max(0,Math.min(N,g.start+shift))]:[...g.word].reduce((sum,b,i)=>{const p=g.start+shift+i,k=alphabet.indexOf(b);return sum+(p>=0&&p<N&&k>=0?sampleCosts[4*p+k]:0);},0);
   const clamp=(v:number)=>Math.max(0,Math.min(N,v)),byLength=new Map<string,{cost:number;a:Geometry;b:Geometry}>();
   for(const a of geometries[h][0])for(const b of geometries[h][1]){
    geometryPairs++;
    if(a.start+a.removed>b.start)continue;
    const i=clamp(a.start),j=clamp(a.start+a.inserted),k=clamp(b.start+a.delta),l=clamp(b.start+a.delta+b.inserted),p=prefix.get(a.delta)!,q=prefix.get(a.delta+b.delta)!;
    const cost=controlCost+prefix.get(0)![i]+wordCost(a,0)+p[k]-p[j]+wordCost(b,a.delta)+q[N]-q[l];
    const key=a.delta+'/'+b.delta;if(cost>=lambda-1e-8||byLength.has(key)&&cost>=byLength.get(key)!.cost)continue;byLength.set(key,{cost,a,b});
   }
   for(const {a,b}of [...byLength.values()].sort((a,b)=>a.cost-b.cost).slice(0,128)){const word=(start:number,length:number,first=false)=>Array.from({length},(_,i)=>alphabet[bestBases[start+i]]||'A').join(''),events:PathEdit[]=[];if(a.removed||a.inserted)events.push({start:a.start,end:a.start+a.removed,inserted:a.word??word(a.start,a.inserted,true)});if(b.removed||b.inserted)events.push({start:b.start,end:b.start+b.removed,inserted:b.word??word(b.start+a.delta,b.inserted)});if(add(h,events))added++;}
  }
  audit?.('column-pricing',performance.now()-pricingStarted,{iteration,added,firstGeometries,admittedFirstGeometries,geometryPairs});
  if(!added){progress('No pricing improvement; error '+violation);return null;}
 }}finally{model.dispose();}
 return null;
}
