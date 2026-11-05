import {SimplexFit,phaseNeighbor} from './sanger-signal';
export interface LatentInsertionFit {fit:SimplexFit;columns:number[][];unknown:boolean[];knownTrainingRows:number;probabilities:{position:number;distribution:number[]}[]}
/** Infer mixture weights from known flanks, never from arbitrary N=25% peaks.
 * Coincident unknown inserted bases receive a pooled probability distribution;
 * this does not assign a unique inserted haplotype to each molecule. */
export function fitLatentInsertions(columns:number[][],sequences:string[],eligible:boolean[],positions:number[],observed:number[][],training:number[],lag:number,fitter:(c:number[][],y:number[])=>SimplexFit):LatentInsertionFit{
 const n=positions.length,unknown=positions.map(p=>sequences.some((s,h)=>eligible[h]&&(s[p]==='N'||s[p-1]==='N'||s[p+1]==='N'))),known=training.filter(k=>!unknown[k]);
 if(!unknown.some(Boolean)||known.length<Math.max(100,columns.length*5)||unknown.filter(Boolean).length>n*.4){return {fit:fitter(columns.map(c=>training.flatMap(k=>c.slice(k*4,k*4+4))),training.flatMap(k=>observed[k])),columns,unknown:positions.map(()=>false),knownTrainingRows:training.length,probabilities:[]};}
 const fit=fitter(columns.map(c=>known.flatMap(k=>c.slice(k*4,k*4+4))),known.flatMap(k=>observed[k])),effective=columns.map(c=>[...c]),probabilities:{position:number;distribution:number[]}[]=[],byPosition=new Map<number,number[]>();
 const project=(v:number[],mass:number)=>{const sorted=[...v].sort((a,b)=>b-a);let sum=0,threshold=0;for(let k=0;k<4;k++){sum+=sorted[k];threshold=(sum-mass)/(k+1);if(k===3||sorted[k+1]<=threshold)break;}return v.map(x=>Math.max(0,x-threshold));};
 const order=Array.from({length:n},(_,k)=>k);if(lag<0)order.reverse();
 for(const k of order){
  const p=positions[k],knownSignal=[0,0,0,0];let missing=0;
  for(let h=0;h<columns.length;h++){
   const current=eligible[h]&&sequences[h][p]==='N',previous=eligible[h]&&sequences[h][phaseNeighbor(p,lag)]==='N',prior=byPosition.get(phaseNeighbor(p,lag))||[.25,.25,.25,.25];
   for(let b=0;b<4;b++){effective[h][4*k+b]=columns[h][4*k+b]-(current?(1-Math.abs(lag))/4:0)+(previous?Math.abs(lag)*prior[b]:0);knownSignal[b]+=fit.weights[h]*effective[h][4*k+b];}if(current)missing+=fit.weights[h]*(1-Math.abs(lag));
  }
  if(missing<=1e-9)continue;const total=knownSignal.reduce((s,v)=>s+v,0)+missing,distribution=project(observed[k].map((y,b)=>y*total-knownSignal[b]),missing).map(v=>v/missing);byPosition.set(p,distribution);probabilities.push({position:p,distribution});
  for(let h=0;h<columns.length;h++)if(eligible[h]&&sequences[h][p]==='N')for(let b=0;b<4;b++)effective[h][4*k+b]+=(1-Math.abs(lag))*distribution[b];
 }
 return {fit,columns:effective,unknown,knownTrainingRows:known.length,probabilities};
}
