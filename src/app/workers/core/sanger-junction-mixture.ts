import {SangerAnnotation} from '../../models/sanger-analysis';
import {fitJunctionPaths,PathEvent,pathAnnotations} from './sanger-junction-paths';
import {leftEquivalentDeletionStart} from '../../utils/sanger-placement';
import {controlContextProbability,controlKnownContextEmission,ControlResponseProfile} from './sanger-control-response';
import {reverseComplement} from '../../utils/biology.utils';
import {fitOrderedOnsets} from './sanger-onset-fit';
import {transitionProbabilityEmission} from './sanger-transition-response';
import {SimplexFit,mixPhase,phaseNeighbor} from './sanger-signal';
import {fitEmissionMixture,emissionVector} from './sanger-emission';

export interface JunctionReference {contextResponse?:ControlResponseProfile;sequence:string;cut:number;cuts?:number[];donors?:number[][];transitionResponse?:number[][];baseYields?:number[]}
export interface JunctionModel {reference:number;size:number;sizes?:number[];path?:PathEvent[];start:number;inserted:string;probabilities?:number[][];origins?:{start:number;weight:number;probabilities?:number[][]}[]}
export interface JunctionFit extends SimplexFit {searchComplete?:boolean;candidateCount?:number;objectiveHistory?:number[];models:JunctionModel[];columns:number[][];iterations:number;duplicateLengths:number;objective:number;predictColumns:(positions:number[])=>number[][]}
export interface JunctionRead {positions:number[];rowIds:number[];references:JunctionReference[];mappings:{reverse:boolean;offset:number}[];lag:number;theta:number[]}
export interface JunctionInput {reads?:JunctionRead[];classes?:{reference:number;size:number;sizes?:number[]}[];references:JunctionReference[];positions:number[];observed:number[];maxDeletion:number;maxInsertion:number;radius:number;lag:number;theta:number[];positionWeights?:number[];fit:(columns:number[][],observed:number[],weights?:number[])=>SimplexFit}
const bases='ACGT';
/** A physical deletion is [start,start+length), on either side of the cut.
 * The radius bounds distance of that interval from the cut, not its direction. */
export function junctionStarts(cut:number,size:number,radius:number,length:number):number[]{
 const lo=Math.max(0,cut-Math.max(0,-size)-radius),hi=Math.min(cut+radius,length-Math.max(0,-size));
 return Array.from({length:Math.max(0,hi-lo+1)},(_,i)=>lo+i);
}
/** Count one candidate per distinct deleted DNA, retaining an allowed origin. */
export function moleculeJunctionStarts(sequence:string,cut:number,size:number,radius:number):number[]{
 const starts=junctionStarts(cut,size,radius,sequence.length);if(size>=0)return starts;
 const unique=new Map<number,number>();for(const start of starts){const canonical=leftEquivalentDeletionStart(sequence,start,-size);if(!unique.has(canonical))unique.set(canonical,start);}return [...unique.values()];
}
export function junctionBase(reference:string,model:JunctionModel,p:number):string {
 if(p<model.start)return reference[p]||'';
 if(model.size<0)return reference[p-model.size]||'';
 return p<model.start+model.size?model.inserted[p-model.start]||'N':reference[p-model.size]||'';
}
/** Exact DNA aliases over the measured span, including neighbouring contexts.
 * Compare physical origins as well as fitted composite columns: different
 * origin mixtures must not hide an identical-DNA alternative length in repeats.
 * Including gaps between measured cycles makes this check conservative. */
export function aliasedDeletionMass(references:JunctionReference[],models:JunctionModel[],weights:number[],positions:number[],maxDeletion:number,radius:number):number {
 const lo=Math.max(0,Math.min(...positions)-1),hi=Math.max(...positions)+2;
 let mass=0;
 references.forEach((ref,h)=>{
  const active=models.flatMap((m,i)=>m.reference===h&&!m.path&&m.size<0&&weights[i]>.01?[i]:[]);if(!active.length)return;
  const signature=(start:number,length:number)=>(ref.sequence.slice(0,start)+ref.sequence.slice(start+length)).slice(lo,hi);
  const lengths=new Map<string,number|undefined>();
  for(let length=1;length<=Math.min(maxDeletion,ref.sequence.length);length++)for(const start of moleculeJunctionStarts(ref.sequence,ref.cut,-length,radius)){
   const key=signature(start,length);if(!lengths.has(key))lengths.set(key,length);else if(lengths.get(key)!==length)lengths.set(key,undefined);
  }
  for(const i of active){const m=models[i],origins=m.origins?.length?m.origins:[{start:m.start,weight:1}],total=origins.reduce((s,o)=>s+o.weight,0);
   for(const o of origins)if(lengths.has(signature(o.start,-m.size))&&lengths.get(signature(o.start,-m.size))===undefined)mass+=weights[i]*o.weight/total;
  }
 });
 return mass;
}
/** WT and edited molecules are evaluated on the SAME measured cycles. Neither
 * the requested indel range nor possible insertion origins mask observations. */
export function fitJunctionMixture(input:JunctionInput):JunctionFit {
 if(input.reads||input.references.some(r=>(r.cuts?.length||1)>1))return fitJunctionPaths(input);
 const {references:refs,positions,observed,maxDeletion,maxInsertion,radius,lag,theta,positionWeights,fit:fitter}=input;
 if(positions.length*4!==observed.length||!positions.length)throw new Error('Invalid junction-fit observations.');
 const rowAt=new Map(positions.map((p,i)=>[p,i]));
 // Empirical control responses retain instrument-specific dye shapes for single
 // loci. A changed junction or an unmeasured donor uses a WT triplet response,
 // learned from the control, rather than excluding that sample observation.
 const responses=refs.map(ref=>{
  const groups=new Map<string,number[][]>();
  ref.donors?.forEach((v,p)=>{if(!v.some(x=>x>0))return;for(const key of [ref.sequence.slice(p-1,p+2),ref.sequence[p]]){if(!key)continue;const g=groups.get(key)||[];g.push(v);groups.set(key,g);}});
  return new Map([...groups].map(([key,values])=>[key,Array.from({length:4},(_,b)=>{const a=values.map(v=>v[b]).sort((a,b)=>a-b);return a[Math.floor(a.length/2)];})]));
 });
 const vector=(m:JunctionModel,p:number):number[]=>{
  const ref=refs[m.reference];
  const probabilities=(at:number):number[]=>{
   if(m.size>0&&at>=m.start&&at<m.start+m.size&&m.probabilities)return m.probabilities[at-m.start];
   const base=junctionBase(ref.sequence,m,at);return Array.from(bases,b=>base==='N'?.25:Number(b===base));
  };
  if(ref.contextResponse){
   const one=(at:number)=>{
    const profile=ref.contextResponse!,wildcard=profile.sites?.[at]?.[-1];
    if(wildcard&&[-1,0,1].map(d=>junctionBase(ref.sequence,m,at+d)).join('')===ref.sequence.slice(Math.max(0,at-1),at+2))return wildcard;
    if(!m.probabilities||at+1<m.start||at-1>=m.start+m.size)return controlKnownContextEmission(junctionBase(ref.sequence,m,at-1),junctionBase(ref.sequence,m,at),junctionBase(ref.sequence,m,at+1),profile,at);
    return controlContextProbability(probabilities(at-1),probabilities(at),probabilities(at+1),profile,at);
   };
   return lag===0?one(p):mixPhase(one(p),one(phaseNeighbor(p,lag)),lag);
  }
  if(!ref.donors){
   if(ref.transitionResponse){const one=(at:number)=>transitionProbabilityEmission(probabilities(at-1),probabilities(at),probabilities(at+1),theta,ref.transitionResponse!,ref.baseYields||[1,1,1,1]),a=one(p),b=one(phaseNeighbor(p,lag));return mixPhase(a,b,lag);}
   const yieldAt=(at:number)=>{const left=probabilities(at-1),current=probabilities(at),right=probabilities(at+1);if(theta.length>8){let yieldValue=0;for(let l=0;l<4;l++)if(left[l])for(let b=0;b<4;b++)if(current[b])for(let r=0;r<4;r++)if(right[r])yieldValue+=left[l]*current[b]*right[r]*Math.exp((theta[l]||0)+(theta[4+r]||0)+(theta[8+16*l+4*b+r]||0));return yieldValue||1;}const a=left.reduce((s,v,b)=>v?s+v*Math.exp(theta[b]||0):s,0),b=right.reduce((s,v,b)=>v?s+v*Math.exp(theta[4+b]||0):s,0);return (a||1)*(b||1);};
   const current=probabilities(p),previous=probabilities(phaseNeighbor(p,lag)),a=yieldAt(p),b=yieldAt(phaseNeighbor(p,lag));
   return mixPhase(current.map(v=>v*a),previous.map(v=>v*b),lag);
  }
  const one=(at:number)=>{
   const base=junctionBase(ref.sequence,m,at);if(!base)return [0,0,0,0];
   if(m.size>0&&at>=m.start&&at<m.start+m.size&&m.probabilities)return probabilities(at).reduce((sum,w,b)=>{const response=responses[m.reference].get(bases[b])||Array.from(bases,x=>Number(x===bases[b])),total=response.reduce((s,v)=>s+v,0);return sum.map((v,k)=>v+w*response[k]/Math.max(1e-12,total));},[0,0,0,0]);
   const donor=at<m.start?at:m.size<0?at-m.size:at>=m.start+m.size?at-m.size:-1;
   const context=[-1,0,1].map(d=>junctionBase(ref.sequence,m,at+d)).join('');
   if(donor>=0&&ref.donors![donor]?.some(v=>v>0)&&context===ref.sequence.slice(donor-1,donor+2))return ref.donors![donor];
   const response=responses[m.reference].get(context)||responses[m.reference].get(base)||Array.from(bases,b=>Number(b===base));
   const sum=response.reduce((s,v)=>s+v,0);return response.map(v=>sum?v/sum:0);
  };
  const current=one(p),previous=one(phaseNeighbor(p,lag));return mixPhase(current,previous,lag);
 };
 const column=(m:JunctionModel)=>positions.flatMap(p=>vector(m,p));
 const models:JunctionModel[]=[];
 refs.forEach((ref,h)=>{
  for(let size=-maxDeletion;size<=maxInsertion;size++){
   if(input.classes&&!input.classes.some(c=>c.reference===h&&c.size===size))continue;
   const start=size===0?ref.cut:Math.max(0,ref.cut-Math.floor(Math.max(0,-size)/2));
   const inserted=size>0?'N'.repeat(size):'';
   models.push({reference:h,size,start,inserted,probabilities:size>0?Array.from({length:size},()=>[.25,.25,.25,.25]):undefined});
  }
 });
 interface Origin {start:number;weight:number;probabilities?:number[][]}
 interface Group {profiled?:boolean;indices:number[];starts:number[];origins:Origin[];rows:number[]}
 const groups:Group[]=[],byKey=new Map<string,Group>();
 models.forEach((m,i)=>{
  if(!m.size)return;
  const ref=refs[m.reference],starts=moleculeJunctionStarts(ref.sequence,ref.cut,m.size,radius);if(!starts.length)return;const lo=starts[0]-4,hi=starts[starts.length-1]+Math.max(0,-m.size)+5;
  // When the local WT sequence and read-axis coordinates are identical,
  // fluorescence cannot phase an insertion/onset to a distant homeolog SNP.
  // Estimate the observable origin/base marginal once, not an arbitrary hard
  // origin per homeolog. Different local reference sequences remain separate.
  const key=[m.size,ref.cut,ref.donors?m.reference:ref.sequence.slice(Math.max(0,lo),hi)].join(':');
  let group=byKey.get(key);
  if(!group){
   const centre=starts.reduce((best,p)=>Math.abs(p-m.start)<Math.abs(best-m.start)?p:best,starts[0]);
   const origins=starts.map(start=>({start,weight:Number(start===centre),probabilities:m.size>0?Array.from({length:m.size},()=>[.25,.25,.25,.25]):undefined}));
   const low=starts[0]-2,high=starts[starts.length-1]+Math.max(0,m.size)+3,rows:number[]=[];for(let p=low;p<high;p++){const at=rowAt.get(p);if(at!==undefined)rows.push(at);}
   group={indices:[],starts,origins,rows};byKey.set(key,group);groups.push(group);
  }
  group.indices.push(i);
 });
 const mixedColumn=(i:number,g:Group,previous?:number[])=>{
  const m=models[i],first=g.origins.find(o=>o.weight>0)!,out=previous?[...previous]:column({...m,start:first.start,probabilities:first.probabilities});
  for(const row of g.rows){const v=[0,0,0,0];for(const o of g.origins)if(o.weight>1e-12){const at=vector({...m,start:o.start,probabilities:o.probabilities},positions[row]);for(let b=0;b<4;b++)v[b]+=o.weight*at[b];}for(let b=0;b<4;b++)out[4*row+b]=v[b];}
  return out;
 };
 let columns=models.map(column);for(const g of groups)for(const i of g.indices)columns[i]=mixedColumn(i,g,columns[i]);
 let fit=fitter(columns,observed,positionWeights),iterations=0;
 const error=(cols:number[][],y:number[],weights:number[],precision?:number[])=>{let result=0;for(let row=0;row<y.length/4;row++){const v=[0,0,0,0];for(let i=0;i<cols.length;i++)if(weights[i])for(let b=0;b<4;b++)v[b]+=weights[i]*cols[i][4*row+b];const total=v.reduce((s,x)=>s+x,0);for(let b=0;b<4;b++)result+=(precision?.[row]??1)*(y[4*row+b]-(total?v[b]/total:0))**2;}return result;};
 const objective=(cols:number[][],weights:number[])=>error(cols,observed,weights,positionWeights);
 let best=objective(columns,fit.weights);const objectiveHistory=[best];
 for(let round=0;round<12;round++){
  iterations++;
  const sum=Array(observed.length).fill(0);for(let a=0;a<columns.length;a++)if(fit.weights[a])for(let k=0;k<sum.length;k++)sum[k]+=fit.weights[a]*columns[a][k];
  for(const g of groups){
   const mass=g.indices.reduce((sum,i)=>sum+fit.weights[i],0);if(mass<1e-8||!g.rows.length)continue;
   const representative=models[g.indices[0]],other=g.rows.flatMap(row=>Array.from({length:4},(_,b)=>sum[4*row+b]-g.indices.reduce((s,i)=>s+fit.weights[i]*columns[i][4*row+b],0)));
   const variants=g.origins.map(o=>g.rows.flatMap(row=>vector({...representative,start:o.start,probabilities:o.probabilities},positions[row])));
   const localY=g.rows.flatMap(i=>observed.slice(4*i,4*i+4)),localWeights=g.rows.map(i=>positionWeights?.[i]??1);
   const refine=(o:number,contribution:number,current:number[],sweeps:number)=>{
    const origin=g.origins[o];if(contribution<1e-8)return;
    for(let sweep=0;sweep<sweeps;sweep++){
     let delta=0;
     for(let j=0;j<representative.size;j++){
      const centre=origin.start+j,ids=g.rows.flatMap((row,k)=>positions[row]>=centre-2&&positions[row]<=centre+2?[k]:[]);if(!ids.length)continue;
      const trials=Array.from({length:4},(_,b)=>({...representative,start:origin.start,probabilities:origin.probabilities!.map((v,k)=>k===j?Array.from({length:4},(_,x)=>Number(x===b)):v)}));
      const basis=trials.map(trial=>ids.flatMap(k=>{const v=vector(trial,positions[g.rows[k]]);return v.map((x,b)=>current[4*k+b]+contribution*(x-variants[o][4*k+b]));}));
      const y=ids.flatMap(k=>localY.slice(4*k,4*k+4)),w=ids.map(k=>localWeights[k]),previous=origin.probabilities![j],estimate=fitEmissionMixture(basis,y,w,24);
      if(error(basis,y,estimate.weights,w)>error(basis,y,previous,w)+1e-12)continue;
      delta=Math.max(delta,...estimate.weights.map((v,b)=>Math.abs(v-previous[b])));origin.probabilities=origin.probabilities!.map((v,k)=>k===j?estimate.weights:v);
      for(const k of ids){const v=vector({...representative,start:origin.start,probabilities:origin.probabilities},positions[g.rows[k]]);for(let b=0;b<4;b++){current[4*k+b]+=contribution*(v[b]-variants[o][4*k+b]);variants[o][4*k+b]=v[b];}}
     }
     if(delta<1e-6)break;
    }
   };
   // Profile every insertion origin before assigning mass: an origin with an
   // initially unknown inserted sequence must not lose solely to its warm start.
   if(representative.size>0&&!g.profiled){
    for(let o=0;o<g.origins.length;o++)refine(o,mass,other.map((v,k)=>v+mass*variants[o][k]),5);
    g.profiled=true;
   }
   const originColumns=variants.map(c=>c.map((v,k)=>other[k]+mass*v)),previousOriginWeights=g.origins.map(o=>o.weight),ordered=null;
   const originWeights=ordered||fitEmissionMixture(originColumns,localY,localWeights,24).weights;
   if(error(originColumns,localY,originWeights,localWeights)<=error(originColumns,localY,previousOriginWeights,localWeights)+1e-12)g.origins.forEach((o,i)=>o.weight=originWeights[i]);
   if(representative.size>0){const current=other.map((v,k)=>v+mass*variants.reduce((s,c,i)=>s+g.origins[i].weight*c[k],0));for(let o=0;o<g.origins.length;o++)refine(o,mass*g.origins[o].weight,current,3);}
   for(const i of g.indices){const next=mixedColumn(i,g,columns[i]);for(const row of g.rows)for(let b=0;b<4;b++)sum[4*row+b]+=fit.weights[i]*(next[4*row+b]-columns[i][4*row+b]);columns[i]=next;}
  }
  const next=fitter(columns,observed,positionWeights);if(objective(columns,next.weights)<=objective(columns,fit.weights)+1e-12)fit=next;const value=objective(columns,fit.weights);objectiveHistory.push(value);
  if(best-value<1e-8){best=value;break;}best=value;
 }
 for(const g of groups){const largest=g.origins.reduce((a,b)=>a.weight>b.weight?a:b);for(const i of g.indices)models[i]={...models[i],start:largest.start,probabilities:largest.probabilities,origins:g.origins.filter(o=>o.weight>1e-7).map(o=>({...o}))};}
 let duplicateLengths=0;
 for(let a=0;a<models.length;a++)if(models[a].size!==0){const wt=models.findIndex(m=>m.reference===models[a].reference&&m.size===0);if(columns[a].every((v,k)=>Math.abs(v-columns[wt][k])<1e-8))duplicateLengths++;}
const predictColumns=(rows:number[])=>models.map(m=>rows.flatMap(p=>{if(!m.origins)return vector(m,p);const v=[0,0,0,0],mass=m.origins.reduce((s,o)=>s+o.weight,0);for(const o of m.origins){const x=vector({...m,start:o.start,probabilities:o.probabilities},p);for(let b=0;b<4;b++)v[b]+=o.weight*x[b]/mass;}return v;}));
 const mean=observed.reduce((s,v)=>s+v,0)/observed.length,variance=observed.reduce((s,v)=>s+(v-mean)**2,0),residual=error(columns,observed,fit.weights);
 return {...fit,r2:variance?1-residual/variance:0,rmse:Math.sqrt(residual/observed.length),models,columns,iterations,duplicateLengths,objective:best,objectiveHistory,predictColumns};
}

/** Render the actual fitted intervals, rather than a deletion forced to one
 * side of the cut. These remain model alternatives, not phased read alleles. */
export function junctionAnnotationComponents(models:JunctionModel[],weights:number[],origin:number,referenceLength:number,reverse=false):SangerAnnotation['components']{
 return models.flatMap<SangerAnnotation['components'][number]>((m,i)=>{
  if(weights[i]<=1e-10)return [];
  if(m.path)return pathAnnotations(m,weights[i],origin,referenceLength,reverse);
  const origins=m.origins?.length?m.origins:[{start:m.start,weight:1}],total=origins.reduce((s,o)=>s+o.weight,0);
  return origins.flatMap(o=>{if(o.weight<=1e-10)return [];let a=origin+o.start,b=a+Math.max(0,-m.size),inserted=m.inserted;if(reverse){[a,b]=[referenceLength-b,referenceLength-a];inserted=reverseComplement(inserted);}return [{size:m.size,fraction:weights[i]*o.weight/total,events:m.size?[{target:1,start:a,end:b,inserted:Math.max(0,m.size),insertedSequence:inserted,insertedBaseProbabilities:o.probabilities?(reverse?[...o.probabilities].reverse().map(v=>[...v].reverse()):o.probabilities.map(v=>[...v])):undefined}]:[]}];});
 });
}
