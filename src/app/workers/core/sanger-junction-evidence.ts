import {pathKey,isWildType} from './sanger-junction-paths';
import {emissionFamilies} from './sanger-signal';
import {fitJunctionMixture,JunctionFit} from './sanger-junction-mixture';
import {fitEmissionMixture} from './sanger-emission';

export interface JunctionPredictionEvidence {validationR2?:number[];foldGains:number[];lengthFoldGains:number[];passed:boolean;lengthPassed:boolean}
type ForwardInput=Omit<Parameters<typeof fitJunctionMixture>[0],'fit'|'classes'>;
const key=(m:{reference:number;size:number;sizes?:number[]})=>`${m.reference}:${pathKey(m)}`;
/** Relearn the SAME physical junction model on complementary contiguous cycle
 * blocks. Insertion bases and origin allocations are never copied from the
 * full-data fit. Validate selected lengths against a training-only null refit.
 * This is a conditional predictive guard, not a probability of biological truth.
 * All cycles still enter the final point estimate and one validation fold.
 */
export function junctionPredictionEvidence(input:ForwardInput,full:JunctionFit,targets:number[]):Map<number,JunctionPredictionEvidence>{
 const answer=new Map(targets.map(i=>[i,{validationR2:[] as number[],foldGains:[] as number[],lengthFoldGains:[] as number[],passed:false,lengthPassed:false}]));
 if(!targets.length)return answer;
 const families=emissionFamilies(full.columns,1e-6);
 const active=full.models.flatMap((m,i)=>full.weights[i]>1e-8||isWildType(m)?[m]:[]).map(m=>({reference:m.reference,size:m.size,sizes:m.sizes}));
 const block=Math.max(2,Math.ceil(Math.cbrt(input.positions.length)));
 const phase=input.positions.map(p=>Math.floor(p/block)%2);
 const loss=(columns:number[][],observed:number[],weights:number[],precision:number[])=>{
  let value=0;for(let p=0;p<observed.length/4;p++){
   const predicted=[0,0,0,0];for(let i=0;i<columns.length;i++)if(weights[i])for(let b=0;b<4;b++)predicted[b]+=weights[i]*columns[i][4*p+b];
   const total=predicted.reduce((s,v)=>s+v,0);for(let b=0;b<4;b++)value+=precision[p]*(observed[4*p+b]-(total?predicted[b]/total:0))**2;
  }return value;
 };
 for(let fold=0;fold<2;fold++){
  const train=input.positions.flatMap((_,i)=>phase[i]!==fold?[i]:[]),valid=input.positions.flatMap((_,i)=>phase[i]===fold?[i]:[]);
  if(train.length<30||valid.length<30)continue;
  const y=train.flatMap(i=>input.observed.slice(4*i,4*i+4)),weights=train.map(i=>input.positionWeights?.[i]??1);
  const fitted=fitJunctionMixture({...input,classes:active,positions:train.map(i=>input.positions[i]),observed:y,positionWeights:weights,fit:(c,y,w)=>fitEmissionMixture(c,y,w,24)});
  const validationColumns=fitted.predictColumns(valid.map(i=>input.positions[i])),validation=valid.flatMap(i=>input.observed.slice(4*i,4*i+4)),precision=valid.map(i=>input.positionWeights?.[i]??1),fullLoss=loss(validationColumns,validation,fitted.weights,precision);
  const mean=validation.reduce((s,v)=>s+v,0)/validation.length,variance=validation.reduce((s,v,k)=>s+precision[k>>2]*(v-mean)**2,0),validationR2=variance?1-fullLoss/variance:0;for(const e of answer.values())e.validationR2.push(validationR2);
  const lengthLosses=new Map<string,number>();
  for(const index of targets){
   const id=fitted.models.findIndex(m=>key(m)===key(full.models[index]));
   if(id<0){answer.get(index)!.foldGains.push(0);answer.get(index)!.lengthFoldGains.push(0);continue;}
   const family=families.find(g=>g.includes(index))||[index],exclude=new Set(family.filter(i=>!isWildType(full.models[i])).map(i=>key(full.models[i])));
   const keep=fitted.models.flatMap((m,i)=>exclude.has(key(m))?[]:[i]);
   const nullFit=fitEmissionMixture(keep.map(i=>fitted.columns[i]),y,weights,24);
   const nullLoss=loss(keep.map(i=>validationColumns[i]),validation,nullFit.weights,precision);
   answer.get(index)!.foldGains.push(nullLoss-fullLoss);
   const size=pathKey(full.models[index]);
   if(!lengthLosses.has(size)){
    const same=full.models.flatMap((m,i)=>pathKey(m)===size?[i]:[]),aliases=new Set(families.filter(g=>g.some(i=>same.includes(i))).flat()),excluded=new Set([...aliases].filter(i=>!isWildType(full.models[i])).map(i=>key(full.models[i]))),keepLength=fitted.models.flatMap((m,i)=>excluded.has(key(m))?[]:[i]),nullLength=fitEmissionMixture(keepLength.map(i=>fitted.columns[i]),y,weights,24);
    lengthLosses.set(size,loss(keepLength.map(i=>validationColumns[i]),validation,nullLength.weights,precision));
   }
   answer.get(index)!.lengthFoldGains.push(lengthLosses.get(size)!-fullLoss);
  }
 }
 for(const evidence of answer.values()){evidence.passed=evidence.foldGains.length===2&&evidence.foldGains.every(v=>v>1e-8*input.positions.length);evidence.lengthPassed=evidence.lengthFoldGains.length===2&&evidence.lengthFoldGains.every(v=>v>1e-8*input.positions.length);}
 return answer;
}
