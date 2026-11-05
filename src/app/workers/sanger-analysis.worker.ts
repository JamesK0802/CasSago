/// <reference lib="webworker" />
import{certifyNativeBenchmark}from'../utils/sanger-native-signal-certificate';
import{assessSangerBenchmark}from'../utils/sanger-benchmark-assessment';
import {analyzeAutomaticSanger} from './core/sanger-automatic';
let optimizer:Promise<any>|undefined;
const localOptimizer=()=>optimizer??=import(/* @vite-ignore */new URL('vendor/highs-1.15.4/highs.mjs',self.location.href).href).then(module=>module.default());
addEventListener('message',async({data})=>{
 try{
  if(data.action==='benchmark'){
   postMessage({type:'progress',stage:'Checking ground truth after inference'});
   const solver=await localOptimizer(),native=certifyNativeBenchmark(data.truth,data.controlBytes,data.sampleBytes,data.result,solver),assessment=assessSangerBenchmark(data.truth,data.references,data.result,data.threshold,native);
   postMessage({type:'benchmark',assessment});return;
  }
  if(data.action!=='automatic')throw Error('Unsupported Sanger request. Use automatic analysis.');
  const needsOptimizer=data.groupRecovery&&data.references?.some((r:any)=>r.cuts?.length===2)&&data.maxDeletion<=10&&data.maxInsertion<=10;
  if(needsOptimizer)postMessage({type:'progress',stage:'Loading local sequence optimizer'});
  const solver=needsOptimizer?await localOptimizer():undefined;
  const result=analyzeAutomaticSanger(data,(stage:string)=>postMessage({type:'progress',stage}),solver);
  postMessage({type:'result',result});
 }catch(error){postMessage({type:'error',message:error instanceof Error?error.message:String(error)});}
});
