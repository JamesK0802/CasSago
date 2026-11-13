import {retrieveMarginalTruth} from './core/sanger-marginal-truth-retrieval';
/// <reference lib="webworker" />
import{certifyNativeBenchmark}from'../utils/sanger-native-signal-certificate';
import{assessSangerBenchmark}from'../utils/sanger-benchmark-assessment';
import {analyzeAutomaticSanger} from './core/sanger-automatic';
import {marginalSessions} from './core/sanger-marginal-engine';
import {marginalAlternativePage} from './core/sanger-marginal-alternatives';
import type {SangerAnalysisResult} from '../models/sanger-analysis';
let lastResult:SangerAnalysisResult|undefined;
let optimizer:Promise<any>|undefined;
const localOptimizer=()=>optimizer??=import(/* @vite-ignore */new URL('vendor/highs-1.15.4/highs.mjs',self.location.href).href).then(async module=>{const engine=await module.default(),conic=await import(/* @vite-ignore */new URL('vendor/clarabel-faer-0.11.1/clarabel.mjs',self.location.href).href);await conic.default();engine.conic=conic;return engine;});
addEventListener('message',async({data})=>{
 try{
  if(data.action==='marginal-truth'){if(!lastResult)throw Error('Analyze the original traces before checking GT retrieval.');const retrieval=retrieveMarginalTruth(lastResult,data.annotations,s=>postMessage({type:'marginal-truth-progress',requestId:data.requestId,stage:s}));postMessage({type:'marginal-truth',retrieval,requestId:data.requestId});return;}
  if(data.action==='marginal-alternatives'){
   if(!lastResult)throw Error('Analyze the original traces before requesting alternatives.');
   const page=marginalAlternativePage(lastResult,data.groupId,data.cursor,s=>postMessage({type:'marginal-progress',requestId:data.requestId,stage:s}));
   postMessage({type:'marginal-alternatives',page,requestId:data.requestId});return;
  }
  if(data.action==='benchmark'){
   postMessage({type:'progress',stage:'Checking ground truth after inference'});
   const native=!data.result.marginalGroups&&data.result.diagnostics?.nativeRawTracesVerified?certifyNativeBenchmark(data.truth,data.controlBytes,data.sampleBytes,data.result,await localOptimizer()):null,assessment=assessSangerBenchmark(data.truth,data.references,data.result,data.threshold,native);
   postMessage({type:'benchmark',assessment});return;
  }
  if(data.action!=='automatic')throw Error('Unsupported Sanger request. Use automatic analysis.');
  if(data.groupRecovery&&data.primer)postMessage({type:'progress',stage:'Loading common sequence optimizer'});
  const solver=data.groupRecovery&&data.primer?await localOptimizer():undefined;
  // Only the active result can serve lazy queries. Release the preceding
  // physical model before constructing its replacement, even if UI rows remain.
  if(lastResult)marginalSessions.delete(lastResult);
  lastResult=undefined;
  const result=analyzeAutomaticSanger(data,(stage:string)=>postMessage({type:'progress',stage}),solver);
  lastResult=result;postMessage({type:'result',result});
 }catch(error){postMessage({type:data.action==='marginal-truth'?'marginal-truth-error':'error',requestId:data.requestId,message:error instanceof Error?error.message:String(error)});}
});
