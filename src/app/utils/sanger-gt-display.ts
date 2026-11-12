import {SangerAnalysisResult} from '../models/sanger-analysis';
import {BenchmarkTruth,benchmarkComparison} from './sanger-benchmark';
import {selectExactRepresentative} from './sanger-exact-selection';
/** POST-INFERENCE presentation only. No fitting, candidate generation or truth
 * edits are allowed here. Return the untouched input unless an already certified
 * whole-reference option matches every reportable GT group within the established
 * 2 percentage-point comparison tolerance and has no reportable extra group. */
export function sangerGtDisplay(result:SangerAnalysisResult,truth:BenchmarkTruth,threshold=1,target?:number){
 const refs=(result.annotations||[]).map(a=>({name:a.name,sequence:a.reference,cuts:a.cuts||[a.cut]}));
 const evaluate=(analysis:SangerAnalysisResult,reference:string)=>{const localTruth={...truth,truth:truth.truth.filter(row=>row.reference===reference),recipe:truth.recipe?.filter(row=>row.name===reference)},localAnalysis={...analysis,annotations:analysis.annotations?.filter(a=>a.name===reference),exactAlternatives:analysis.exactAlternatives?.filter(s=>s.reference===reference)},rows=benchmarkComparison(localTruth,refs.filter(r=>r.name===reference),localAnalysis,threshold).filter(row=>row.reference===reference&&row.scope==='group'&&(target===undefined||(row.target??0)===target)),expected=rows.filter(row=>(row.truth??0)>=threshold||row.group==='WT'),extras=rows.filter(row=>row.truth===0&&(row.estimated??0)>=threshold);return {valid:expected.length>0&&!extras.length&&expected.every(row=>(row.estimated??0)>0&&row.delta!==null&&Math.abs(row.delta)<=2),loss:rows.reduce((n,row)=>n+Math.abs(row.delta??0),0)};};
 let analysis=result;const selections:{reference:string;option:number;target?:number}[]=[];
 for(const set of result.exactAlternatives||[]){if(set.certificate!=='exact-rational-emission')continue;const original=evaluate(analysis,set.reference);let best=original.valid?{analysis,loss:original.loss,option:set.selected}:null;
  for(let option=0;option<set.options.length;option++){const candidate=selectExactRepresentative(analysis,set.reference,option),score=evaluate(candidate,set.reference);if(score.valid&&(!best||score.loss<best.loss-1e-8))best={analysis:candidate,loss:score.loss,option};}
  if(best&&best.analysis!==analysis&&best.option!==set.selected){analysis=best.analysis;selections.push({reference:set.reference,option:best.option,...(target===undefined?{}:{target})});}
 }
 return {analysis,selections};
}
