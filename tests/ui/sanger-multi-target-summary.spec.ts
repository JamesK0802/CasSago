import '@angular/compiler';
import {it,expect} from 'vitest';
import {SangerAnalysisPanelComponent} from '../../src/app/pages/sanger-analysis-panel/sanger-analysis-panel.component';
import {benchmarkMetrics,parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
it('keeps target-local WT, substitution and indel metrics separate on the same molecule axis',()=>{
 const p=new SangerAnalysisPanelComponent({markForCheck(){}}as any),components=[{size:0,pathSizes:[0,0],fraction:.6,events:[]},{size:0,pathSizes:[0,0],fraction:.2,events:[{target:1,start:150,end:151,inserted:1,insertedSequence:'C'}]},{size:-1,pathSizes:[0,-1],fraction:.2,events:[{target:2,start:300,end:301,inserted:0,insertedSequence:''}]}];
 p.result={referenceName:'R',cutPosition:150,diagnostics:{algorithmVersion:'sanger-sequence-groups-v3'},pooled:{amplicons:[{name:'R',share:1}]},annotations:[{name:'R',reference:'A'.repeat(500),cut:150,cuts:[150,300],share:1,estimable:true,sequenceGrouping:true,components}]}as any;
 expect(p.resultSummary[0]).toMatchObject({wt:.8,inFrame:0,outOfFrame:0});expect(p.resultSummary[1]).toMatchObject({wt:.8,inFrame:0,outOfFrame:.2});expect(p.substitutionFor(0,0)).toBe(.2);expect(p.substitutionFor(0,1)).toBe(0);expect(p.spectrumComponents).toEqual([{size:0,fraction:1}]);
 p.selectResultTarget(1);expect(p.spectrumComponents).toEqual([{size:-1,fraction:.2},{size:0,fraction:.8}]);
 p.result!.annotations![0].estimable=false;expect(p.resultSummary[0].wt).toBeNull();expect(p.substitutionFraction).toBeNull();
});
it('compares only the selected target with truth and leaves CSV substitution attribution unresolved',()=>{
 const base={reference:'R',sourceReference:'R',share:100,sequence:'',inputSequence:''},truth=parseBenchmarkTruth([{...base,allele:'WT',percent:60,totalPercent:60,pattern:'0/0',description:'WT',events:[]},{...base,allele:'Sub',percent:20,totalPercent:20,pattern:'0/0',description:'A→C',substitution:true,events:[{kind:'substitution',target:0}]},{...base,allele:'Del',percent:20,totalPercent:20,pattern:'0/-1',description:'T2 deletion',events:[{kind:'deletion',target:1}]}]);
 expect(benchmarkMetrics(truth,'R',0,true,true)).toMatchObject({wt:.8,substitution:.2,inFrame:0,outOfFrame:0});expect(benchmarkMetrics(truth,'R',1,true,true)).toMatchObject({wt:.8,substitution:0,outOfFrame:.2});truth.truth[1].events=[];expect(benchmarkMetrics(truth,'R',0,true,true)).toBeNull();truth.truth[1].substitution=false;truth.truth[1].description='No edit';truth.truth[2].events=[];expect(benchmarkMetrics(truth,'R',1,true,true)?.outOfFrame).toBe(.2);
});
