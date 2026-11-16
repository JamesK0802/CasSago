import '@angular/compiler';
import {it,expect} from 'vitest';
import {SangerAnalysisPanelComponent} from '../../src/app/pages/sanger-analysis-panel/sanger-analysis-panel.component';
import type {SangerAnalysisResult,SangerMarginalGroup} from '../../src/app/models/sanger-analysis';
function panel(){
 const p=new SangerAnalysisPanelComponent({markForCheck(){}}as any),g:SangerMarginalGroup={id:'zero',reference:'R',target:1,size:0,fraction:1,fractionRange:[1,1],baselineReadOffset:0,bases:[{position:2,estimate:[0,.8,.2,0],lower:[0,.79,.19,0],upper:[0,.81,.21,0]}],unmodifiedRange:[.79,.81],substitutionRange:[.19,.21],scope:'read-cycle-by-length',status:'bounded'};
 p.result={annotations:[{name:'R',reference:'AACCAA',cut:1,cuts:[1],direction:'forward',share:1,estimable:true,components:[],marginalGroups:[g]}],marginalGroups:[g]}as unknown as SangerAnalysisResult;return p;
}
it('does not present a quantified class-0 SNP fraction as the all-length Sub result',()=>{
 const p=panel();expect(p.paneSubstitution(false,0,0)).toBeNull();expect(p.fractionLabel(0,0,'substitution',.2)).toBe('—');expect(p.fractionLabel(0,0,'wt',null)).toBe('80.00%');
});
it('uses the actual all-length Sub envelope even when the class-0 range differs',()=>{
 const p=panel();p.result!.substitutionSummaries=[{reference:'R',target:1,scope:'all-length-groups',method:'joint-observable-departure-envelope',fractionRange:[.0484,.0677],zeroFeasible:false}];expect(p.fractionLabel(0,0,'substitution',null)).toBe('4.84–6.77%');expect(p.paneSubstitution(false,0,0)).toBeNull();
});
