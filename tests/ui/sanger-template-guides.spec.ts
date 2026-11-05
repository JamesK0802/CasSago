import '@angular/compiler';
import {expect,it} from 'vitest';
import {parseAutomaticSangerTemplate} from '../../src/app/utils/sanger-automatic-template';
import {resolveSangerReference} from '../../src/app/utils/sanger-design.utils';
import {reverseComplement} from '../../src/app/utils/biology.utils';
import {SangerAnalysisPanelComponent} from '../../src/app/pages/sanger-analysis-panel/sanger-analysis-panel.component';
const guide='GACCTAGTGCATGACCTGAC',sequence='A'.repeat(100)+guide+'AGG'+'C'.repeat(350);
function template(reverse=false){const seq=reverse?reverseComplement(sequence):sequence,cuts=resolveSangerReference('R',seq,guide,'guides').cuts;return {format:'casmango-sanger-analysis-v1',amplicons:[{name:'R',sequence:seq,cuts,guides:[guide]}]};}
it('loads design gRNA sequences rather than numeric coordinates on both strands',()=>{for(const reverse of [false,true]){const t=template(reverse),r=parseAutomaticSangerTemplate(t).references[0];expect(r.mode).toBe('guides');expect(r.targets).toBe(guide);expect(resolveSangerReference(r.name,r.sequence,r.targets,r.mode).cuts).toEqual(t.amplicons[0].cuts);}});
it('retains legacy coordinate-only templates and rejects inconsistent guide metadata',()=>{const t=template();delete(t.amplicons[0]as any).guides;expect(parseAutomaticSangerTemplate(t).references[0].mode).toBe('cuts');const bad=template();bad.amplicons[0].cuts[0]++;expect(()=>parseAutomaticSangerTemplate(bad)).toThrow('disagree');});
it('collapses Advanced when a new template is applied',()=>{const p=new SangerAnalysisPanelComponent({markForCheck(){}}as any);p.advancedOpen=true;p.applyTemplate(template());expect(p.advancedOpen).toBe(false);expect(p.references[0].mode).toBe('guides');});
it('shows the selected fit as point percentages while alternatives remain separate',()=>{const p=new SangerAnalysisPanelComponent({markForCheck(){}}as any);p.result={annotations:[{name:'R'}],diagnostics:{substitutionLinkage:[{reference:'R',target:1,blocks:2,fraction:.1,wtRange:[.8,.9],modelSequence:'ACGT'}]}}as any;expect(p.fractionLabel(0,0,'wt',.9)).toBe('90.00%');expect(p.fractionLabel(0,0,'substitution',.1)).toBe('10.00%');expect(p.identifiedFraction(0,0,.9)).toBe(.9);expect(p.fractionLabel(0,1,'wt',.7)).toBe('70.00%');});
