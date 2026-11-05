import '@angular/compiler';
import {afterEach,expect,it,vi} from 'vitest';
import {readFileSync} from 'node:fs';
import {SangerDesignPageComponent} from '../../src/app/pages/sanger-design-page/sanger-design-page.component';
import {guideCut} from '../../src/app/utils/sanger-design.utils';
import {designSangerPrimers} from '../../src/app/workers/core/sanger-primer-designer';
import {makeAnalysisTemplate} from '../../src/app/utils/sanger-design-export';
const create=()=>new SangerDesignPageComponent({markForCheck:vi.fn()}as any);
afterEach(()=>vi.unstubAllGlobals());
it('limits each gene to two named guides and preserves the remaining target on removal',()=>{const p=create(),row=p.rows[0];p.addTarget(row);p.targetRows(row)[1].name='Distal';p.targetRows(row)[1].guide='A'.repeat(20);p.addTarget(row);expect(p.targetRows(row)).toHaveLength(2);p.removeTarget(row,0);expect(p.targetRows(row)[0].name).toBe('Distal');p.addTarget(row);expect(new Set(p.targetRows(row).map(t=>t.name)).size).toBe(2);p.removeTarget(row,1);p.removeTarget(row,0);expect(p.targetRows(row)).toHaveLength(1);});
it('resolves two guides in their entered order and passes names and automatic settings to the design worker',()=>{
 const workers:any[]=[];vi.stubGlobal('Worker',class{postMessage=vi.fn();terminate=vi.fn();constructor(){workers.push(this);}});
 const p=create(),a='GACCTAGTGCATGACCTGAC',b='CTGACGATCGTAGCTAGTCA';p.rows=[{name:'GeneA',sequence:'ACGT'.repeat(80)+a+'AGG'+'TGCA'.repeat(80)+b+'TGG'+'ACGT'.repeat(80),targets:'',mode:'guides',targetRows:[{name:'Proximal',guide:a},{name:'Distal',guide:b}]}];
 p.settings.platform='ngs';p.settings.direction='reverse';p.settings.indelRange=31;p.search();expect(p.error).toBe('');
 const payload=workers[0].postMessage.mock.calls[0][0];expect(payload.references[0].cuts).toEqual([337,680]);expect(payload.references[0].targetNames).toEqual(['Proximal','Distal']);expect(payload.settings).toMatchObject({platform:'sanger',direction:'both',workflow:'standard-pcr',readEnd:1000,indelRange:10,requirePartner:true});p.ngOnDestroy();
});
it('keeps the example as valid named gRNAs and preserves names in the selected PCR template',()=>{
 const p=create();p.loadExample();for(const row of p.rows){const ref=(p as any).designReference(row);expect(ref.cuts[0]).toBe(guideCut(row.sequence,p.targetRows(row)[0].guide));expect(ref.targetNames).toEqual(['Target 1','Target 2']);}
 const row=p.rows[0],ref=(p as any).designReference(row);ref.targetNames=['FAD2-site','Distal'];const result=designSangerPrimers([ref],{...p.settings,workflow:'standard-pcr',primerMin:20,primerMax:20});const template=makeAnalysisTemplate(result,result.candidates[0]);expect(template.amplicons[0].targetNames).toEqual(['FAD2-site','Distal']);
},20000);
it('removes the three import buttons and moves secondary settings into a collapsed advanced section',()=>{
 const html=readFileSync('src/app/pages/sanger-design-page/sanger-design-page.component.html','utf8');for(const text of ['Import FASTA','Use analysis references','Load design','Sequencing platform','Read direction','Indel allowance'])expect(html).not.toContain(text);
 expect(html).toContain('Gene name');expect(html).toContain('Gene sequence');expect(html).toContain('Target name');expect(html).toContain('Target sequence');expect(html).toContain('<details class="advanced-design"><summary>Advanced settings</summary>');expect(html.indexOf('Minimum target offset')).toBeGreaterThan(html.indexOf('Advanced settings'));
});

it('keeps primer-quality constraints internal and exposes only essential advanced controls',()=>{const p=create();p.settings.primerNm=999;p.changed();expect(p.settings).toMatchObject({primerMin:20,primerMax:20,tmMin:52,tmMax:68,gcMin:30,gcMax:70,sodiumMm:50,primerNm:250});const html=readFileSync('src/app/pages/sanger-design-page/sanger-design-page.component.html','utf8');for(const name of ['settings.primerNm','settings.sodiumMm','settings.maxTmDifference','Specificity method','Near-binding filter'])expect(html).not.toContain(name.startsWith('settings.')?'[(ngModel)]="'+name+'"':name);expect(html).toContain('class="gene-fields compact-pair"');expect(html).toContain('rows="1"');expect(html).toContain('aria-label="Sanger read range"');expect(html).not.toContain('aria-label="PCR product length range"');});
it('runs the same single-reference design settings for both amplification choices while retaining the selected label',()=>{const p=create();p.rows=[p.blank(1)];p.setAmplification('shared');const a={...p.settings};p.setAmplification('specific');expect(p.amplification).toBe('specific');expect(p.settings).toEqual(a);p.addReference();expect(p.settings.workflow).toBe('homeolog-pcr');p.setAmplification('shared');expect(p.settings.workflow).toBe('shared-sanger');});

it('keeps an explicit primer length range editable, with a 20-nt preset',()=>{const p=create();expect([p.settings.primerMin,p.settings.primerMax]).toEqual([20,20]);p.settings.primerMin=19;p.settings.primerMax=23;p.changed();expect([p.settings.primerMin,p.settings.primerMax]).toEqual([19,23]);});

it('ignores obsolete PCR product range input while preserving the user Sanger range',()=>{const p=create();p.settings.productMin=999;p.settings.productMax=10000;p.settings.readStart=60;p.settings.readEnd=900;p.changed();expect([p.settings.productMin,p.settings.productMax]).toEqual([150,1500]);expect([p.settings.readStart,p.settings.readEnd]).toEqual([60,900]);});

it('discards cached references and settings and never stores changed inputs',()=>{const removeItem=vi.fn(),setItem=vi.fn();vi.stubGlobal('sessionStorage',{getItem:()=>JSON.stringify({rows:[{name:'Cached',sequence:'A'.repeat(900),targets:'250',mode:'cuts'}],settings:{tmMin:1}}),removeItem,setItem});const p=create();expect(p.rows.every(r=>r.sequence==='')).toBe(true);expect(p.settings.tmMin).toBe(52);p.rows[0].sequence='ACGT';p.settings.tmMin=55;p.changed();expect(setItem).not.toHaveBeenCalled();const fresh=create();expect(fresh.rows[0].sequence).toBe('');expect(fresh.settings.tmMin).toBe(52);expect(removeItem).toHaveBeenCalledWith('casmango.sanger-design.draft.v2');});

it('keeps one shared chroma preference while selecting different candidates',()=>{const p=create(),a={id:'a'}as any,b={id:'b'}as any;expect(p.showChroma).toBe(false);p.toggleCandidateTrace(a);expect(p.showChroma).toBe(true);p.selectCandidate(b);expect(p.showChroma).toBe(true);p.toggleCandidateTrace(b);expect(p.showChroma).toBe(false);p.selectCandidate(a);expect(p.showChroma).toBe(false);});
