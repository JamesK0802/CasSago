import '@angular/compiler';import {describe,it,expect,vi,afterEach} from 'vitest';
import {GrnaDesignPageComponent} from '../../src/app/pages/grna-design-page/grna-design-page.component';
import {SangerDesignPageComponent} from '../../src/app/pages/sanger-design-page/sanger-design-page.component';
import {DesignHandoffService} from '../../src/app/services/design-handoff.service';
import {designGuides} from '../../src/app/workers/core/grna-designer';
import {GUIDE_DEFAULTS} from '../../src/app/models/grna-design';
const g='GACCTAGTGCATGACCTGAC',refs=()=>[{name:'Single',sequence:'N'.repeat(200)+g+'AGG'+'N'.repeat(300),selected:true,regionStart:null,regionEnd:null}];
afterEach(()=>vi.unstubAllGlobals());
describe('gRNA to primer workflow isolation',()=>{
 it('transfers a single guide with explicit cut boundaries to standard PCR without writing analysis state',()=>{
  const handoff=new DesignHandoffService(),page=new GrnaDesignPageComponent({markForCheck:vi.fn()} as any,handoff);page.result=designGuides(refs(),[],GUIDE_DEFAULTS);page.selected=page.result.candidates.find(c=>c.sequence===g)!;
  const emit=vi.spyOn(page.primerRequested,'emit');page.toPrimers();expect(emit).toHaveBeenCalledOnce();const primer=new SangerDesignPageComponent({markForCheck:vi.fn()} as any,handoff);
  expect(primer.settings.workflow).toBe('standard-pcr');expect(primer.rows).toHaveLength(1);expect(primer.targetRows(primer.rows[0])[0]).toMatchObject({guide:g,cut:217});expect(primer.analysisDesignRows[0].targets).toBe('217');expect(primer.analysisDesignRows[0].mode).toBe('cuts');expect(handoff.take()).toBeNull();expect('state' in primer).toBe(false);
 });
 it('transfers all intended shared references in target order',()=>{
  const handoff=new DesignHandoffService(),page=new GrnaDesignPageComponent({markForCheck:vi.fn()} as any,handoff),r=refs();page.result=designGuides([r[0],{...r[0],name:'Other',sequence:'N'.repeat(10)+r[0].sequence}],[],{...GUIDE_DEFAULTS,mode:'shared'});page.selected=page.result.candidates.find(c=>c.sequence===g)!;page.toPrimers();const next=handoff.take()!;expect(next.workflow).toBe('shared-sanger');expect(next.references.map(r=>r.cuts)).toEqual([[217],[227]]);
 });
 it('terminates background search and invalidates stale results on edit',()=>{
  const workers:any[]=[];vi.stubGlobal('Worker',class{terminate=vi.fn();postMessage=vi.fn();constructor(){workers.push(this);}});const page=new GrnaDesignPageComponent({markForCheck:vi.fn()} as any,new DesignHandoffService());page.references=refs();page.search();expect(page.busy).toBe(true);page.ngOnDestroy();expect(workers[0].terminate).toHaveBeenCalledOnce();page.result={} as any;page.changed();expect(page.result).toBeNull();
 });
});
