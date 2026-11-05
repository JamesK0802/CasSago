import '@angular/compiler';
import { afterEach,describe,expect,it,vi } from 'vitest';
import { SangerDesignPageComponent } from '../../src/app/pages/sanger-design-page/sanger-design-page.component';
import { sangerDesignExample, resolveSangerReference } from '../../src/app/utils/sanger-design.utils';
import { designSangerPrimers } from '../../src/app/workers/core/sanger-primer-designer';
import { SANGER_DEFAULTS } from '../../src/app/models/sanger-design';

afterEach(()=>vi.unstubAllGlobals());
const create=()=>new SangerDesignPageComponent({markForCheck:vi.fn()} as any);
describe('Sanger design page lifecycle',()=>{
  it('runs independently without the CasMango application state',()=>{const page=create();expect(page.rows.length).toBeGreaterThan(0);expect('state' in page).toBe(false);});
  it('clears stale results when input or settings change',()=>{
    const page=create();page.result={} as any;page.selected={} as any;page.changed();expect(page.result).toBeNull();expect(page.selected).toBeNull();
  });
  it('maps a specific pair to its own reference rather than the first reference',()=>{
    const page=create();const refs=sangerDesignExample();
    page.result=designSangerPrimers(refs,{...SANGER_DEFAULTS,workflow:'homeolog-pcr',specificityStrategy:'mismatch',primerMin:20,primerMax:20});
    page.selectCandidate(page.result.candidates.find(c=>c.referenceIndex===2)!);
    expect(page.reference(0).name).toBe('Demo-D');expect(page.refIndex(0)).toBe(2);
    expect(page.mapX(refs[2].cuts[0],0)).toBeGreaterThan(0);
    page.sequenceFocus='primer';expect(page.sequenceCells(0).filter(c=>c.kind==='primer')).toHaveLength(page.selected!.sequence.length);
    page.sequenceFocus='partner';expect(page.sequenceCells(0).filter(c=>c.kind==='partner')).toHaveLength(page.selected!.partner!.sequence.length);
    page.sequenceFocus='target';expect(page.sequenceCells(0).filter(c=>c.cut)).toHaveLength(1);
    expect(page.status).toMatch(/3\/3 references/);
    const download=vi.spyOn(page as any,'download').mockImplementation(()=>{});page.exportAmplicons();
    const hit=page.selected!.sites[0],other=page.selected!.partner!.sites[0];
    expect(download.mock.calls[0][1]).toContain(refs[2].sequence.slice(Math.min(hit.start,other.start),Math.max(hit.end,other.end)));
    expect(download.mock.calls[0][1]).toContain('Demo-D');

  });
  it('aligns full-reference maps by target boundary instead of FASTA origin',()=>{
    const page=create();page.result=designSangerPrimers(sangerDesignExample(),{...SANGER_DEFAULTS,primerMin:20,primerMax:20});
    page.selectCandidate(page.result.candidates[0]);
    expect(page.mapX(page.reference(0).cuts[0],0)).toBeCloseTo(page.mapX(page.reference(2).cuts[0],2),8);
    expect(page.mapX(page.selected!.sites[0].start,0)).not.toBe(page.mapX(page.selected!.sites[2].start,2));
  });
  it('terminates work on cancel and component destruction',()=>{
    const workers:any[]=[];
    vi.stubGlobal('Worker',class{terminate=vi.fn();postMessage=vi.fn();constructor(){workers.push(this);}});
    const page=create();page.loadExample();page.search();expect(page.busy).toBe(true);page.cancel();expect(workers[0].terminate).toHaveBeenCalledOnce();expect(page.busy).toBe(false);
    page.search();page.ngOnDestroy();expect(workers[1].terminate).toHaveBeenCalledOnce();
  });
  it('receives worker results and reports empty outcomes without a successful status',()=>{
    const workers:any[]=[];vi.stubGlobal('Worker',class{terminate=vi.fn();postMessage=vi.fn();constructor(){workers.push(this);}});
    const page=create();page.loadExample();page.search();workers[0].onmessage({data:{type:'result',result:{candidates:[],stats:{commonUnique:0},complete:true,settings:SANGER_DEFAULTS,references:sangerDesignExample()}}});
    expect(page.busy).toBe(false);expect(page.status).toBe('No candidate under these settings');expect(page.noResultReason).toMatch(/exact site/);
  });
});

describe('primer option browsing and handoff',()=>{
 it('loads candidates in tens and clears the list on input changes',()=>{const page=create();page.result={candidates:Array.from({length:37},(_,i)=>({id:String(i)}))} as any;expect(page.visibleCandidates.length).toBe(10);page.loadMore();expect(page.visibleCandidates.length).toBe(20);page.loadMore();page.loadMore();expect(page.visibleCandidates.length).toBe(37);page.changed();expect(page.visibleCount).toBe(10);});
 it('opens the analysis tab with selected PCR products instead of full gene sequences',()=>{const page=create(),ref=sangerDesignExample()[0];page.result=designSangerPrimers([ref],{...SANGER_DEFAULTS,workflow:'standard-pcr',primerMin:20,primerMax:20});page.selected=page.result.candidates[0];const applyTemplate=vi.fn();page.analysisPanel={applyTemplate} as any;page.useInAnalysis();expect(page.activeSection).toBe('analysis');expect(applyTemplate).toHaveBeenCalledOnce();const v=applyTemplate.mock.calls[0][0];expect(v.amplicons[0].sequence.length).toBe(page.selected.partner!.products[0]);expect(v.amplicons[0].sequence.length).toBeLessThan(ref.sequence.length);});
});

 it('loads three demo references with two distinct, resolvable guide targets each',()=>{const page=create();page.loadExample();expect(page.rows.map(r=>r.name)).toEqual(['Demo-A','Demo-B','Demo-D']);for(const row of page.rows){expect(row.targetRows).toHaveLength(2);const resolved=resolveSangerReference(row.name,row.sequence,row.targets,'guides');expect(resolved.cuts).toHaveLength(2);expect(resolved.cuts[1]-resolved.cuts[0]).toBe(250);} });
