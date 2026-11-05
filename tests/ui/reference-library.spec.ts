import '@angular/compiler';import {describe,it,expect,vi} from 'vitest';
import {ReferenceSourceComponent} from '../../src/app/components/reference-library/reference-source.component';import {ReferenceLibraryComponent} from '../../src/app/components/reference-library/reference-library.component';
import {GrnaDesignPageComponent} from '../../src/app/pages/grna-design-page/grna-design-page.component';import {SangerDesignPageComponent} from '../../src/app/pages/sanger-design-page/sanger-design-page.component';import {DesignHandoffService} from '../../src/app/services/design-handoff.service';import {referenceChecksum} from '../../src/app/utils/reference-checksum';
const cd={markForCheck:vi.fn()},source:any={provider:'ensembl',division:'EnsemblPlants',species:'triticum_aestivum',speciesLabel:'Wheat',geneId:'GeneA',geneVersion:1,assembly:'IWGSC',chromosome:'4A',geneStart:10,geneEnd:13,strand:1,sequenceStart:10,sequenceEnd:13,upstream:0,downstream:0,retrievedAt:'2026-10-28',sha256:referenceChecksum('ACGT'),viewUrl:'https://rest.ensembl.org/lookup/id/GeneA'};
const reference=()=>({name:'GeneA',sequence:'ACGT',source:{...source}});
describe('reference import and explicit update workflow',()=>{
 it('stages a changed source without replacing user data until Apply is clicked',async()=>{
  const next={...reference(),sequence:'ACGTA',source:{...source,sha256:referenceChecksum('ACGTA'),sequenceEnd:14}},api={latest:vi.fn().mockResolvedValue(next)},c=new ReferenceSourceComponent(api as any,cd as any);c.source=source;c.sequence='ACGT';c.ngOnChanges();const emit=vi.spyOn(c.updated,'emit');await c.check();expect(c.pending).toEqual(next);expect(c.sequence).toBe('ACGT');expect(emit).not.toHaveBeenCalled();c.apply();expect(emit).toHaveBeenCalledWith(next);
 });
 it('reports unchanged sources and detects local edits without calling them source updates',async()=>{
  const c=new ReferenceSourceComponent({latest:vi.fn().mockResolvedValue(reference())} as any,cd as any);c.source=source;c.sequence='ACGA';c.ngOnChanges();expect(c.modified).toBe(true);await c.check();expect(c.pending).toBeNull();expect(c.message).toContain('unchanged');
 });
 it('drops stale asynchronous updates when the row changes',async()=>{
  let resolve:(v:any)=>void=()=>{};const c=new ReferenceSourceComponent({latest:vi.fn(()=>new Promise(r=>resolve=r))} as any,cd as any);c.source=source;c.sequence='ACGT';c.ngOnChanges();const pending=c.check();c.sequence='CCCC';c.ngOnChanges();resolve({...reference(),sequence:'TTTT'});await pending;expect(c.pending).toBeNull();expect(c.message).toBe('');
 });
 it('cancels a closed library before late responses can emit imported sequences',async()=>{
  let resolve:(v:any)=>void=()=>{};const api={importGene:vi.fn(()=>new Promise(r=>resolve=r))},c=new ReferenceLibraryComponent(api as any,cd as any);c.species=[{name:source.species,display_name:'Wheat',assembly:'IWGSC'}];c.speciesName=source.species;c.genes=[{id:'GeneA'} as any];c.selected.add('GeneA');const emit=vi.spyOn(c.imported,'emit'),request=c.importSelected();c.cancel();resolve(reference());await request;expect(emit).not.toHaveBeenCalled();
 });
 it('requires explicit species selection after a filter hides the selected cultivar',()=>{
  const c=new ReferenceLibraryComponent({} as any,cd as any);c.species=[{name:'triticum_aestivum',display_name:'Wheat',assembly:'IWGSC'}];c.speciesName='triticum_aestivum';c.speciesFilter='human';c.filterChanged();expect(c.speciesName).toBe('');expect(c.choice).toBeUndefined();
 });
 it('preserves provenance when importing and handing a guide to primer design',()=>{
  const service=new DesignHandoffService(),c=new GrnaDesignPageComponent(cd as any,service);c.references=[c.blank(1)];c.importLibrary([reference()]);expect(c.references[0].source?.sha256).toBe(source.sha256);
  c.result={references:c.references} as any;c.selected={sequence:'ACGT',sites:[{reference:0,cut:2}]} as any;c.toPrimers();const primer=new SangerDesignPageComponent(cd as any,service);expect(primer.rows[0].source?.geneId).toBe('GeneA');expect(primer.rows[0].targets).toBe('2');
 });
 it('invalidates old target positions and search results on explicit reference replacement',()=>{
  const c=new GrnaDesignPageComponent(cd as any,new DesignHandoffService());c.references=[{...c.blank(1),...reference(),name:'My gene',regionStart:2,regionEnd:3}];c.result={} as any;c.updateLibrary(0,{...reference(),sequence:'AAAA'});expect(c.references[0].regionStart).toBeNull();expect(c.references[0].name).toBe('My gene');expect(c.result).toBeNull();
  const primer=new SangerDesignPageComponent(cd as any);primer.rows=[{...reference(),targets:'2',mode:'cuts'}];primer.result={} as any;primer.updateLibrary(0,{...reference(),sequence:'AAAA'});expect(primer.rows[0].targets).toBe('');expect(primer.result).toBeNull();
 });
 it('selects single or shared mode to match a fresh library import',()=>{
  const c=new GrnaDesignPageComponent(cd as any,new DesignHandoffService());c.settings.mode='shared';c.references=[c.blank(1)];c.importLibrary([reference()]);expect(c.settings.mode).toBe('single');
  c.references=[c.blank(1)];c.importLibrary([reference(),{...reference(),name:'GeneB',source:{...source,geneId:'GeneB'}}]);expect(c.settings.mode).toBe('shared');expect(c.references.every(r=>r.selected)).toBe(true);
 });
 it('does not overwrite a running design or append duplicate imported identifiers',()=>{
  const c=new GrnaDesignPageComponent(cd as any,new DesignHandoffService());c.references=[{...c.blank(1),...reference()}];c.importLibrary([reference()]);expect(c.references).toHaveLength(1);expect(c.error).toContain('already present');c.references[0].name='Renamed locally';c.importLibrary([reference()]);expect(c.references).toHaveLength(1);c.busy=true;c.importLibrary([{...reference(),name:'B'}]);expect(c.references).toHaveLength(1);
 });
});
