import {describe,it,expect,vi,afterEach} from 'vitest';import {createHash} from 'node:crypto';
import {ReferenceLibraryService} from '../../src/app/services/reference-library.service';import {referenceChecksum} from '../../src/app/utils/reference-checksum';
const gene={id:'GeneA',species:'triticum_aestivum',object_type:'Gene',assembly_name:'IWGSC',seq_region_name:'4A',start:100,end:199,strand:-1 as const,version:2};
const choice={name:'triticum_aestivum',display_name:'Triticum aestivum',assembly:'IWGSC',strain:'Chinese Spring',release:116};
const respond=(data:any,status=200)=>({ok:status===200,status,json:async()=>data});
afterEach(()=>vi.unstubAllGlobals());
describe('live reference source validation',()=>{
 it('uses standard SHA-256 independent of HTTP/WebCrypto and normalizes DNA whitespace',()=>{
  for(const seq of ['', 'ACGT', 'ACGTN'.repeat(20000)])expect(referenceChecksum(seq)).toBe(createHash('sha256').update(seq).digest('hex'));
  expect(referenceChecksum('ac gt\nN')).toBe(referenceChecksum('ACGTN'));
 });
 it('fetches current species metadata without credentials or a local sequence database',async()=>{
  const request=vi.fn().mockResolvedValue(respond({species:[choice]}));vi.stubGlobal('fetch',request);const result=await new ReferenceLibraryService().species('EnsemblPlants');expect(result).toHaveLength(1);expect(request.mock.calls[0][0]).toContain('division=EnsemblPlants');expect(request.mock.calls[0][1]).toMatchObject({cache:'no-store',credentials:'omit'});
 });
 it('pins species, assembly and reverse-strand flanks in the actual sequence request',async()=>{
  const sequence='ACGT'.repeat(40),request=vi.fn().mockResolvedValueOnce(respond(gene)).mockResolvedValueOnce(respond({assembly_name:'IWGSC',top_level_region:[{name:'4A',length:500}]})).mockResolvedValueOnce(respond({seq:sequence,id:'chromosome:IWGSC:4A:60:219:-1'}));vi.stubGlobal('fetch',request);
  const out=await new ReferenceLibraryService().importGene(gene,choice,'EnsemblPlants',20,40);expect(out.sequence).toBe(sequence);expect(out.source).toMatchObject({geneId:'GeneA',geneVersion:2,sequenceStart:60,sequenceEnd:219,strand:-1,upstream:20,downstream:40,assembly:'IWGSC',species:choice.name});expect(out.source.sha256).toBe(referenceChecksum(sequence));
  const region=new URL(request.mock.calls[2][0]);expect(decodeURIComponent(region.pathname)).toContain('/4A:60:219:-1');expect(region.searchParams.get('coord_system_version')).toBe('IWGSC');expect(new URL(request.mock.calls[0][0]).searchParams.get('species')).toBe(choice.name);
 });
 it('clips flanks to chromosome bounds instead of claiming unavailable sequence',async()=>{
  vi.stubGlobal('fetch',vi.fn().mockResolvedValueOnce(respond({...gene,start:3,end:9,strand:1})).mockResolvedValueOnce(respond({assembly_name:'IWGSC',top_level_region:[{name:'4A',length:12}]})).mockResolvedValueOnce(respond({seq:'ACGTACGTACGT',id:'chromosome:IWGSC:4A:1:12:1'})));
  const out=await new ReferenceLibraryService().importGene(gene,choice,'EnsemblPlants',20,40);expect(out.source.sequenceStart).toBe(1);expect(out.source.sequenceEnd).toBe(12);
 });
 it('rejects another cultivar with a reused gene ID',async()=>{
  vi.stubGlobal('fetch',vi.fn().mockResolvedValue(respond({...gene,species:'triticum_aestivum_alchemy'})));await expect(new ReferenceLibraryService().lookup(choice.name,'GeneA')).rejects.toThrow(/different species/);
 });
 it('rejects assembly and coordinate mismatches before importing any reference',async()=>{
  vi.stubGlobal('fetch',vi.fn().mockResolvedValueOnce(respond(gene)).mockResolvedValueOnce(respond({assembly_name:'OTHER'})));await expect(new ReferenceLibraryService().importGene(gene,choice,'EnsemblPlants',0,0)).rejects.toThrow(/assembly disagree/);
  vi.stubGlobal('fetch',vi.fn().mockResolvedValueOnce(respond(gene)).mockResolvedValueOnce(respond({assembly_name:'IWGSC',top_level_region:[{name:'4A',length:500}]})).mockResolvedValueOnce(respond({seq:'A'.repeat(100),id:'chromosome:IWGSC:4B:100:199:-1'})));await expect(new ReferenceLibraryService().importGene(gene,choice,'EnsemblPlants',0,0)).rejects.toThrow(/coordinates\/assembly/);
 });
 it('rejects oversized genes before downloading DNA and restricts input IDs',async()=>{
  const request=vi.fn().mockResolvedValueOnce(respond({...gene,end:600100})).mockResolvedValueOnce(respond({assembly_name:'IWGSC',top_level_region:[{name:'4A',length:1000000}]}));vi.stubGlobal('fetch',request);await expect(new ReferenceLibraryService().importGene(gene,choice,'EnsemblPlants',0,0)).rejects.toThrow(/500 kb/);expect(request).toHaveBeenCalledTimes(2);
  await expect(new ReferenceLibraryService().lookup(choice.name,'../secret')).rejects.toThrow(/valid species/);
 });
 it('resolves gene symbols to species-qualified gene lookups, never guessing among transcripts',async()=>{
  const request=vi.fn().mockResolvedValueOnce(respond([{id:'TranscriptX',type:'transcript'},{id:'GeneA',type:'gene'}])).mockResolvedValueOnce(respond(gene));vi.stubGlobal('fetch',request);const out=await new ReferenceLibraryService().search(choice.name,'PDS','symbol');expect(out.map(g=>g.id)).toEqual(['GeneA']);expect(request.mock.calls[0][0]).toContain('object_type=gene');expect(request.mock.calls[1][0]).toContain('species=triticum_aestivum');
 });
 it('labels a direct-symbol fallback when the linked-gene endpoint is unavailable',async()=>{
  const request=vi.fn().mockResolvedValueOnce(respond({},503)).mockResolvedValueOnce(respond(gene));vi.stubGlobal('fetch',request);const found=await new ReferenceLibraryService().search(choice.name,'PDS','symbol');expect(found[0].id).toBe('GeneA');expect(found[0].lookupNote).toContain('listing unavailable');expect(request.mock.calls[1][0]).toContain('/lookup/symbol/');
 });
 it('does not bypass provider rate limits with a fallback request',async()=>{
  const request=vi.fn().mockResolvedValue(respond({},429));vi.stubGlobal('fetch',request);await expect(new ReferenceLibraryService().search(choice.name,'PDS','symbol')).rejects.toThrow(/request limit/);expect(request).toHaveBeenCalledTimes(1);
 });
 it('does not silently return partial success for a multi-ID request',async()=>{
  vi.stubGlobal('fetch',vi.fn().mockResolvedValueOnce(respond(gene)).mockResolvedValueOnce(respond({},404)));await expect(new ReferenceLibraryService().search(choice.name,'GeneA GeneMissing','id')).rejects.toThrow(/404/);
 });
 it('reports provider rate limiting and supports cancellation',async()=>{
  vi.stubGlobal('fetch',vi.fn().mockResolvedValue(respond({},429)));await expect(new ReferenceLibraryService().lookup(choice.name,'GeneA')).rejects.toThrow(/request limit/);
  vi.stubGlobal('fetch',vi.fn((_url,options)=>new Promise((_ok,reject)=>{options.signal.addEventListener('abort',()=>reject(new Error('aborted')));})));const controller=new AbortController(),p=new ReferenceLibraryService().lookup(choice.name,'GeneA',controller.signal);controller.abort();await expect(p).rejects.toThrow(/cancelled/);
 });
});
