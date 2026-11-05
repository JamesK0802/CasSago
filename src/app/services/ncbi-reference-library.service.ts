import {Injectable} from '@angular/core';
import {LibraryGene,LibrarySpecies,LibraryReference,ReferenceSource,EnsemblDivision} from '../models/reference-source';
import {referenceChecksum} from '../utils/reference-checksum';
const BASE='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
const endpoint=(method:string,params:Record<string,string|number>)=>BASE+method+'.fcgi?'+new URLSearchParams(Object.fromEntries(Object.entries({tool:'CasMango',...params}).map(([k,v])=>[k,String(v)])));
@Injectable({providedIn:'root'})
export class NcbiReferenceLibraryService {
 private nextRequest=0;
 private async request(method:string,params:Record<string,string|number>,signal?:AbortSignal,text=false):Promise<any>{
  const wait=Math.max(0,this.nextRequest-Date.now());this.nextRequest=Date.now()+wait+350;
  if(wait)await new Promise<void>(resolve=>setTimeout(resolve,wait));
  if(signal?.aborted)throw Error('Request cancelled.');
  const controller=new AbortController(),cancel=()=>controller.abort(),timer=setTimeout(cancel,30000);signal?.addEventListener('abort',cancel,{once:true});
  try{const response=await fetch(endpoint(method,params),{signal:controller.signal,credentials:'omit',cache:'no-store'});if(!response.ok)throw Error(response.status===429?'NCBI request limit reached. Retry shortly.':`NCBI returned ${response.status}.`);const data=text?await response.text():await response.json();if(data?.error||data?.esearchresult?.ERROR)throw Error(String(data.error||data.esearchresult.ERROR));return data;}
  catch(e){if(controller.signal.aborted)throw Error(signal?.aborted?'Request cancelled.':'NCBI timed out. Retry or choose another source.');if(e instanceof TypeError)throw Error('Could not reach NCBI. Retry or choose another source.');throw e;}
  finally{clearTimeout(timer);signal?.removeEventListener('abort',cancel);}
 }
 async genome(name:string,signal?:AbortSignal):Promise<LibrarySpecies>{
  const scientific=name.trim().replace(/_/g,' ');if(!/^[A-Za-z][A-Za-z .()-]{1,100}$/.test(scientific))throw Error('Enter a common or scientific species name.');
  const search=await this.request('esearch',{db:'taxonomy',term:scientific+'[Scientific Name]',retmode:'json',retmax:2},signal),ids=search.esearchresult?.idlist;
  if(!Array.isArray(ids)||ids.length!==1||Number(search.esearchresult.count)!==1)throw Error('NCBI species not uniquely identified. Enter its full scientific name.');
  const data=await this.request('esummary',{db:'taxonomy',id:ids[0],retmode:'json'},signal),tax=data.result?.[ids[0]];
  if(String(tax?.scientificname||'').toLowerCase()!==scientific.toLowerCase()||String(tax.taxid??tax.uid)!==String(ids[0]))throw Error('NCBI returned unverified species metadata.');
  return {name:String(ids[0]),display_name:tax.scientificname,assembly:'RefSeq genomic accession selected per gene'};
 }
 private gene(data:any,taxonomyId:string):LibraryGene{
  if(String(data?.organism?.taxid)!==taxonomyId)throw Error('The NCBI gene belongs to a different species.');
  if(!/^\d+$/.test(String(data.uid))||data.status||data.currentid&&data.currentid!=='0')throw Error('The NCBI gene record is discontinued or replaced. Search its current Gene ID.');
  const loci=(data.genomicinfo||[]).filter((g:any)=>/^(NC|NW|NT)_\d+\.\d+$/.test(g.chraccver)&&Number.isInteger(g.chrstart)&&Number.isInteger(g.chrstop)&&g.chrstart>=0&&g.chrstop>=0);
  if(loci.length!==1)throw Error('This NCBI gene lacks one unambiguous RefSeq genomic locus. Use a verified accession/FASTA or Ensembl.');
  const g=loci[0];return {id:String(data.uid),display_name:data.name,description:data.description,species:taxonomyId,taxonomyId,assembly_name:'RefSeq '+g.chraccver,seq_region_name:g.chraccver,start:Math.min(g.chrstart,g.chrstop)+1,end:Math.max(g.chrstart,g.chrstop)+1,strand:g.chrstart<=g.chrstop?1:-1,object_type:'Gene'};
 }
 async lookup(taxonomyId:string,id:string,signal?:AbortSignal):Promise<LibraryGene>{
  if(!/^\d+$/.test(id)||!/^\d+$/.test(taxonomyId))throw Error('NCBI requires a numeric Gene ID and taxonomy ID. Ensembl IDs require the Ensembl source.');
  const data=await this.request('esummary',{db:'gene',id,retmode:'json'},signal);const gene=this.gene(data.result?.[id],taxonomyId);if(gene.id!==id)throw Error('NCBI returned a different Gene ID.');return gene;
 }
 async search(taxonomyId:string,query:string,mode:'id'|'symbol',signal?:AbortSignal):Promise<LibraryGene[]>{
  if(!/^\d+$/.test(taxonomyId))throw Error('Choose a verified NCBI species.');
  let ids:string[];
  if(mode==='id'){ids=[...new Set(query.trim().split(/[\s,;]+/).filter(Boolean))];if(!ids.length||ids.length>8||ids.some(id=>!/^\d+$/.test(id)))throw Error('Enter 1–8 numeric NCBI Gene IDs, or choose Ensembl for Ensembl IDs.');}
  else {if(!/^[A-Za-z0-9_. -]{1,100}$/.test(query.trim()))throw Error('Enter a gene name or numeric NCBI Gene ID.');const data=await this.request('esearch',{db:'gene',term:`${query.trim()}[Gene Name] AND ${taxonomyId}[Taxonomy ID] AND alive[prop]`,retmode:'json',retmax:9},signal);if(Number(data.esearchresult?.count)>8)throw Error('More than eight genes match. Use numeric NCBI Gene IDs.');ids=data.esearchresult?.idlist;if(!Array.isArray(ids))throw Error('NCBI returned an invalid gene search.');}
  if(!ids.length)return [];const data=await this.request('esummary',{db:'gene',id:ids.join(','),retmode:'json'},signal);return ids.map(id=>{const gene=this.gene(data.result?.[id],taxonomyId);if(gene.id!==id)throw Error('NCBI returned a different Gene ID.');return gene;});
 }
 async importGene(gene:Pick<LibraryGene,'id'>,choice:LibrarySpecies,division:EnsemblDivision,upstream:number,downstream:number,signal?:AbortSignal):Promise<LibraryReference>{
  if(![upstream,downstream].every(n=>Number.isInteger(n)&&n>=0&&n<=20000))throw Error('Use whole-number flanks within 0–20,000 bp.');
  const current=await this.lookup(choice.name,gene.id,signal),data=await this.request('esummary',{db:'nuccore',id:current.seq_region_name,retmode:'json'},signal),uids=data.result?.uids,chrom=Array.isArray(uids)&&uids.length===1?data.result[uids[0]]:null;
  if(chrom?.accessionversion!==current.seq_region_name||!Number.isInteger(chrom.slen)||chrom.slen<current.end||Number(chrom.taxid)!==Number(choice.name))throw Error('NCBI could not verify the genomic accession, species and bounds.');
  const start=Math.max(1,current.start-(current.strand===1?upstream:downstream)),end=Math.min(chrom.slen,current.end+(current.strand===1?downstream:upstream));
  if(end-start+1>500000)throw Error('The gene + flanks exceeds 500 kb. Import a smaller verified region manually.');
  const params={db:'nuccore',id:current.seq_region_name,seq_start:start,seq_stop:end,strand:current.strand===1?1:2,rettype:'fasta',retmode:'text'},fasta=String(await this.request('efetch',params,signal,true)).trim(),lines=fasta.split(/\r?\n/),header=lines.shift()||'',sequence=lines.join('').replace(/\s/g,'').toUpperCase();
  const expected='>'+current.seq_region_name+':'+(current.strand===1?`${start}-${end}`:`c${end}-${start}`);
  if(header.split(/\s/)[0]!==expected||!/^[ACGTN]+$/.test(sequence)||sequence.length!==end-start+1)throw Error('NCBI returned a different accession, interval, strand or sequence length. Nothing was imported.');
  return {name:current.display_name||current.id,sequence,source:{provider:'ncbi',division,species:choice.name,speciesLabel:choice.display_name,geneId:current.id,geneVersion:null,assembly:current.assembly_name,chromosome:current.seq_region_name,geneStart:current.start,geneEnd:current.end,strand:current.strand,sequenceStart:start,sequenceEnd:end,upstream,downstream,retrievedAt:new Date().toISOString(),sha256:referenceChecksum(sequence),lookupUrl:endpoint('esummary',{db:'gene',id:current.id,retmode:'json'}),sequenceUrl:endpoint('efetch',params),viewUrl:'https://www.ncbi.nlm.nih.gov/gene/'+current.id}};
 }
 latest(source:ReferenceSource,signal?:AbortSignal){return this.importGene({id:source.geneId},{name:source.species,display_name:source.speciesLabel,assembly:source.assembly},source.division,source.upstream,source.downstream,signal);}
}
