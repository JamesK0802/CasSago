import {Component,Input,Output,EventEmitter,ChangeDetectorRef,OnDestroy} from '@angular/core';import {CommonModule} from '@angular/common';import {FormsModule} from '@angular/forms';
import {NcbiReferenceLibraryService} from '../../services/ncbi-reference-library.service';
import {ReferenceLibraryService} from '../../services/reference-library.service';import {EnsemblDivision,LibrarySpecies,LibraryGene,LibraryReference} from '../../models/reference-source';
@Component({selector:'app-reference-library',standalone:true,imports:[CommonModule,FormsModule],templateUrl:'./reference-library.component.html',styleUrl:'./reference-library.component.css'})
export class ReferenceLibraryComponent implements OnDestroy{
 @Input() simpleSearch=false;
 provider:'ensembl'|'ncbi'='ensembl';
 get activeApi(){return this.provider==='ncbi'?this.ncbi:this.api;}
 providerChanged(){this.manualSpecies='';this.species=[];this.speciesName='';this.collectionAuto=true;this.inputsChanged();}
 private manualSpecies='';
 speciesInput='Wheat';searchType:'auto'|'id'|'symbol'='auto';collectionAuto=true;
 readonly popularSpecies=['Wheat','Rice','Arabidopsis','Maize','Tomato','Canola','Tobacco','Human','Mouse','Zebrafish'];
 @Input() disabled=false;@Output() imported=new EventEmitter<LibraryReference[]>();open=false;division:EnsemblDivision='EnsemblPlants';species:LibrarySpecies[]=[];speciesName='triticum_aestivum';speciesFilter='';query='';queryMode:'id'|'symbol'='id';upstream=1500;downstream=1500;genes:LibraryGene[]=[];selected=new Set<string>();loading=false;stage='';error='';private controller:AbortController|null=null;private epoch=0;
 constructor(private api:ReferenceLibraryService,private cd:ChangeDetectorRef,private ncbi:NcbiReferenceLibraryService=new NcbiReferenceLibraryService()){}
 toggle(){this.open=!this.open;if(this.open&&!this.simpleSearch)void this.loadSpecies();else this.cancel();}
 cancel(){this.epoch++;this.controller?.abort();this.controller=null;this.loading=false;this.stage='';}
 ngOnDestroy(){this.cancel();}
 private begin(stage:string){this.cancel();this.loading=true;this.stage=stage;this.error='';this.controller=new AbortController();return {epoch:this.epoch,signal:this.controller.signal};}
 private done(epoch:number){if(epoch===this.epoch){this.loading=false;this.stage='';this.cd.markForCheck();}}
 async loadSpecies(){const request=this.begin('Loading current species catalogue…');this.genes=[];this.selected.clear();try{const list=await this.api.species(this.division,request.signal);if(request.epoch!==this.epoch)return;this.species=list;if(!list.some(s=>s.name===this.speciesName))this.speciesName=list.find(s=>s.name==='homo_sapiens')?.name||list[0]?.name||'';}catch(e){if(request.epoch===this.epoch)this.error=e instanceof Error?e.message:String(e);}finally{this.done(request.epoch);}}
 get filteredSpecies(){const q=this.speciesFilter.toLowerCase();return this.species.filter(s=>[s.name,s.display_name,s.strain,s.assembly,...(s.aliases||[])].some(v=>v?.toLowerCase().includes(q)));}
 filterChanged(){if(!this.filteredSpecies.some(s=>s.name===this.speciesName)){this.speciesName='';this.inputsChanged();}}
 get choice(){return this.species.find(s=>s.name===this.speciesName);}
 inputsChanged(){this.cancel();this.genes=[];this.selected.clear();this.error='';}
 async simpleSearchGenes(){
  const name=this.speciesInput.trim().toLowerCase().replace(/_/g,' '),query=this.query.trim();if(!name||!query){this.error='Enter a species and gene name or ID.';return;}
  const aliases:Record<string,[string,EnsemblDivision]>={wheat:['triticum_aestivum','EnsemblPlants'],rice:['oryza_sativa','EnsemblPlants'],arabidopsis:['arabidopsis_thaliana','EnsemblPlants'],maize:['zea_mays','EnsemblPlants'],tomato:['solanum_lycopersicum','EnsemblPlants'],canola:['brassica_napus','EnsemblPlants'],tobacco:['nicotiana_tabacum','EnsemblPlants'],human:['homo_sapiens','EnsemblVertebrates'],mouse:['mus_musculus','EnsemblVertebrates'],zebrafish:['danio_rerio','EnsemblVertebrates']};
  const preferred=aliases[name]||Object.values(aliases).find(([id])=>id.replace(/_/g,' ')===name);if(this.collectionAuto&&preferred)this.division=preferred[1];
  const request=this.begin('Finding species…');this.genes=[];this.selected.clear();
  try{
   const direct=this.manualSpecies||preferred?.[0]||(/^[a-z]+ [a-z]+$/.test(name)?name.replace(/ /g,'_'):'');
   let resolved=false;
   if(direct){try{const species=await this.activeApi.genome(direct,request.signal);if(request.epoch!==this.epoch)return;if(this.provider==='ensembl'&&!this.collectionAuto&&species.division!==this.division)throw Error('Species belongs to a different collection. Choose that collection in Advanced.');if(species.division)this.division=species.division as EnsemblDivision;this.species=[species];this.speciesName=species.name;this.speciesInput=species.display_name;resolved=true;}catch(e){if(!/Ensembl returned (400|404)/.test(String(e)))throw e;}}
   if(!resolved&&this.provider==='ncbi'){const species=await this.ncbi.genome(name,request.signal);if(request.epoch!==this.epoch)return;this.species=[species];this.speciesName=species.name;this.speciesInput=species.display_name;resolved=true;}
   if(!resolved){
   let list=await this.api.species(this.division,request.signal);if(request.epoch!==this.epoch)return;
   const matching=(items:LibrarySpecies[])=>{const exact=this.manualSpecies?items.find(s=>s.name===this.manualSpecies):preferred?items.find(s=>s.name===preferred[0]):items.find(s=>[s.name.replace(/_/g,' '),s.display_name.toLowerCase(),...(s.aliases||[]).map(a=>a.toLowerCase())].includes(name));return exact?[exact]:items.filter(s=>[s.name,s.display_name,s.strain,...(s.aliases||[])].some(v=>v?.toLowerCase().replace(/_/g,' ').includes(name)));};
   let matches=matching(list);
   if(!matches.length&&this.collectionAuto){const divisions:EnsemblDivision[]=['EnsemblPlants','EnsemblVertebrates','EnsemblMetazoa','EnsemblFungi','EnsemblProtists'];const others=divisions.filter(d=>d!==this.division),settled=await Promise.allSettled(others.map(async d=>(await this.api.species(d,request.signal)).map(s=>({...s,division:d}))));if(request.epoch!==this.epoch)return;const found=settled.flatMap(r=>r.status==='fulfilled'?r.value:[]);matches=matching(found);if(matches.length===1){this.division=matches[0].division as EnsemblDivision;list=found.filter(s=>s.division===this.division);}else list=[...list,...found];}
   this.species=list;
   if(matches.length!==1){this.speciesFilter=this.speciesInput;throw Error(matches.length?'Several species/cultivars match. Choose one in Advanced, or enter the full scientific name.':'Species not found. Use a scientific name or choose its collection in Advanced.');}
   this.speciesName=matches[0].name;this.speciesInput=matches[0].display_name;
   }
   const mode=this.searchType==='auto'?((this.provider==='ncbi'?/^[0-9\s,;]+$/.test(query):/[0-9]{3,}/.test(query))||query.split(/[,;]+/).length>1?'id':'symbol'):this.searchType;
   this.stage='Searching genes…';this.cd.markForCheck();
   let genes:LibraryGene[];
   try{genes=await this.activeApi.search(this.speciesName,query,mode,request.signal);}catch(e){if(this.searchType!=='auto'||mode!=='id'||request.signal.aborted||/limit|timed out|reach|cancel/i.test(String(e)))throw e;genes=await this.activeApi.search(this.speciesName,query,'symbol',request.signal);}
   if(request.epoch!==this.epoch)return;this.genes=genes;genes.forEach(g=>this.selected.add(g.id));if(!genes.length)this.error='No matching gene. Try its gene symbol or a stable gene ID.';
  }catch(e){if(request.epoch===this.epoch)this.error=e instanceof Error?e.message:String(e);}finally{this.done(request.epoch);}
 }
 speciesTextChanged(){this.manualSpecies='';this.inputsChanged();}
 chooseSimpleSpecies(){this.manualSpecies=this.speciesName;this.collectionAuto=false;this.speciesInput=this.choice?.display_name||this.speciesInput;this.inputsChanged();}
 async search(){if(!this.choice){this.error='Choose a species from the live catalogue first.';return;}const request=this.begin('Looking up gene identifiers…');this.genes=[];this.selected.clear();try{const genes=await this.api.search(this.speciesName,this.query,this.queryMode,request.signal);if(request.epoch!==this.epoch)return;this.genes=genes;genes.forEach(g=>this.selected.add(g.id));if(!genes.length)this.error='No matching gene in this species. Try its stable gene ID.';}catch(e){if(request.epoch===this.epoch)this.error=e instanceof Error?e.message:String(e);}finally{this.done(request.epoch);}}
 select(id:string){this.selected.has(id)?this.selected.delete(id):this.selected.add(id);}
 async importSelected(){const choice=this.choice,genes=this.genes.filter(g=>this.selected.has(g.id));if(!choice||!genes.length)return;const request=this.begin('Fetching assembly-verified genomic DNA…');const imported:LibraryReference[]=[];try{for(const [i,gene] of genes.entries()){this.stage=`Fetching ${gene.id} (${i+1}/${genes.length})…`;this.cd.markForCheck();const reference=await this.activeApi.importGene(gene,choice,this.division,this.upstream,this.downstream,request.signal);if(request.epoch!==this.epoch)return;imported.push(reference);}this.imported.emit(imported);this.open=false;}catch(e){if(request.epoch===this.epoch)this.error=e instanceof Error?e.message:String(e);}finally{this.done(request.epoch);}}
}
