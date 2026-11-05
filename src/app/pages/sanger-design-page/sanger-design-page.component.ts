import {SangerSynthBuilderComponent} from '../../components/sanger-synth-builder/sanger-synth-builder.component';
import {SynthDataset} from '../../utils/sanger-synthetic-benchmark';
import {SangerPrimerDifferenceComponent} from '../../components/sanger-primer-difference/sanger-primer-difference.component';
import {SangerAmpliconAlignmentComponent} from '../../components/sanger-amplicon-alignment/sanger-amplicon-alignment.component';
import {SangerPrimerMapComponent} from '../../components/sanger-primer-map/sanger-primer-map.component';
import {selectedAmplicons,makeAnalysisTemplate,DesignedAmplicon} from '../../utils/sanger-design-export';
import {ReferenceSourceComponent} from '../../components/reference-library/reference-source.component';
import {LibraryReference,ReferenceSource} from '../../models/reference-source';
import { AfterViewChecked, ChangeDetectorRef, Component, ElementRef, OnDestroy, Optional, ViewChild } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { SANGER_DEFAULTS, SangerCandidate, SangerDesignResult, SangerDesignSettings } from '../../models/sanger-design';
import { parseSangerFasta, resolveSangerReference, sangerDesignExample } from '../../utils/sanger-design.utils';
import { validateDesign } from '../../workers/core/sanger-primer-designer';
import { SangerAnalysisPanelComponent } from '../sanger-analysis-panel/sanger-analysis-panel.component';
import { DesignHandoffService } from '../../services/design-handoff.service';

interface NamedTarget {name:string;guide:string;cut?:number}
interface ReferenceRow { name: string; sequence: string; targets: string; mode: 'cuts' | 'guides';source?:ReferenceSource;targetRows?:NamedTarget[] }
@Component({
  selector:'app-sanger-design-page', standalone:true, imports:[CommonModule,FormsModule,SangerAnalysisPanelComponent,SangerPrimerMapComponent,SangerAmpliconAlignmentComponent,SangerPrimerDifferenceComponent,SangerSynthBuilderComponent,ReferenceSourceComponent],
  templateUrl:'./sanger-design-page.component.html', styleUrls:['./sanger-design-page.component.css','../sanger-ui.css'],
})
export class SangerDesignPageComponent implements OnDestroy, AfterViewChecked {
  @ViewChild(SangerAnalysisPanelComponent) analysisPanel?:SangerAnalysisPanelComponent;
  visibleCount=10;
  showMap=false;
  showSynth=false;
  showChroma=false;
  @ViewChild('sequenceScroll') sequenceScroll?: ElementRef<HTMLElement>;
  private sequenceCenterPending = false;
  activeSection:'design'|'analysis'='design';
  private defaults():SangerDesignSettings{return {...SANGER_DEFAULTS,platform:'sanger',direction:'both',readEnd:1000,requirePartner:true,primerMin:20,primerMax:20};}
  settings: SangerDesignSettings = this.defaults();
  amplification:'shared'|'specific'='shared';
  rows: ReferenceRow[] = [this.blank(1),this.blank(2)];
  result: SangerDesignResult | null = null;
  selected: SangerCandidate | null = null;
  mapTarget = 0; sequenceFocus: 'target' | 'primer' | 'partner' = 'target';
  busy = false; error = ''; progress = 0; stage = ''; notice = '';
  private worker: Worker | null = null;
  readonly colours = ['#008c71','#347cb6','#c28a16','#a357a0','#df7750','#727cc0','#708f45','#8e6e5e'];
  constructor(private changeDetector: ChangeDetectorRef, @Optional() handoff?:DesignHandoffService) {
    try {sessionStorage.removeItem('casmango.sanger-design.draft.v1');sessionStorage.removeItem('casmango.sanger-design.draft.v2');} catch { /* Inputs stay in memory only. */ }
    const incoming=handoff?.take();
    if(incoming){this.rows=incoming.references.map(r=>({...r,targets:r.cuts.join(', '),mode:'cuts'}));this.settings={...this.defaults(),workflow:incoming.workflow};this.amplification=incoming.workflow==='homeolog-pcr'?'specific':'shared';this.rows.forEach(row=>{if(/^[ACGT]{20}$/i.test(incoming.guide))row.targetRows=[{name:'Target 1',guide:incoming.guide,cut:Number(row.targets.split(',')[0])}];});this.changed();this.notice='Selected gRNA copied: '+incoming.guide+'. Set the required window and check primer coverage.';}
    this.rows.forEach(row=>this.targetRows(row));this.normalizeSettings();
  }
  blank(n: number): ReferenceRow { return {name:`Gene ${n}`,sequence:'',targets:'',mode:'guides',targetRows:[{name:'Target 1',guide:''}]}; }
  targetRows(row:ReferenceRow):NamedTarget[]{
    if(!row.targetRows?.length){const tokens=row.targets.trim().split(/[\s,;]+/).filter(Boolean);row.targetRows=tokens.length?tokens.map((value,i)=>({name:'Target '+(i+1),guide:row.mode==='guides'?value:'',...(row.mode==='cuts'?{cut:Number(value)}:{})})):[{name:'Target 1',guide:''}];}
    return row.targetRows;
  }
  addTarget(row:ReferenceRow){if(this.targetRows(row).length>=2)return;const names=new Set(row.targetRows!.map(t=>t.name));let n=1;while(names.has('Target '+n))n++;row.targetRows!.push({name:'Target '+n,guide:''});this.changed();}
  removeTarget(row:ReferenceRow,index:number){if(this.targetRows(row).length<=1)return;row.targetRows!.splice(index,1);this.changed();}
  referenceChanged(row:ReferenceRow){for(const target of this.targetRows(row))if(target.guide.trim())delete target.cut;this.changed();}
  targetChanged(row:ReferenceRow,target:NamedTarget){delete target.cut;this.changed();}
  private normalizeSettings(){const fixed=['productMin','productMax','sodiumMm','primerNm','maxTmDifference','specificityMismatches','specificBoth','specificityStrategy','threePrimeWindow'] as const;for(const key of fixed)(this.settings as any)[key]=SANGER_DEFAULTS[key];this.settings.platform='sanger';this.settings.direction='both';this.settings.requirePartner=true;this.settings.indelRange=SANGER_DEFAULTS.indelRange;this.settings.readEnd=Math.min(1000,this.settings.readEnd);this.settings.workflow=this.rows.length===1?'standard-pcr':this.amplification==='specific'?'homeolog-pcr':'shared-sanger';}
  private syncTargets(row:ReferenceRow){const targets=this.targetRows(row);row.mode=targets.every(t=>!t.guide.trim()&&Number.isInteger(t.cut))?'cuts':'guides';row.targets=targets.map(t=>t.guide.trim()||String(t.cut??'')).join(', ');}
  private designReference(row:ReferenceRow){
    const targets=this.targetRows(row);if(targets.length<1||targets.length>2)throw Error(row.name+': use one or two targets.');
    if(targets.some(t=>!t.name.trim())||new Set(targets.map(t=>t.name.trim())).size!==targets.length)throw Error(row.name+': give each target a unique name.');
    const cuts=targets.map(t=>resolveSangerReference(row.name,row.sequence,Number.isInteger(t.cut)?String(t.cut):t.guide.trim(),Number.isInteger(t.cut)?'cuts':'guides').cuts[0]);
    if(new Set(cuts).size!==cuts.length)throw Error(row.name+': targets resolve to the same cut position.');
    return {...resolveSangerReference(row.name,row.sequence,cuts.join(', '),'cuts'),source:row.source,targetNames:targets.map(t=>t.name.trim()),guides:targets.map(t=>t.guide.trim())};
  }
  get analysisDesignRows(){return this.rows.map(row=>{try{return {...row,targets:this.designReference(row).cuts.join(', '),mode:'cuts' as const};}catch{return row;}});}
  changed() {
    this.rows.forEach(row=>this.syncTargets(row));this.normalizeSettings();
    this.result = null; this.selected = null;this.visibleCount=10; this.mapTarget = 0; this.error = ''; this.notice = '';
  }
  addReference() { if (this.rows.length < 8) { this.rows.push(this.blank(this.rows.length+1)); this.changed(); } }
  removeReference(i: number) { this.rows.splice(i,1); this.changed(); }
  loadExample() {
    const guides=['GACCTAGTGCATGACCTGAC','TGCAGATCGTACCTAGACGA'];
    this.rows = sangerDesignExample(2).map(r => {
      let sequence=r.sequence;
      r.cuts.forEach((cut,t)=>{sequence=sequence.slice(0,cut-17)+guides[t]+'AGG'+sequence.slice(cut+6);});
      return {name:r.name,sequence,targets:guides.join(', '),mode:'guides',targetRows:guides.map((guide,t)=>({name:'Target '+(t+1),guide}))};
    });
    this.settings = this.defaults();this.amplification='shared';this.changed();
    this.notice = 'Constructed example · shared binding motif with different target distances. Not experimental data.';
  }
  importLibrary(imported:LibraryReference[]){if(this.busy){this.error='Finish or cancel the current search before importing.';return;}const existing=this.rows.filter(r=>r.sequence.trim());if(existing.length+imported.length>8){this.error='Keep at most eight references. Remove an input before importing.';return;}if(imported.some(r=>existing.some(e=>e.name===r.name||(e.source?.geneId===r.source.geneId&&e.source?.species===r.source.species)))){this.error='This gene is already present. Use Check latest on that reference.';return;}this.rows=[...existing,...imported.map(r=>({...r,targets:'',mode:'cuts' as const}))];if(this.rows.length===1)this.settings.workflow='standard-pcr';else if(this.settings.workflow==='standard-pcr')this.settings.workflow='shared-sanger';this.changed();this.notice='Genomic DNA imported. Add one or two named gRNAs before primer design.';}
  updateLibrary(i:number,next:LibraryReference){if(this.busy)return;this.rows[i]={...this.rows[i],...next,name:this.rows[i].name,targets:'',mode:'guides',targetRows:[{name:'Target 1',guide:''}]};this.changed();this.notice='Updated reference applied. Its target coordinates and previous primer results were cleared.';}
  async upload(event: Event, type: 'fasta' | 'project') {
    const input = event.target as HTMLInputElement, file = input.files?.[0]; if (!file) return;
    try {
      if (file.size > 20000000) throw new Error('Use an input smaller than 20 MB.');
      const text = await file.text();
      if (type === 'fasta') {
        const refs = parseSangerFasta(text);
        if (refs.length < 1 || refs.length > 8) throw new Error('Import 1–8 FASTA references.');
        if(refs.length===1)this.settings.workflow='standard-pcr';else if(this.settings.workflow==='standard-pcr')this.settings.workflow='shared-sanger';
        this.rows = refs.map(r => ({...r,targets:'',mode:'cuts'}));this.changed();
      } else {
        const parsed = JSON.parse(text);
        if (parsed.format !== 'casmango-sanger-design-v1' || !Array.isArray(parsed.references)) throw new Error('Choose a CasSago design JSON export.');
        const settings = {...SANGER_DEFAULTS,...parsed.settings};validateDesign(parsed.references,settings);
        this.rows = parsed.references.map((r:any) => ({name:r.name,sequence:r.sequence,targets:r.cuts.join(', '),mode:'cuts',source:r.source}));
        this.settings = settings;this.changed();
      }
    } catch(e) { this.error = e instanceof Error ? e.message : String(e); }
    input.value = '';this.changeDetector.markForCheck();
  }
  search() {
    this.error = ''; this.notice = ''; this.result = null; this.selected = null;
    try {
      this.normalizeSettings();const references = this.rows.map(r=>this.designReference(r));
      validateDesign(references,this.settings);
      if (typeof Worker === 'undefined') throw new Error('This browser does not support background workers. Open the page in a current browser.');
      this.busy = true;this.progress = 0;this.stage = 'Preparing reference search';
      this.worker = new Worker(new URL('../../workers/sanger-design.worker',import.meta.url),{type:'module'});
      this.worker.onmessage = ({data}) => {
        if (data.type === 'progress') {this.progress = Math.round(data.fraction*100);this.stage = data.stage;}
        else if (data.type === 'result') {
          this.result = data.result;this.visibleCount=10;this.selected = this.result?.candidates[0] || null;this.mapTarget=0;this.sequenceFocus='target';this.centerSequence();this.stopWorker();
        } else {this.error = data.message;this.stopWorker();}
        this.changeDetector.markForCheck();
      };
      this.worker.onerror = () => {this.error = 'The design worker could not complete. Reduce reference length and try again.';this.stopWorker();this.changeDetector.markForCheck();};
      this.worker.postMessage({references,settings:{...this.settings}});
    } catch(e) {this.error = e instanceof Error ? e.message : String(e);this.stopWorker();}
  }
  cancel() {this.stopWorker();this.stage = '';this.notice = 'Search cancelled.';}
  private stopWorker() {this.worker?.terminate();this.worker = null;this.busy = false;}
  ngOnDestroy() {this.stopWorker();}
  mismatchLabel(screen?:{minimumMismatches:number;minimumExact?:boolean} | null){return !screen?'—':screen.minimumExact===false?'Not resolved':String(screen.minimumMismatches);}
  get primerScreens(){const c=this.selected;return c?.primerSafety?[{label:'Sequencing primer',screen:c.primerSafety.primary},...(c.primerSafety.partner?[{label:'Opposite PCR primer',screen:c.primerSafety.partner}]:[])]:[];}
  get visibleCandidates(){return this.result?.candidates.slice(0,this.visibleCount)||[];}
  loadMore(){this.visibleCount+=10;}
  setPlatform(value:'sanger'|'ngs'){this.settings.platform=value;if(value==='ngs'&&this.settings.workflow==='shared-sanger')this.settings.workflow='shared-pcr';this.changed();}
  setAmplification(value:'shared'|'specific'|'single'){this.amplification=value==='specific'?'specific':'shared';this.changed();}
  get sharedMode(){return this.settings.workflow==='shared-pcr'||this.settings.workflow==='shared-sanger';}
  setOffsetDesign(enabled:boolean){this.settings.workflow=enabled?'shared-sanger':'shared-pcr';this.changed();}
  get analysisAmplicons():DesignedAmplicon[]{try{return this.result&&this.selected?selectedAmplicons(this.result,this.selected):[];}catch{return [];}}
  get templateIssue():string{try{if(this.result&&this.selected)makeAnalysisTemplate(this.result,this.selected);return '';}catch(e){return e instanceof Error?e.message:String(e);}}
  useGeneratedSynth(data:SynthDataset){try{if(!this.analysisPanel)return;this.analysisPanel.applySynthetic(data.template,data.control,data.edit,{recipe:data.recipe,truth:data.truth,mixing:data.mixing,cycles:data.cycles});this.activeSection='analysis';}catch(e){this.error=e instanceof Error?e.message:String(e);}}
  useInAnalysis(){try{if(!this.result||!this.selected||!this.analysisPanel)return;this.analysisPanel.applyTemplate(makeAnalysisTemplate(this.result,this.selected));this.activeSection='analysis';}catch(e){this.error=e instanceof Error?e.message:String(e);}}
  exportAnalysisTemplate(){try{if(this.result&&this.selected)this.download('CasSago_profile.json',JSON.stringify(makeAnalysisTemplate(this.result,this.selected),null,2),'application/json');}catch(e){this.error=e instanceof Error?e.message:String(e);}}
  exportTargets(){const rows=[['Reference','Target','Cut position (after base N)'],...this.analysisAmplicons.flatMap(a=>a.cuts.map((cut,i)=>[a.name,a.targetNames?.[i]||'T'+(i+1),String(cut)]))];this.download('CasSago_selected_targets.csv',rows.map(row=>row.map(v=>'"'+v.replace(/"/g,'""')+'"').join(',')).join('\r\n'),'text/csv');}
  get status(): string {
    if (!this.result) return '';
    if (this.result.settings.workflow === 'homeolog-pcr') {
      const covered=this.result.coverage?.filter(r=>r.candidates>0).length || 0;
      return `${this.result.candidates.length} pairs · ${covered}/${this.result.references.length} references`;
    }
    if (this.result.candidates.length) return `${this.result.candidates.length} candidate${this.result.candidates.length === 1 ? '' : 's'} found`;
    return this.result.complete ? 'No candidate under these settings' : 'No candidate in the checked shortlist';
  }
  get noResultReason(): string {
    const st = this.result?.stats;if (!st) return '';
    if(this.result!.settings.workflow==='standard-pcr')return 'No PCR pair passed the required window, usable read range and oligo/product filters.';
    if (this.result!.settings.workflow === 'homeolog-pcr') return 'No separate PCR pair passed the window, product and selected mismatch-specificity screen.';
    if(this.result!.settings.workflow==='shared-pcr')return 'No common PCR pair passed the required window, platform and oligo/product filters.';
    if (!st.commonUnique) return 'No quality-screened primer has exactly one exact site per reference in a shared orientation.';
    if (!st.geometryPassed) return 'Shared primers exist, but none covers every required target window within the usable Sanger range.';
    if (!st.offsetPassed) return 'Shared primers cover the targets, but their target distances do not meet the minimum separation.';
    if(st.separabilityRejected)return 'Primer/product candidates were excluded because their observable WT signals cannot be distinguished. Change the primer position, reference span or usable read range.';
    return 'Offset primers exist, but no checked candidate has a qualifying common PCR partner. Try sequencing-primer-only mode or adjust PCR limits.';
  }
  selectCandidate(c: SangerCandidate) {this.selected=c;this.mapTarget=0;this.sequenceFocus='target';this.centerSequence();}
  toggleCandidateMap(c:SangerCandidate){const same=this.selected?.id===c.id;this.selectCandidate(c);this.showMap=same?!this.showMap:true;}
  toggleCandidateTrace(c:SangerCandidate){this.showMap=true;const same=this.selected?.id===c.id;this.selectCandidate(c);this.showChroma=same?!this.showChroma:true;}
  centerSequence() {this.sequenceCenterPending=true;}
  ngAfterViewChecked() {
    if (!this.sequenceCenterPending || !this.sequenceScroll || !this.selected) return;
    this.sequenceCenterPending=false;
    const el=this.sequenceScroll.nativeElement, bounds=this.sequenceBounds(0);
    el.scrollLeft=Math.max(0,(bounds[1]-bounds[0])*13/2-(el.clientWidth-106)/2);
  }

  refIndex(i:number): number {return this.selected?.referenceIndex ?? i;}
  reference(i:number) {return this.result!.references[this.refIndex(i)];}
  get specificMode(): boolean {return (this.result?.settings || this.settings).workflow==='homeolog-pcr';}
  // Target-relative coordinates make different genomic FASTA origins irrelevant.
  origin(i:number):number {return this.reference(i).cuts[this.mapTarget] ?? this.reference(i).cuts[0];}
  mapBounds():[number,number] {
    const values:number[]=[];
    this.selected?.sites.forEach((hit,i)=>{
      const o=this.origin(i);values.push(hit.start-o,hit.end-o);
      for (const cut of this.reference(i).cuts) values.push(cut-o-this.result!.settings.window/2,cut-o+this.result!.settings.window/2);
      const p=this.selected?.partner?.sites[i];if(p) values.push(p.start-o,p.end-o);
    });
    return [Math.min(...values,0)-50,Math.max(...values,0)+50];
  }
  mapX(bp:number,i:number):number {const [lo,hi]=this.mapBounds();return 100*(bp-this.origin(i)-lo)/(hi-lo);}
  mapWidth(bp:number):number {const [lo,hi]=this.mapBounds();return 100*bp/(hi-lo);}
  mapTick(which:'start'|'end'):number {return Math.round(this.mapBounds()[which==='start'?0:1]);}
  windowStart(i:number,t:number):number {return this.reference(i).cuts[t]-Math.floor(this.result!.settings.window/2);}
  sequenceBounds(i:number):[number,number] {
    const ref=this.reference(i),hit=this.selected!.sites[i],partner=this.selected!.partner?.sites[i];
    const center=this.sequenceFocus==='primer'?(hit.start+hit.end)/2:this.sequenceFocus==='partner'&&partner?(partner.start+partner.end)/2:this.origin(i);
    const width=this.sequenceFocus==='target'?this.result!.settings.window+80:100;
    const start=Math.floor(center-width/2);return [start,start+width];
  }
  sequenceCells(i:number):{base:string;position:number;kind:string;cut:boolean}[] {
    const [lo,hi]=this.sequenceBounds(i),ref=this.reference(i),hit=this.selected!.sites[i],p=this.selected!.partner?.sites[i];
    return Array.from({length:hi-lo},(_,j)=>{
      const position=lo+j,base=ref.sequence[position] || '·';
      const kind=position>=hit.start&&position<hit.end?'primer':p&&position>=p.start&&position<p.end?'partner':ref.cuts.some(c=>position>=c-Math.floor(this.result!.settings.window/2)&&position<c+Math.ceil(this.result!.settings.window/2))?'window':'';
      return {base,position,kind,cut:ref.cuts.includes(position)};
    });
  }
  sequenceRangeLabel(i:number):string {const [lo,hi]=this.sequenceBounds(i);return `${Math.max(1,lo+1).toLocaleString()}–${Math.min(this.reference(i).sequence.length,hi).toLocaleString()}`;}
  trackBase(_:number,cell:{position:number}) {return cell.position;}
  extensionPreview(i: number): string {return this.selected!.sites[i].extension.slice(this.result!.settings.readStart,this.result!.settings.readStart+90);}
  exportJson() {
    if (!this.result) return;
    this.download('CasSago_design.json',JSON.stringify({format:'casmango-sanger-design-v1',createdAt:new Date().toISOString(),
      scope:'Input-reference primer screening only; no genome-wide specificity or trace-deconvolution validation.',...this.result},null,2),'application/json');
  }
  exportCsv() {
    if (!this.result) return;
    const rows: (string|number)[][] = [['Candidate','Sequencing primer (5 to 3)','Direction','Tm Na-only (C)','GC (%)','Reference','Binding start (1-based)','Binding end (inclusive)','Target offsets after primer (bp)','Min target separation (bp)','WT trace differences','PCR partner (5 to 3)','Product length (bp)','Workflow','Primary specific','Partner specific','Mismatch screen','Specificity strategy','3-prime variants','Status','Farthest cut after sequencing primer (bp)','Primary minimum unintended-site mismatches','Primary sites with <=3 mismatches','Partner minimum unintended-site mismatches','Partner sites with <=3 mismatches']];
    for (const c of this.result.candidates) c.sites.forEach((hit,i) => rows.push([c.id,c.sequence,c.direction,c.tm.toFixed(2),c.gc.toFixed(1),hit.reference,hit.start+1,hit.end,hit.cutOffsets.join('; '),c.referenceIndex !== undefined ? '' : c.minOffset,c.referenceIndex !== undefined ? '' : c.minTraceDifferences,c.partner?.sequence || '',c.partner?.products[i] || '',this.result!.settings.workflow,c.specificity ? String(c.specificity.primary) : '',c.specificity ? String(c.specificity.partner) : '',c.specificity?.mismatches ?? '',c.specificity?.strategy || '',[...(c.specificity?.primaryVariants||[]),...(c.specificity?.partnerVariants||[])].map(v=>`${v.background}:3prime-${v.fromThreePrime-1}:${v.targetBase}/${v.backgroundBase}@${v.position}`).join('; '), 'Experimental candidate; input-reference screen only',c.targetDistance??Math.max(...c.sites.flatMap(h=>h.cutOffsets)),this.mismatchLabel(c.primerSafety?.primary),c.primerSafety?.primary.nearSiteCount??'',this.mismatchLabel(c.primerSafety?.partner),c.primerSafety?.partner?.nearSiteCount??'']));
    this.download('CasSago_primers.csv',rows.map(row => row.map(v => `"${String(v).replace(/"/g,'""')}"`).join(',')).join('\r\n'),'text/csv');
  }
  exportSelected() {
    const c = this.selected;if (!c) return;
    this.download('CasSago_selected_primers.fasta',`>${this.specificMode || this.result!.settings.workflow==='standard-pcr' ? this.reference(0).name.replace(/[^a-zA-Z0-9_-]/g,'_') : 'Shared'}_sequencing_primer_${c.direction}\n${c.sequence}\n`+(c.partner ? `>${this.specificMode || this.result!.settings.workflow==='standard-pcr' ? 'Opposite' : 'Common'}_PCR_partner_not_for_Sanger_mixture\n${c.partner.sequence}\n` : ''),'text/plain');
  }
  exportAmplicons() {
    const c=this.selected;if (!c?.partner) return;
    const records=this.analysisAmplicons.map(a=>`>${a.name}\n${a.sequence}`);
    this.download('CasSago_selected_amplicons.fasta',records.join('\n')+'\n','text/plain');
  }
  private download(name: string,text: string,type: string) {
    const url = URL.createObjectURL(new Blob([text],{type}));const a = document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(() => URL.revokeObjectURL(url),1000);
  }
}
