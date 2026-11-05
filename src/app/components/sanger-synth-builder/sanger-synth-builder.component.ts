import {initSynthProfiles,targetWT,SynthTargetProfile} from '../../utils/sanger-synth-profiles';
import {Component,Input,Output,EventEmitter,OnChanges,OnDestroy,ViewChild,ElementRef,ChangeDetectorRef} from '@angular/core';
import {CommonModule} from '@angular/common';
import {FormsModule} from '@angular/forms';
import {zipSync,strToU8} from 'fflate';
import {SangerCandidate,SangerDesignResult} from '../../models/sanger-design';
import {SynthReference,SynthAllele,SynthEvent,SynthDataset,synthReferences,synthWT,generateSangerSynth,synthTruthCsv,synthPickerContext,randomSynthEvent,randomSynthPercent,synthDescription} from '../../utils/sanger-synthetic-benchmark';
@Component({selector:'app-sanger-synth-builder',standalone:true,imports:[CommonModule,FormsModule],templateUrl:'./sanger-synth-builder.component.html',styleUrl:'./sanger-synth-builder.component.css'})
export class SangerSynthBuilderComponent implements OnChanges,OnDestroy {
 @Input() references:{name:string;targets:string;mode:string}[]=[];
 @Input() result:SangerDesignResult|null=null;@Input() selected:SangerCandidate|null=null;
 @Output() analyze=new EventEmitter<SynthDataset>();
 @ViewChild('picker') picker?:ElementRef<HTMLDialogElement>;
 refs:SynthReference[]=[];data:SynthDataset|null=null;error='';stage='';busy=false;
 pickerContext:ReturnType<typeof synthPickerContext>|null=null;
 picking:{ref:SynthReference;event:SynthEvent}|null=null;selection:{start:number;end:number}|null=null;
 private nextId=1;
 constructor(private cd:ChangeDetectorRef){}
 ngOnChanges(){if(this.picker?.nativeElement.open)this.picker.nativeElement.close();this.picking=null;this.selection=null;this.stop();this.refs=this.result&&this.selected?.partner?synthReferences(this.result,this.selected):[];this.refs.forEach(initSynthProfiles);this.invalidate();}
 ngOnDestroy(){this.stop();}
 referenceName(name:string){return this.refs.find(r=>r.name===name)?.sourceReference||this.data?.truth.find(t=>t.reference===name)?.sourceReference||name;}
 wt(profile:SynthTargetProfile){return targetWT(profile);}
 get shareTotal(){return this.refs.reduce((s,r)=>s+Number(r.share),0);}
 get combinationRows(){return (this.data?.truth||[]).map(row=>{const ref=this.refs.find(r=>r.name===row.reference)!;return {reference:row.sourceReference,percent:row.percent,samplePercent:row.totalPercent,targets:(ref.profiles||[]).map(profile=>{const events=row.events.filter(e=>e.target===profile.target),index=profile.alleles.findIndex(g=>JSON.stringify(g.events)===JSON.stringify(events));return {group:events.length?'G'+(row.targetStates?.[profile.target]??index+1):'WT',description:events.map(e=>e.kind==='substitution'?e.before+'→'+e.sequence:(e.offset>=0?'+':'')+e.offset+' bp · '+(e.kind==='insertion'?'+'+e.sequence:'−'+e.length+' bp')).join('; ')};})};});}
 get profileTableTargets(){return this.refs[0]?.profiles?.map(p=>p.name)||[];}
 get profileRows(){return this.refs.flatMap(ref=>(ref.profiles||[]).flatMap(profile=>[{reference:this.referenceName(ref.name),target:profile.name,group:'WT',percent:this.wt(profile),samplePercent:ref.share*this.wt(profile)/100,description:'WT'},...profile.alleles.map((a,i)=>({reference:this.referenceName(ref.name),target:profile.name,group:'Group '+(i+1),percent:a.percent,samplePercent:ref.share*a.percent/100,description:synthDescription(a.events)}))]));}
 invalidate(){this.stop();this.data=null;this.error='';this.stage='';}
 addAllele(ref:SynthReference,profile:SynthTargetProfile){try{const event=randomSynthEvent(ref,Math.random,profile.target);profile.alleles.push({id:this.nextId++,percent:randomSynthPercent(profile.alleles.map(a=>a.percent),Math.random,this.refs.flatMap(r=>(r.profiles||[]).flatMap(p=>p.alleles.map(a=>a.percent)))),events:[event]});this.invalidate();}catch(e){this.error=e instanceof Error?e.message:String(e);}}
 removeAllele(profile:SynthTargetProfile,allele:SynthAllele){profile.alleles=profile.alleles.filter(a=>a!==allele);this.invalidate();}
 newEvent():SynthEvent{return {kind:'deletion',target:0,offset:0,length:1,sequence:''};}
 addEvent(allele:SynthAllele){allele.events.push(this.newEvent());this.invalidate();}
 removeEvent(allele:SynthAllele,event:SynthEvent){allele.events=allele.events.filter(e=>e!==event);this.invalidate();}
 eventType(event:SynthEvent){delete event.selectionStart;delete event.before;event.sequence='';this.invalidate();}
 validSub(event:SynthEvent){return !!event.before&&/^[ACGT]+$/.test(event.sequence.toUpperCase())&&event.sequence.length===event.before.length&&event.sequence.toUpperCase()!==event.before;}
 openPicker(ref:SynthReference,event:SynthEvent){this.picking={ref,event};this.selection=null;const raw=this.references.find(r=>r.name===ref.sourceReference),guides=raw?.mode==='guides'?raw.targets.split(/[\s,;]+/).filter(Boolean):[];this.pickerContext=synthPickerContext(ref,event.target,guides);this.picker?.nativeElement.showModal();}
 captureSelection(root:HTMLElement){const selection=window.getSelection();if(!selection?.rangeCount)return;const range=selection.getRangeAt(0),first=(range.startContainer.nodeType===Node.TEXT_NODE?range.startContainer.parentElement:range.startContainer as HTMLElement)?.closest<HTMLElement>('[data-base]'),last=(range.endContainer.nodeType===Node.TEXT_NODE?range.endContainer.parentElement:range.endContainer as HTMLElement)?.closest<HTMLElement>('[data-base]');if(!first||!last||!root.contains(first)||!root.contains(last))return;const start=Number(first.dataset['base'])+Number(range.startOffset>0),end=Number(last.dataset['base'])+Number(range.endOffset>0);this.selection=end>start?{start,end}:null;}
 useSelection(){if(!this.picking||!this.selection)return;const{ref,event}=this.picking;event.selectionStart=this.selection.start;event.before=ref.sequence.slice(this.selection.start,this.selection.end);event.offset=this.selection.start-ref.cuts[event.target];this.invalidate();this.picker?.nativeElement.close();}
 generate(){this.error='';try{if(!this.result||!this.selected)throw Error('Select a primer pair.');this.data=generateSangerSynth(this.result,this.selected,this.refs);}catch(e){this.data=null;this.error=e instanceof Error?e.message:String(e);}}
 downloadFiles(){if(!this.data)return;const d=this.data,truth={format:'casmango-sanger-synth-v1',signal:'Ideal Gaussian peaks; no added noise or instrument response',cycles:d.cycles,recipe:d.recipe,mixing:d.mixing,truth:d.truth};this.save('CasSago_SYNTHETIC.zip',zipSync({'SYNTHETIC_control.ab1':new Uint8Array(d.control),'SYNTHETIC_edit.ab1':new Uint8Array(d.edit),'truth.csv':strToU8(synthTruthCsv(d)),'target-profiles.csv':strToU8([['Reference','Target','Group','Within-target %','Sample %','Edit'],...this.profileRows.map(r=>[r.reference,r.target,r.group,r.percent,r.samplePercent,r.description])].map(row=>row.map(v=>'"'+String(v).replace(/"/g,'""')+'"').join(',')).join('\r\n')),'truth.json':strToU8(JSON.stringify(truth,null,2)),'analysis-template.json':strToU8(JSON.stringify(d.template,null,2))}),'application/zip');}
 stop(){this.busy=false;}
 private save(name:string,data:string|Uint8Array,type:string){const blob=new Blob([typeof data==='string'?data:new Uint8Array(data).buffer],{type}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),30000);}
}
