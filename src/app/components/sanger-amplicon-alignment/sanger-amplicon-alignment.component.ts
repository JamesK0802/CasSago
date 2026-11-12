import {Component,Input,OnChanges} from '@angular/core';
import {CommonModule} from '@angular/common';
import {SangerCandidate,SangerDesignResult} from '../../models/sanger-design';
import {selectedAmplicons} from '../../utils/sanger-design-export';
import {reverseComplement} from '../../utils/biology.utils';
interface Segment {text:string;kind:string;position:number;cut:boolean}
interface AlignedProduct {name:string;length:number;reversed:boolean;cuts:number[];primerEnd:number;padding:number;guides:{start:number;end:number;name:string}[];segments:Segment[]}
@Component({selector:'app-sanger-amplicon-alignment',standalone:true,imports:[CommonModule],templateUrl:'./sanger-amplicon-alignment.component.html',styleUrl:'./sanger-amplicon-alignment.component.css'})
export class SangerAmpliconAlignmentComponent implements OnChanges {
 @Input() result:SangerDesignResult|null=null;
 @Input() selected:SangerCandidate|null=null;
 @Input() references:{name:string;targets:string;mode:string}[]=[];
 rows:AlignedProduct[]=[];width=0;ticks:{position:number;label:string}[]=[];
 ngOnChanges(){
  this.rows=[];this.width=0;this.ticks=[];
  if(!this.result||!this.selected?.partner)return;
  const c=this.selected,reverse=c.direction==='reverse';
  // Actual predicted PCR products, oriented to extension from the chosen sequencing primer.
  this.rows=selectedAmplicons(this.result,c).map(a=>{
   const sequence=reverse?reverseComplement(a.sequence):a.sequence;
   const cuts=a.cuts.map(cut=>reverse?sequence.length-cut:cut).sort((a,b)=>a-b);
   const primary=c.sequence.length,partner=(reverse?a.forwardPrimer:a.reversePrimer)?.length||c.partner!.sequence.length;
   const guides:AlignedProduct['guides']=[],raw=this.references.find(r=>r.name===a.sourceReference);
   if(raw?.mode==='guides')raw.targets.toUpperCase().split(/[\s,;]+/).filter(Boolean).forEach((g,i)=>{for(const query of new Set([g,reverseComplement(g)])){let at=sequence.indexOf(query);while(at>=0){if(cuts.some(cut=>at<=cut&&at+query.length>=cut))guides.push({start:at,end:at+query.length,name:a.targetNames?.[i]||'gRNA '+(i+1)});at=sequence.indexOf(query,at+1);}}});
   const segments:Segment[]=[];
   for(let i=0;i<sequence.length;i++){
    const kind=i<primary?'primer':i>=sequence.length-partner?'partner':guides.some(g=>i>=g.start&&i<g.end)?'guide':cuts.some(cut=>i>=cut-Math.floor(this.result!.settings.window/2)&&i<cut+Math.ceil(this.result!.settings.window/2))?'window':'';
    const cut=cuts.includes(i),last=segments.at(-1);
    if(last&&last.kind===kind&&!cut)last.text+=sequence[i];else segments.push({text:sequence[i],kind,cut,position:i});
   }
   return {name:a.sourceReference,length:sequence.length,reversed:reverse,cuts,primerEnd:primary,padding:0,guides,segments};
  });
  const anchor=Math.max(0,...this.rows.map(r=>r.primerEnd));
  this.rows.forEach(r=>r.padding=anchor-r.primerEnd);
  this.width=Math.max(0,...this.rows.map(r=>r.padding+r.length));
  this.ticks=Array.from({length:Math.floor(this.width/50)+1},(_,i)=>({position:i*50,label:String(i*50-anchor)}));
 }
}
