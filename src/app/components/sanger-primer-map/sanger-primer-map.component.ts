import {SangerWTPreviewComponent} from '../sanger-wt-preview/sanger-wt-preview.component';
import {SangerPrimerDifferenceComponent} from '../sanger-primer-difference/sanger-primer-difference.component';
import {Component,Input,Output,EventEmitter,OnChanges,HostListener} from '@angular/core';
import {CommonModule} from '@angular/common';
import {oppositePrimer,SangerCandidate,SangerDesignResult,PrimerSite} from '../../models/sanger-design';
import {reverseComplement} from '../../utils/biology.utils';
import {oppositePrimerColour} from '../../utils/sanger-primer-colours';
interface Track {name:string;primary:PrimerSite;partner?:PrimerSite;length:number;cuts:{name:string;offset:number}[];guides:{start:number;end:number}[]}
interface CandidateMap {candidate:SangerCandidate;index:number;tracks:Track[];low:number;high:number;ticks:number[]}
@Component({selector:'app-sanger-primer-map',standalone:true,imports:[CommonModule,SangerPrimerDifferenceComponent,SangerWTPreviewComponent],templateUrl:'./sanger-primer-map.component.html',styleUrl:'./sanger-primer-map.component.css'})
export class SangerPrimerMapComponent implements OnChanges {
 @Input() compact=false;
 @Input() showOpposite=true;
 @Input() result:SangerDesignResult|null=null;
 @Input() candidates:SangerCandidate[]=[];
 @Input() selected:SangerCandidate|null=null;
 @Input() references:{name:string;targets:string;mode:string}[]=[];
 @Output() candidateSelected=new EventEmitter<SangerCandidate>();
 readonly Math=Math;
 @Input() showChroma=false;
 @Output() chromaRequested=new EventEmitter<SangerCandidate>();
 private lastKey='';private lastResult:SangerDesignResult|null=null;
 groups:CandidateMap[]=[];
 popup:{id:string;label:string;sequence:string;location?:string}|null=null;
 hover:{id:string;reference:string;kind:'primary'|'partner';sequence:string;cells:{base:string;mismatch:boolean}[];left:number;top:number}|null=null;
 ngOnChanges(){
  const key=this.candidates.map(c=>c.id).join('|')+';'+this.selected?.id+';'+this.showOpposite+';'+this.references.map(r=>r.name+':'+r.targets+':'+r.mode).join('|');if(this.lastResult===this.result&&this.lastKey===key)return;this.lastResult=this.result;this.lastKey=key;
  this.hover=null;
  if(!this.result){this.groups=[];return;}
  this.groups=this.candidates.map((candidate,index)=>{
   const tracks=candidate.sites.map(primary=>{
    const ref=this.result!.references.find(r=>r.name===primary.reference)!,partner=this.showOpposite?candidate.partner?.sites.find(p=>p.reference===primary.reference):undefined;
    const raw=this.references.find(r=>r.name===ref.name),guides:Track['guides']=[];
    if(raw?.mode==='guides')for(const g of raw.targets.toUpperCase().split(/[\s,;]+/).filter(Boolean))for(const seq of new Set([g,reverseComplement(g)])){
     let at=ref.sequence.indexOf(seq);while(at>=0){if(ref.cuts.some(c=>at<=c&&at+seq.length>=c))guides.push({start:this.offset(at,primary),end:this.offset(at+seq.length,primary)});at=ref.sequence.indexOf(seq,at+1);}
    }
    return {name:ref.name,primary,partner,length:partner?Math.max(primary.end,partner.end)-Math.min(primary.start,partner.start):Math.min(this.result!.settings.readEnd,primary.direction==='forward'?ref.sequence.length-primary.end:primary.start)+candidate.sequence.length,cuts:ref.cuts.map((c,i)=>({name:ref.targetNames?.[i]||'T'+(i+1),offset:this.offset(c,primary)})),guides};
   });
   const low=-candidate.sequence.length,high=Math.max(1,...tracks.map(t=>t.length-candidate.sequence.length));
   const step=high>1000?250:100,ticks=Array.from({length:Math.floor(high/step)+1},(_,i)=>i*step);
   const originalIndex=this.result!.candidates.findIndex(c=>c.id===candidate.id);return {candidate,index:originalIndex<0?index:originalIndex,tracks,low,high,ticks};
  });
  if(this.popup&&(!this.groups.some(g=>g.candidate.id===this.popup!.id)||this.popup.id!==this.selected?.id))this.popup=null;
 }
 offset(bp:number,primary:PrimerSite){return primary.direction==='forward'?bp-primary.end:primary.start-bp;}
 choose(candidate:SangerCandidate){this.popup=null;this.candidateSelected.emit(candidate);}
 toggleTrace(candidate:SangerCandidate){this.chromaRequested.emit(candidate);}
 trackGroup(_:number,group:CandidateMap){return group.candidate.id;}
 trackReference(_:number,track:Track){return track.name;}
 @HostListener('document:pointerdown',['$event']) dismissPopup(event:PointerEvent){const el=event.target as Element|null;if(!el?.closest('.binding,.primer-popup'))this.popup=null;}
 showPrimer(group:CandidateMap,kind:'primary'|'partner',track?:Track){const seq=kind==='primary'?group.candidate.sequence:oppositePrimer(group.candidate,track?group.candidate.sites.indexOf(track.primary):0)?.sequence;if(!seq)return;this.choose(group.candidate);this.popup={id:group.candidate.id,label:kind==='primary'?'Sequencing primer · 5′ → 3′':'Opposite PCR primer · 5′ → 3′',sequence:seq,...(track?{location:track.name+' · '+((kind==='primary'?track.primary:track.partner)!.start+1)+'–'+(kind==='primary'?track.primary:track.partner)!.end+' · '+(kind==='primary'?track.primary:track.partner)!.direction}:{})};}
 primerCells(group:CandidateMap,track:Track,kind:'primary'|'partner'){
  const sequence=kind==='primary'?group.candidate.sequence:this.oppositeSequence(group,track),positions=kind==='primary'?track.primary.mismatchPositions||[]:[];
  return [...sequence].map((base,i)=>({base,mismatch:positions.includes(i+1)}));
 }
 hoverPrimer(group:CandidateMap,track:Track,kind:'primary'|'partner',index:number){
  const sequence=kind==='primary'?group.candidate.sequence:this.oppositeSequence(group,track),site=kind==='primary'?track.primary:track.partner!,ref=this.result!.references.find(r=>r.name===track.name)!;
  const binding=ref.sequence.slice(site.start,site.end).toUpperCase(),oriented=site.direction==='forward'?binding:reverseComplement(binding);
  const cells=[...oriented].map((base,i)=>({base,mismatch:base!==sequence[i]})),offset=kind==='primary'?group.low:this.start(site,track.primary);
  this.hover={id:group.candidate.id,reference:track.name,kind,sequence,cells,left:Math.max(210,Math.min(990-cells.length*7.3,this.x(offset,group)+70)),top:this.y(index)-37};
 }
 hidePrimer(){this.hover=null;}
 mismatchLow(position:number,group:CandidateMap){return group.low+position-1;}
 oppositeSequence(group:CandidateMap,track:Track){return oppositePrimer(group.candidate,group.candidate.sites.indexOf(track.primary))?.sequence||'';}
 oppositeColour(group:CandidateMap,track:Track){return oppositePrimerColour(group.candidate.partner,group.candidate.sites.indexOf(track.primary));}
 x(offset:number,group:CandidateMap){return 140+820*(offset-group.low)/(group.high-group.low);}
 y(index:number){return 75+index*34;}
 height(group:CandidateMap){return 59+group.tracks.length*34;}
 start(site:PrimerSite,primary:PrimerSite){return Math.min(this.offset(site.start,primary),this.offset(site.end,primary));}
 end(site:PrimerSite,primary:PrimerSite){return Math.max(this.offset(site.start,primary),this.offset(site.end,primary));}
 windowLow(offset:number){return offset-Math.floor((this.result?.settings.window||90)/2);}
 windowHigh(offset:number){return offset+Math.ceil((this.result?.settings.window||90)/2);}
}
