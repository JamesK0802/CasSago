import {Component,Input,OnChanges} from '@angular/core';
import {CommonModule} from '@angular/common';
import {SangerCandidate,SangerDesignResult} from '../../models/sanger-design';
import {TracePlot} from '../../models/sanger-analysis';
import {SangerChromatogramComponent} from '../sanger-chromatogram/sanger-chromatogram.component';
import {selectedAmplicons} from '../../utils/sanger-design-export';
import {reverseComplement} from '../../utils/biology.utils';
import {idealWTTrace} from './ideal-wt-trace';
@Component({selector:'app-sanger-wt-preview',standalone:true,imports:[CommonModule,SangerChromatogramComponent],templateUrl:'./sanger-wt-preview.component.html',styleUrl:'./sanger-wt-preview.component.css'})
export class SangerWTPreviewComponent implements OnChanges {
 @Input() result:SangerDesignResult|null=null;@Input() selected:SangerCandidate|null=null;
 references:{name:string;sequence:string}[]=[];choice=-1;plot:TracePlot|null=null;start=0;end=0;
 ngOnChanges(){
  this.references=[];this.plot=null;this.choice=-1;
  if(!this.result||!this.selected?.partner)return;
  const c=this.selected;
  this.references=selectedAmplicons(this.result,c).map(a=>({name:a.sourceReference,sequence:(c.direction==='reverse'?reverseComplement(a.sequence):a.sequence).slice(c.sequence.length)}));
  const cuts=c.sites.flatMap(h=>h.cutOffsets),margin=Math.min(75,Math.ceil(this.result.settings.window/2));
  const usableEnd=Math.min(this.result.settings.readEnd,...this.references.map(r=>r.sequence.length))-50;
  this.start=Math.max(50,this.result.settings.readStart,Math.min(...cuts)-margin);
  this.end=Math.min(usableEnd,Math.max(...cuts)+margin,this.start+400);
  this.update();
 }
 choose(index:number){this.choice=index;this.update();}
 private update(){this.plot=idealWTTrace((this.choice<0?this.references:[this.references[this.choice]]).map(r=>r.sequence),this.start,this.end);}
 get title(){return this.choice<0&&this.references.length>1?'Mixed WT · ideal signal':(this.references[this.choice<0?0:this.choice]?.name||'WT')+' · ideal signal';}
 get caption(){return (this.choice<0&&this.references.length>1?'Equal contributions from '+this.references.length+' references. ':'')+'No indels or added noise. Schematic peaks; not measured data or a prediction of experimental peak heights. Mixed-base letters use IUPAC codes.';}
}
