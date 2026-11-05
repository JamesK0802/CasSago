import {Component,Input,OnChanges,AfterViewInit,OnDestroy,ViewChild,ElementRef} from '@angular/core';
import {Chart} from 'chart.js/auto';
import {SangerAnalysisResult} from '../../models/sanger-analysis';
/** The CRISPR dashboard chart cards, fed by signal fractions rather than read counts. */
@Component({selector:'app-sanger-result-charts',standalone:true,template:`<div class="charts-section"><div class="chart-card pie"><canvas #composition aria-label="Signal composition"></canvas></div><div class="chart-card"><canvas #spectrum aria-label="Indel size spectrum"></canvas></div><div class="chart-card donut"><canvas #shares aria-label="Indel editing"></canvas></div></div>`,styles:[`:host{display:block}.charts-section{display:grid;grid-template-columns:1fr 1.4fr 1fr;gap:16px;margin-bottom:24px}.chart-card{height:230px;min-width:0;padding:16px;position:relative}canvas{max-width:100%}@media(max-width:850px){.charts-section{grid-template-columns:1fr}.chart-card{height:220px}}`]})
export class SangerResultChartsComponent implements OnChanges,AfterViewInit,OnDestroy {
 @Input() substitution:number|null=null;
 @Input() profile:{share:number;wt:number|null;inFrame:number|null;outOfFrame:number|null}|null=null;@Input() components:{size:number;fraction:number}[]|null=null;@Input() targetLabel='Target 1';@Input() targetIndex=0;@Input() zeroNetUnresolved=false;
 @Input() result:SangerAnalysisResult|null=null;@Input() selected=0;
 @ViewChild('composition') composition?:ElementRef<HTMLCanvasElement>;@ViewChild('spectrum') spectrum?:ElementRef<HTMLCanvasElement>;@ViewChild('shares') shares?:ElementRef<HTMLCanvasElement>;
 private ready=false;private charts:Chart[]=[];
 ngAfterViewInit(){this.ready=true;this.draw();}ngOnChanges(){if(this.ready)this.draw();}ngOnDestroy(){this.charts.forEach(c=>c.destroy());}
 private draw(){this.charts.forEach(c=>c.destroy());this.charts=[];const r=this.result;if(!r?.mixture||!this.composition||!this.spectrum||!this.shares)return;const p=this.profile||(r.pooled?.amplicons[this.selected]||r.mixture),components=this.components||(r.pooled?.amplicons[this.selected]?.components||r.mixture.components);
  const options={responsive:true,maintainAspectRatio:false,animation:false as const,plugins:{legend:{position:'bottom' as const,labels:{boxWidth:10,font:{size:10}}}}};
  const a=r.annotations?.[this.selected],sequenceGroups=r.diagnostics?.algorithmVersion?.startsWith('sanger-sequence-groups-v'),sub=this.substitution??(sequenceGroups&&a?.estimable?a.components.filter(c=>{const es=(c.events||[]).filter(e=>e.target===this.targetIndex+1||e.target===0);return !!es.length&&es.every(e=>e.inserted===e.end-e.start);}).reduce((sum,c)=>sum+c.fraction,0):0);
  const compositionLabels=this.zeroNetUnresolved?['In-frame','Out-of-frame','WT / substitution unresolved']:null,compositionData=this.zeroNetUnresolved?[p.inFrame??0,p.outOfFrame??0,Math.max(0,1-(p.inFrame??0)-(p.outOfFrame??0))]:null;
  this.charts.push(new Chart(this.composition.nativeElement,{type:'doughnut',data:{labels:compositionLabels||[sequenceGroups?'Unedited':'0 bp','In-frame','Out-of-frame',...(sequenceGroups?['Substitution']:[])],datasets:[{data:(compositionData||[p.wt??0,p.inFrame??0,p.outOfFrame??0,...(sequenceGroups?[sub]:[])]).map(v=>v*100),backgroundColor:this.zeroNetUnresolved?['#0072B2','#E69F00','#E6C857']:['#009E73','#0072B2','#E69F00','#CC79A7'],borderWidth:0}]},options}));
  this.charts.push(new Chart(this.spectrum.nativeElement,{
   type:'bar',
   data:{labels:components.map(c=>(c.size>0?'+':'')+c.size),datasets:[{label:'Signal %',data:components.map(c=>c.fraction*100),backgroundColor:components.map(c=>c.size===0?'#009E73':c.size%3===0?'#0072B2':'#E69F00')}]},
   options:{...options,plugins:{legend:{display:false}},scales:{
    y:{beginAtZero:true,title:{display:true,text:'Signal %',font:{size:10}}},
    x:{title:{display:true,text:'Indel (bp)',font:{size:10}},ticks:{font:{size:9}}}
   }}
  }));
  this.charts.push(new Chart(this.shares.nativeElement,{type:'doughnut',data:{labels:this.zeroNetUnresolved?['Indel edited','WT / substitution unresolved']:['Indel edited',sequenceGroups?'Unedited':'0 bp',...(sequenceGroups?['Substitution']:[])],datasets:[{data:(this.zeroNetUnresolved?[(p.inFrame??0)+(p.outOfFrame??0),Math.max(0,1-(p.inFrame??0)-(p.outOfFrame??0))]:[(p.inFrame??0)+(p.outOfFrame??0),p.wt??0,...(sequenceGroups?[sub]:[])]).map(v=>v*100),backgroundColor:this.zeroNetUnresolved?['#E69F00','#E6C857']:['#E69F00','#009E73','#CC79A7'],borderWidth:0}]},options}));
 }
}
