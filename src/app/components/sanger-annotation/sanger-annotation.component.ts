import {exactDisplayView,exactFractionLabel} from '../../utils/sanger-exact-display';
import {projectTargetAnnotation,projectTargetTruth} from '../../utils/sanger-target-profile';
import {BenchmarkTruth} from '../../utils/sanger-benchmark';
import {sangerTruthGroups,TruthGroupView} from '../../utils/sanger-truth-groups';
import {SynthComparison} from '../../utils/sanger-synthetic-benchmark';
import {benchmarkPattern} from '../../utils/sanger-benchmark';
import {SangerChromatogramComponent} from '../sanger-chromatogram/sanger-chromatogram.component';import {rawSeparatedGroupTrace} from '../sanger-chromatogram/sanger-raw-group-trace';import {SangerAnalysisResult,SangerTrace} from '../../models/sanger-analysis';
import {Component,Input,Output,EventEmitter,OnChanges,ViewChildren,QueryList,ElementRef,AfterViewChecked,SimpleChanges} from '@angular/core';import {CommonModule} from '@angular/common';
import {SangerAnnotation} from '../../models/sanger-analysis';import {sangerAnnotationViews,SangerAnnotationView} from '../../utils/sanger-annotation';
@Component({selector:'app-sanger-annotation',standalone:true,imports:[CommonModule,SangerChromatogramComponent],templateUrl:'./sanger-annotation.component.html',styleUrl:'./sanger-annotation.component.css'})
export class SangerAnnotationComponent implements OnChanges,AfterViewChecked{
 @Input() truth:BenchmarkTruth|null=null;truthViews:TruthGroupView[][]=[];
 @Input() equivalentReferences:string[]=[];
 @Input() benchmark:SynthComparison[]=[];@Input() threshold=1;showSmall=false;
 get exactSet(){return this.analysis?.exactAlternatives?.find(s=>s.reference===this.active?.name&&s.targets.includes(this.targetIndex+1));}
 @Output() exactSelectionChange=new EventEmitter<{reference:string;option:number}>();
 selectExact(event:Event){if(this.exactSet)this.exactSelectionChange.emit({reference:this.exactSet.reference,option:Number((event.target as HTMLSelectElement).value)});}
 isExactAmbiguous(row:{id?:string;sequenceKey?:string}){if(row.id?.startsWith('exact-insertion:'))return true;return !!row.sequenceKey&&!![...(this.exactSet?.ambiguousSequences[String(this.targetIndex+1)]||[]),...(this.exactSet?.ambiguousPlacements?.[String(this.targetIndex+1)]||[])].includes(row.sequenceKey);}
 private labelSet:typeof this.exactSet;private labelTarget=-1;private optionLabels:string[]=[];
 exactLabel(index:number){const set=this.exactSet;if(set!==this.labelSet||this.labelTarget!==this.targetIndex){this.labelSet=set;this.labelTarget=this.targetIndex;const reference=this.annotations[this.selectedIndex]?.reference||'';this.optionLabels=(set?.options||[]).map((option,i)=>{const groups=new Map<string,number>();for(const c of option.components){const events=(c.events||[]).filter(e=>e.target===this.targetIndex+1||e.target===0),label=!events.length?'WT':events.map(e=>e.start===e.end?'+'+e.insertedSequence:!e.inserted?'−'+(e.end-e.start)+' bp @'+(e.start+1):reference.slice(e.start,e.end)+'→'+e.insertedSequence).join(' / ');groups.set(label,(groups.get(label)||0)+c.fraction);}return 'Combination '+(i+1)+' · '+[...groups].map(([label,fraction])=>label+' '+(fraction*100).toFixed(2)+'%').join(' · ');});}return this.optionLabels[index]||'';}
 get linkage(){if(this.exactSet)return [];return this.analysis?.diagnostics?.substitutionLinkage?.filter(a=>a.reference===this.active?.name&&a.target===this.targetIndex+1)||[];}
 get linkageWtRange(){const a=this.linkage;return a.length?`${(Math.min(...a.map(x=>x.wtRange[0]))*100).toFixed(2)}–${(Math.max(...a.map(x=>x.wtRange[1]))*100).toFixed(2)}%`:'';}
 get wholeSequenceFit(){return this.analysis?.diagnostics?.algorithmVersion?.startsWith('sanger-sequence-groups-v');}
 get cutoff(){return Math.max(0,Math.min(100,Number(this.threshold)||0))/100;}
 showGroup(row:{edited:boolean;fraction:number;fractionRange?:[number,number]}){return !row.edited||(row.fractionRange?.[1]??row.fraction)+1e-12>=this.cutoff;}
 get smallGroups(){return this.active?.rows.filter(r=>!this.showGroup(r))||[];}
 get smallMass(){return this.smallGroups.reduce((sum,r)=>sum+r.fraction,0);}
 truthGroup(row:{size:number;patternLabel?:string;sequenceKey?:string}){return this.benchmark.find(t=>t.reference===this.active?.name&&t.scope==='group'&&t.pattern===benchmarkPattern(row)&&(t.sequenceScope==='sequence'?!!row.sequenceKey&&t.sequenceKey===row.sequenceKey:true));}
 truthGroupLabel(t:SynthComparison){if(t.truth===null)return 'GT not attributable';if(t.sequenceScope==='pattern'&&this.active?.rows.filter(r=>benchmarkPattern(r)===t.pattern).length>1)return `GT size total ${t.truth.toFixed(2)}%`;return `GT ${t.truth.toFixed(2)}%${t.delta===null?'':` · Δ ${t.delta>=0?'+':''}${t.delta.toFixed(2)} pp`}`;}
 get missingTruth(){return this.benchmark.filter(t=>t.reference===this.active?.name&&t.scope==='group'&&(t.truth??0)>0&&!this.active?.rows.some(r=>benchmarkPattern(r)===t.pattern&&(t.sequenceScope==='sequence'?!!r.sequenceKey&&r.sequenceKey===t.sequenceKey:true)));}
 @Input() sampleTrace:SangerTrace|null=null;@Input() analysis:SangerAnalysisResult|null=null;expanded=new Set<string>();private traceCache=new Map<string,ReturnType<typeof rawSeparatedGroupTrace>>();
 groupKey(row:{size:number;patternLabel?:string;id?:string}){return this.selectedIndex+':'+this.targetIndex+':'+(row.id??row.patternLabel??row.size);}
 toggleTrace(row:{size:number;patternLabel?:string;id?:string}){const key=this.groupKey(row);if(this.expanded.has(key))this.expanded.delete(key);else this.expanded.add(key);}
 groupTrace(row:NonNullable<typeof this.active>['rows'][number]){const key=this.groupKey(row);if(!this.traceCache.has(key))this.traceCache.set(key,this.analysis?rawSeparatedGroupTrace(this.analysis,this.sampleTrace,this.selectedIndex,row,this.window,this.targetIndex):{plot:null,note:'The original raw signal is unavailable.'});return this.traceCache.get(key)!;}
 @Input() targetIndex=0;@Input() targetLabel='';@Input() annotations:SangerAnnotation[]=[];@Input() window=200;@Input() selectedIndex=0;@Output() selectedIndexChange=new EventEmitter<number>();@ViewChildren('annotationScroll') scrolls?:QueryList<ElementRef<HTMLElement>>;views:SangerAnnotationView[]=[];private center=true;
 private referenceTokenCache=new Map<string,{type:'equal';val:string}[]>();
 referenceTokens(reference:string){if(!this.referenceTokenCache.has(reference))this.referenceTokenCache.set(reference,[{type:'equal',val:reference}]);return this.referenceTokenCache.get(reference)!;}
 get truthActive(){return this.truthViews[this.selectedIndex]||[];}
 truthLabel(row:SangerAnnotationView['rows'][number]){const gt=this.truthActive.find(t=>t.row.sequenceKey&&t.row.sequenceKey===row.sequenceKey);return gt?`GT ${gt.label} · ${(gt.fraction*100).toFixed(2)}% · Exact`:'No exact GT group';}
 truthMatches(index:number){const gt=this.truthActive[index],matches=(this.active?.rows||[]).flatMap((row,i)=>gt.row.sequenceKey&&row.sequenceKey===gt.row.sequenceKey?[{row,index:i+1}]:[]);return matches.length?`Group ${matches.map(m=>m.index).join(' + ')} · ${(matches.reduce((sum,m)=>sum+m.row.fraction,0)*100).toFixed(2)}%`:'No exact estimated group';}
 truthDetected(gt:TruthGroupView){return !!this.active?.estimable&&(this.equivalentReferences.includes(this.active.name)||!!gt.row.sequenceKey&&this.active.rows.some(row=>row.sequenceKey===gt.row.sequenceKey&&this.showGroup(row)));}
 truthDetectionLabel(gt:TruthGroupView){const rows=(this.active?.rows||[]).filter(row=>!!gt.row.sequenceKey&&row.sequenceKey===gt.row.sequenceKey&&this.showGroup(row));if(!this.active?.estimable)return 'Not detected: fit QC review';if(!rows.length)return this.equivalentReferences.includes(this.active?.name||'')?'Exactly signal-equivalent; individual sequence assignment is unresolved':'Not detected: no matching estimated sequence above the threshold';const fraction=rows.reduce((sum,row)=>sum+row.fraction,0)*100;return 'Detected sequence · estimated '+fraction.toFixed(2)+'% · Δ '+(fraction-gt.fraction*100).toFixed(2)+' pp';}
 fractionText(row:SangerAnnotationView['rows'][number]){return exactFractionLabel(row);}
 badge(row:SangerAnnotationView['rows'][number]){if(!row.edited)return 'WT';const events=row.representative?.events||[];return row.size===0&&events.length&&events.every(e=>e.inserted===e.end-e.start)?'Substitution':(row.size>0?'+':'')+row.size+' bp';}
 private syncing=false;
 syncScroll(event:Event){if(this.syncing)return;this.syncing=true;const source=event.target as HTMLElement;this.scrolls?.forEach(ref=>{if(ref.nativeElement!==source)ref.nativeElement.scrollLeft=source.scrollLeft;});this.syncing=false;}
 get active(){return this.views[Math.min(Math.max(0,this.selectedIndex),Math.max(0,this.views.length-1))];}
 choose(i:number){this.selectedIndex=i;this.selectedIndexChange.emit(i);this.center=true;}
 private viewCache=new Map<string,{views:SangerAnnotationView[];truthViews:TruthGroupView[][]}>();
 ngOnChanges(changes?:SimpleChanges){if(!changes||changes['analysis']||changes['annotations']||changes['window']||changes['truth'])this.viewCache.clear();if(!changes||changes['analysis']||changes['sampleTrace']||changes['window'])this.traceCache.clear();this.rebuild();}
 rebuild(){const key=this.targetIndex+':'+this.window;const cached=this.viewCache.get(key);if(cached){this.views=cached.views;this.truthViews=cached.truthViews;this.center=true;return;}this.views=this.annotations.map(a=>exactDisplayView(a,this.targetIndex,Math.floor(this.window/2),this.analysis?.exactAlternatives?.find(s=>s.reference===a.name&&s.targets.includes(this.targetIndex+1))));this.truthViews=this.annotations.map(a=>sangerTruthGroups(this.truth&&(a.cuts?.length||1)>1?projectTargetTruth(this.truth,a,this.targetIndex):this.truth,a,Math.floor(this.window/2),a.cuts?.[this.targetIndex]??a.cut));this.truthViews=this.truthViews.map((truth,index)=>{const order=this.views[index].rows,rank=(gt:TruthGroupView)=>{const row=order.findIndex(row=>!!gt.row.sequenceKey&&row.sequenceKey===gt.row.sequenceKey);return row<0?Infinity:row;};return truth.map((gt,original)=>({gt,original})).sort((a,b)=>rank(a.gt)-rank(b.gt)||a.original-b.original).map(t=>t.gt);});this.viewCache.set(key,{views:this.views,truthViews:this.truthViews});this.center=true;}
 ngAfterViewChecked(){if(!this.center||!this.scrolls?.length)return;this.center=false;this.scrolls.forEach((ref,i)=>{const el=ref.nativeElement,cut=this.active?.cutIndex||0;el.scrollLeft=Math.max(0,150+15+cut*13-el.clientWidth/2);});}
}
