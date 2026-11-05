import {Component,Input,OnChanges} from '@angular/core';
import {CommonModule} from '@angular/common';
import {SangerCandidate,SangerDesignResult,PrimerSite,SangerReference} from '../../models/sanger-design';
import {reverseComplement} from '../../utils/biology.utils';
interface Comparison {mismatches:number|null;terminal:number|null;description:string}
interface PairComparison {from:string;to:string;primary:Comparison;partner:Comparison|null}
let serial=0;
export function correspondingPrimerComparison(sequence:string,site:PrimerSite,source:SangerReference,target:SangerReference):Comparison {
 const start=target.cuts[0]+site.start-source.cuts[0],end=start+sequence.length;
 if(start<0||end>target.sequence.length)return{mismatches:null,terminal:null,description:'Corresponding interval is outside the supplied reference'};
 const dna=target.sequence.slice(start,end),query=site.direction==='reverse'?reverseComplement(dna):dna;
 const description=target.name+': '+(start+1)+'–'+end+' · '+site.direction+' · '+query;
 if(/[^ACGT]/.test(query))return {mismatches:null,terminal:null,description:description+' · unresolved bases'};
 const count=(a:string,b:string)=>[...a].filter((base,i)=>base!==b[i]).length;
 return {mismatches:count(sequence,query),terminal:count(sequence.slice(-3),query.slice(-3)),description};
}
@Component({selector:'app-sanger-primer-difference',standalone:true,imports:[CommonModule],templateUrl:'./sanger-primer-difference.component.html',styleUrl:'./sanger-primer-difference.component.css'})
export class SangerPrimerDifferenceComponent implements OnChanges {
 @Input() candidate:SangerCandidate|null=null;@Input() result:SangerDesignResult|null=null;
 readonly id='primer-difference-'+(++serial);pairs:PairComparison[]=[];
 label(screen?:{minimumMismatches:number;minimumExact?:boolean}|null){return !screen?'—':screen.minimumExact===false?'Unresolved':screen.minimumMismatches+' bp';}
 ngOnChanges(){this.pairs=[];if(!this.result||!this.candidate)return;const c=this.candidate;
  for(const site of c.sites){const source=this.result.references.find(r=>r.name===site.reference);if(!source)continue;for(const target of this.result.references){if(target.name===source.name)continue;const partner=c.partner?.sites.find(p=>p.reference===source.name);this.pairs.push({from:source.name,to:target.name,primary:correspondingPrimerComparison(c.sequence,site,source,target),partner:partner&&c.partner?correspondingPrimerComparison(c.partner.sequence,partner,source,target):null});}}
 }
}
