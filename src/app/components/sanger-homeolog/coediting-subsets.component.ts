import { Component, ElementRef, EventEmitter, HostListener, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { CoeditingMetric, ComparisonRow, homeologSubsets } from './sanger-homeolog-presentation';

@Component({
  selector:'app-coediting-subsets', standalone:true, imports:[CommonModule],
  template:`<div class="pickers"><div class="picker" *ngFor="let k of counts"><button type="button" class="arrow" [attr.aria-label]="k + '/' + rows.length + ' combinations'" [attr.aria-expanded]="openCount === k" (click)="openCount = openCount === k ? null : k">{{ count === null ? k + '/' + rows.length + ' ' : '' }}▾</button><div class="menu" *ngIf="openCount === k"><button type="button" *ngFor="let item of combinations(k).items" [class.selected]="isSelected(item.names)" [attr.aria-pressed]="isSelected(item.names)" (click)="toggle(item.names)"><span>{{item.names.join(' + ')}}</span><strong>{{item.percent | number:'1.1-1'}}%</strong></button><button type="button" *ngIf="combinations(k).hasMore" (click)="limit = limit + 50">More</button></div></div></div>`,
  styles:[`:host{display:block}.pickers{display:flex;gap:6px;flex-wrap:wrap}.picker{position:relative}.arrow{border:0;background:transparent;color:#64748b;padding:2px 5px;cursor:pointer;font:inherit;font-size:11px}.menu{position:absolute;top:100%;left:0;z-index:20;width:280px;max-width:calc(100vw - 50px);max-height:240px;overflow:auto;background:#fff;border:1px solid #dbe1e8;border-radius:8px;box-shadow:0 8px 24px #0f172a24;padding:5px}.menu button{display:flex;width:100%;justify-content:space-between;gap:10px;border:0;background:#fff;border-radius:5px;padding:8px;color:#475569;font:inherit;font-size:11px;cursor:pointer;text-align:left}.menu button.selected{background:#fffbeb;color:#92400e}.menu button span{overflow-wrap:anywhere}.menu strong{white-space:nowrap}button:focus-visible{outline:2px solid #f59e0b;outline-offset:1px}`]
})
export class CoeditingSubsetsComponent {
  @Input() rows: ComparisonRow[]=[];
  @Input() categories: CoeditingMetric=['inFrame','outOfFrame'];
  @Input() count:number|null=null;
  @Input() selected:readonly string[][]=[];
  @Output() selectedChange=new EventEmitter<string[][]>();
  openCount:number|null=null;
  limit=50;
  constructor(private element:ElementRef<HTMLElement>) {}
  get counts() { return this.count === null ? Array.from({length:Math.max(0,this.rows.length-1)},(_,i)=>i+1) : [this.count]; }
  combinations(count:number) { return homeologSubsets(this.rows,this.categories,count,this.limit); }
  private key(names:readonly string[]) { return JSON.stringify([...names].sort()); }
  isSelected(names:readonly string[]) { return this.selected.some(item=>this.key(item)===this.key(names)); }
  toggle(names:string[]) { this.selectedChange.emit(this.isSelected(names) ? this.selected.filter(item=>this.key(item)!==this.key(names)).map(item=>[...item]) : [...this.selected.map(item=>[...item]),[...names]]); }
  @HostListener('document:click',['$event']) closeOutside(event:MouseEvent) { if(!this.element.nativeElement.contains(event.target as Node)) this.openCount=null; }
}
