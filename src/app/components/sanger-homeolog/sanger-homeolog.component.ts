import { Component, EventEmitter, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';

import { ComparisonRow, homeologCombination, FIGURE_COLORS, FigureSource, CoeditingMetric, coeditingMetricLabel } from './sanger-homeolog-presentation';
import { CoeditingBasisComponent } from './coediting-basis.component';
import { CoeditingSubsetsComponent } from './coediting-subsets.component';

@Component({
  selector: 'app-sanger-homeolog', standalone: true, imports: [CommonModule, FormsModule, CoeditingBasisComponent, CoeditingSubsetsComponent],
  template: `
    <div class="comparison" id="guide-homeo-comparison" [attr.data-selection-count]="selectedNames.length">
      <div class="selectors"><label>Group<select aria-label="Homeolog group" [ngModel]="activeGroup" (ngModelChange)="setGroup($event)"><option *ngFor="let name of groupNames" [value]="name">{{ name }}</option></select></label>
      <label>Target<select aria-label="Comparison target" [ngModel]="activeTarget" (ngModelChange)="target = $event"><option *ngFor="let id of targetIds" [value]="id">{{ id }}</option></select></label></div>
      <div class="selection-heading"><strong>Compare references</strong><div><button type="button" id="guide-homeo-select-all" (click)="chooseAll()">All</button><button type="button" (click)="selection = []">Clear</button></div></div>
      <div class="reference-selection"><div class="reference-picker" id="guide-homeo-reference-picker">
        <button type="button" *ngFor="let gene of groupGenes; let i = index" [attr.data-reference]="gene.gene"
          [attr.id]="i === 1 ? 'guide-homeo-second-reference' : null" [attr.aria-pressed]="selectedNames.includes(gene.gene)"
          [class.selected]="selectedNames.includes(gene.gene)" (click)="toggle(gene.gene)"><span>{{ selectedNames.includes(gene.gene) ? '✓' : '+' }}</span> {{ gene.gene }}</button>
      </div><app-coediting-basis [selected]="coeditingMetric" (selectedChange)="coeditingMetric = $event; coeditingMetricChange.emit($event)" /></div>
      <ng-container *ngIf="rows.length; else empty">
        <section class="combination-analysis">
          <h3>Co-editing model <small>· {{ metricLabel }} · independence assumed</small></h3>
          <ng-container *ngIf="combination as combo; else incomplete">
            <div class="model">
              <div class="model-totals"><span>All {{ rows.length }} selected: <strong>{{ format(combo.all) }}</strong></span><span>At least one: <strong>{{ format(combo.any) }}</strong></span></div>
              <div class="model-row" *ngFor="let item of combo.distribution; trackBy: trackCount"><span style="display:flex;align-items:center;min-width:70px">{{ item.count }} / {{ rows.length }}<app-coediting-subsets *ngIf="item.count > 0 && item.count < rows.length" [count]="item.count" [rows]="rows" [categories]="coeditingMetric" [selected]="coeditingSubsets" (selectedChange)="coeditingSubsets = $event; coeditingSubsetsChange.emit($event)" /></span><div class="model-track"><i [style.width.%]="item.percent"></i></div><strong>{{ format(item.percent) }}</strong></div>
            </div>
          </ng-container>
          <ng-template #incomplete><p class="note" *ngIf="coeditingMetric.length">{{ rows.length < 2 ? 'Select at least two references.' : 'Every selected reference needs evaluable signal for this target.' }}</p></ng-template>
        </section>
        <section class="profile"><h3>Editing profile</h3>
          <div class="legend"><span *ngFor="let name of labels; let i = index"><i [style.background]="colors[i]"></i>{{ name }}</span></div>
          <div *ngFor="let row of rows" class="profile-row">
            <div class="row-heading"><button type="button" (click)="inspect.emit(row.name)" title="Show this reference in Standard">{{ row.name }} ↗</button><small>{{row.share|number:'1.1-1'}}% signal share</small><strong>{{ row.wt === null ? 'Not evaluated' : (row.wt | number:'1.1-1') + '% WT' }}</strong></div>
            <div class="stack" [attr.aria-label]="row.name + ': ' + (row.wt === null ? 'not evaluated' : row.wt + ' percent WT')">
              <span *ngFor="let value of [row.wt, row.substitution, row.inFrame, row.outOfFrame]; let i = index" [style.width.%]="value ?? 0" [style.background]="colors[i]" [title]="labels[i] + ': ' + (value ?? 0).toFixed(2) + '%'"></span>
            </div>
          </div>
        </section>
        <details open><summary>Signal summary</summary>
          <div class="summary-scroll"><table><thead><tr><th>Reference</th><th>WT</th><th>Sub.</th><th>In-frame</th><th>Out-of-frame</th></tr></thead>
          <tbody><tr *ngFor="let row of rows"><th>{{ row.name }}</th><td>{{ format(row.wt) }}</td><td>{{ format(row.substitution) }}</td><td>{{ format(row.inFrame) }}</td><td>{{ format(row.outOfFrame) }}</td></tr></tbody></table></div>
        </details>
        <details open><summary>Major mutation groups by reference</summary>
          <div *ngFor="let row of rows" class="allele-summary"><strong>{{ row.name }}</strong>
            <div *ngFor="let group of row.target?.top_groups?.slice(0, 5)"><span>{{ group.classification }} · {{ group.net_indel > 0 ? '+' : '' }}{{ group.net_indel }} bp</span><span>{{ group.read_pct | number:'1.1-1' }}%</span></div>
            <small *ngIf="!row.target?.top_groups?.length">No groups available</small>
          </div>
        </details>
      </ng-container>
      <ng-template #empty><div class="empty">Select references to compare.</div></ng-template>
    </div>`,
  styles: [`
    .reference-selection{display:flex;align-items:center;justify-content:space-between;gap:12px;flex-wrap:wrap}.reference-selection .reference-picker{flex:1;min-width:0}.reference-selection app-coediting-basis{margin-left:auto}
    .balance{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px}.balance>div{display:grid;align-content:start;gap:6px;padding:12px;background:#f8fafc;border-radius:8px;overflow-wrap:anywhere}.balance small{font-size:11px;color:#64748b}.balance strong{font-size:14px}.balance span{font-size:13px}.model-totals{display:flex;flex-wrap:wrap;gap:16px;font-size:12px;margin:14px 0}.model-row{display:flex;align-items:center;gap:10px;font-size:11px;margin:9px 0}.model-row>span{min-width:45px}.model-row>strong{min-width:50px;text-align:right}.model-track{flex:1;background:#f1f5f9;height:14px;border-radius:3px;overflow:hidden}.model-track i{display:block;height:100%;background:#0072b2}.combination-analysis>label{max-width:220px;margin-bottom:18px}
    :host{display:block;min-width:0}.comparison{padding:18px}.selectors{display:flex;gap:12px}label{display:grid;gap:5px;font-size:11px;color:var(--color-text-secondary,#64748b)}.selectors label{flex:1;min-width:0}
    button,input,select{font:inherit;color:inherit;border:1px solid var(--color-border,#dbe1e8);background:var(--color-surface,#fff);border-radius:7px;padding:7px 10px;font-size:12px}button{cursor:pointer}button:focus-visible,input:focus-visible,select:focus-visible{outline:2px solid #f59e0b;outline-offset:2px}select{width:100%;min-width:0}
    .selection-heading{display:flex;justify-content:space-between;align-items:center;margin:20px 0 10px;font-size:12px}.selection-heading button{border:0;padding:4px 7px;background:transparent}.reference-picker{display:flex;flex-wrap:wrap;gap:7px;max-height:150px;overflow:auto;padding:3px}.reference-picker .selected{border-color:#e9b34b;background:#fffbeb;color:#92400e}.context,.note{font-size:11px;line-height:1.6;color:var(--color-text-secondary,#64748b)}.context{margin:14px 0}
    .overview{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:9px}.overview>div{border:1px solid var(--color-border,#e2e8f0);border-radius:10px;padding:12px}.overview strong{font-size:23px;display:block}.overview small{font-size:10px;color:var(--color-text-secondary,#64748b)}.threshold{display:flex;align-items:center;margin:12px 0;gap:7px}.threshold input{width:65px}
    section{margin-top:22px}h3{font-size:14px;margin:0 0 14px}.legend{display:flex;flex-wrap:wrap;gap:10px;font-size:10px;color:var(--color-text-secondary,#64748b)}.legend i{display:inline-block;width:8px;height:8px;border-radius:2px;margin-right:4px}.profile-row{margin-top:18px}.row-heading{display:flex;align-items:center;gap:10px;font-size:11px;flex-wrap:wrap}.row-heading button{border:0;padding:0;background:transparent;font-weight:700}.row-heading small{color:#64748b}.row-heading strong{margin-left:auto}.stack{height:24px;background:var(--color-surface-alt,#f1f5f9);display:flex;margin-top:7px;border-radius:4px;overflow:hidden}.stack>span{height:100%}
    .summary-scroll{overflow:auto}table{width:100%;font-size:11px;border-collapse:collapse}th,td{padding:10px 8px;border-bottom:1px solid var(--color-border,#e2e8f0);text-align:right;white-space:nowrap}th:first-child{text-align:left}thead{color:#64748b}details{margin-top:18px;border-top:1px solid var(--color-border,#e2e8f0);padding-top:13px}summary{font-size:12px;font-weight:600;cursor:pointer}.allele-summary{margin:15px 0;font-size:12px}.allele-summary>div{display:flex;justify-content:space-between;padding:6px 0;gap:12px;font-size:11px}.empty{padding:50px 12px;text-align:center;line-height:1.8;color:#64748b}
  `]
})
export class SangerHomeologComponent {
  @Input() profiles: ComparisonRow[] = [];
  @Input() scopeLabel = '';
  @Input() coeditingMetric: CoeditingMetric = ['inFrame', 'outOfFrame'];
  @Output() coeditingMetricChange = new EventEmitter<CoeditingMetric>();
  @Input() coeditingSubsets: string[][] = [];
  @Output() coeditingSubsetsChange = new EventEmitter<string[][]>();
  get metricLabel() { return coeditingMetricLabel(this.coeditingMetric); }
  trackCount(_: number, item: {count: number}) { return item.count; }
  @Output() inspect = new EventEmitter<string>();
  group = ''; target = ''; selection: string[] | null = null;
  colors = FIGURE_COLORS.okabe;
  labels = ['WT', 'Substitution', 'In-frame indel', 'Out-of-frame indel'];
  private genesCache?:{source:ComparisonRow[];genes:{gene:string}[]};
  get groupGenes(){if(this.genesCache?.source!==this.profiles)this.genesCache={source:this.profiles,genes:[...new Set(this.profiles.map(row=>row.name))].map(gene=>({gene}))};return this.genesCache.genes;}
  get groupNames(){return ['References'];}
  get activeGroup(){return 'References';}
  get targetIds(){return [...new Set(this.profiles.map(row=>row.targetId))];}
  get activeTarget(){return this.targetIds.includes(this.target)?this.target:this.targetIds[0]||'';}
  get selectedNames(){return [...new Set(this.profiles.map(row=>row.name))].filter(name=>this.selection===null||this.selection.includes(name));}
  get rows(){return this.profiles.filter(row=>row.targetId===this.activeTarget&&this.selectedNames.includes(row.name));}
  get combination() { return homeologCombination(this.rows, this.coeditingMetric); }
  get figureSource(): FigureSource { return { label: `${this.scopeLabel} · ${this.activeGroup} · ${this.activeTarget}`, rows: this.rows }; }
  setGroup(value: string) { this.group = value; this.target = ''; this.selection = null; }
  chooseAll() { this.selection = null; }
  toggle(name: string) { const selected = this.selectedNames; this.selection = selected.includes(name) ? selected.filter(n => n !== name) : [...selected, name]; }
  format(value: number | null) { return value === null ? 'N/A' : value.toFixed(1) + '%'; }
}
