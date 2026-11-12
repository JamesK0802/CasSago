import { Component, EventEmitter, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { CoeditingCategory, CoeditingMetric, COEDITING_METRICS } from './sanger-homeolog-presentation';

@Component({
  selector: 'app-coediting-basis', standalone: true, imports: [CommonModule],
  template: `<div class="basis" role="group" aria-label="Co-editing categories"><button type="button" *ngFor="let category of categories" [attr.aria-pressed]="selected.includes(category.value)" [class.selected]="selected.includes(category.value)" (click)="toggle(category.value)">{{ category.label }}</button></div>`,
  styles: [`:host{display:block;max-width:100%}.basis{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:7px;width:345px;max-width:100%}button{box-sizing:border-box;height:34px;min-width:0;padding:7px 8px;border:1px solid var(--color-border,#dbe1e8);border-radius:7px;background:var(--color-surface,#fff);color:var(--color-text-secondary,#64748b);font:inherit;font-size:12px;line-height:18px;white-space:nowrap;cursor:pointer}button.selected{border-color:#e9b34b;background:#fffbeb;color:#92400e}button:focus-visible{outline:2px solid #f59e0b;outline-offset:2px}`]
})
export class CoeditingBasisComponent {
  @Input() selected: CoeditingMetric = ['inFrame', 'outOfFrame'];
  @Output() selectedChange = new EventEmitter<CoeditingMetric>();
  categories = COEDITING_METRICS;
  toggle(value: CoeditingCategory) {
    const next = this.selected.includes(value) ? this.selected.filter(item => item !== value)
      : this.categories.filter(item => item.value === value || this.selected.includes(item.value)).map(item => item.value);
    this.selectedChange.emit(next);
  }
}
