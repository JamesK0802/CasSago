export interface ComparisonRow {name:string;share:number;targetId:string;reads:number;wt:number|null;substitution:number|null;inFrame:number|null;outOfFrame:number|null;edited:number|null;target?:{top_groups?:{classification:string;net_indel:number;read_pct:number}[]}}
export interface FigureSource {label:string;rows:ComparisonRow[]}
export type CoeditingCategory = 'wt' | 'inFrame' | 'outOfFrame';
export type CoeditingMetric = readonly CoeditingCategory[];
export const COEDITING_METRICS: { value: CoeditingCategory; label: string }[] = [
  { value: 'wt', label: 'WT' },
  { value: 'inFrame', label: 'In-frame' },
  { value: 'outOfFrame', label: 'Out-of-frame' }
];
export function coeditingMetricLabel(metric: CoeditingMetric = ['inFrame', 'outOfFrame']): string {
  return COEDITING_METRICS.filter(item => metric.includes(item.value)).map(item => item.label).join(' + ') || '—';
}
function coeditingRates(rows: ComparisonRow[], metric: CoeditingMetric) {
  const categories = [...new Set(metric)];
  if (!categories.length || rows.length < 2 || rows.some(row => row.edited === null || categories.some(category => row[category] === null))) return null;
  return rows.map(row => ({name:row.name, rate:Math.max(0,Math.min(100,categories.reduce((sum,category)=>sum+row[category]!,0)))}));
}
export interface CoeditingSubset { names: string[]; percent: number; }
export function homeologSubset(rows: ComparisonRow[], metric: CoeditingMetric, names: readonly string[]): CoeditingSubset | null {
  const rates = coeditingRates(rows,metric), selected = new Set(names);
  if (!rates || !selected.size || selected.size >= rows.length || selected.size !== names.length || names.some(name=>!rates.some(row=>row.name===name))) return null;
  return { names:rates.filter(row=>selected.has(row.name)).map(row=>row.name), percent:rates.reduce((p,row)=>p*(selected.has(row.name) ? row.rate/100 : 1-row.rate/100),1)*100 };
}
export function homeologSubsets(rows: ComparisonRow[], metric: CoeditingMetric, count: number, limit = 50) {
  const items: CoeditingSubset[] = [];
  if (!coeditingRates(rows,metric) || count < 1 || count >= rows.length) return {items,hasMore:false};
  const visit = (start:number, names:string[]) => {
    if (items.length > limit) return;
    if (names.length === count) { const item=homeologSubset(rows,metric,names); if(item) items.push(item); return; }
    for(let i=start;i<=rows.length-(count-names.length);i++) { visit(i+1,[...names,rows[i].name]); if(items.length>limit) break; }
  };
  visit(0,[]);
  return {items:items.slice(0,limit),hasMore:items.length>limit};
}
export function homeologCombination(rows: ComparisonRow[], metric: CoeditingMetric = ['inFrame', 'outOfFrame']) {
  const values = coeditingRates(rows,metric);
  if (!values) return null;
  const rates = values.map(value => Math.max(0, Math.min(1, value.rate! / 100)));
  let distribution = [1];
  for (const p of rates) {
    const next = Array(distribution.length + 1).fill(0);
    distribution.forEach((value, k) => { next[k] += value * (1 - p); next[k + 1] += value * p; });
    distribution = next;
  }
  const ordered = [...values].sort((a, b) => a.rate! - b.rate!);
  return { weakest: ordered[0], strongest: ordered[ordered.length - 1],
    gap: ordered[ordered.length - 1].rate! - ordered[0].rate!,
    all: distribution[distribution.length - 1] * 100, any: (1 - distribution[0]) * 100,
    distribution: distribution.map((p, count) => ({ count, percent: p * 100 })) };
}

export const FIGURE_COLORS={okabe:['#009E73','#CC79A7','#0072B2','#E69F00']};
