import {describe,it,expect} from 'vitest';import {sangerAnnotationViews} from '../../src/app/utils/sanger-annotation';import {SangerAnnotation} from '../../src/app/models/sanger-analysis';
const input=(direction:'forward'|'reverse'='forward'):SangerAnnotation=>({name:'A',reference:'ACGTTCAACGTGAGTC',cut:8,direction,share:.4,estimable:true,components:[{size:0,fraction:.6},{size:-2,fraction:.25},{size:2,fraction:.145},{size:1,fraction:.005}]});
describe('Sanger net-indel model annotation',()=>{
 it('uses main-annotation boundary coordinates and preserves reference bases',()=>{const a=input(),view=sangerAnnotationViews([a],2,5)[0];expect(view.reference).toBe(a.reference.slice(3,13));expect(view.cutIndex).toBe(5);expect(view.rows[0].classification).toBe('No net indel model');expect(view.rows[0].tokens.map(t=>t.val).join('')).toBe(view.reference);});
 it('does not fabricate a cut-anchored deletion or insertion when only length is known',()=>{
  for(const direction of ['forward','reverse'] as const){const view=sangerAnnotationViews([input(direction)],0,8)[0];for(const row of view.rows.filter(r=>r.size!==0)){expect(row.tokens.every(t=>t.type==='equal')).toBe(true);expect(row.positionLabel).toBe('Position not estimated');expect(row.placements).toEqual([]);}}
 });
 it('renders a supplied physical interval in either direction without relocating it to the cut',()=>{
  for(const direction of ['forward','reverse'] as const){const a={...input(direction),components:[{size:-2,fraction:1,events:[{target:1,start:3,end:5,inserted:0}]}]},row=sangerAnnotationViews([a],0,8)[0].rows[0];expect(row.tokens.find(t=>t.type==='delete')?.val).toBe(a.reference.slice(3,5));expect(row.placements[0].location).toBe('Deletion bases 4–5');expect(row.positionLabel).toContain('Best fit:');}
 });
 it('keeps all fitted denominators and does not mutate results when thresholding',()=>{const a=input(),before=JSON.stringify(a),v=sangerAnnotationViews([a],20)[0];expect(v.rows.map(r=>r.fraction)).toEqual([.6,.25]);expect(v.hiddenFraction).toBeCloseTo(.15);expect(JSON.stringify(a)).toBe(before);});
 it('displays low-share fitted groups with the uncertainty flag retained',()=>{const a={...input(),share:.01,estimable:false},v=sangerAnnotationViews([a],0)[0];expect(v.rows).toHaveLength(4);expect(v.estimable).toBe(false);});
});

it('renders both zero-net joint edits and suppresses arbitrary unresolved junctions',()=>{const a:any={name:'A',reference:'A'.repeat(300),cut:80,cuts:[80,180],direction:'forward',share:1,estimable:true,components:[{size:0,fraction:.5,events:[{target:1,start:80,end:81,inserted:0},{target:2,start:180,end:180,inserted:1}]},{size:-100,fraction:.5,events:[],unresolved:true}]};const view=sangerAnnotationViews([a],0)[0];expect(view.cutIndices).toHaveLength(2);expect(view.rows[0].classification).toBe('Net-zero multi-target edit');expect(view.rows[0].tokens.some(t=>t.type==='insert')).toBe(true);expect(view.rows[0].tokens.some(t=>t.type==='delete')).toBe(true);expect(view.rows[1].classification).toBe('Breakpoint unresolved');expect(view.rows[1].tokens.every(t=>t.type==='equal')).toBe(true);});

it('shows QC-flagged estimates by default and allows optional hiding without altering the fit',()=>{const a={name:'A',reference:'ACGT'.repeat(50),cut:70,direction:'forward' as const,share:.4,estimable:false,components:[{size:-1,fraction:1}]};expect(sangerAnnotationViews([a],2)[0].rows).toHaveLength(1);expect(sangerAnnotationViews([a],2,30,false)[0].rows).toHaveLength(0);expect(sangerAnnotationViews([{...a,share:.001}],2)[0].rows).toHaveLength(1);expect(sangerAnnotationViews([{...a,share:0}],2)[0].rows).toHaveLength(0);expect(a.components[0].fraction).toBe(1);});

it('groups same-length placements before display thresholding and preserves all allocations',()=>{
 const a={...input(),components:[{size:0,fraction:.97},{size:-2,fraction:.012,events:[{target:1,start:4,end:6,inserted:0}]},{size:-2,fraction:.018,events:[{target:1,start:7,end:9,inserted:0}]}]},before=JSON.stringify(a),row=sangerAnnotationViews([a],2,8)[0].rows.find(r=>r.size===-2)!;
 expect(row.fraction).toBeCloseTo(.03);expect(row.placements).toHaveLength(2);expect(row.placements.reduce((s,p)=>s+p.fraction,0)).toBeCloseTo(.03);expect(row.tokens.find(t=>t.type==='delete')?.val).toBe(a.reference.slice(7,9));expect(row.positionLabel).toBe('Best fit: Deletion bases 8–9');expect(JSON.stringify(a)).toBe(before);
});
it('lowering the display threshold never upgrades unresolved length evidence to a real edit',()=>{
 const support={status:'unresolved' as const,score:.5,requiredScore:4,scope:'reference-length' as const},a={...input(),components:[{size:0,fraction:.995},{size:1,fraction:.005,support,events:[{target:1,start:6,end:6,inserted:1}]}]};expect(sangerAnnotationViews([a],2)[0].rows).toHaveLength(1);const row=sangerAnnotationViews([a],0)[0].rows[1];expect(row.support).toEqual(support);expect(row.positionLabel).toContain('Best fit:');expect(row.placements[0].tokens.find(t=>t.type==='insert')?.val).toBe('N');
});

it('shows291 as the representative site while retaining the whole25% length fraction and source alternatives',()=>{
 const a:SangerAnnotation={name:'A',reference:'ACGT'.repeat(100),cut:290,direction:'forward',share:1,estimable:true,components:[{size:-1,fraction:.2498,events:[{target:1,start:290,end:291,inserted:0}]},{size:-1,fraction:.0001,events:[{target:1,start:293,end:294,inserted:0}]},{size:-1,fraction:.0001,events:[{target:1,start:295,end:296,inserted:0}]}]},before=JSON.stringify(a),row=sangerAnnotationViews([a],0,30)[0].rows[0];
 expect(row.positionLabel).toBe('Best fit: Deletion base 291');expect(row.fraction).toBeCloseTo(.25);expect(row.representative!.events[0].start).toBe(290);expect(row.tokens.filter(t=>t.type==='delete').map(t=>t.val).join('')).toHaveLength(1);expect(row.placements).toHaveLength(3);expect(JSON.stringify(a)).toBe(before);
});
it('chooses a deterministic representative on ties and reports the tie',()=>{
 const components=[{size:-1,fraction:.5,events:[{target:1,start:7,end:8,inserted:0}]},{size:-1,fraction:.5,events:[{target:1,start:4,end:5,inserted:0}]}];for(const cs of[components,[...components].reverse()]){const row=sangerAnnotationViews([{...input(),components:cs}],0)[0].rows[0];expect(row.positionLabel).toBe('Best fit: Deletion base 4 (tied best fit)');}
});
