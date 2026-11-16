import '@angular/compiler';
import {it,expect} from 'vitest';
import {SangerAnnotationComponent} from '../../src/app/components/sanger-annotation/sanger-annotation.component';
import {exactDisplayView} from '../../src/app/utils/sanger-exact-display';
import {marginalTruthStatus,marginalSequenceResolved} from '../../src/app/utils/sanger-marginal-display-status';
import type {SangerAnnotation,SangerMarginalGroup} from '../../src/app/models/sanger-analysis';

const group:SangerMarginalGroup={id:'family',reference:'ref',target:1,size:2,fraction:.1234,fractionRange:[0,1],status:'estimated',scope:'read-cycle-by-length',boundMethod:'deferred-report-envelope',representative:{events:[{target:1,start:10,end:10,inserted:2,insertedSequence:'CG'}],witnessFraction:.07,source:'fitted-marginal-proposal'},bases:[{position:9,estimate:[0,1,0,0],lower:[0,0,0,0],upper:[1,1,1,1]},{position:10,estimate:[.2,.8,0,0],lower:[0,0,0,0],upper:[1,1,1,1]}]};
const annotation:SangerAnnotation={name:'ref',reference:'ACGT'.repeat(10),cut:10,cuts:[10],direction:'forward',share:1,estimable:true,components:[],marginalGroups:[group]};

it('shows an actual representative and point percentages without presenting uncomputed intervals or a uniqueness claim',()=>{
 const row=exactDisplayView(annotation,0,10).rows[0],ui=new SangerAnnotationComponent();
 expect(row.tokens).toContainEqual({type:'insert',val:'CG'});expect(row.fractionRange).toBeUndefined();expect(ui.fractionText(row)).toBe('12.34');
 expect(ui.ambiguousBases(row).map(b=>b.position)).toEqual([10]);expect(ui.baseFitText(.8)).toBe('80.0%');expect(marginalSequenceResolved(group)).toBe(false);
 const truth={...annotation,components:[{size:2,fraction:.1234,events:group.representative!.events}]};
 expect(marginalTruthStatus(truth,group,row).status).toBe('mismatch');
 const requested:any[]=[];ui.boundsRequested.subscribe(v=>requested.push(v));ui.requestBounds(row);expect(requested).toEqual([{groupId:'family'}]);
});

it('preserves an estimated single-site WT/Sub partition and never promotes its placeholder ranges when parent bounds arrive',()=>{
 const family:SangerMarginalGroup={...group,size:0,fraction:1,displayComponents:[{id:'family::base:1',rangeStatus:'estimated',fraction:.8,fractionRange:[0,1],events:[],bases:group.bases,mixture:false,selection:{marker:{position:10,bases:[1]}}},{id:'family::base:2',rangeStatus:'estimated',fraction:.2,fractionRange:[0,1],events:[{target:1,start:10,end:11,inserted:1,insertedSequence:'G'}],bases:group.bases,mixture:false,selection:{marker:{position:10,bases:[2]}}}]};
 for(const status of ['estimated','bounded']as const){const rows=exactDisplayView({...annotation,marginalGroups:[{...family,status}]},0,10).rows,ui=new SangerAnnotationComponent();expect(rows.map(r=>ui.badge(r))).toEqual(['WT','Substitution']);expect(rows.map(r=>ui.fractionText(r))).toEqual(['80.00','20.00']);expect(rows.every(r=>r.marginal!.status==='estimated'&&r.fractionRange===undefined)).toBe(true);}
});

it('does not hide a shared allocation merely because its fitted base happens to be 100 percent',()=>{
 const shared={...group,bases:[{position:10,estimate:[1,0,0,0],lower:[0,0,0,0],upper:[1,1,1,1],sharedWith:['other-family']}]},row=exactDisplayView({...annotation,marginalGroups:[shared]},0,10).rows[0],ui=new SangerAnnotationComponent();
 expect(ui.ambiguousBases(row)).toHaveLength(1);expect(ui.baseFitText(row.marginal!.bases[0].estimate[0])).toBe('100.0%');expect(ui.sharedBaseTitle(row.marginal!.bases[0])).toContain('other-family');expect(row.fractionRange).toBeUndefined();expect(marginalSequenceResolved(shared)).toBe(false);
 const bounded={...shared,status:'bounded' as const,bases:[{...shared.bases[0],lower:[1,0,0,0],upper:[1,0,0,0]}]},boundedRow=exactDisplayView({...annotation,marginalGroups:[bounded]},0,10).rows[0];expect(ui.ambiguousBases(boundedRow)).toHaveLength(0);
});
