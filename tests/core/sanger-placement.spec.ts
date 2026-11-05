import{it,expect}from'vitest';
import{canonicalSangerEvents}from'../../src/app/utils/sanger-placement';
import{sangerAnnotationViews}from'../../src/app/utils/sanger-annotation';
it('combines equivalent deletion mass before selecting the most likely edited molecule',()=>{
 const reference='CCAAAAGT',truth=[2,3,4].map(start=>({size:-1,fraction:.18,events:[{target:1,start,end:start+1,inserted:0}]})),wrong={size:-1,fraction:.46,events:[{target:1,start:6,end:7,inserted:0}]},a={name:'A',reference,cut:4,direction:'forward' as const,share:1,estimable:true,components:[...truth,wrong]},before=JSON.stringify(a),row=sangerAnnotationViews([a],0,30)[0].rows[0];
 expect(row.representative!.events[0].start).toBe(2);expect(row.representative!.fraction).toBeCloseTo(.54);expect(row.fraction).toBeCloseTo(.54);expect(row.placements).toHaveLength(1);expect(sangerAnnotationViews([a],0,30)[0].rows).toHaveLength(2);expect(JSON.stringify(a)).toBe(before);
});
it('canonicalizes only edits preserving the final DNA, including known inserted-base rotation',()=>{
 const reference='CCATATGG',apply=(e:any)=>reference.slice(0,e.start)+(e.insertedSequence||'')+reference.slice(e.end);
 for(const e of[{target:1,start:4,end:6,inserted:0},{target:1,start:4,end:4,inserted:2,insertedSequence:'AT'}]){const normalized=canonicalSangerEvents(reference,[e])[0];expect(apply(normalized)).toBe(apply(e));}
 const unknown={target:1,start:4,end:4,inserted:2,insertedSequence:'NN'};expect(canonicalSangerEvents(reference,[unknown])[0]).toEqual(unknown);
});
