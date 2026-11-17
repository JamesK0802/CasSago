import{it,expect}from'vitest';
import{summarizeMarginalPoint}from'../../src/app/workers/core/sanger-marginal-point-report';
import{substitutionBlockDisplay}from'../../src/app/workers/core/sanger-substitution-block-display';
import{marginalGroupBounds}from'../../src/app/workers/core/sanger-marginal-group-bounds';
import{marginalSessions}from'../../src/app/workers/core/sanger-marginal-state';
import{marginalAlternativePage}from'../../src/app/workers/core/sanger-marginal-alternatives';
import loadCommon from'../probes/load-common-optimizer';
function fixture(amounts:number[],positions=amounts.map((_,i)=>97+i)){
 const point=[1,0,0,0],program:any={upper:[1,0,0,0],cost:[0,0,0,0],rows:[{terms:[[0,1]],lower:1,upper:1}]},bases=new Map(),flow=new Map(),departures=new Map();
 for(let i=0;i<amounts.length;i++){const w=point.length,s=w+1;point.push(1-amounts[i],amounts[i]);program.upper.push(1,1);program.cost.push(0,1);program.rows.push({terms:[[w,1],[s,1],[0,-1]],lower:0,upper:0},{terms:[[s,1]],lower:amounts[i],upper:amounts[i]});bases.set(positions[i],[[],[[w,1]],[[s,1]],[]]);flow.set(positions[i],{unchanged:[{base:1,column:w,cap:[[0,1]]}],changed:[1,2,s,3]});departures.set(positions[i],[[s,1]]);}
 const budgetRow=program.rows.length;program.rows.push({terms:program.cost.flatMap((c:number,i:number)=>c?[[i,c]]:[]),lower:0,upper:amounts.reduce((a,b)=>a+b,0)});program.dualHint={budgetRow,shadow:[]};
 const edge={reference:0,target:0,incoming:0,outgoing:0,size:0,start:100,mass:0,free:new Map()},raw={id:'zero',reference:0,target:0,size:0,mass:[[0,1]],bases,flow,departures,edges:[edge],fraction:1,fractionRange:[0,1]},ref={name:'R',sequence:'C'.repeat(220),offset:0,reverse:false,cuts:[100],input:{name:'R',sequence:'C'.repeat(220),cut:100,cuts:[100]}},target={program,groups:[raw],point,positions},model:any={program,refs:[ref],shares:[1],noise:0,drafts:[target]};return{point,program,raw,ref,target,model};
}
function conserved(group:any){for(const p of group.bases)for(let b=0;b<4;b++)expect(group.displayComponents.reduce((v:number,c:any)=>v+c.fraction*c.bases.find((q:any)=>q.position===p.position).estimate[b],0)).toBeCloseTo(group.fraction*p.estimate[b],10);}
it('splits the user example10% at-3..0 and20% at1..3 without a nested long10% profile',()=>{
 const f=fixture([.1,.1,.1,.1,.2,.2,.2]),before=JSON.stringify(f.program),g=summarizeMarginalPoint(f.model,f.point,2)[0],rows=g.displayComponents!;expect(rows.map(c=>c.fraction)).toEqual([.7,.1,.2]);expect(rows[0].events).toEqual([]);expect(rows[1].events.map(e=>e.start)).toEqual([97,98,99,100]);expect(rows[2].events.map(e=>e.start)).toEqual([101,102,103]);expect(rows.some(c=>c.events.length===7)).toBe(false);expect(g.unmodifiedEstimate).toBeCloseTo(.7,12);conserved(g);expect(JSON.stringify(f.program)).toBe(before);
});
it('merges only contiguous same-fraction sites and keeps separated same-fraction blocks distinct',()=>{
 const a=fixture([.1,.1,.1]),one=summarizeMarginalPoint(a.model,a.point,2)[0];expect(one.displayComponents!.map(c=>c.fraction)).toEqual([.9,.1]);const b=fixture([.1,.1,.1,.1],[97,98,102,103]),two=summarizeMarginalPoint(b.model,b.point,2)[0];expect(two.displayComponents!.map(c=>c.fraction)).toEqual([.8,.1,.1]);conserved(one);conserved(two);
});
it('retains a full mixture when required overlap prevents a disjoint representative partition',()=>{
 const f=fixture([.8,.8],[97,103]),g=summarizeMarginalPoint(f.model,f.point,2)[0];expect(g.displayComponents).toBeUndefined();expect(g.bases.map(p=>p.estimate[2])).toEqual([.8,.8]);
});
it('conserves small within-block fraction differences rather than rounding away their bases',()=>{
 const f=fixture([.1,.10001,.2,.2]),g=summarizeMarginalPoint({...f.model,noise:.001},f.point,2)[0];for(const[i,v]of[.69999,.10001,.2].entries())expect(g.displayComponents![i].fraction).toBeCloseTo(v,12);conserved(g);expect(g.displayComponents![1].bases[0].estimate[2]).toBeLessThan(1);
});
it('keeps reverse-orientation block markers on the original read axis',()=>{
 const f=fixture([.1,.1,.2,.2]),g=summarizeMarginalPoint(f.model,f.point,2)[0];const rows=substitutionBlockDisplay({...g,bases:g.bases.slice().reverse()},()=>1,p=>219-p,true,0,.02)!.components;for(const c of rows){if(c.selection.marker){expect(c.selection.marker.bases).toEqual([1]);expect([119,121]).toContain(c.selection.marker.position);}}conserved({...g,bases:g.bases.slice().reverse(),displayComponents:rows});
});
it('bounds WT as an all-WT subset and retains overlap in selected-block bounds and10-item alternatives',async()=>{
 const f=fixture([.1,.1,.1,.1,.2,.2,.2]),groups=summarizeMarginalPoint(f.model,f.point,2),parent=groups[0],wt=parent.displayComponents![0],block=parent.displayComponents![1],result:any={marginalGroups:groups,warnings:[]},before=JSON.stringify(f.program),o=await loadCommon();marginalSessions.set(result,{refs:[f.ref],shares:[1],targets:[f.target],optimizer:o,pages:new Map()}as any);
 expect(marginalGroupBounds(result,wt.id,()=>{},0).complete).toBe(false);expect(marginalGroupBounds(result,wt.id).complete).toBe(true);expect(wt.fractionRange[0]).toBeCloseTo(0,5);expect(wt.fractionRange[1]).toBeCloseTo(.8,5);expect(wt.bases.every(p=>p.lower[1]===1&&p.upper[1]===1)).toBe(true);
 expect(marginalGroupBounds(result,block.id).complete).toBe(true);const other=block.bases.find(p=>p.position===101)!;expect(other.lower[2]).toBe(0);expect(other.upper[2]).toBe(1);const page=marginalAlternativePage(result,wt.id,0);expect(page.candidates).toHaveLength(1);expect(page.candidates[0].events).toEqual([]);expect(page.exhausted).toBe(true);expect(JSON.stringify(f.program)).toBe(before);conserved(parent);marginalSessions.delete(result);
});

it('keeps raw-identical overlapping populations as genuine alternatives to the chosen block representatives',async()=>{
 const{readFileSync}=await import('node:fs'),{writeSynthTrace,synthAlleleSequence}=await import('../../src/app/utils/sanger-synthetic-benchmark'),{parseSangerTrace}=await import('../../src/app/utils/sanger-trace'),f=JSON.parse(readFileSync(new URL('../fixtures/substitution-blocks.json',import.meta.url),'utf8')),ref:any={name:'portable',sourceReference:'portable',sequence:f.reference,cuts:[f.cut],targetNames:['T1'],share:100,primerLength:f.primerLength,partnerLength:20,alleles:[]};
 const traces=new Map<string,ReturnType<typeof parseSangerTrace>>();for(const c of f.cases){const parts=c.parts.map((p:any)=>{const events=p.offsets.map((offset:number)=>{const at=f.cut+offset,before=f.reference[at];return{kind:'substitution',target:0,offset,selectionStart:at,length:1,before,sequence:'ACGT'[('ACGT'.indexOf(before)+1)%4]};}),dna=events.length?synthAlleleSequence(ref,events):f.reference;return{sequence:dna.slice(f.primerLength),weight:p.percent/100};});traces.set(c.id,parseSangerTrace(writeSynthTrace(parts,f.cycles,c.id),c.id));}
 for(const[a,b]of f.rawIdenticalPairs){const x=traces.get(a)!,y=traces.get(b)!;for(const base of['A','C','G','T']as const)expect(x.channels[base]).toEqual(y.channels[base]);expect(x.peaks).toEqual(y.peaks);expect(x.sequence).toBe(y.sequence);}
});
