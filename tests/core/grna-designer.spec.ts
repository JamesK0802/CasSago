import {describe,it,expect} from 'vitest';
import {scanSpCas9,designGuides,validateGuideInput} from '../../src/app/workers/core/grna-designer';
import {GUIDE_DEFAULTS,GuideReference} from '../../src/app/models/grna-design';
import {reverseComplement} from '../../src/app/utils/biology.utils';
const g='GACCTAGTGCATGACCTGAC';const ref=(name:string,sequence:string,selected=true):GuideReference=>({name,sequence,selected,regionStart:null,regionEnd:null});
const fixture=(spacer=g,pam='AGG')=>'N'.repeat(50)+spacer+pam+'N'.repeat(50);
describe('SpCas9 guide candidate geometry and supplied-reference screen',()=>{
 it('reports correct input-strand coordinates and boundary on both orientations',()=>{
  const plus=scanSpCas9(fixture(),0).find(h=>h.sequence===g)!;expect(plus).toMatchObject({start:50,end:70,pamStart:70,pamEnd:73,strand:'+',cut:67,pam:'AGG'});
  const minus=scanSpCas9('N'.repeat(50)+reverseComplement(g+'AGG')+'N'.repeat(50),0).find(h=>h.sequence===g)!;
  expect(minus).toMatchObject({start:53,end:73,pamStart:50,pamEnd:53,strand:'-',cut:56,pam:'AGG'});
 });
 it('finds an identical guide with different PAMs and orientations in selected references',()=>{
  const r=[ref('A',fixture()),ref('B',fixture(g,'TGG')),ref('D',reverseComplement(fixture()))];const out=designGuides(r,[],{...GUIDE_DEFAULTS,mode:'shared'}),c=out.candidates.find(c=>c.sequence===g)!;
  expect(c.sites).toHaveLength(3);expect(c.otherHits).toEqual([0,0,0]);expect(out.complete).toBe(true);
 });
 it('counts near sites with canonical PAM using disjoint-seed retrieval, without counting intended sites',()=>{
  const mutated=(n:number)=>[...g].map((c,i)=>i<n?({A:'C',C:'G',G:'T',T:'A'}[c]!):c).join('');
  const bg=[0,1,2,3,4].map(n=>fixture(mutated(n))).join('N'.repeat(30))+fixture(mutated(1),'AGA')+fixture(g.slice(0,4)+'N'+g.slice(5));
  for(let m=0;m<=3;m++){const out=designGuides([ref('Target',fixture())],[{name:'Background',sequence:bg}],{...GUIDE_DEFAULTS,mismatches:m}),c=out.candidates.find(c=>c.sequence===g)!;expect(c.otherHits).toEqual(Array(m+1).fill(1));expect(c.exampleOtherSites.every(h=>h.reference===1)).toBe(true);}
 });
 it('treats unselected references as backgrounds and reference-specific mode rejects cross-binding',()=>{
  const refs=[ref('A',fixture()),ref('B',fixture(),false)];const out=designGuides(refs,[],GUIDE_DEFAULTS);expect(out.candidates.find(c=>c.sequence===g)!.otherHits[0]).toBe(1);
  refs[1].selected=true;expect(designGuides(refs,[],{...GUIDE_DEFAULTS,mode:'specific'}).candidates.some(c=>c.sequence===g)).toBe(false);
 });
 it('returns positive reference-specific candidates when the input targets differ sufficiently',()=>{
  const other='CTGA'+g.slice(4),out=designGuides([ref('A',fixture()),ref('B',fixture(other))],[],{...GUIDE_DEFAULTS,mode:'specific'});
  expect(out.candidates.some(c=>c.sequence===g&&c.sites[0].reference===0)).toBe(true);
  expect(out.candidates.some(c=>c.sequence===other&&c.sites[0].reference===1)).toBe(true);
  expect(out.candidates.every(c=>c.otherHits.every(n=>n===0))).toBe(true);
 });
 it('requires exactly one selected target in single mode and two in shared mode',()=>{
  expect(()=>validateGuideInput([ref('A',fixture()),ref('B',fixture())],[],GUIDE_DEFAULTS)).toThrow(/exactly one/);
  expect(()=>validateGuideInput([ref('A',fixture())],[],{...GUIDE_DEFAULTS,mode:'shared'})).toThrow(/at least two/);
 });
 it('filters cut region by boundary and reports bounded shortlist honestly',()=>{
  const r=ref('A',fixture());r.regionStart=67;r.regionEnd=67;expect(designGuides([r],[],GUIDE_DEFAULTS).candidates.some(c=>c.sequence===g)).toBe(true);
  r.regionStart=68;r.regionEnd=70;expect(designGuides([r],[],GUIDE_DEFAULTS).candidates.some(c=>c.sequence===g)).toBe(false);
  r.regionStart=r.regionEnd=null;r.sequence=fixture()+fixture('ACGTAGCTAGCTGACGTAGC');const out=designGuides([r],[],{...GUIDE_DEFAULTS,maxCandidates:1});expect(out.complete).toBe(false);expect(out.screened).toBe(1);expect(out.eligible).toBeGreaterThan(1);
 });
 it('never calls repeated exact targets unique or infers a cut from ambiguous spacer/PAM',()=>{
  const r=ref('A',fixture()+fixture());expect(designGuides([r],[],GUIDE_DEFAULTS).candidates.some(c=>c.sequence===g)).toBe(false);
  expect(scanSpCas9(fixture(g,'NGG'),0).some(h=>h.sequence===g)).toBe(false);
 });
 it('filters poly-T only when requested and checks invalid input bounds',()=>{
  const r=ref('A',fixture('TTTTGCGCGCGCACGTAGCA'));const s={...GUIDE_DEFAULTS,gcMin:0,gcMax:100};expect(designGuides([r],[],s).candidates.some(c=>c.polyT)).toBe(true);expect(designGuides([r],[],{...s,excludePolyT:true}).candidates.some(c=>c.polyT)).toBe(false);
  r.regionEnd=999;expect(()=>validateGuideInput([r],[],s)).toThrow(/bounds/);r.regionEnd=null;r.sequence='ACGTX';expect(()=>validateGuideInput([r],[],s)).toThrow(/DNA/);
 });
});
