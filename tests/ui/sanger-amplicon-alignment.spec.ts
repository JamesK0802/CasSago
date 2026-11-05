import '@angular/compiler';
import {it,expect} from 'vitest';
import {SangerAmpliconAlignmentComponent} from '../../src/app/components/sanger-amplicon-alignment/sanger-amplicon-alignment.component';
import {reverseComplement} from '../../src/app/utils/biology.utils';
function fixture(reverse=false){
 const primary='ACGTACGTACGTACGTACGA',partner='TTGCAATTGCAATTGCAATT';
 const refs=[60,90].map((middle,i)=>({name:'Ref'+i,sequence:'N'.repeat(30)+(reverse?partner:primary)+'G'.repeat(middle)+reverseComplement(reverse?primary:partner)+'N'.repeat(40),cuts:[30+20+middle/2]}));
 const sites=refs.map((r,i)=>({reference:r.name,start:reverse?r.sequence.length-60:30,end:reverse?r.sequence.length-40:50,direction:reverse?'reverse':'forward',cutOffsets:[],extension:''}));
 const partners=refs.map(r=>({reference:r.name,start:reverse?30:r.sequence.length-60,end:reverse?50:r.sequence.length-40,direction:reverse?'forward':'reverse',cutOffsets:[],extension:''}));
 const candidate:any={id:'one',sequence:primary,direction:reverse?'reverse':'forward',sites,partner:{sequence:partner,sites:partners}};
 const result:any={references:refs,settings:{window:20},candidates:[candidate]};
 const p=new SangerAmpliconAlignmentComponent();p.result=result;p.selected=candidate;p.ngOnChanges();return {p,result,candidate,primary,partner};
}
it.each([false,true])('trims to both PCR primers and aligns in sequencing orientation, reverse=%s',reverse=>{
 const {p,result,primary,partner}=fixture(reverse),before=JSON.stringify(result);
 expect(p.rows.map(r=>r.length)).toEqual([100,130]);
 for(const row of p.rows){const sequence=row.segments.map(s=>s.text).join('');expect(sequence.startsWith(primary)).toBe(true);expect(sequence.endsWith(reverseComplement(partner))).toBe(true);expect(sequence).not.toContain('N');expect(row.reversed).toBe(reverse);expect(row.padding).toBe(0);expect(row.primerEnd).toBe(20);expect(row.segments.some(s=>s.cut)).toBe(true);}
 expect(p.rows[1].cuts[0]-p.rows[0].cuts[0]).toBe(15);expect(JSON.stringify(result)).toBe(before);
});
it('updates all products together when a different candidate pair is selected',()=>{const {p,candidate}=fixture();p.selected={...candidate,id:'two',direction:'reverse',sequence:reverseComplement(candidate.partner.sequence),sites:candidate.partner.sites,partner:{...candidate.partner,sequence:reverseComplement(candidate.sequence),sites:candidate.sites}};p.ngOnChanges();expect(p.rows.every(r=>r.reversed)).toBe(true);p.selected=null;p.ngOnChanges();expect(p.rows).toEqual([]);});
it.each([false,true])('marks supplied gRNAs in the sequenced strand and keeps the DNA rows primer-aligned, reverse=%s',reverse=>{const {p}=fixture(reverse);p.references=p.result!.references.map(r=>({name:r.name,targets:'G'.repeat(20),mode:'guides'}));p.ngOnChanges();expect(p.rows.every(r=>r.guides.length>0&&r.segments.some(s=>s.kind==='guide'))).toBe(true);expect(p.rows[0].padding+p.rows[0].primerEnd).toBe(p.rows[1].padding+p.rows[1].primerEnd);expect(p.width).toBe(Math.max(...p.rows.map(r=>r.padding+r.length)));});
