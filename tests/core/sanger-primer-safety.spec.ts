import {expect,it} from 'vitest';
import {screenPrimerBinding,usablePrimerRead,comparePrimerSafety} from '../../src/app/workers/core/sanger-primer-safety';
import {reverseComplement} from '../../src/app/utils/biology.utils';
import {SANGER_DEFAULTS} from '../../src/app/models/sanger-design';
const oligo='ACCTGATCGTACCGTAGCTA';
it('excludes intended positions in every homeolog but detects a three-mismatch alternate site on the reverse strand',()=>{
 const near='TCCTGATCGCACCGTAGCTC',refs=[{name:'A',sequence:oligo+'GGGGGG',cuts:[20]},{name:'B',sequence:oligo+'GGGGGG'+reverseComplement(near),cuts:[20]}];
 const result=screenPrimerBinding(oligo,refs,refs.map(r=>({reference:r.name,start:0,end:20,direction:'forward'})));
 expect(result.minimumMismatches).toBe(3);expect(result.nearSiteCount).toBe(1);expect(result.sites[0]).toMatchObject({reference:'B',start:26,end:46,direction:'reverse',mismatches:3,terminalMismatches:1});
 expect(JSON.stringify(refs)).toContain(oligo);
});
it('does not exclude an unintended exact copy just because another site is intended',()=>{
 const ref={name:'B',sequence:oligo+'GGGGGG'+oligo,cuts:[20]};const result=screenPrimerBinding(oligo,[ref],[{reference:'B',start:0,end:20,direction:'forward'}]);expect(result.minimumMismatches).toBe(0);expect(result.sites.some(s=>s.start===26)).toBe(true);
});
it('agrees with exhaustive scanning on both strands across independent substitution patterns',()=>{
 let seed=92714;const random=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed;};
 for(let replicate=0;replicate<24;replicate++){
  let dna=Array.from({length:180},()=> 'ACGT'[random()>>>30]).join('');const mutated=[...oligo],positions=new Set<number>();while(positions.size<replicate%9)positions.add(random()%oligo.length);for(const p of positions)mutated[p]='ACGT'[('ACGT'.indexOf(mutated[p])+1)%4];dna=dna.slice(0,62)+mutated.join('')+dna.slice(82);
  const ref={name:'R',sequence:dna,cuts:[100]},expected:any[]=[],distances:number[]=[];
  for(const [direction,s]of [['forward',dna],['reverse',reverseComplement(dna)]]as const)for(let at=0;at<=s.length-oligo.length;at++){const mismatches=[...oligo].filter((b,j)=>b!==s[at+j]).length;distances.push(mismatches);if(mismatches<=3)expected.push({direction,start:direction==='forward'?at:s.length-at-oligo.length,mismatches});}
  const actual=screenPrimerBinding(oligo,[ref],[]);expect(actual.nearSiteCount).toBe(expected.length);expect(actual.minimumMismatches).toBe(Math.min(...distances));for(const row of expected)expect(actual.sites).toContainEqual(expect.objectContaining(row));
 }
});
it('reports unresolved N sequence instead of counting N as a known mismatch',()=>{const r=screenPrimerBinding(oligo,[{name:'R',sequence:oligo.slice(0,10)+'N'+oligo.slice(11),cuts:[10]}],[]);expect(r.references[0].hasUnknownBases).toBe(true);expect(r.sites).toHaveLength(0);});
it('reserves both nominal read ends including a shorter physical product',()=>{expect(usablePrimerRead({...SANGER_DEFAULTS,readStart:0,readEnd:700},520)).toEqual({start:50,end:470});expect(usablePrimerRead({...SANGER_DEFAULTS,readStart:100,readEnd:700},900)).toEqual({start:100,end:650});});
it('ranks fewer unintended matches ahead of target proximity without treating intended common binding as risk',()=>{const candidate=(m:number,n:number)=>({primerSafety:{primary:{minimumMismatches:m,nearSiteCount:n},partner:null}}as any);expect(comparePrimerSafety(candidate(4,0),candidate(1,1))).toBeLessThan(0);expect(comparePrimerSafety(candidate(3,1),candidate(3,5))).toBeLessThan(0);});
