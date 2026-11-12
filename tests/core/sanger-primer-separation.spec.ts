import {it,expect} from 'vitest';import {screenSharedSangerWT} from '../../src/app/workers/core/sanger-primer-separation';import {PrimerSite} from '../../src/app/models/sanger-design';
import {minOffset} from '../../src/app/workers/core/sanger-primer-designer';
const site=(extension:string,cut:number):PrimerSite=>({reference:String(cut),start:0,end:20,direction:'forward',cutOffsets:[cut],extension});
const cuts=(...offsets:number[]):PrimerSite=>({...site('ACGT'.repeat(200),offsets[0]),cutOffsets:offsets});
it('compares different target indices across genes on the same read axis',()=>{
 expect(minOffset([cuts(150,400),cuts(400,650)])).toBe(0);
 expect(minOffset([cuts(150,400),cuts(425,650)])).toBe(25);
 expect(minOffset([cuts(150,400),cuts(430,650)])).toBe(30);
});
it('also compares targets within a reference and unequal target counts',()=>{
 expect(minOffset([cuts(150,170),cuts(300)])).toBe(20);
 expect(minOffset([cuts(150),cuts(300,160,500)])).toBe(10);
 expect(minOffset([cuts(300,150)])).toBe(150);
 expect(minOffset([cuts(150)])).toBe(0);
});
it('is independent of reference order and reverse-read target ordering',()=>{
 expect(minOffset([cuts(400,150),cuts(455,205),cuts(510,260)])).toBe(55);
 expect(minOffset([cuts(510,260),cuts(455,205),cuts(400,150)])).toBe(55);
});
it('excludes identical observable products even when target offsets differ',()=>{const s='ACGT'.repeat(150),r=screenSharedSangerWT([site(s,150),site(s,220)],{readStart:35,readEnd:500});expect(r.passed).toBe(false);expect(r.reason).toMatch(/identical/);});
it('does not count differences outside the sequencing range as usable evidence',()=>{const prefix='ACGT'.repeat(120),r=screenSharedSangerWT([site(prefix+'AAAA',150),site(prefix+'CCCC',220)],{readStart:35,readEnd:400});expect(r.passed).toBe(false);});
it('excludes recombinant WT dictionaries with an explicit mixture nullspace',()=>{const prefix='A'.repeat(50),suffix='G'.repeat(100),sites=['AC','AT','GC','GT'].map((x,i)=>site(prefix+x[0]+'A'.repeat(20)+x[1]+suffix,100+i*30));expect(screenSharedSangerWT(sites,{readStart:35,readEnd:170}).passed).toBe(false);});
it('retains distinguishable products without promising edited-allele identifiability',()=>{const a='ACGT'.repeat(150),b=a.slice(0,100)+'TGCATG'+a.slice(106),r=screenSharedSangerWT([site(a,150),site(b,220)],{readStart:35,readEnd:500});expect(r.passed).toBe(true);expect(r.minimumDifferences).toBeGreaterThan(0);expect(r.reason).toMatch(/not guaranteed/);});
it('removes rejected products from the actual shared-primer candidate list and records their reasons',async()=>{const {designSangerPrimers}=await import('../../src/app/workers/core/sanger-primer-designer'),{sangerDesignExample}=await import('../../src/app/utils/sanger-design.utils'),{SANGER_DEFAULTS}=await import('../../src/app/models/sanger-design');const original=sangerDesignExample()[0],result=designSangerPrimers([original,{...original,name:'Duplicate',cuts:original.cuts.map(c=>c+50)}],{...SANGER_DEFAULTS,primerMin:20,primerMax:20,requirePartner:false});expect(result.candidates).toHaveLength(0);expect(result.stats.separabilityRejected).toBeGreaterThan(0);expect(result.excludedCandidates!.every(x=>!!x.reason)).toBe(true);},30000);
