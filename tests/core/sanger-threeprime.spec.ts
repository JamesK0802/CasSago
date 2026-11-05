import { describe,it,expect } from 'vitest';
import { SANGER_DEFAULTS } from '../../src/app/models/sanger-design';
import { sangerDesignExample } from '../../src/app/utils/sanger-design.utils';
import { reverseComplement } from '../../src/app/utils/biology.utils';
import { designSangerPrimers,hasTerminalCompatibleSite,homoeologVariantMaps } from '../../src/app/workers/core/sanger-primer-designer';
const settings={...SANGER_DEFAULTS,workflow:'homeolog-pcr' as const,specificityStrategy:'three-prime' as const,primerMin:20,primerMax:20,tmMin:48,tmMax:72};
const refs=()=>{const original=sangerDesignExample()[0],at=original.sequence.indexOf('CCGGCGGGTACGACCATAAC')+19;return ['A','B','D'].map((name,i)=>({...original,name,sequence:original.sequence.slice(0,at)+'CTG'[i]+original.sequence.slice(at+1)}));};
describe('aligned 3-prime homeolog PCR design',()=>{
  it('deliberately ranks diagnostic 3-prime SNPs, allowing a shared opposite primer',()=>{
    const r=refs(),out=designSangerPrimers(r,settings);expect(out.candidates.length).toBeGreaterThan(0);
    const c=out.candidates.find(c=>c.referenceIndex===0&&c.sequence==='CCGGCGGGTACGACCATAAC')!;expect(c).toBeDefined();expect(c.specificity!.primary).toBe(true);expect(c.specificity!.strategy).toBe('three-prime');expect(c.specificity!.primaryVariants).toHaveLength(2);
    for(const v of c.specificity!.primaryVariants!){expect(v.fromThreePrime).toBe(1);expect(v.targetBase).toBe('C');expect(v.backgroundBase).not.toBe('C');}
    expect(c.partner).not.toBeNull();
  });
  it('reports reverse-primer 3-prime bases in oligo orientation',()=>{
    const r=refs().map(v=>({...v,sequence:reverseComplement(v.sequence),cuts:v.cuts.map(c=>v.sequence.length-c)}));
    const out=designSangerPrimers(r,{...settings,direction:'reverse'});const c=out.candidates.find(c=>c.referenceIndex===0&&c.sequence==='CCGGCGGGTACGACCATAAC')!;
    expect(c).toBeDefined();expect(c.specificity!.primaryVariants![0].targetBase).toBe('C');expect(c.specificity!.primaryVariants![0].fromThreePrime).toBe(1);
  });
  it('does not call a 5-prime-only SNP a diagnostic 3-prime end',()=>{
    const original=sangerDesignExample()[0],at=original.sequence.indexOf('CCGGCGGGTACGACCATAAC');const r=['C','A','T'].map((base,i)=>({...original,name:'Ref'+i,sequence:original.sequence.slice(0,at)+base+original.sequence.slice(at+1)}));
    const out=designSangerPrimers(r,settings);expect(out.candidates.some(c=>c.referenceIndex===0&&c.sequence==='CCGGCGGGTACGACCATAAC')).toBe(false);
  });
  it('screens additional near sites with compatible 3-prime ends',()=>{
    const oligo='CCGGCGGGTACGACCATAAC',fivePrimeMismatch='ACGGCGGGTACGACCATAAC',terminalMismatch='CCGGCGGGTACGACCATAAT';
    expect(hasTerminalCompatibleSite(oligo,[fivePrimeMismatch],2,3)).toBe(true);expect(hasTerminalCompatibleSite(oligo,[terminalMismatch],2,3)).toBe(false);
    const r=refs();r[1].sequence+=fivePrimeMismatch;expect(designSangerPrimers(r,settings).candidates.some(c=>c.referenceIndex===0&&c.sequence===oligo)).toBe(false);
  });
  it('requires evidence against every other homeolog at the target-anchored alignment',()=>{
    const r=refs(),maps=homoeologVariantMaps(r,0);expect(maps).toHaveLength(2);expect(maps.every(m=>m.identity>.99)).toBe(true);
    r[2].sequence=r[0].sequence;const out=designSangerPrimers(r,settings);expect(out.candidates.some(c=>c.referenceIndex===0)).toBe(false);
  });
});
