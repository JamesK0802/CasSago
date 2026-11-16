import {expect,it} from 'vitest';
import {assertDesignedCutSpacing} from '../probes/designed-cut-spacing';
const primer='ACGTCGTAGCTAGGATCCGA';
const ref=(name:string,cuts:number[])=>({name,sequence:primer+'A'.repeat(800),cuts:cuts.map(c=>c+primer.length)});
it('accepts the inclusive design boundary and records every reference/target coordinate',()=>{
 const result=assertDesignedCutSpacing([ref('G1',[150,400]),ref('G2',[180,430])],primer);
 expect(result.minimum).toBe(30);expect(result.cuts.map(c=>c.offset)).toEqual([150,180,400,430]);
});
it('rejects a cross-target collision before invoking the analysis engine',()=>{
 expect(()=>assertDesignedCutSpacing([ref('G1',[150,400]),ref('G2',[400,650])],primer)).toThrow(/minimum is 0/);
 expect(()=>assertDesignedCutSpacing([ref('G1',[150,400]),ref('G2',[429,650])],primer)).toThrow(/minimum is 29/);
 expect(()=>assertDesignedCutSpacing([ref('G1',[150,170]),ref('G2',[300])],primer)).toThrow(/minimum is 20/);
});
