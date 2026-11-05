import {expect,it} from 'vitest';
import {substitutionLinkageAliases} from '../../src/app/workers/core/sanger-substitution-linkage';
const reference='ACGTACGTACGTACGT',sequence='ACGAACGTTCGTACGT';
const column=(s:string)=>s.split('').flatMap(c=>'ACGT'.split('').map(b=>Number(b===c)));
it('detects the exact WT90/linked10 versus WT80/separate10+10 signal ambiguity without truth',()=>{
 const refs=[{sequence:reference,cuts:[8]}],chosen=[{reference:0,sequence:reference,events:[],fraction:.9},{reference:0,sequence,events:[{start:3,end:9,inserted:sequence.slice(3,9)}],fraction:.1}];
 const aliases=substitutionLinkageAliases(refs,chosen,column);expect(aliases).toHaveLength(1);expect(aliases[0].blocks).toBe(2);expect(aliases[0].wtRange[0]).toBeCloseTo(.8);expect(aliases[0].wtRange[1]).toBeCloseTo(.9);
});
it('also detects adjacent equal-frequency substitutions when emission columns are identical',()=>{
 const s=reference.slice(0,3)+'GG'+reference.slice(5);expect(substitutionLinkageAliases([{sequence:reference,cuts:[8]}],[{reference:0,sequence:reference,events:[],fraction:.9},{reference:0,sequence:s,events:[{start:3,end:5,inserted:'GG'}],fraction:.1}],column)).toHaveLength(1);
});
it('rejects an alias when context-sensitive emissions distinguish the recipes',()=>{
 const context=(s:string)=>[...column(s),Number(s===sequence)];expect(substitutionLinkageAliases([{sequence:reference,cuts:[8]}],[{reference:0,sequence:reference,events:[],fraction:.9},{reference:0,sequence,events:[{start:3,end:9,inserted:sequence.slice(3,9)}],fraction:.1}],context)).toEqual([]);
});
it('does not call a one-bit context difference an exact alias',()=>{
 const context=(s:string)=>[...column(s),s===sequence?Number.MIN_VALUE:0];expect(substitutionLinkageAliases([{sequence:reference,cuts:[8]}],[{reference:0,sequence:reference,events:[],fraction:.9},{reference:0,sequence,events:[{start:3,end:9,inserted:sequence.slice(3,9)}],fraction:.1}],context)).toEqual([]);
});
