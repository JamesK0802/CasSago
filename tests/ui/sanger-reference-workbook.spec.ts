import {expect,it} from 'vitest';
import {referenceWorkbook,readReferenceWorkbook} from '../../src/app/utils/sanger-reference-workbook';
it('round-trips CasMango-compatible Excel reference configuration',async()=>{const refs=[{name:'Gene A',sequence:'ACGT'.repeat(50),targets:'ACGTACGTACGTACGTACGT',mode:'guides' as const,uiTargets:[{name:'Target 1',value:'ACGTACGTACGTACGTACGT'}]}],buffer=await referenceWorkbook(refs,true),loaded=await readReferenceWorkbook(buffer as unknown as ArrayBuffer);expect(loaded).toEqual(refs);});
