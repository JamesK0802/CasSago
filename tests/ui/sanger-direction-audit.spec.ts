import {it,expect} from 'vitest';
import {writeFileSync} from 'node:fs';
import {sangerDesignExample} from '../../src/app/utils/sanger-design.utils';
import {exactBindingSites,designSangerPrimers} from '../../src/app/workers/core/sanger-primer-designer';
import {SANGER_DEFAULTS} from '../../src/app/models/sanger-design';
it('audits demo direction geometry without changing filters',()=>{const refs=sangerDesignExample(),sequence='CAGCGACTCCTTGAGTGACG';const offsets=refs.map(r=>{const h=exactBindingSites(sequence,r).find(h=>h.direction==='reverse')!;return h.start-r.cuts[0];});expect(new Set(offsets).size).toBe(1);const result=designSangerPrimers(refs,{...SANGER_DEFAULTS,readEnd:1000,primerMin:20,primerMax:20});const report={reverseSharedPrimer:sequence,reverseOffsets:offsets,minimumReverseOffset:0,returned:result.candidates.map(c=>({direction:c.direction,sequence:c.sequence,minOffset:c.minOffset})),stats:result.stats};writeFileSync('../direction-audit.json',JSON.stringify(report,null,2));expect(result.candidates.length).toBeGreaterThan(0);});
