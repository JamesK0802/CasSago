import '@angular/compiler';
import {it,expect} from 'vitest';
import {profileFixture} from '../fixtures/sanger-profile-fixture';
import {generateSangerSynth} from '../../src/app/utils/sanger-synthetic-benchmark';
import {analyzeAutomaticSanger} from '../../src/app/workers/core/sanger-automatic';
import {benchmarkComparison,parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
it('fits molecules with edits at both targets and evaluates their two marginal profiles separately',()=>{
 const f=profileFixture(),d=generateSangerSynth(f.result,f.candidate,f.refs,9033),result=analyzeAutomaticSanger({...d.input,groupRecovery:true,maxDeletion:3,maxInsertion:2}),rows=benchmarkComparison(parseBenchmarkTruth({truth:d.truth,recipe:d.recipe}),d.template.amplicons,result).filter(r=>r.scope==='group'&&(r.truth??0)>0);
 expect(result.diagnostics?.algorithmVersion).toBe('sanger-sequence-groups-v4');expect(result.annotations?.some(a=>a.components.some(c=>(c.events?.length||0)>1))).toBe(true);expect(rows).toHaveLength(4);for(const row of rows)expect(Math.abs(row.delta??100),row.group+' T'+row.target).toBeLessThan(2);
},120000);
