import '@angular/compiler';
import {expect,it} from 'vitest';
import {fiveGroupIdentifiabilityFixture,fiveIndelSequenceIdentifiabilityFixture} from '../fixtures/sanger-identifiability-fixture';
import {generateSangerSynth} from '../../src/app/utils/sanger-synthetic-benchmark';
it('cannot uniquely identify five groups and their percentages from identical input traces',()=>{
 const f=fiveGroupIdentifiabilityFixture();for(const recipe of [f.refs,f.alternate]){const frequencies=recipe.flatMap(r=>r.profiles!.flatMap(p=>p.alleles.map(a=>a.percent)));expect(new Set(frequencies).size).toBe(30);expect(recipe).toHaveLength(3);expect(recipe.every(r=>r.profiles!.length===2&&r.profiles!.every(p=>p.alleles.length===5&&p.alleles.every(a=>a.percent>0)))).toBe(true);}
 for(const seed of [20261104,20261105,20261106]){const a=generateSangerSynth(f.result,f.candidate,f.refs,seed),b=generateSangerSynth(f.result,f.candidate,f.alternate,seed);expect(a.recipe).not.toEqual(b.recipe);expect(a.template).toEqual(b.template);expect(new Uint8Array(a.control)).toEqual(new Uint8Array(b.control));expect(new Uint8Array(a.edit)).toEqual(new Uint8Array(b.edit));}
});

it('admits disjoint insertion sequences under the same five-indel-group observation',()=>{const f=fiveIndelSequenceIdentifiabilityFixture(),a=generateSangerSynth(f.result,f.candidate,f.refs,20261104),b=generateSangerSynth(f.result,f.candidate,f.alternate,20261104);expect(new Uint8Array(a.edit)).toEqual(new Uint8Array(b.edit));const strings=(recipe:typeof f.refs)=>recipe[0].profiles![0].alleles.filter(g=>g.events[0].kind==='insertion').map(g=>g.events[0].sequence);expect(strings(f.refs).some(s=>strings(f.alternate).includes(s))).toBe(false);for(const recipe of [f.refs,f.alternate]){const values=recipe.flatMap(r=>r.profiles!.flatMap(p=>p.alleles.map(a=>a.percent)));expect(new Set(values).size).toBe(30);}});

import {gramCondition} from '../../src/app/workers/core/sanger-signal';
it('does not mistake full rank of a selected basis for uniqueness across alternative supports',()=>{
 const encode=(s:string)=>s.split('').flatMap(c=>'ACGT'.split('').map(b=>Number(b===c))),a=['AAA','ACC','CAC','CCA'].map(encode),b=['AAC','ACA','CAA','CCC'].map(encode);
 expect(gramCondition(a).rank).toBe(4);expect(gramCondition(b).rank).toBe(4);expect(gramCondition([...a,...b]).rank).toBe(4);
});
