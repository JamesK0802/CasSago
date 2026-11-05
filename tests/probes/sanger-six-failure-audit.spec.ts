import {it} from 'vitest';
import {readFileSync,writeFileSync} from 'node:fs';
import {benchmarkCompatibilityDiagnostics,parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
(process.env.CASSAGO_RESEARCH_AUDIT==='1'?it:it.skip)('audits saved fixed fit without changing any inferred calls',()=>{
const path='../benchmarks/point-display-20261104/',saved=JSON.parse(readFileSync(path+'result-final.json','utf8'));
const data=parseBenchmarkTruth(saved.syntheticTruth),refs=saved.input.references.map((r:any)=>({name:r.name,sequence:r.sequence,cuts:r.cuts}));
const start=Date.now(),rows=benchmarkCompatibilityDiagnostics(data,refs,saved.result);
const certificates=[...new Set(rows.map(r=>r.reference))].map(reference=>({reference,certificate:rows.find(r=>r.reference===reference)?.profileCertificate??null}));
writeFileSync(path+'six-failure-whole-profile-audit.json',JSON.stringify({elapsedMs:Date.now()-start,certificates,note:'Post-hoc compatibility with a fixed blind fit, allows 2 pp marginal tolerance and up to 1 pp off-profile mass. Not genotype recovery or proof of uniqueness.'},null,2));
console.log(JSON.stringify(certificates));
},300000);
