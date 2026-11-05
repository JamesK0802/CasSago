import '@angular/compiler';import{it,expect}from'vitest';import{readFileSync}from'node:fs';
import loadHighs from'../../public/vendor/highs-1.15.4/highs.mjs';
import{parseSangerTrace}from'../../src/app/utils/sanger-trace';import{analyzeConditionalGroups}from'../../src/app/workers/core/sanger-conditional-groups';import{automaticSangerSettings}from'../../src/app/workers/core/sanger-automatic';import{parseBenchmarkTruth}from'../../src/app/utils/sanger-benchmark';import{assessSangerBenchmark}from'../../src/app/utils/sanger-benchmark-assessment';import{certifyNativeBenchmark}from'../../src/app/utils/sanger-native-signal-certificate';
it('directly recovers the near-equal substitution/insertion failure without truth-derived inference inputs',async()=>{
 const base=new URL('../fixtures/release-v42/near-frequency/',import.meta.url),template=JSON.parse(readFileSync(new URL('template.json',base),'utf8')),bytes=(n:string)=>{const b=readFileSync(new URL(n,base));return b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength);},controlBytes=bytes('control.ab1'),sampleBytes=bytes('edit.ab1'),input={references:template.amplicons.map((r:any)=>({...r,cut:r.cuts[0]})),primer:template.primer,control:parseSangerTrace(controlBytes,'control'),sample:parseSangerTrace(sampleBytes,'edit'),maxDeletion:10,maxInsertion:10},solver=await loadHighs();
 const result=analyzeConditionalGroups({amplicons:input.references,primer:input.primer,control:input.control,sample:input.sample,ampliconConfirmed:true,settings:{...automaticSangerSettings(input),readEnd:1000,editStartRadius:3,maxIndel:10}},()=>{},solver);
 // Ground truth enters only after the production engine has returned.
 const truth=parseBenchmarkTruth(readFileSync(new URL('truth.json',base),'utf8')),native=certifyNativeBenchmark(truth,controlBytes,sampleBytes,result,solver),assessment=assessSangerBenchmark(truth,input.references,result,1,native);
 expect(result.diagnostics?.algorithmVersion).toBe('sanger-sequence-groups-v42-near-frequency-phase');expect(assessment.pass).toBe(true);expect(assessment.unexplained).toBe(0);
 const first=assessment.rows.filter(r=>r.reference==='Demo-A_6001-6790'&&r.target===0&&r.scope==='group'&&(r.truth??0)>0);
 expect(first.length).toBe(6);for(const r of first){expect(r.estimated).toBeGreaterThan(0);expect(Math.abs(r.delta!)).toBeLessThanOrEqual(2);}
},600000);
