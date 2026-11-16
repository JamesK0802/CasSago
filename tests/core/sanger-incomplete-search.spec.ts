import '@angular/compiler';
import{it,expect,beforeAll}from'vitest';
import loadHighs from'../../public/vendor/highs-1.15.4/highs.mjs';
import{profileFixture}from'../fixtures/sanger-profile-fixture';
import{generateSangerSynth}from'../../src/app/utils/sanger-synthetic-benchmark';
import{analyzeAutomaticSanger}from'../../src/app/workers/core/sanger-automatic';
function input(){const f=profileFixture(1);f.refs[0].profiles![0].alleles=[{id:1,percent:50,events:[{kind:'deletion',target:0,offset:0,length:1,sequence:''}]}];f.refs[0].profiles![1].alleles=[];const data=generateSangerSynth(f.result,f.candidate,f.refs,19072);return{...data.input,maxDeletion:1,maxInsertion:1,groupRecovery:true};}
let optimizer:any;beforeAll(async()=>{optimizer=await loadHighs();});
it('withholds group certification after a truncated search even when later fitting returns a spectrum',()=>{
 let interrupted=false;const stages:string[]=[];
 const result=analyzeAutomaticSanger(input(),stage=>{stages.push(stage);if(!interrupted&&stage.startsWith('Proposing target 1 groups')){interrupted=true;throw Object.assign(Error('Injected proposal work budget exhaustion'),{name:'SangerWorkBudgetExceeded'});}},optimizer);
 expect(interrupted).toBe(true);expect(stages.some(s=>s.startsWith('Certifying selectable sequence alternatives'))).toBe(false);expect(result.exactAlternatives).toEqual([]);expect(result.diagnostics?.['initialization']).toMatchObject({status:'unavailable',candidates:0});expect(result.diagnostics?.['measuredValidation']).toMatchObject({profileSearchComplete:false});expect(result.mixture?.passesChecks).toBe(false);expect(result.annotations!.every(a=>!a.estimable)).toBe(true);expect(result.warnings?.some(s=>s.includes('reconstruction did not complete'))).toBe(true);
},120000);
it('keeps completed common-engine reconstruction eligible for validation',()=>{const result=analyzeAutomaticSanger(input(),()=>{},optimizer);expect(result.diagnostics?.['measuredValidation']).toMatchObject({profileSearchComplete:true});expect(result.mixture?.passesChecks).toBe(true);},120000);
