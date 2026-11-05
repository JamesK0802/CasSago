import loadHighs from'../../public/vendor/highs-1.15.4/highs.mjs';
import{certifyNativeBenchmark}from'../../src/app/utils/sanger-native-signal-certificate';
import '@angular/compiler';
import {it,expect} from 'vitest';
import {randomBytes,createHash} from 'node:crypto';
import {readFileSync,mkdirSync,writeFileSync,readdirSync,existsSync} from 'node:fs';
import {profileFixture} from '../fixtures/sanger-profile-fixture';
import {randomSynthEvent,randomSynthPercent,generateSangerSynth,synthAlleleSequence} from '../../src/app/utils/sanger-synthetic-benchmark';
import {analyzeConditionalGroups} from '../../src/app/workers/core/sanger-conditional-groups';
import {automaticSangerSettings} from '../../src/app/workers/core/sanger-automatic';
import {benchmarkComparison,parseBenchmarkTruth} from '../../src/app/utils/sanger-benchmark';
import {assessSangerBenchmark} from '../../src/app/utils/sanger-benchmark-assessment';
const run=process.env.CASSAGO_PROSPECTIVE==='1'?it:it.skip;
run('prospectively validates fresh blind datasets against a frozen five-consecutive-pass gate',async()=>{
 const folder=process.env.CASSAGO_PROSPECTIVE_FOLDER!;if(!folder)throw Error('A new output folder is required.');
 const optimizer=await loadHighs();
 if(existsSync(folder+'/manifest.json'))throw Error('Use a fresh output folder; previous trials are immutable.');
 const vendorHash=createHash('sha256').update(readFileSync('public/vendor/highs-1.15.4/highs.mjs')).update(readFileSync('public/vendor/highs-1.15.4/highs.wasm')).digest('hex');
 const count=Number(process.env.CASSAGO_PROSPECTIVE_COUNT||5),files=readdirSync('src/app/workers/core').filter(f=>f.endsWith('.ts')).sort(),engineHash=()=>createHash('sha256').update(files.map(f=>f+'\n'+readFileSync('src/app/workers/core/'+f,'utf8')).join('\n')).digest('hex'),hash=engineHash(),seeds=Array.from({length:count},()=>randomBytes(4).readUInt32LE()),report:any={createdAt:new Date().toISOString(),optimizerVersion:'highs-1.15.4',optimizerHash:vendorHash,engineHash:hash,seeds,route:'Production conditional engine; an exception or failed QC is a failed case. General fallback and browser UI are not evaluated.',gate:{consecutive:5,editedGroups:30,maxDeltaPp:2,wtMaxDeltaPp:2,extraGroupsAtLeast1pp:0,requiresSignalQc:true,optionMembershipCredit:false,nativeControlAndEditByteIdentityCredit:true,exactWholeSignalEquivalenceCredit:true,offProfileBudgetPp:0},generatorHash:createHash('sha256').update(readFileSync('src/app/utils/sanger-synthetic-benchmark.ts')).update(readFileSync('src/app/utils/sanger-synth-profiles.ts')).digest('hex'),evaluationHash:createHash('sha256').update(readFileSync('src/app/utils/sanger-benchmark-assessment.ts')).update(readFileSync('src/app/utils/sanger-signal-equivalence-benchmark.ts')).update(readFileSync('src/app/utils/sanger-native-signal-certificate.ts')).digest('hex'),referencePanel:'Current three-reference two-target demo. New edit recipes, unequal percentages and molecular pairing per case; genomic references are held fixed.',cases:[],streak:0,maxStreak:0};
 mkdirSync(folder,{recursive:true});writeFileSync(folder+'/manifest.json',JSON.stringify(report,null,2));
 for(const [index,initialSeed] of seeds.entries()){
  expect(engineHash()).toBe(hash);let seed=initialSeed;const random=()=>((seed=(Math.imul(seed,1664525)+1013904223)>>>0)/4294967296),fixture=profileFixture(3),used:number[]=[];
  for(const ref of fixture.refs)for(const profile of ref.profiles!){profile.alleles=[];const sequences=new Set<string>();while(profile.alleles.length<5){const event=randomSynthEvent(ref,random,profile.target),sequence=synthAlleleSequence(ref,[event]);if(sequences.has(sequence))continue;sequences.add(sequence);const percent=randomSynthPercent(profile.alleles.map(a=>a.percent),random,used);used.push(percent);profile.alleles.push({id:profile.alleles.length+1,percent,events:[event]});}}
  const data=generateSangerSynth(fixture.result,fixture.candidate,fixture.refs,initialSeed^0xa17e963b),path=folder+'/case-'+(index+1);mkdirSync(path,{recursive:true});
  // Freeze generated files before calling inference. Truth is not an engine input.
  writeFileSync(path+'/control.ab1',new Uint8Array(data.control));writeFileSync(path+'/edit.ab1',new Uint8Array(data.edit));writeFileSync(path+'/template.json',JSON.stringify(data.template));writeFileSync(path+'/truth.json',JSON.stringify({recipe:data.recipe,truth:data.truth}));
  const stages:any[]=[],start=Date.now();let item:any={case:index+1,seed:initialSeed,pass:false};
  try{
   const result=analyzeConditionalGroups({amplicons:data.input.references,primer:data.input.primer!,sample:data.input.sample,control:data.input.control,ampliconConfirmed:true,settings:{...automaticSangerSettings(data.input),readEnd:1000,editStartRadius:3,maxIndel:10}},stage=>{stages.push({stage,ms:Date.now()-start});writeFileSync(path+'/progress.json',JSON.stringify(stages));},optimizer);
   // Evaluation begins only after blind inference has returned.
   const rows=benchmarkComparison(parseBenchmarkTruth({recipe:data.recipe,truth:data.truth}),data.input.references.map(r=>({...r,cuts:r.cuts!})),result),edited=rows.filter(r=>r.scope==='group'&&(r.truth??0)>0&&r.group!=='WT'),correct=(r:typeof rows[number])=>(r.estimated??0)>0&&Math.abs(r.delta!)<=2,wt=rows.filter(r=>r.group==='WT'),extras=rows.filter(r=>r.scope==='group'&&r.truth===0&&(r.estimated??0)>=1);
   item={...item,elapsedMs:Date.now()-start,algorithm:result.diagnostics?.algorithmVersion,signalQc:result.mixture?.passesChecks,correct:edited.filter(correct).length,total:edited.length,wtCorrect:wt.filter(correct).length,wtTotal:wt.length,extras:extras.length,byTarget:[0,1].map(t=>({target:t+1,correct:edited.filter(r=>r.target===t&&correct(r)).length,total:edited.filter(r=>r.target===t).length})),failures:edited.filter(r=>!correct(r)),pass:edited.length===30&&edited.every(correct)&&wt.length===6&&wt.every(correct)&&extras.length===0&&result.mixture?.passesChecks===true};
   const certificateStart=Date.now(),benchmarkTruth=parseBenchmarkTruth({recipe:data.recipe,truth:data.truth}),nativeCertificate=certifyNativeBenchmark(benchmarkTruth,data.control,data.edit,result,optimizer),assessment=assessSangerBenchmark(parseBenchmarkTruth({recipe:data.recipe,truth:data.truth}),data.input.references.map(r=>({...r,cuts:r.cuts!})),result,1,nativeCertificate);
   item={...item,pointRecovered:assessment.pointRecovered,signalEquivalent:assessment.signalEquivalent,unexplained:assessment.unexplained,unexplainedGroups:assessment.unexplainedGroups,wtAccepted:assessment.wtAccepted,falseExtraGroups:assessment.falseExtraGroups,referenceAllocationPassed:assessment.referenceAllocationPassed,certificateMs:Date.now()-certificateStart,elapsedMs:Date.now()-start,pass:assessment.pass&&assessment.total===30&&assessment.wtTotal===6};
   writeFileSync(path+'/native-equivalence-certificate.json',JSON.stringify(nativeCertificate));writeFileSync(path+'/signal-equivalence-certificates.json',JSON.stringify(assessment.certificates));writeFileSync(path+'/result.json',JSON.stringify(result));writeFileSync(path+'/comparison.json',JSON.stringify(rows,null,2));
  }catch(error){item={...item,elapsedMs:Date.now()-start,error:String(error)};}
  report.streak=item.pass?report.streak+1:0;report.maxStreak=Math.max(report.maxStreak,report.streak);report.cases.push(item);writeFileSync(path+'/summary.json',JSON.stringify(item,null,2));writeFileSync(folder+'/manifest.json',JSON.stringify(report,null,2));
 }
 expect(engineHash()).toBe(hash);report.finishedAt=new Date().toISOString();report.gatePassed=report.streak>=5;writeFileSync(folder+'/manifest.json',JSON.stringify(report,null,2));if(process.env.CASSAGO_ENFORCE_GATE==='1')expect(report.gatePassed).toBe(true);
},1800000);
