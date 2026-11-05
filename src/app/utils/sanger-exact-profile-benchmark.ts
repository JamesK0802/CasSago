import{projectTargetTruth,targetEvents}from'./sanger-target-profile';
import type{BenchmarkTruth}from'./sanger-benchmark';
import type{SangerAnalysisResult}from'../models/sanger-analysis';
import type{SynthComparison}from'./sanger-synthetic-benchmark';
import{sangerEditedSequence}from'./sanger-sequence-groups';
import{reverseComplement}from'./biology.utils';
import{exactBoundedMixture,exactFromNumber,exactNumber,exactSum,exactWeightedEquality,exactSubtract,exactScale,exactProduct}from'../workers/core/exact-mixture-fiber';
function localSequence(reference:string,sequence:string,cut:number){
 const size=sequence.length-reference.length;if(!size){let start=0,end=reference.length;while(start<end&&reference[start]===sequence[start])start++;if(start===end)return true;while(end>start&&reference[end-1]===sequence[end-1])end--;return start>=cut-3&&start<=cut+3&&end-start<=10;}
 for(let start=Math.max(0,cut-3);start<=cut+3;start++){const end=start+Math.max(0,-size),inserted=size>0?sequence.slice(start,start+size):'';if(reference.slice(0,start)+inserted+reference.slice(end)===sequence)return true;}return false;
}
/** Post-inference membership of a whole profile in the fixed signal fiber.
 * This never generates inferred alleles, selects a dropdown option or changes
 * point estimates. Finite displayed representatives need not include every
 * compatible point of a continuous mixture polytope. */
export function compareExactProfiles(data:BenchmarkTruth,rows:SynthComparison[],result:SangerAnalysisResult,tolerancePp=2):SynthComparison[]{
 const d=result.diagnostics;if(!d||!/^sanger-sequence-groups-v(?:5-offset|6-conditional)/.test(d?.algorithmVersion||'')||d.contextEmission?.used||!d.contextEmission?.theta.every(x=>x===0)||d.sampleCalibration?.lag!==0)return rows;
 const certificates=new Map<string,NonNullable<SynthComparison['profileCertificate']>>();
 for(const a of result.annotations||[]){
  if((a.cuts?.length||1)===2){const certificate=twoTargetCertificate(data,a,result,tolerancePp);if(certificate)certificates.set(a.name,certificate);continue;}
  if((a.cuts?.length||1)!==1||!a.share)continue;
  const truth=data.truth.filter(t=>t.reference===a.name);if(!truth.length||truth.length>16||truth.some(t=>!t.inputSequence||/[^ACGT]/.test(t.inputSequence)))continue;
  const oriented=(s:string)=>a.direction==='reverse'?reverseComplement(s):s,set=result.exactAlternatives?.find(s=>s.reference===a.name),base=set?.options[0],components=base?.components||a.components;if(components.length>24||components.some(c=>c.unresolved))continue;
  const weights=base?.exactFractions||components.map(c=>exactFromNumber(c.fraction*a.share)),totalMass=exactSum(weights),cap=exactScale(totalMass,1,100),massAt=(percent:number)=>exactScale(exactProduct(exactFromNumber(percent),totalMass),1,100),share=exactNumber(totalMass),gt=truth.map(t=>oriented(t.inputSequence)),fitted=components.map(c=>oriented(sangerEditedSequence(a.reference,c.events||[]))),extras=[...new Set(fitted.filter(s=>!gt.includes(s)))],sequences=[...gt,...extras];if(new Set(gt).size!==gt.length||gt.some(s=>!localSequence(oriented(a.reference),s,a.direction==='reverse'?a.reference.length-a.cut:a.cut)))continue;
  const length=Math.max(...sequences.map(s=>s.length));if(length>4000)continue;
  const column=(s:string)=>Array.from({length},(_,p)=>Array.from('ACGT',b=>Number(s[p]===b))).flat(),columns=sequences.map(column),source=fitted.map((s,i)=>({column:column(s),weight:weights[i]})),bounds=[...truth.map(t=>({lower:massAt(Math.max(0,t.percent-tolerancePp)),upper:massAt(Math.min(100,t.percent+tolerancePp))})),...extras.map(()=>({lower:exactFromNumber(0),upper:exactFromNumber(.01*share)}))],witness=exactBoundedMixture(columns,source,bounds);
  if(!witness)continue;const offProfilePp=witness.slice(gt.length).reduce((sum,w)=>sum+exactNumber(w),0)/share*100;if(offProfilePp>1||!exactWeightedEquality(source,columns.map((column,i)=>({column,weight:witness[i]}))))continue;
  certificates.set(a.name,{kind:'exact-complete-emission',maxDeltaPp:Math.max(...truth.map((t,i)=>Math.abs(exactNumber(witness[i])/share*100-t.percent))),offProfilePp});
 }
 return rows.map(row=>certificates.has(row.reference)?{...row,profileCertificate:certificates.get(row.reference)}:row);
}

/** Auxiliary bounded variables encode target marginals. Their defining rows
 * conserve each complete target partition, so the augmented mass constraint is
 * redundant with physical mass; one witness couples both targets together. */
function twoTargetCertificate(data:BenchmarkTruth,a:NonNullable<SangerAnalysisResult['annotations']>[number],result:SangerAnalysisResult,tolerancePp:number):NonNullable<SynthComparison['profileCertificate']>|null{
 if(a.direction!=='forward')return null;
 const truth=data.truth.filter(t=>t.reference===a.name);if(!truth.length||truth.length>64||truth.some(t=>!t.inputSequence||/[^ACGT]/.test(t.inputSequence)))return null;
 const convert=(t:typeof truth[number])=>t.events.map(e=>{const start=Number(e.kind==='substitution'?e.selectionStart:a.cuts![e.target]+e.offset),end=start+(e.kind==='insertion'?0:e.length),insertedSequence=e.kind==='deletion'?'':e.sequence;return{target:e.target+1,start,end,inserted:insertedSequence.length,insertedSequence};}),truthComponents=truth.map(t=>({events:convert(t)}));
 if(truth.some((t,i)=>sangerEditedSequence(a.reference,truthComponents[i].events)!==t.inputSequence))return null;
 const projected=[0,1].map(t=>projectTargetTruth(data,a,t));if(projected.some((p,t)=>p.truth.some(row=>!localSequence(a.reference,row.inputSequence,a.cuts![t]))))return null;
 const base=result.exactAlternatives?.find(s=>s.reference===a.name)?.options[0],source=base?.components||a.components;if(source.some(c=>c.unresolved))return null;
 const weights=base?.exactFractions||source.map(c=>exactFromNumber(c.fraction*a.share)),totalMass=exactSum(weights),cap=exactScale(totalMass,1,100),massAt=(percent:number)=>exactScale(exactProduct(exactFromNumber(percent),totalMass),1,100),share=exactNumber(totalMass),gt=truth.map(t=>t.inputSequence),fitted=source.map(c=>sangerEditedSequence(a.reference,c.events||[])),extras=[...new Set(fitted.filter(s=>!gt.includes(s)))],sequences=[...gt,...extras],components=[...truthComponents,...extras.map(sequence=>source.find(c=>sangerEditedSequence(a.reference,c.events||[])===sequence)!)];if(new Set(gt).size!==gt.length)return null;
 const buckets:[number,string,number][]=[];for(const [target,p]of projected.entries()){for(const row of p.truth)buckets.push([target,row.inputSequence,row.percent]);buckets.push([target,'unmatched',0]);}buckets.push([-1,'unmatched',0],[-1,'full-truth',100]);if(sequences.length+buckets.length>128)return null;
 const known=[0,1].map(t=>new Set(buckets.filter(b=>b[0]===t&&b[1]!=='unmatched').map(b=>b[1]))),matches=(events:NonNullable<typeof a.components[number]['events']>,bucket:typeof buckets[number])=>{if(bucket[0]===-1){const inTruth=gt.includes(sangerEditedSequence(a.reference,events));return bucket[1]==='full-truth'?inTruth:!inTruth;}const key=sangerEditedSequence(a.reference,targetEvents(events,bucket[0]));return bucket[1]==='unmatched'?!known[bucket[0]].has(key):key===bucket[1];},length=Math.max(...sequences.map(s=>s.length));if(length>4000)return null;
 const column=(s:string)=>Array.from({length},(_,p)=>Array.from('ACGT',b=>Number(s[p]===b))).flat(),physical=sequences.map(column),columns=physical.map((c,i)=>[...c,...buckets.map(b=>Number(matches(components[i].events||[],b)))]),bounds=sequences.map((_,i)=>({lower:exactFromNumber(0),upper:i<gt.length?totalMass:cap})),original=fitted.map((s,i)=>({column:column(s),weight:weights[i]})),sourceColumns=fitted.map((s,i)=>({column:[...column(s),...buckets.map(b=>Number(matches(source[i].events||[],b)))],weight:weights[i]}));
 for(const [j,bucket]of buckets.entries()){const aux=[...Array(4*length).fill(0),...buckets.map((_,k)=>j===k?-1:0)],mass=exactSum(source.flatMap((c,i)=>matches(c.events||[],bucket)?[weights[i]]:[]));columns.push(aux);bounds.push({lower:massAt(Math.max(0,bucket[2]-tolerancePp)),upper:bucket[1]==='unmatched'?cap:massAt(Math.min(100,bucket[2]+tolerancePp))});sourceColumns.push({column:aux,weight:mass});}
 const witness=exactBoundedMixture(columns,sourceColumns,bounds);if(!witness)return null;const offMass=exactSum(witness.slice(gt.length,sequences.length)),offProfilePp=Math.min(1,exactNumber(offMass)/share*100);if(BigInt(exactSubtract(offMass,cap).n)>0n||!exactWeightedEquality(original,physical.map((column,i)=>({column,weight:witness[i]}))))return null;
 return{kind:'exact-complete-emission',offProfilePp,maxDeltaPp:Math.max(...buckets.flatMap((b,i)=>b[1]==='unmatched'||b[0]===-1?[]:[Math.abs(exactNumber(witness[sequences.length+i])/share*100-b[2])]))};
}
