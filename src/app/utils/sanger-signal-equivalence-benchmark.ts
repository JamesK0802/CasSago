import type {BenchmarkTruth} from './sanger-benchmark';
import type {SangerAnalysisResult,SangerAnnotation} from '../models/sanger-analysis';
import {sangerEditedSequence} from './sanger-sequence-groups';
import {projectTargetTruth,targetEvents} from './sanger-target-profile';
import {reverseComplement} from './biology.utils';
import {exactBoundedMixture,exactFromNumber,exactNumber,exactSum,exactProduct,exactScale,exactSubtract,exactWeightedEquality} from '../workers/core/exact-mixture-fiber';
import type {ExactFraction} from '../workers/core/exact-mixture-fiber';

export interface SignalEquivalenceCertificate {
 reference:string;
 kind:'exact-complete-emission';
 maxDeltaPp:number;
 offProfilePp:0;
 /** One coherent whole-molecule witness, not independent matches per target. */
 truthWeights:ExactFraction[];
 sourceWeights:ExactFraction[];
 comparedBases:number;
}
/** Post-inference scoring only. A certificate must reproduce EVERY modeled
 * cycle/channel exactly using ONLY truth molecules. No extra allele budget,
 * approximate equality, added inference candidates or selected UI option.
 * Null means uncertified, not a claim that a different solution cannot exist. */
export function certifySignalEquivalence(data:BenchmarkTruth,result:SangerAnalysisResult,tolerancePp=2):SignalEquivalenceCertificate[]{
 const d=result.diagnostics;
 if(!d||d.contextEmission?.used||!d.contextEmission?.theta?.every(x=>x===0)||d.sampleCalibration?.lag!==0)return[];
 return(result.annotations||[]).flatMap(a=>{const c=certifyReference(data,result,a,tolerancePp);return c?[c]:[];});
}
function certifyReference(data:BenchmarkTruth,result:SangerAnalysisResult,a:SangerAnnotation,tolerancePp:number):SignalEquivalenceCertificate|null {
 const truth=data.truth.filter(t=>t.reference===a.name),source=a.components;
 if(!truth.length||truth.length>64||!source.length||source.some(c=>c.unresolved)||truth.some(t=>!t.inputSequence||/[^ACGT]/.test(t.inputSequence)))return null;
 const cuts=a.cuts||[a.cut],convert=(events:typeof truth[number]['events'])=>events.map(e=>{const start=Number(e.kind==='substitution'?e.selectionStart:cuts[e.target]+e.offset),end=start+(e.kind==='insertion'?0:e.kind==='deletion'?e.length:(e.before||'').length),word=e.kind==='deletion'?'':e.sequence;return{target:e.target+1,start,end,inserted:word.length,insertedSequence:word};}),truthEvents=truth.map(t=>convert(t.events));
 if(truth.some((t,i)=>sangerEditedSequence(a.reference,truthEvents[i])!==t.inputSequence))return null;
 const oriented=(s:string)=>a.direction==='reverse'?reverseComplement(s):s,sourceSequences=source.map(c=>oriented(sangerEditedSequence(a.reference,c.events||[]))),truthSequences=truth.map(t=>oriented(t.inputSequence));
 // Prefer original fitted masses to a round trip through within-ref fractions.
 const base=result.exactAlternatives?.find(s=>s.reference===a.name)?.options[0],joint=result.joint?.groups.filter(g=>g.reference===a.name),weights=base?.exactFractions.length===source.length?base.exactFractions:joint?.length===source.length&&joint.every((g,i)=>oriented(sangerEditedSequence(a.reference,g.events))===sourceSequences[i])?joint.map(g=>exactFromNumber(g.fraction)):source.map(c=>exactFromNumber(c.fraction*a.share)),total=exactSum(weights),massAt=(percent:number)=>exactScale(exactProduct(exactFromNumber(percent),total),1,100),length=Math.max(...truthSequences.map(s=>s.length),...sourceSequences.map(s=>s.length))+2;
 if(length>4000||new Set(truthSequences).size!==truth.length||sourceSequences.some(s=>/[^ACGT]/.test(s)))return null;
 const column=(s:string)=>Array.from({length},(_,p)=>Array.from('ACGT',b=>Number(s[p]===b))).flat(),physical=truthSequences.map(column),original=sourceSequences.map((s,i)=>({column:column(s),weight:weights[i]}));
 type Bucket={target:number;sequence:string;percent:number};const buckets:Bucket[]=[],known=new Map<number,Set<string>>();
 for(let target=0;target<cuts.length;target++){
  const projected=projectTargetTruth(data,a,target);known.set(target,new Set(projected.truth.map(t=>t.inputSequence)));
  for(const t of projected.truth)buckets.push({target,sequence:t.inputSequence,percent:t.percent});
  // These partition rows make the auxiliary total-mass equation redundant.
  buckets.push({target,sequence:'unmatched',percent:0});
 }
 const matches=(events:NonNullable<typeof source[number]['events']>,b:Bucket)=>{const local=sangerEditedSequence(a.reference,targetEvents(events,b.target));return b.sequence==='unmatched'?!known.get(b.target)!.has(local):local===b.sequence;},columns=physical.map((c,i)=>[...c,...buckets.map(b=>Number(matches(truthEvents[i],b)))]),bounds=truth.map(()=>({lower:exactFromNumber(0),upper:total})),rhs=original.map((row,i)=>({column:[...row.column,...buckets.map(b=>Number(matches(source[i].events||[],b)))],weight:row.weight}));
 for(const [j,b]of buckets.entries()){
  const auxiliary=[...Array(4*length).fill(0),...buckets.map((_,i)=>i===j?-1:0)],mass=exactSum(source.flatMap((c,i)=>matches(c.events||[],b)?[weights[i]]:[]));columns.push(auxiliary);rhs.push({column:auxiliary,weight:mass});bounds.push({lower:b.sequence==='unmatched'?exactFromNumber(0):massAt(Math.max(0,b.percent-tolerancePp)),upper:b.sequence==='unmatched'?exactFromNumber(0):massAt(Math.min(100,b.percent+tolerancePp))});
 }
 // Identical coefficient rows are the same exact equation, not extra evidence.
 const rowIds:number[]=[],seenRows=new Set<string>();for(let p=0;p<columns[0].length;p++){const key=columns.map(c=>c[p]).join(',')+'|'+rhs.map(r=>r.column[p]).join(',');if(!seenRows.has(key)){seenRows.add(key);rowIds.push(p);}}
 const witness=exactBoundedMixture(columns.map(c=>rowIds.map(p=>c[p])),rhs.map(r=>({...r,column:rowIds.map(p=>r.column[p])})),bounds);
 if(!witness)return null;const physicalWeights=witness.slice(0,truth.length);
 if(exactSubtract(exactSum(physicalWeights),total).n!=='0'||!exactWeightedEquality(original,physical.map((c,i)=>({column:c,weight:physicalWeights[i]}))))return null;
 const maxDeltaPp=Math.max(...buckets.map((b,j)=>b.sequence==='unmatched'?0:Math.abs(exactNumber(witness[truth.length+j])/exactNumber(total)*100-b.percent)));
 return{reference:a.name,kind:'exact-complete-emission',maxDeltaPp,offProfilePp:0,truthWeights:physicalWeights,sourceWeights:weights,comparedBases:length};
}
