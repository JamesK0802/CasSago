import {SynthTargetProfile,mixTargetProfiles} from './sanger-synth-profiles';
import {referenceChecksum} from './reference-checksum';
import {SangerCandidate,SangerDesignResult} from '../models/sanger-design';
import {SangerAnalysisResult,SangerJointEvent,TRACE_BASES} from '../models/sanger-analysis';
import {selectedAmplicons,makeAnalysisTemplate,SangerAnalysisTemplate} from './sanger-design-export';
import {reverseComplement} from './biology.utils';
import {writeAbif} from './sanger-analysis-demo';
import {parseSangerTrace} from './sanger-trace';
import {AutomaticSangerInput} from '../workers/core/sanger-automatic';
import {sangerAnnotationViews,SangerAnnotationView} from './sanger-annotation';
export interface SynthEvent {kind:'insertion'|'deletion'|'substitution';target:number;offset:number;sequence:string;length:number;selectionStart?:number;before?:string}
export interface SynthAllele {id:number;percent:number;events:SynthEvent[]}
export interface SynthReference {name:string;sourceReference?:string;sequence:string;cuts:number[];targetNames:string[];share:number;primerLength:number;partnerLength:number;alleles:SynthAllele[];profiles?:SynthTargetProfile[]}
export interface SynthTruth {targetStates?:number[];reference:string;sourceReference:string;allele:string;share:number;percent:number;totalPercent:number;pattern:string;description:string;events:SynthEvent[];sequence:string;inputSequence:string;substitution:boolean}
export interface SynthDataset {control:ArrayBuffer;edit:ArrayBuffer;truth:SynthTruth[];template:SangerAnalysisTemplate;cycles:number;recipe:SynthReference[];mixing?:{method:'random-target-pairing';seed:number;populationPerReference:number};input:AutomaticSangerInput}
export interface SynthComparison {profileCertificate?:{kind:'exact-complete-emission';maxDeltaPp:number;offProfilePp:number};optionMatch?:{indices:number[];sequenceIndices:number[];withinTolerance:boolean;range:[number,number];profileOptions:number[];tolerancePp:number};target?:number;reference:string;group:string;truth:number|null;estimated:number|null;delta:number|null;recovery:string;placement:string;pattern?:string;scope?:'reference'|'group';sequenceKey?:string;sequenceScope?:'sequence'|'pattern'}
const dna=(s:string)=>s.replace(/\s/g,'').toUpperCase();
export function synthReferences(result:SangerDesignResult,c:SangerCandidate):SynthReference[]{const a=selectedAmplicons(result,c);return a.map((r,i)=>({name:r.name,sourceReference:r.sourceReference,sequence:c.direction==='reverse'?reverseComplement(r.sequence):r.sequence,cuts:r.cuts.map(cut=>c.direction==='reverse'?r.sequence.length-cut:cut),targetNames:r.targetNames||r.cuts.map((_,t)=>'T'+(t+1)),share:i===a.length-1?100-Math.round(10000/a.length)/100*(a.length-1):Math.round(10000/a.length)/100,primerLength:c.sequence.length,partnerLength:c.partner!.sequence.length,alleles:[]}));}
export function synthWT(ref:SynthReference){return 100-ref.alleles.reduce((s,a)=>s+Number(a.percent),0);}
/** Editable random defaults in sequencing orientation; never alter either primer site. */
export function randomSynthEvent(ref:SynthReference,random:()=>number=Math.random,onlyTarget?:number):SynthEvent {
 const pick=<T>(items:T[])=>items[Math.min(items.length-1,Math.floor(random()*items.length))];
 const end=ref.sequence.length-ref.partnerLength;
 const choices=(['insertion','deletion','substitution']as const).map(kind=>({kind,targets:ref.cuts.map((cut,target)=>({target,positions:Array.from({length:7},(_,i)=>cut+i-3).filter(start=>start>=ref.primerLength&&start<=end&&(kind==='insertion'||start<end&&/^[ACGT]$/.test(ref.sequence[start])))})).filter(t=>t.positions.length&&(onlyTarget===undefined||t.target===onlyTarget))})).filter(c=>c.targets.length);
 if(!choices.length)throw Error(ref.name+': no editable bases within ±3 bp of a cut.');
 const {kind,targets}=pick(choices),{target,positions}=pick(targets),start=pick(positions),offset=start-ref.cuts[target];
 let maxLength=kind==='insertion'?10:Math.min(10,end-start);
 if(kind==='substitution'){let valid=0;while(valid<maxLength&&/^[ACGT]$/.test(ref.sequence[start+valid]))valid++;maxLength=valid;}
 const length=1+Math.floor(random()*maxLength);
 const sequence=kind==='deletion'?'':Array.from({length},(_,i)=>pick('ACGT'.split('').filter(base=>kind!=='substitution'||base!==ref.sequence[start+i]))).join('');
 return {kind,target,offset,length,sequence,...(kind==='substitution'?{selectionStart:start,before:ref.sequence.slice(start,start+length)}:{})};
}

/** Random editable default, drawn without replacing existing or manually edited
 * percentages. Integer hundredths keep WT subtraction and finite mixing exact. */
export function randomSynthPercent(existing:number[],random:()=>number=Math.random,avoid:number[]=[]):number {
 const remaining=Math.floor((100-existing.reduce((sum,value)=>sum+Number(value),0))*100+1e-6);
 if(remaining<1)throw Error('No WT percentage remains. Reduce an existing group first.');
 const upper=Math.min(1700,remaining),lower=Math.min(300,upper),used=new Set([...existing,...avoid].map(value=>Math.round(value*100))),choices:number[]=[];
 for(let value=lower;value<=upper;value++)if(!used.has(value))choices.push(value);
 if(!choices.length)for(let value=1;value<=remaining;value++)if(!used.has(value))choices.push(value);
 if(!choices.length)throw Error('No distinct positive percentage remains. Adjust existing groups first.');
 return choices[Math.min(choices.length-1,Math.floor(random()*choices.length))]/100;
}

function resolved(ref:SynthReference,event:SynthEvent){
 const sequence=dna(event.sequence||''),start=event.kind==='substitution'?event.selectionStart:ref.cuts[event.target]+Number(event.offset);
 if(!Number.isInteger(start)||start===undefined)throw Error(ref.name+': select a valid edit position.');
 const removed=event.kind==='insertion'?0:event.kind==='deletion'?Number(event.length):(event.before||'').length;
 if(event.kind==='insertion'&&(!/^[ACGT]+$/.test(sequence)||sequence.length>30))throw Error('Insertion requires 1–30 A/C/G/T bases.');
 if(event.kind==='deletion'&&(!Number.isInteger(removed)||removed<1||removed>500))throw Error('Deletion length must be 1–500 bp.');
 if(event.kind==='substitution'&&(!event.before||sequence.length!==event.before.length||!/^[ACGT]+$/.test(sequence)||sequence===event.before||ref.sequence.slice(start,start+removed)!==event.before))throw Error('Select the original bases and enter a different A/C/G/T sequence of the same length.');
 if(start<ref.primerLength||start+removed>ref.sequence.length-ref.partnerLength)throw Error('Edits must stay between the two primer binding sites.');
 return {start,end:start+removed,inserted:event.kind==='deletion'?'':sequence,event};
}
export function synthAlleleSequence(ref:SynthReference,events:SynthEvent[]):string {
 if(!events.length)throw Error('Add at least one edit to each variant allele.');
 const edits=events.map(e=>resolved(ref,e)).sort((a,b)=>a.start-b.start);
 for(let i=1;i<edits.length;i++)if(edits[i].start<edits[i-1].end||edits[i].start===edits[i-1].start)throw Error('Edits on one allele must not overlap or share an insertion boundary.');
 let out=ref.sequence;for(const e of edits.reverse())out=out.slice(0,e.start)+e.inserted+out.slice(e.end);return out;
}
export function synthPattern(ref:SynthReference,events:SynthEvent[]){const sizes=ref.cuts.map(()=>0);for(const e of events){if(e.kind!=='substitution')sizes[e.target]+=(e.kind==='insertion'?dna(e.sequence).length:-Number(e.length));}return sizes.join('/');}
export function synthDescription(events:SynthEvent[]){return events.length?events.map(e=>e.kind==='substitution'?`${e.before}→${dna(e.sequence)} at amplicon ${Number(e.selectionStart)+1}`:`T${e.target+1} ${e.kind} ${e.offset>=0?'+':''}${e.offset}: ${e.kind==='insertion'?dna(e.sequence):e.length+' bp'}`).join('; '):'WT';}
export function writeSynthTrace(parts:{sequence:string;weight:number}[],cycles:number,name:string){const peaks=Array.from({length:cycles},(_,i)=>14+i*10),channels=TRACE_BASES.map(()=>new Uint16Array(cycles*10+30).fill(2)),calls:string[]=[];
 for(let i=0;i<cycles;i++){const values=TRACE_BASES.map(b=>parts.reduce((s,p)=>s+(p.sequence[i]===b?p.weight:0),0));calls.push(TRACE_BASES[values.indexOf(Math.max(...values))]);for(let b=0;b<4;b++)for(let d=-5;d<=5;d++)channels[b][peaks[i]+d]+=Math.round(1200*values[b]*Math.exp(-d*d/4));}
 return writeAbif(calls.join(''),peaks,channels,'SYNTHETIC_'+name);
}
export function generateSangerSynth(result:SangerDesignResult,c:SangerCandidate,refs:SynthReference[],mixSeed?:number):SynthDataset {
 const seed=(mixSeed??Math.floor(Math.random()*4294967296))>>>0;let randomState=seed;const random=()=>{randomState=(Math.imul(randomState,1664525)+1013904223)>>>0;return randomState/4294967296;};
 if(refs.some(r=>!Number.isFinite(r.share)||r.share<0)||Math.abs(refs.reduce((s,r)=>s+r.share,0)-100)>.001)throw Error('Reference shares must total 100%.');
 const template=makeAnalysisTemplate(result,c),products=selectedAmplicons(result,c);
 if(refs.length!==products.length||refs.some((r,i)=>r.name!==products[i].name))throw Error('The recipe no longer matches the selected primer pair.');
 const truth:SynthTruth[]=[],controlParts:{sequence:string;weight:number}[]=[],editParts:{sequence:string;weight:number}[]=[];let maxDeletion=10,maxInsertion=10;
 refs.forEach((ref,i)=>{
  const wt=ref.profiles?100:synthWT(ref);if(!Number.isFinite(wt)||wt<-.000001||ref.alleles.some(a=>!Number.isFinite(a.percent)||a.percent<0||a.percent>100))throw Error(ref.name+': variant percentages must total no more than 100%.');
  if(/[^ACGT]/.test(ref.sequence))throw Error('Resolve N or other ambiguous bases before generating synthetic traces.');
  const expected=c.direction==='reverse'?reverseComplement(products[i].sequence):products[i].sequence;
  if(ref.sequence!==expected)throw Error('Reference sequence cannot be edited directly.');
  const add=(allele:string,percent:number,events:SynthEvent[],targetStates?:number[])=>{if(percent<=0||ref.share===0)return;const sequence=events.length?synthAlleleSequence(ref,events):ref.sequence;for(const e of events){if(e.kind==='deletion')maxDeletion=Math.max(maxDeletion,e.length);if(e.kind==='insertion')maxInsertion=Math.max(maxInsertion,dna(e.sequence).length);}const inputSequence=c.direction==='reverse'?reverseComplement(sequence):sequence;truth.push({targetStates,reference:ref.name,sourceReference:products[i].sourceReference,allele,share:ref.share,percent,totalPercent:ref.share*percent/100,pattern:synthPattern(ref,events),description:synthDescription(events),events:structuredClone(events),sequence,inputSequence,substitution:events.some(e=>e.kind==='substitution')});editParts.push({sequence:sequence.slice(ref.primerLength),weight:ref.share/100*percent/100});};
  controlParts.push({sequence:ref.sequence.slice(ref.primerLength),weight:ref.share/100});if(ref.profiles){for(const group of mixTargetProfiles(ref,random))add(group.allele,group.percent,group.events,group.targetStates);}else{add('WT',Math.max(0,wt),[]);ref.alleles.forEach((a,j)=>{if(a.events.length!==1)throw Error('Each group must contain exactly one edit.');add('Group '+(j+1),a.percent,a.events);});}
 });
 const active=[...controlParts,...editParts].filter(p=>p.weight>0),cycles=Math.min(result.settings.readEnd,...active.map(p=>p.sequence.length))-20;
 if(cycles<150)throw Error('The edited products leave fewer than 150 complete sequencing cycles.');
 if(c.sites.some(h=>h.cutOffsets.some(cut=>cut+Math.ceil(result.settings.window/2)>cycles)))throw Error('An edited product ends before the full target window. Use a longer product or a shorter deletion.');
 const control=writeSynthTrace(controlParts,cycles,'control'),edit=writeSynthTrace(editParts,cycles,'edit');template.settings={...template.settings,maxDeletion,maxInsertion,maxIndel:Math.min(30,Math.max(maxDeletion,maxInsertion)),analysisPolicy:'automatic'};
 const input:AutomaticSangerInput={references:template.amplicons.map(a=>({name:a.name,sequence:a.sequence,cut:a.cuts[0],cuts:a.cuts})),primer:template.primer,control:{...parseSangerTrace(control,'SYNTHETIC_control.ab1'),synthetic:true},sample:{...parseSangerTrace(edit,'SYNTHETIC_edit.ab1'),synthetic:true},maxDeletion,maxInsertion,synthetic:true};
 return {control,edit,truth,template,cycles,recipe:structuredClone(refs),mixing:refs.some(r=>r.profiles)?{method:'random-target-pairing',seed,populationPerReference:10000}:undefined,input};
}
export function synthTruthCsv(data:SynthDataset){const rows=[['Reference','Source reference','Group','Reference %','Within-reference %','Total sample %','Target size pattern','Edits'],...data.truth.map(t=>[t.reference,t.sourceReference,t.allele,t.share,t.percent,t.totalPercent,t.pattern,t.description])];return rows.map(r=>r.map(v=>'"'+String(v).replace(/"/g,'""')+'"').join(',')).join('\r\n');}
export function compareSangerSynth(data:{truth:SynthTruth[];template:{amplicons:{name:string;sequence:string;cuts:number[]}[]}},result:SangerAnalysisResult,displayViews?:SangerAnnotationView[]):SynthComparison[]{
 const out:SynthComparison[]=[],views=displayViews||sangerAnnotationViews(result.annotations||[],0,60,true);
 for(const ref of data.template.amplicons){
  const truth=data.truth.filter(t=>t.reference===ref.name),view=views.find(a=>a.name===ref.name),count=ref.cuts.length,zero=ref.cuts.map(()=>0).join('/'),withSub=truth.some(t=>t.substitution),sequenceMode=truth.every(t=>!!t.inputSequence),groups=new Map<string,typeof truth>();
  const sequenceGroups=result.diagnostics?.algorithmVersion?.startsWith('sanger-sequence-groups-v');
  const expectedShare=truth[0]?.share??0,estimatedShare=view?view.share*100:result.pooled?.amplicons.find(a=>a.name===ref.name)?.share!==undefined?result.pooled.amplicons.find(a=>a.name===ref.name)!.share*100:data.template.amplicons.length===1?100:0;
  out.push({reference:ref.name,group:'Reference share',scope:'reference',truth:expectedShare,estimated:estimatedShare,delta:estimatedShare-expectedShare,recovery:'Reference allocation',placement:'—'});
  for(const row of truth){const key=row.substitution&&!sequenceGroups?'sub:'+row.description:sequenceMode?row.pattern+':'+referenceChecksum(row.inputSequence):row.pattern;groups.set(key,[...(groups.get(key)||[]),row]);}
  type Inferred={fraction:number;pattern:string;sequenceKey?:string;unphased:boolean;unknown:boolean};const inferred=new Map<string,Inferred>();
  for(const row of view?.rows||[]){const pattern=row.patternLabel?row.patternLabel.split(' · ').map(p=>Number(p.split(': ')[1].replace(' bp',''))).join('/'):count===1?String(row.size):null;if(pattern===null)continue;const unknown=!row.sequenceKey,key=sequenceMode?pattern+':'+(row.sequenceKey||'unresolved'):pattern,old=inferred.get(key);if(old){old.fraction+=row.fraction*100;old.unphased||=row.sequenceState==='unphased';old.unknown||=unknown;}else inferred.set(key,{fraction:row.fraction*100,pattern,sequenceKey:row.sequenceKey,unphased:row.sequenceState==='unphased',unknown});}
  for(const [key,rows] of groups){
   const pattern=rows[0].pattern,expected=rows.reduce((sum,t)=>sum+t.percent,0),unresolved=key.startsWith('sub:')||!sequenceGroups&&withSub&&pattern===zero,actual=unresolved?undefined:inferred.get(key),unknownMass=sequenceMode&&[...inferred.values()].some(r=>r.pattern===pattern&&r.unknown&&r.fraction>0),estimated=unresolved||!actual&&unknownMass?null:actual?.fraction??0;
   const placement=unresolved?'Substitution attribution not resolved by current engine':!actual?unknownMass?'Sequence not attributable':'Not recovered':!sequenceMode?'Truth sequence not provided':actual.unphased?'Sequence candidate matched; linkage unresolved':'Sequence matched';
   out.push({reference:ref.name,group:rows.map(r=>r.description).join(' | '),scope:'group',pattern,sequenceKey:sequenceMode&&!unresolved?referenceChecksum(rows[0].inputSequence):undefined,sequenceScope:sequenceMode?'sequence':'pattern',truth:expected,estimated,delta:estimated===null?null:estimated-expected,recovery:estimated===null?'Not attributable':actual?.unphased?'Candidate group':estimated>0?'Recovered group':'Missed group',placement});
   if(!unresolved)inferred.delete(key);
  }
  for(const row of inferred.values()){const unresolved=row.unknown||!sequenceGroups&&withSub&&row.pattern===zero;out.push({reference:ref.name,group:'Estimated '+row.pattern+' bp',scope:'group',pattern:row.pattern,sequenceKey:row.sequenceKey,sequenceScope:sequenceMode?'sequence':'pattern',truth:unresolved?null:0,estimated:row.fraction,delta:unresolved?null:row.fraction,recovery:unresolved?'Unresolved sequences':'Extra group',placement:row.unknown?'Sequence not attributable':row.unphased?'Unphased sequence candidate':'Estimated sequence available'});}
 }
 return out;
}

export function synthPickerContext(ref:SynthReference,target:number,guides:string[]=[]) {
 const cut=ref.cuts[target],start=Math.max(ref.primerLength,cut-60),end=Math.min(ref.sequence.length-ref.partnerLength,cut+60);
 const hits:{start:number;end:number;target:number}[]=[];
 guides.forEach((g,t)=>{g=dna(g);if(!g)return;for(const query of new Set([g,reverseComplement(g)])){let at=ref.sequence.indexOf(query,start);while(at>=0&&at<end){if(at<=ref.cuts[t]&&at+query.length>=ref.cuts[t])hits.push({start:at,end:at+query.length,target:t});at=ref.sequence.indexOf(query,at+1);}}});
 const cells=Array.from({length:Math.max(0,end-start)},(_,i)=>{const index=start+i;return {index,base:ref.sequence[index],relative:index-cut,cut:ref.cuts.includes(index),cutName:ref.targetNames[ref.cuts.indexOf(index)]||'',guide:hits.some(g=>index>=g.start&&index<g.end)};});
 const rows=Array.from({length:Math.ceil(cells.length/60)},(_,i)=>({relative:cells[i*60].relative,cells:cells.slice(i*60,(i+1)*60)}));
 return {start,end,cut,rows,hasGuide:hits.length>0};
}
