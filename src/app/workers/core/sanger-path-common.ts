import {normalizeMoleculeEdits} from './sanger-edit-normalize';
import {validateSangerProductBounds} from '../../utils/sanger-product-bounds';
import {repairSignalTraceCycles} from './sanger-cycle';
import {PooledAmplicon,PooledPreflight,SangerAnalysisSettings,SangerAnalysisResult,SangerJointEvent,SangerTrace,TRACE_BASES} from '../../models/sanger-analysis';
import {exactBindingSites} from './sanger-primer-designer';
import {reverseComplement} from '../../utils/biology.utils';
import {registerTrace} from './sanger-pooled';
import {contrastPeaks,integratedPeaks,calibrateSignal,correctGain,gramCondition,rawSignalScore} from './sanger-signal';
import {emissionVector,fitContextEmission,fitNormalizedEmissionMixture as fitEmissionMixtureOriginal} from './sanger-emission';
import {plotTrace} from './sanger-analysis';
const DNA='ACGT';
const fitEmissionMixture=(columns:number[][],observed:number[])=>fitEmissionMixtureOriginal(columns,observed);
export interface PathEdit {start:number;end:number;inserted:string}
export interface DiscoveredPath {sequence:string;events:PathEdit[];score:number}
/** Affine match/substitution/insertion/deletion alignment of a signal profile.
 * Only training-position signal enters profile; unobserved inserted bases are N.
 * This is a bounded Viterbi-style optimization, not a posterior-probability HMM. */
export function discoverSignalPath(reference:string,profile:number[][],anchor:number,end:number,cuts:number[],options:{maxDeletion:number;maxInsertion:number;alpha:number;local:boolean;training:boolean[];startRadius?:number}):DiscoveredPath|null{
 const n=end-anchor,m=Math.min(reference.length-anchor,end-anchor+options.maxDeletion),cols=m+1,cells=(n+1)*cols,NEG=-1e20;
 if(n<40||m<40||cells>2_000_000)throw Error('Alignment grid exceeds the browser budget. Shorten the read region.');
 const M=new Float32Array(cells).fill(NEG),I=new Float32Array(cells).fill(NEG),D=new Float32Array(cells).fill(NEG),pm=new Uint8Array(cells),pi=new Uint8Array(cells),pd=new Uint8Array(cells),mb=new Uint8Array(cells),ib=new Uint8Array(cells),il=new Uint8Array(cells),dl=new Uint16Array(cells);
 const inLocal=(p:number)=>cuts.some(c=>Math.abs(p-c)<=25),penalty=(p:number)=>{const distance=Math.min(...cuts.map(c=>Math.abs(p-c)));return 2.8+1.7*Math.min(1,distance/60)+options.alpha*2;};
 const extend=.015,insertExtend=.25,subPenalty=Math.log(3)+3;
 const insertionAllowed=(p:number)=>options.startRadius===undefined||cuts.some(c=>Math.abs(p-c)<=options.startRadius!),deletionAllowed=(start:number,end:number)=>options.startRadius===undefined||cuts.some(c=>start<=c+options.startRadius!&&end>=c-options.startRadius!);
 M[0]=0;for(let j=1;j<=Math.min(m,options.maxDeletion);j++){if(options.local&&!inLocal(anchor+j-1))break;D[j]=-penalty(anchor)-extend*(j-1);pd[j]=j===1?0:2;dl[j]=j;}
 let terminalScore=NEG,terminalI=n;
 for(let i=1;i<=n;i++){
  const row=i*cols,prev=(i-1)*cols,p=anchor+i-1,observed=profile[p]||[0,0,0,0],trained=options.training[p]||false;
  const maxScore=Math.max(...observed),prob=observed.map(v=>Math.exp(8*(v-maxScore))),sumProb=prob.reduce((a,b)=>a+b,0),logEmission=prob.map(v=>Math.log(4*v/sumProb));
  let bestBase=0;for(let b=1;b<4;b++)if(observed[b]>observed[bestBase])bestBase=b;
  const lo=Math.max(0,i-options.maxInsertion),hi=Math.min(m,i+options.maxDeletion);
  for(let j=lo;j<=hi;j++){
   const at=row+j,refP=anchor+j-1,allowed=!options.local||inLocal(refP);
   if(j>0){const before=prev+j-1;let state=0,value=M[before];if(I[before]>value){state=1;value=I[before];}if(D[before]>value&&deletionAllowed(anchor+j-1-dl[before],anchor+j-1)){state=2;value=D[before];}
    const base=DNA.indexOf(reference[refP]);let emitted=base>=0?logEmission[base]:0,chosen=base;
    if(trained&&allowed&&logEmission[bestBase]-subPenalty>emitted){emitted=logEmission[bestBase]-subPenalty;chosen=bestBase;}
    M[at]=value+emitted;pm[at]=state;mb[at]=chosen<0?4:chosen;
   }
   if(allowed&&j<m){const before=prev+j,open=penalty(anchor+j);let value=insertionAllowed(anchor+j)?M[before]-open:NEG,state=0,length=1;
    if(il[before]<options.maxInsertion&&I[before]-insertExtend>value){value=I[before]-insertExtend;state=1;length=il[before]+1;}
    I[at]=value;pi[at]=state;ib[at]=trained?bestBase:4;il[at]=length;
   }
   if(j>0&&allowed){const before=at-1,open=penalty(refP);let value=M[before]-open,state=0;let length=1;if(dl[before]<options.maxDeletion&&D[before]-extend>value){value=D[before]-extend;state=2;length=dl[before]+1;}D[at]=value;pd[at]=state;dl[at]=length;}
  }
  // A molecule can terminate before other pooled products; it emits no further signal.
  if(m===reference.length-anchor&&m>=lo&&m<=hi&&M[row+m]>terminalScore){terminalScore=M[row+m];terminalI=i;}
 }
 let i=n,j=0,state=0,best=NEG;
 for(let k=Math.max(0,n-options.maxInsertion);k<=Math.min(m,n+options.maxDeletion);k++)for(let st=0;st<2;st++){const v=(st?I:M)[n*cols+k];if(v>best){best=v;j=k;state=st;}}
 if(terminalScore>best){best=terminalScore;i=terminalI;j=m;state=0;}
 if(best<NEG/2)return null;const finalJ=j,ops:{kind:string;p:number;base:string}[]=[];
 while(i>0||j>0){const at=i*cols+j;if(state===0){if(i===0||j===0)return null;ops.push({kind:'M',p:anchor+j-1,base:DNA[mb[at]]||'N'});state=pm[at];i--;j--;}else if(state===1){if(i===0)return null;ops.push({kind:'I',p:anchor+j,base:DNA[ib[at]]||'N'});state=pi[at];i--;}else{if(j===0)return null;ops.push({kind:'D',p:anchor+j-1,base:''});state=pd[at];j--;}}
 ops.reverse();let sequence=reference.slice(0,anchor),events:PathEdit[]=[],active:PathEdit|undefined;
 const close=()=>{if(active){events.push(active);active=undefined;}};
 for(const op of ops){const equal=op.kind==='M'&&op.base===reference[op.p];if(equal){close();sequence+=op.base;continue;}if(!active)active={start:op.p,end:op.p,inserted:''};if(op.kind!=='I')active.end=op.p+1;if(op.kind!=='D'){active.inserted+=op.base;sequence+=op.base;}}
 close();sequence+=reference.slice(anchor+finalJ);
 // Preserve compound edits separated only by a short exact bridge.
 const merged:PathEdit[]=[];for(const e of events){const last=merged.at(-1);if(last&&e.start-last.end<=3){last.inserted+=reference.slice(last.end,e.start)+e.inserted;last.end=e.end;}else merged.push({...e});}
 if(merged.length>8)return null;return {sequence,events:merged,score:best};
}
/** Insertion site must be cut-proximal. A deletion can extend in either
 * direction, so its interval must intersect a cut-proximal initiation zone.
 * Equal-length substitutions do not have an indel origin. */
export function eventsRespectIndelOrigin(events:PathEdit[],cuts:number[],radius:number|undefined):boolean{
 return radius===undefined||events.every(e=>e.inserted.length===e.end-e.start||cuts.some(c=>e.start<=c+radius&&e.end>=c-radius));
}
export function materialize(reference:string,events:PathEdit[]):string|null {let at=0,result='';for(const e of [...events].sort((a,b)=>a.start-b.start)){if(e.start<at||e.end>reference.length||e.start<0)return null;result+=reference.slice(at,e.start)+e.inserted;at=e.end;}return result+reference.slice(at);}
export interface Prepared {name:string;sequence:string;cuts:number[];input:PooledAmplicon;reverse:boolean;offset:number}
export interface Model {reference:number;sequence:string;events:PathEdit[];aliases:{reference:number;events:PathEdit[]}[]}
export const failure=(reason:string):PooledPreflight=>({passed:false,reason,positions:[],modelCount:0,rank:0,condition:null,minimumDifferences:0,minimumCutOffset:0,amplicons:[]});
export function prepare(amplicons:PooledAmplicon[],primer:string,s:SangerAnalysisSettings){
 if(!/^[ACGT]{12,40}$/.test(primer)||amplicons.length<1||amplicons.length>4||new Set(amplicons.map(a=>a.name)).size!==amplicons.length)throw Error('Use 1–4 uniquely named submitted amplicons and the actual sequencing primer.');
 if(s.method!=='mixture'||s.signalModel==='classic')throw Error('Signal-guided alignment requires calibrated trace decomposition.');
 if(!Number.isInteger(s.window)||s.window<60||s.window>500||!Number.isInteger(s.readStart)||!Number.isInteger(s.readEnd)||s.readStart<0||s.readEnd>1000||s.readEnd<=s.readStart)throw Error('Invalid signal-guided window or usable read range.');
 if(s.editStartRadius!==undefined&&(!Number.isInteger(s.editStartRadius)||s.editStartRadius<0||s.editStartRadius>100))throw Error('Indel origin radius must be 0–100 bp, or unrestricted.');
 const maxDeletion=s.largeDeletions?(s.maxDeletion??200):s.maxIndel,maxInsertion=s.maxInsertion??30;
 if(!Number.isInteger(maxDeletion)||maxDeletion<1||maxDeletion>500||!Number.isInteger(maxInsertion)||maxInsertion<1||maxInsertion>60)throw Error('Use deletion limit 1–500 and insertion limit 1–60 bp.');
 const report=failure(''),refs:Prepared[]=amplicons.map(a=>{
  validateSangerProductBounds(a);const cuts=a.cuts??[a.cut];if(!a.name.trim()||!/^[ACGT]+$/.test(a.sequence)||a.sequence.length>20000||cuts.length<1||cuts.length>2||new Set(cuts).size!==cuts.length||cuts.some(c=>!Number.isInteger(c)||c<=0||c>=a.sequence.length))throw Error('Use one or two distinct cuts within each full submitted product.');
  const hits=exactBindingSites(primer,{name:a.name,sequence:a.sequence,cuts});if(hits.length!==1)throw Error(a.name+': sequencing primer must bind exactly once.');
  const h=hits[0],reverse=h.direction==='reverse',offset=reverse?a.sequence.length-h.start:h.end,sequence=(reverse?reverseComplement(a.sequence):a.sequence).slice(offset),positions=cuts.map(c=>(reverse?a.sequence.length-c:c)-offset).sort((a,b)=>a-b);
  if(positions[0]<s.readStart+100||positions.at(-1)!+30>s.readEnd)throw Error(a.name+': preserve ≥100 usable upstream bases and ≥30 downstream bases inside the Sanger range.');
  const windows=positions.map(c=>({cutOffset:c,windowStart:Math.max(s.readStart,c-Math.floor(s.window/2)),windowEnd:Math.min(s.readEnd,c+Math.ceil(s.window/2))}));
  report.amplicons.push({name:a.name,direction:h.direction,primerStart:h.start+1,primerEnd:h.end,cutOffset:positions[0],windowStart:windows[0].windowStart,windowEnd:windows.at(-1)!.windowEnd,targets:windows,maximumWindow:Math.min(500,Math.max(0,2*Math.min(...positions.map(c=>Math.min(c-s.readStart,s.readEnd-c))))),maximumIndependentWindow:positions.length===2?positions[1]-positions[0]:undefined,mergedWindows:windows.length===2&&windows[0].windowEnd>windows[1].windowStart});
  return {name:a.name,sequence,cuts:positions,input:a,reverse,offset};
 });
 const anchor=s.readStart+84,end=s.readEnd;
 const positions=Array.from({length:end-anchor},(_,i)=>i+anchor);if(positions.length<120)throw Error('Too little shared signal for signal-guided path inference.');
 const columns=refs.map(r=>positions.flatMap(p=>TRACE_BASES.map(b=>Number(r.sequence[p]===b))));
 let min=positions.length;for(let i=0;i<refs.length;i++)for(let j=0;j<i;j++)min=Math.min(min,positions.filter(p=>refs[i].sequence[p]!==refs[j].sequence[p]).length);
 const geometry=gramCondition(columns);if(min<3||geometry.rank<refs.length||geometry.condition===null||geometry.condition>100)throw Error('WT references are indistinguishable on the sequencing-primer axis.');
 report.passed=true;report.reason='WT primer geometry passed; edit paths are discovered from training signal and require independent QC.';report.positions=positions;report.modelCount=refs.length;report.minimumDifferences=min;report.rank=geometry.rank;report.condition=geometry.condition;report.joint={mode:s.jointEngine||'signal-path',candidateCount:refs.length,ambiguousClasses:0,coverageExcluded:0,maxDeletion};
 return {refs,report,anchor,end,maxDeletion,maxInsertion};
}
export function checkSignalPathSeparation(amplicons:PooledAmplicon[],primer:string,s:SangerAnalysisSettings){try{return prepare(amplicons,primer,s).report;}catch(e){return failure(e instanceof Error?e.message:String(e));}}
/** Large isolated peak-spacing jumps can mimic a biological phase shift. */
export function assertSignalTraceGrid(trace:SangerTrace):void{
 const spacing=Array.from({length:trace.peaks.length-1},(_,i)=>trace.peaks[i+1]-trace.peaks[i]);
 for(let i=8;i<spacing.length-8;i++){const local=spacing.slice(i-7,i+8).sort((a,b)=>a-b),median=local[7];if(median>0&&(spacing[i]>1.8*median||spacing[i]<.55*median))throw Error(trace.name+': irregular base-call cycle spacing. Re-basecall or inspect the chromatogram before interpreting a phase shift as an indel.');}
}
export interface SignalPathInput {amplicons:PooledAmplicon[];primer:string;ampliconConfirmed:boolean;settings:SangerAnalysisSettings;control:SangerTrace;sample:SangerTrace;synthetic?:boolean}
