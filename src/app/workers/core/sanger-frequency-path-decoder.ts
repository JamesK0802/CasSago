import {frequencySubsetIndex} from'./sanger-subset-index';
import {assignSequenceBases,SangerLinearOptimizer} from'./sanger-discrete-phase';
import{Prepared,PathEdit,materialize}from'./sanger-path-common';import{normalizeMoleculeEdits}from'./sanger-edit-normalize';import{SangerFitCache}from'./sanger-fit-cache';
interface Atom {size:number;fraction:number}
export interface FrequencyRefinementMemo{context?:string;states:Map<string,Atom[]>}
interface Family extends Atom {reference:number}
interface Allocation {cost:number;code:number}
const bases='ACGT';
function subsetSums(atoms:Atom[]){if(2**atoms.length>1048576)throw Error("Frequency subset work budget exceeded; unresolved signal");const sums=[0];for(const a of atoms){const n=sums.length;for(let i=0;i<n;i++)sums.push(sums[i]+a.fraction);}return sums;}
const nearest=(values:number[],v:number)=>values.reduce((best,x)=>Math.min(best,Math.abs(v-x)),Infinity);
function refineAtoms(input:Atom[],levels:number[],tolerance:number){
 let atoms=input.map(a=>({...a}));const families=[...new Set(atoms.map(a=>a.size))],anchors=new Map(families.map(size=>[size,atoms.findIndex(a=>a.size===size)])),variables=atoms.flatMap((a,i)=>anchors.get(a.size)===i?[]:[i]);if(!variables.length)return atoms;
 const totals=new Map(families.map(size=>[size,atoms.filter(a=>a.size===size).reduce((s,a)=>s+a.fraction,0)]));
 for(let iteration=0;iteration<12;iteration++){const sums=subsetSums(atoms),lookup=frequencySubsetIndex(sums),n=variables.length,matrix=Array.from({length:n},()=>Array(n+1).fill(0));let count=0;
  for(const value of levels){const mask=lookup.mask(value);if(Math.abs(sums[mask]-value)>2*tolerance)continue;const row=variables.map(i=>Number(!!(mask&(1<<i)))-Number(!!(mask&(1<<anchors.get(atoms[i].size)!)))),base=families.reduce((s,size)=>s+((mask&(1<<anchors.get(size)!))?totals.get(size)!:0),0),rhs=value-base;for(let i=0;i<n;i++){matrix[i][n]+=row[i]*rhs;for(let j=0;j<n;j++)matrix[i][j]+=row[i]*row[j];}count++;}if(count<n)break;
  let valid=true;for(let c=0;c<n;c++){let pivot=c;for(let r=c+1;r<n;r++)if(Math.abs(matrix[r][c])>Math.abs(matrix[pivot][c]))pivot=r;if(Math.abs(matrix[pivot][c])<1e-10){valid=false;break;}[matrix[c],matrix[pivot]]=[matrix[pivot],matrix[c]];const d=matrix[c][c];for(let j=c;j<=n;j++)matrix[c][j]/=d;for(let r=0;r<n;r++)if(r!==c){const f=matrix[r][c];for(let j=c;j<=n;j++)matrix[r][j]-=f*matrix[c][j];}}
  if(!valid)break;const next=atoms.map(a=>({...a}));variables.forEach((id,i)=>next[id].fraction=matrix[i][n]);for(const size of families)next[anchors.get(size)!].fraction=totals.get(size)!-variables.filter(i=>atoms[i].size===size).reduce((s,i)=>s+next[i].fraction,0);if(next.some(a=>a.fraction<=0))break;const change=Math.max(...next.map((a,i)=>Math.abs(a.fraction-atoms[i].fraction)));atoms=next;if(change<1e-9)break;
 }
 return atoms;
}
/** Decode every group count with the same squared-error integer model.
 * Origin/base variables represent inserted words without enumerating 4^length
 * strings. Group count never selects a different solver or loss function. */
function decodePaths(ref:Prepared,atoms:Atom[],y:number[][],radius:number,span:number,maxVariants:number,distinctGroups=false,optimizer?:SangerLinearOptimizer){
 if(!optimizer)throw Error('Sequence profile decoding requires the common integer optimizer.');
 const cut=ref.cuts[0],groups=atoms.map(a=>({events:[{start:cut,end:cut+(a.size<0?-a.size:a.size===0?span:0),inserted:a.size>0?'A'.repeat(a.size):a.size===0?ref.sequence.slice(cut,cut+span):''}],size:a.size,sequence:ref.sequence,editable:true,weights:new Map([[0,a.fraction]])}));
 // A zero-net replacement may equal WT. One unrestricted model therefore
 // already contains every former choice of which frequency slot is WT.
 // Solving it once preserves that entire union without repeated WT branches.
 const events=assignSequenceBases(ref.sequence,cut,groups,Array.from({length:y.length},(_,j)=>cut-radius+j),y,optimizer,undefined,120,maxVariants);
 if(!events)return undefined;const sequences=events.map(es=>materialize(ref.sequence,es)!);if(distinctGroups&&new Set(sequences).size!==sequences.length)return undefined;
 let cost=0;for(let j=0;j<y.length;j++)for(let b=0;b<4;b++)cost+=(sequences.reduce((n,seq,i)=>n+Number(seq[cut-radius+j]===bases[b])*atoms[i].fraction,0)-y[j][b])**2;
 return{cost,groups:events.map((events,i)=>({events,fraction:atoms[i].fraction}))};
}

/** Frequency refinement never merges groups because their percentages are
 * alike. Observational levels propose splits of a length class; the complete
 * sequence model and raw reconstruction must independently validate them. */
export function decodeFrequencyPaths(refs:Prepared[],families:Family[],positions:number[],observed:number[],noise:number,radius=3,span=10,maxVariants=Infinity,maxStates=10000,progress:(stage:string)=>void=()=>{},collectCandidates=false,requestedCounts?:Record<string,number>,optimizer?:SangerLinearOptimizer,minimumFraction=.01,memo?:FrequencyRefinementMemo){
 const index=new Map(positions.map((p,i)=>[p,i])),groups:{reference:number;events:PathEdit[];fraction:number}[]=[],details:any[]=[],candidates:{reference:number;groups:{events:PathEdit[];fraction:number}[];rmse:number;maximum:number}[]=[],candidateKeys=new Set<string>();
 for(let h=0;h<refs.length;h++){
  progress(`Decoding offset edits · reference ${h+1} / ${refs.length}`);
  const ref=refs[h],cut=ref.cuts[0],begin=cut-radius,end=cut+radius+Math.max(span,...families.filter(f=>f.reference===h).map(f=>f.size)),own=families.filter(f=>f.reference===h&&f.fraction>0),share=own.reduce((s,a)=>s+a.fraction,0),tolerance=Math.max(1e-5,4*noise),minimum=share*minimumFraction;
  const other=(p:number)=>{const v=[0,0,0,0];for(const f of families)if(f.reference!==h&&f.fraction>0){const r=refs[f.reference],c=r.cuts[0],b=bases.indexOf(r.sequence[p<c-radius?p:p-f.size]);if(b>=0)v[b]+=f.fraction;}return v;};
  const y=Array.from({length:end-begin},(_,j)=>{const p=begin+j,k=index.get(p);if(k===undefined)throw Error('Missing local trace coverage.');for(const [i,r]of refs.entries())if(i!==h&&p>=r.cuts[0]-radius&&p<r.cuts[0]+radius+Math.max(span,10))throw Error('Local windows overlap.');const bg=other(p);return bg.map((v,b)=>observed[4*k+b]-v);}),values=y.flat().filter(v=>v>tolerance&&v<share-tolerance).sort((a,b)=>a-b),clusters:number[][]=[];
  for(const v of values){const last=clusters.at(-1);if(last&&v-last[0]<=tolerance/2)last.push(v);else clusters.push([v]);}const levels=clusters.map(c=>c.reduce((a,b)=>a+b,0)/c.length),key=(a:Atom[])=>a.map(x=>x.size+':'+x.fraction.toFixed(7)).sort().join('/'),ordered=(a:Atom[])=>a.sort((a,b)=>a.size-b.size||b.fraction-a.fraction);
  interface State {atoms:Atom[];missing:number[];error:number}
  const context=JSON.stringify([ref.sequence,cut,radius,span,tolerance,values]),rawSeen=new Set<string>(),refined=memo?.context===context?memo.states:new Map<string,Atom[]>();if(memo&&memo.context!==context){memo.context=context;memo.states=refined;}
  const refineCached=(atoms:Atom[],raw:string)=>{let result=refined.get(raw);if(!result){result=refineAtoms(atoms,values,tolerance);if(refined.size>=20000)refined.delete(refined.keys().next().value!);refined.set(raw,result.map(a=>({...a})));}return result.map(a=>({...a}));};

  const queue:State[]=[],seen=new Set<string>();const enqueue=(atoms:Atom[])=>{if(requestedCounts&&atoms.some(a=>atoms.filter(b=>b.size===a.size).length>(requestedCounts[String(a.size)]||0)))return;atoms=ordered(atoms);const raw=JSON.stringify(atoms.map(a=>[a.size,a.fraction]));if(rawSeen.has(raw))return;rawSeen.add(raw);atoms=ordered(refineCached(atoms,raw));const k=key(atoms);if(seen.has(k))return;seen.add(k);const sums=subsetSums(atoms),lookup=frequencySubsetIndex(sums),errors=levels.map(v=>lookup.distance(v));queue.push({atoms,missing:levels.filter((_,i)=>errors[i]>tolerance),error:errors.reduce((s,e)=>s+e*e,0)});};enqueue(own.map(({size,fraction})=>({size,fraction})));let accepted:{events:PathEdit[];fraction:number}[]|undefined,states=0,decoded=0,bestError=Infinity,bestState:State|undefined;
  while(queue.length&&states++<maxStates){if(states%32===0)progress("Frequency allocations · "+states+" states · "+queue.length+" pending");queue.sort((a,b)=>requestedCounts?a.missing.length-b.missing.length||a.error-b.error||b.atoms.length-a.atoms.length:a.atoms.length-b.atoms.length||a.missing.length-b.missing.length||a.error-b.error);const state=queue.shift()!;if(!bestState||state.missing.length<bestState.missing.length||state.missing.length===bestState.missing.length&&state.error<bestState.error)bestState=state;
   if((!requestedCounts||Object.entries(requestedCounts).every(([size,count])=>state.atoms.filter(a=>a.size===Number(size)).length===count))&&(collectCandidates||!state.missing.length||state.atoms.length===maxVariants+1)){const candidate=decodePaths(ref,state.atoms,y,radius,span,maxVariants,!!requestedCounts,optimizer);decoded++;if(candidate){const merged=new Map<string,{events:PathEdit[];sequence:string}>();for(const g of candidate.groups){const sequence=materialize(ref.sequence,g.events)!;if(!merged.has(sequence))merged.set(sequence,{events:g.events,sequence});}const models=[...merged.values()],rows=positions.flatMap((p,i)=>refs.every((r,k)=>k===h||p<r.cuts[0]-radius||p>=r.cuts[0]+radius+Math.max(span,10))?[i]:[]),columns=models.map(m=>rows.flatMap(i=>[...bases].map(b=>Number(m.sequence[positions[i]]===b)))),target=rows.flatMap(i=>{const bg=other(positions[i]);return bg.map((v,b)=>(observed[4*i+b]-v)/share);}),fit=new SangerFitCache(columns,target).fit(models.map((_,i)=>i)),weights=fit.weights.map(w=>w*share);let maximum=0;for(let j=0;j<y.length;j++){const p=begin+j,v=[0,0,0,0];models.forEach((m,i)=>{const b=bases.indexOf(m.sequence[p]);if(b>=0)v[b]+=weights[i];});maximum=Math.max(maximum,Math.sqrt(v.reduce((s,x,b)=>s+(x-y[j][b])**2,0)/4));}bestError=Math.min(bestError,maximum);if(collectCandidates||maximum<=tolerance){accepted=models.flatMap((m,i)=>(collectCandidates||weights[i]>1e-9)?[{events:m.events,fraction:weights[i]}]:[]);if(!collectCandidates)break;const key=models.map(m=>m.sequence).sort().join('|');if(!candidateKeys.has(key)){candidateKeys.add(key);candidates.push({reference:h,groups:accepted,rmse:fit.rmse,maximum});}}}}
   if(state.atoms.length>=maxVariants+1||requestedCounts&&Object.entries(requestedCounts).every(([size,count])=>state.atoms.filter(a=>a.size===Number(size)).length===count))continue;const targets=state.missing.length?(requestedCounts?state.missing:[state.missing[0]]):levels;
   for(const value of targets)for(let i=0;i<state.atoms.length;i++){const atom=state.atoms[i];if(requestedCounts&&state.atoms.filter(a=>a.size===atom.size).length>=(requestedCounts[String(atom.size)]||0))continue;const others=state.atoms.filter((_,j)=>j!==i);for(const sum of subsetSums(others)){const u=value-sum;if(u<minimum||atom.fraction-u<minimum)continue;enqueue([...others,{size:atom.size,fraction:u},{size:atom.size,fraction:atom.fraction-u}]);}}
  }
  details.push({reference:ref.name,states,decoded,bestError,accepted:!collectCandidates&&!!accepted,proposalOnly:collectCandidates});if(!accepted)throw Error('Frequency/path decoding unresolved for '+ref.name+' after '+states+' states; best local error '+bestError+' '+JSON.stringify({tolerance,levels,own,bestState}));groups.push(...accepted.map(g=>({...g,reference:h})));
 }
 return{groups,details,candidates};
}
