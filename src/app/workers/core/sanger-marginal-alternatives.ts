import {originBaseQuantities} from './sanger-origin-base-support';
import type {SangerAnalysisResult,SangerMarginalCandidate,SangerMarginalPage} from '../../models/sanger-analysis';
import {sangerSequenceKey} from '../../utils/sanger-sequence-groups';
import {marginalSessions,marginalInputEvent,type MarginalEdge,type MarginalSession,type MarginalInternalGroup,type MarginalTargetState} from './sanger-marginal-engine';
import {MarginalLinearSolver,LP_INFINITY,type MarginalProgram,type LinearTerm,validMarginalSolution} from './sanger-marginal-lp';
import {materialize} from './sanger-path-common';
const DNA='ACGT',RESOLUTION=1e-7,PAGE_SIZE=10;
interface Node {edge:number;word:string;bound:number;order:number;nextBase?:number;probed?:boolean}
const better=(a:Node,b:Node)=>Math.abs(a.bound-b.bound)>2e-7?a.bound>b.bound:a.word.length!==b.word.length?a.word.length>b.word.length:a.order<b.order;
class Frontier {items:Node[]=[];push(x:Node){let i=this.items.push(x)-1;while(i){const p=(i-1)>>1;if(!better(this.items[i],this.items[p]))break;[this.items[i],this.items[p]]=[this.items[p],this.items[i]];i=p;}}pop(){const first=this.items[0],last=this.items.pop();if(this.items.length&&last){this.items[0]=last;let i=0;for(;;){const a=2*i+1,b=a+1;if(a>=this.items.length)break;const child=b<this.items.length&&better(this.items[b],this.items[a])?b:a;if(!better(this.items[child],this.items[i]))break;[this.items[i],this.items[child]]=[this.items[child],this.items[i]];i=child;}}return first;}}
interface Search {embeddingRows:Map<string,number[]>;embeddingEdges:Set<number>;selection?:import('../../models/sanger-analysis').SangerMarginalDisplayComponent['selection'];onlySub?:boolean;session:MarginalSession;target:MarginalTargetState;group:MarginalInternalGroup;program:MarginalProgram;solver:MarginalLinearSolver;z:number;rootIndex:number;counter:number;frontier:Frontier;seedFrontier:Frontier;seedInitialized:boolean;seedDone:boolean;seedPaused?:Node;seedPoint:number[];proposalIndex:number;proposalSolver?:MarginalLinearSolver;coordinates:number[];free:number[][];baseRows:Map<number,number[]>;edgeRows:number[];active:Set<number>;candidates:SangerMarginalCandidate[];seen:Set<string>;witnesses:Map<string,{x:number[];active:number[]}>;paused?:Node;exhausted:boolean}
function createSearch(session:MarginalSession,target:MarginalTargetState,group:MarginalInternalGroup):Search{
 const p:MarginalProgram={...target.program,upper:target.program.upper.slice(),cost:target.program.cost.map(()=>0),rows:target.program.rows.slice()},z=p.upper.length;p.upper.push(Math.min(1,group.fractionRange[1]));p.cost.push(0);
 const coordinates=[...group.bases.keys()].sort((a,b)=>a-b),baseRows=new Map<number,number[]>(),edgeRows:number[]=[],add=(terms:LinearTerm[])=>{const id=p.rows.length;p.rows.push({terms,lower:-LP_INFINITY,upper:LP_INFINITY});return id;};
 for(const q of coordinates)baseRows.set(q,group.bases.get(q)!.map(f=>add([...f,[z,-1]])));
 for(const e of group.edges)edgeRows.push(add([[e.mass,1],[z,-1]]));
 const free=group.edges.map(e=>coordinates.filter(q=>Math.abs(q-e.incoming-session.refs[e.reference].cuts[e.target])<=45||q>=e.start+e.incoming&&q<e.start+e.incoming+Math.max(0,e.size)));
 for(const coordinates of free)coordinates.sort((a,b)=>{const certainty=(q:number)=>Math.max(...group.bases.get(q)!.map(f=>f.reduce((n,[i,w])=>n+w*target.point[i],0)));return certainty(a)-certainty(b)||a-b;});
 const solver=new MarginalLinearSolver(p,session.optimizer);return{embeddingRows:new Map(),embeddingEdges:new Set(),session,target,group,program:p,solver,z,rootIndex:0,counter:0,frontier:new Frontier(),seedFrontier:new Frontier(),seedInitialized:false,seedDone:false,seedPoint:target.point,proposalIndex:0,coordinates,free,baseRows,edgeRows,active:new Set(),candidates:[],seen:new Set(),witnesses:new Map(),exhausted:false};
}
function fixedBase(s:Search,e:MarginalEdge,q:number){const ref=s.session.refs[e.reference],at=q<e.start+e.incoming?q-e.incoming:q-e.outgoing;return DNA.indexOf(ref.sequence[at]);}
/** Nonnegative quantities remaining after extracting one real allele from
 * this particular origin. Hall slacks preserve off-diagonal substitution
 * transport; group-wise base marginals alone do not prove an origin allele. */
function embeddingQuantities(s:Search,e:MarginalEdge,q:number,b:number){return originBaseQuantities(s.group,e,q,b,s.session.refs[e.reference].sequence);}
function alleleBaseMass(s:Search,e:MarginalEdge,q:number,b:number,x:number[]){return Math.min(x[e.mass],...embeddingQuantities(s,e,q,b).map(form=>form.reduce((n,[i,w])=>n+w*x[i],0)));}
function ensureEmbeddingRows(s:Search,edge:number){
 if(s.embeddingEdges.has(edge))return;const rows:MarginalProgram['rows']=[],offset=s.program.rows.length;
 for(const q of s.coordinates)for(let b=0;b<4;b++){const ids:number[]=[];for(const quantity of embeddingQuantities(s,s.group.edges[edge],q,b)){const coefficients=new Map<number,number>();for(const[i,w]of [...quantity,[s.z,-1]as LinearTerm])coefficients.set(i,(coefficients.get(i)||0)+w);ids.push(offset+rows.length);rows.push({terms:[...coefficients].filter(([,w])=>w!==0),lower:-LP_INFINITY,upper:LP_INFINITY});}s.embeddingRows.set(edge+':'+q+':'+b,ids);}
 s.solver.addRows(rows);s.embeddingEdges.add(edge);
}
function activateNode(s:Search,node:Node){
 ensureEmbeddingRows(s,node.edge);
 const e=s.group.edges[node.edge],free=s.free[node.edge],freeSet=new Set(free),active=new Set<number>([s.edgeRows[node.edge]]);
 for(const q of s.coordinates)if(!freeSet.has(q)){const b=fixedBase(s,e,q);if(b>=0){active.add(s.baseRows.get(q)![b]);for(const row of s.embeddingRows.get(node.edge+':'+q+':'+b)||[])active.add(row);}}
 for(let i=0;i<node.word.length;i++){const b=DNA.indexOf(node.word[i]);active.add(s.baseRows.get(free[i])![b]);for(const row of s.embeddingRows.get(node.edge+':'+free[i]+':'+b)||[])active.add(row);}
 const changes=[...new Set([...s.active,...active])].filter(i=>s.active.has(i)!==active.has(i)).sort((a,b)=>a-b);if(changes.length)s.solver.setRowBounds(changes,changes.map(i=>active.has(i)?0:-LP_INFINITY),changes.map(()=>LP_INFINITY));s.active=active;
}
function solveNode(s:Search,node:Node){activateNode(s,node);return s.solver.feasibleMaximum([[s.z,1]]);
}
function candidate(s:Search,node:Node,x:number[]):SangerMarginalCandidate|null{
 const e=s.group.edges[node.edge],ref=s.session.refs[e.reference],word=new Map(s.free[node.edge].map((p,i)=>[p,node.word[i]])),inserted=Array.from({length:Math.max(0,e.size)},(_,i)=>word.get(e.start+e.incoming+i)!).join(''),events:{start:number;end:number;inserted:string}[]=[];
 if(e.size)events.push({start:e.start,end:e.start+Math.max(0,-e.size),inserted});
 for(const[p,base]of word){if(p>=e.start+e.incoming&&p<e.start+e.incoming+Math.max(0,e.size))continue;const at=p<e.start+e.incoming?p-e.incoming:p-e.outgoing;if(base!==ref.sequence[at])events.push({start:at,end:at+1,inserted:base});}
 events.sort((a,b)=>a.start-b.start||a.end-b.end);if(!materialize(ref.sequence,events))throw Error('Inconsistent selected-group sequence');
 if(s.onlySub&&!events.length)return null;
 if(s.selection?.marker){const m=s.selection.marker,b=word.has(m.position)?DNA.indexOf(word.get(m.position)!):fixedBase(s,e,m.position);if(!m.bases.includes(b))return null;}
 if(s.selection?.origins&&!s.selection.origins.some(o=>o.edge===node.edge&&(!o.marker||DNA.indexOf(word.get(o.marker.position)||DNA[fixedBase(s,e,o.marker.position)])===o.marker.base)))return null;
 const inputEvents=events.map(event=>marginalInputEvent(ref,event.start,event.end,event.inserted,e.target)),key=sangerSequenceKey(ref.input.sequence,inputEvents);if(s.seen.has(key))return null;
 s.seen.add(key);const id=s.group.id+':candidate:'+s.candidates.length,companionChanges=[];
 for(const g of s.target.groups){if(g.id===s.group.id)continue;const mass=g.mass.reduce((n,[i,w])=>n+w*x[i],0),fraction=mass/s.session.shares[g.reference];let basesChanged=0;for(const forms of g.bases.values())if(forms.some(f=>Math.abs(f.reduce((n,[i,w])=>n+w*x[i],0)-f.reduce((n,[i,w])=>n+w*s.target.point[i],0))>1e-6))basesChanged++;if(Math.abs(mass-g.fraction)>1e-6||basesChanged)companionChanges.push({groupId:g.id,fraction,basesChanged});}
 s.witnesses.set(id,{x,active:[...s.active]});while(s.witnesses.size>20)s.witnesses.delete(s.witnesses.keys().next().value!);return{id,events:inputEvents,sequenceKey:key,fractionRange:[0,Math.min(1,node.bound/s.session.shares[e.reference])],fractionRangeKind:'conservative-outer-envelope',witnessFraction:x[s.z]/s.session.shares[e.reference],certificate:'measured-marginal-feasibility',companionChanges};
}
// Enumerate full words supported by the already fitted feasible allocation
// before requesting new allocations. This is initialization of the SAME model,
// shared by all inputs. No genotype list is constructed during basic inference.
// The full joint-LP frontier is retained for allocations outside this seed.
function seedBound(s:Search,node:Node){const e=s.group.edges[node.edge],free=s.free[node.edge],chosen=new Map(free.slice(0,node.word.length).map((q,i)=>[q,DNA.indexOf(node.word[i])])),freeSet=new Set(free),x=s.seedPoint,available=(q:number,b:number)=>alleleBaseMass(s,e,q,b,x);let bound=x[e.mass];for(const q of s.coordinates){const b=chosen.get(q),v=b!==undefined?available(q,b):freeSet.has(q)?Math.max(...[0,1,2,3].map(b=>available(q,b))):available(q,fixedBase(s,e,q));bound=Math.min(bound,v);}return Math.max(0,bound);}

/** Lazy best-bound prefix enumeration for ONE selected length family. Each
 * complete sequence has a feasible witness in the SAME measured joint model.
 * Cursor indexes emitted unique DNA sequences; it is never a total-search cap.
 * A time slice retains the frontier and returns incomplete, not "exhausted". */
export function marginalAlternativePage(result:SangerAnalysisResult,groupId:string,cursor=0,progress:(s:string)=>void=()=>{},timeSliceSeconds=20):SangerMarginalPage{
 const started=performance.now();if(!Number.isInteger(cursor)||cursor<0)throw Error('Invalid alternative cursor.');const session=marginalSessions.get(result);if(!session)throw Error('The measured analysis session is unavailable. Analyze the original traces before requesting alternatives.');const child=result.marginalGroups?.flatMap(g=>g.displayComponents||[]).find(c=>c.id===groupId),parentId=child?groupId.split('::profile:')[0]:groupId.endsWith('::sub')?groupId.slice(0,-5):groupId,target=session.targets.find(t=>t.groups.some(g=>g.id===parentId)),group=target?.groups.find(g=>g.id===parentId);if(!target||!group)throw Error('Select an existing reference/target length group.');
 let s=session.pages.get(groupId)as Search|undefined;if(!s){s=createSearch(session,target,group);s.onlySub=groupId.endsWith('::sub');s.selection=child?.selection;session.pages.set(groupId,s);}if(cursor>s.candidates.length)throw Error('Alternative cursor is beyond the computed frontier.');const wanted=cursor+PAGE_SIZE;let incomplete=false,message:string|undefined;
 const expired=()=>performance.now()-started>=timeSliceSeconds*1000;
 try{
  while(s.candidates.length<wanted&&!s.exhausted){
   if(expired()){incomplete=true;break;}
   if(!s.seedInitialized){for(let edge=0;edge<s.group.edges.length;edge++){const node={edge,word:'',bound:0,order:s.counter++};node.bound=seedBound(s,node);if(node.bound>RESOLUTION)s.seedFrontier.push(node);}s.seedInitialized=true;}
   if(!s.seedDone){const node=s.seedPaused||s.seedFrontier.pop();s.seedPaused=node;if(!node){s.seedDone=true;continue;}if(node.word.length===s.free[node.edge].length){const x=s.seedPoint.slice();x[s.z]=node.bound;activateNode(s,node);if(!validMarginalSolution(s.program,x))throw Error('Seed full-word feasibility check failed');const c=candidate(s,{...node,bound:s.group.fractionRange[1]},x);if(c)s.candidates.push(c);s.seedPaused=undefined;continue;}for(let b=node.nextBase||0;b<4;b++){node.nextBase=b;const child={edge:node.edge,word:node.word+DNA[b],bound:0,order:s.counter++};child.bound=seedBound(s,child);if(child.bound>RESOLUTION)s.seedFrontier.push(child);node.nextBase=b+1;}s.seedPaused=undefined;continue;}
   // Reuse marginal-bound objectives to propose further feasible allocations.
   // This changes objective costs in the original joint LP, rather than pinning
   // hundreds of full-word constraints during every proposal. Each allocation
   // is consumed lazily and every emitted full word is checked. The exhaustive
   // prefix frontier below is still retained after these proposals.
   if(s.proposalIndex<s.coordinates.length*4){if(!s.proposalSolver)s.proposalSolver=new MarginalLinearSolver({...s.target.program,upper:s.target.program.upper.slice(),cost:s.target.program.cost.map(()=>0),rows:s.target.program.rows.slice()},s.session.optimizer);const q=s.coordinates[Math.floor(s.proposalIndex/4)],b=s.proposalIndex%4,allocation=s.proposalSolver.solve(s.group.bases.get(q)![b],true);s.proposalIndex++;s.seedPoint=allocation.x;s.seedInitialized=false;s.seedDone=false;continue;}
   if(s.rootIndex<s.group.edges.length){const node={edge:s.rootIndex,word:'',bound:0,order:s.counter++};s.paused=node;const solved=solveNode(s,node);if(solved&&solved.value>RESOLUTION){node.bound=solved.value;s.frontier.push(node);}s.rootIndex++;s.paused=undefined;continue;}
   const node=s.paused||s.frontier.pop();s.paused=undefined;if(!node){s.exhausted=true;break;}s.paused=node;
   if(node.word.length===s.free[node.edge].length){const solved=solveNode(s,node);if(solved&&solved.value>RESOLUTION){const c=candidate(s,node,solved.x);if(c)s.candidates.push(c);}s.paused=undefined;continue;}
   // Propose a complete word directly from a feasible allocation, without
   // claiming that its largest possible mass has been optimized. The safe
   // witness mass is the minimum supported chosen-base/geometry flow. Every
   // full-word constraint is checked, and the relaxed prefix bound remains a
   // conservative upper bound. Prefix expansion is retained without pruning.
   if(!node.probed){const supported=solveNode(s,node);if(!supported){s.paused=undefined;continue;}const free=s.free[node.edge],word=node.word+free.slice(node.word.length).map(q=>{const values=[0,1,2,3].map(b=>alleleBaseMass(s,s.group.edges[node.edge],q,b,supported.x));return DNA[values.indexOf(Math.max(...values))];}).join(''),leaf={...node,word},x=supported.x.slice();let safe=x[s.group.edges[node.edge].mass];const chosen=new Map(free.map((q,i)=>[q,word[i]]));for(const q of s.coordinates){const b=chosen.has(q)?DNA.indexOf(chosen.get(q)!):fixedBase(s,s.group.edges[node.edge],q);if(b>=0)safe=Math.min(safe,alleleBaseMass(s,s.group.edges[node.edge],q,b,x));}x[s.z]=Math.max(0,safe);activateNode(s,leaf);if(safe>RESOLUTION){if(!validMarginalSolution(s.program,x))throw Error('Full-word candidate feasibility check failed');const c=candidate(s,leaf,x);if(c)s.candidates.push(c);}node.probed=true;if(s.candidates.length>=wanted)break;if(expired()){incomplete=true;break;}}
   for(let b=node.nextBase||0;b<4;b++){node.nextBase=b;const child={edge:node.edge,word:node.word+DNA[b],bound:0,order:s.counter++},solved=solveNode(s,child);if(solved&&solved.value>RESOLUTION){child.bound=solved.value;s.frontier.push(child);}node.nextBase=b+1;if(expired()){incomplete=true;break;}}
   if(node.nextBase===4)s.paused=undefined;else if(incomplete)break;
   if(s.counter%32===0)progress('Exploring selected group only · '+s.candidates.length+' sequences · '+s.frontier.items.length+' pending prefixes');
  }
 }catch(error){incomplete=true;message=error instanceof Error?error.message:String(error);}
 const candidates=s.candidates.slice(cursor,wanted),next=cursor+candidates.length,complete=candidates.length===PAGE_SIZE||s.exhausted;return{groupId,cursor,previousCursor:cursor?Math.max(0,cursor-PAGE_SIZE):null,nextCursor:s.exhausted&&next>=s.candidates.length?null:next,pageSize:10,resolution:RESOLUTION,candidates,complete:complete&&!incomplete,exhausted:s.exhausted,seconds:(performance.now()-started)/1000,message:message||(incomplete?'Search is incomplete; the frontier is retained. Continue to calculate the next candidates.':undefined)};
}
/** Diagnostic access for feasibility regression tests; never used by GT scoring. */
export function verifyFittedMarginalRepresentatives(result:SangerAnalysisResult){
 const session=marginalSessions.get(result);if(!session)throw Error('Live marginal session required');
 return(result.marginalGroups||[]).filter(g=>g.fraction>=.02).map(display=>{
  const target=session.targets.find(t=>t.groups.some(g=>g.id===display.id))!,group=target.groups.find(g=>g.id===display.id)!,s=createSearch(session,target,group),expected=display.representative;
  try{if(!expected)return{groupId:display.id,valid:false};const ref=session.refs[group.reference],key=sangerSequenceKey(ref.input.sequence,expected.events),x=[...target.point,expected.witnessFraction*session.shares[group.reference]];
   for(let edge=0;edge<group.edges.length;edge++){
    const word=s.free[edge].map(q=>{const v=group.bases.get(q)!.map(f=>f.reduce((n,[i,w])=>n+w*target.point[i],0));return DNA[v.indexOf(Math.max(...v))];}).join(''),node={edge,word,bound:x[s.z],order:0};activateNode(s,node);
    if(!validMarginalSolution(s.program,x))continue;const c=candidate(s,node,x);if(c?.sequenceKey===key)return{groupId:display.id,valid:true,witnessFraction:expected.witnessFraction};
   }return{groupId:display.id,valid:false};
  }finally{s.solver.dispose();}
 });
}
export function marginalCandidateWitness(result:SangerAnalysisResult,groupId:string,candidateId:string){const s=marginalSessions.get(result)?.pages.get(groupId)as Search|undefined,w=s?.witnesses.get(candidateId);return s&&w?{program:{...s.program,rows:[...s.target.program.rows,...w.active.map(i=>({terms:s.program.rows[i].terms,lower:0,upper:LP_INFINITY}))]},x:w.x}:null;}
