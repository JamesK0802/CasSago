import {conditionalZeroWitnesses} from './sanger-conditional-zero-witness';
import {expandInsertionTransports} from './sanger-insertion-transport';
import type {SangerMarginalGroup,SangerMarginalDisplayComponent,SangerJointEvent} from '../../models/sanger-analysis';
import type {MarginalInternalGroup,MarginalTargetState,MarginalSession,MarginalEdge} from './sanger-marginal-engine';
import {marginalInputEvent} from './sanger-marginal-engine';
import {MarginalLinearSolver,type LinearTerm} from './sanger-marginal-lp';
import {sangerEditedSequence} from '../../utils/sanger-sequence-groups';

const DNA='ACGT',EPS=1e-7;
const sum=(terms:LinearTerm[],x:number[])=>terms.reduce((n,[i,w])=>n+w*x[i],0);
interface Atom {edge:MarginalEdge;mass:number;word:string;marker?:{q:number;b:number}}
interface Part {terms:LinearTerm[];atoms:Atom[];events:SangerJointEvent[];mixed:boolean;marker?:{q:number;bases:number[]};range?:[number,number]}
/** Compact DISPLAY partition of the existing fitted allocation, never an
 * exhaustive haplotype search. No input, measurement, envelope or fit changes.
 * Multiple variable insertion sites remain a mixture (no invented phasing).
 * Substituted indels are split only in the one-variable-site case below.
 */
export function withinLengthDisplay(session:MarginalSession,target:MarginalTargetState,g:MarginalInternalGroup,display:SangerMarginalGroup):SangerMarginalDisplayComponent[]|undefined{
 if(!g.size||g.fraction<=EPS)return;
 const x=target.point,ref=session.refs[g.reference],share=session.shares[g.reference],event=(a:number,b:number,s:string)=>marginalInputEvent(ref,a,b,s,g.target);
 const active=g.edges.filter(e=>x[e.mass]>EPS),parts=new Map<string,Part>();
 // One genuinely mixed output site has no inter-site phase ambiguity. Its
 // nucleotide marginal itself defines a disjoint, exhaustive partition.
 const mixed=[...g.bases].filter(([,forms])=>forms.filter(f=>sum(f,x)>EPS).length>1);
 if(active.length===1&&mixed.length===1&&display.representative){
  const e=active[0],[q,forms]=mixed[0],insertion=q>=e.start+e.incoming&&q<e.start+e.incoming+Math.max(0,e.size),at=q<e.start+e.incoming?q-e.incoming:q-e.outgoing;
  const rows:Part[]=[];
  for(let b=0;b<4;b++)if(sum(forms[b],x)>EPS){const es=display.representative.events.map(e=>({...e}));if(insertion){const indel=es.find(e=>e.inserted-(e.end-e.start)===g.size);if(!indel?.insertedSequence)continue;const j=ref.reverse?e.size-1-(q-e.start-e.incoming):q-e.start-e.incoming,base=ref.reverse?DNA[3-b]:DNA[b];indel.insertedSequence=indel.insertedSequence.slice(0,j)+base+indel.insertedSequence.slice(j+1);}else{const change=event(at,at+1,DNA[b]);for(let i=es.length-1;i>=0;i--)if(es[i].start===change.start&&es[i].end===change.end&&es[i].inserted===1)es.splice(i,1);if(DNA[b]!==ref.sequence[at])es.push(change);}
   rows.push({terms:forms[b],atoms:[],events:es,mixed:false});
  }
  if(rows.length>1){const zero=conditionalZeroWitnesses(target.program,x,g,rows.map(p=>p.terms),session.optimizer);let optional=false;try{optional=rows.some((_,i)=>zero.possible(i));}finally{zero.dispose();}if(!optional){const solver=new MarginalLinearSolver(target.program,session.optimizer);try{for(const p of rows)p.range=solver.bounds(p.terms);if(rows.every(p=>p.range![0]>EPS)){
     const answer=rows.map((p,i)=>component(p,i,true,q,forms.findIndex(f=>f===p.terms))),other=Array.from({length:4},(_,b)=>b).filter(b=>!rows.some(p=>p.terms===forms[b]));
     if(other.length){const terms=other.flatMap(b=>forms[b]),range=solver.bounds(terms);if(range[1]>EPS)answer.push(component({terms,atoms:[],events:display.representative.events,mixed:true,marker:{q,bases:other},range},answer.length,true,q,-1));}
     return answer;
    }}finally{solver.dispose();}}}
 }
 // For a substitution-free fitted allocation, geometry and inserted-base
 // flows already give the component proportions, without fitting haplotypes.
 if([...g.flow?.values()||[]].some(f=>f.changed.some(i=>x[i]>1e-8)))return;
 // Two varying output sites do not identify insertion haplotypes. Do not
 // split them into geometry-dependent words that could imply false linkage.
 if(g.size>0&&mixed.length>1)return;
 if(g.sharedInsertions&&active.length>1){const expanded=expandInsertionTransports(target,g);return withinLengthDisplay(session,expanded.target,expanded.group,display);}
 for(const e of active){
  const free=[...e.free].sort((a,b)=>a[0]-b[0]),variable=free.filter(([,ids])=>ids.filter(i=>x[i]>EPS).length>1),word=free.map(([,ids])=>DNA[ids.map(i=>x[i]).indexOf(Math.max(...ids.map(i=>x[i]))) ]).join('');
  const branches=variable.length===1?variable[0][1].flatMap((id,b)=>x[id]>EPS?[{terms:[[id,1]]as LinearTerm[],mass:x[id],marker:{q:variable[0][0],b},word:word.slice(0,variable[0][0]-e.start-e.incoming)+DNA[b]+word.slice(variable[0][0]-e.start-e.incoming+1)}]:[]):[{terms:[[e.mass,1]]as LinearTerm[],mass:x[e.mass],word,marker:undefined}];
  for(const branch of branches){const es=[event(e.start,e.start+Math.max(0,-e.size),branch.word)],key=sangerEditedSequence(ref.input.sequence,es)+(variable.length>1?'::mixture':'');let p=parts.get(key);if(!p){p={terms:[],atoms:[],events:es,mixed:variable.length>1};parts.set(key,p);}p.terms.push(...branch.terms);p.atoms.push({edge:e,mass:branch.mass,word:branch.word,marker:branch.marker});}
 }
 if(parts.size<2)return;
 const rows=[...parts.values()];let solver:MarginalLinearSolver|undefined;const exact=()=>solver??=new MarginalLinearSolver(target.program,session.optimizer);
 // Complete the partition with a compact residual, including zero fitted
 // allocations that remain possible. They are never discarded or enumerated.
 const restAtoms:Atom[]=[],rest:LinearTerm[]=[];
 for(const e of g.edges){const taken=rows.flatMap(p=>p.atoms).filter(a=>a.edge===e);if(taken.some(a=>!a.marker))continue;
  if(!taken.length){rest.push([e.mass,1]);restAtoms.push({edge:e,mass:x[e.mass],word:'A'.repeat(Math.max(0,e.size))});continue;}
  const q=taken[0].marker!.q;for(let b=0;b<4;b++)if(!taken.some(a=>a.marker!.b===b)){const id=e.free.get(q)![b];rest.push([id,1]);restAtoms.push({edge:e,mass:x[id],word:'A'.repeat(Math.max(0,e.size)),marker:{q,b}});}
 }
 if(rest.length)rows.push({terms:rest,atoms:restAtoms,events:display.representative?.events||rows[0].events,mixed:true});
 const zero=conditionalZeroWitnesses(target.program,x,g,rows.map(p=>p.terms),session.optimizer);
 try{
  const supported:Part[]=[];
  for(let i=0;i<rows.length;i++){const p=rows[i];if(p.mixed)continue;const lo=zero.possible(i)?0:Math.max(0,exact().solve(p.terms).value-1e-7);if(lo>EPS){p.range=[lo,Math.max(0,exact().solve(p.terms,true).value+1e-7)];supported.push(p);}}
  if(!supported.length)return;
  const uncertain=rows.filter(p=>!supported.includes(p));
  if(uncertain.length){const terms=uncertain.flatMap(p=>p.terms),range=exact().bounds(terms);if(range[1]>EPS){const best=uncertain.reduce((a,b)=>sum(a.terms,x)>sum(b.terms,x)?a:b);supported.push({terms,atoms:uncertain.flatMap(p=>p.atoms),events:best.events,mixed:true,range});}}
  if(supported.length<2)return;
  return supported.map((p,i)=>component(p,i));
 }finally{zero.dispose();solver?.dispose();}

 function component(p:Part,i:number,markerOnly=false,qMarker=0,bMarker=0):SangerMarginalDisplayComponent{
  const mass=sum(p.terms,x),range=p.range!,fhi=display.fractionRange[1],lo=range[0]/share,hi=range[1]/share;
  const bases=display.bases.map(parent=>{
   const q=ref.reverse?ref.input.sequence.length-ref.offset-parent.position-1:parent.position-ref.offset,values=[0,0,0,0];
   if(markerOnly){const b=ref.reverse?3-bMarker:bMarker;for(let j=0;j<4;j++)values[j]=q===qMarker&&bMarker>=0?Number(j===b):parent.estimate[j];}
   else{for(const a of p.atoms){if(a.mass<=EPS)continue;const e=a.edge,vars=e.free.get(q);if(vars){if(a.marker?.q===q)values[a.marker.b]+=a.mass;else for(let b=0;b<4;b++)values[b]+=a.mass*x[vars[b]]/x[e.mass];}else{const at=q<e.start+e.incoming?q-e.incoming:q-e.outgoing,b=DNA.indexOf(ref.sequence[at]);if(b>=0)values[b]+=a.mass;}}for(let b=0;b<4;b++)values[b]/=Math.max(EPS,mass);if(ref.reverse)values.reverse();}
   const lower=parent.lower.map(v=>Math.max(0,Math.min(1,1-fhi*(1-v)/Math.max(EPS,lo)))),upper=parent.upper.map(v=>Math.max(0,Math.min(1,fhi*v/Math.max(EPS,lo))));
   if(markerOnly&&q===qMarker&&bMarker>=0){for(let b=0;b<4;b++)lower[b]=upper[b]=values[b];}
   if(mass<=EPS){for(let b=0;b<4;b++){values[b]=0;lower[b]=0;upper[b]=1;}}
   return{position:parent.position,estimate:values,lower,upper};
  });
  return{id:g.id+'::profile:'+i,fraction:mass/share,fractionRange:[lo,Math.min(1,hi)],bases,events:p.events,mixture:p.mixed,selection:markerOnly?{marker:{position:qMarker,bases:bMarker>=0?[bMarker]:p.marker!.bases}}:{origins:p.atoms.map(a=>({edge:g.edges.indexOf(a.edge),marker:a.marker?{position:a.marker.q,base:a.marker.b}:undefined}))}};
 }
}
