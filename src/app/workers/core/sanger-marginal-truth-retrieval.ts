import {originContextContains} from './sanger-origin-prefix';
import type {SangerAnalysisResult,SangerAnnotation,SangerTruthRetrieval} from '../../models/sanger-analysis';
import {reverseComplement} from '../../utils/biology.utils';
import {sangerEditedSequence,sangerSequenceKey} from '../../utils/sanger-sequence-groups';
import {targetEvents} from '../../utils/sanger-target-profile';
import {marginalSessions,type MarginalInternalGroup,type MarginalTargetState,type MarginalSession} from './sanger-marginal-engine';
import {MarginalLinearSolver,type MarginalProgram,type LinearTerm} from './sanger-marginal-lp';

interface Profile{sequenceKey:string;fraction:number;sequence:string}
interface Family{group:MarginalInternalGroup;reference:string;target:number;share:number;profiles:Profile[];truthMass:number}
const DNA='ACGT';
/** POST-INFERENCE hypothesis check only. Freeze group creation/ranking first.
 * Reconstruct the complete requested within-family allele mixture, including
 * its proportions, rather than accepting independent coordinate bounds. */
export function retrieveMarginalTruth(result:SangerAnalysisResult,truth:SangerAnnotation[],progress:(s:string)=>void=()=>{}):SangerTruthRetrieval{
 const began=performance.now(),session=marginalSessions.get(result);if(!session)throw Error('Analyze the original traces before checking GT retrieval.');const checks:SangerTruthRetrieval['checks']=[];
 for(const target of session.targets){
  const families:Family[]=[];
  for(const g of target.groups){const ref=session.refs[g.reference],t=ref.reverse?ref.cuts.length-g.target-1:g.target,a=truth.find(a=>a.name===ref.name);if(!a)continue;const map=new Map<string,Profile>();
   for(const c of a.components){const es=targetEvents(c.events||[],t),size=c.pathSizes?.[t]??es.reduce((n,e)=>n+e.inserted-(e.end-e.start),0);if(size!==g.size||c.fraction<=0||c.unresolved)continue;const key=sangerSequenceKey(a.reference,es),old=map.get(key);if(old){old.fraction+=c.fraction;continue;}const dna=sangerEditedSequence(a.reference,es),forward=ref.reverse?reverseComplement(dna):dna;map.set(key,{sequenceKey:key,fraction:c.fraction,sequence:forward.slice(ref.offset)});}
   if(map.size)families.push({group:g,reference:ref.name,target:t+1,share:session.shares[g.reference],profiles:[...map.values()],truthMass:[...map.values()].reduce((n,p)=>n+p.fraction,0)});
  }
  if(!families.length)continue;
  progress('Checking GT group and proportion retrieval · target '+families[0].target);
  const all=solveFamilies(session,target,families);
  if(all.status==='retrieved'){checks.push(...all.checks);continue;}
  for(const f of families){const answer=families.length===1?all:solveFamilies(session,target,[f]);checks.push(...answer.checks);}
 }
 return{method:'joint-family-profile-retrieval',seconds:(performance.now()-began)/1000,checks};
}

function solveFamilies(session:MarginalSession,target:MarginalTargetState,families:Family[]){
 const began=performance.now(),p:MarginalProgram={...target.program,upper:target.program.upper.slice(),cost:target.program.cost.map(()=>0),rows:target.program.rows.slice()},atoms=new Map<Family,number[][]>();
 const row=(terms:LinearTerm[],lower=0,upper=lower)=>{const m=new Map<number,number>();for(const[i,v]of terms)m.set(i,(m.get(i)||0)+v);p.rows.push({terms:[...m].filter(([,v])=>Math.abs(v)>1e-15),lower,upper});};
 let impossible=false;
 for(const f of families){const g=f.group,ref=session.refs[g.reference];if(Math.abs(g.fraction/f.share-f.truthMass)>.02+1e-10)impossible=true;const byVariable=new Map<number,LinearTerm[]>(),owned=new Set<number>();
  for(const e of g.edges){owned.add(e.mass);for(const vars of e.free.values())vars.forEach(i=>owned.add(i));}for(const flow of g.flow?.values()||[]){flow.unchanged.forEach(c=>owned.add(c.column));flow.changed.forEach(i=>owned.add(i));}
  const profileAtoms:number[][]=[];
  for(const profile of f.profiles){const ids:number[]=[];
   for(const e of g.edges){const vector=new Map<number,number>([[e.mass,1]]);let valid=true;
    for(const[q,vars]of e.free){const b=DNA.indexOf(profile.sequence[q-e.incoming]);if(b<0){valid=false;break;}vector.set(vars[b],1);}
    if(!valid)continue;
    for(const q of g.bases.keys()){const b=DNA.indexOf(profile.sequence[q-e.incoming]),inserted=q>=e.start+e.incoming&&q<e.start+e.incoming+Math.max(0,e.size),at=q<e.start+e.incoming?q-e.incoming:q-e.outgoing,base=DNA.indexOf(ref.sequence[at]),flow=g.flow?.get(q);if(b<0){valid=false;break;}
     if(inserted||b===base){if(flow){const atom=inserted?e.free.get(q)![b]:e.mass,context=flow.unchanged.find(c=>c.base===b&&(inserted?c.cap.some(([i])=>i===atom):originContextContains(c.cap,e)));if(!context){valid=false;break;}vector.set(context.column,1);}}
     else{if(!flow){valid=false;break;}vector.set(flow.changed[b],1);}
    }
    if(!valid)continue;const id=p.upper.length;p.upper.push(f.share);p.cost.push(0);ids.push(id);for(const[i,v]of vector){const terms=byVariable.get(i)||[];terms.push([id,v]);byVariable.set(i,terms);}
   }
   if(!ids.length)impossible=true;profileAtoms.push(ids);const ratio=profile.fraction/f.truthMass;row([...ids.map(i=>[i,1]as LinearTerm),...g.mass.map(([i,v])=>[i,-v*ratio]as LinearTerm)]);row(ids.map(i=>[i,1]),Math.max(1e-8,(profile.fraction-.02)*f.share),Math.min(f.share,(profile.fraction+.02)*f.share));
  }
  atoms.set(f,profileAtoms);for(const i of owned)row([[i,1],...(byVariable.get(i)||[]).map(([j,v])=>[j,-v]as LinearTerm)]);
 }
 let x:number[]|undefined,status:'retrieved'|'not-retrieved'|'incomplete'='not-retrieved',reason:string|undefined;
 if(!impossible){let solver:MarginalLinearSolver|undefined;try{solver=new MarginalLinearSolver(p,session.optimizer);const answer=solver.feasibleMaximum([]);if(answer){x=answer.x;status='retrieved';}}catch(error){status='incomplete';reason=String(error);}finally{solver?.dispose();}}
 const seconds=(performance.now()-began)/1000;
 const checks=families.map(f=>({reference:f.reference,target:f.target,familyId:f.group.id,status,seconds,reason,profiles:f.profiles.map((profile,k)=>({sequenceKey:profile.sequenceKey,expectedFraction:profile.fraction,estimatedFraction:x?atoms.get(f)![k].reduce((n,i)=>n+x![i],0)/f.share:undefined}))}));
 return{status,checks};
}
