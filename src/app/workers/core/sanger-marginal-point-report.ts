import type {SangerMarginalGroup,SangerMarginalDisplayComponent} from '../../models/sanger-analysis';
import type {MarginalModel} from './sanger-marginal-model';
import {fittedMarginalRepresentative} from './sanger-marginal-representative';
import {marginalInputEvent} from './sanger-marginal-state';
import type {LinearTerm} from './sanger-marginal-lp';
const DNA='ACGT';
/** No optimization or candidate enumeration. These are estimates from the
 * existing jointly fitted allocation, not bounds or phased DNA abundances.
 * [0,1] is a conservative internal placeholder, explicitly not reportable until
 * a selected-group range query finishes. It never confers GT acceptance.
 */
export function summarizeMarginalPoint(model:MarginalModel,point:number[],displayThreshold=2):SangerMarginalGroup[]{
 const evaluate=(terms:LinearTerm[])=>terms.reduce((n,[i,w])=>n+w*point[i],0);
 return model.drafts.flatMap(target=>{
  target.point=point;
  // Reporting metadata ONLY: structural base alternatives shared by displayed
  // families. This is not an uncertainty bound or a proof of nonidentifiability.
  // In particular, a fitted 100% must not hide a potentially shared allocation.
  // Low/hidden hypotheses are still present in the original fit and lazy domain.
  const potentials=new Map<string,Map<number,Set<number>>>();
  for(const g of target.groups){
   if(Math.max(0,evaluate(g.mass))/model.shares[g.reference]+1e-12<displayThreshold/100)continue;
   const ref=model.refs[g.reference],byPosition=new Map<number,Set<number>>();
   for(const q of g.bases.keys()){
    const bases=new Set<number>();for(const e of g.edges){if(model.program.upper[e.mass]<=1e-7)continue;if(e.free.has(q)){for(let b=0;b<4;b++)bases.add(b);}else{const at=q<e.start+e.incoming?q-e.incoming:q-e.outgoing,b=DNA.indexOf(ref.sequence[at]);if(b>=0)bases.add(b);}}
    byPosition.set(q,bases);
   }
   potentials.set(g.id,byPosition);
  }
  return target.groups.map(g=>{
   const ref=model.refs[g.reference],share=model.shares[g.reference],mass=Math.max(0,evaluate(g.mass));
   g.fraction=mass;g.fractionRange=[0,share];
   const group:SangerMarginalGroup={id:g.id,reference:ref.name,target:ref.reverse?ref.cuts.length-g.target:g.target+1,size:g.size,fraction:mass/share,fractionRange:[0,1],scope:'read-cycle-by-length',status:'estimated',boundMethod:'deferred-report-envelope',bases:[...g.bases].map(([q,forms])=>({position:ref.reverse?ref.input.sequence.length-ref.offset-q-1:ref.offset+q,estimate:(ref.reverse?forms.slice().reverse():forms).map(f=>mass?Math.max(0,evaluate(f))/mass:0),lower:[0,0,0,0],upper:[1,1,1,1]})),representative:fittedMarginalRepresentative(g,point,ref,share,(a,b,s)=>marginalInputEvent(ref,a,b,s,g.target))};
   const own=potentials.get(g.id);
   if(own)for(const base of group.bases){const q=ref.reverse?ref.input.sequence.length-ref.offset-base.position-1:base.position-ref.offset,possible=own.get(q)!;if(possible.size<2)continue;const peers=[...potentials].filter(([id,other])=>id!==g.id&&[...possible].filter(b=>other.get(q)?.has(b)).length>=2).map(([id])=>id);if(peers.length)base.sharedWith=peers;}
   const incoming=[...new Set(g.edges.map(e=>e.incoming))];
   if(g.size===0&&incoming.length===1){
    group.baselineReadOffset=ref.reverse?incoming[0]:-incoming[0];
    // Partition at one reportable SNP marker. Minor other departures must not
    // suppress a clear marker split, nor be silently declared phased: retain
    // their full conditional estimates on BOTH mixture rows. This conserves
    // every position marginal without claiming pure WT/SNP haplotypes.
    const baseline=(position:number)=>DNA.indexOf(ref.input.sequence[position+group.baselineReadOffset!]);
    const mixed=group.bases.filter(p=>{const b=baseline(p.position);return b>=0&&1-p.estimate[b]>1e-5;});
    // This fitted-marginal WT lower amount exists even when several major
    // SNP sites prevent a marker partition. Missing linkage must not turn
    // the engine's WT lower amount into an implicit zero via a missing field.
    group.unmodifiedEstimate=Math.max(0,group.fraction*(1-mixed.reduce((n,p)=>n+1-p.estimate[baseline(p.position)],0)));
    const reportable=mixed.filter(p=>group.fraction*(1-p.estimate[baseline(p.position)])+1e-12>=displayThreshold/100),site=reportable.length===1?reportable[0]:mixed.length===1?mixed[0]:undefined;
    if(site){
     const mixture=mixed.length>1;
     group.displayComponents=site.estimate.flatMap((v,b):SangerMarginalDisplayComponent[]=>{
      if(v<=1e-7)return[];
      const events=group.bases.flatMap(p=>{const base=p.position===site.position?b:p.estimate.indexOf(Math.max(...p.estimate)),position=p.position+group.baselineReadOffset!;return base===baseline(p.position)?[]:[{target:group.target,start:position,end:position+1,inserted:1,insertedSequence:DNA[base]}];}).sort((a,b)=>a.start-b.start);
      return[{rangeStatus:'estimated',id:g.id+'::base:'+b,fraction:group.fraction*v,fractionRange:[0,1],events,mixture,selection:{marker:{position:ref.reverse?ref.input.sequence.length-ref.offset-site.position-1:site.position-ref.offset,bases:[ref.reverse?3-b:b]}},bases:group.bases.map(p=>({...p,estimate:p.position===site.position?[0,1,2,3].map(j=>Number(j===b)):p.estimate}))}];
     });
    }
   }
   return group;
  });
 });
}
