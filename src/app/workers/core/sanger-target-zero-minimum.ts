import {conditionalZeroWitnesses} from './sanger-conditional-zero-witness';
import type {MarginalInternalGroup} from './sanger-marginal-engine';
import type {MarginalProgram,LinearTerm} from './sanger-marginal-lp';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
type Product={u:number;cap:LinearTerm[];column:number};
/** Restricted target solves propose witnesses only. The existing verifier lifts
 * each candidate into EVERY original row, cone and box. A nonnegative mass
 * attained at zero has the original global lower bound zero. Failed proposals
 * leave the unrestricted bound query intact; no physical hypothesis is removed.
 */
export function targetZeroMinimumWitnesses(program:MarginalProgram,anchor:number[],groups:MarginalInternalGroup[],products:readonly Product[],optimizer:SangerLinearOptimizer){
 if(!groups.length)throw Error('A target witness needs at least one group.');
 const reference=groups[0].reference,target=groups[0].target;
 if(groups.some(g=>g.reference!==reference||g.target!==target))throw Error('A target witness must contain exactly one reference and target.');
 const ids=new Set<number>();
 for(const g of groups){
  for(const[i]of g.mass)ids.add(i);
  for(const forms of g.bases.values())for(const form of forms)for(const[i]of form)ids.add(i);
  for(const e of g.edges){ids.add(e.mass);for(const i of e.prefixColumns||[])ids.add(i);for(const vars of e.free.values())for(const i of vars)ids.add(i);}
  for(const flow of g.flow?.values()||[]){for(const i of flow.changed)ids.add(i);for(const c of flow.unchanged){ids.add(c.column);for(const[i]of c.cap)ids.add(i);}}
 }
 const donors=new Set(products.filter(p=>p.cap.some(([i,w])=>w!==0&&ids.has(i))).map(p=>p.column));
 for(const i of donors)ids.add(i);
 for(const p of products)if(donors.has(p.column)||p.cap.some(([i,w])=>w!==0&&ids.has(i)))ids.add(p.u);
 const witness=conditionalZeroWitnesses(program,anchor,groups[0],groups.map(g=>g.mass),optimizer,{additionalColumns:[...ids],seconds:60});
 return{...witness,possible(group:MarginalInternalGroup){const index=groups.indexOf(group);return index>=0&&witness.possible(index);}};
}
