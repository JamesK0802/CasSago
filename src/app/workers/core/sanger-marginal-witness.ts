import type {MarginalInternalGroup} from './sanger-marginal-engine';
import type {MarginalProgram,LinearTerm} from './sanger-marginal-lp';
/** Prove coordinate endpoints with sparse full-model feasible witnesses.
 * A successful witness proves a 0/1 bound; failure simply retains the LP query.
 * No row is removed, no envelope is widened, and no genotype is enumerated. */
export function marginalEndpointWitnesses(program:MarginalProgram,point:number[]){
 const adjacency=program.upper.map(()=>[]as LinearTerm[]),values=program.rows.map((r,k)=>{for(const[i,w]of r.terms)adjacency[i].push([k,w]);return r.terms.reduce((n,[i,w])=>n+w*point[i],0);});
 return(g:MarginalInternalGroup,q:number)=>{const flow=g.flow?.get(q);if(!flow)return[false,false,false,false];const mass=g.mass.reduce((n,[i,w])=>n+w*point[i],0);if(mass<=1e-8)return[false,false,false,false];return Array.from({length:4},(_,base)=>{const changes=new Map<number,number>();let kept=0;for(const c of flow.unchanged){const v=c.base===base?c.cap.reduce((n,[i,w])=>n+w*point[i],0):0;kept+=v;changes.set(c.column,v);}flow.changed.forEach((col,b)=>changes.set(col,b===base?Math.max(0,mass-kept):0));const delta=new Map<number,number>();for(const[i,v]of changes){if(v< -1e-8||v>program.upper[i]+1e-8)return false;const d=v-point[i];if(!d)continue;for(const[k,w]of adjacency[i])delta.set(k,(delta.get(k)||0)+w*d);}for(const[k,d]of delta){const v=values[k]+d,r=program.rows[k];if(v<r.lower-1e-8||v>r.upper+1e-8)return false;}return true;});};
}
