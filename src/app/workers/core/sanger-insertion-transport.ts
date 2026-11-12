import type {MarginalInternalGroup,MarginalTargetState} from './sanger-marginal-engine';
import type {LinearTerm} from './sanger-marginal-lp';
/** Restore latent origin/colour transports for a selected display query.
 * Basic inference retains their exact projection. For any nonnegative origin
 * masses and colour totals with equal sums, a transport exists. This explicit
 * lift is used when a query needs a particular origin's colour allocation. */
export function expandInsertionTransports(target:MarginalTargetState,group:MarginalInternalGroup){
 const p={...target.program,upper:target.program.upper.slice(),cost:target.program.cost.slice(),rows:target.program.rows.map(r=>({...r}))},point=target.point.slice(),edges=group.edges.map(e=>({...e,free:new Map<number,number[]>()})),byCycle=new Map<number,{flow:number[];edges:number[]}>();
 group.edges.forEach((e,i)=>{for(const[q,flow]of e.free){let f=byCycle.get(q);if(!f){f={flow,edges:[]};byCycle.set(q,f);}f.edges.push(i);}});
 for(const[q,f]of byCycle){const total=f.edges.reduce((n,i)=>n+Math.max(0,point[edges[i].mass]),0),columns:number[][]=[];
  for(const i of f.edges){const e=edges[i],mass=Math.max(0,point[e.mass]),vars=f.flow.map(v=>{const id=p.upper.length;p.upper.push(p.upper[e.mass]);p.cost.push(0);point.push(total?Math.max(0,point[v])*mass/total:0);return id;});e.free.set(q,vars);columns.push(vars);p.rows.push({terms:[...vars.map(i=>[i,1]as LinearTerm),[e.mass,-1]],lower:0,upper:0});}
  for(let b=0;b<4;b++)p.rows.push({terms:[...columns.map(v=>[v[b],1]as LinearTerm),[f.flow[b],-1]],lower:0,upper:0});
 }
 return{target:{...target,program:p,point},group:{...group,sharedInsertions:false,edges}};
}
