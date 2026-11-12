import type {MarginalInternalGroup,MarginalEdge} from './sanger-marginal-engine';
import type {LinearTerm,MarginalProgram} from './sanger-marginal-lp';
/** Exact sparse factorization of repeated origin sums. Original origin masses
 * remain explicit. Prefix recurrences replace repeated dense incidence rows;
 * every original origin allocation has one unique lifted prefix allocation. */
export function originPrefixes(group:MarginalInternalGroup,p:MarginalProgram,variable:(upper:number)=>number,row:(terms:LinearTerm[],lower:number,upper:number)=>unknown){
 const incoming=[...new Set(group.edges.map(e=>e.incoming))];
 const buckets=incoming.map(a=>{const edges=group.edges.filter(e=>e.incoming===a).sort((x,y)=>x.start-y.start),columns:number[]=[];for(let j=0;j<edges.length;j++){const column=j?variable(p.upper[edges[j].mass]):edges[j].mass;if(j)row([[column,1],[columns[j-1],-1],[edges[j].mass,-1]],0,0);columns.push(column);}edges.forEach((e,j)=>e.prefixColumns=columns.slice(j));return{incoming:a,outgoing:edges[0].outgoing,edges,total:columns.at(-1)!,at(threshold:number):LinearTerm[]{let lo=0,hi=edges.length;while(lo<hi){const m=(lo+hi)>>1;if(edges[m].start<=threshold)lo=m+1;else hi=m;}return lo?[[columns[lo-1],1]]:[];}};});
 group.mass=buckets.map(b=>[b.total,1]);return buckets;
}
export function originContextContains(terms:LinearTerm[],edge:MarginalEdge){return terms.reduce((n,[i,w])=>n+(i===edge.mass||edge.prefixColumns?.includes(i)?w:0),0)>.5;}
