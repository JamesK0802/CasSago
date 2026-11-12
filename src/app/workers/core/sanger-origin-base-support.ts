import {originContextContains} from './sanger-origin-prefix';
import type {MarginalInternalGroup,MarginalEdge} from './sanger-marginal-engine';
import type {LinearTerm} from './sanger-marginal-lp';
/** Exact nonnegative residual quantities required to extract a base from one
 * origin, including off-diagonal substitution Hall constraints. The shared
 * read-cycle base marginal is not an origin-specific allele certificate. */
export function originBaseQuantities(group:MarginalInternalGroup,e:MarginalEdge,q:number,b:number,sequence:string):LinearTerm[][]{
 const flow=group.flow?.get(q),inserted=q>=e.start+e.incoming&&q<e.start+e.incoming+Math.max(0,e.size);
 if(inserted){const free=e.free.get(q);if(!free)return[[]];if(!flow)return[[[free[b],1]]];const context=flow.unchanged.find(c=>c.base===b&&c.cap.some(([i,w])=>i===free[b]&&w>0));return context?[[[free[b],1]],[[context.column,1]]]:[[]];}
 const source='ACGT'.indexOf(sequence[q<e.start+e.incoming?q-e.incoming:q-e.outgoing]);if(!flow)return b===source?[]:[[]];const context=flow.unchanged.find(c=>c.base===source&&originContextContains(c.cap,e));if(!context)return[[]];if(b===source)return[[[context.column,1]]];
 const quantities:LinearTerm[][]=[[...context.cap,[context.column,-1]],[[flow.changed[b],1]]];for(let a=0;a<4;a++)if(a!==source&&a!==b)quantities.push([...flow.unchanged.filter(c=>c.base===a).flatMap(c=>[...c.cap.map(([i,w])=>[i,-w]as LinearTerm),[c.column,1]as LinearTerm]),...flow.changed.flatMap((i,c)=>c!==a?[[i,1]as LinearTerm]:[])]);return quantities;
}
export function originBaseMass(group:MarginalInternalGroup,e:MarginalEdge,q:number,b:number,sequence:string,x:number[]){return Math.min(x[e.mass],...originBaseQuantities(group,e,q,b,sequence).map(f=>f.reduce((n,[i,w])=>n+w*x[i],0)));}
