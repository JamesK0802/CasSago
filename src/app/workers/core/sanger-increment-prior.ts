import type {LinearTerm} from './sanger-marginal-lp';
interface Edge {reference:number;target:number;incoming:number;size:number;mass:number;prefixColumns?:number[];transition?:number}
/** A soft representative prior, never a physical constraint: conditional
 * increment profiles prefer independence from the preceding phase. Aggregate
 * conditional differences are stabilized by measured noise resolution.
 * Aggregate positional aliases FIRST, so their number cannot influence the prior.
 * The returned auxiliary identities preserve every original allocation. */
export function incrementPrior(edges:Edge[],states:{size:number;column:number}[][][],seed:{column:number;value:number}[],shares:number[],firstAuxiliary:number,noise:number){
 if(!(Number.isFinite(noise)&&noise>0))throw Error('A positive measured noise scale is required.');
 const transitions=new Map<number,Edge>();for(const e of edges){const id=e.transition??e.prefixColumns?.at(-1)??e.mass,old=transitions.get(id);if(old&&(old.reference!==e.reference||old.target!==e.target||old.incoming!==e.incoming||old.size!==e.size))throw Error('Inconsistent aggregate transition identity.');transitions.set(id,e);}
 const cohorts=new Map<string,{reference:number;terms:LinearTerm[];entries:{id:number;edge:Edge}[]}>();for(const[id,edge]of transitions){const key=[edge.reference,edge.target,edge.size].join(':'),c=cohorts.get(key)||{reference:edge.reference,terms:[],entries:[]};c.terms.push([id,1]);c.entries.push({id,edge});cohorts.set(key,c);}
 const values=new Map(seed.map(s=>[s.column,s.value])),upper:number[]=[],rows:{terms:LinearTerm[];lower:number;upper:number}[]=[],forms:LinearTerm[][]=[];
 for(const c of cohorts.values()){const z=firstAuxiliary+upper.length;upper.push(shares[c.reference]);rows.push({terms:[[z,1],...c.terms.map(([i,w])=>[i,-w]as LinearTerm)],lower:0,upper:0});for(const{id,edge:e}of c.entries){let probability=1;if(e.target){const state=states[e.target-1][e.reference].find(s=>s.size===e.incoming);if(!state||!values.has(state.column))throw Error('Missing measured incoming-phase estimate.');probability=shares[e.reference]>0?values.get(state.column)!/shares[e.reference]:0;}const resolution=shares[e.reference]>0?noise/shares[e.reference]:noise,weight=shares[c.reference]>0?1/Math.sqrt(shares[c.reference]*(Math.max(0,probability)+resolution)):0;forms.push([[id,weight],[z,-probability*weight]]);}}
 return{upper,rows,forms,transitions:transitions.size};
}
