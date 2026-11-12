import {LP_INFINITY,type MarginalProgram,type LinearTerm} from './sanger-marginal-lp';
/** Convex outer transport of independently calibrated WT photon uncertainty.
 * WT photons are conserved exactly at every shared-colour donor position.
 * McCormick envelopes retain every physical product of allele mass and unknown
 * donor response. This uncertainty carries no biological substitution cost. */
export function emissionUncertainty(program:MarginalProgram,shares:number[],sequences:string[],responses:number[][][],bounds:[number,number][][],blocks:Map<number,Set<number>>){
 const groups=new Map<string,Map<number,{column:number;lo:number;range:number}>>(),products=new Map<string,LinearTerm[]>(),seeds:{u:number;cap:LinearTerm[];column:number}[]=[],variable=(upper:number)=>{const i=program.upper.length;program.upper.push(upper);program.cost.push(0);return i;};
 const donor=(h:number,q:number)=>{
  const base=sequences[h][q],key=q+':'+base;if(groups.has(key))return groups.get(key)!.get(h);
  const columns=new Map<number,{column:number;lo:number;range:number}>(),terms:LinearTerm[]=[];let constant=0;
  for(let j=0;j<shares.length;j++)if(sequences[j][q]===base){const[lo,hi]=bounds[j]?.[q]||[0,0],range=hi-lo;if(!(range>1e-12))continue;const column=variable(range),mass=shares[j]*(responses[j][q]||[]).reduce((a,b)=>a+b,0);columns.set(j,{column,lo,range});terms.push([column,mass]);constant-=mass*lo;}
  if(terms.length)program.rows.push({terms,lower:constant,upper:constant});groups.set(key,columns);return columns.get(h);
 };
 const correction=(cap:LinearTerm[],h:number,q:number,p:number):LinearTerm[]=>{
  const d=donor(h,q);if(!d||!cap.length)return[];const key=h+':'+q+':'+JSON.stringify(cap),known=products.get(key);if(known)return known;const share=shares[h],u=variable(share*d.range),rows=program.rows;
  rows.push({terms:[[u,1],...cap.map(([i,w])=>[i,-w*d.range]as LinearTerm)],lower:-LP_INFINITY,upper:0},{terms:[[u,1],[d.column,-share]],lower:-LP_INFINITY,upper:0},{terms:[[u,1],...cap.map(([i,w])=>[i,-w*d.range]as LinearTerm),[d.column,-share]],lower:-share*d.range,upper:LP_INFINITY});
  const set=blocks.get(p)||new Set<number>();set.add(u);blocks.set(p,set);const result:LinearTerm[]=[...cap.map(([i,w])=>[i,w*d.lo]as LinearTerm),[u,1]];products.set(key,result);seeds.push({u,cap,column:d.column});return result;
 };
 return Object.assign(correction,{seed:(x:number[])=>{for(const group of groups.values())for(const d of group.values())x[d.column]=-d.lo;for(const s of seeds)x[s.u]=s.cap.reduce((v,[i,w])=>v+w*x[i],0)*x[s.column];}});
}
