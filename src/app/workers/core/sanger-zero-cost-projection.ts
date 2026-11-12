import type {MarginalProgram,LinearTerm} from './sanger-marginal-lp';
/** Exact substitution on the nonnegative objective's zero-cost face. The
 * supplied identities are implied by the original balance/capacity equations
 * when positive-cost columns are zero. Every lifted answer is still checked
 * against the original program by the caller. No candidate is discarded. */
export function zeroCostProjection(p:MarginalProgram){
 if(!p.zeroCostSubstitutions?.length||p.cost.some(c=>c<0))return;
 const substitution=new Map(p.zeroCostSubstitutions.map(s=>[s.column,s.terms])),retained=p.upper.flatMap((u,i)=>u>0&&p.cost[i]===0&&!substitution.has(i)?[i]:[]),indices=new Map(retained.map((i,j)=>[i,j])),cache=new Map<number,LinearTerm[]>(),visiting=new Set<number>();
 const form=(terms:LinearTerm[]):LinearTerm[]=>{const result=new Map<number,number>();for(const[i,w]of terms)for(const[j,v]of expand(i))result.set(j,(result.get(j)||0)+w*v);return[...result].filter(([,v])=>v!==0);};
 const expand=(i:number):LinearTerm[]=>{const found=cache.get(i);if(found)return found;if(visiting.has(i))throw Error('Cyclic zero-cost substitution.');visiting.add(i);const next=p.cost[i]>0||p.upper[i]===0?[]:substitution.has(i)?form(substitution.get(i)!):[[indices.get(i)!,1]as LinearTerm];visiting.delete(i);cache.set(i,next);return next;};
 // Canonicalize one expanded row at a time. Keeping a second complete
 // expanded-row array alongside the canonical model doubles the peak heap.
 const unique=new Map<string,MarginalProgram['rows'][number]>();
 const add=(input:LinearTerm[],lo:number,hi:number)=>{if(!input.length&&lo<=0&&hi>=0)return;const sorted=input.slice().sort((a,b)=>a[0]-b[0]),scale=sorted[0]?.[1]||1,terms=sorted.map(([i,w])=>[i,w/scale]as LinearTerm),lower=(scale>0?lo:hi)/scale,upper=(scale>0?hi:lo)/scale,key=terms.map(([i,w])=>i+':'+w).join(','),previous=unique.get(key);if(previous){previous.lower=Math.max(previous.lower,lower);previous.upper=Math.min(previous.upper,upper);}else unique.set(key,{terms,lower,upper});};
 for(const row of p.rows)add(form(row.terms),row.lower,row.upper);
 for(const[column]of substitution)add(expand(column),0,p.upper[column]);
 const program:MarginalProgram={upper:retained.map(i=>p.upper[i]),cost:retained.map(()=>0),rows:[...unique.values()],quadratics:p.quadratics?.map(q=>({...q,terms:q.terms.map(form)}))};
 return{program,indices,lift:(x:number[])=>p.upper.map((_,i)=>expand(i).reduce((s,[j,w])=>s+w*x[j],0))};
}
