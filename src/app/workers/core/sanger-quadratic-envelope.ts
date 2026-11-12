import type {LinearTerm} from './sanger-marginal-lp';
export interface QuadraticEnvelope {name:string;terms:LinearTerm[][];upper:number}
export function evaluateQuadratic(envelope:QuadraticEnvelope,x:number[]){let value=0;for(const row of envelope.terms){const residual=row.reduce((n,[i,w])=>n+w*x[i],0);value+=residual*residual;}return value;}
/** A supporting plane at the radial boundary of a convex sum of squared
 * raw-signal residuals. Cauchy–Schwarz gives g(x)·y <= 2 sqrt(Q(x) U)
 * for every Q(y)<=U; this is tighter than a tangent at an exterior point.
 * Every feasible point satisfies this plane. Returning an LP optimum only
 * after checking the original quadratic gives a global bound certificate. */
export function quadraticCut(envelope:QuadraticEnvelope,x:number[]){const gradient=new Map<number,number>();let value=0;for(const row of envelope.terms){const r=row.reduce((n,[i,w])=>n+w*x[i],0);value+=r*r;for(const[i,w]of row)gradient.set(i,(gradient.get(i)||0)+2*r*w);}return{value,terms:[...gradient].filter(([,w])=>w!==0)as LinearTerm[],upper:2*Math.sqrt(Math.max(0,envelope.upper*value))};}
/** Necessary tangent planes, seeded from pure length states. This is only a
 * lossless numerical outer approximation; the original squared-residual
 * envelope is still verified before any result is accepted. */
export function seedQuadraticTangents(envelope:QuadraticEnvelope,seed:number[],groups:{columns:number[];mass:number}[]){
 const initial=envelope.terms.map(row=>row.reduce((n,[i,w])=>n+w*seed[i],0)),baseline=initial.reduce((n,v)=>n+v*v,0),point=seed.slice(),cuts:{terms:LinearTerm[];upper:number}[]=[];
 for(const group of groups){const members=new Set(group.columns);for(const chosen of group.columns){let a=0,b=0;for(let k=0;k<envelope.terms.length;k++){const delta=envelope.terms[k].reduce((n,[i,w])=>n+(members.has(i)?w*((i===chosen?group.mass:0)-seed[i]):0),0);a+=delta*delta;b+=initial[k]*delta;}if(baseline+2*b+a<=envelope.upper||a===0)continue;const extra=Math.max(0,envelope.upper-baseline),root=Math.sqrt(b*b+a*extra),alpha=b>=0?extra/(b+root):(-b+root)/a;for(const i of group.columns)point[i]=seed[i]+alpha*((i===chosen?group.mass:0)-seed[i]);const cut=quadraticCut(envelope,point);cuts.push({terms:cut.terms,upper:cut.upper+1e-10*(1+envelope.upper)});for(const i of group.columns)point[i]=seed[i];}}
 return cuts;
}
