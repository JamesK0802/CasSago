import type {LinearTerm} from './sanger-marginal-lp';
/** Immutable sparse row storage. Keep every index and IEEE-754 coefficient;
 * avoid retaining a separate JS array object for every coefficient pair.
 * Array-producing operations return ordinary mutable arrays, as before. */
class PackedTerms {
 private indices:Uint32Array;
 private weights:Float64Array;
 readonly length:number;
 constructor(terms:LinearTerm[]){this.length=terms.length;this.indices=new Uint32Array(this.length);this.weights=new Float64Array(this.length);terms.forEach(([i,w],k)=>{if(!Number.isSafeInteger(i)||i<0||i>0xffffffff)throw Error('Invalid sparse column index.');this.indices[k]=i;this.weights[k]=w;});}
 private term(k:number):LinearTerm{return[this.indices[k],this.weights[k]];}
 dot(x:number[]){let sum=0;for(let k=0;k<this.length;k++)sum+=this.weights[k]*x[this.indices[k]];return sum;}
 *[Symbol.iterator](){for(let k=0;k<this.length;k++)yield this.term(k);}
 map<T>(fn:(v:LinearTerm,k:number,a:LinearTerm[])=>T):T[]{const out:T[]=[];for(let k=0;k<this.length;k++)out.push(fn(this.term(k),k,this as unknown as LinearTerm[]));return out;}
 flatMap<T>(fn:(v:LinearTerm,k:number,a:LinearTerm[])=>T|T[]):T[]{const out:T[]=[];for(let k=0;k<this.length;k++){const value=fn(this.term(k),k,this as unknown as LinearTerm[]);if(Array.isArray(value))out.push(...value);else out.push(value);}return out;}
 reduce<T>(fn:(v:T,term:LinearTerm,k:number)=>T,initial:T){let value=initial;for(let k=0;k<this.length;k++)value=fn(value,this.term(k),k);return value;}
 forEach(fn:(v:LinearTerm,k:number)=>void){for(let k=0;k<this.length;k++)fn(this.term(k),k);}
 every(fn:(v:LinearTerm,k:number)=>boolean){for(let k=0;k<this.length;k++)if(!fn(this.term(k),k))return false;return true;}
 some(fn:(v:LinearTerm,k:number)=>boolean){for(let k=0;k<this.length;k++)if(fn(this.term(k),k))return true;return false;}
 filter(fn:(v:LinearTerm,k:number)=>boolean){const out:LinearTerm[]=[];for(let k=0;k<this.length;k++){const value=this.term(k);if(fn(value,k))out.push(value);}return out;}
 slice(start=0,end=this.length){return this.map(v=>v).slice(start,end);}
 at(k:number){if(k<0)k+=this.length;return k>=0&&k<this.length?this.term(k):undefined;}
 toJSON(){return this.map(v=>v);}
}
export function packLinearTerms(terms:LinearTerm[],minimumLength=64):LinearTerm[]{return terms.length<minimumLength?terms:new PackedTerms(terms)as unknown as LinearTerm[];}
/** Keep original accumulation order while avoiding a callback per coefficient. */
export function linearTermValue(terms:LinearTerm[],x:number[]){if(terms instanceof PackedTerms)return terms.dot(x);let sum=0;for(let k=0;k<terms.length;k++){const term=terms[k];sum+=term[1]*x[term[0]];}return sum;}
