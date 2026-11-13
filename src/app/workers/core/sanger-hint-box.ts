import type {LinearTerm} from './sanger-marginal-lp';
/** Exact partial minimization over omitted boxes for repeated dual hints.
 * Bounds and the hint belong to one immutable bound-certificate workspace.
 * Sparse query changes need not copy the entire global coefficient map. */
export function prepareHintBox(image:Map<number,number>,upper:readonly number[],geometric:ReadonlyMap<number,number>){
 const geometry:LinearTerm[]=[];let box=0,compensation=0,magnitude=0;for(const[i,w]of image){if(geometric.has(i)){geometry.push([i,w]);continue;}const value=Math.min(0,-w)*upper[i],adjusted=value-compensation,sum=box+adjusted;compensation=(sum-box)-adjusted;box=sum;magnitude+=Math.abs(value);}
 return(reduced:Map<number,number>,alpha:number)=>{if(!Number.isFinite(alpha)||alpha<0)throw Error('A nonnegative finite dual-hint multiplier is required.');const coefficients=new Map<number,number>();let constant=alpha*box,scale=alpha*magnitude;for(const[i,w]of reduced){if(geometric.has(i)){coefficients.set(i,w);continue;}const before=alpha*Math.min(0,-(image.get(i)||0))*upper[i],after=Math.min(0,w-alpha*(image.get(i)||0))*upper[i];constant+=after-before;scale+=Math.abs(before)+Math.abs(after);}if(alpha)for(const[i,w]of geometry)coefficients.set(i,(coefficients.get(i)||0)-alpha*w);return{coefficients,constant,magnitude:scale};};
}
