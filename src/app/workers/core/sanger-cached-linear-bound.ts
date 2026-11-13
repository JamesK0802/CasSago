import {marginalViolation,type MarginalProgram,type LinearTerm} from './sanger-marginal-lp';
interface Solver {solve(terms:LinearTerm[],maximum?:boolean):{x:number[];value:number};objectiveLowerBound():number|undefined}
/** For a fixed polytope, L(c0)+min_box((c-c0)x) is a valid global lower
 * bound for c. Reuse only exact objectives or an attained primal/dual gap at
 * the ORIGINAL 1e-7 tolerance. This is a proof transfer, not rounded matching. */
export function cachedLinearBound(p:MarginalProgram,solver:Solver){
 let domain=p.rows.map(r=>({terms:r.terms,lower:r.lower,upper:r.upper})),upper=p.upper.slice();const cache:{cost:Float64Array;point:number[];lower:number;pivot:number}[]=[],stats={hits:0,solves:0,transfers:0};
 const unchanged=()=>upper.length===p.upper.length&&upper.every((v,i)=>v===p.upper[i])&&domain.length===p.rows.length&&domain.every((r,k)=>r.terms===p.rows[k].terms&&r.lower===p.rows[k].lower&&r.upper===p.rows[k].upper);
 return{stats,lower(terms:LinearTerm[],maximum=false){
  if(!unchanged()){cache.length=0;domain=p.rows.map(r=>({terms:r.terms,lower:r.lower,upper:r.upper}));upper=p.upper.slice();}
  const cost=new Float64Array(upper.length);for(const[i,w]of terms){if(!Number.isSafeInteger(i)||i<0||i>=upper.length||!Number.isFinite(w))throw Error('Invalid cached objective coefficient.');cost[i]+=(maximum?-1:1)*w;}
  if(p.immutableRows&&!p.quadratics?.length&&!p.quadraticCost?.length)for(const entry of cache){const ratio=entry.cost[entry.pivot]?cost[entry.pivot]/entry.cost[entry.pivot]:0;
   for(const scale of ratio>0&&Number.isFinite(ratio)&&ratio!==1?[1,ratio]:[1]){let penalty=0,compensation=0,magnitude=Math.abs(scale*entry.lower),value=0,identical=true;
    for(let i=0;i<cost.length;i++){const difference=cost[i]-scale*entry.cost[i];if(difference)identical=false;const v=Math.min(0,difference)*upper[i],adjusted=v-compensation,sum=penalty+adjusted;compensation=(sum-penalty)-adjusted;penalty=sum;magnitude+=Math.abs(v);value+=cost[i]*entry.point[i];}
    if(identical&&scale===1){stats.hits++;return entry.lower;}
    const bound=scale*entry.lower+penalty-128*Number.EPSILON*(1+magnitude),gap=value-bound;if(Number.isFinite(bound)&&gap>=-1e-7&&gap<=1e-7){stats.hits++;stats.transfers++;return bound;}
   }
  }
  stats.solves++;const answer=solver.solve(terms,maximum),lower=solver.objectiveLowerBound();
  if(lower===undefined||!Number.isFinite(lower))throw Error('Geometry bound requires an independently certified objective lower bound.');
  if(p.immutableRows&&!p.quadratics?.length&&!p.quadraticCost?.length&&Number.isFinite(lower)&&marginalViolation(p,answer.x)<=1e-7){let pivot=0;for(let i=1;i<cost.length;i++)if(Math.abs(cost[i])>Math.abs(cost[pivot]))pivot=i;cache.unshift({cost,point:answer.x.slice(),lower,pivot});if(cache.length>8)cache.pop();}return lower;
 }};
}
