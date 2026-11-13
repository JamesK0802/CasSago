import{LP_INFINITY,nativeMarginalProgram,type MarginalProgram,type LinearTerm}from'./sanger-marginal-lp';
import type{SangerLinearOptimizer}from'./sanger-discrete-phase';

/** A lower envelope, NEVER an attained optimizer result. Any finite signed
 * row multipliers give a valid Lagrangian bound after every original column's
 * box correction. No primal feasibility, solver status or dual-gap claim is
 * made. This is only for the optional linear geometry relaxation used to
 * strengthen an outer bound; inference/optimizer acceptance remains unchanged.
 */
export function linearDualOuterLower(p:MarginalProgram,cost:ArrayLike<number>,dual:ArrayLike<number>):number|undefined{
 if(p.quadratics?.length||p.quadraticCost?.length||cost.length!==p.upper.length||dual.length!==p.rows.length)return;
 const reduced=Float64Array.from(cost),compensation=new Float64Array(cost.length);
 let lower=0,correction=0,magnitude=0,magnitudeCorrection=0;
 const addMagnitude=(v:number)=>{const a=v-magnitudeCorrection,s=magnitude+a;magnitudeCorrection=(s-magnitude)-a;magnitude=s;};
 const accumulate=(value:number)=>{const adjusted=value-correction,sum=lower+adjusted;correction=(sum-lower)-adjusted;lower=sum;};
 for(let i=0;i<cost.length;i++){const u=p.upper[i],c=cost[i];if(!Number.isFinite(c)||!Number.isFinite(u)||u<0||u>=LP_INFINITY)return;addMagnitude(Math.abs(c)*u);}
 for(let k=0;k<p.rows.length;k++){
  const r=p.rows[k];if(!Number.isFinite(r.lower)||!Number.isFinite(r.upper)||r.lower>r.upper)return;for(const[i,w]of r.terms)if(!Number.isSafeInteger(i)||i<0||i>=cost.length||!Number.isFinite(w))return;let y=dual[k];if(!Number.isFinite(y))return;
  if(y>0&&r.lower<=-LP_INFINITY||y<0&&r.upper>=LP_INFINITY)y=0;
  if(!y)continue;const constant=y*(y>0?r.lower:r.upper);accumulate(constant);addMagnitude(Math.abs(constant));
  for(const[i,w]of r.terms){if(!Number.isSafeInteger(i)||i<0||i>=cost.length||!Number.isFinite(w))return;const product=-y*w,adjusted=product-compensation[i],sum=reduced[i]+adjusted;compensation[i]=(sum-reduced[i])-adjusted;reduced[i]=sum;addMagnitude(Math.abs(product)*p.upper[i]);}
 }
 for(let i=0;i<cost.length;i++){const value=Math.min(0,reduced[i])*p.upper[i];accumulate(value);addMagnitude(Math.abs(value));}
 // Include input products and all boxed coefficient magnitudes in the guard.
 // Compensated accumulation retains small signed coefficients, without
 // rounding objective keys or dropping physical terms.
 const bound=lower-256*Number.EPSILON*(1+magnitude);return Number.isFinite(bound)?bound:undefined;
}

/** Short, warm LP proposals; only the independently derived original-row/box
 * lower envelope escapes. Failed/incomplete proposals never become fits.
 */
export function cachedLinearDualOuter(p:MarginalProgram,optimizer:SangerLinearOptimizer,seconds=1){
 let native:any,previous:Float64Array|undefined,domain=p.rows.map(r=>({...r})),upper=p.upper.slice();
 const cache=new Map<string,number>(),stats={hits:0,solves:0,partial:0,failures:0,seconds:0,lastFailure:undefined as string|undefined};
 const unchanged=()=>upper.length===p.upper.length&&upper.every((u,i)=>u===p.upper[i])&&domain.length===p.rows.length&&domain.every((r,k)=>r.terms===p.rows[k].terms&&r.lower===p.rows[k].lower&&r.upper===p.rows[k].upper);
 const release=()=>{const old=native;native=undefined;previous=undefined;try{old?.dispose();}catch(error){stats.lastFailure=String(error);}};
 return{stats,lower(terms:LinearTerm[],maximum=false){
  if(!p.immutableRows||!unchanged()){release();cache.clear();domain=p.rows.map(r=>({...r}));upper=p.upper.slice();}
  const cost=new Float64Array(p.upper.length);let objectiveError=0,errorCompensation=0;
  for(const[i,w]of terms){
   if(!Number.isSafeInteger(i)||i<0||i>=cost.length||!Number.isFinite(w))throw Error('Invalid outer objective.');
   const value=(maximum?-1:1)*w,old=cost[i],sum=old+value;if(!Number.isFinite(sum))throw Error('Invalid coalesced outer objective.');
   // TwoSum exposes the exact rounding residual of this coefficient addition.
   // Correct the returned lower bound even if duplicate terms cancel after
   // rounding, rather than certifying a falsely zero original objective.
   const part=sum-old,error=(old-(sum-part))+(value-part);cost[i]=sum;
   if(error){const u=p.upper[i];if(!Number.isFinite(u)||u<0||u>=LP_INFINITY)throw Error('Invalid objective box.');let weighted=Math.abs(error)*u;if(!Number.isFinite(weighted))throw Error('Invalid objective roundoff.');if(!weighted&&u>0)weighted=Number.MIN_VALUE;const adjusted=weighted-errorCompensation,next=objectiveError+adjusted;errorCompensation=(next-objectiveError)-adjusted;objectiveError=next;}
  }
  if(objectiveError)objectiveError+=64*Number.EPSILON*(1+objectiveError);
  if(!Number.isFinite(objectiveError))throw Error('Invalid objective roundoff.');
  const outward=(lower:number)=>objectiveError?lower-objectiveError-64*Number.EPSILON*(1+Math.abs(lower)+objectiveError):lower;
  if(p.quadratics?.length||p.quadraticCost?.length)throw Error('Linear outer envelope requires a linear relaxation.');
  const key=JSON.stringify(Array.from(cost).flatMap((w,i)=>w?[[i,w]]:[])),cached=cache.get(key);if(cached!==undefined){stats.hits++;return outward(cached);}
  // The identically zero objective has the universal lower bound zero.
  // Solving a large feasibility LP cannot strengthen this certificate.
  if(cost.every(c=>c===0))return outward(0);
  // The zero-multiplier certificate is always available for finite boxes.
  let lower=linearDualOuterLower(p,cost,new Float64Array(p.rows.length));if(lower===undefined)throw Error('Invalid finite-box outer program.');
  const start=performance.now();stats.solves++;
  try{
   if(!native){native=optimizer.createModel!();native.options.set({output_flag:false,solver:'simplex',time_limit:seconds,primal_feasibility_tolerance:1e-9,dual_feasibility_tolerance:1e-9});native.passModel(nativeMarginalProgram({...p,cost:Array.from(cost)}));}
   else{const changed=Array.from(cost).flatMap((c,i)=>c!==previous![i]?[i]:[]);if(changed.length)native.changeColsCost({kind:'set',indices:changed},changed.map(i=>cost[i]));}
   previous=cost;native.zeroAllClocks();const status=native.run().modelStatus;if(status!==7)stats.partial++;
   const bound=linearDualOuterLower(p,cost,native.getSolution().rowDual as ArrayLike<number>);if(bound!==undefined)lower=Math.max(lower,bound);
  }catch(error){stats.failures++;stats.lastFailure=String(error);release();}
  finally{stats.seconds+=(performance.now()-start)/1000;}
  cache.set(key,lower);if(cache.size>16)cache.delete(cache.keys().next().value!);return outward(lower);
 },dispose(){release();cache.clear();}};
}
