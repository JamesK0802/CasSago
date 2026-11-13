import {marginalViolation,LP_INFINITY,OptimizationPending,type MarginalProgram,type LinearTerm} from './sanger-marginal-lp';
export interface SangerConicRuntime {solve:(pointers:Uint32Array,indices:Uint32Array,values:Float64Array,cost:Float64Array,aPointers:Uint32Array,aIndices:Uint32Array,aValues:Float64Array,b:Float64Array,n:number,m:number,cones:string,settings:string)=>any}
/** The same feasible set, represented directly by nonnegative and Euclidean
 * cones. Independently validate the original primal and a repaired conic dual
 * certificate, including all original variable boxes. Solver status alone is
 * never sufficient. */
export class VerifiedConicSolver {
 private timeLimit=300;
 private iterationLimit=300;
 private tolerance=1e-8;
 private certifiedLower?:number;
 private retryableCertificate=false;
 objectiveLowerBound(){return this.certifiedLower;}
 setTimeLimit(seconds:number){this.timeLimit=seconds;}
 private compiled?:ReturnType<VerifiedConicSolver['compile']>;
 private dual:number[]=[];
 private planes:{terms:LinearTerm[];upper:number}[]=[];
 readonly stats={queries:0,seconds:0,iterations:0,maximumViolation:0,maximumGap:0,certificationRetries:0,maximumRejectedGap:0,scalingRetries:0,statusIndependentCertificates:0,boxCertificates:0};
 constructor(private p:MarginalProgram,private runtime:SangerConicRuntime){}
 invalidate(){this.compiled=undefined;}
 private compile(){
  const p=this.p,columns=Int32Array.from(p.upper.flatMap((u,i)=>u===0?[]:[i])),index=new Int32Array(p.upper.length).fill(-1);columns.forEach((i,j)=>index[i]=j);const n=columns.length,upper=Float64Array.from(columns,i=>p.upper[i]),equalities:number[]=[],linear:number[]=[],sizes=(p.quadratics||[]).map(q=>q.terms.length+1),counts=new Uint32Array(n);let finiteUpper=0;
  // Keep row indices and signs, rather than duplicating millions of sparse
  // coefficient pairs for negative inequalities, boxes and cone residuals.
  p.rows.forEach((r,k)=>{const sides=r.lower===r.upper?1:Number(r.lower>-LP_INFINITY)+Number(r.upper<LP_INFINITY);if(r.lower===r.upper)equalities.push(k);else{if(r.lower>-LP_INFINITY)linear.push(2*k+1);if(r.upper<LP_INFINITY)linear.push(2*k);}for(const[i,w]of r.terms)if(w&&index[i]>=0)counts[index[i]]+=sides;});
  for(let i=0;i<n;i++){counts[i]++;if(upper[i]<LP_INFINITY){counts[i]++;finiteUpper++;}}
  for(const q of p.quadratics||[])for(const terms of q.terms)for(const[i,w]of terms)if(w&&index[i]>=0)counts[index[i]]++;
  const zero=equalities.length,nonneg=linear.length+n+finiteUpper,rowCount=zero+nonneg+sizes.reduce((a,b)=>a+b,0),starts=new Uint32Array(n+1);for(let i=0;i<n;i++)starts[i+1]=starts[i]+counts[i];
  const indices=new Uint32Array(starts[n]),values=new Float64Array(starts[n]),cursor=starts.slice(),b=new Float64Array(rowCount),originalRows=new Int32Array(rowCount).fill(-1),shadowSigns=new Int8Array(rowCount),lowerBoxes=new Int32Array(n),upperBoxes=new Int32Array(n).fill(-1);let k=0;
  const put=(terms:LinearTerm[],scale:number,rhs:number,original:number,shadow:number)=>{b[k]=rhs;originalRows[k]=original;shadowSigns[k]=shadow;for(const[i,w]of terms)if(w&&index[i]>=0){const j=cursor[index[i]]++;indices[j]=k;values[j]=scale*w;}k++;};
  for(const row of equalities)put(p.rows[row].terms,1,p.rows[row].lower,row,-1);
  for(const code of linear){const row=Math.floor(code/2),lower=code%2===1,r=p.rows[row];put(r.terms,lower?-1:1,lower?-r.lower:r.upper,row,lower?1:-1);}
  for(let i=0;i<n;i++){lowerBoxes[i]=k;let j=cursor[i]++;indices[j]=k;values[j]=-1;k++;if(upper[i]<LP_INFINITY){upperBoxes[i]=k;j=cursor[i]++;indices[j]=k;values[j]=1;b[k]=upper[i];k++;}}
  for(const q of p.quadratics||[]){b[k++]=Math.sqrt(q.upper);for(const terms of q.terms)put(terms,-1,0,-1,0);}
  const hessian=new Map<number,Map<number,number>>();for(const form of p.quadraticCost||[]){const merged=new Map<number,number>();for(const[i,w]of form)if(index[i]>=0)merged.set(index[i],(merged.get(index[i])||0)+w);const terms=[...merged].filter(([,w])=>w!==0).sort((a,b)=>a[0]-b[0]);for(let a=0;a<terms.length;a++)for(let b=a;b<terms.length;b++){const[i,v]=terms[a],[j,w]=terms[b],col=hessian.get(j)||new Map<number,number>();col.set(i,(col.get(i)||0)+2*v*w);hessian.set(j,col);}}const pStarts=new Uint32Array(n+1);for(let i=0;i<n;i++)pStarts[i+1]=pStarts[i]+[...(hessian.get(i)?.values()||[])].filter(w=>w!==0).length;const pIndices=new Uint32Array(pStarts[n]),pValues=new Float64Array(pStarts[n]);for(let j=0;j<n;j++){let k=pStarts[j];for(const[i,w]of [...(hessian.get(j)||[])].sort((a,b)=>a[0]-b[0]))if(w!==0){pIndices[k]=i;pValues[k++]=w;}}
  return{pStarts,pIndices,pValues,n,columns,index,upper,rowCount,originalRows,shadowSigns,lowerBoxes,upperBoxes,starts,indices,values,b,zero,nonneg,sizes,cones:JSON.stringify({zero,nonneg,soc:sizes,exp:0,power:[]})};
 }
 solve(terms?:LinearTerm[],maximize=false):{x:number[];value:number}{
  this.certifiedLower=undefined;if(maximize&&this.p.quadraticCost?.length)throw Error('Convex quadratic objectives require minimization.');
  for(;;){try{return this.attempt(terms,maximize);}catch(error){if(this.tolerance<=1e-12||!this.retryableCertificate)throw error;this.tolerance=Math.max(1e-12,this.tolerance/100);this.stats.certificationRetries++;}}
 }
 private attempt(terms?:LinearTerm[],maximize=false){
  this.retryableCertificate=false;
  if(this.timeLimit===0)throw new OptimizationPending();
  const start=performance.now(),m=this.compiled??=this.compile(),cost=new Float64Array(m.n);if(terms){for(const[i,w]of terms)if(m.index[i]>=0)cost[m.index[i]]+=(maximize?-1:1)*w;}else m.columns.forEach((i,j)=>cost[j]=this.p.cost[i]);if(!m.n){const x=this.p.upper.map(()=>0),violation=marginalViolation(this.p,x);if(violation>1e-7)throw Error('Zero-box model violates original constraints.');this.certifiedLower=0;this.dual=this.p.rows.map(()=>0);this.planes=[];return{x,value:0};}
  let result:any,lastTimeLimit=this.timeLimit;
  // Positive objective scaling leaves the entire feasible set and optimizer
  // unchanged. A numerical failure may be sensitive to that scale. Convert
  // dual multipliers back BEFORE checking the original certificate below.
  // Retries share the existing time slice; they never relax a tolerance.
  for(const scale of [1,10,.1]){
   if(scale!==1&&!(result?.status==='unknown'&&!(result.z?.length)))break;
   const remaining=scale===1?this.timeLimit:this.timeLimit>0?this.timeLimit-(performance.now()-start)/1000:this.timeLimit;if(this.timeLimit>0&&remaining<=0)throw new OptimizationPending();
   lastTimeLimit=remaining;if(scale!==1)this.stats.scalingRetries++;
   const answer=this.runtime.solve(m.pStarts,m.pIndices,scale===1?m.pValues:Float64Array.from(m.pValues,v=>v*scale),scale===1?cost:Float64Array.from(cost,v=>v*scale),m.starts,m.indices,m.values,m.b,m.n,m.rowCount,m.cones,JSON.stringify({verbose:false,max_iter:this.iterationLimit,time_limit:remaining,tol_gap_abs:this.tolerance,tol_gap_rel:this.tolerance}));
   result=scale===1||!answer.z?answer:{...answer,z:Array.from(answer.z as ArrayLike<number>,v=>v/scale)};
  }
  if(result.status==='max_iterations'&&((result.x?.length||0)!==m.n||(result.z?.length||0)!==m.rowCount)&&(result.iterations>=this.iterationLimit||lastTimeLimit>0&&Number.isFinite(lastTimeLimit)&&result.solve_time>=lastTimeLimit*.8)){if(result.iterations>=this.iterationLimit)this.iterationLimit*=2;throw new OptimizationPending();}
  const raw=result.x||[],x=this.p.upper.map(()=>0);m.columns.forEach((i,j)=>x[i]=raw[j]);const violation=raw.length===m.n?marginalViolation(this.p,x):Infinity,y=Float64Array.from(result.z||[]);if(y.length!==m.rowCount)throw Error('Conic optimization returned no dual witness: '+result.status);
  for(let k=m.zero;k<m.zero+m.nonneg;k++)y[k]=Math.max(0,y[k]);let offset=m.zero+m.nonneg;for(const size of m.sizes){let norm=0;for(let j=1;j<size;j++)norm=Math.hypot(norm,y[offset+j]);y[offset]=Math.max(y[offset],norm);offset+=size;}
  const gradient=cost.slice();let quadraticValue=0;for(const form of this.p.quadraticCost||[]){const value=form.reduce((s,[i,w])=>s+w*x[i],0);quadraticValue+=value*value;for(const[i,w]of form)if(m.index[i]>=0)gradient[m.index[i]]+=2*w*value;}
  // Remove stationarity roundoff using existing bound multipliers before
  // paying a full box penalty. Reducing a lower-bound multiplier costs
  // nothing; reducing an upper-bound multiplier improves the dual bound.
  // Both operations preserve dual-cone feasibility and the original model.
  for(let i=0;i<m.n;i++){let residual=gradient[i];for(let j=m.starts[i];j<m.starts[i+1];j++)residual+=m.values[j]*y[m.indices[j]];if(residual<0)y[m.lowerBoxes[i]]-=Math.min(-residual,y[m.lowerBoxes[i]]);else if(m.upperBoxes[i]>=0)y[m.upperBoxes[i]]-=Math.min(residual,y[m.upperBoxes[i]]);}
  this.planes=[];offset=m.zero+m.nonneg;for(const q of this.p.quadratics||[]){const coefficients=new Map<number,number>();q.terms.forEach((row,k)=>{for(const[i,w]of row)coefficients.set(i,(coefficients.get(i)||0)-y[offset+k+1]*w);});const upper=Math.sqrt(q.upper)*y[offset];this.planes.push({terms:[...coefficients].filter(([,w])=>w!==0),upper:upper+128*Number.EPSILON*(1+Math.abs(upper))});offset+=q.terms.length+1;}
  let lower=-quadraticValue,magnitude=Math.abs(quadraticValue);for(let k=0;k<y.length;k++){const v=-m.b[k]*y[k];lower+=v;magnitude+=Math.abs(v);}for(let i=0;i<m.n;i++){let reduced=gradient[i];for(let j=m.starts[i];j<m.starts[i+1];j++)reduced+=m.values[j]*y[m.indices[j]];if(reduced<0){const v=reduced*m.upper[i];lower+=v;magnitude+=Math.abs(v);}}
  lower-=128*Number.EPSILON*(1+magnitude);
  // Squared objective forms are nonnegative. The original variable boxes
  // supply another global bound, independently of redundant row multipliers.
  let boxLower=0,correction=0,boxMagnitude=0;for(let i=0;i<m.n;i++)if(cost[i]<0){const term=cost[i]*m.upper[i],adjusted=term-correction,sum=boxLower+adjusted;correction=(sum-boxLower)-adjusted;boxLower=sum;boxMagnitude+=Math.abs(term);}
  if(boxMagnitude)boxLower-=64*Number.EPSILON*(1+boxMagnitude);
  const boxCertified=Number.isFinite(boxLower)&&boxLower>lower;if(boxCertified){lower=boxLower;this.planes=[];}
  const value=quadraticValue+cost.reduce((s,v,i)=>s+v*(x[m.columns[i]]||0),0),gap=value-lower;this.stats.queries++;this.stats.seconds+=(performance.now()-start)/1000;this.stats.iterations+=result.iterations||0;
  if(!Number.isFinite(violation)||violation>1e-7||!Number.isFinite(gap)||gap< -1e-7||gap>1e-7){
   if(result.status==='max_iterations'&&(result.iterations>=this.iterationLimit||lastTimeLimit>0&&Number.isFinite(lastTimeLimit)&&result.solve_time>=lastTimeLimit*.8)){if(result.iterations>=this.iterationLimit)this.iterationLimit*=2;throw new OptimizationPending();}
   this.retryableCertificate=Number.isFinite(violation)&&violation<=1e-7&&Number.isFinite(gap);this.stats.maximumRejectedGap=Math.max(this.stats.maximumRejectedGap,Math.abs(gap));throw Error('Conic optimization did not certify the original model: '+result.status+'; violation='+violation+'; primal-dual gap='+gap);
  }
  // The original feasible primal and independently repaired global dual are
  // the proof. A backend termination label cannot override that certificate.
  if(result.status!=='optimal')this.stats.statusIndependentCertificates++;
  if(boxCertified)this.stats.boxCertificates++;
  this.stats.maximumViolation=Math.max(this.stats.maximumViolation,violation);this.stats.maximumGap=Math.max(this.stats.maximumGap,gap);
  this.certifiedLower=lower;this.dual=this.p.rows.map(()=>0);for(let k=0;k<m.rowCount;k++)if(!boxCertified&&m.originalRows[k]>=0)this.dual[m.originalRows[k]]+=m.shadowSigns[k]*y[k];return{x,value:maximize?-value:value};
 }
 rowDuals(){return this.dual.slice();}
 supportPlanes(){return this.planes.map(p=>({terms:p.terms.map(t=>[...t]as LinearTerm),upper:p.upper}));}
 dispose(){this.compiled=undefined;this.dual=[];this.planes=[];this.certifiedLower=undefined;}
}
