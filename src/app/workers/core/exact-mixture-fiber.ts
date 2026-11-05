/** Exact rational arithmetic for equivalence certificates. Floating tolerances
 * are deliberately not used to decide whether two mixtures are identical. */
export interface ExactFraction {n:string;d:string}
interface Q {n:bigint;d:bigint}
const gcd=(a:bigint,b:bigint):bigint=>{a=a<0n?-a:a;while(b){const t=a%b;a=b;b=t;}return a;};
const q=(n:bigint,d=1n):Q=>{if(!n)return {n:0n,d:1n};if(d<0n){n=-n;d=-d;}const g=gcd(n,d);return {n:n/g,d:d/g};};
const add=(a:Q,b:Q)=>q(a.n*b.d+b.n*a.d,a.d*b.d),sub=(a:Q,b:Q)=>q(a.n*b.d-b.n*a.d,a.d*b.d),mul=(a:Q,b:Q)=>q(a.n*b.n,a.d*b.d),div=(a:Q,b:Q)=>q(a.n*b.d,a.d*b.n),cmp=(a:Q,b:Q)=>a.n*b.d-b.n*a.d;
const numberQ=(v:number):Q=>{if(!Number.isFinite(v))throw Error('Non-finite mixture coefficient.');if(!v)return q(0n);const view=new DataView(new ArrayBuffer(8));view.setFloat64(0,v);const bits=view.getBigUint64(0),exp=Number((bits>>52n)&2047n),mantissa=(bits&((1n<<52n)-1n))+(exp?1n<<52n:0n),power=(exp||1)-1023-52,sign=bits>>63n?-1n:1n;return power>=0?q(sign*(mantissa<<BigInt(power))):q(sign*mantissa,1n<<BigInt(-power));};
const encode=(v:Q):ExactFraction=>({n:String(v.n),d:String(v.d)}),decode=(v:ExactFraction)=>q(BigInt(v.n),BigInt(v.d));
export const exactNumber=(v:ExactFraction)=>Number(v.n)/Number(v.d);
export const exactProduct=(a:ExactFraction,b:ExactFraction)=>encode(mul(decode(a),decode(b)));
export const exactScale=(v:ExactFraction,numerator:number,denominator:number)=>encode(mul(decode(v),q(BigInt(numerator),BigInt(denominator))));
export const exactSubtract=(a:ExactFraction,b:ExactFraction)=>encode(sub(decode(a),decode(b)));
export const exactMin=(a:ExactFraction,b:ExactFraction)=>cmp(decode(a),decode(b))<=0?a:b;
export const exactSum=(values:ExactFraction[])=>encode(values.reduce((sum,v)=>add(sum,decode(v)),q(0n)));
export const exactFromNumber=(v:number)=>encode(numberQ(v));
export function exactWeightedEquality(a:{weight:ExactFraction;column:number[]}[],b:{weight:ExactFraction;column:number[]}[]):boolean{
 if(!a.length||!b.length||a[0].column.length!==b[0].column.length||[...a,...b].some(c=>c.column.length!==a[0].column.length))return false;
 const weightsA=a.map(c=>decode(c.weight)),weightsB=b.map(c=>decode(c.weight)),mass=(weights:Q[])=>weights.reduce(add,q(0n));if(cmp(mass(weightsA),mass(weightsB))!==0n)return false;
 const verified=new Set<string>();
 for(let k=0;k<a[0].column.length;k++){
  const common=a[0].column[k];if(a.every(c=>c.column[k]===common)&&b.every(c=>c.column[k]===common))continue;
  // Repeated channel patterns have exactly the same rational equation. Cache
  // the equation, never a tolerance or rounded fluorescence value.
  const pattern=JSON.stringify([...a.map(c=>c.column[k]),...b.map(c=>c.column[k])]);if(verified.has(pattern))continue;
  const sum=(rows:typeof a,weights:Q[])=>rows.reduce((s,c,i)=>add(s,mul(weights[i],numberQ(c.column[k]))),q(0n));if(cmp(sum(a,weightsA),sum(b,weightsB))!==0n)return false;verified.add(pattern);
 }
 return true;
}

export interface ExactHaplotype {sequence:string;weight:number;exactWeight:ExactFraction}
export interface HaplotypeRepresentatives {options:ExactHaplotype[][];verticesComplete:boolean;}
/** Feasible vertices of the same per-position marginal constraints. For large
 * tracts, deterministic quantile couplings give certified representatives.
 * The caller must additionally certify equality under its full emission model. */
export function exactHaplotypeRepresentatives(input:{sequence:string;weight:number}[],limit=16):HaplotypeRepresentatives{
 if(!input.length)return {options:[],verticesComplete:true};const length=input[0].sequence.length;
 if(!length||input.some(x=>x.sequence.length!==length||/[^ACGT]/.test(x.sequence)||x.weight<0))return {options:[],verticesComplete:false};
 const merged=new Map<string,Q>();for(const x of input)merged.set(x.sequence,add(merged.get(x.sequence)||q(0n),numberQ(x.weight)));
 const base=[...merged].filter(([,w])=>w.n>0n).sort(([a],[b])=>a.localeCompare(b)),mass=base.reduce((s,[,w])=>add(s,w),q(0n)),alphabet=Array.from({length},(_,p)=>[...new Set(base.map(([s])=>s[p]))].sort());
 const rowLabels=alphabet.flatMap((letters,p)=>letters.slice(0,-1).map(base=>({p,base}))),margins=alphabet.map((letters,p)=>new Map(letters.map(letter=>[letter,base.reduce((sum,[s,w])=>s[p]===letter?add(sum,w):sum,q(0n))]))),features=(word:string)=>[1,...rowLabels.map(({p,base})=>Number(word[p]===base))],rhs=[mass,...rowLabels.map(({p,base})=>margins[p].get(base)!)],rank=rhs.length;
 const toRows=(rows:[string,Q][])=>rows.filter(([,w])=>w.n>0n).sort(([a],[b])=>a.localeCompare(b)).map(([sequence,w])=>({sequence,weight:Number(w.n)/Number(w.d),exactWeight:encode(w)})),options:ExactHaplotype[][]=[toRows(base)],keys=new Set([JSON.stringify(options[0].map(x=>[x.sequence,x.exactWeight]))]);
 const keep=(rows:[string,Q][])=>{if(options.length>=limit)return;const result=toRows(rows),key=JSON.stringify(result.map(x=>[x.sequence,x.exactWeight]));if(!keys.has(key)){keys.add(key);options.push(result);}};
 let count=1;for(const a of alphabet){count*=a.length;if(count>16)break;}let verticesComplete=false;
 if(count<=16&&count>rank){let words=[''];for(const letters of alphabet)words=words.flatMap(prefix=>letters.map(b=>prefix+b));let combinations=1;for(let k=0;k<rank;k++)combinations=combinations*(words.length-k)/(k+1);
  if(combinations<=5000){let visited=0;const solve=(ids:number[])=>{const matrix=rhs.map((value,r)=>[...ids.map(id=>q(BigInt(features(words[id])[r]))),value]);for(let col=0;col<rank;col++){const pivot=matrix.findIndex((row,r)=>r>=col&&row[col].n!==0n);if(pivot<0)return;[matrix[col],matrix[pivot]]=[matrix[pivot],matrix[col]];const d=matrix[col][col];matrix[col]=matrix[col].map(v=>div(v,d));for(let r=0;r<rank;r++)if(r!==col){const scale=matrix[r][col];if(scale.n)matrix[r]=matrix[r].map((v,k)=>sub(v,mul(scale,matrix[col][k])));}}const values=matrix.map(row=>row[rank]);if(values.every(v=>v.n>=0n))keep(ids.map((id,i)=>[words[id],values[i]]));};
   const visit=(ids:number[],start:number)=>{if(options.length>=limit)return;if(ids.length===rank){visited++;solve(ids);return;}for(let i=start;i<=words.length-(rank-ids.length);i++)visit([...ids,i],i+1);};visit([],0);verticesComplete=visited===Math.round(combinations);
  }
 }
 // These couplings preserve marginals exactly without enumerating 4^length.
 if(!verticesComplete&&options.length<limit){const variable=alphabet.flatMap((a,p)=>a.length>1?[p]:[]),patterns=[new Set<number>(),new Set(variable),new Set(variable.filter((_,i)=>i%2===0)),...variable.map(p=>new Set([p]))];for(const reversed of patterns){if(options.length>=limit)break;const intervals=alphabet.map((letters,p)=>{let at=q(0n);return (reversed.has(p)?[...letters].reverse():letters).map(letter=>{const start=at;at=add(at,margins[p].get(letter)!);return {letter,start,end:at};});}),breaks=new Map<string,Q>();for(const rows of intervals)for(const x of rows)for(const point of [x.start,x.end])breaks.set(point.n+'/'+point.d,point);const sorted=[...breaks.values()].sort((a,b)=>cmp(a,b)<0n?-1:cmp(a,b)>0n?1:0),rows:[string,Q][]=[];for(let i=1;i<sorted.length;i++){const mid=div(add(sorted[i-1],sorted[i]),q(2n)),word=intervals.map(entries=>entries.find(x=>cmp(mid,x.start)>=0n&&cmp(mid,x.end)<0n)!.letter).join('');rows.push([word,sub(sorted[i],sorted[i-1])]);}keep(rows);}}
 return {options,verticesComplete};
}

/** Find a bounded mixture, then certify it with rational arithmetic. The
 * floating simplex only proposes active bounds; it can never certify equality.
 * A null result means no certificate was found, not a proof of infeasibility. */
export function exactBoundedMixture(columns:number[][],source:{weight:ExactFraction;column:number[]}[],bounds:{lower:ExactFraction;upper:ExactFraction}[]):ExactFraction[]|null{
 const n=columns.length;if(!n||bounds.length!==n||!source.length||columns.some(c=>c.length!==columns[0].length)||source.some(c=>c.column.length!==columns[0].length))return null;
 const lower=bounds.map(b=>decode(b.lower)),upper=bounds.map(b=>decode(b.upper));if(lower.some((l,i)=>l.n<0n||cmp(l,upper[i])>0n))return null;
 const equations:Q[][]=[],seen=new Set<string>();
 const keep=(coeff:number[],rhs:Q)=>{const key=JSON.stringify(coeff)+'='+rhs.n+'/'+rhs.d;if(!seen.has(key)){seen.add(key);equations.push([...coeff.map(numberQ),rhs]);}};
 keep(Array(n).fill(1),source.reduce((s,r)=>add(s,decode(r.weight)),q(0n)));
 for(let p=0;p<columns[0].length;p++)keep(columns.map(c=>c[p]),source.reduce((s,r)=>add(s,mul(decode(r.weight),numberQ(r.column[p]))),q(0n)));
 // Reduce redundant cycle/channel constraints before asking for a vertex.
 const basis:Q[][]=[],pivots:number[]=[];
 const insert=(input:Q[],required:boolean)=>{const row=input.slice();for(let i=0;i<basis.length;i++){const f=row[pivots[i]];if(f.n)for(let j=0;j<=n;j++)row[j]=sub(row[j],mul(f,basis[i][j]));}const pivot=row.findIndex((x,j)=>j<n&&x.n!==0n);if(pivot<0)return !required||row[n].n===0n;const scale=row[pivot];for(let j=0;j<=n;j++)row[j]=div(row[j],scale);for(let i=0;i<basis.length;i++){const f=basis[i][pivot];if(f.n)for(let j=0;j<=n;j++)basis[i][j]=sub(basis[i][j],mul(f,row[j]));}basis.push(row);pivots.push(pivot);return true;};
 for(const row of equations)if(!insert(row,true))return null;
 const reduced=basis.map(row=>row.slice()),reducedPivots=[...pivots],validate=(result:Q[])=>{if(result.some((x,i)=>cmp(x,lower[i])<0n||cmp(x,upper[i])>0n))return null;for(const row of equations)if(cmp(row.slice(0,n).reduce((s,x,i)=>add(s,mul(x,result[i])),q(0n)),row[n])!==0n)return null;return result.map(encode);};
 const rationalFallback=()=>{const free=Array.from({length:n},(_,i)=>i).filter(i=>!reducedPivots.includes(i)),constraints:Q[][]=[],rhs:Q[]=[];for(const [r,pivot]of reducedPivots.entries()){const row=reduced[r],coeff=free.map(j=>row[j]);constraints.push(coeff,coeff.map(v=>q(-v.n,v.d)));rhs.push(sub(row[n],lower[pivot]),sub(upper[pivot],row[n]));}for(const [j,id]of free.entries()){const row=free.map((_,k)=>q(j===k?1n:0n));constraints.push(row,row.map(v=>q(-v.n,v.d)));rhs.push(upper[id],q(-lower[id].n,lower[id].d));}const solution=free.length?rationalFeasibleSimplex(constraints,rhs):[];if(!solution)return null;const result=Array.from({length:n},()=>q(0n));free.forEach((id,i)=>result[id]=solution[i]);reduced.forEach((row,i)=>result[reducedPivots[i]]=sub(row[n],free.reduce((sum,id,k)=>add(sum,mul(row[id],solution[k])),q(0n))));return validate(result);};
 const value=(v:Q)=>Number(v.n)/Number(v.d),A:number[][]=[],b:number[]=[];
 for(const row of basis){const coefficients=row.slice(0,n).map(value),rhs=value(sub(row[n],row.slice(0,n).reduce((s,x,i)=>add(s,mul(x,lower[i])),q(0n))));A.push(coefficients,coefficients.map(v=>-v));b.push(rhs,-rhs);}
 for(let i=0;i<n;i++){A.push(Array.from({length:n},(_,j)=>Number(i===j)));b.push(value(sub(upper[i],lower[i])));}
 const proposal=feasibleSimplex(A,b);if(!proposal)return rationalFallback();
 const active=bounds.flatMap((_,i)=>[{i,v:lower[i],distance:Math.abs(proposal[i])},{i,v:upper[i],distance:Math.abs(proposal[i]-value(sub(upper[i],lower[i])))}]).sort((a,b)=>a.distance-b.distance);
 for(const bound of active){if(basis.length===n)break;if(bound.distance>1e-7)break;const row=Array.from({length:n+1},(_,j)=>j===n?bound.v:q(j===bound.i?1n:0n));insert(row,false);}
 if(basis.length!==n)return rationalFallback();const result=Array.from({length:n},()=>q(0n));basis.forEach((row,i)=>result[pivots[i]]=row[n]);return validate(result)||rationalFallback();
}
/** Phase-I simplex for a nonnegative feasibility proposal, Ax <= b. */
function feasibleSimplex(A:number[][],b:number[]):number[]|null{
 const m=b.length,n=A[0].length,B=Array.from({length:m},(_,i)=>n+i),N=Array.from({length:n+1},(_,i)=>i<n?i:-1),D=Array.from({length:m+2},()=>Array(n+2).fill(0));
 for(let i=0;i<m;i++){for(let j=0;j<n;j++)D[i][j]=A[i][j];D[i][n]=-1;D[i][n+1]=b[i];}D[m+1][n]=1;
 const pivot=(r:number,s:number)=>{const inv=1/D[r][s];for(let i=0;i<m+2;i++)if(i!==r)for(let j=0;j<n+2;j++)if(j!==s)D[i][j]-=D[r][j]*D[i][s]*inv;for(let j=0;j<n+2;j++)if(j!==s)D[r][j]*=inv;for(let i=0;i<m+2;i++)if(i!==r)D[i][s]*=-inv;D[r][s]=inv;[B[r],N[s]]=[N[s],B[r]];};
 let r=0;for(let i=1;i<m;i++)if(D[i][n+1]<D[r][n+1])r=i;
 if(D[r][n+1]<-1e-10){pivot(r,n);let done=false;for(let iteration=0;iteration<10000;iteration++){let s=0;for(let j=1;j<=n;j++)if(D[m+1][j]<D[m+1][s]-1e-10||Math.abs(D[m+1][j]-D[m+1][s])<=1e-10&&N[j]<N[s])s=j;if(D[m+1][s]>=-1e-10){done=true;break;}r=-1;for(let i=0;i<m;i++)if(D[i][s]>1e-10&&(r<0||D[i][n+1]/D[i][s]<D[r][n+1]/D[r][s]-1e-10||Math.abs(D[i][n+1]/D[i][s]-D[r][n+1]/D[r][s])<=1e-10&&B[i]<B[r]))r=i;if(r<0)return null;pivot(r,s);}if(!done||Math.abs(D[m+1][n+1])>1e-8)return null;}
 const x=Array(n).fill(0);for(let i=0;i<m;i++)if(B[i]>=0&&B[i]<n)x[B[i]]=D[i][n+1];return x;
}

/** Exact phase-I fallback on the reduced free-variable polytope. Bland's rule
 * avoids degenerate cycling. Exhaustion returns no certificate. */
function rationalFeasibleSimplex(A:Q[][],b:Q[]):Q[]|null{
 const m=b.length,n=A[0]?.length||0,B=Array.from({length:m},(_,i)=>n+i),N=Array.from({length:n+1},(_,i)=>i<n?i:-1),D=Array.from({length:m+1},()=>Array.from({length:n+2},()=>q(0n)));
 for(let i=0;i<m;i++){for(let j=0;j<n;j++)D[i][j]=A[i][j];D[i][n]=q(-1n);D[i][n+1]=b[i];}D[m][n]=q(1n);
 const pivot=(r:number,s:number)=>{const d=D[r][s];for(let i=0;i<=m;i++)if(i!==r&&D[i][s].n){const f=div(D[i][s],d);for(let j=0;j<n+2;j++)if(j!==s&&D[r][j].n)D[i][j]=sub(D[i][j],mul(f,D[r][j]));D[i][s]=q(-f.n,f.d);}for(let j=0;j<n+2;j++)if(j!==s)D[r][j]=div(D[r][j],d);D[r][s]=div(q(1n),d);[B[r],N[s]]=[N[s],B[r]];};
 let r=0;for(let i=1;i<m;i++)if(cmp(D[i][n+1],D[r][n+1])<0n)r=i;
 if(D[r][n+1].n<0n){pivot(r,n);let done=false;for(let iteration=0;iteration<5000;iteration++){let s=-1;for(let j=0;j<=n;j++)if(D[m][j].n<0n&&(s<0||N[j]<N[s]))s=j;if(s<0){done=true;break;}r=-1;for(let i=0;i<m;i++)if(D[i][s].n>0n&&(r<0||cmp(div(D[i][n+1],D[i][s]),div(D[r][n+1],D[r][s]))<0n||cmp(div(D[i][n+1],D[i][s]),div(D[r][n+1],D[r][s]))===0n&&B[i]<B[r]))r=i;if(r<0)return null;pivot(r,s);}if(!done||D[m][n+1].n!==0n)return null;}
 const x=Array.from({length:n},()=>q(0n));for(let i=0;i<m;i++)if(B[i]>=0&&B[i]<n)x[B[i]]=D[i][n+1];return x;
}
