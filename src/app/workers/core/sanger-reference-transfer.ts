/** Paired unedited-prefix detector transfer, with a dye gain and common
 * read-position trend for each submitted reference. Fits absolute photons
 * before normalization; the same loop handles one and multiple references.
 * Only pre-cut training rows are read. Validation rows are never fitted. */
export function referenceTransfer(responses:number[][][],shares:number[],sampleAreas:number[][],positions:number[],sampleTrim:number,scale:number,initial:{gains:number[];slope:number;logScale:number;first:number;span:number}){
 const hCount=shares.length,n=5*hCount,first=initial.first,span=initial.span,at=(p:number)=>2*(p-first)/span-1,rows:{p:number;b:number;y:number;weight:number;source:number[]}[]=[];
 for(const p of positions){const area=sampleAreas[p-sampleTrim];if(!area)continue;const total=area.reduce((a,b)=>a+b,0),c=[0,1,2,3].map(b=>responses.reduce((s,r,h)=>s+shares[h]*r[p][b],0)),ctotal=c.reduce((a,b)=>a+b,0);for(let b=0;b<4;b++){const cv=c[b]/ctotal,sv=area[b]/total;if(cv<.05||sv<.05)continue;rows.push({p,b,y:area[b]/scale,weight:cv*sv,source:responses.map((r,h)=>shares[h]*r[p][b])});}}
 let parameters=shares.flatMap(()=>[...initial.gains.map(v=>Math.log(v)+initial.logScale),initial.slope]),robust=rows.map(r=>r.weight);
 const one=(r:typeof rows[number],theta:number[])=>r.source.map((v,h)=>v*Math.exp(theta[5*h+r.b]+theta[5*h+4]*at(r.p))),loss=(theta:number[])=>rows.reduce((s,r,k)=>{const predicted=one(r,theta).reduce((a,b)=>a+b,0),error=Math.log(Math.max(1e-30,predicted)/r.y);return s+robust[k]*error*error;},0);
 for(let outer=0;outer<5;outer++){
  for(let step=0;step<40;step++){
   const matrix=Array.from({length:n},(_,i)=>Array.from({length:n+1},(_,j)=>j===i?1e-8:0));
   rows.forEach((r,k)=>{const v=one(r,parameters),predicted=v.reduce((a,b)=>a+b,0),error=Math.log(Math.max(1e-30,predicted)/r.y),jacobian=Array(n).fill(0);for(let h=0;h<hCount;h++){jacobian[5*h+r.b]=v[h]/Math.max(1e-30,predicted);jacobian[5*h+4]=jacobian[5*h+r.b]*at(r.p);}for(let i=0;i<n;i++){for(let j=0;j<n;j++)matrix[i][j]+=robust[k]*jacobian[i]*jacobian[j];matrix[i][n]-=robust[k]*jacobian[i]*error;}});
   for(let i=0;i<n;i++){let pivot=i;for(let j=i+1;j<n;j++)if(Math.abs(matrix[j][i])>Math.abs(matrix[pivot][i]))pivot=j;[matrix[i],matrix[pivot]]=[matrix[pivot],matrix[i]];const divisor=matrix[i][i];for(let j=i;j<=n;j++)matrix[i][j]/=divisor;for(let r=0;r<n;r++)if(r!==i){const d=matrix[r][i];for(let j=i;j<=n;j++)matrix[r][j]-=d*matrix[i][j];}}
   const delta=matrix.map(r=>r[n]);if(Math.max(...delta.map(Math.abs))<1e-9)break;const before=loss(parameters);let alpha=1,accepted=false;while(alpha>=1/65536){const next=parameters.map((v,j)=>v+alpha*delta[j]);if(loss(next)<before){parameters=next;accepted=true;break;}alpha*=.5;}if(!accepted)break;
  }
  const errors=rows.map(r=>Math.abs(Math.log(Math.max(1e-30,one(r,parameters).reduce((a,b)=>a+b,0))/r.y))).sort((a,b)=>a-b),noise=Math.max(.01,1.4826*(errors[Math.floor(errors.length/2)]||0));robust=rows.map(r=>r.weight*Math.min(1,1.5*noise/Math.max(1e-12,Math.abs(Math.log(Math.max(1e-30,one(r,parameters).reduce((a,b)=>a+b,0))/r.y)))));
 }
 const residualLogVariance=rows.reduce((s,r)=>{const error=Math.log(Math.max(1e-30,one(r,parameters).reduce((a,b)=>a+b,0))/r.y);return s+r.weight*error*error;},0)/Math.max(1e-30,rows.reduce((s,r)=>s+r.weight,0));
 const quantizationLogVariance=rows.reduce((s,r)=>s+r.weight*(1/Math.max(1,r.y*scale)**2+1/Math.max(1,r.source.reduce((a,b)=>a+b,0)*scale)**2),0)/Math.max(1e-30,rows.reduce((s,r)=>s+r.weight,0)),logVariance=Math.max(residualLogVariance,quantizationLogVariance);
 return{first,span,parameters,logVariance,quantizationLogVariance,trainingRows:rows.length,apply(h:number,p:number,dye:number[]){return dye.map((v,b)=>v*Math.exp(parameters[5*h+b]-Math.log(initial.gains[b])+parameters[5*h+4]*at(p)));}};
}
