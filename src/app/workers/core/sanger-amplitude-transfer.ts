/** Detector colour gains learned from paired WT-prefix photon amplitudes.
 * A common linear read-position trend is a nuisance parameter. No edited
 * residual, allele, genotype count or ground truth is used. */
export function amplitudeTransfer(control:{vectors:number[][];strength:number[]},sample:{vectors:number[][];strength:number[]},positions:number[],controlTrim:number,sampleTrim:number){
 const first=positions[0],span=Math.max(1,positions.at(-1)!-first),rows:{x:number[];y:number;weight:number}[]=[],counts=[0,0,0,0];
 for(const p of positions){const c=p-controlTrim,s=p-sampleTrim;if(!control.strength[c]||!sample.strength[s])continue;for(let b=0;b<4;b++){const cv=control.vectors[c][b],sv=sample.vectors[s][b];if(cv<.05||sv<.05)continue;const x=Array(5).fill(0);x[b]=1;x[4]=2*(p-first)/span-1;rows.push({x,y:Math.log(sample.strength[s]*sv/(control.strength[c]*cv)),weight:cv*sv});counts[b]++;}}
 let fit=Array(5).fill(0),weights=rows.map(r=>r.weight);
 for(let round=0;round<5;round++){const a=Array.from({length:5},(_,i)=>Array.from({length:6},(_,j)=>j===i?1e-8:0));for(let k=0;k<rows.length;k++){const r=rows[k],w=weights[k];for(let i=0;i<5;i++){for(let j=0;j<5;j++)a[i][j]+=w*r.x[i]*r.x[j];a[i][5]+=w*r.x[i]*r.y;}}
  for(let i=0;i<5;i++){let pivot=i;for(let j=i+1;j<5;j++)if(Math.abs(a[j][i])>Math.abs(a[pivot][i]))pivot=j;[a[i],a[pivot]]=[a[pivot],a[i]];const d=a[i][i];for(let j=i;j<6;j++)a[i][j]/=d;for(let r=0;r<5;r++)if(r!==i){const d=a[r][i];for(let j=i;j<6;j++)a[r][j]-=d*a[i][j];}}fit=a.map(r=>r[5]);const errors=rows.map(r=>Math.abs(r.y-r.x.reduce((n,v,i)=>n+v*fit[i],0))).sort((a,b)=>a-b),scale=Math.max(.01,1.4826*(errors[Math.floor(errors.length/2)]||0));weights=rows.map((r,i)=>r.weight*Math.min(1,1.5*scale/Math.max(1e-12,Math.abs(r.y-r.x.reduce((n,v,i)=>n+v*fit[i],0)))));
 }
 const mean=fit.slice(0,4).reduce((n,v)=>n+v,0)/4;return{gains:fit.slice(0,4).map(v=>Math.exp(v-mean)),slope:fit[4],logScale:mean,first,span,counts};
}
