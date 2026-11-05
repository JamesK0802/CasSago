/** Joint proposals for contiguous substitution molecules. Input is fluorescence
 * and inferred candidates only. No truth recipes or group percentages enter. */
export function substitutionProposals(columns:number[][],wt:number[],residual:number[],mass:number,limit=8):number[][]{
 const sparse=columns.map(c=>{const at:number[]=[],values:number[]=[];let rhs=0,norm=0;for(let k=0;k<c.length;k++){const v=c[k]-wt[k];if(Math.abs(v)<1e-10)continue;at.push(k);values.push(v);rhs+=v*residual[k];norm+=v*v;}return{at,values,rhs,norm};});
 const dot=(i:number,j:number)=>{const a=sparse[i],b=sparse[j];let x=0,y=0,sum=0;while(x<a.at.length&&y<b.at.length){if(a.at[x]===b.at[y]){sum+=a.values[x++]*b.values[y++];}else if(a.at[x]<b.at[y])x++;else y++;}return sum;};
 const best:{ids:number[];gain:number}[]=[];let floor=-Infinity;
 const add=(ids:number[],gain:number)=>{if(gain<=floor)return;best.push({ids,gain});best.sort((a,b)=>b.gain-a.gain);if(best.length>limit)best.pop();if(best.length===limit)floor=best.at(-1)!.gain;};
 for(let i=0;i<sparse.length;i++){const a=sparse[i];if(!a.norm||a.rhs<=0)continue;const u=Math.min(mass,a.rhs/a.norm);add([i],2*u*a.rhs-u*u*a.norm);
  for(let j=0;j<i;j++){const b=sparse[j];if(!b.norm||b.rhs<=0)continue;const cross=dot(i,j),det=a.norm*b.norm-cross*cross;if(det<1e-10)continue;let x=(a.rhs*b.norm-b.rhs*cross)/det,y=(b.rhs*a.norm-a.rhs*cross)/det;if(x<=0||y<=0)continue;if(x+y>mass){const curvature=a.norm+b.norm-2*cross;x=Math.max(0,Math.min(mass,curvature>0?(a.rhs-b.rhs+mass*(b.norm-cross))/curvature:mass/2));y=mass-x;}add([i,j],2*(x*a.rhs+y*b.rhs)-x*x*a.norm-y*y*b.norm-2*x*y*cross);
  }
 }
 return best.map(b=>b.ids);
}
