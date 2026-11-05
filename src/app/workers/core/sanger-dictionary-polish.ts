/** Joint sparse dictionary refit used as an initialization, not a claim of
 * unique allele phasing. Family totals come from the measured-flank fit;
 * the caller subsequently frees every retained weight in the full model. */
export function polishSignalDictionary(columns:number[][],observed:number[],families:string[],selected:number[],initial:number[],sweeps=160){
 const weights=columns.map(()=>0);selected.forEach((id,j)=>weights[id]=initial[j]);
 const groups=new Map<string,{ids:number[];anchor:number}>();for(let id=0;id<columns.length;id++){const key=families[id];let group=groups.get(key);if(!group){group={ids:[],anchor:id};groups.set(key,group);}group.ids.push(id);if(weights[id]>weights[group.anchor])group.anchor=id;}
 const activeGroups=[...groups.values()].filter(g=>g.ids.some(id=>weights[id]>0));
 const sparse=new Map<number,{at:number[];values:number[]}>();for(const g of activeGroups)for(const id of g.ids){const at:number[]=[],values:number[]=[];for(let k=0;k<observed.length;k++){const d=columns[id][k]-columns[g.anchor][k];if(Math.abs(d)>1e-10){at.push(k);values.push(d);}}sparse.set(id,{at,values});}
 for(let outer=0;outer<3;outer++){
  const prediction=observed.map(()=>0);for(let id=0;id<weights.length;id++)if(weights[id]>0)for(let k=0;k<observed.length;k++)prediction[k]+=weights[id]*columns[id][k];
  const totals=Array.from({length:observed.length/4},(_,p)=>Math.max(1e-8,prediction[4*p]+prediction[4*p+1]+prediction[4*p+2]+prediction[4*p+3])),residual=observed.map((v,k)=>v*totals[k>>2]-prediction[k]),rowWeights=observed.map((_,k)=>1/totals[k>>2]**2);
  const gradient=(id:number)=>{const c=sparse.get(id)!;let g=0;for(let j=0;j<c.at.length;j++){const k=c.at[j];g+=c.values[j]*residual[k]*rowWeights[k];}return g;};
  for(let sweep=0;sweep<sweeps;sweep++){let maximumGain=0;
   for(const group of activeGroups){let best=group.anchor,worst=group.ids.find(id=>weights[id]>1e-12)!,bestScore=-Infinity,worstScore=Infinity;
    for(const id of group.ids){const g=gradient(id);if(g>bestScore){best=id;bestScore=g;}if(weights[id]>1e-12&&g<worstScore){worst=id;worstScore=g;}}
    if(best===worst||bestScore-worstScore<1e-9)continue;
    const a=sparse.get(best)!,b=sparse.get(worst)!,at:number[]=[],values:number[]=[];let i=0,j=0,norm=0,rhs=0;
    while(i<a.at.length||j<b.at.length){const k=Math.min(a.at[i]??Infinity,b.at[j]??Infinity),v=(a.at[i]===k?a.values[i++]:0)-(b.at[j]===k?b.values[j++]:0);if(!v)continue;at.push(k);values.push(v);norm+=v*v*rowWeights[k];rhs+=v*residual[k]*rowWeights[k];}
    if(norm<1e-15)continue;const step=Math.max(0,Math.min(weights[worst],rhs/norm));if(step<1e-12)continue;weights[best]+=step;weights[worst]-=step;for(let k=0;k<at.length;k++)residual[at[k]]-=step*values[k];maximumGain=Math.max(maximumGain,2*step*rhs-step*step*norm);
   }
   if(maximumGain<1e-11)break;
  }
 }
 const ids=weights.flatMap((w,id)=>w>=.0005?[id]:[]);return {ids,weights:ids.map(id=>weights[id])};
}
