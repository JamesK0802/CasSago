import type{SangerMarginalGroup,SangerMarginalDisplayComponent,SangerJointEvent}from'../../models/sanger-analysis';
interface Block{sites:number[];width:number;minimum:number;marker:number}
/** Representative partition ONLY. Adjacent substitutions share one block only
 * at the same measured fraction resolution. Different fraction levels get
 * disjoint representative populations, never a nested long mutation plus a
 * remainder. Every original fitted base marginal and hypothesis is retained.
 * If disjoint blocks cannot fit the family, retain the unsplit mixture.
 */
export function substitutionBlockDisplay(group:SangerMarginalGroup,baseline:(position:number)=>number,toRead:(position:number)=>number,reverse:boolean,noiseWithinReference:number,threshold:number){
 const sites=group.bases,amounts=sites.map(p=>p.estimate.reduce((n,v,b)=>n+(b===baseline(p.position)?0:Math.max(0,v)),0)),resolution=Math.max(1e-7,noiseWithinReference),blocks:Block[]=[];
 if(sites.filter((_,i)=>group.fraction*amounts[i]+1e-12>=threshold).length<2)return;
 for(let i=0;i<sites.length;i++){
  if(!(amounts[i]>0))continue;const amount=group.fraction*amounts[i],last=blocks.at(-1),previous=last?.sites.at(-1),adjacent=previous!==undefined&&Math.abs(sites[previous].position-sites[i].position)===1;
  if(last&&adjacent&&Math.max(last.width,amount)-Math.min(last.minimum,amount)<=resolution){last.sites.push(i);last.minimum=Math.min(last.minimum,amount);if(amount>last.width){last.width=amount;last.marker=i;}}
  else blocks.push({sites:[i],width:amount,minimum:amount,marker:i});
 }
 const total=blocks.reduce((n,b)=>n+b.width,0);if(total>group.fraction+1e-10)return;
 const components:SangerMarginalDisplayComponent[]=[],pure=(i:number)=>[0,1,2,3].map(b=>Number(b===baseline(sites[i].position))),event=(i:number,b:number):SangerJointEvent=>({target:group.target,start:sites[i].position+group.baselineReadOffset!,end:sites[i].position+group.baselineReadOffset!+1,inserted:1,insertedSequence:'ACGT'[b]});
 const add=(id:string,fraction:number,estimates:number[][],selection:SangerMarginalDisplayComponent['selection'])=>{
  if(!(fraction>0))return;const events=estimates.flatMap((v,i)=>{const b=v.indexOf(Math.max(...v));return b===baseline(sites[i].position)?[]:[event(i,b)];}).sort((a,b)=>a.start-b.start);
  components.push({id,rangeStatus:'estimated',fraction,fractionRange:[0,1],events,mixture:!selection.unmodified,selection,bases:sites.map((p,i)=>({...p,estimate:estimates[i]}))});
 };
 add(group.id+'::block-wt',Math.max(0,group.fraction-total),sites.map((_,i)=>pure(i)),{unmodified:true});
 for(const block of blocks){
  const marker=sites[block.marker],wt=baseline(marker.position);let lo=0;
  for(let b=0;b<4;b++)if(b!==wt){const weight=group.fraction*Math.max(0,marker.estimate[b]),hi=lo+weight;if(weight>0){
    const estimates=sites.map((_,i)=>{if(!block.sites.includes(i))return pure(i);const values=[0,0,0,0];let at=0;for(let a=0;a<4;a++)if(a!==baseline(sites[i].position)){const end=at+group.fraction*Math.max(0,sites[i].estimate[a]);values[a]=Math.max(0,Math.min(hi,end)-Math.max(lo,at))/weight;at=end;}values[baseline(sites[i].position)]=Math.max(0,1-values.reduce((s,v)=>s+v,0));return values;});
    add(group.id+'::block:'+marker.position+':'+b,weight,estimates,{marker:{position:toRead(marker.position),bases:[reverse?3-b:b]}});
   }lo=hi;}
 }
 // A display partition must conserve the original signal; otherwise use the
 // unchanged parent mixture. Numerical fit tolerances are never data pruning.
 const error=Math.max(0,...sites.flatMap((p,i)=>p.estimate.map((v,b)=>Math.abs(group.fraction*v-components.reduce((n,c)=>n+c.fraction*c.bases[i].estimate[b],0)))));
 if(error>1e-7||Math.abs(components.reduce((n,c)=>n+c.fraction,0)-group.fraction)>1e-7)return;
 return{components,unmodifiedEstimate:Math.max(0,group.fraction-total)};
}
