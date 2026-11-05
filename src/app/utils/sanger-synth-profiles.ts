import type {SynthReference,SynthAllele,SynthEvent} from './sanger-synthetic-benchmark';
export interface SynthTargetProfile {target:number;name:string;alleles:SynthAllele[]}
export function initSynthProfiles(ref:SynthReference){ref.profiles=ref.cuts.map((_,target)=>({target,name:ref.targetNames[target]||'Target '+(target+1),alleles:[]}));return ref.profiles;}
export function targetWT(profile:SynthTargetProfile){return 100-profile.alleles.reduce((s,a)=>s+Number(a.percent),0);}
/** Finite random pairing preserves every target marginal exactly to 0.01%.
 * It is an explicit simulation coupling, not inferred biological linkage. */
export function mixTargetProfiles(ref:SynthReference,random:()=>number){
 const profiles=ref.profiles!;if(profiles.length!==ref.cuts.length||profiles.some((p,i)=>p.target!==i))throw Error(ref.name+': give each target one profile.');
 const population=10000,axes=profiles.map(profile=>{
  if(!Number.isFinite(targetWT(profile))||targetWT(profile)<-1e-8)throw Error(ref.name+' '+profile.name+': group percentages must total at most 100%.');
  const values=new Uint16Array(population);let at=0;
  profile.alleles.forEach((allele,i)=>{const percent=Number(allele.percent),n=Math.round(percent*100);if(!Number.isFinite(percent)||percent<0||percent>100||Math.abs(n-percent*100)>1e-6||allele.events.length!==1||allele.events[0].target!==profile.target)throw Error('Each target group needs one edit and a percentage with at most two decimals.');values.fill(i+1,at,at+n);at+=n;});
  for(let j=population-1;j>0;j--){const k=Math.min(j,Math.floor(random()*(j+1))),v=values[j];values[j]=values[k];values[k]=v;}return values;
 });
 const counts=new Map<string,number>();for(let i=0;i<population;i++){const key=axes.map(a=>a[i]).join('/');counts.set(key,(counts.get(key)||0)+1);}
 return [...counts].sort((a,b)=>a[0].localeCompare(b[0],undefined,{numeric:true})).map(([key,count])=>{const choices=key.split('/').map(Number),events:SynthEvent[]=choices.flatMap((group,t)=>group?structuredClone(profiles[t].alleles[group-1].events):[]);return {allele:events.length?choices.map((group,t)=>'T'+(t+1)+' '+(group?'Group '+group:'WT')).join(' · '):'WT',percent:count/100,events,targetStates:choices};});
}
