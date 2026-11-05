import type{JunctionInput,JunctionFit,JunctionModel}from'./sanger-junction-mixture';
import{junctionStarts}from'./sanger-junction-mixture';
import{jointPathColumn,PathEvent}from'./sanger-junction-paths';
const bases='ACGT';
function molecule(sequence:string,events:PathEvent[]){let out=sequence;for(const e of [...events].reverse()){const o=e.origins[0],insert=e.size>0?o.probabilities!.map(v=>bases[v.indexOf(Math.max(...v))]).join(''):'';out=out.slice(0,o.start)+insert+out.slice(o.start+Math.max(0,-e.size));}return out;}
/** Search physical recombinations that change target-size marginals but leave
 * every measured four-channel emission exactly unchanged. This is an explicit
 * non-identifiability witness, not a posterior probability or exhaustive proof
 * of identifiability when none is found. Point estimates are never withheld. */
export function findPathAliases(input:JunctionInput,fit:JunctionFit){
 const reports:{reference:number;target:number;shiftableJointFraction:number}[]=[];
 input.references.forEach((ref,h)=>{
  const cuts=ref.cuts||[ref.cut];if(cuts.length!==2)return;const first=cuts[0]<cuts[1]?0:1,last=1-first;
  const hard=fit.models.flatMap((m,index)=>{
   if(m.reference!==h||fit.weights[index]<=1e-8||!m.sizes)return[];
   let mass=fit.weights[index];const path=(m.path||[]).map(e=>{const total=e.origins.reduce((s,o)=>s+o.weight,0),o=e.origins.reduce((a,b)=>a.weight>b.weight?a:b);mass*=o.weight/Math.max(1e-12,total);const probabilities=o.probabilities?.map(v=>{const best=v.indexOf(Math.max(...v));mass*=v[best];return v.map((_,b)=>Number(b===best));});return{target:e.target,size:e.size,origins:[{start:o.start,weight:1,probabilities}]};}),model={...m,path},early=path.filter(e=>e.target===first+1),baseline=molecule(ref.sequence,early);
   return[{model,index,mass,early,baseline,sequence:molecule(ref.sequence,path)}];
  });
  const find=(row:typeof hard[number],sequence:string):JunctionModel|null=>{
   const size=sequence.length-row.baseline.length;if(size< -input.maxDeletion||size>input.maxInsertion)return null;
   let event:PathEvent|undefined;
   if(!size){if(sequence!==row.baseline)return null;}
   else{
    const prefix=row.model.sizes![first];
    for(const start of junctionStarts(cuts[last],size,input.radius,ref.sequence.length)){
     if(row.early.some(e=>start<e.origins[0].start+Math.max(0,-e.size)))continue;
     const p=start+prefix,insert=size>0?sequence.slice(p,p+size):'';
     if(row.baseline.slice(0,p)+insert+row.baseline.slice(p+Math.max(0,-size))!==sequence)continue;
     event={target:last+1,size,origins:[{start,weight:1,probabilities:size>0?Array.from(insert,b=>Array.from(bases,v=>Number(v===b))):undefined}]};break;
    }
    if(!event)return null;
   }
   const sizes=[...row.model.sizes!];sizes[last]=size;
   return{...row.model,sizes,size:sizes.reduce((s,v)=>s+v,0),path:[...row.early,...(event?[event]:[])]};
  };
  let shifted=0;
  for(let a=0;a<hard.length;a++)for(let b=0;b<a;b++){
   const A=hard[a],B=hard[b],mass=Math.min(A.mass,B.mass);if(mass<=shifted||A.model.sizes![first]===B.model.sizes![first])continue;
   const firstEnd=Math.max(...[A,B].flatMap(row=>row.early.map(e=>e.origins[0].start+Math.max(0,e.size))) ,0);
   const lo=Math.max(firstEnd+3,cuts[last]-input.maxDeletion-input.radius-input.maxInsertion),hi=Math.min(A.sequence.length,B.sequence.length,cuts[last]+input.radius+2*input.maxInsertion);
   for(let p=lo;p<=hi;p++){
    const C=find(A,A.sequence.slice(0,p)+B.sequence.slice(p)),D=find(B,B.sequence.slice(0,p)+A.sequence.slice(p));if(!C||!D)continue;
    const before=[A.model.sizes![last],B.model.sizes![last]].sort((a,b)=>a-b).join('/'),after=[C.sizes![last],D.sizes![last]].sort((a,b)=>a-b).join('/');if(before===after)continue;
    const oldA=jointPathColumn(input,A.model),oldB=jointPathColumn(input,B.model),newA=jointPathColumn(input,C),newB=jointPathColumn(input,D);
    if(oldA.every((v,k)=>Math.abs(v+oldB[k]-newA[k]-newB[k])<1e-9)){shifted=mass;break;}
   }
  }
  if(shifted>0)reports.push({reference:h,target:last+1,shiftableJointFraction:Math.min(1,2*shifted)});
 });return reports;
}
