import{materialize,type Prepared,type PathEdit}from'./sanger-path-common';
import{normalizeMoleculeEdits}from'./sanger-edit-normalize';
import type{SangerLinearOptimizer}from'./sanger-discrete-phase';
type Term=[number,number];
interface Variable{cost:number;upper:number;integer:number}
/** One bounded-noise, minimum-support model for every profile size. Abundances,
 * legal origins and insertion words are optimized together. Capacity grows by
 * net-length class, never by enumerating frequency partitions/count allocations. */
export function fitSparseFrequencyProfile(ref:Prepared,families:{size:number;fraction:number}[],positions:number[],observed:number[],noise:number,budget:number,progress:(s:string)=>void,_threshold:number,optimizer?:SangerLinearOptimizer){
 if(!optimizer?.createModel)throw Error('Initialize the common sequence optimizer.');
 const started=performance.now(),alphabet='ACGT',cut=ref.cuts[0],span=10,left=cut-3,right=cut+3+Math.max(span,...families.map(f=>f.size)),ps=positions.filter(p=>p>=left&&p<right),lookup=new Map(positions.map((p,i)=>[p,i])),target=ps.flatMap(p=>observed.slice(4*lookup.get(p)!,4*lookup.get(p)!+4)),tolerance=Math.max(1e-5,3*noise),INF=1e30,classes=families.filter(f=>f.fraction>0),details:any[]=[];
 for(let capacity=1;;capacity++){
  const variables:Variable[]=[],constraints:{terms:Term[];lower:number;upper:number}[]=[],prediction=ps.map(()=>Array.from({length:4},()=>[]as Term[])),groups:{family:number;size:number;active:number;origins:{start:number;bit:number;mass:number}[];bases:number[][]}[]=[],wt:{family:number;mass:number}[]=[],classMass:Term[][]=classes.map(()=>[]);
  const variable=(cost=0,upper=1,integer=0)=>{const id=variables.length;variables.push({cost,upper,integer});return id;},row=(terms:Term[],lower:number,upper:number)=>constraints.push({terms,lower,upper});
  classes.forEach((family,f)=>{
   if(family.size===0){const mass=variable(0,family.fraction);wt.push({family:f,mass});classMass[f].push([mass,1]);ps.forEach((p,i)=>prediction[i][alphabet.indexOf(ref.sequence[p])].push([mass,1]));}
   const previous:number[][]=[];
   for(let slot=0;slot<capacity;slot++){
    const active=variable(1,1,1),length=family.size>0?family.size:family.size===0?span:0,bases=Array.from({length},()=>{const bits=Array.from({length:4},()=>variable(0,1,1));row([...bits.map(id=>[id,1]as Term),[active,-1]],0,0);return bits;}),origins:{start:number;bit:number;mass:number}[]=[];
    for(let start=cut-3;start<=cut+3;start++){
     if(start<0||start+Math.max(0,-family.size)>ref.sequence.length)continue;if(family.size<0&&start>cut-3&&ref.sequence[start-1]===ref.sequence[start-1-family.size])continue;
     const bit=variable(0,1,1),mass=variable(0,family.fraction);origins.push({start,bit,mass});classMass[f].push([mass,1]);row([[mass,1],[bit,-family.fraction]],-INF,0);
     if(length&&start>cut-3){const forbidden=alphabet.indexOf(ref.sequence[family.size>0?start-1:start+length-1]);row([[bit,1],[bases[length-1][forbidden],1]],-INF,1);}
     if(family.size===0){const matches=bases.map((bs,k)=>[bs[alphabet.indexOf(ref.sequence[start+k])],1]as Term);row([...matches,[bit,span]],-INF,2*span-1);}
     const cells=new Map<number,number[]>();
     for(let k=0;k<length;k++){const flow=Array.from({length:4},(_,b)=>{const id=variable(0,family.fraction);row([[id,1],[bases[k][b],-family.fraction]],-INF,0);return id;});cells.set(start+k,flow);row([...flow.map(id=>[id,1]as Term),[mass,-1]],0,0);}
     ps.forEach((p,i)=>{const flow=cells.get(p);if(flow)flow.forEach((id,b)=>prediction[i][b].push([id,1]));else{const at=p<start?p:p-family.size,b=alphabet.indexOf(ref.sequence[at]);if(b>=0)prediction[i][b].push([mass,1]);}});
    }
    row([...origins.map(o=>[o.bit,1]as Term),[active,-1]],0,0);
    // Slots are unlabeled; sorting their FREE abundances removes permutations
    // without discarding a single allele/frequency association.
    const mass=origins.map(o=>o.mass);if(previous.length)row([...previous.at(-1)!.map(id=>[id,1]as Term),...mass.map(id=>[id,-1]as Term)],0,INF);previous.push(mass);groups.push({family:f,size:family.size,active,origins,bases});
   }
  });
  classes.forEach((f,i)=>row(classMass[i],f.fraction,f.fraction));
  const secondary=.01/(target.length*tolerance);prediction.forEach((p,i)=>p.forEach((terms,b)=>{const y=target[4*i+b],e=variable(secondary,tolerance);row(terms,y-tolerance,y+tolerance);row([...terms,[e,1]],y,INF);row([...terms,[e,-1]],-INF,y);}));
  const starts=[0],indices:number[]=[],values:number[]=[],lower:number[]=[],upper:number[]=[];for(const r of constraints){for(const[id,w]of r.terms){indices.push(id);values.push(w);}starts.push(indices.length);lower.push(r.lower);upper.push(r.upper);}
  const model:any=optimizer.createModel();let solution:number[]=[],status:any,bound:number|undefined;
  progress('Solving common sparse profiles · '+capacity+' slots per length class · '+variables.length+' variables');
  try{model.options.set({output_flag:false,time_limit:540,mip_abs_gap:.05,mip_rel_gap:0,mip_feasibility_tolerance:1e-7,primal_feasibility_tolerance:1e-8});model.passModel({numCols:variables.length,numRows:constraints.length,colCost:variables.map(v=>v.cost),colLower:variables.map(()=>0),colUpper:variables.map(v=>v.upper),integrality:variables.map(v=>v.integer),rowLower:lower,rowUpper:upper,matrix:{format:'csr',numCols:variables.length,numRows:constraints.length,starts,indices,values}});let lastUpdate=-Infinity;status=model.run({6:(event:any)=>{const now=performance.now();if(now-lastUpdate<2000)return;lastUpdate=now;const d=event.data,format=(value:number)=>Number.isFinite(value)?value.toFixed(3):'pending';progress('Optimizing common sequence profiles · '+((now-started)/1000).toFixed(1)+' seconds · '+d.mip_node_count+' nodes · lower bound '+format(d.mip_dual_bound)+' · best '+format(d.mip_primal_bound));}}).modelStatus;solution=Array.from(model.getSolution().colValue as ArrayLike<number>);try{bound=model.info.get('mip_dual_bound');}catch{}}finally{model.dispose();}
  progress('Common sparse profile status · capacity '+capacity+' · status '+status+' · lower bound '+bound+' · '+((performance.now()-started)/1000).toFixed(3)+' seconds');
  details.push({capacity,status,variables:variables.length,constraints:constraints.length,bound,seconds:(performance.now()-started)/1000});
  if(status===8)continue;
  const valid=solution.length===variables.length&&solution.every(Number.isFinite)&&variables.every((v,i)=>!v.integer||Math.abs(solution[i]-Math.round(solution[i]))<1e-6)&&constraints.every(r=>{const value=r.terms.reduce((n,[id,w])=>n+w*solution[id],0);return value>=r.lower-1e-7&&value<=r.upper+1e-7;});if(!valid)throw Error('Sparse profile search incomplete: '+status);
  const activeCount=groups.reduce((n,g)=>n+Math.round(solution[g.active]),0),supportOptimal=status===7||Number.isFinite(bound)&&activeCount-Math.ceil(bound!-.01)<1;if(!supportOptimal)throw Error('Sparse allele support not certified; '+activeCount+' groups, lower bound '+bound);
  const profile:{events:PathEdit[];fraction:number}[]=wt.flatMap(w=>solution[w.mass]>0?[{events:[],fraction:solution[w.mass]}]:[]);
  for(const g of groups){const fraction=g.origins.reduce((n,o)=>n+solution[o.mass],0);if(!(fraction>0))continue;const origin=g.origins.reduce((a,b)=>solution[a.bit]>=solution[b.bit]?a:b),word=g.bases.map(bs=>alphabet[bs.reduce((best,id,b)=>solution[id]>solution[bs[best]]?b:best,0)]).join(''),events=normalizeMoleculeEdits(ref.sequence,[{start:origin.start,end:origin.start+(g.size<0?-g.size:g.size===0?span:0),inserted:word}]);if(!materialize(ref.sequence,events))throw Error('Invalid sparse profile sequence');profile.push({events,fraction});}
  const sequences=profile.map(g=>materialize(ref.sequence,g.events)!),errors=ps.map((p,i)=>Math.sqrt([...'ACGT'].reduce((n,b,k)=>n+(sequences.reduce((v,s,j)=>v+profile[j].fraction*Number(s[p]===b),0)-target[4*i+k])**2,0)/4)),maximum=Math.max(...errors);if(maximum>tolerance+1e-7)throw Error('Sparse profile reconstruction mismatch '+maximum);return{options:[profile],scores:[{maximum,rmse:Math.sqrt(errors.reduce((n,e)=>n+e*e,0)/errors.length)}],details};
 }
}
