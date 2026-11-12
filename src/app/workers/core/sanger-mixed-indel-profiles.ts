import {Prepared,PathEdit,materialize} from './sanger-path-common';
import {normalizeMoleculeEdits} from './sanger-edit-normalize';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
/** Jointly infer frequencies and bases. Deletions have a finite legal geometry
 * dictionary; insertion/replacement words use binary bases and exact binary x
 * continuous linearization. Model order grows from signal, never a fixed K. */
export function mixedIndelProfiles(ref:Prepared,families:{size:number;fraction:number}[],positions:number[],observed:number[],noise:number,budget:number,progress:(s:string)=>void,_displayThresholdPercent:number,optimizer?:SangerLinearOptimizer,proposalOnly=false){
 if(!optimizer)throw Error('Initialize the discrete optimizer before profile inference.');
 const bases='ACGT',share=families.reduce((n,f)=>n+f.fraction,0),y=observed.map(v=>v/share),minimum=Math.max(1e-8,2*noise/share),epsilon=Math.max(1e-5,4*noise/share),cut=ref.cuts[0],q=Array.from({length:16},(_,i)=>cut-3+i),qIndex=new Map(q.map((p,i)=>[p,i]));
 const variable=families.filter(f=>f.fraction>0&&f.size>=0),initial=variable.map(f=>f.size===0?0:1),queue=[initial],visited=new Set<string>();let trials=0;const started=Date.now(),proposals:{groups:{events:PathEdit[];fraction:number}[];score:{maximum:number;rmse:number};details:any}[]=[];
 interface Group {size:number;name:string;fixed?:PathEdit[];weight:number;fixedWeight?:number;events?:PathEdit[]}
 while(queue.length&&trials<budget){
  if(proposalOnly&&trials>0&&Date.now()-started>30000)break;
  const counts=queue.shift()!,key=counts.join('/');if(visited.has(key))continue;visited.add(key);trials++;
  progress('Joint sequence/frequency model · '+trials+' · variable groups '+counts.join('/'));
  const groups:Group[]=[],constraints:string[]=[],bounds:string[]=[],binary:string[]=[],cost:string[]=[],weightsBySize=new Map<number,string[]>();let row=0;
  const term=(coefficient:number,name:string)=>(coefficient<0?' - ':' + ')+Math.abs(coefficient)+' '+name;
  const add=(terms:[number,string][],op:string,rhs:number)=>{constraints.push(' r'+row+++': '+(terms.filter(([c])=>c!==0).map(([c,n])=>term(c,n)).join('')||'0 w0')+' '+op+' '+rhs);};
  const register=(size:number,fixed?:PathEdit[])=>{const name='w'+groups.length,g={size,name,fixed,weight:0};groups.push(g);weightsBySize.set(size,[...(weightsBySize.get(size)||[]),name]);bounds.push(' 0 <= '+name+' <= 1');return g;};
  // Legal deletions have unknown continuous weights, not frequency atoms.
  for(const f of families.filter(f=>f.size<0&&f.fraction>0))for(let start=cut-3;start<=cut+3;start++){
   if(start<0||start-f.size>ref.sequence.length)continue;const g=register(f.size,[{start,end:start-f.size,inserted:''}]),z='d'+groups.length;binary.push(z);add([[1,g.name],[-1,z]],'<=',0);if(minimum>0)add([[1,g.name],[-(proposalOnly?Math.min(minimum,f.fraction/share):minimum),z]],'>=',0);cost.push('0.000000001 '+z);
  }
  if(families.some(f=>f.size===0))register(0,[]);
  variable.forEach((f,i)=>{for(let n=0;n<counts[i];n++){const g=register(f.size) as Group;if(f.size>0&&counts[i]===1){g.fixedWeight=f.fraction/share;add([[1,g.name]],"=",g.fixedWeight);}}});
  // Symmetry breaking: equal-length unknown alleles ordered by frequency.
  for(const f of variable){const own=groups.filter(g=>g.size===f.size&&!g.fixed);for(let i=1;i<own.length;i++)add([[1,own[i-1].name],[-1,own[i].name]],'>=',0);}
  for(const f of families){const names=weightsBySize.get(f.size)||[];if(!names.length)continue;add(names.map(n=>[1,n]),'>=',Math.max(0,f.fraction/share-epsilon));add(names.map(n=>[1,n]),'<=',f.fraction/share+epsilon);}
  add(groups.map(g=>[1,g.name]),'=',1);
  const unknown=groups.flatMap((g,j)=>!g.fixed?[j]:[]);
  for(const j of unknown){const g=groups[j],origins=Array.from({length:7},(_,i)=>'z'+j+'_'+i);binary.push(...origins);add(origins.map(n=>[1,n]),'=',1);if(minimum>0)add([[1,g.name]],'>=',proposalOnly?Math.min(minimum,(families.find(f=>f.size===g.size)?.fraction||0)/share):minimum);
   for(let k=0;k<q.length;k++){
    const bs=Array.from({length:4},(_,b)=>'b'+j+'_'+k+'_'+b);binary.push(...bs);add(bs.map(n=>[1,n]),'=',1);
    for(let i=0;i<7;i++){const start=cut-3+i;if(q[k]>=start&&q[k]<start+(g.size>0?g.size:10))continue;const b=bases.indexOf(ref.sequence[q[k]<start?q[k]:q[k]-g.size]);if(b<0)throw Error('Incomplete profile reference coverage.');add([[1,origins[i]],[-1,bs[b]]],'<=',0);}
    if(g.fixedWeight===undefined)for(let b=0;b<4;b++){const u='u'+j+'_'+k+'_'+b;bounds.push(' 0 <= '+u+' <= 1');add([[1,u],[-1,g.name]],'<=',0);add([[1,u],[-1,bs[b]]],'<=',0);add([[1,u],[-1,g.name],[-1,bs[b]]],'>=',-1);}
   }
  }
  const fixedSequences=groups.map(g=>g.fixed?materialize(ref.sequence,g.fixed):null);
  for(let i=0;i<positions.length;i++)for(let b=0;b<4;b++){
   const terms:[number,string][]=[],p=positions[i];for(let j=0;j<groups.length;j++){const g=groups[j];if(g.fixed){if(fixedSequences[j]?.[p]===bases[b])terms.push([1,g.name]);}else{const k=qIndex.get(p);if(k!==undefined)terms.push(g.fixedWeight===undefined?[1,'u'+j+'_'+k+'_'+b]:[g.fixedWeight,'b'+j+'_'+k+'_'+b]);else if(ref.sequence[p<q[0]?p:p-g.size]===bases[b])terms.push([1,g.name]);}}
   const e='e'+i+'_'+b;bounds.push(' 0 <= '+e+' <= '+(proposalOnly?1:epsilon));cost.push(e);add([...terms,[-1,e]],'<=',y[4*i+b]);add([...terms,[1,e]],'>=',y[4*i+b]);
  }
  const program=['Minimize',' obj: '+cost.join(' + '),'Subject To',...constraints,'Bounds',...bounds,'Binary',...binary,'End'].join('\n');
  const answer=optimizer.solve(program,{output_flag:false,time_limit:proposalOnly?2:10,mip_rel_gap:1e-7,mip_abs_gap:1e-9,mip_feasibility_tolerance:1e-9,primal_feasibility_tolerance:1e-9});
  const read=(n:string)=>answer.Columns[n]?.Primal??NaN;
  let valid=groups.every(g=>Number.isFinite(read(g.name))&&read(g.name)>=-1e-8)&&binary.every(n=>Number.isFinite(read(n))&&Math.abs(read(n)-Math.round(read(n)))<1e-7);
  if(valid){
   for(let j=0;j<groups.length;j++){const g=groups[j];g.weight=Math.max(0,read(g.name));if(g.fixed)g.events=g.fixed;else{let origin=0;for(let i=1;i<7;i++)if(read('z'+j+'_'+i)>read('z'+j+'_'+origin))origin=i;const start=cut-3+origin,length=g.size>0?g.size:10,word=Array.from({length},(_,i)=>{const k=qIndex.get(start+i)!;let b=0;for(let c=1;c<4;c++)if(read('b'+j+'_'+k+'_'+c)>read('b'+j+'_'+k+'_'+b))b=c;return bases[b];}).join('');g.events=normalizeMoleculeEdits(ref.sequence,[{start,end:start+(g.size===0?10:0),inserted:word}]);}}
   const selected=groups.filter(g=>g.weight>1e-8),sequences=selected.map(g=>materialize(ref.sequence,g.events!)!);let maximum=0,error=0;
   for(let i=0;i<positions.length;i++){let e=0;for(let b=0;b<4;b++){const d=selected.reduce((n,g,j)=>n+g.weight*Number(sequences[j]?.[positions[i]]===bases[b]),0)-y[4*i+b];if(!proposalOnly&&Math.abs(d)>epsilon+1e-7)valid=false;e+=d*d;}error+=e;maximum=Math.max(maximum,share*Math.sqrt(e/4));}
   if(Math.abs(selected.reduce((n,g)=>n+g.weight,0)-1)>1e-7||!proposalOnly&&selected.some(g=>g.events!.length&&g.weight+1e-8<minimum))valid=false;
   if(valid){const groups=selected.map(g=>({events:g.events!,fraction:g.weight*share})),score={maximum,rmse:share*Math.sqrt(error/y.length)},details={method:'joint-frequency-and-sequence',trials,counts,status:answer.Status,proposalOnly,exhaustive:false};if(maximum<=Math.max(1e-5,4*noise)+1e-7)return{options:[groups],scores:[score],details:[details]};if(proposalOnly)proposals.push({groups,score,details});}
  }
  // Increase only data-supported edited slots. Threshold gives a physical mass
  // bound; at zero threshold, the computational trial budget alone bounds work.
  for(let i=0;i<counts.length;i++){const f=variable[i];if(minimum>0&&(counts[i]+1)*minimum>f.fraction/share+epsilon)continue;const next=[...counts];next[i]++;if(!visited.has(next.join('/')))queue.push(next);}
 }
 if(proposalOnly&&proposals.length){proposals.sort((a,b)=>a.score.rmse-b.score.rmse);const selected=proposals.slice(0,4);return{options:selected.map(p=>p.groups),scores:selected.map(p=>p.score),details:selected.map(p=>p.details)};}
 throw Error('Joint profile discovery exhausted its computational budget without a signal-valid solution.');
}
