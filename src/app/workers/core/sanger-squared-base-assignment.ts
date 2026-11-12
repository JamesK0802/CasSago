import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import {normalizeMoleculeEdits} from './sanger-edit-normalize';
import type {PathEdit} from './sanger-path-common';
interface Group{events:PathEdit[];size:number;sequence:string;weights:Map<number,number>;editable?:boolean}
interface Observation{terms:[number,string][];lower:number;upper:number;error:string;points:Set<number>}
/** Legal origins and bases under the SAME squared fluorescence objective as
 * the complete mixture fit. Supporting planes give valid global lower bounds;
 * adding the integer incumbent's exact prediction closes the convex envelope.
 * The same supporting-plane algorithm is used at every group count.
 * Every returned incumbent satisfies the complete discrete edit geometry. */
export function assignSquaredSequenceBases(sequence:string,cut:number,groups:Group[],positions:number[],observed:number[][],optimizer:SangerLinearOptimizer,seconds=10,intervals?:{lower:number[][];upper:number[][];scale?:number},distinctGroups=true,maxEditedGroups=Infinity){
 const began=performance.now(),bases='ACGT',q=Array.from({length:6+Math.max(10,...groups.map(g=>g.size))},(_,i)=>cut-3+i),index=new Map(q.map((p,i)=>[p,i])),constraints:string[]=[],binary:string[]=[],checks:{terms:[number,string][];op:string;rhs:number}[]=[],rows:Observation[]=[];let serial=0,constantObjective=0;const lossScale=1e6;
 const term=(coefficient:number,name:string)=>(coefficient<0?'- ':'+ ')+Math.abs(coefficient)+' '+name,expression=(terms:[number,string][])=>terms.filter(([w])=>w!==0).map(([w,n])=>term(w,n)).join(' '),add=(terms:[number,string][],op:string,rhs:number)=>{checks.push({terms,op,rhs});constraints.push(' r'+serial+++': '+expression(terms)+' '+op+' '+rhs);},bit=(g:number,k:number,b:number)=>`b${g}_${k}_${b}`,origin=(g:number,i:number)=>`z${g}_${i}`;
 groups.forEach((g,j)=>{if(!(g.editable??g.events.length>0))return;const origins=Array.from({length:7},(_,i)=>origin(j,i));binary.push(...origins);add(origins.map(n=>[1,n]),'=',1);q.forEach((p,k)=>{const bits=Array.from({length:4},(_,b)=>bit(j,k,b));binary.push(...bits);add(bits.map(n=>[1,n]),'=',1);for(let i=0;i<7;i++){const start=cut-3+i;if(start<=p&&p<start+(g.size>0?g.size:g.size===0?10:0))continue;const at=p<start?p:p-g.size,b=bases.indexOf(sequence[at]);if(b<0)throw Error('Incomplete reference coverage in squared base assignment.');add([[1,origins[i]],[-1,bits[b]]],'<=',0);}});});
 if(distinctGroups)for(let a=0;a<groups.length;a++)for(let b=a+1;b<groups.length;b++)if(groups[a].size===groups[b].size){
  const ae=groups[a].editable??groups[a].events.length>0,be=groups[b].editable??groups[b].events.length>0;
  if(!ae&&!be){if(groups[a].sequence===groups[b].sequence)return null;continue;}
  if(ae&&be){const differences:[number,string][]=[];for(let k=0;k<q.length;k++)for(let c=0;c<4;c++){const x=bit(a,k,c),y=bit(b,k,c),d=`d${a}_${b}_${k}_${c}`;differences.push([1,d]);add([[1,d],[-1,x],[1,y]],'>=',0);add([[1,d],[1,x],[-1,y]],'>=',0);add([[1,d],[-1,x],[-1,y]],'<=',0);add([[1,d],[1,x],[1,y]],'<=',2);}add(differences,'>=',2);}
  else{const editable=ae?a:b,fixed=ae?b:a,terms:[number,string][]=q.map((p,k)=>[1,bit(editable,k,bases.indexOf(groups[fixed].sequence[p]))]);add(terms,'<=',q.length-1);}
 }
 if(Number.isFinite(maxEditedGroups)){
  const wtTerms:[number,string][]=[],fixedWt=groups.filter(g=>!(g.editable??g.events.length>0)&&!g.events.length).length;
  groups.forEach((g,j)=>{if(g.size!==0||!(g.editable??g.events.length>0))return;const wt='wt'+j;binary.push(wt);wtTerms.push([1,wt]);const matches=q.map((p,k)=>bit(j,k,bases.indexOf(sequence[p])));for(const b of matches)add([[1,wt],[-1,b]],'<=',0);add([[1,wt],...matches.map(b=>[-1,b]as[number,string])],'>=',1-q.length);});
  const needed=groups.length-maxEditedGroups-fixedWt;if(needed>wtTerms.length)return null;if(needed>0)add(wtTerms,'>=',needed);
 }
 const distance=(value:number,lower:number,upper:number)=>value<lower?value-lower:value>upper?value-upper:0;
 positions.forEach((p,i)=>{for(let b=0;b<4;b++){const sums=new Map<string,number>();let fixed=0;groups.forEach((g,j)=>{for(const[shift,w]of g.weights){if(!w)continue;const at=p-shift,k=index.get(at);if((g.editable??g.events.length>0)&&k!==undefined){const name=bit(j,k,b);sums.set(name,(sums.get(name)||0)+w);}else if((g.events.length?sequence[at<q[0]?at:at-g.size]:g.sequence[at])===bases[b])fixed+=w;}});const terms=[...sums].map(([n,w])=>[w,n]as[number,string]),lower=(intervals?.lower[i][b]??observed[i][b])-fixed,upper=(intervals?.upper[i][b]??observed[i][b])-fixed;
  if(!terms.length){constantObjective+=distance(0,lower,upper)**2;continue;}
  const total=terms.reduce((n,[w])=>n+w,0),clamp=(x:number)=>Math.max(0,Math.min(total,x)),points=new Set<number>([0,total,clamp(lower),clamp(upper)]);
  for(const f of[-.125,-.03125,.03125,.125]){points.add(clamp(lower+total*f));points.add(clamp(upper+total*f));}
  rows.push({terms,lower,upper,error:'e'+rows.length,points});}});

 if(!binary.length||!rows.length)return null;
 const readEvents=(read:(name:string)=>number)=>groups.map((g,j)=>{if(!(g.editable??g.events.length>0))return[];let choice=0;for(let i=1;i<7;i++)if(read(origin(j,i))>read(origin(j,choice)))choice=i;const start=cut-3+choice,length=g.size>0?g.size:g.size===0?10:0,word=Array.from({length},(_,i)=>{const k=index.get(start+i)!;let dye=0;for(let b=1;b<4;b++)if(read(bit(j,k,b))>read(bit(j,k,dye)))dye=b;return bases[dye];}).join('');return normalizeMoleculeEdits(sequence,[{start,end:start+(g.size<0?-g.size:g.size===0?10:0),inserted:word}]);});
 let best:{events:PathEdit[][];objective:number;lowerBound:number;optimal:boolean;rounds:number}|null=null,lowerBound=constantObjective;
 for(let round=0;round<64;round++){
  const remaining=seconds-(performance.now()-began)/1000;if(remaining<=0)break;const planes:string[]=[];let id=0;
  for(const row of rows)for(const p of row.points){const error=distance(p,row.lower,row.upper),slope=lossScale*2*error;planes.push(' p'+id+++': '+expression([...row.terms.map(([w,n])=>[slope*w,n]as[number,string]),[-1,row.error]])+' <= '+(lossScale*(2*error*p-error*error)));}
  const program=['Minimize',' obj: '+rows.map(r=>r.error).join(' + '),'Subject To',...constraints,...planes,'Binary',...binary,'End'].join('\n'),answer=optimizer.solve(program,{output_flag:false,time_limit:remaining,mip_abs_gap:1e-12,mip_rel_gap:1e-8,mip_feasibility_tolerance:1e-9,primal_feasibility_tolerance:1e-9})as ReturnType<SangerLinearOptimizer['solve']>&{ObjectiveValue?:number},read=(name:string)=>answer.Columns?.[name]?.Primal||0;
  if(!answer.Columns||binary.some(n=>!Number.isFinite(answer.Columns[n]?.Primal)||Math.abs(read(n)-Math.round(read(n)))>1e-7))break;
  if(checks.some(c=>{const v=c.terms.reduce((n,[w,k])=>n+w*(binary.includes(k)?Math.round(read(k)):read(k)),0);return(c.op==='='?Math.abs(v-c.rhs):c.op==='<='?v-c.rhs:c.rhs-v)>1e-7;}))break;
  const predictions=rows.map(r=>r.terms.reduce((n,[w,k])=>n+w*Math.round(read(k)),0)),objective=constantObjective+rows.reduce((n,r,i)=>n+distance(predictions[i],r.lower,r.upper)**2,0);if(answer.Status==='Optimal'&&Number.isFinite(answer.ObjectiveValue))lowerBound=Math.max(lowerBound,answer.ObjectiveValue!/lossScale+constantObjective);
  const optimal=answer.Status==='Optimal'&&objective-lowerBound<=Math.max(1e-11,objective*1e-8);if(!best||objective<best.objective)best={events:readEvents(n=>Math.round(read(n))),objective,lowerBound,optimal,rounds:round+1};if(optimal){best!.optimal=true;best!.lowerBound=lowerBound;break;}
  let added=0;rows.forEach((r,i)=>{if(!r.points.has(predictions[i])){r.points.add(predictions[i]);added++;}});if(!added)break;
 }
 return best;
}
