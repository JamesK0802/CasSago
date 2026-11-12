import {alleleComplexity} from './sanger-allele-complexity';
import type {Prepared,PathEdit} from './sanger-path-common';
interface Geometry {start:number;removed:number;inserted:number;delta:number;complexity?:number}
/** Price whole molecules against a development residual. All legal insertion
 * words are represented implicitly: the best base is selected at each inserted
 * position. No recipe, truth, file label or reporting cutoff enters pricing. */
export function priceMolecules(refs:Prepared[],maxDeletion:number,maxInsertion:number,positions:number[],residual:number[],folds:number[],means:number[][],control:number[][],lag:number,penalty:number){
 const alphabet='ACGT',N=Math.max(...positions)+maxDeletion*2+maxInsertion*2+5,cost=Array.from({length:N},()=>[0,0,0,0]),count=folds.filter(f=>f!==2).length;
 for(let k=0;k<residual.length;k++)if(folds[k]!==2){const p=positions[Math.floor(k/4)],b=k%4,v=-residual[k]/count;cost[p][b]+=(1-Math.abs(lag))*v;const q=p+(lag<0?1:-1);if(q>=0&&q<N)cost[q][b]+=Math.abs(lag)*v;}
 const baseCost=cost.map(row=>means.map(emission=>row.reduce((n,v,b)=>n+v*emission[b],0))),best=baseCost.map(row=>row.indexOf(Math.min(...row))),wordPrefix=new Float64Array(N+1);for(let p=0;p<N;p++)wordPrefix[p+1]=wordPrefix[p]+baseCost[p][best[p]];
 const result:{reference:number;events:PathEdit[];cost:number}[]=[],clamp=(p:number)=>Math.max(0,Math.min(N,p));
 refs.forEach((ref,h)=>{
  const prefix=new Map<number,Float64Array>();for(let shift=-maxDeletion*ref.cuts.length;shift<=maxInsertion*ref.cuts.length;shift++){const sums=new Float64Array(N+1);for(let p=0;p<N;p++){const donor=p-shift,b=alphabet.indexOf(ref.sequence[donor]),emission=refs.length===1?control[donor]:undefined;sums[p+1]=sums[p]+(emission?cost[p].reduce((n,v,k)=>n+v*emission[k],0):b>=0?baseCost[p][b]:0);}prefix.set(shift,sums);}
  const geometries=ref.cuts.map(cut=>{const rows:Geometry[]=[{start:cut,removed:0,inserted:0,delta:0}];for(let size=1;size<=maxDeletion;size++)for(let start=cut-size-3;start<=cut+3;start++)if(start>=0&&start+size<=ref.sequence.length)rows.push({start,removed:size,inserted:0,delta:-size});for(let size=1;size<=maxInsertion;size++)for(let start=cut-3;start<=cut+3;start++)if(start>=0)rows.push({start,removed:0,inserted:size,delta:size});return rows;});
  for(const rows of geometries)for(const g of rows)g.complexity=alleleComplexity([{start:g.start,end:g.start+g.removed,inserted:'N'.repeat(g.inserted)}],ref.cuts);
  const options=new Map<string,{geometries:Geometry[];cost:number}>(),consider=(gs:Geometry[])=>{let shift=0,from=0,value=0;for(const g of gs){const a=clamp(g.start+shift),b=clamp(g.start+shift+g.inserted),p=prefix.get(shift)!;value+=p[a]-p[from]+wordPrefix[b]-wordPrefix[a];if(g.removed||g.inserted)value+=penalty*g.complexity!;shift+=g.delta;from=b;}value+=prefix.get(shift)![N]-prefix.get(shift)![from];const key=gs.map(g=>g.delta).join('/'),old=options.get(key);if(!old||value<old.cost)options.set(key,{geometries:gs,cost:value});};
  if(geometries.length===1)for(const a of geometries[0])consider([a]);else for(const a of geometries[0])for(const b of geometries[1])if(a.start+a.removed<=b.start)consider([a,b]);
  if(ref.cuts.length===2)for(let start=Math.max(0,ref.cuts[0]-maxDeletion);start<=ref.cuts[0]+3;start++)for(let end=ref.cuts[1]-3;end<=Math.min(ref.sequence.length,ref.cuts[1]+maxDeletion);end++)if(end>start&&end-start<=maxDeletion){const g:Geometry={start,removed:end-start,inserted:0,delta:start-end};g.complexity=alleleComplexity([{start,end,inserted:''}],ref.cuts);consider([g]);}
  for(const option of [...options.values()].sort((a,b)=>a.cost-b.cost).slice(0,32)){let shift=0;const events:PathEdit[]=[];for(const g of option.geometries){if(g.removed||g.inserted)events.push({start:g.start,end:g.start+g.removed,inserted:Array.from({length:g.inserted},(_,i)=>alphabet[best[clamp(g.start+shift+i)]??0]).join('')});shift+=g.delta;}result.push({reference:h,events,cost:option.cost});}
 });
 return result;
}
