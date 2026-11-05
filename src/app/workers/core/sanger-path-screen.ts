import type{JunctionInput,JunctionFit,JunctionModel,JunctionRead}from'./sanger-junction-mixture';
import{jointPathColumn,pathColumn,mapPathModel}from'./sanger-junction-paths';
import{controlKnownContextEmission}from'./sanger-control-response';
import{mixPhase,phaseNeighbor}from'./sanger-signal';
const bases='ACGT';
/** Exact directional scores for canonical two-event columns, using prefix/suffix
 * sums. Every requested size is scanned; full columns are retained only for the
 * working set. Memory grows with shift states, not their Cartesian product. */
export function screenPathCandidates(input:JunctionInput,current:JunctionFit,range:number[],make:(h:number,sizes:number[])=>JunctionModel,limit=8){
 const key=(m:JunctionModel)=>m.reference+':'+m.sizes!.join('/'),known=new Set(current.models.map(key));
 const gradients=input.positions.map((_,p)=>{const v=[0,0,0,0];current.columns.forEach((c,i)=>{if(current.weights[i])for(let b=0;b<4;b++)v[b]+=current.weights[i]*c[4*p+b];});const total=v.reduce((s,v)=>s+v,0),pi=v.map(v=>v/Math.max(1e-12,total)),residual=pi.map((v,b)=>input.observed[4*p+b]-v),mean=residual.reduce((s,v,b)=>s+v*pi[b],0);return residual.map(v=>(input.positionWeights?.[p]??1)*(v-mean)/Math.max(1e-12,total));});
 const reads:JunctionRead[]=input.reads||[{references:input.references,positions:input.positions,rowIds:input.positions,mappings:input.references.map(()=>({reverse:false,offset:0})),lag:input.lag,theta:input.theta}];
 const rowAt=new Map(input.positions.map((p,i)=>[p,i])),best:{m:JunctionModel;gain:number}[]=[];
 const retain=(m:JunctionModel,gain:number)=>{if(best.length<limit||gain>best.at(-1)!.gain){best.push({m,gain});best.sort((a,b)=>b.gain-a.gain);if(best.length>limit)best.pop();}};
 input.references.forEach((ref,h)=>{
  const axes=reads.map(read=>{
   const records=read.positions.flatMap((p,j)=>rowAt.has(read.rowIds[j])?[{p,row:rowAt.get(read.rowIds[j])!}]:[]).sort((a,b)=>a.p-b.p),positions=records.map(v=>v.p),firstCache=new Map<number,Float64Array>(),lastCache=new Map<number,Float64Array>();
   const cuts=read.references[h].cuts||[read.references[h].cut],first=cuts.indexOf(Math.min(...cuts)),last=1-first;
   const prefix=(column:number[])=>{const result=new Float64Array(positions.length+1);for(let i=0;i<positions.length;i++)result[i+1]=result[i]+gradients[records[i].row].reduce((s,g,b)=>s+g*column[4*i+b],0);return result;};
   const firstScore=(size:number)=>{let value=firstCache.get(size);if(!value){const sizes=cuts.map(()=>0);sizes[first]=size;value=prefix(pathColumn(read.references[h],mapPathModel(make(h,sizes),read.mappings[h]),positions,read.lag,read.theta));firstCache.set(size,value);}return value;};
   const lastScore=(size:number)=>{let value=lastCache.get(size);if(!value){const r=read.references[h],one=(p:number)=>{const q=p-size;return r.contextResponse?controlKnownContextEmission(r.sequence[q-1]||'',r.sequence[q]||'',r.sequence[q+1]||'',r.contextResponse,p):Array.from(bases,b=>Number(r.sequence[q]===b));};value=prefix(positions.flatMap(p=>read.lag?mixPhase(one(p),one(phaseNeighbor(p,read.lag)),read.lag):one(p)));lastCache.set(size,value);}return value;};
   const lower=(p:number)=>{let lo=0,hi=positions.length;while(lo<hi){const m=(lo+hi)>>1;if(positions[m]<p)lo=m+1;else hi=m;}return lo;};
   return{read,records,positions,first,last,firstScore,lastScore,lower};
  });
  const scan=(sizes:number[])=>{
   const m=make(h,sizes);if(known.has(key(m))||m.path?.some(e=>!e.origins.length))return;
   if((ref.cuts?.length||1)===1){const column=jointPathColumn(input,m);retain(m,gradients.reduce((s,g,p)=>s+g.reduce((s,v,b)=>s+v*column[4*p+b],0),0));return;}
   const events=m.path||[];if(events.length===2&&events[0].origins[0].start+Math.max(0,-events[0].size)>events[1].origins[0].start)return;
   let gain=0;
   for(const axis of axes){
    const first=axis.firstScore(sizes[axis.first]);if(!sizes[axis.last]){gain+=first.at(-1)!;continue;}
    const mapped=mapPathModel(m,axis.read.mappings[h]),event=mapped.path!.find(e=>e.target===axis.last+1)!,q=event.origins[0].start+sizes[axis.first],lo=axis.lower(q-2),hi=axis.lower(q+Math.max(0,event.size)+3),last=axis.lastScore(m.size);
    gain+=first[lo]+last.at(-1)!-last[hi];
    const column=pathColumn(axis.read.references[h],mapped,axis.positions.slice(lo,hi),axis.read.lag,axis.read.theta);
    for(let j=lo;j<hi;j++)for(let b=0;b<4;b++)gain+=gradients[axis.records[j].row][b]*column[4*(j-lo)+b];
   }
   retain(m,gain);
  };
  if((ref.cuts?.length||1)===1)range.forEach(a=>scan([a]));else for(const a of range)for(const b of range)scan([a,b]);
 });return best;
}
