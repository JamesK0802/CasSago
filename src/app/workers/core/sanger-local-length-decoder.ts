import{Prepared,PathEdit,materialize}from'./sanger-path-common';import{normalizeMoleculeEdits}from'./sanger-edit-normalize';
interface Family {reference:number;size:number;fraction:number}
interface Line {b:number;c:number;code:number;base:number;start?:number}
interface Piece {lo:number;hi:number;b:number;c:number}
interface Cycle {active:number[];fixed:number[];y:number[];wt:number;no:{c:number;code:number};lines:Line[];pieces:Piece[]}
function piecesFor(lines:Line[],baseline:number,maximum:number):Piece[]{
 const sorted=[...lines].sort((a,b)=>b.b-a.b||a.c-b.c),hull:(Line&{start:number})[]=[];
 for(const line of sorted){if(hull.at(-1)?.b===line.b)continue;let start=-Infinity;while(hull.length){const last=hull.at(-1)!;start=(last.c-line.c)/(line.b-last.b);if(start>last.start)break;hull.pop();}hull.push({...line,start:hull.length?start:-Infinity});}
 const out:Piece[]=[];for(let i=0;i<hull.length;i++){const line=hull[i],d=line.c-baseline,disc=line.b*line.b-8*d;if(disc<=0)continue;const root=Math.sqrt(disc),right=(-line.b+root)/4,left=right>0?d/(2*right):(-line.b-root)/4,lo=Math.max(0,line.start,left),hi=Math.min(maximum,hull[i+1]?.start??Infinity,right);if(hi>lo)out.push({lo,hi,b:line.b,c:d});}return out;
}
/** Exact continuous one-substitution-weight minimization for fixed indel
 * origins. The per-cycle ALT assignment has a quadratic cost; a lower envelope
 * of lines avoids a percentage grid or an exponential inserted-string list. */
function bestSubstitution(cycles:Cycle[],radius:number,span:number,maximum:number){
 let best={gain:0,weight:0,start:0};
 for(let start=0;start<=2*radius;start++){const events:{x:number;a:number;b:number;c:number}[]=[];for(const p of cycles.slice(start,start+span))for(const q of p.pieces){events.push({x:q.lo,a:2,b:q.b,c:q.c},{x:q.hi,a:-2,b:-q.b,c:-q.c});}events.push({x:0,a:0,b:0,c:0},{x:maximum,a:0,b:0,c:0});events.sort((a,b)=>a.x-b.x);let a=0,b=0,c=0;
  for(let i=0;i<events.length;){const x=events[i].x;while(i<events.length&&events[i].x===x){a+=events[i].a;b+=events[i].b;c+=events[i].c;i++;}const end=events[i]?.x??maximum;if(end<x)continue;const weight=a>0?Math.max(x,Math.min(end,-b/(2*a))):x,gain=a*weight*weight+b*weight+c;if(gain<best.gain)best={gain,weight,start};}
 }
 return best;
}
/** Minimal local hypothesis: one allele per nonzero length and at most one
 * substitution allele. Other references' offset windows must not overlap.
 * This is an experimental unit-emission decoder; residual QC must reject a
 * sample needing a richer hypothesis rather than certify it as recovered. */
export function decodeLocalLengths(refs:Prepared[],families:Family[],positions:number[],observed:number[],radius=3,subSpan=10){
 if(refs.some(r=>r.cuts.length!==1))throw Error('Local length decoder currently supports one target per reference.');
 const index=new Map(positions.map((p,i)=>[p,i])),result:{reference:number;events:PathEdit[];fraction:number}[]=[],details:any[]=[];
 for(let h=0;h<refs.length;h++){
  const ref=refs[h],cut=ref.cuts[0],own=families.filter(f=>f.reference===h),zero=own.find(f=>!f.size)?.fraction??0,indels=own.filter(f=>f.size!==0&&f.fraction>0),extent=Math.max(subSpan,...indels.filter(f=>f.size>0).map(f=>f.size)),begin=cut-radius,end=cut+radius+extent;
  if(indels.length>6)throw Error('Local length decoder hypothesis exceeds six nonzero length classes.');
  const observations=Array.from({length:end-begin},(_,i)=>{const p=begin+i,k=index.get(p);if(k===undefined)throw Error('Local editing window is not covered by the fitted trace.');const y=observed.slice(4*k,4*k+4);for(const f of families)if(f.reference!==h&&f.fraction>0){const r=refs[f.reference],c=r.cuts[0];if(p>=c-radius&&p<c+radius+Math.max(subSpan,10))throw Error('Offset edit windows overlap.');const base=r.sequence[p<c-radius?p:p-f.size],b='ACGT'.indexOf(base);if(b>=0)y[b]-=f.fraction;}return y;});
  let best:{score:number;origins:number[];cycles:Cycle[];sub:ReturnType<typeof bestSubstitution>}|undefined,visited=0;
  const evaluate=(origins:number[])=>{
   const cycles:Cycle[]=observations.map((y,j)=>{const p=begin+j,wt='ACGT'.indexOf(ref.sequence[p]),fixed=[0,0,0,0],active:number[]=[];fixed[wt]+=zero;for(let k=0;k<indels.length;k++){const f=indels[k],start=origins[k];if(f.size>0&&p>=start&&p<start+f.size)active.push(k);else{const pos=p<start?p:p-f.size,b='ACGT'.indexOf(ref.sequence[pos]);if(b>=0)fixed[b]+=f.fraction;}}
    let no={c:Infinity,code:0};const lines:Line[]=[];for(let code=0;code<4**active.length;code++){const v=[...fixed];for(let k=0;k<active.length;k++)v[(code>>(2*k))&3]+=indels[active[k]].fraction;const residual=y.map((x,b)=>x-v[b]),c=residual.reduce((s,x)=>s+x*x,0);if(c<no.c)no={c,code};for(let b=0;b<4;b++)if(b!==wt)lines.push({b:-2*(residual[b]-residual[wt]),c,code,base:b});}return{active,fixed,y,wt,no,lines,pieces:piecesFor(lines,no.c,zero)};
   });
   const sub=bestSubstitution(cycles,radius,subSpan,zero),score=cycles.reduce((s,p)=>s+p.no.c,0)+sub.gain;visited++;if(!best||score<best.score)best={score,origins:[...origins],cycles,sub};
  };
  const enumerate=(origins:number[])=>{if(origins.length===indels.length){evaluate(origins);return;}for(let start=cut-radius;start<=cut+radius;start++)enumerate([...origins,start]);};enumerate([]);
  const chosen=best!,words=indels.map(f=>f.size>0?Array(f.size).fill('A'):[]),subStart=begin+chosen.sub.start,subWord=ref.sequence.slice(subStart,subStart+subSpan).split('');
  chosen.cycles.forEach((cycle,j)=>{let code=cycle.no.code,base=cycle.wt;if(j>=chosen.sub.start&&j<chosen.sub.start+subSpan&&chosen.sub.weight>0){let value=cycle.no.c;for(const line of cycle.lines){const v=line.c+line.b*chosen.sub.weight+2*chosen.sub.weight**2;if(v<value){value=v;code=line.code;base=line.base;}}}for(let k=0;k<cycle.active.length;k++){const id=cycle.active[k];words[id][begin+j-chosen.origins[id]]='ACGT'[(code>>(2*k))&3];}if(j>=chosen.sub.start&&j<chosen.sub.start+subSpan)subWord[j-chosen.sub.start]='ACGT'[base];});
  result.push({reference:h,events:[],fraction:zero-chosen.sub.weight});for(let k=0;k<indels.length;k++){const f=indels[k],start=chosen.origins[k],events=normalizeMoleculeEdits(ref.sequence,[{start,end:f.size<0?start-f.size:start,inserted:words[k].join('')}]);result.push({reference:h,events,fraction:f.fraction});}const subEvents=normalizeMoleculeEdits(ref.sequence,[{start:subStart,end:subStart+subSpan,inserted:subWord.join('')}]);if(chosen.sub.weight>0&&subEvents.length)result.push({reference:h,events:subEvents,fraction:chosen.sub.weight});details.push({reference:ref.name,score:chosen.score,origins:chosen.origins,substitutionFraction:chosen.sub.weight,visited});
 }
 return{groups:result,details};
}
