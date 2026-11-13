import {tightenAuxiliaryBoxes} from './sanger-linear-boxes';
import {LP_INFINITY,OptimizationPending,marginalViolation,type MarginalProgram,type LinearTerm} from './sanger-marginal-lp';
interface FlatSolver {solve(terms?:LinearTerm[],maximize?:boolean):{x:number[];value:number};objectiveLowerBound?:()=>number|undefined;addRow(terms:LinearTerm[],lower:number,upper:number):void;setRowBounds(indices:number[],lower:number[],upper:number[]):void;rowDuals():number[];useDualForBoundChanges():void;dispose():void}
type Factory=(p:MarginalProgram)=>FlatSolver;
interface LocalRow {original?:number;link?:boolean;terms:LinearTerm[];master:LinearTerm[];lower:number;upper:number}
const compact=(terms:LinearTerm[])=>{const m=new Map<number,number>();for(const[i,v]of terms)m.set(i,(m.get(i)||0)+v);return[...m].filter(([,v])=>v!==0).sort((a,b)=>a[0]-b[0])as LinearTerm[];};
const dot=(terms:LinearTerm[],x:number[])=>terms.reduce((n,[i,w])=>n+w*x[i],0);
interface Link {column:number;minimum:number;terms:LinearTerm[]}
/** Exact LP decomposition. Local phase-I duals supply necessary supporting
 * cuts to the master. A result is returned only when its assembled witness
 * satisfies EVERY original row and attains the master objective bound.
 * No original variable, measurement, origin or constraint is discarded. */
export class FactorizedMarginalSolver {
 static observe?: (event:Record<string,unknown>)=>void;
 private columns:number[][];
 private masterColumns:number[];
 private masterMap:Map<number,number>;
 private owner:Map<number,{block:number;column:number}>;
 private base:MarginalProgram;
 private rows:LocalRow[][];
 private links=new Map<string,Link>();
 private elastic:ElasticBlock[];
 private cuts:{terms:LinearTerm[];lower:number;upper:number}[]=[];
 private cutKeys=new Set<string>();
 readonly stats={queries:0,iterations:0,cuts:0,blocks:0,seededCuts:0};
 constructor(readonly program:MarginalProgram,private flat:Factory){
  this.columns=(program.localBlocks||[]).map(v=>[...new Set(v)].filter(i=>program.upper[i]>0).sort((a,b)=>a-b)).filter(v=>v.length);this.owner=new Map();this.columns.forEach((v,block)=>v.forEach((id,column)=>{if(this.owner.has(id))throw Error('Overlapping local LP blocks.');this.owner.set(id,{block,column});}));
  this.masterColumns=program.upper.flatMap((_,i)=>program.upper[i]===0||this.owner.has(i)?[]:[i]);this.masterMap=new Map(this.masterColumns.map((id,j)=>[id,j]));this.base={upper:this.masterColumns.map(i=>program.upper[i]),cost:this.masterColumns.map(()=>0),rows:[]};this.rows=this.columns.map(()=>[]);
  if(program.quadratics)this.base.quadratics=program.quadratics.map(q=>({...q,terms:q.terms.map(row=>row.filter(([i,w])=>w!==0&&program.upper[i]>0).map(([i,w])=>{const id=this.masterMap.get(i);if(id===undefined)throw Error('A quadratic envelope crosses local LP blocks.');return[id,w]as LinearTerm;}))}));
  if(program.originPricing)this.base.originPricing={columns:program.originPricing.columns.filter(i=>this.masterMap.has(i)).map(i=>this.masterMap.get(i)!),seed:program.originPricing.seed.filter(i=>this.masterMap.has(i)).map(i=>this.masterMap.get(i)!)};
  for(const[original,row]of program.rows.entries()){const master:LinearTerm[]=[],byBlock=new Map<number,LinearTerm[]>();for(const[i,w]of row.terms){if(program.upper[i]===0)continue;const place=this.owner.get(i);if(place){const terms=byBlock.get(place.block)||[];terms.push([place.column,w]);byBlock.set(place.block,terms);}else master.push([this.masterMap.get(i)!,w]);}
   if(byBlock.size===1){const[block,terms]=[...byBlock][0];this.rows[block].push({original,terms:compact(terms),master:compact(master),lower:row.lower,upper:row.upper});}
   else{let constant=0;for(const[block,terms]of byBlock){const linked=this.link(block,terms,this.base,this.rows,this.links);master.push([linked.link.column,linked.factor]);constant+=linked.factor*linked.link.minimum;}this.base.rows.push({terms:compact(master),lower:row.lower<=-LP_INFINITY?-LP_INFINITY:row.lower-constant,upper:row.upper>=LP_INFINITY?LP_INFINITY:row.upper-constant});}
  }
  // The original default cost receives stable shared links. After a departure
  // budget is added these are the very same algebraic block contributions.
  const costLinks=new Map<number,{column:number;factor:number}>(),costBlocks=new Map<number,LinearTerm[]>();program.cost.forEach((w,i)=>{if(!w)return;const place=this.owner.get(i);if(place){const terms=costBlocks.get(place.block)||[];terms.push([place.column,w]);costBlocks.set(place.block,terms);}});for(const[block,terms]of costBlocks){const value=this.link(block,terms,this.base,this.rows,this.links);costLinks.set(block,{column:value.link.column,factor:value.factor});}
  // The hint certifies the original departure budget, even when a later
  // representative projection has a different optimization objective.
  if(program.dualHint?.shadow.length){costLinks.clear();const budgetBlocks=new Map<number,LinearTerm[]>();for(const[i,w]of program.rows[program.dualHint.budgetRow]?.terms||[]){const place=this.owner.get(i);if(place){const terms=budgetBlocks.get(place.block)||[];terms.push([place.column,w]);budgetBlocks.set(place.block,terms);}}for(const[b,terms]of budgetBlocks){const linked=this.link(b,terms,this.base,this.rows,this.links);costLinks.set(b,{column:linked.link.column,factor:linked.factor});}}
  tightenAuxiliaryBoxes(this.base,this.masterColumns.length);
  this.base.cost=this.base.upper.map(()=>0);this.elastic=this.rows.map((rows,b)=>new ElasticBlock(this.columns[b].map(i=>program.upper[i]),rows,flat));this.stats.blocks=this.columns.length;
  // Reuse the full original fit's local row duals as necessary cost cuts.
  // Box correction in seedCostCut makes every cut valid independently of
  // the hint's accuracy; no fitted allocation is fixed or removed.
  if(program.dualHint?.shadow.length)for(const[b,link]of costLinks){const cut=this.elastic[b].seedCostCut(program.dualHint.shadow,link.column,link.factor,this.base.upper),key=cut.terms.map(([i,w])=>i+':'+w.toPrecision(12)).join(',')+'<'+cut.upper.toPrecision(12);if(!this.cutKeys.has(key)){this.cutKeys.add(key);this.cuts.push({...cut,lower:-LP_INFINITY});this.stats.cuts++;this.stats.seededCuts++;}}
  FactorizedMarginalSolver.observe?.({phase:'constructed',...this.stats,masterVariables:this.base.upper.length,masterRows:this.base.rows.length});
 }
 private link(block:number,input:LinearTerm[],master:MarginalProgram,rows:LocalRow[][],registry:Map<string,Link>){
  const terms=compact(input),factor=terms.reduce((a,[,w])=>Math.abs(w)>Math.abs(a)?w:a,terms[0][1]),normalized=terms.map(([i,w])=>[i,w/factor]as LinearTerm),key=block+':'+normalized.map(([i,w])=>i+':'+w).join(','),existing=registry.get(key);if(existing)return{link:existing,factor};
  let minimum=0,maximum=0;for(const[i,w]of normalized){const u=this.program.upper[this.columns[block][i]];if(w<0)minimum+=w*u;else maximum+=w*u;}
  const column=master.upper.length,link={column,minimum,terms:normalized};master.upper.push(maximum-minimum);master.cost.push(0);registry.set(key,link);rows[block].push({link:true,terms:normalized,master:[[column,-1]],lower:minimum,upper:minimum});return{link,factor};
 }
 solve(terms?:LinearTerm[],maximize=false):{x:number[];value:number}{const steps=this.steps(terms,maximize);for(;;){const next=steps.next();if(next.done)return next.value;}}
 *steps(terms?:LinearTerm[],maximize=false):Generator<void,{x:number[];value:number},unknown>{
  this.stats.queries++;const requested=terms||this.program.cost.flatMap((v,i)=>v?[[i,v]as LinearTerm]:[]),cost=compact(requested.map(([i,w])=>[i,(maximize?-1:1)*w]as LinearTerm)),master:MarginalProgram={...this.base,upper:this.base.upper.slice(),cost:this.base.cost.map(()=>0),rows:[...this.base.rows,...this.cuts]},rows=this.rows.map(v=>v.slice()),registry=new Map(this.links),queryBlocks=new Set<number>(),byBlock=new Map<number,LinearTerm[]>();let constant=0;
  for(const[i,w]of cost){if(this.program.upper[i]===0)continue;const place=this.owner.get(i);if(place){const v=byBlock.get(place.block)||[];v.push([place.column,w]);byBlock.set(place.block,v);}else master.cost[this.masterMap.get(i)!]+=w;}
  for(const[block,local]of byBlock){const count=rows[block].length,linked=this.link(block,local,master,rows,registry);master.cost[linked.link.column]+=linked.factor;constant+=linked.factor*linked.link.minimum;if(rows[block].length!==count)queryBlocks.add(block);}
  const elastic=this.elastic.slice();for(const b of queryBlocks)elastic[b]=new ElasticBlock(this.columns[b].map(i=>this.program.upper[i]),rows[b],this.flat);
  tightenAuxiliaryBoxes(master,this.masterColumns.length);
  const originalMasterRows=master.rows.length,solver=this.flat(master),queryKeys=new Set<string>();solver.useDualForBoundChanges();
  try{for(;;){this.stats.iterations++;FactorizedMarginalSolver.observe?.({phase:'master',...this.stats});let answer:{x:number[];value:number};for(;;){try{answer=solver.solve();break;}catch(error){if(error instanceof OptimizationPending){yield;continue;}throw error;}}const point=this.program.upper.map(()=>0);this.masterColumns.forEach((i,j)=>point[i]=answer.x[j]);const pending:{terms:LinearTerm[];lower:number;upper:number;common:boolean;key:string}[]=[];
    for(let b=0;b<elastic.length;b++){let local:ReturnType<ElasticBlock['solve']>;for(;;){try{local=elastic[b].solve(answer.x,master.upper);break;}catch(error){if(error instanceof OptimizationPending){yield;continue;}throw error;}}this.columns[b].forEach((i,j)=>point[i]=local.x[j]);if(local.cut&&local.cut.violation>1e-10){const r=local.cut,key=r.terms.map(([i,w])=>i+':'+w.toPrecision(12)).join(',')+'<'+r.upper.toPrecision(12),common=r.terms.every(([i])=>i<this.base.upper.length);if(!queryKeys.has(key)&&!(common&&this.cutKeys.has(key))){queryKeys.add(key);pending.push({terms:r.terms,lower:-LP_INFINITY,upper:r.upper,common,key});}}if(b%16===15)yield;}
    const objective=dot(cost,point),gap=objective-(solver.objectiveLowerBound?.()??answer.value)-constant,violation=marginalViolation(this.program,point);FactorizedMarginalSolver.observe?.({phase:'verify',...this.stats,pendingCuts:pending.length,gap,violation});if(violation<=1e-7&&Math.abs(gap)<=1e-7){for(const r of master.rows.slice(originalMasterRows)){if(r.lower> -LP_INFINITY||r.terms.some(([i])=>i>=this.base.upper.length))continue;const key=r.terms.map(([i,w])=>i+':'+w.toPrecision(12)).join(',')+'<'+r.upper.toPrecision(12);if(!this.cutKeys.has(key)){this.cutKeys.add(key);this.cuts.push(r);this.stats.cuts++;}}return{x:point,value:dot(requested,point)};}
    if(!pending.length)throw Error('Factorized LP could not certify its original primal/dual witness: violation='+violation+'; gap='+gap);
    for(const r of pending){solver.addRow(r.terms,r.lower,r.upper);this.stats.cuts++;if(r.common){const cut={terms:r.terms,lower:r.lower,upper:r.upper};this.cuts.push(cut);this.cutKeys.add(r.key);}}
    yield;
   }}finally{solver.dispose();for(const block of elastic)block.releaseWorkspace();for(const b of queryBlocks)elastic[b].dispose();}
 }
 dispose(){for(const block of this.elastic)block.dispose();}
}
class ElasticBlock {
 private program?:MarginalProgram;
 private solver?:FlatSolver;
 private sides:{row:LocalRow;lower:boolean;soft:boolean}[]=[];
 private secondary?:ElasticBlock;
 private lastLower:number[]=[];
 private lastUpper:number[]=[];
 private cached?:{x:number[];dual:number[]};
 constructor(private upper:number[],private rows:LocalRow[],private flat:Factory,private allElastic=false){
  for(const row of rows)for(const lower of [true,false])if(lower?row.lower>-LP_INFINITY:row.upper<LP_INFINITY)this.sides.push({row,lower,soft:allElastic||!!row.link});
 }
 private materialize(lower:number[],upper:number[]){
  const p:MarginalProgram={upper:this.upper.slice(),cost:this.upper.map(()=>0),rows:[]};
  this.sides.forEach((side,k)=>{const terms=side.row.terms.slice();if(side.soft){const slack=p.upper.length;p.upper.push(LP_INFINITY);p.cost.push(1);terms.push([slack,side.lower?1:-1]);}p.rows.push({terms,lower:lower[k],upper:upper[k]});});
  return this.program=p;
 }
 solve(master:number[],masterUpper:number[]):{x:number[];cut:{terms:LinearTerm[];upper:number;violation:number}}{
  const lower=this.sides.map(s=>s.lower?s.row.lower-dot(s.row.master,master):-LP_INFINITY),upper=this.sides.map(s=>s.lower?LP_INFINITY:s.row.upper-dot(s.row.master,master)),changed=lower.flatMap((v,i)=>v!==this.lastLower[i]||upper[i]!==this.lastUpper[i]?[i]:[]);
  if(changed.length){if(this.solver)this.solver.setRowBounds(changed,changed.map(i=>lower[i]),changed.map(i=>upper[i]));this.cached=undefined;this.lastLower=lower;this.lastUpper=upper;}
  if(!this.cached){try{if(!this.solver){this.solver=this.flat(this.materialize(lower,upper));this.solver.useDualForBoundChanges();}const answer=this.solver.solve();this.cached={x:answer.x.slice(0,this.upper.length),dual:this.solver.rowDuals()};this.releaseWorkspace();}catch(error){if(error instanceof OptimizationPending)throw error;this.releaseWorkspace();if(!this.allElastic&&/did not complete: 8(?:;|$)/.test(String(error))){this.secondary??=new ElasticBlock(this.upper,this.rows,this.flat,true);return this.secondary.solve(master,masterUpper);}throw error;}}
  const cut=this.dualCut(this.cached.dual,masterUpper);return{x:this.cached.x,cut:{...cut,violation:dot(cut.terms,master)-cut.upper}};
 }
 seedCostCut(shadow:number[],column:number,factor:number,masterUpper:number[]){
  const scale=1/Math.max(1,Math.abs(factor)),dual=this.sides.map(side=>{const y=side.row.original!==undefined?(shadow[side.row.original]||0)*scale:side.row.link&&side.row.master.length===1&&side.row.master[0][0]===column?-factor*scale:0;return side.lower?Math.max(0,y):Math.min(0,y);});return this.dualCut(dual,masterUpper);
 }
 private dualCut(dual:number[],masterUpper:number[]){
  const reduced=this.upper.map(()=>0),gradient:LinearTerm[]=[];let intercept=0,magnitude=0;
  for(let k=0;k<this.sides.length;k++){const{row,lower,soft}=this.sides[k],y=lower?Math.max(0,Math.min(soft?1:Infinity,dual[k])):Math.max(soft?-1:-Infinity,Math.min(0,dual[k]));if(!Number.isFinite(y))throw Error('Invalid local LP dual certificate.');const v=y*(lower?row.lower:row.upper);intercept+=v;magnitude+=Math.abs(v);for(const[i,w]of row.terms)reduced[i]-=y*w;for(const[i,w]of row.master)gradient.push([i,-y*w]);}
  for(let i=0;i<reduced.length;i++)if(reduced[i]<0){const v=reduced[i]*this.upper[i];intercept+=v;magnitude+=Math.abs(v);}
  const terms=compact(gradient),guard=64*Number.EPSILON*(1+magnitude+terms.reduce((n,[i,w])=>n+Math.abs(w)*masterUpper[i],0)),upperBound=-intercept+guard;
  return{terms,upper:upperBound};
 }
 /** Keep verified primal/dual vectors, not every native factorization. A
  * pending solve retains its one workspace until resumed or cancelled. */
 releaseWorkspace(){this.solver?.dispose();this.solver=undefined;this.program=undefined;this.secondary?.releaseWorkspace();}
 dispose(){this.releaseWorkspace();this.secondary?.dispose();}
}
