import type {MarginalTargetState,MarginalSession} from './sanger-marginal-engine';
import {validMarginalSolution,type LinearTerm,LP_INFINITY} from './sanger-marginal-lp';

/** A compact joint hypothesis in the existing measured model. All origins
 * and inserted bases are chosen together. This is a proposal, not a restriction
 * on the uncertainty model or the number of biological alleles. Infeasibility
 * or an incomplete search never removes any of the original feasible results. */
export function proposeJointFamilies(session:MarginalSession,target:MarginalTargetState,seconds=5){
 const start=performance.now(),p=target.program,upper=p.upper.slice(),cost=p.cost.map(()=>0),integrality=p.upper.map(()=>0),rows=p.rows.map(r=>({...r,terms:r.terms.map(t=>[...t]as LinearTerm)}));
 const variable=(bound=1,integer=1,objective=0)=>{const id=upper.length;upper.push(bound);integrality.push(integer);cost.push(objective);return id;},row=(terms:LinearTerm[],lo:number,hi=lo)=>{const merged=new Map<number,number>();for(const[i,w]of terms)merged.set(i,(merged.get(i)||0)+w);rows.push({terms:[...merged].filter(([,v])=>Math.abs(v)>1e-15),lower:lo,upper:hi});};
 const departures=p.cost.flatMap((v,i)=>v?[[i,v]as LinearTerm]:[]),minimum=departures.reduce((n,[i,v])=>n+v*target.point[i],0);row(departures,0,minimum+1e-7);
 for(const g of target.groups){const mass=g.mass.reduce((n,[i,v])=>n+v*target.point[i],0);row(g.mass,mass);if(mass<=1e-8)continue;
  const origins:number[]=[],insertions=new Map<number,{flow:number[];origins:number[]}>();for(const e of g.edges){const z=variable();origins.push(z);row([[e.mass,1],[z,-mass]],0);for(const flow of e.free.values()){let f=insertions.get(flow[0]);if(!f){f={flow,origins:[]};insertions.set(flow[0],f);}f.origins.push(z);}}
  for(const f of insertions.values()){const bits=f.flow.map(()=>variable());row([...bits.map(i=>[i,1]as LinearTerm),...f.origins.map(i=>[i,-1]as LinearTerm)],0);f.flow.forEach((v,b)=>row([[v,1],[bits[b],-mass]],0));}
  row(origins.map(i=>[i,1]),1);
  if(minimum+1e-7<mass){for(const f of g.flow?.values()||[])for(const i of f.changed)upper[i]=0;}
  else for(const forms of g.bases.values()){const bits=forms.map(()=>variable());row(bits.map(i=>[i,1]),1);forms.forEach((f,b)=>row([...f,[bits[b],-mass]],0));}
 }
 if(target.measurements){for(const m of target.measurements){const den=Math.max(1e-12,m.total.reduce((n,[i,w])=>n+w*target.point[i],0));for(let b=0;b<4;b++){const terms=[...m.forms[b],...m.total.map(([i,w])=>[i,-w*m.observed[b]]as LinearTerm)],e=variable(10,0,1/den);row([...terms,[e,-1]],-LP_INFINITY,0);row([...terms,[e,1]],0,LP_INFINITY);}}}else{const measurementRows=p.rows.slice(-(4*target.positions.length+1),-1);for(const r of measurementRows){const y=(r.lower+r.upper)/2,e=variable(1,0,1);row([...r.terms,[e,-1]],-LP_INFINITY,y);row([...r.terms,[e,1]],y,LP_INFINITY);}}
 const starts=[0],indices:number[]=[],values:number[]=[];for(const r of rows){for(const[i,v]of r.terms){indices.push(i);values.push(v);}starts.push(indices.length);}
 const model:any=session.optimizer.createModel!();let status=-1,x:number[]=[],gap:number|undefined;
 try{model.options.set({output_flag:false,time_limit:seconds,mip_rel_gap:1e-6,mip_feasibility_tolerance:1e-8,primal_feasibility_tolerance:1e-9});model.passModel({numCols:upper.length,numRows:rows.length,colCost:cost,colLower:upper.map(()=>0),colUpper:upper,integrality,rowLower:rows.map(r=>r.lower),rowUpper:rows.map(r=>r.upper),matrix:{format:'csr',numCols:upper.length,numRows:rows.length,starts,indices,values}});status=model.run().modelStatus;x=Array.from(model.getSolution().colValue as ArrayLike<number>);try{gap=model.info.get('mip_gap');}catch{}}
 finally{model.dispose();}
 const valid=x.length===upper.length&&validMarginalSolution({upper,cost,rows},x)&&integrality.every((v,i)=>!v||Math.abs(x[i]-Math.round(x[i]))<=1e-7);
 return{point:valid?x.slice(0,p.upper.length):undefined,seconds:(performance.now()-start)/1000,status,valid,optimal:status===7,gap,variables:upper.length,binaries:integrality.reduce((n,v)=>n+v,0)};
}
