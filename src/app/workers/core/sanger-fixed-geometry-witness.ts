import {tryMarginalWitness} from './sanger-witness-proposal';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
import {LP_INFINITY,validMarginalSolution,type LinearTerm,type MarginalProgram} from './sanger-marginal-lp';
/** Find a witness on an existing geometry, without claiming an optimum for
 * that restricted problem. Only a witness checked in the complete original
 * model is returned. Failure leaves the unrestricted endpoint solve intact. */
export function fixedGeometryWitness(program:MarginalProgram,anchor:number[],freeColumns:number[],objective:LinearTerm[],value:number,optimizer:SangerLinearOptimizer):number[]|undefined{
 const columns=[...new Set(freeColumns)].sort((a,b)=>a-b),index=new Map(columns.map((i,j)=>[i,j]));
 if((program.quadratics||[]).some(q=>q.terms.some(row=>row.some(([i,w])=>w!==0&&index.has(i)))))return;
 const project=(terms:LinearTerm[])=>{const local:LinearTerm[]=[];let constant=0;for(const[i,w]of terms){const j=index.get(i);if(j===undefined)constant+=w*anchor[i];else local.push([j,w]);}return{terms:local,constant};},p:MarginalProgram={upper:columns.map(i=>program.upper[i]),cost:columns.map(()=>0),rows:[]};
 for(const row of program.rows){const r=project(row.terms),lower=row.lower<=-LP_INFINITY?-LP_INFINITY:row.lower-r.constant,upper=row.upper>=LP_INFINITY?LP_INFINITY:row.upper-r.constant;if(!r.terms.length){if(lower>1e-7||upper< -1e-7)return;continue;}p.rows.push({terms:r.terms,lower,upper});}
 p.localBlocks=program.localBlocks?.map(block=>block.filter(i=>index.has(i)).map(i=>index.get(i)!)).filter(block=>block.length);
 const goal=project(objective);p.rows.push({terms:goal.terms,lower:value-goal.constant,upper:value-goal.constant});
 const answer=tryMarginalWitness(p,optimizer,10);if(!answer)return;const x=anchor.slice();columns.forEach((i,j)=>x[i]=answer.x[j]);if(validMarginalSolution(program,x)&&Math.abs(objective.reduce((n,[i,w])=>n+w*x[i],0)-value)<=1e-7)return x;
 return;
}
