import {SangerAnalysisResult} from '../models/sanger-analysis';
import {exactDisplayView} from './sanger-exact-display';
import {SangerAnnotationRow} from './sanger-annotation';
export const sangerGroupKey=(reference:string,target:number,id:string)=>JSON.stringify([reference,target,id]);
export function sangerCurationView(result:SangerAnalysisResult,index:number,target:number,excluded:readonly string[],bounds:number|{left:number;right:number}=45){
 const annotation=result.annotations?.[index];if(!annotation)return null;
 const view=exactDisplayView(annotation,target,bounds,result.exactAlternatives?.find(s=>s.reference===annotation.name));
 const retained=view.rows.filter(row=>!excluded.includes(sangerGroupKey(annotation.name,target,row.id))),mass=retained.reduce((sum,row)=>sum+row.fraction,0);
 const rows=retained.map(row=>({...row,fraction:mass?row.fraction/mass:0}));
 const sum=(test:(row:SangerAnnotationRow)=>boolean)=>annotation.estimable&&mass?rows.reduce((n,row)=>n+(test(row)?row.fraction:0),0):null;
 const substitution=(row:SangerAnnotationRow)=>row.edited&&!!row.representative?.events.length&&row.representative.events.every(e=>e.inserted===e.end-e.start);
 return {view,rows,mass,wt:sum(row=>!row.edited),substitution:sum(substitution),inFrame:sum(row=>row.edited&&!substitution(row)&&row.size%3===0),outOfFrame:sum(row=>row.edited&&row.size%3!==0)};
}
