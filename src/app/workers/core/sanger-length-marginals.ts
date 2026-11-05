import {Prepared,materialize}from'./sanger-path-common';import{SangerFitCache}from'./sanger-fit-cache';
/** Flanking observations identify net-length mass without committing to edit
 * placement or inserted strings. Local edit windows (including emission
 * neighbours) are excluded; no genotype recipe is part of this input. */
export function fitLengthMarginals(refs:Prepared[],positions:number[],observed:number[],emit:(s:string,p:number)=>number[],maxDeletion:number,maxInsertion:number,radius:number){
 const localSpan=radius+Math.max(maxDeletion,maxInsertion,10)+2,rows=positions.flatMap((p,i)=>refs.every(r=>r.cuts.every(c=>p<c-radius-2||p>c+localSpan))?[i]:[]);
 if(rows.length<80)throw Error('Insufficient invariant flanking signal for length marginals.');
 const families=refs.flatMap((r,h)=>Array.from({length:maxDeletion+maxInsertion+1},(_,i)=>{const size=i-maxDeletion,cut=r.cuts[0],events=size>0?[{start:cut,end:cut,inserted:'A'.repeat(size)}]:size<0?[{start:cut,end:cut-size,inserted:''}]:[],sequence=materialize(r.sequence,events)!;return{reference:h,size,sequence};})),columns=families.map(f=>rows.flatMap(i=>emit(f.sequence,positions[i]))),y=rows.flatMap(i=>observed.slice(4*i,4*i+4)),fit=new SangerFitCache(columns,y).fit(families.map((_,i)=>i));
 return{families:families.map((f,i)=>({reference:f.reference,size:f.size,fraction:fit.weights[i]})),rmse:fit.rmse,rows:rows.map(i=>positions[i]),converged:fit.converged};
}
