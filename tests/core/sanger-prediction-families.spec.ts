import{it,expect}from'vitest';
import{junctionPredictionEvidence}from'../../src/app/workers/core/sanger-junction-evidence';
import{junctionBase,JunctionFit}from'../../src/app/workers/core/sanger-junction-mixture';
it('tests an observable edit family instead of erasing all mutually substitutable deletion lengths',()=>{
 const sequence='ACGT'.repeat(150),cut=180,positions=Array.from({length:380},(_,p)=>p+10),models=[{reference:0,size:0,start:cut,inserted:''},{reference:0,size:-1,start:cut,inserted:''},{reference:0,size:-5,start:cut,inserted:''}],columns=models.map(m=>positions.flatMap(p=>Array.from('ACGT',b=>Number(junctionBase(sequence,m,p)===b)))),weights=[.5,.25,.25],observed=columns[0].map((_,k)=>columns.reduce((s,c,i)=>s+weights[i]*c[k],0)),full:JunctionFit={models,columns,weights,r2:1,rmse:0,converged:true,iterations:0,duplicateLengths:0,objective:0,predictColumns:rows=>models.map(m=>rows.flatMap(p=>Array.from('ACGT',b=>Number(junctionBase(sequence,m,p)===b))))},evidence=junctionPredictionEvidence({references:[{sequence,cut}],positions,observed,maxDeletion:5,maxInsertion:1,radius:0,lag:0,theta:Array(8).fill(0)},full,[1,2]);
 expect(columns[1]).toEqual(columns[2]);for(const i of[1,2]){expect(evidence.get(i)!.foldGains.every(g=>g>1)).toBe(true);expect(evidence.get(i)!.lengthFoldGains.every(g=>g>1)).toBe(true);}
});
