import '@angular/compiler';import{it}from'vitest';import{readFileSync,writeFileSync}from'node:fs';import{parseSangerTrace}from'../../src/app/utils/sanger-trace';import{sangerEditedSequence}from'../../src/app/utils/sanger-sequence-groups';
(process.env.CASSAGO_BIN_COLLISION==='1'?it:it.skip)('audits repeated emission masks against distinct original ADC bins without truth',()=>{
 const folder=process.env.CASSAGO_RAW_INPUT!,result=JSON.parse(readFileSync(process.env.CASSAGO_RAW_RESULT!,'utf8')),read=(n:string)=>{const b=readFileSync(folder+'/'+n);return parseSangerTrace(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),n);},sample=read('edit.ab1'),seqs=result.annotations.flatMap((a:any,h:number)=>a.components.map((c:any)=>sangerEditedSequence(a.reference,c.events).slice(result.pooled.preflight.amplicons[h].primerEnd))),map=new Map<string,any>();
 for(let p=84;p<Math.min(sample.peaks.length-3,...seqs.map((s:string)=>s.length));p++)for(const b of'ACGT')for(let d=0;d<=4;d++){
 const mask=seqs.flatMap((s:string,i:number)=>s[p]===b?[i]:[]).join(','),key=mask+'/'+d,value=sample.channels[b as'A'][sample.peaks[p]+d],old=map.get(key),row={p,b,d,value};if(old){old.min=Math.min(old.min,value);old.max=Math.max(old.max,value);old.rows.push(row);}else map.set(key,{mask,d,min:value,max:value,rows:[row]});
 }
 const conflicts=[...map.values()].filter(r=>r.max>r.min);writeFileSync(process.env.CASSAGO_RAW_OUTPUT!,JSON.stringify({noTruth:true,conflicts:conflicts.map(r=>({...r,rows:r.rows.slice(0,30)}))},null,2));
});
