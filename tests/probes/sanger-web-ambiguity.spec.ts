import {it,expect} from 'vitest';
import {readFileSync,writeFileSync} from 'node:fs';
import {writeSynthTrace} from '../../src/app/utils/sanger-synthetic-benchmark';
const run=process.env.CASSAGO_WEB_AMBIGUITY==='1'?it:it.skip;
run('checks an actual same-signal alternative independently of blind inference',()=>{
 const root='../benchmarks/web-3ref-2target',data=JSON.parse(readFileSync(root+'/inputs/truth.json','utf8')),rows=data.truth,ref=data.recipe.find((r:any)=>r.sourceReference==='Demo-B'),parts=rows.map((r:any)=>({sequence:r.sequence.slice(ref.primerLength),weight:r.share/100*r.percent/100}));
 const indices=rows.flatMap((r:any,i:number)=>r.sourceReference==='Demo-B'&&r.events.length===1&&r.events[0].kind==='insertion'&&r.events[0].target===0&&r.events[0].sequence.length===6?[i]:[]);expect(indices).toHaveLength(2);
 const [a,b]=indices,ea=rows[a].events[0],eb=rows[b].events[0],pa=ref.cuts[0]+ea.offset-ref.primerLength,pb=ref.cuts[0]+eb.offset-ref.primerLength,at=Math.max(pa,pb);expect(at).toBeLessThan(Math.min(pa+6,pb+6));
 const alternative=parts.map((r:any)=>({...r})),sa=parts[a].sequence,sb=parts[b].sequence;expect(sa[at]).not.toBe(sb[at]);alternative[a].sequence=sa.slice(0,at)+sb[at]+sa.slice(at+1);alternative[b].sequence=sb.slice(0,at)+sa[at]+sb.slice(at+1);
 const original=writeSynthTrace(parts,data.cycles,'edit'),other=writeSynthTrace(alternative,data.cycles,'edit');expect(new Uint8Array(original)).toEqual(new Uint8Array(readFileSync(root+'/inputs/SYNTHETIC_edit.ab1')));expect(new Uint8Array(other)).toEqual(new Uint8Array(original));
 writeFileSync(root+'/same-signal-proof.json',JSON.stringify({reference:ref.name,groups:indices.map(i=>rows[i].allele),original:indices.map(i=>rows[i].events[0]),alternativeInsertedSequences:indices.map((i,j)=>{const event=rows[i].events[0],position=at-(j?pb:pa);return event.sequence.slice(0,position)+(j?sa[at]:sb[at])+event.sequence.slice(position+1);}),identicalAB1Bytes:true,scope:'Evaluation-only construction. Neither recipe nor alternatives are passed to inference.'},null,2));
});
