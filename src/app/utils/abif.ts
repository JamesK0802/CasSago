import { SangerTrace, TRACE_BASES, TraceBase } from '../models/sanger-analysis';
interface Entry {type:number;size:number;count:number;bytes:number;offset:number}
/** ABIF big-endian directory reader. Only analyzed sequencing traces are used. */
export function parseAbif(buffer:ArrayBuffer,name='trace.ab1'):SangerTrace {
  const view=new DataView(buffer),bytes=new Uint8Array(buffer);
  const decoder=new TextDecoder();
  const text=(at:number,n:number)=>decoder.decode(bytes.subarray(at,at+n));
  const bounds=(at:number,n:number)=>{if(!Number.isSafeInteger(at)||!Number.isSafeInteger(n)||at<0||n<0||at+n>buffer.byteLength)throw new Error('Truncated or invalid ABIF data offset.');};
  bounds(0,34);if(text(0,4)!=='ABIF') throw new Error('Choose an ABIF Sanger .ab1/.abi file. FASTQ, FASTA and SCF are not supported here.');
  const rootCount=view.getUint32(18),rootOffset=view.getUint32(26);
  if(rootCount<1||rootCount>4096)throw new Error('Invalid ABIF directory length.');
  bounds(rootOffset,rootCount*28);
  const entries=new Map<string,Entry>();
  for(let i=0;i<rootCount;i++) {
    const at=rootOffset+i*28,key=text(at,4)+view.getUint32(at+4);
    const entry:Entry={type:view.getUint16(at+8),size:view.getUint16(at+10),count:view.getUint32(at+12),bytes:view.getUint32(at+16),offset:view.getUint32(at+20)};
    if(entry.bytes<=4)entry.offset=at+20;
    bounds(entry.offset,entry.bytes);
    if(entry.count*entry.size>entry.bytes)throw new Error(`Invalid ABIF element size: ${key}.`);
    entries.set(key,entry);
  }
  const get=(tag:string)=>{const e=entries.get(tag);if(!e)throw new Error(`Missing ${tag}: an analyzed, basecalled Sanger trace is required.`);return e;};
  const chars=(e:Entry)=>text(e.offset,e.bytes).replace(/\0+$/,'');
  const pair=entries.has('PBAS2')&&entries.has('PLOC2')?2:1;
  const baseEntry=get('PBAS'+pair);if(baseEntry.size!==1)throw new Error('Unsupported ABIF base-call type.');
  const sequence=chars(baseEntry).toUpperCase();
  if(sequence.length<40||sequence.length>3000||/[^ACGTRYSWKMBDHVN]/.test(sequence))throw new Error('Use a basecalled Sanger read of 40–3,000 bases.');
  const order=chars(get('FWO_1')).slice(0,4);
  if(order.length!==4||new Set(order).size!==4||/[^ACGT]/.test(order))throw new Error('ABIF dye order FWO_1 is missing or invalid.');
  const numeric=(e:Entry):number[]=>{
    if(!((e.size===1&&[1,2].includes(e.type))||(e.size===2&&[3,4].includes(e.type))||(e.size===4&&e.type===5)))throw new Error('Unsupported ABIF numeric array.');
    return Array.from({length:e.count},(_,i)=>e.size===1?view.getUint8(e.offset+i):e.size===2?(e.type===4?view.getInt16(e.offset+2*i):view.getUint16(e.offset+2*i)):view.getInt32(e.offset+4*i));
  };
  const channels={} as Record<TraceBase,Float32Array>;
  for(let i=0;i<4;i++) {
    const e=get('DATA'+(9+i));if(e.count<40||e.count>1000000)throw new Error('ABIF trace has an unsupported number of data points.');
    channels[order[i] as TraceBase]=Float32Array.from(numeric(e),v=>Math.max(0,v));
  }
  const length=channels.A.length;
  if(TRACE_BASES.some(b=>channels[b].length!==length))throw new Error('ABIF channel lengths do not match.');
  const positions=numeric(get('PLOC'+pair));
  if(positions.length!==sequence.length||positions.some((p,i)=>p<0||p>=length||(i>0&&p<=positions[i-1])))throw new Error('ABIF peak positions are incomplete or not strictly increasing.');
  const q=entries.get('PCON'+pair);const quality=q?Uint8Array.from(numeric(q)):null;
  if(quality&&quality.length!==sequence.length)throw new Error('ABIF quality values do not match base calls.');
  const sample=entries.get('SMPL1');
  let sampleName=sample?chars(sample):name;
  if(sample?.type===18)sampleName=sampleName.slice(1,1+bytes[sample.offset]);
  return {name,sampleName,sequence,peaks:Uint32Array.from(positions),quality,channels};
}
