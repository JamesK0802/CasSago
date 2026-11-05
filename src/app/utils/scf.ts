import {SangerTrace,TRACE_BASES,TraceBase} from '../models/sanger-analysis';
/** Staden SCF 2.x/3.x. V3 stores channel-major second differences modulo 2^bits.
 * https://staden.sourceforge.net/manual/formats_unix_8.html */
export function parseScf(buffer:ArrayBuffer,name='trace.scf'):SangerTrace{
 const v=new DataView(buffer),bytes=new Uint8Array(buffer),text=(p:number,n:number)=>new TextDecoder().decode(bytes.subarray(p,p+n));
 const bounds=(p:number,n:number)=>{if(!Number.isSafeInteger(p)||!Number.isSafeInteger(n)||p<128||n<0||p+n>bytes.length)throw Error('Truncated or invalid SCF data offset.');};
 if(bytes.length<128)throw Error('Truncated SCF header.');if(text(0,4)!=='.scf')throw Error('Invalid SCF signature.');
 const samples=v.getUint32(4),sampleAt=v.getUint32(8),bases=v.getUint32(12),baseAt=v.getUint32(24),version=Number(text(36,4)),size=v.getUint32(40);
 if(!(version>=2&&version<4)||![1,2].includes(size))throw Error('Supported SCF versions are 2.x/3.x with 8- or 16-bit samples.');
 if(samples<40||samples>1000000||bases<40||bases>3000)throw Error('Use a basecalled SCF read of 40–3,000 bases and at most 1,000,000 scan points.');
 bounds(sampleAt,4*samples*size);bounds(baseAt,12*bases);
 const channels={} as Record<TraceBase,Float32Array>,mod=size===1?256:65536,read=(p:number)=>size===1?v.getUint8(p):v.getUint16(p);
 TRACE_BASES.forEach((b,k)=>{const values=Array.from({length:samples},(_,i)=>read(sampleAt+size*(version>=3?k*samples+i:4*i+k)));if(version>=3)for(let pass=0;pass<2;pass++){let sum=0;for(let i=0;i<samples;i++){sum=(sum+values[i])%mod;values[i]=sum;}}channels[b]=Float32Array.from(values);});
 const peaks=new Uint32Array(bases),quality=new Uint8Array(bases);let sequence='';
 for(let i=0;i<bases;i++){
  const peak=v.getUint32(baseAt+(version>=3?4*i:12*i)),base=text(baseAt+(version>=3?8*bases+i:12*i+8),1).toUpperCase().replace('-','N');
  if(!/^[ACGTRYSWKMBDHVN]$/.test(base))throw Error('Invalid SCF base call.');
  if(peak>=samples||(i>0&&peak<=peaks[i-1]))throw Error('SCF peak positions are incomplete or not strictly increasing.');
  peaks[i]=peak;sequence+=base;const k=TRACE_BASES.indexOf(base as TraceBase);quality[i]=k<0?0:v.getUint8(baseAt+(version>=3?4*bases+k*bases+i:12*i+4+k));
 }
 const n=v.getUint32(28),at=v.getUint32(32);let sampleName=name;if(n){bounds(at,n);sampleName=text(at,n).match(/(?:^|\n)NAME=([^\n\0]+)/)?.[1]||name;}
 return {name,sampleName,sequence,peaks,quality:quality.some(q=>q>0)?quality:null,channels};
}
