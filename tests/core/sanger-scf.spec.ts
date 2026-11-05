import {describe,it,expect} from 'vitest';
import {readFileSync} from 'node:fs';
import {parseScf} from '../../src/app/utils/scf';
import {parseAbif} from '../../src/app/utils/abif';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
import {TRACE_BASES} from '../../src/app/models/sanger-analysis';
function fixture(name:string){const b=readFileSync(new URL('../fixtures/sanger/public-sangerseqR/'+name,import.meta.url));return b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength);}
function encoded(version:number,size:number){
 const points=300,bases=60,at=128,baseAt=at+4*points*size,buf=new ArrayBuffer(baseAt+12*bases),v=new DataView(buf),b=new Uint8Array(buf),txt=(p:number,s:string)=>b.set(new TextEncoder().encode(s),p),mod=size===1?256:65536;
 txt(0,'.scf');v.setUint32(4,points);v.setUint32(8,at);v.setUint32(12,bases);v.setUint32(24,baseAt);txt(36,version+'.00');v.setUint32(40,size);
 const signals=TRACE_BASES.map((_,k)=>Array.from({length:points},(_,i)=>(i*i+31*k+9*i*k)%mod));
 signals.forEach((values,k)=>values.forEach((value,i)=>{const q=version===3?(value-2*(values[i-1]||0)+(values[i-2]||0)+2*mod)%mod:value,p=at+size*(version===3?k*points+i:4*i+k);size===1?v.setUint8(p,q):v.setUint16(p,q);}));
 for(let i=0;i<bases;i++){v.setUint32(baseAt+(version===3?4*i:12*i),2+4*i);for(let k=0;k<4;k++)v.setUint8(baseAt+(version===3?4*bases+k*bases+i:12*i+4+k),k===i%4?30:0);txt(baseAt+(version===3?8*bases+i:12*i+8),TRACE_BASES[i%4]);}
 return {buf,signals};
}
describe('SCF chromatogram input',()=>{
 it.each([[2,1],[2,2],[3,1],[3,2]])('reads v%s %s-byte channel layouts and modulo double deltas',(version,size)=>{const {buf,signals}=encoded(version,size),r=parseScf(buf);TRACE_BASES.forEach((b,k)=>expect([...r.channels[b]]).toEqual(signals[k]));expect([...r.quality!]).toEqual(Array(60).fill(30));expect([...r.peaks]).toEqual(Array.from({length:60},(_,i)=>2+4*i));});
 it('matches all measured channels and positions in independently encoded real AB1 and SCF files',()=>{const a=parseAbif(fixture('heterozygous.ab1')),s=parseScf(fixture('heterozygous.scf'));expect(s.sequence).toBe(a.sequence);expect(s.peaks).toEqual(a.peaks);for(const b of TRACE_BASES)expect(s.channels[b]).toEqual(a.channels[b]);expect(a.quality!.some(q=>q>0)).toBe(true);expect(s.quality).toBeNull();});
 it('accepts actual WT SCF with absent quality without inventing confidence',()=>{const t=parseSangerTrace(fixture('homozygous.scf'),'WT.scf');expect(t.sequence.length).toBe(722);expect(t.quality).toBeNull();});
 it('rejects truncated, oversized, unsupported and nonmonotonic traces',()=>{expect(()=>parseScf(new ArrayBuffer(127))).toThrow(/Truncated/);const {buf}=encoded(3,2);expect(()=>parseScf(buf.slice(0,-1))).toThrow(/offset/);const bad=buf.slice(0);new DataView(bad).setUint32(4,1000001);expect(()=>parseScf(bad)).toThrow(/scan points/);const version=buf.slice(0);new Uint8Array(version).set(new TextEncoder().encode('4.00'),36);expect(()=>parseScf(version)).toThrow(/versions/);const peaks=buf.slice(0),v=new DataView(peaks);v.setUint32(v.getUint32(24)+4,2);expect(()=>parseScf(peaks)).toThrow(/increasing/);expect(()=>parseSangerTrace(new TextEncoder().encode('>FASTA\nACGT').buffer)).toThrow(/chromatogram/);});
});
