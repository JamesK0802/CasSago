import {parseAbif} from './abif';
import {parseScf} from './scf';
/** Detect by content, not by extension. FASTA/FASTQ cannot supply a chromatogram. */
export function parseSangerTrace(buffer:ArrayBuffer,name='trace.ab1'){
 const signature=new TextDecoder().decode(new Uint8Array(buffer,0,Math.min(4,buffer.byteLength)));
 if(signature==='ABIF')return parseAbif(buffer,name);
 if(signature==='.scf')return parseScf(buffer,name);
 throw Error('Choose a basecalled Sanger AB1/ABI or SCF chromatogram.');
}
