import type {TracePlot,TraceBase} from '../../models/sanger-analysis';
const BASES='ACGT';
const IUPAC:Record<string,string>={A:'A',C:'C',G:'G',T:'T',AC:'M',AG:'R',AT:'W',CG:'S',CT:'Y',GT:'K',ACG:'V',ACT:'H',AGT:'D',CGT:'B',ACGT:'N'};
/** Schematic, deterministic base-colour peaks. No noise, sequencing response or fitted data. */
export function idealWTTrace(sequences:string[],start:number,end:number):TracePlot|null {
 if(!sequences.length||end<=start)return null;
 const n=end-start,width=n*10+32,paths={} as TracePlot['paths'];
 for(const base of BASES){const points:string[]=[];
  for(let j=0;j<=n*5;j++){
   const cycle=start+j/5-.5;let signal=0;
   for(let p=Math.max(0,Math.floor(cycle)-1);p<=Math.ceil(cycle)+1;p++){
    const fraction=sequences.filter(s=>s[p]===base).length/sequences.length;
    signal+=fraction*Math.exp(-.5*((cycle-p)/.19)**2);
   }
   points.push(`${j?'L':'M'}${(16+j*2).toFixed(1)},${(124-96*signal).toFixed(2)}`);
  }
  paths[base as TraceBase]=points.join(' ');
 }
 const letters=Array.from({length:n},(_,i)=>{const p=start+i,seen=[...new Set(sequences.map(s=>s[p]).filter(b=>BASES.includes(b||'?')))].sort().join('');return{x:21+i*10,base:sequences.some(s=>!s[p]||s[p]==='N')?'N':IUPAC[seen]||'N',offset:p};});
 return {name:'Ideal WT chromatogram',start,end,width,cutX:0,cutIndex:0,paths,letters};
}
