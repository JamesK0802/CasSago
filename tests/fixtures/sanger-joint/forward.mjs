/** Independent forward simulator. No inference model, regression or fit-position imports.
 * Molecules are complete edited DNA strings. Electrophoresis amplitudes are summed
 * on a fragment-length axis and convolved with asymmetric dye-dependent peaks.
 * Numerical presets are engineering stress conditions, not fitted instrument data.
 */
export const BASES='ACGT';
export function random(seed){let x=seed>>>0;return ()=>((x=(Math.imul(x,1664525)+1013904223)>>>0)+.5)/4294967296;}
export function editAllele(sequence,cut,size,insertion='AGTCA',junction=0){const at=cut+junction;return size<0?sequence.slice(0,at)+sequence.slice(at-size):sequence.slice(0,at)+insertion.repeat(Math.ceil(size/insertion.length)).slice(0,size)+sequence.slice(at);}
export const PRESETS={
 ideal:{noise:0,gain:0,context:0,bleed:0,phase:0,warp:0,decay:.05},
 moderate:{noise:.004,gain:.06,context:.08,bleed:.006,phase:.012,warp:.2,decay:.35},
 stress:{noise:.035,gain:.35,context:.4,bleed:.06,phase:.12,warp:1.1,decay:.8},
};
export function simulate(molecules,{seed,trim=12,length=700,preset='moderate',name='synthetic'}){
 const p=PRESETS[preset],rng=random(seed),gaussian=()=>Math.sqrt(-2*Math.log(rng()))*Math.cos(2*Math.PI*rng());
 const n=Math.min(length,...molecules.map(m=>m.sequence.length))-trim-3;
 if(n<80)throw new Error('Too little molecule sequence for a chromatogram');
 const centers=Array.from({length:n},(_,j)=>35+j*11.7+.0008*j*j+p.warp*Math.sin(j*.067));
 const points=Math.ceil(centers[n-1]+45),channels=Array.from({length:4},()=>new Float64Array(points)),gain=Array.from({length:4},()=>1+p.gain*(rng()*2-1));
 const weights=molecules.map(m=>m.weight),sum=weights.reduce((a,b)=>a+b,0);
 for(let j=0;j<n;j++){
  const amplitudes=[0,0,0,0];
  for(let h=0;h<molecules.length;h++){
   const m=molecules[h],at=j+trim,base=BASES.indexOf(m.sequence[at]);if(base<0)continue;
   // Sequence-context yield is generated from DNA, independently of regression.
   const context=(BASES.indexOf(m.sequence[at-1]||'A')+BASES.indexOf(m.sequence[at+1]||'T')-3)/3;
   const amount=weights[h]/sum*(1+p.context*context);amplitudes[base]+=amount*(1-p.phase);
   const lag=BASES.indexOf(m.sequence[at-1]||m.sequence[at]);amplitudes[lag]+=amount*p.phase;
  }
  const decay=Math.exp(-p.decay*j/n),height=12000*decay*(1+.04*Math.sin(j*.19));
  for(let b=0;b<4;b++){
   const a=(amplitudes[b]*(1-p.bleed)+amplitudes.reduce((s,x,i)=>s+(i===b?0:x),0)*p.bleed/3)*gain[b]*height;
   const center=centers[j]+p.warp*(b-1.5)*.22,width=1.45+j/n*.55;
   for(let t=Math.max(0,Math.floor(center-10));t<Math.min(points,center+14);t++){const z=(t-center)/(width*(t>center?1.17:1));channels[b][t]+=a*Math.exp(-z*z/2);}
  }
 }
 const quantized=channels.map((c,b)=>Uint16Array.from(c,(v,t)=>Math.max(0,Math.min(65535,Math.round(v+12000*(p.noise*.4*(1+Math.sin(t*.003+b))+p.noise*.35*gaussian()))))));
 // Recover local maxima from total signal around the migration grid; do not use
 // inference reference bases. This assumes no missing/collapsed called peaks.
 const peaks=centers.map(center=>{let best=Math.round(center),height=-1;for(let t=Math.round(center)-3;t<=Math.round(center)+3;t++){const h=quantized.reduce((s,c)=>s+c[t],0);if(h>height){height=h;best=t;}}return best;});
 const quality=[],calls=peaks.map(t=>{const v=quantized.map(c=>c[t]),total=v.reduce((a,b)=>a+b,0),best=Math.max(...v);quality.push(Math.max(2,Math.min(50,Math.round(-10*Math.log10(Math.max(.00001,1-best/total))))));return BASES[v.indexOf(best)];}).join('');
 return {buffer:writeAB1(calls,peaks,quantized,quality,name),calls,peaks,parameters:{...p,seed,trim,calledBases:n,preset,quality:'Synthetic dominant-peak confidence; not calibrated Phred',basecaller:'One local maximum per migration-grid cycle; missing/collapsed cycles not simulated'}};
}
/** Standalone ABIF writer; FWO_ GATC, processed DATA9–12, PBAS/PLOC/PCON. */
function writeAB1(sequence,peaks,channels,quality,name){
 const enc=new TextEncoder(),shorts=v=>{const b=new Uint8Array(v.length*2),d=new DataView(b.buffer);v.forEach((x,i)=>d.setUint16(i*2,x));return b;};
 const order='GATC',fields=[['FWO_',1,2,1,enc.encode(order)],['PBAS',2,2,1,enc.encode(sequence)],['PLOC',2,3,2,shorts(peaks)],['PCON',2,1,1,Uint8Array.from(quality)],['SMPL',1,19,1,enc.encode(name+'\0')],...Array.from(order,(b,i)=>['DATA',9+i,3,2,shorts(channels[BASES.indexOf(b)])])];
 const root=128;let offset=root+fields.length*28;const bytes=new Uint8Array(offset+fields.reduce((s,f)=>s+(f[4].length>4?f[4].length:0),0)),view=new DataView(bytes.buffer);
 bytes.set(enc.encode('ABIF'));view.setUint16(4,101);bytes.set(enc.encode('tdir'),6);view.setUint16(14,1023);view.setUint16(16,28);view.setUint32(18,fields.length);view.setUint32(22,fields.length*28);view.setUint32(26,root);
 fields.forEach(([tag,num,type,size,data],i)=>{const at=root+28*i;bytes.set(enc.encode(tag),at);view.setUint32(at+4,num);view.setUint16(at+8,type);view.setUint16(at+10,size);view.setUint32(at+12,data.length/size);view.setUint32(at+16,data.length);if(data.length<=4)bytes.set(data,at+20);else{view.setUint32(at+20,offset);bytes.set(data,offset);offset+=data.length;}});return bytes.buffer;
}
