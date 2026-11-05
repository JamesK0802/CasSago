import { TRACE_BASES } from '../models/sanger-analysis';
/** Constructed ABIF fixture: explicitly synthetic, exercises the real file reader. */
export function sangerAnalysisDemo():{reference:string;cut:number;control:ArrayBuffer;sample:ArrayBuffer} {
  let seed=718223;const dna=(n:number)=>Array.from({length:n},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return 'ACGT'[seed>>>30];}).join('');
  const sequence=dna(650),cut=230,flank=dna(60);
  const make=(mixed:boolean)=>{
    const count=sequence.length,peaks=Array.from({length:count},(_,i)=>14+i*10),length=count*10+30;
    const channels=TRACE_BASES.map(()=>new Uint16Array(length));
    for(let i=0;i<count;i++) {
      const mixture=mixed&&i>=cut?[[sequence[i],.6],[sequence[Math.min(count-1,i+1)],.25],[sequence[Math.max(0,i-1)],.15]] as const:[[sequence[i],1]] as const;
      for(const [base,weight] of mixture){const channel=TRACE_BASES.indexOf(base as any);for(let d=-5;d<=5;d++){const at=peaks[i]+d;channels[channel][at]+=Math.round(weight*1100*Math.exp(-d*d/4));}}
    }
    for(const c of channels)for(let j=0;j<c.length;j++)c[j]+=2;
    return writeAbif(sequence,peaks,channels,mixed?'Synthetic_sample_60WT_25del1_15ins1':'Synthetic_control');
  };
  return {reference:flank+sequence+dna(60),cut:flank.length+cut,control:make(false),sample:make(true)};
}
/** Minimal standards-shaped ABIF writer used only for the labelled demo. */
export function writeAbif(sequence:string,peaks:number[],channels:Uint16Array[],name:string,order='ACGT',quality=40):ArrayBuffer {
  const encoder=new TextEncoder(),shorts=(values:ArrayLike<number>)=>{const out=new Uint8Array(values.length*2),v=new DataView(out.buffer);Array.from(values).forEach((x,i)=>v.setUint16(i*2,x));return out;};
  const fields=[{tag:'FWO_',number:1,type:2,size:1,data:encoder.encode(order)},
    {tag:'PBAS',number:2,type:2,size:1,data:encoder.encode(sequence)},
    {tag:'PLOC',number:2,type:3,size:2,data:shorts(peaks)},
    {tag:'PCON',number:2,type:2,size:1,data:new Uint8Array(sequence.length).fill(quality)},
    {tag:'SMPL',number:1,type:19,size:1,data:encoder.encode(name+'\0')},
    ...Array.from(order,(base,i)=>({tag:'DATA',number:9+i,type:3,size:2,data:shorts(channels[TRACE_BASES.indexOf(base as any)])})),
    {tag:'PBAS',number:1,type:2,size:1,data:encoder.encode(sequence)},
    {tag:'PLOC',number:1,type:3,size:2,data:shorts(peaks)},
    {tag:'PCON',number:1,type:2,size:1,data:new Uint8Array(sequence.length).fill(quality)},];
  const rootOffset=128;let offset=rootOffset+28*fields.length;
  const total=offset+fields.reduce((n,f)=>n+(f.data.length>4?f.data.length:0),0),buffer=new ArrayBuffer(total),bytes=new Uint8Array(buffer),view=new DataView(buffer);
  bytes.set(encoder.encode('ABIF'));view.setUint16(4,101);bytes.set(encoder.encode('tdir'),6);view.setUint16(14,1023);view.setUint16(16,28);view.setUint32(18,fields.length);view.setUint32(22,fields.length*28);view.setUint32(26,rootOffset);
  fields.forEach((f,i)=>{const at=rootOffset+i*28;bytes.set(encoder.encode(f.tag),at);view.setUint32(at+4,f.number);view.setUint16(at+8,f.type);view.setUint16(at+10,f.size);view.setUint32(at+12,f.data.length/f.size);view.setUint32(at+16,f.data.length);
    if(f.data.length<=4)bytes.set(f.data,at+20);else{view.setUint32(at+20,offset);bytes.set(f.data,offset);offset+=f.data.length;}
  });
  return buffer;
}

/** Known two-amplicon mixture; constructed traces, not biological validation. */
export function sangerPooledDemo() {
  let seed=3921881;const dna=(n:number)=>Array.from({length:n},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return 'ACGT'[seed>>>30];}).join('');
  const primer='CCGGCGGGTACGACCATAAC',first=dna(680),second=dna(60)+first.slice(0,620),cuts=[230,290],refs=[first,second],trim=12;
  const make=(edited:boolean)=>{
    const count=580,peaks=Array.from({length:count},(_,i)=>14+i*10),channels=TRACE_BASES.map(()=>new Uint16Array(count*10+30).fill(2)),calls:string[]=[];
    for(let i=0;i<count;i++){
      const p=i+trim,values=[0,0,0,0];
      refs.forEach((ref,h)=>{
        const parts=edited&&p>=cuts[h]?(h===0?[{size:0,f:.6},{size:-1,f:.25},{size:1,f:.15}]:[{size:0,f:.4},{size:-3,f:.6}]):[{size:0,f:1}];
        for(const part of parts){const base=ref[p-part.size];values[TRACE_BASES.indexOf(base as any)]+=(h===0?.55:.45)*part.f;}
      });
      calls.push(TRACE_BASES[values.indexOf(Math.max(...values))]);
      values.forEach((value,c)=>{for(let d=-5;d<=5;d++)channels[c][peaks[i]+d]+=Math.round(value*1200*Math.exp(-d*d/4));});
    }
    return writeAbif(calls.join(''),peaks,channels,edited?'Synthetic_pooled_sample':'Synthetic_pooled_control');
  };
  return {primer,amplicons:refs.map((sequence,i)=>({name:`Demo-${'AB'[i]}`,sequence:primer+sequence,cut:primer.length+cuts[i]})),control:make(false),sample:make(true)};
}

/** Two nearby substitutions with three unchanged bases; consensus-only fixture. */
export function sangerLocalConsensusDemo(){
 const original=sangerAnalysisDemo(),reference=original.reference,cut=original.cut,positions=[cut-2,cut+2],sequence=Array.from(reference).map((b,i)=>positions.includes(i)?(b==='A'?'C':'A'):b).join('');
 const peaks=Array.from({length:sequence.length},(_,i)=>14+i*10),channels=TRACE_BASES.map(()=>new Uint16Array(sequence.length*10+30).fill(2));
 Array.from(sequence).forEach((b,i)=>{const c=TRACE_BASES.indexOf(b as any);if(c>=0)for(let d=-5;d<=5;d++)channels[c][peaks[i]+d]+=Math.round(1100*Math.exp(-d*d/4));});
 return {reference,cut,positions,sample:writeAbif(sequence,peaks,channels,'Synthetic_local_substitutions')};
}

/** Two targets, including a compensating −1/+1 allele. Entire molecules are
 * synthesized independently of the inference proposal generator. */
export function sangerJointDemo(){
 let seed=9201;const dna=()=>Array.from({length:900},()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return 'ACGT'[seed>>>30];}).join('');
 const primer='ACGTTAGCTACGATCGATGC',refs=[dna(),dna()],cuts=[[260,370],[305,415]];
 const edited=[refs[0].slice(0,260)+refs[0].slice(261,370)+'T'+refs[0].slice(370),refs[1].slice(0,305)+refs[1].slice(307)];
 const make=(sample:boolean)=>{const molecules=sample?[{sequence:refs[0],weight:.3},{sequence:edited[0],weight:.3},{sequence:refs[1],weight:.2},{sequence:edited[1],weight:.2}]:[{sequence:refs[0],weight:.6},{sequence:refs[1],weight:.4}];
  const n=720,peaks=Array.from({length:n},(_,i)=>14+i*10),channels=TRACE_BASES.map(()=>new Uint16Array(n*10+30).fill(2)),calls:string[]=[];
  for(let i=0;i<n;i++){const values=TRACE_BASES.map(b=>molecules.reduce((sum,m)=>sum+(m.sequence[i+12]===b?m.weight:0),0));calls.push(TRACE_BASES[values.indexOf(Math.max(...values))]);values.forEach((v,b)=>{for(let d=-5;d<=5;d++)channels[b][peaks[i]+d]+=Math.round(v*12000*Math.exp(-d*d/4));});}
  return writeAbif(calls.join(''),peaks,channels,sample?'Synthetic_joint_sample':'Synthetic_joint_control');
 };
 return {primer,amplicons:refs.map((sequence,i)=>({name:`Joint-${'AB'[i]}`,sequence:primer+sequence,cut:primer.length+cuts[i][0],cuts:cuts[i].map(c=>c+primer.length)})),control:make(false),sample:make(true)};
}
