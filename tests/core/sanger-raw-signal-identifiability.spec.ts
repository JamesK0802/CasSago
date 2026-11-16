import {it,expect} from 'vitest';
import {readFileSync} from 'node:fs';
import {parseSangerTrace} from '../../src/app/utils/sanger-trace';
import {writeSynthTrace} from '../../src/app/utils/sanger-synthetic-benchmark';
const root=new URL('../fixtures/identical-chromatograms/',import.meta.url);
const manifest=JSON.parse(readFileSync(new URL('manifest.json',root),'utf8'));
function parse(bytes:Uint8Array){return parseSangerTrace(bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength)as ArrayBuffer,'synthetic.ab1');}
for(const c of manifest.cases){
 it(c.id+': different finite DNA populations produce identical raw signals',()=>{
  const folder=new URL(c.id+'/',root),baseline=parse(readFileSync(new URL('baseline.ab1',folder)));let originalMass:Record<string,number>|undefined,originalPopulation:string|undefined;
  for(const name of c.names){
   const trace=parse(readFileSync(new URL(name+'.ab1',folder))),population=JSON.parse(readFileSync(new URL(name+'-population.json',folder),'utf8'))as{mass:number;sequence:string;size:number;nominalSequence?:string}[],masses:Record<string,number>={};
   for(const p of population){expect(Number.isInteger(p.mass)).toBe(true);expect(p.mass).toBeGreaterThan(0);masses[p.size]=(masses[p.size]||0)+p.mass;}
   expect(Object.values(masses).reduce((n,v)=>n+v,0)).toBe(10000);
   if(!originalMass){originalMass=masses;originalPopulation=JSON.stringify(population.map(p=>[p.sequence,p.mass]));}else{expect(masses).toEqual(originalMass);expect(JSON.stringify(population.map(p=>[p.sequence,p.mass]))).not.toBe(originalPopulation);}
   if(c.globalSubMass!==undefined){const carriers:Record<string,number>={};for(const p of population){expect(p.nominalSequence).toBeDefined();if(p.sequence!==p.nominalSequence)carriers[p.size]=(carriers[p.size]||0)+p.mass;}expect(Object.values(carriers).reduce((n,v)=>n+v,0)).toBe(c.globalSubMass);expect(carriers).toEqual(c.subCarrierMassByName[name]);}
   const reconstructed=parse(new Uint8Array(writeSynthTrace(population.map(p=>({sequence:p.sequence.slice(c.primerLength),weight:p.mass/10000})),trace.sequence.length,'reconstructed')));
   for(const b of ['A','C','G','T']as const){expect(trace.channels[b]).toEqual(baseline.channels[b]);expect(reconstructed.channels[b]).toEqual(trace.channels[b]);}
   expect(trace.sequence).toBe(baseline.sequence);expect(trace.peaks).toEqual(baseline.peaks);
  }
 });
}
