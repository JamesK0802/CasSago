import {reverseComplement} from '../../utils/biology.utils';
import type {PrimerSite, SangerReference, SangerDesignSettings, SangerCandidate, PrimerBindingScreen} from '../../models/sanger-design';
/** Design guard only: measured read quality still needs checking after sequencing. */
export const SANGER_EDGE_RESERVE = 50;
export function usablePrimerRead(settings:Pick<SangerDesignSettings,'readStart'|'readEnd'>,length:number) {
  return {start:Math.max(SANGER_EDGE_RESERVE,settings.readStart),end:Math.min(settings.readEnd,length)-SANGER_EDGE_RESERVE};
}
/** Enumerate ungapped <=3 substitution sites on both strands. Four disjoint
 * seeds make the search exhaustive at that bound. Intended sites are excluded
 * by reference, coordinates AND orientation, not by similar sequence alone.
 * This is sequence similarity, not a binding-probability prediction. */
export function screenPrimerBinding(oligo:string,refs:SangerReference[],intended:Pick<PrimerSite,'reference'|'start'|'end'|'direction'>[]):PrimerBindingScreen {
  if(!/^[ACGT]{4,}$/.test(oligo))throw Error('Primer screen requires an A/C/G/T oligo.');
  const limit=3;let count=0,minimum=oligo.length+1;
  const sites:PrimerBindingScreen['sites']=[];
  const references=refs.map(ref=>{
    let nearest=oligo.length+1,nearSiteCount=0;
    for(const direction of ['forward','reverse'] as const){
      const dna=direction==='forward'?ref.sequence:reverseComplement(ref.sequence),visited=new Set<number>();
      for(let part=0;part<=limit;part++){
        const lo=Math.floor(part*oligo.length/4),hi=Math.floor((part+1)*oligo.length/4),seed=oligo.slice(lo,hi);
        for(let hit=dna.indexOf(seed);hit>=0;hit=dna.indexOf(seed,hit+1)){
          const at=hit-lo;if(at<0||at+oligo.length>dna.length||visited.has(at))continue;visited.add(at);
          const start=direction==='forward'?at:dna.length-at-oligo.length,end=start+oligo.length;
          if(intended.some(s=>s.reference===ref.name&&s.start===start&&s.end===end&&s.direction===direction))continue;
          let mismatches=0,terminalMismatches=0,unknown=false;
          for(let j=0;j<oligo.length;j++){
            if(!'ACGT'.includes(dna[at+j])){unknown=true;break;}
            if(dna[at+j]!==oligo[j]){mismatches++;if(j>=oligo.length-3)terminalMismatches++;}
            if(mismatches>limit)break;
          }
          if(unknown||mismatches>limit)continue;
          count++;nearSiteCount++;nearest=Math.min(nearest,mismatches);minimum=Math.min(minimum,mismatches);
          const row={reference:ref.name,start,end,direction,mismatches,terminalMismatches};sites.push(row);
          sites.sort((a,b)=>a.mismatches-b.mismatches||a.terminalMismatches-b.terminalMismatches||a.reference.localeCompare(b.reference)||a.start-b.start);
          if(sites.length>12)sites.pop();
        }
      }
    }
    // If no near site was found, obtain the exact nearest known sequence by
    // branch-and-bound Hamming scanning. No binding efficiency is inferred.
    if(nearest>3)for(const direction of ['forward','reverse']as const){
      const dna=direction==='forward'?ref.sequence:reverseComplement(ref.sequence);
      for(let at=0;at+oligo.length<=dna.length;at++){
        const start=direction==='forward'?at:dna.length-at-oligo.length,end=start+oligo.length;
        if(intended.some(s=>s.reference===ref.name&&s.start===start&&s.end===end&&s.direction===direction))continue;
        let errors=0,j=0;
        for(;j<oligo.length;j++){const b=dna[at+j];if(!'ACGT'.includes(b))break;if(b!==oligo[j])errors++;if(errors>=nearest)break;}
        if(j===oligo.length)nearest=errors;
      }
    }
    minimum=Math.min(minimum,nearest);
    return {reference:ref.name,minimumMismatches:nearest<=oligo.length?nearest:4,minimumExact:nearest<=oligo.length,nearSiteCount,hasUnknownBases:ref.sequence.includes('N')};
  });
  return {searchLimit:3,minimumMismatches:minimum<=oligo.length?minimum:4,minimumExact:minimum<=oligo.length,nearSiteCount:count,sites,references};
}
/** Apply only to the existing bounded design shortlist; no edited trace used. */
export function annotatePrimerSafety(candidates:SangerCandidate[],refs:SangerReference[],settings:SangerDesignSettings,progress:(stage:string)=>void=()=>{}){
  if(settings.platform==='ngs')return;
  const cache=new Map<string,PrimerBindingScreen>();
  const check=(sequence:string,sites:PrimerSite[])=>{const key=sequence+'|'+sites.map(s=>`${s.reference}:${s.start}:${s.end}:${s.direction}`).join('|');let value=cache.get(key);if(!value){value=screenPrimerBinding(sequence,refs,sites);cache.set(key,value);}return value;};
  for(const [index,c]of candidates.entries()){if(index%10===0)progress(`Checking unintended primer sites ${index+1}/${candidates.length}`);c.primerSafety={primary:check(c.sequence,c.sites),partner:c.partner?check(c.partner.sequence,c.partner.sites):null};}
}
export function comparePrimerSafety(a:SangerCandidate,b:SangerCandidate):number {
  const ap=a.primerSafety,bp=b.primerSafety;if(!ap||!bp)return 0;
  return bp.primary.minimumMismatches-ap.primary.minimumMismatches || ap.primary.nearSiteCount-bp.primary.nearSiteCount || (bp.partner?.minimumMismatches??4)-(ap.partner?.minimumMismatches??4) || (ap.partner?.nearSiteCount??0)-(bp.partner?.nearSiteCount??0);
}
