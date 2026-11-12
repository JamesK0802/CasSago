import {usablePrimerRead,annotatePrimerSafety,comparePrimerSafety} from './sanger-primer-safety';
import {screenSharedSangerWT} from './sanger-primer-separation';
import { alignTraceSequence } from './sanger-analysis';
import { reverseComplement } from '../../utils/biology.utils';
import {designProgressiveSharedPrimers} from './sanger-shared-primer-designer';
import { CommonPartner, PrimerProperties, PrimerSite, SangerCandidate, SangerDesignResult,
  SangerDesignSettings, SangerReference, ThreePrimeVariant } from '../../models/sanger-design';

// SantaLucia (1998), unified DNA/DNA nearest-neighbour parameters.
const NN: Record<string, [number, number]> = {
  AA: [-7.9,-22.2], TT: [-7.9,-22.2], AT: [-7.2,-20.4], TA: [-7.2,-21.3],
  CA: [-8.5,-22.7], TG: [-8.5,-22.7], GT: [-8.4,-22.4], AC: [-8.4,-22.4],
  CT: [-7.8,-21.0], AG: [-7.8,-21.0], GA: [-8.2,-22.2], TC: [-8.2,-22.2],
  CG: [-10.6,-27.2], GC: [-9.8,-24.4], GG: [-8.0,-19.9], CC: [-8.0,-19.9],
};
export function nearestNeighbourTm(seq: string, sodiumMm: number, primerNm: number): number {
  if (!/^[ACGT]{2,}$/.test(seq) || sodiumMm <= 0 || primerNm <= 0) return NaN;
  let h = .2, s = -5.7;
  for (let i = 0; i < seq.length - 1; i++) { const p = NN[seq.slice(i,i+2)]; h += p[0]; s += p[1]; }
  for (const base of [seq[0],seq[seq.length-1]]) if (base === 'A' || base === 'T') { h += 2.2; s += 6.9; }
  const symmetric = seq === reverseComplement(seq);
  if (symmetric) s -= 1.4;
  s += .368 * (seq.length - 1) * Math.log(sodiumMm / 1000);
  return 1000 * h / (s + 1.987 * Math.log(primerNm * 1e-9 / (symmetric ? 1 : 4))) - 273.15;
}

/** Contiguous complement screens; deliberately not presented as thermodynamic ΔG. */
export function complementarity(a: string, b: string): { any: number; end: number } {
  const rc = reverseComplement(b);
  let any = 0, end = 0;
  for (let shift = -rc.length + 1; shift < a.length; shift++) {
    let run = 0;
    for (let i = 0; i < a.length; i++) {
      const j = i - shift;
      if (j >= 0 && j < rc.length && a[i] === rc[j]) {
        run++; any = Math.max(any,run);
        if (i === a.length-1 || j-run+1 === 0) end = Math.max(end,run);
      } else run = 0;
    }
  }
  return { any, end };
}
const COMPLEMENT:Record<string,string>={A:'T',T:'A',C:'G',G:'C'};
function hairpinStem(seq: string): number {
  let best = 0;
  for (let i = 0; i < seq.length; i++) for (let j = i + 4; j < seq.length; j++) {
    let k = 0;
    while (i+k < j-k-3 && seq[i+k] === COMPLEMENT[seq[j-k]]) k++;
    best = Math.max(best,k);
  }
  return best;
}
export function properties(seq: string, s: SangerDesignSettings): PrimerProperties | null {
  if (/[^ACGT]/.test(seq) || /([ACGT])\1{4}/.test(seq)) return null;
  const gc = 100 * (seq.match(/[GC]/g)?.length || 0) / seq.length;
  const tm = nearestNeighbourTm(seq,s.sodiumMm,s.primerNm);
  if (gc < s.gcMin || gc > s.gcMax || tm < s.tmMin || tm > s.tmMax) return null;
  const comp = complementarity(seq,seq), hairpin = hairpinStem(seq);
  if (comp.end >= 5 || comp.any >= 9 || hairpin >= 6) return null;
  return { sequence: seq, tm, gc, selfEnd: comp.end, selfAny: comp.any, hairpinStem: hairpin };
}

type Direction = 'forward' | 'reverse';
/** Enumerate target-local oligos, then scan the entire supplied sequence for
 * those candidates. This avoids indexing every genomic k-mer in a long reference. */
export function localOligos(ref: SangerReference, k: number, direction: Direction, reach: number): Set<string> {
  const dna = direction === 'forward' ? ref.sequence : reverseComplement(ref.sequence);
  const cuts = ref.cuts.map(c => direction === 'forward' ? c : dna.length-c);
  const lo = Math.max(0,Math.min(...cuts)-reach-k);
  const hi = Math.min(dna.length-k,Math.max(...cuts)+reach);
  const out = new Set<string>();
  for (let at=lo;at<=hi;at++) {
    const seq=dna.slice(at,at+k);if (!/[^ACGT]/.test(seq)) out.add(seq);
  }
  return out;
}
function candidateIndex(dna: string, wanted: Set<string>, k: number): Map<string,number[]> {
  const out = new Map<string,number[]>();
  for (let at=0;at<=dna.length-k;at++) {
    const seq=dna.slice(at,at+k);
    if (!wanted.has(seq)) continue;
    const positions=out.get(seq);
    // Only uniqueness is needed. Keep at most two positions, including repeats.
    if (!positions) out.set(seq,[at]);else if (positions.length<2) positions.push(at);
  }
  return out;
}
export function exactBindingSites(seq: string, ref: SangerReference): { start: number; end: number; direction: Direction }[] {
  const hits: { start: number; end: number; direction: Direction }[] = [];
  for (const direction of ['forward','reverse'] as const) {
    const dna = direction === 'forward' ? ref.sequence : reverseComplement(ref.sequence);
    for (let at = dna.indexOf(seq); at >= 0; at = dna.indexOf(seq,at+1)) {
      hits.push({ direction, start: direction === 'forward' ? at : dna.length-at-seq.length,
        end: direction === 'forward' ? at+seq.length : dna.length-at });
    }
  }
  return hits;
}
export function site(ref: SangerReference, start: number, k: number, direction: Direction, s: SangerDesignSettings): PrimerSite {
  const dna = direction === 'forward' ? ref.sequence : reverseComplement(ref.sequence);
  return { reference: ref.name, direction,
    start: direction === 'forward' ? start : dna.length-start-k,
    end: direction === 'forward' ? start+k : dna.length-start,
    cutOffsets: ref.cuts.map(c => (direction === 'forward' ? c : dna.length-c) - start-k),
    extension: dna.slice(start+k,start+k+s.readEnd) };
}
export function windows(ref: SangerReference, s: SangerDesignSettings): [number,number][] {
  return ref.cuts.map(c => [c-Math.floor(s.window/2),c+Math.ceil(s.window/2)]);
}
export function primerTargetDistance(c: SangerCandidate): number {return Math.max(...c.sites.flatMap(hit=>hit.cutOffsets));}
export function comparePrimerDistance(a:SangerCandidate,b:SangerCandidate):number {return primerTargetDistance(a)-primerTargetDistance(b) || Math.abs(a.tm-60)-Math.abs(b.tm-60) || (Math.max(...(a.partner?.products||[0]))-Math.max(...(b.partner?.products||[0]))) || a.id.localeCompare(b.id);}
export function coversRead(sites: PrimerSite[], refs: SangerReference[], s: SangerDesignSettings): boolean {
  // Quantitative Sanger needs an intact upstream calibration interval.
  // Match the template handoff contract even for a single PCR product.
  if(s.platform!=='ngs'&&s.workflow!=='shared-pcr'&&sites.some(hit=>hit.cutOffsets.some(c=>c<usablePrimerRead(s,hit.extension.length).start+72)))return false;
  return sites.every((hit,i) => windows(refs[i],s).every(([lo,hi]) => {
    const start = hit.direction === 'forward' ? lo-hit.end : hit.start-hi;
    const end = hit.direction === 'forward' ? hi-hit.end : hit.start-lo;
    return lo >= 0 && hi <= refs[i].sequence.length && (s.platform==='ngs' ? start>=0 : start-s.indelRange >= usablePrimerRead(s,hit.extension.length).start && end+s.indelRange <= usablePrimerRead(s,hit.extension.length).end);
  }));
}
export function minOffset(sites: PrimerSite[]): number {
  // All ref/target cuts share the sequencing-primer read axis. A different
  // gene, guide or target index does not exempt a pair from the spacing rule.
  // Keep duplicates: two distinct targets at the same read cycle have gap 0.
  const cuts=sites.flatMap(s=>s.cutOffsets).sort((a,b)=>a-b);
  let gap = cuts.length<2 ? 0 : Infinity;
  for(let i=1;i<cuts.length;i++)gap=Math.min(gap,cuts[i]-cuts[i-1]);
  return gap;
}
export function traceDistance(sites: PrimerSite[], s: SangerDesignSettings): [number,number] {
  let differences = sites.length<2 ? 0 : Infinity, compared = sites.length<2 ? Math.max(0,Math.min(s.readEnd,sites[0].extension.length)-s.readStart) : Infinity;
  for (let a = 0; a < sites.length; a++) for (let b = a+1; b < sites.length; b++) {
    const range=s.platform==='ngs'?{start:s.readStart,end:Math.min(s.readEnd,sites[a].extension.length,sites[b].extension.length)}:usablePrimerRead(s,Math.min(sites[a].extension.length,sites[b].extension.length)),end=range.end;
    let diff = 0, observed = 0;
    for (let i = range.start; i < end; i++) if (/[ACGT]/.test(sites[a].extension[i]) && /[ACGT]/.test(sites[b].extension[i])) {
      observed++; if (sites[a].extension[i] !== sites[b].extension[i]) diff++;
    }
    differences = Math.min(differences,diff); compared = Math.min(compared,observed);
  }
  return [differences,compared];
}

export function validateDesign(refs: SangerReference[], s: SangerDesignSettings): void {
  if (s.workflow==='standard-pcr' ? refs.length!==1 : refs.length<2 || refs.length>8) throw new Error(s.workflow==='standard-pcr'?'Use one reference for single-reference PCR.':'Use 2–8 references for shared or homeolog-specific design.');
  if (new Set(refs.map(r => r.name)).size !== refs.length || refs.some(r => !r.name)) throw new Error('Reference names must be non-empty and unique.');
  if (refs.some(r => !r.sequence || /[^ACGTN]/.test(r.sequence))) throw new Error('References must contain A, C, G, T or N.');
  if (refs.reduce((n,r) => n+r.sequence.length,0) > 8000000 || refs.some(r => r.sequence.length > 2000000)) throw new Error('Use full reference sequences up to 2,000,000 bp each and 8,000,000 bp total.');
  if (!refs[0].cuts.length || refs[0].cuts.length > 8 || refs.some(r => r.cuts.length !== refs[0].cuts.length)) throw new Error('Give each reference the same number of corresponding targets (1–8), in the same target order.');
  if (refs.some(r => r.cuts.some(c => !Number.isInteger(c) || c < 0 || c > r.sequence.length))) throw new Error('Invalid cut coordinates.');
  for (const [key,value] of Object.entries(s)) if (typeof value === 'number' && !Number.isFinite(value)) throw new Error(`Invalid ${key}.`);
  const ints = [s.primerMin,s.primerMax,s.window,s.readStart,s.readEnd,s.minOffset,s.indelRange,s.productMin,s.productMax];
  if (ints.some(v => !Number.isInteger(v))) throw new Error('Lengths, positions and offsets must be whole numbers.');
  if (s.primerMin < 16 || s.primerMax > 30 || s.primerMax < s.primerMin) throw new Error('Primer length must be within 16–30 nt.');
  if (s.window < 10 || s.window > 2000 || s.readStart < 0 || s.readEnd > 1500 || s.readEnd <= s.readStart) throw new Error('Check required window and usable read range (maximum 1,500 bp).');
  if (s.minOffset < 1 || s.indelRange < 0 || s.indelRange > 100) throw new Error('Use a positive offset and an indel range of 0–100 bp.');
  if (s.sodiumMm < 1 || s.sodiumMm > 1000 || s.primerNm < 1 || s.primerNm > 10000) throw new Error('Check Na+ and oligo concentration.');
  if (s.gcMin < 0 || s.gcMax > 100 || s.gcMin > s.gcMax || s.tmMin > s.tmMax || s.tmMin < 0 || s.tmMax > 100) throw new Error('Check GC and Tm ranges.');
  if (s.productMin < 1 || s.productMax > 20000 || s.productMin > s.productMax || s.maxTmDifference < 0) throw new Error('Check PCR product and paired-primer Tm limits.');
  if (!['forward','reverse','both'].includes(s.direction)) throw new Error('Invalid reading direction.');
  if (!['shared-sanger','shared-pcr','homeolog-pcr','standard-pcr'].includes(s.workflow)) throw new Error('Invalid design workflow.');
  if (!Number.isInteger(s.specificityMismatches) || s.specificityMismatches < 0 || s.specificityMismatches > 2 || typeof s.specificBoth !== 'boolean') throw new Error('Use 0–2 mismatches for the input-reference specificity screen.');
  if (s.workflow === 'homeolog-pcr' && refs.some(r=>r.sequence.includes('N'))) throw new Error('Resolve N bases before the homeolog-specific mismatch screen (A/C/G/T references required).');
  if (!['three-prime','mismatch'].includes(s.specificityStrategy) || ![1,3].includes(s.threePrimeWindow)) throw new Error('Choose a 3′-variant or mismatch-only specificity strategy and a terminal window of 1 or 3 bases.');
  if(s.platform!==undefined&&!['sanger','ngs'].includes(s.platform))throw new Error('Invalid sequencing platform.');
  if(s.platform==='ngs'&&s.workflow==='shared-sanger')throw new Error('Offset-separated shared Sanger requires the Sanger platform.');
  if (typeof s.requirePartner !== 'boolean') throw new Error('Invalid PCR partner setting.');
}

/** Prefer a usable opposite sequencing read for multi-target designs, without
 * adding a design mode or rejecting otherwise valid one-direction PCR pairs. */
export function oppositeReadReady(sites:PrimerSite[],refs:SangerReference[],s:SangerDesignSettings):boolean {
 if(!refs.some(r=>r.cuts.length>1)||!coversRead(sites,refs,s))return false;
 if(s.workflow==='shared-sanger'&&minOffset(sites)<=2*s.indelRange)return false;
 return true;
}
function partnersFor(c: SangerCandidate, refs: SangerReference[], s: SangerDesignSettings,
  pool: { props: PrimerProperties; sites: PrimerSite[]; direction: Direction }[]): CommonPartner[] {
  const direction = c.direction === 'forward' ? 'reverse' : 'forward';
  const ranked:{partner:CommonPartner;score:number}[]=[];
  for (const item of pool) {
    if (item.direction !== direction) continue;
    const props = item.props;
    if (Math.abs(props.tm-c.tm) > s.maxTmDifference) continue;
    const pair = complementarity(c.sequence,props.sequence);
    if (pair.end >= 5 || pair.any >= 9) continue;
    const sites: PrimerSite[] = [], products: number[] = [];
    for (let i = 0; i < refs.length; i++) {
      const hit = item.sites[i], bound = windows(refs[i],s);
      if (direction === 'reverse' ? hit.start < Math.max(...bound.map(w => w[1]))+(s.platform==='ngs'?0:s.indelRange) : hit.end > Math.min(...bound.map(w => w[0]))-(s.platform==='ngs'?0:s.indelRange)) break;
      const product = c.direction === 'forward' ? hit.end-c.sites[i].start : c.sites[i].end-hit.start;
      if (product < s.productMin || product > s.productMax) break;
      sites.push(hit); products.push(product);
    }
    if (sites.length !== refs.length) continue;
    if(s.workflow==='shared-sanger'){const requiredEnd=Math.max(...c.sites.flatMap(hit=>hit.cutOffsets))+Math.ceil(s.window/2)+s.indelRange;if(Math.min(...products.map(product=>product-c.sequence.length))<requiredEnd)continue;}
    const score = Math.abs(props.tm-c.tm)*10 + Math.max(...products)/1000 + pair.end;
    ranked.push({partner:{...props,sites:sites.map((hit,i)=>({...hit,extension:hit.extension.slice(0,products[i]-props.sequence.length)})),products,pairEnd:pair.end},score});
  }
  return ranked.sort((a,b)=>Number(oppositeReadReady(b.partner.sites,refs,s))-Number(oppositeReadReady(a.partner.sites,refs,s))||a.score-b.score||a.partner.sequence.localeCompare(b.partner.sequence)).slice(0,3).map(r=>r.partner);
}

export function designSangerPrimers(refs: SangerReference[], s: SangerDesignSettings,
  progress: (fraction: number, stage: string) => void = () => {}): SangerDesignResult {
  validateDesign(refs,s);
  if (s.workflow === 'homeolog-pcr') return designSpecificPairs(refs,s,progress);
  if(s.workflow==='shared-sanger' && s.platform!=='ngs' && s.sharedPrimerPolicy==='progressive')return designProgressiveSharedPrimers(refs,s,progress);
  const stats = { enumerated:0,qualityPassed:0,commonUnique:0,geometryPassed:0,offsetPassed:0,partnerChecked:0 };
  const candidates: SangerCandidate[] = [];
  const canPair=refs.every(r=>r.sequence.length>=s.productMin);
  const sharedPool: { props: PrimerProperties; sites: PrimerSite[]; direction: Direction }[] = [];
  const forward = refs.map(r => r.sequence), reverse = refs.map(r => reverseComplement(r.sequence));
  for (let k = s.primerMin; k <= s.primerMax; k++) {
    progress(.65*(k-s.primerMin)/(s.primerMax-s.primerMin+1),`Searching ${k}-nt shared primers`);
    const wanted = new Set<string>();
    const quality = new Map<string,PrimerProperties>();
    for (const direction of ['forward','reverse'] as const) for (const seq of localOligos(refs[0],k,direction,canPair ? Math.max(s.readEnd,s.productMax) : s.readEnd)) {
      if (wanted.has(seq)) continue;
      const props=properties(seq,s); if (props) {wanted.add(seq);quality.set(seq,props);}
    }
    const fi = forward.map(d => candidateIndex(d,wanted,k)), ri = reverse.map(d => candidateIndex(d,wanted,k));
    for (const direction of ['forward','reverse'] as const) {
      const intended = direction === 'forward' ? fi : ri, opposite = direction === 'forward' ? ri : fi;
      for (const [seq,positions] of intended[0]) {
        stats.enumerated++;
        if (positions.length !== 1 || opposite[0].has(seq)) continue;
        const props = quality.get(seq)!; stats.qualityPassed++;
        const hits: PrimerSite[] = [];
        for (let i = 0; i < refs.length; i++) {
          const list = intended[i].get(seq);
          if (!list || list.length !== 1 || opposite[i].has(seq)) break;
          hits.push(site(refs[i],list[0],k,direction,s));
        }
        if (hits.length !== refs.length) continue; stats.commonUnique++;
        sharedPool.push({ props,sites:hits,direction });
        if (s.direction !== 'both' && s.direction !== direction) continue;
        if (!coversRead(hits,refs,s)) continue; stats.geometryPassed++;
        const offset = minOffset(hits); if (s.workflow==='shared-sanger' && offset < s.minOffset) continue; stats.offsetPassed++;
        const [diff,bases] = traceDistance(hits,s);
        candidates.push({ ...props,id:`${direction}-${seq}`,direction,sites:hits,minOffset:offset,
          minTraceDifferences:diff,traceComparedBases:bases,indelOffsetOverlap:s.workflow==='shared-sanger' && offset <= 2*s.indelRange,partner:null });
      }
    }
  }
  candidates.sort(comparePrimerDistance);
  const returned: SangerCandidate[] = [],excludedCandidates:{id:string;reason:string}[]=[];
  const screen=(c:SangerCandidate)=>{if(s.workflow!=='shared-sanger'||s.platform==='ngs')return true;const range=usablePrimerRead(s,Math.min(...c.sites.map(h=>h.extension.length)));const check=screenSharedSangerWT(c.sites,{readStart:range.start,readEnd:range.end});if(!check.passed){excludedCandidates.push({id:c.id,reason:check.reason});return false;}c.wtSeparation={minimumDifferences:check.minimumDifferences,comparedBases:check.comparedBases,condition:check.condition,warnings:check.warnings};return true;};
  // Bounded primary search, with multiple opposite-primer choices per primary.
  const limit = Math.min(candidates.length,600);
  for (let i = 0; i < limit; i++) {
    progress(.65+.34*i/Math.max(1,limit),`Checking PCR partner ${i+1}/${limit}`);
    const c = candidates[i]; stats.partnerChecked++;
    const partners = canPair ? partnersFor(c,refs,s,sharedPool) : [];
    for(const partner of partners){
      const paired={...c,id:c.id+'-'+partner.sequence,partner,targetDistance:primerTargetDistance(c),sites:c.sites.map((hit,j)=>({...hit,extension:hit.extension.slice(0,partner.products[j]-c.sequence.length)}))};
      [paired.minTraceDifferences,paired.traceComparedBases]=traceDistance(paired.sites,s);if(coversRead(paired.sites,refs,s)&&screen(paired))returned.push(paired);
    }
    if(!partners.length&&!s.requirePartner&&s.workflow==='shared-sanger'){const unpaired={...c,targetDistance:primerTargetDistance(c)};if(screen(unpaired))returned.push(unpaired);}
    if(returned.length>=300&&i+1<limit&&primerTargetDistance(candidates[i+1])>primerTargetDistance(c))break;
  }
  annotatePrimerSafety(returned,refs,s,stage=>progress(.99,stage));
  returned.sort((a,b)=>comparePrimerSafety(a,b)|| (s.platform!=='ngs'&&s.workflow==='shared-sanger'?b.minTraceDifferences-a.minTraceDifferences:0)||Number(!!b.partner&&oppositeReadReady(b.partner.sites,refs,s))-Number(!!a.partner&&oppositeReadReady(a.partner.sites,refs,s))||comparePrimerDistance(a,b));const capped=returned.length>300;
  progress(1,'Complete');
  return { candidates:returned.slice(0,300), stats:{...stats,separabilityRejected:excludedCandidates.length},excludedCandidates, complete:!capped && stats.partnerChecked === candidates.length,
    references:refs,settings:{...s} };
}

/** Pigeonhole seed search for ungapped matches with <=m substitutions.
 * m+1 disjoint seeds guarantee at least one exact seed in every such match. */
export function nearBindingSiteCount(oligo: string, strands: string[], mismatches: number, limit=2): number {
  if (!/^[ACGT]+$/.test(oligo) || mismatches<0 || mismatches>2) throw new Error('Invalid mismatch search.');
  let count=0;
  for (const dna of strands) {
    const visited=new Set<number>();
    for (let part=0;part<=mismatches;part++) {
      const lo=Math.floor(part*oligo.length/(mismatches+1));
      const hi=Math.floor((part+1)*oligo.length/(mismatches+1));
      const seed=oligo.slice(lo,hi);
      for (let hit=dna.indexOf(seed);hit>=0;hit=dna.indexOf(seed,hit+1)) {
        const at=hit-lo;
        if (at<0 || at+oligo.length>dna.length || visited.has(at)) continue;
        visited.add(at);
        let errors=0;
        for (let j=0;j<oligo.length;j++) {if (dna[at+j]!==oligo[j]) errors++;if (errors>mismatches) break;}
        if (errors<=mismatches && ++count>=limit) return count;
      }
    }
  }
  return count;
}

interface VariantMap {background:string;identity:number;positions:Map<number,{base:string;position:number}>}
/** Pairwise homoeolog alignment anchored at the configured corresponding target.
 * SNP evidence is from aligned genomic bases, not arbitrary nearest oligo hits. */
export function homoeologVariantMaps(refs:SangerReference[],own:number,reach=1500):VariantMap[] {
  const target=refs[own],center=target.cuts[0],radius=Math.min(1500,reach)+100;
  const left=Math.max(0,center-radius),right=Math.min(target.sequence.length,center+radius),a=target.sequence.slice(left,right);
  return refs.flatMap((other,b)=>{
    if(b===own)return [];
    const lo=Math.max(0,other.cuts[0]-radius),hi=Math.min(other.sequence.length,other.cuts[0]+radius),raw=other.sequence.slice(lo,hi);
    let best:{identity:number;score:number;reverse:boolean;columns:ReturnType<typeof alignTraceSequence>['columns']}|null=null;
    for(const reverse of [false,true]) {
      const aligned=alignTraceSequence(a,reverse?reverseComplement(raw):raw);
      const anchored=aligned.columns.some(c=>c.ref!==null&&c.query!==null&&Math.abs(c.ref+left-center)<=2&&Math.abs((reverse?hi-c.query-1:lo+c.query)-other.cuts[0])<=8);
      if(!anchored||aligned.identity<.7)continue;
      if(!best||aligned.score>best.score)best={...aligned,reverse};
    }
    const positions=new Map<number,{base:string;position:number}>();
    if(best)for(const c of best.columns)if(c.ref!==null&&c.query!==null&&/[ACGT]/.test(c.a)&&/[ACGT]/.test(c.b)&&c.a!==c.b)positions.set(left+c.ref,{base:c.b,position:best.reverse?hi-c.query-1:lo+c.query});
    return [{background:other.name,identity:best?.identity||0,positions}];
  });
}
function diagnosticVariants(hit:PrimerSite,maps:VariantMap[],terminalWindow:number):ThreePrimeVariant[]|null {
  const variants:ThreePrimeVariant[]=[];
  for(const bg of maps){
    let found:ThreePrimeVariant|null=null;
    for(let j=0;j<terminalWindow;j++){
      const pos=hit.direction==='forward'?hit.end-1-j:hit.start+j,other=bg.positions.get(pos);
      if(!other)continue;
      const base=hit.direction==='forward'?other.base:reverseComplement(other.base);
      // The primer oligo base is assigned by the caller (below) for the same 3′ offset.
      found={background:bg.background,position:pos+1,backgroundPosition:other.position+1,targetBase:'',backgroundBase:base,fromThreePrime:j+1,alignmentIdentity:bg.identity};break;
    }
    if(!found)return null;variants.push(found);
  }
  return variants;
}
/** Reject any near-binding background site whose terminal bases still match. */
export function hasTerminalCompatibleSite(oligo:string,strands:string[],mismatches:number,terminalWindow:number):boolean {
  const suffix=oligo.slice(-terminalWindow);
  for(const dna of strands)for(let hit=dna.indexOf(suffix);hit>=0;hit=dna.indexOf(suffix,hit+1)){
    const at=hit-(oligo.length-terminalWindow);if(at<0||at+oligo.length>dna.length)continue;
    let errors=0;for(let j=0;j<oligo.length;j++){if(dna[at+j]!==oligo[j])errors++;if(errors>mismatches)break;}
    if(errors<=mismatches)return true;
  }
  return false;
}

function designSpecificPairs(refs: SangerReference[], s: SangerDesignSettings,
  progress: (fraction:number,stage:string)=>void): SangerDesignResult {
  const stats={enumerated:0,qualityPassed:0,commonUnique:0,geometryPassed:0,offsetPassed:0,partnerChecked:0};
  const candidates:SangerCandidate[]=[];
  const coverage:{reference:string;candidates:number}[]=[];
  const strands=refs.map(r=>[r.sequence,reverseComplement(r.sequence)]);
  let complete=true;
  const hitCache=new Map<string,number>();
  const count=(seq:string,i:number)=>{
    const key=`${i}:${seq}`;
    if (!hitCache.has(key)) hitCache.set(key,nearBindingSiteCount(seq,strands[i],s.specificityMismatches));
    return hitCache.get(key)!;
  };
  for (let r=0;r<refs.length;r++) {
    progress(r/refs.length,`Designing separate PCR pairs for ${refs[r].name}`);
    const ref=refs[r], bounds=windows(ref,s).map(([lo,hi])=>[lo-(s.platform==='ngs'?0:s.indelRange),hi+(s.platform==='ngs'?0:s.indelRange)] as [number,number]);
    const variantMaps=s.specificityStrategy==='three-prime'?homoeologVariantMaps(refs,r,s.productMax):[];
    if (bounds.some(([lo,hi])=>lo<0 || hi>ref.sequence.length)) {coverage.push({reference:ref.name,candidates:0});continue;}
    const pool:{props:PrimerProperties;hit:PrimerSite;variants:ThreePrimeVariant[]|null}[]=[];
    // Enumerate only possible flanks. Specificity is still checked over the full inputs.
    for (const direction of ['forward','reverse'] as const) {
      const dna=strands[r][direction==='forward'?0:1];
      const oriented=bounds.map(([lo,hi])=>direction==='forward'?[lo,hi]:[dna.length-hi,dna.length-lo]);
      const first=Math.min(...oriented.map(w=>w[0])), last=Math.max(...oriented.map(w=>w[1]));
      for (let k=s.primerMin;k<=s.primerMax;k++) {
        const begin=Math.max(0,last-s.productMax), end=Math.min(dna.length-k,first-k);
        for (let at=begin;at<=end;at++) {
          stats.enumerated++;
          const props=properties(dna.slice(at,at+k),s);if (!props) continue;stats.qualityPassed++;
          const hit=site(ref,at,k,direction,s),variants=s.specificityStrategy==='three-prime'?diagnosticVariants(hit,variantMaps,s.threePrimeWindow):null;
          if(variants)for(const variant of variants)variant.targetBase=props.sequence[k-variant.fromThreePrime];
          pool.push({props,hit,variants});
        }
      }
    }
    const rank=(a:typeof pool[number],b:typeof pool[number])=>Math.max(...a.hit.cutOffsets)-Math.max(...b.hit.cutOffsets) || (s.specificityStrategy==='three-prime'?Number(!!b.variants)-Number(!!a.variants):0) || Math.abs(a.props.tm-60)-Math.abs(b.props.tm-60) || a.props.selfEnd-b.props.selfEnd || a.hit.start-b.hit.start || a.props.sequence.localeCompare(b.props.sequence);
    const byDirection={forward:pool.filter(p=>p.hit.direction==='forward').sort(rank),reverse:pool.filter(p=>p.hit.direction==='reverse').sort(rank)};
    const ownUnique=(p:typeof pool[number])=>count(p.props.sequence,r)===1;
    const specificCache=new Map<string,boolean>();
    const specific=(p:typeof pool[number])=>{
      if(specificCache.has(p.props.sequence))return specificCache.get(p.props.sequence)!;
      const passed=s.specificityStrategy==='three-prime'?!!p.variants&&refs.every((_,i)=>i===r||!hasTerminalCompatibleSite(p.props.sequence,strands[i],s.specificityMismatches,s.threePrimeWindow)):refs.every((_,i)=>i===r||count(p.props.sequence,i)===0);
      specificCache.set(p.props.sequence,passed);return passed;
    };
    const found:SangerCandidate[]=[];
    for (const direction of ['forward','reverse'] as const) {
      if (s.direction!=='both' && s.direction!==direction) continue;
      const primaries=byDirection[direction].filter(p=>coversRead([p.hit],[ref],s)&&(!s.specificBoth||s.specificityStrategy!=='three-prime'||!!p.variants));
      const partners=byDirection[direction==='forward'?'reverse':'forward'];
      stats.geometryPassed+=primaries.length;
      if (primaries.length>240 || partners.length>240) complete=false;
      const allPartners=partners.slice(0,240),diagnosticPartners=partners.filter(q=>!!q.variants).slice(0,240);
      if(partners.filter(q=>!!q.variants).length>240)complete=false;
      for (const p of primaries.slice(0,240)) {
        if (!ownUnique(p)) continue;stats.commonUnique++;
        const primarySpecific=specific(p);
        if (s.specificBoth && !primarySpecific) continue;
        let partnerOptions=0;
        const partnerShortlist=s.specificityStrategy==='three-prime'&&(s.specificBoth||!primarySpecific)?diagnosticPartners:allPartners;
        for (const q of partnerShortlist) {
          stats.partnerChecked++;
          if (Math.abs(p.props.tm-q.props.tm)>s.maxTmDifference) continue;
          const product=direction==='forward'?q.hit.end-p.hit.start:p.hit.end-q.hit.start;
          if (product<s.productMin || product>s.productMax) continue;
          const comp=complementarity(p.props.sequence,q.props.sequence);
          if (comp.end>=5 || comp.any>=9 || !ownUnique(q)) continue;
          const partnerSpecific=specific(q);
          if (s.specificBoth?!partnerSpecific:!primarySpecific&&!partnerSpecific) continue;
          if(s.platform!=='ngs'&&!coversRead([{...p.hit,extension:p.hit.extension.slice(0,product-p.props.sequence.length)}],[ref],s))continue;
          stats.offsetPassed++;
          found.push({...p.props,id:`specific-${r}-${direction}-${p.props.sequence}-${q.props.sequence}`,direction,
            referenceIndex:r,sites:[{...p.hit,extension:p.hit.extension.slice(0,product-p.props.sequence.length)}],
            minOffset:0,minTraceDifferences:0,traceComparedBases:0,indelOffsetOverlap:false,
            specificity:{primary:primarySpecific,partner:partnerSpecific,mismatches:s.specificityMismatches,strategy:s.specificityStrategy,terminalWindow:s.threePrimeWindow,primaryVariants:primarySpecific?p.variants||[]:[],partnerVariants:partnerSpecific?q.variants||[]:[]},
            partner:{...q.props,sites:[q.hit],products:[product],pairEnd:comp.end}});
          if(++partnerOptions>=3)break;
        }
        if (found.length>=100){complete=false;break;}
      }

    }
    found.sort(comparePrimerDistance);if(found.length>100)complete=false;const selected=found.slice(0,100);selected.forEach(c=>c.targetDistance=primerTargetDistance(c));
    candidates.push(...selected);coverage.push({reference:ref.name,candidates:selected.length});
  }
  progress(1,'Complete');
  annotatePrimerSafety(candidates,refs,s,stage=>progress(.99,stage));candidates.sort((a,b)=>comparePrimerSafety(a,b)||comparePrimerDistance(a,b));return {candidates,coverage,stats,complete,references:refs,settings:{...s}};
}
