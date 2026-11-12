import {findGrnaCutSite,resolveMatchedCutSite} from '../workers/core/classifier';
import { reverseComplement } from './biology.utils';
import { SangerReference } from '../models/sanger-design';

export function cleanReference(value: string): string {
  const seq = value.replace(/\s/g, '').toUpperCase();
  if (!seq || /[^ACGTN]/.test(seq)) throw new Error('References must contain A, C, G, T or N only.');
  return seq;
}

export function parseSangerFasta(text: string): { name: string; sequence: string }[] {
  const entries: { name: string; sequence: string }[] = [];
  let name = '', lines: string[] = [];
  const flush = () => { if (name) entries.push({ name, sequence: cleanReference(lines.join('')) }); };
  for (const line of text.replace(/^\uFEFF/, '').split(/\r?\n/)) {
    if (!line.trim()) continue;
    if (line.startsWith('>')) {
      flush(); name = line.slice(1).trim(); lines = [];
      if (!name) throw new Error('Each FASTA record needs a name.');
    } else {
      if (!name) throw new Error('Use FASTA headers beginning with >.');
      lines.push(line.trim());
    }
  }
  flush();
  if (!entries.length) throw new Error('No FASTA records found.');
  if (new Set(entries.map(e => e.name)).size !== entries.length) throw new Error('Reference names must be unique.');
  return entries;
}

/** Same matched-target cut convention as CRISPR Analysis. PAM-confirmed
 * sites take priority; a unique exact target is also usable without a PAM. */
export function guideTarget(sequence:string,input:string):{cut:number;direction:'forward'|'reverse'} {
 const guide=input.replace(/\s/g,'').toUpperCase();
 if(!/^[ACGT]{17,30}$/.test(guide))throw Error('Use a 17–30 nt target sequence, or enter cut coordinates.');
 const hits:{start:number;direction:'forward'|'reverse'}[]=[];
 for(const direction of ['forward','reverse']as const){const query=direction==='forward'?guide:reverseComplement(guide);for(let start=sequence.indexOf(query);start>=0;start=sequence.indexOf(query,start+1))hits.push({start,direction});}
 const matches=hits.flatMap(h=>{const cut=resolveMatchedCutSite(sequence,h.start,guide.length,h.direction);return cut.pam_found?[{cut:cut.cut_site,direction:h.direction}]:[];}),confirmed=[...new Set(matches.map(h=>h.cut))];
 if(confirmed.length===1)return matches[0];
 if(confirmed.length>1||hits.length>1)throw Error('Target has multiple matching sites; use explicit cut coordinates.');
 if(!hits.length)throw Error('No exact target sequence found in either strand.');
 return {cut:findGrnaCutSite(sequence,guide).cut_site,direction:hits[0].direction};
}
export function guideCut(sequence:string,input:string):number{return guideTarget(sequence,input).cut;}

/** Align input strands using an actual shared, configured guide site. Never
 * infer orientation from unrelated targets or mutate the supplied references. */
export function orientSangerReferences(refs:SangerReference[]):SangerReference[]{
 if(refs.length<2)return refs;
 const first=refs[0],anchors=(first.guides||[]).filter(Boolean).map(g=>g.replace(/\s/g,'').toUpperCase());
 return refs.map((ref,i)=>{
  if(!i)return ref;
  for(const guide of anchors){
   let a:ReturnType<typeof guideTarget>,b:ReturnType<typeof guideTarget>;
   try{a=guideTarget(first.sequence,guide);b=guideTarget(ref.sequence,guide);}catch{continue;}
   if(!first.cuts.includes(a.cut)||!ref.cuts.includes(b.cut))continue;
   if(a.direction===b.direction)return ref;
   return {...ref,sequence:reverseComplement(ref.sequence),cuts:ref.cuts.map(c=>ref.sequence.length-c),inputReversed:!ref.inputReversed};
  }
  return ref;
 });
}

export function resolveSangerReference(name: string, dna: string, targets: string, mode: 'cuts' | 'guides'): SangerReference {
  const sequence = cleanReference(dna);
  const tokens = targets.trim().split(/[\s,;]+/).filter(Boolean);
  if (!tokens.length) throw new Error(`${name || 'Reference'}: add at least one target.`);
  const cuts = mode === 'guides' ? tokens.map(g => guideCut(sequence, g)) : tokens.map(t => {
    if (!/^\d+$/.test(t)) throw new Error('Cut positions must be whole numbers.');
    return Number(t);
  });
  if (cuts.some(c => c < 0 || c > sequence.length)) throw new Error(`${name}: cut position is outside the reference.`);
  if (new Set(cuts).size !== cuts.length) throw new Error(`${name}: duplicate target positions.`);
 return { name: name.trim(), sequence, cuts,...(mode==='guides'?{guides:tokens}: {}) };
}

/** Reproducible constructed geometry example, not biological validation data. */
export function sangerDesignExample(targetCount: 1 | 2 = 1): SangerReference[] {
  let seed = 19481;
  const dna = (n: number) => Array.from({ length: n }, () => {
    seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
    return 'ACGT'[seed >>> 30];
  }).join('');
  const primer = 'CCGGCGGGTACGACCATAAC';
  const partner = 'CAGCGACTCCTTGAGTGACG';
  const shared = dna(targetCount === 2 ? 700 : 450);
  return [0, 50, 105].map((offset, i) => {
    const upstream = dna(6000+i*1700) + primer + dna(50 + offset);
    return { name: `Demo-${'ABD'[i]}`, sequence: upstream + shared + reverseComplement(partner) + dna(6500), cuts: (targetCount === 2 ? [100, 350] : [100]).map(cut => upstream.length + cut) };
  });
}
