import { SequenceFeature, Topology } from '../models/sequence.model';

export function reverseComplement(seq: string): string {
  const complementMap: Record<string, string> = {
    A: 'T', T: 'A', U: 'A',
    C: 'G', G: 'C',
    R: 'Y', Y: 'R', // R = A/G, Y = C/T
    K: 'M', M: 'K', // K = G/T, M = A/C
    S: 'S', W: 'W', // S = C/G, W = A/T
    B: 'V', V: 'B', // B = C/G/T, V = A/C/G
    D: 'H', H: 'D', // D = A/G/T, H = A/C/T
    N: 'N',
    a: 't', t: 'a', u: 'a',
    c: 'g', g: 'c',
    r: 'y', y: 'r',
    k: 'm', m: 'k',
    s: 's', w: 'w',
    b: 'v', v: 'b',
    d: 'h', h: 'd',
    n: 'n',
    '-': '-'
  };

  let rc = '';
  for (let i = seq.length - 1; i >= 0; i--) {
    const char = seq[i];
    rc += complementMap[char] || char;
  }
  return rc;
}

export function complement(seq: string): string {
  const complementMap: Record<string, string> = {
    A: 'T', T: 'A', U: 'A',
    C: 'G', G: 'C',
    R: 'Y', Y: 'R', // R = A/G, Y = C/T
    K: 'M', M: 'K', // K = G/T, M = A/C
    S: 'S', W: 'W', // S = C/G, W = A/T
    B: 'V', V: 'B', // B = C/G/T, V = A/C/G
    D: 'H', H: 'D', // D = A/G/T, H = A/C/T
    N: 'N',
    a: 't', t: 'a', u: 'a',
    c: 'g', g: 'c',
    r: 'y', y: 'r',
    k: 'm', m: 'k',
    s: 's', w: 'w',
    b: 'v', v: 'b',
    d: 'h', h: 'd',
    n: 'n',
    '-': '-'
  };

  let c = '';
  for (let i = 0; i < seq.length; i++) {
    const char = seq[i];
    c += complementMap[char] || char;
  }
  return c;
}

export function getGCContent(seq: string): number {
  if (!seq || seq.length === 0) return 0;
  const gcCount = (seq.match(/[GCgcSscg]/g) || []).length;
  return (gcCount / seq.length) * 100;
}

const GENETIC_CODE: Record<string, string> = {
  ATA:'I', ATC:'I', ATT:'I', ATG:'M',
  ACA:'T', ACC:'T', ACG:'T', ACT:'T',
  AAC:'N', AAT:'N', AAA:'K', AAG:'K',
  AGC:'S', AGT:'S', AGA:'R', AGG:'R',
  CTA:'L', CTC:'L', CTG:'L', CTT:'L',
  CCA:'P', CCC:'P', CCG:'P', CCT:'P',
  CAC:'H', CAT:'H', CAA:'Q', CAG:'Q',
  CGA:'R', CGC:'R', CGG:'R', CGT:'R',
  GTA:'V', GTC:'V', GTG:'V', GTT:'V',
  GCA:'A', GCC:'A', GCG:'A', GCT:'A',
  GAC:'D', GAT:'D', GAA:'E', GAG:'E',
  GGA:'G', GGC:'G', GGG:'G', GGT:'G',
  TCA:'S', TCC:'S', TCG:'S', TCT:'S',
  TTC:'F', TTT:'F', TTA:'L', TTG:'L',
  TAC:'Y', TAT:'Y', TAA:'*', TAG:'*',
  TGC:'C', TGT:'C', TGA:'*', TGG:'W',
};

export function translateDNA(seq: string): string {
  const cleanSeq = seq.replace(/[^ATCGU]/gi, '').replace(/U/gi, 'T').toUpperCase();
  let protein = '';
  for (let i = 0; i < cleanSeq.length - 2; i += 3) {
    const codon = cleanSeq.substring(i, i + 3);
    protein += GENETIC_CODE[codon] || 'X';
  }
  return protein;
}

export interface ORF {
  start: number;
  end: number;
  frame: number;
  strand: 1 | -1;
  translation: string;
}

export function findORFs(seq: string, minAminoAcids: number = 30): ORF[] {
  const orfs: ORF[] = [];
  const seqUpper = seq.replace(/U/gi, 'T').toUpperCase();
  const rcSeq = reverseComplement(seqUpper);
  const minLen = minAminoAcids * 3;

  const findInFrame = (dna: string, strand: 1 | -1) => {
    for (let frame = 0; frame < 3; frame++) {
      let currentStart = -1;
      for (let i = frame; i <= dna.length - 3; i += 3) {
        const codon = dna.substring(i, i + 3);
        if (codon === 'ATG' && currentStart === -1) {
          currentStart = i;
        } else if ((codon === 'TAA' || codon === 'TAG' || codon === 'TGA') && currentStart !== -1) {
          if (i + 3 - currentStart >= minLen) {
            let start = currentStart;
            let end = i + 3;
            if (strand === -1) {
              const rcStart = dna.length - end;
              const rcEnd = dna.length - start;
              start = rcStart;
              end = rcEnd;
            }
            const translation = translateDNA(dna.substring(currentStart, i + 3));
            orfs.push({ start, end, frame: strand === 1 ? frame + 1 : -(frame + 1), strand, translation });
          }
          currentStart = -1; // Reset to find next ORF in same frame
        }
      }
    }
  };

  findInFrame(seqUpper, 1);
  findInFrame(rcSeq, -1);

  return orfs;
}

export function calculateTm(seq: string): number {
  const clean = seq.toUpperCase().replace(/[^ATCG]/g, '');
  if (clean.length < 14) {
    const wS = (clean.match(/[AT]/g) || []).length;
    const sS = (clean.match(/[CG]/g) || []).length;
    return (wS * 2) + (sS * 4);
  } else {
    const cg = (seq.match(/[CG]/g) || []).length;
    return 64.9 + 41 * (cg - 16.4) / seq.length;
  }
}

export function shiftFeatures(
  features: SequenceFeature[],
  editStart: number,
  deleteCount: number,
  insertLength: number,
  originalSeqLen: number,
  topology: Topology
): SequenceFeature[] {
  const diff = insertLength - deleteCount;
  if (diff === 0 && deleteCount === 0) return features;

  const newFeatures: SequenceFeature[] = [];

  for (const feat of features) {
    let newStart = feat.start;
    let newEnd = feat.end;
    let deleted = false;

    // Simplify wraparound edits for circular for now
    
    // Shift start
    if (feat.start >= editStart + deleteCount) {
      newStart += diff;
    } else if (feat.start >= editStart && feat.start < editStart + deleteCount) {
      newStart = editStart; // truncated start
    }

    // Shift end
    if (feat.end > editStart + deleteCount) {
      newEnd += diff;
    } else if (feat.end > editStart && feat.end <= editStart + deleteCount) {
      newEnd = editStart; // truncated end
    }

    if (newStart >= newEnd) {
      deleted = true; // feature completely wiped out
    }

    if (!deleted) {
      newFeatures.push({ ...feat, start: newStart, end: newEnd });
    }
  }

  return newFeatures;
}
