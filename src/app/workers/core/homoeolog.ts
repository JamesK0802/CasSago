import { findGrnaCutSite, reverseComplement } from './classifier';
import { ReferenceMatcher as SequenceMatcher } from './reference-comparison';

export interface HomoeologGuideMatch {
  matched: boolean;
  exact: boolean;
  mismatches: number;
  identity: number;
  strand: 'forward' | 'reverse' | 'unknown';
  grnaStart: number;
  grnaEnd: number;
  cutSite: number;
  pam: string;
  pamFound: boolean;
  referenceSpacer: string;
  analysisGuide: string;
  error?: string;
}

const clean = (value: unknown) => String(value ?? '').replace(/\s+/g, '').toUpperCase();

function mismatchCount(a: string, b: string): number {
  let count = 0;
  for (let i = 0; i < Math.min(a.length, b.length); i++) if (a[i] !== b[i]) count++;
  return count + Math.abs(a.length - b.length);
}

/** Locate an exact guide or the best PAM-aware homoeologous protospacer. */
export function findHomoeologGuideSite(reference: string, guide: string, maxMismatchRate = 0.25): HomoeologGuideMatch {
  const ref = clean(reference);
  const query = clean(guide);
  if (!ref || !query || query.length > ref.length) {
    return {
      matched: false, exact: false, mismatches: query.length, identity: 0, strand: 'unknown',
      grnaStart: -1, grnaEnd: -1, cutSite: -1, pam: 'N/A', pamFound: false,
      referenceSpacer: '', analysisGuide: query, error: 'Guide or reference is empty'
    };
  }

  const exact = findGrnaCutSite(ref, query);
  if (exact.grna_start >= 0) {
    const spacer = ref.slice(exact.grna_start, exact.grna_end);
    return {
      matched: true, exact: true, mismatches: 0, identity: 100,
      strand: exact.strand as 'forward' | 'reverse', grnaStart: exact.grna_start,
      grnaEnd: exact.grna_end, cutSite: exact.cut_site, pam: exact.pam,
      pamFound: exact.pam_found, referenceSpacer: spacer,
      analysisGuide: exact.strand === 'reverse' ? reverseComplement(spacer) : spacer
    };
  }

  const rc = reverseComplement(query);
  const candidates: HomoeologGuideMatch[] = [];
  for (let pos = 0; pos <= ref.length - query.length; pos++) {
    const spacer = ref.slice(pos, pos + query.length);

    const forwardPam = ref.slice(pos + query.length, pos + query.length + 3);
    const forwardPamFound = forwardPam.length === 3 && forwardPam[1] === 'G' && forwardPam[2] === 'G';
    const forwardMismatches = mismatchCount(spacer, query);
    candidates.push({
      matched: true, exact: false, mismatches: forwardMismatches,
      identity: ((query.length - forwardMismatches) / query.length) * 100,
      strand: 'forward', grnaStart: pos, grnaEnd: pos + query.length,
      cutSite: pos + query.length - 3, pam: forwardPamFound ? forwardPam : 'NOT_FOUND',
      pamFound: forwardPamFound, referenceSpacer: spacer, analysisGuide: spacer
    });

    const reversePam = ref.slice(pos - 3, pos);
    const reversePamFound = reversePam.length === 3 && reversePam[0] === 'C' && reversePam[1] === 'C';
    const reverseMismatches = mismatchCount(spacer, rc);
    candidates.push({
      matched: true, exact: false, mismatches: reverseMismatches,
      identity: ((query.length - reverseMismatches) / query.length) * 100,
      strand: 'reverse', grnaStart: pos, grnaEnd: pos + query.length,
      cutSite: pos + 3, pam: reversePamFound ? reversePam : 'NOT_FOUND',
      pamFound: reversePamFound, referenceSpacer: spacer,
      analysisGuide: reverseComplement(spacer)
    });
  }

  candidates.sort((a, b) =>
    a.mismatches - b.mismatches || Number(b.pamFound) - Number(a.pamFound) || a.grnaStart - b.grnaStart
  );
  const best = candidates[0];
  const maxMismatches = Math.max(1, Math.floor(query.length * maxMismatchRate));
  if (!best || best.mismatches > maxMismatches) {
    return {
      matched: false, exact: false, mismatches: best?.mismatches ?? query.length,
      identity: best?.identity ?? 0, strand: 'unknown', grnaStart: -1, grnaEnd: -1,
      cutSite: -1, pam: 'N/A', pamFound: false, referenceSpacer: '', analysisGuide: query,
      error: `No homologous guide site within ${maxMismatches} mismatches`
    };
  }
  best.identity = Math.round(best.identity * 10) / 10;
  return best;
}

/** Project a coordinate from a similar reference onto the anchor reference. */
export function projectHomoeologPosition(reference: string, anchor: string, position: number): number | null {
  const ref = clean(reference);
  const base = clean(anchor);
  if (!ref || !base || position < 0 || position > ref.length) return null;
  if (ref === base) return Math.min(position, base.length);
  const opcodes = new SequenceMatcher(null, ref, base, false).getOpcodes();
  for (const [, i1, i2, j1, j2] of opcodes) {
    if (position < i1 || position > i2) continue;
    const sourceLength = Math.max(1, i2 - i1);
    const targetLength = j2 - j1;
    const ratio = Math.max(0, Math.min(1, (position - i1) / sourceLength));
    return Math.round(j1 + ratio * targetLength);
  }
  return position <= 0 ? 0 : base.length;
}

export function homoeologSimilarity(a: string, b: string): number {
  const left = clean(a);
  const right = clean(b);
  if (!left || !right) return 0;
  if (left === right) return 100;
  return Math.round(new SequenceMatcher(null, left, right, false).ratio() * 1000) / 10;
}
