/**
 * classifier.ts — X-Padding Classification Core for CRISPR Analysis.
 *
 * 1:1 TypeScript port of backend/core/classifier.py.
 *
 * Pipeline:
 * 1. gRNA-based coordinate alignment between read and reference window.
 * 2. Anchor search for precise inner-region extraction (indel detection).
 * 3. Cut-site coverage gate, with observed two-flank evidence for deletions.
 * 4. X-padding for unobserved terminal positions.
 * 5. X-aware anchor comparison (skip X, exact on observed).
 * 6. X-aware alignment scoring for gene classification.
 *
 * Core principle:
 *   Observed bases are evidence.
 *   X bases are unknown separators and ignored regardless of origin.
 *   Cut-site coverage cannot be inferred by a constant offset across an indel.
 *   Two unambiguous observed terminal anchors may establish a spanning molecule.
 */

import { SequenceMatcher, type Opcode } from './sequence-matcher';
import type { QualityScores } from './fastq-parser';

// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

const CUT_SITE_MIN_FLANK = 15;
const ANCHOR_LEN = 15;
const MAX_ANCHOR_ERRORS = 4;

// ─────────────────────────────────────────────────────────────────────────────
// Type Helpers & Basic Utilities
// ─────────────────────────────────────────────────────────────────────────────

function toStr(seq: string | Uint8Array): string {
  if (typeof seq === 'string') return seq;
  return new TextDecoder().decode(seq);
}

const COMP: Record<string, string> = {
  A: 'T', T: 'A', G: 'C', C: 'G', N: 'N', X: 'X',
  a: 't', t: 'a', g: 'c', c: 'g', n: 'n', x: 'x',
};

export function reverseComplement(seq: string): string {
  const s = toStr(seq).toUpperCase();
  let result = '';
  for (let i = s.length - 1; i >= 0; i--) {
    result += COMP[s[i]] || 'N';
  }
  return result;
}

export function avgPhred(qual: QualityScores | null): number {
  if (!qual || qual.length === 0) return 40.0;
  let sum = 0;
  for (let i = 0; i < qual.length; i++) sum += qual[i];
  return sum / qual.length;
}

// ─────────────────────────────────────────────────────────────────────────────
// Reference Optimization
// ─────────────────────────────────────────────────────────────────────────────

export interface CutSiteInfo {
  strand: string;
  grna_start: number;
  grna_end: number;
  cut_site: number;
  pam: string;
  pam_found: boolean;
}

export function resolveMatchedCutSite(
  reference: string,
  matchStart: number,
  matchLength: number,
  strand: 'forward' | 'reverse'
): Pick<CutSiteInfo, 'cut_site' | 'pam' | 'pam_found'> {
  const refUp = toStr(reference).toUpperCase();
  const matchEnd = matchStart + matchLength;
  if (strand === 'forward') {
    const adjacentPam = refUp.slice(matchEnd, matchEnd + 3);
    if (/^[ACGT]GG$/.test(adjacentPam)) {
      return { cut_site: matchEnd - 3, pam: adjacentPam, pam_found: true };
    }
    const includedPam = refUp.slice(matchEnd - 3, matchEnd);
    if (matchLength >= 23 && /^[ACGT]GG$/.test(includedPam)) {
      return { cut_site: matchEnd - 6, pam: includedPam, pam_found: true };
    }
    return { cut_site: matchEnd - 3, pam: 'NOT_FOUND', pam_found: false };
  }

  const adjacentPam = refUp.slice(matchStart - 3, matchStart);
  if (/^CC[ACGT]$/.test(adjacentPam)) {
    return { cut_site: matchStart + 3, pam: adjacentPam, pam_found: true };
  }
  const includedPam = refUp.slice(matchStart, matchStart + 3);
  if (matchLength >= 23 && /^CC[ACGT]$/.test(includedPam)) {
    return { cut_site: matchStart + 6, pam: includedPam, pam_found: true };
  }
  return { cut_site: matchStart + 3, pam: 'NOT_FOUND', pam_found: false };
}

export function findGrnaCutSite(reference: string, grna: string): CutSiteInfo {
  const refUp = toStr(reference).toUpperCase();
  const grnaUp = toStr(grna).toUpperCase();
  const grnaRc = reverseComplement(grnaUp);
  const refLen = refUp.length;
  const grnaLen = grnaUp.length;

  if (grnaLen === 0) {
    return { strand: 'unknown', grna_start: -1, grna_end: -1, cut_site: Math.floor(refLen / 2), pam: 'N/A', pam_found: false };
  }

  // Forward: gRNA + [NGG]
  for (let pos = 0; pos < refLen - grnaLen; pos++) {
    if (refUp.substring(pos, pos + grnaLen) === grnaUp) {
      const cut = resolveMatchedCutSite(refUp, pos, grnaLen, 'forward');
      if (cut.pam_found) {
        return { strand: 'forward', grna_start: pos, grna_end: pos + grnaLen, ...cut };
      }
    }
  }

  // Reverse: [CCN] + gRNA_RC
  for (let pos = 0; pos < refLen - grnaRc.length; pos++) {
    if (refUp.substring(pos, pos + grnaRc.length) === grnaRc) {
      const cut = resolveMatchedCutSite(refUp, pos, grnaRc.length, 'reverse');
      if (cut.pam_found) {
        return { strand: 'reverse', grna_start: pos, grna_end: pos + grnaRc.length, ...cut };
      }
    }
  }

  // A 23+ bp input may already include its PAM. Accept either reference
  // orientation and use the same rule as result annotations and viewers.
  if (grnaLen >= 23) {
    const direct = refUp.indexOf(grnaUp);
    if (direct !== -1) {
      for (const strand of ['forward', 'reverse'] as const) {
        const cut = resolveMatchedCutSite(refUp, direct, grnaLen, strand);
        if (cut.pam_found) return { strand, grna_start: direct, grna_end: direct + grnaLen, ...cut };
      }
    }
    const reverse = refUp.indexOf(grnaRc);
    if (reverse !== -1) {
      for (const strand of ['reverse', 'forward'] as const) {
        const cut = resolveMatchedCutSite(refUp, reverse, grnaLen, strand);
        if (cut.pam_found) return { strand, grna_start: reverse, grna_end: reverse + grnaLen, ...cut };
      }
    }
  }

  // Fallback
  let idx = refUp.indexOf(grnaUp);
  if (idx !== -1) {
    return { strand: 'forward', grna_start: idx, grna_end: idx + grnaLen, cut_site: idx + grnaLen - 3, pam: 'NOT_FOUND', pam_found: false };
  }
  idx = refUp.indexOf(grnaRc);
  if (idx !== -1) {
    return { strand: 'reverse', grna_start: idx, grna_end: idx + grnaRc.length, cut_site: idx + 3, pam: 'NOT_FOUND', pam_found: false };
  }

  return { strand: 'unknown', grna_start: -1, grna_end: -1, cut_site: Math.floor(refLen / 2), pam: 'N/A', pam_found: false };
}

export function getWindowBounds(reference: string, cutSite: number, windowSize: number, leftSize?: number, rightSize?: number): [number, number] {
  if (leftSize !== undefined && rightSize !== undefined) {
    const start = Math.max(0, cutSite - Math.max(0, leftSize));
    const end = Math.min(reference.length, cutSite + Math.max(0, rightSize));
    return [start, end];
  }
  const half = Math.floor(windowSize / 2);
  const start = Math.max(0, cutSite - half);
  const end = Math.min(reference.length, cutSite + half);
  return [start, end];
}

export function cutIndexInWindow(reference: string, cutSite: number, windowSize: number, leftSize?: number, rightSize?: number): number {
  const [start] = getWindowBounds(reference, cutSite, windowSize, leftSize, rightSize);
  return Math.max(0, cutSite - start);
}

export function extractWindow(reference: string, cutSite: number, windowSize: number, leftSize?: number, rightSize?: number): string {
  const [start, end] = getWindowBounds(reference, cutSite, windowSize, leftSize, rightSize);
  return reference.substring(start, end);
}

// ─────────────────────────────────────────────────────────────────────────────
// X-Padding Pipeline
// ─────────────────────────────────────────────────────────────────────────────

function findOffset(seqUp: string, refUp: string, sgrnaSeq: string): number | null {
  if (sgrnaSeq) {
    const sgrnaUp = sgrnaSeq.toUpperCase();
    const sgrnaRc = reverseComplement(sgrnaUp);
    for (const search of [sgrnaUp, sgrnaRc]) {
      const ri = seqUp.indexOf(search);
      if (ri === -1) continue;
      const rr = refUp.indexOf(search);
      if (rr === -1) continue;
      return ri - rr;
    }
  }

  // Left anchor fallback
  const leftAnchor = refUp.substring(0, ANCHOR_LEN);
  const li = seqUp.indexOf(leftAnchor);
  if (li !== -1) return li;

  // Right anchor fallback
  const rightAnchor = refUp.substring(refUp.length - ANCHOR_LEN);
  const ri2 = seqUp.indexOf(rightAnchor);
  if (ri2 !== -1) return ri2 - (refUp.length - ANCHOR_LEN);

  return null;
}

function editDistance(a: string, b: string, maxErrors: number = MAX_ANCHOR_ERRORS): number {
  if (Math.abs(a.length - b.length) > maxErrors) return maxErrors + 1;
  // Only distances within the existing error allowance are relevant. Any path
  // outside this diagonal band already needs more than maxErrors gaps.
  const limit = maxErrors + 1;
  let previous = new Int32Array(b.length + 1).fill(limit);
  let current = new Int32Array(b.length + 1);
  for (let j = 0; j <= Math.min(b.length, maxErrors); j++) previous[j] = j;
  for (let i = 1; i <= a.length; i++) {
    current.fill(limit);
    current[0] = Math.min(i, limit);
    let rowMin = current[0];
    for (let j = Math.max(1, i - maxErrors); j <= Math.min(b.length, i + maxErrors); j++) {
      const substitution = previous[j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1);
      current[j] = Math.min(previous[j] + 1, current[j - 1] + 1, substitution);
      rowMin = Math.min(rowMin, current[j]);
    }
    if (rowMin > maxErrors) return maxErrors + 1;
    [previous, current] = [current, previous];
  }
  return Math.min(previous[b.length], limit);
}

/** Locate an observed terminal anchor without letting SequenceMatcher discard
 * a valid downstream flank after an indel. Exact matching keeps this hot path
 * cheap; the existing alignment fallback handles reads without such a match. */
function exactAnchorPositions(reference: string, anchor: string): number[] {
  const positions: number[] = [];
  if (!anchor.length) return positions;
  let index = reference.indexOf(anchor);
  while (index !== -1) {
    positions.push(index);
    index = reference.indexOf(anchor, index + 1);
  }
  return positions;
}

function boundedAnchorPositions(reference: string, anchor: string): number[] {
  const exact = exactAnchorPositions(reference, anchor);
  if (exact.length || !anchor.length) return exact;
  let bestErrors = MAX_ANCHOR_ERRORS + 1;
  const positions: number[] = [];
  for (let start = 0; start <= reference.length - anchor.length; start++) {
    let errors = 0;
    let run = 0, longestRun = 0;
    for (let i = 0; i < anchor.length; i++) {
      if (reference[start + i] !== anchor[i]) {
        run = 0;
        if (++errors > MAX_ANCHOR_ERRORS) break;
      } else {
        longestRun = Math.max(longestRun, ++run);
      }
    }
    if (errors > MAX_ANCHOR_ERRORS) continue;
    // An approximate terminal coincidence cannot establish new reference
    // coverage without the same 10-bp placement support as the fallback.
    // Ordinary anchor error tolerance on already located windows is unchanged.
    if (longestRun < Math.min(10, anchor.length)) continue;
    if (errors < bestErrors) {
      bestErrors = errors;
      positions.length = 0;
    }
    if (errors === bestErrors) positions.push(start);
  }
  return positions;
}

interface AlignResult {
  fail: string | null;
  observed_read?: string;
  read_window?: string;
  qual_observed?: QualityScores | null;
  left_x?: number;
  right_x?: number;
}

/** Locate the observed endpoint in reference coordinates. Query length is
 * not reference coverage: internal insertions/deletions change their offset.
 * Use the same matching-block fallback as terminal truncation, but require
 * the established 10-bp placement support instead of one-base coincidences. */
function observedWindowEndpoints(reference: string, observation: string): {
  refStart: number; refEnd: number; obsStart: number; obsEnd: number;
} | null {
  const blocks = new SequenceMatcher(null, reference, observation, observation.length < 200)
    .getMatchingBlocks().filter(block => block[2] >= 10);
  if (!blocks.length) return null;
  const first = blocks[0], last = blocks[blocks.length - 1];
  let refStart = first[0], obsStart = first[1];
  let refEnd = last[0] + last[2], obsEnd = last[1] + last[2];
  // Preserve observed terminal mismatches. Dropping them and marking their
  // reference columns X would turn an anchor failure into fabricated missing
  // evidence. Only the genuinely absent remainder receives padding.
  if (refStart <= ANCHOR_LEN && obsStart <= ANCHOR_LEN) {
    const extension = Math.min(refStart, obsStart);
    refStart -= extension; obsStart -= extension;
  }
  const refTail = reference.length - refEnd, obsTail = observation.length - obsEnd;
  if (refTail <= ANCHOR_LEN && obsTail <= ANCHOR_LEN) {
    const extension = Math.min(refTail, obsTail);
    refEnd += extension; obsEnd += extension;
  }
  return { refStart, obsStart, refEnd, obsEnd };
}

/** A failed exact terminal seed is not evidence that the molecule lacks the
 * cut. Search only bounded terminal reference flanks, using exact 10 bp seeds
 * followed by the existing error allowance. Require unique placements on BOTH
 * sides; never infer coverage from one flank or from an X separator. */
function rescueObservedFlanks(
  seq: string, qual: QualityScores | null, reference: string, cut: number
): AlignResult | null {
  const length = 30;
  // Keep short-window assignment unchanged. The bounded 90-bp terminal search
  // is intended to rescue spanning large-indel windows, not relax local calls.
  if (reference.length < 180 || cut < 45 || reference.length - cut < 45) return null;
  type Hit = { ref: number; query: number; end: number; errors: number };
  const findFlank = (starts: number[]): Hit | null => {
    const hits: Hit[] = [];
    for (const ref of starts) {
      const anchor = reference.substring(ref, ref + length);
      if (reference.indexOf(anchor) !== ref || reference.indexOf(anchor, ref + 1) !== -1) continue;
      const candidates = new Set<number>();
      for (let seed = 0; seed < length; seed += 10) {
        const word = anchor.substring(seed, seed + 10);
        let found = seq.indexOf(word);
        // Repetitive seeds cannot supply unique positional evidence.
        let occurrences = 0;
        while (found !== -1 && occurrences++ < 16) {
          for (let shift = -MAX_ANCHOR_ERRORS; shift <= MAX_ANCHOR_ERRORS; shift++) {
            const start = found - seed + shift;
            if (start >= 0 && start + length - MAX_ANCHOR_ERRORS <= seq.length) candidates.add(start);
          }
          found = seq.indexOf(word, found + 1);
        }
        if (found !== -1) return null;
      }
      for (const query of candidates) {
        let bestErrors = MAX_ANCHOR_ERRORS + 1, end = -1;
        for (let size = length - MAX_ANCHOR_ERRORS; size <= length + MAX_ANCHOR_ERRORS; size++) {
          if (query + size > seq.length) continue;
          const errors = editDistance(anchor, seq.substring(query, query + size));
          if (errors < bestErrors || (errors === bestErrors && Math.abs(size - length) < Math.abs(end - query - length))) {
            bestErrors = errors; end = query + size;
          }
        }
        if (bestErrors <= MAX_ANCHOR_ERRORS) hits.push({ ref, query, end, errors: bestErrors });
      }
    }
    if (!hits.length) return null;
    hits.sort((a, b) => a.errors - b.errors || a.ref - b.ref || a.query - b.query);
    const best = hits[0];
    // Nearby starts describe the same noisy flank. Separate placements are
    // ambiguous even if one has a slightly better score.
    if (hits.some(hit => Math.abs((hit.query - hit.ref) - (best.query - best.ref)) > 2 * MAX_ANCHOR_ERRORS)) return null;
    return best;
  };
  const leftStarts: number[] = [], rightStarts: number[] = [];
  for (let distance = 0; distance <= 60; distance += 15) {
    if (distance + length <= cut - CUT_SITE_MIN_FLANK) leftStarts.push(distance);
    const start = reference.length - length - distance;
    if (start >= cut + CUT_SITE_MIN_FLANK) rightStarts.push(start);
  }
  const left = findFlank(leftStarts);
  if (!left) return null;
  const right = findFlank(rightStarts);
  if (!right || left.end > right.query) return null;
  const start = left.query - Math.min(left.ref, left.query);
  const leftX = left.ref - Math.min(left.ref, left.query);
  const refTail = reference.length - right.ref - length;
  const extension = Math.min(refTail, seq.length - right.end);
  const end = right.end + extension;
  const rightX = refTail - extension;
  const observed = seq.substring(start, end);
  return { fail: null, observed_read: observed,
    read_window: 'X'.repeat(leftX) + observed + 'X'.repeat(rightX),
    qual_observed: qual ? qual.slice(start, end) : null, left_x: leftX, right_x: rightX };
}

function alignReadToWindow(
  seq: string,
  qual: QualityScores | null,
  refWindow: string,
  cutIdxInWindow: number,
  sgrnaSeq: string
): AlignResult {
  const seqUp = seq.toUpperCase();
  const refUp = refWindow.toUpperCase();
  const winLen = refWindow.length;
  // A short direct-input window cannot hold two independent 15 bp anchors.
  // Split it into non-overlapping terminal anchors instead of treating the
  // entire window as both anchors (which can never form a valid pair).
  const anchorLen = Math.min(ANCHOR_LEN, Math.max(1, Math.floor(winLen / 2)));

  // Step 1: Find coordinate offset
  const offset = findOffset(seqUp, refUp, sgrnaSeq);
  if (offset === null) return { fail: 'no_alignment' };

  // A constant offset is only a preliminary coverage estimate. A large deletion
  // shifts reference coordinates after its junction and can put this inferred
  // cut outside a molecule that actually contains both observed flanks.
  const cutInRead = offset + cutIdxInWindow;
  const offsetCoversCut = cutInRead - CUT_SITE_MIN_FLANK >= 0 &&
    cutInRead + CUT_SITE_MIN_FLANK <= seq.length;

  // Step 2: Find anchors before rejecting coverage
  const leftAnchor = refUp.substring(0, anchorLen);
  const rightAnchor = refUp.substring(refUp.length - anchorLen);

  // Find all anchor positions
  const leftPositions: number[] = [];
  let idx = seqUp.indexOf(leftAnchor);
  while (idx !== -1) {
    leftPositions.push(idx);
    idx = seqUp.indexOf(leftAnchor, idx + 1);
  }

  const rightPositions: number[] = [];
  idx = seqUp.indexOf(rightAnchor);
  while (idx !== -1) {
    rightPositions.push(idx);
    idx = seqUp.indexOf(rightAnchor, idx + 1);
  }

  const canUseObservedPair = anchorLen === ANCHOR_LEN &&
    leftPositions.length === 1 && rightPositions.length === 1 &&
    cutIdxInWindow >= anchorLen && cutIdxInWindow <= winLen - anchorLen;
  // Do not enter the pair loop for repeated anchors without offset coverage.
  if (!offsetCoversCut && !canUseObservedPair) return { fail: 'no_coverage' };

  // Find best anchor pair
  const refInnerLen = winLen - 2 * anchorLen;
  let bestLi = -1, bestRi = -1, bestDiff = Infinity;
  for (const l of leftPositions) {
    for (const r of rightPositions) {
      if (r >= l + anchorLen) {
        const innerLen = r - l - anchorLen;
        const diff = Math.abs(innerLen - refInnerLen);
        if (diff < bestDiff) {
          bestDiff = diff;
          bestLi = l;
          bestRi = r;
        }
      }
    }
  }

  // Rescue only an unambiguous, ordered pair of *observed* exact anchors that
  // bracket the reference cut. One-sided/truncated evidence and X separators
  // retain the original coverage gate. Reuse the searches above: no new DP.
  const observedPairSpansCut = canUseObservedPair && bestLi !== -1 && bestRi !== -1;
  if (!offsetCoversCut && !observedPairSpansCut) {
    return { fail: 'no_coverage' };
  }

  // Step 4: Build observed_read + read_window
  let observedRead: string;
  let readWindow: string;
  let leftX: number;
  let rightX: number;
  let qualObserved: QualityScores | null;

  if (bestLi !== -1 && bestRi !== -1) {
    // CASE A: Both anchors found
    observedRead = seq.substring(bestLi, bestRi + anchorLen).toUpperCase();
    leftX = 0;
    rightX = 0;
    readWindow = observedRead;
    qualObserved = qual ? qual.slice(bestLi, bestRi + anchorLen) : null;

  } else if (leftPositions.length > 0) {
    // CASE B: Left anchor found, right truncated
    const bestLeft = leftPositions.reduce((best, x) => Math.abs(x - offset) < Math.abs(best - offset) ? x : best);
    const candidateObs = seq.substring(bestLeft).toUpperCase();
    const fixedWindow = seqUp.substring(bestLeft, bestLeft + winLen);
    if (bestLeft + winLen > seq.length || !xawareAnchorCheck(fixedWindow, refUp)) {
      // Right truncated by end of read / X-gap
      const terminalLen = Math.min(ANCHOR_LEN, candidateObs.length);
      const terminalAnchor = candidateObs.substring(candidateObs.length - terminalLen);
      const terminalPositions = boundedAnchorPositions(refUp, terminalAnchor);
      if (terminalPositions.length > 0) {
        const expectedEnd = Math.min(winLen, candidateObs.length);
        const terminalPos = terminalPositions.reduce((best, position) =>
          Math.abs(position + terminalLen - expectedEnd) < Math.abs(best + terminalLen - expectedEnd) ? position : best
        );
        const observedRefEnd = terminalPos + terminalLen;
        if (observedRefEnd < cutIdxInWindow + CUT_SITE_MIN_FLANK) {
          return { fail: 'no_coverage' };
        }
        observedRead = candidateObs;
        leftX = 0;
        rightX = Math.max(0, winLen - observedRefEnd);
        readWindow = observedRead + 'X'.repeat(rightX);
        qualObserved = qual ? qual.slice(bestLeft) : null;
      } else {
      const endpoints = observedWindowEndpoints(refUp, candidateObs);
      if (!endpoints || endpoints.refEnd < cutIdxInWindow + CUT_SITE_MIN_FLANK) {
        return { fail: 'no_coverage' };
      }
      const lastRefEnd = endpoints.refEnd, lastObsEnd = endpoints.obsEnd;

      if (lastRefEnd >= winLen) {
        observedRead = candidateObs.substring(0, lastObsEnd);
        leftX = 0;
        rightX = 0;
        readWindow = observedRead;
        qualObserved = qual ? qual.slice(bestLeft, bestLeft + lastObsEnd) : null;
      } else {
        observedRead = candidateObs.substring(0, lastObsEnd);
        leftX = 0;
        rightX = Math.max(0, winLen - lastRefEnd);
        readWindow = observedRead + 'X'.repeat(rightX);
        qualObserved = qual ? qual.slice(bestLeft, bestLeft + lastObsEnd) : null;
      }
      }
    } else {
      observedRead = seq.substring(bestLeft, Math.min(seq.length, bestLeft + winLen)).toUpperCase();
      leftX = 0;
      rightX = 0;
      readWindow = observedRead;
      qualObserved = qual ? qual.slice(bestLeft, bestLeft + winLen) : null;
    }

  } else if (rightPositions.length > 0) {
    // CASE C: Right anchor found, left truncated
    const bestRight = rightPositions.reduce((best, x) =>
      Math.abs(x - (offset + winLen - anchorLen)) < Math.abs(best - (offset + winLen - anchorLen)) ? x : best
    );
    const endPos = bestRight + anchorLen;
    const winStart = bestRight - (winLen - anchorLen);
    const fixedWindow = seqUp.substring(Math.max(0, winStart), endPos);
    if (winStart < 0 || !xawareAnchorCheck(fixedWindow, refUp)) {
      // Left truncated by start of read
      const candidateObs = seq.substring(0, endPos).toUpperCase();
      const terminalLen = Math.min(ANCHOR_LEN, candidateObs.length);
      const terminalAnchor = candidateObs.substring(0, terminalLen);
      const terminalPositions = boundedAnchorPositions(refUp, terminalAnchor);
      if (terminalPositions.length > 0) {
        const expectedStart = Math.max(0, winLen - candidateObs.length);
        const terminalPos = terminalPositions.reduce((best, position) =>
          Math.abs(position - expectedStart) < Math.abs(best - expectedStart) ? position : best
        );
        if (terminalPos > cutIdxInWindow - CUT_SITE_MIN_FLANK) {
          return { fail: 'no_coverage' };
        }
        observedRead = candidateObs;
        leftX = terminalPos;
        rightX = 0;
        readWindow = 'X'.repeat(leftX) + observedRead;
        qualObserved = qual ? qual.slice(0, endPos) : null;
      } else {
      const endpoints = observedWindowEndpoints(refUp, candidateObs);
      if (!endpoints || endpoints.refStart > cutIdxInWindow - CUT_SITE_MIN_FLANK) {
        return { fail: 'no_coverage' };
      }
      const firstRefStart = endpoints.refStart, firstObsStart = endpoints.obsStart;

      if (firstRefStart <= 0) {
        observedRead = candidateObs.substring(firstObsStart);
        leftX = 0;
        rightX = 0;
        readWindow = observedRead;
        qualObserved = qual ? qual.slice(firstObsStart, endPos) : null;
      } else {
        observedRead = candidateObs.substring(firstObsStart);
        leftX = Math.max(0, firstRefStart);
        rightX = 0;
        readWindow = 'X'.repeat(leftX) + observedRead;
        qualObserved = qual ? qual.slice(firstObsStart, endPos) : null;
      }
      }
    } else {
      observedRead = seq.substring(winStart, endPos).toUpperCase();
      leftX = 0;
      rightX = 0;
      readWindow = observedRead;
      qualObserved = qual ? qual.slice(winStart, endPos) : null;
    }

  } else {
    // CASE D: No anchors found — offset-based extraction
    const rwStart = offset;
    const rwEnd = offset + winLen;
    leftX = Math.max(0, -rwStart);
    rightX = Math.max(0, rwEnd - seq.length);
    const actualStart = Math.max(0, rwStart);
    const actualEnd = Math.min(seq.length, rwEnd);
    observedRead = seq.substring(actualStart, actualEnd).toUpperCase();
    readWindow = 'X'.repeat(leftX) + observedRead + 'X'.repeat(rightX);
    qualObserved = qual ? qual.slice(actualStart, actualEnd) : null;
  }

  return {
    fail: null,
    observed_read: observedRead,
    read_window: readWindow,
    qual_observed: qualObserved,
    left_x: leftX,
    right_x: rightX,
  };
}

function xawareAnchorCheck(readWindow: string, refWindow: string): boolean {
  const rw = readWindow.toUpperCase();
  const rf = refWindow.toUpperCase();
  if (!rw || !rf) return false;

  const leadingX = Math.min((rw.match(/^X+/)?.[0].length || 0), ANCHOR_LEN, rf.length);
  const trailingX = Math.min((rw.match(/X+$/)?.[0].length || 0), ANCHOR_LEN, Math.max(0, rf.length - leadingX));

  const anchorPasses = (refAnchor: string, observedEdge: string): boolean => {
    if (!refAnchor.length) return true;
    // Substitution-only evidence is already a valid edit-distance bound.
    // Avoid nine small DP runs for the common exact/SNP/error case without
    // changing the allowance; indels still use the existing DP fallback.
    if (observedEdge.length >= refAnchor.length) {
      let errors = 0;
      for (let i = 0; i < refAnchor.length; i++) {
        if (refAnchor[i] !== observedEdge[i] && ++errors > MAX_ANCHOR_ERRORS) break;
      }
      if (errors <= MAX_ANCHOR_ERRORS) return true;
    }
    const minLength = Math.max(1, refAnchor.length - MAX_ANCHOR_ERRORS);
    const maxLength = Math.min(observedEdge.length, refAnchor.length + MAX_ANCHOR_ERRORS);
    if (maxLength < minLength) return false;
    let best = MAX_ANCHOR_ERRORS + 1;
    for (let length = minLength; length <= maxLength; length++) {
      best = Math.min(best, editDistance(refAnchor, observedEdge.substring(0, length), MAX_ANCHOR_ERRORS));
    }
    return best <= MAX_ANCHOR_ERRORS;
  };

  const leftRef = rf.substring(leadingX, Math.min(ANCHOR_LEN, rf.length));
  const leftObserved = rw.substring(leadingX, Math.min(rw.length - trailingX, ANCHOR_LEN + MAX_ANCHOR_ERRORS));
  if (!anchorPasses(leftRef, leftObserved)) return false;

  const rightRefEnd = Math.max(0, rf.length - trailingX);
  const rightRef = rf.substring(Math.max(leadingX, rf.length - ANCHOR_LEN), rightRefEnd).split('').reverse().join('');
  const observedEnd = Math.max(leadingX, rw.length - trailingX);
  const rightObserved = rw.substring(Math.max(leadingX, observedEnd - ANCHOR_LEN - MAX_ANCHOR_ERRORS), observedEnd).split('').reverse().join('');
  return anchorPasses(rightRef, rightObserved);
}

/** X is an unobserved separator, never nucleotide evidence. */
function observedSegments(
  seq: string,
  qual: QualityScores | null,
): Array<{ seq: string; qual: QualityScores | null }> {
  const segments: Array<{ seq: string; qual: QualityScores | null }> = [];
  const seqUp = seq.toUpperCase();
  const matcher = /[^X]+/g;
  let match: RegExpExecArray | null;
  while ((match = matcher.exec(seqUp)) !== null) {
    const start = match.index;
    const end = start + match[0].length;
    segments.push({
      seq: match[0],
      qual: qual ? Array.from(qual).slice(start, end) : null,
    });
  }
  return segments;
}

// ─────────────────────────────────────────────────────────────────────────────
// Usability Filter
// ─────────────────────────────────────────────────────────────────────────────

// Map-based memoization (replaces Python's lru_cache)
const usabilityCache = new Map<string, [boolean, string, ReadResult | null]>();
const MAX_USABILITY_CACHE_ENTRIES = 8192;

function qualityCacheKey(qual: QualityScores | null): string {
  if (!qual || qual.length === 0) return '';
  let key = '';
  for (let i = 0; i < qual.length; i++) key += String.fromCharCode(qual[i] + 33);
  return key;
}

export interface ReadResult {
  fail: null;
  observed_read: string;
  read_window: string;
  qual_observed: QualityScores | null;
  left_x: number;
  right_x: number;
  is_rc: boolean;
  // Mutation calling retains both X-separated mates, but classification must
  // compare the located target window, not count off-window mate sequence as
  // an insertion after correcting the gap-scoring directions.
  classification_window?: string;
  seed_recovered?: boolean;
}

export interface ValidatedTargetWindow {
  ref_window: string;
  cut_index_in_window: number;
  result: ReadResult;
}

function isReadUsableCached(
  seq: string,
  qualTuple: QualityScores | null,
  refWindow: string,
  phredThreshold: number,
  sgrnaSeq: string,
  cutIdxInWindow: number
): [boolean, string, ReadResult | null] {
  const winLen = refWindow.length;
  if (cutIdxInWindow < 0) cutIdxInWindow = Math.floor(winLen / 2);

  const candidates: Array<{ res: AlignResult; isRc: boolean }> = [];
  const failReasons = new Set<string>();
  const rescueCandidates: Array<{ seq: string; qual: QualityScores | null; isRc: boolean }> = [];
  const consider = (res: AlignResult, isRc: boolean): void => {
    if (res.fail !== null) { failReasons.add(res.fail); return; }
    if (!xawareAnchorCheck(res.read_window!, refWindow)) { failReasons.add('no_anchor'); return; }
    if (res.qual_observed && res.qual_observed.length > 0 && avgPhred(res.qual_observed) < phredThreshold) { failReasons.add('quality'); return; }
    candidates.push({ res, isRc });
  };

  // Evaluate each contiguous observed segment independently. This makes an
  // Illumina inter-mate X guard equivalent to terminal X padding: X never
  // supplies coverage, anchor, quality, scoring, or mutation evidence.
  for (const segment of observedSegments(seq, qualTuple)) {
    const fwRes = alignReadToWindow(segment.seq, segment.qual, refWindow, cutIdxInWindow, sgrnaSeq);
    const rcSeq = reverseComplement(segment.seq);
    const rcQual = segment.qual ? Array.from(segment.qual).reverse() : null;
    const rcRes = alignReadToWindow(rcSeq, rcQual, refWindow, cutIdxInWindow, sgrnaSeq);

    for (const [res, isRc] of [[fwRes, false], [rcRes, true]] as [AlignResult, boolean][]) {
      consider(res, isRc);
      if (res.fail === 'no_alignment' || res.fail === 'no_coverage') {
        rescueCandidates.push({ seq: isRc ? rcSeq : segment.seq, qual: isRc ? rcQual : segment.qual, isRc });
      }
    }
  }

  // Keep every previously usable molecule and its orientation unchanged.
  // A rescue is needed only after ALL native orientation/segment paths fail.
  if (!candidates.length) {
    for (const candidate of rescueCandidates) {
      const rescued = rescueObservedFlanks(candidate.seq, candidate.qual, refWindow.toUpperCase(), cutIdxInWindow);
      if (rescued) consider(rescued, candidate.isRc);
    }
  }

  if (candidates.length === 0) {
    let reason: string;
    if (failReasons.has('quality')) reason = 'quality';
    else if (failReasons.has('no_anchor')) reason = 'no_anchor';
    else if (failReasons.has('no_coverage')) reason = 'no_coverage';
    else reason = 'no_alignment';
    return [false, reason, null];
  }

  // Prefer the candidate with least missing sequence, then most observed bases.
  const best = candidates.reduce((a, b) => {
    const aMissing = a.res.left_x! + a.res.right_x!;
    const bMissing = b.res.left_x! + b.res.right_x!;
    if (aMissing !== bMissing) return aMissing < bMissing ? a : b;
    return a.res.observed_read!.length >= b.res.observed_read!.length ? a : b;
  });

  const inputSegments = observedSegments(seq, qualTuple);
  const preservePairedObservation = inputSegments.length > 1;
  // The Illumina preprocessor already emits R1 + X + reverse-complemented R2
  // in genomic order. A single repetitive mate can tie in both orientations,
  // so using that mate's tie-break to flip the whole pair is incorrect.
  const orientedComposite = seq.toUpperCase();

  const result: ReadResult = {
    fail: null,
    // For an X-separated paired read, usability is still decided from the
    // strongest observed mate, but mutation calling must receive both mates.
    // Returning only `best.res.read_window` silently discarded the other mate
    // and made the result depend on analysis-window size.
    observed_read: preservePairedObservation
      ? orientedComposite.replace(/X/g, '')
      : best.res.observed_read!,
    read_window: preservePairedObservation
      ? orientedComposite
      : best.res.read_window!,
    qual_observed: preservePairedObservation ? null : (best.res.qual_observed || null),
    left_x: preservePairedObservation ? 0 : best.res.left_x!,
    right_x: preservePairedObservation ? 0 : best.res.right_x!,
    is_rc: best.isRc,
    classification_window: best.res.read_window!,
  };

  return [true, 'ok', result];
}

export function isReadUsable(
  seq: string,
  qual: QualityScores | null,
  refWindow: string,
  phredThreshold: number,
  sgrnaSeq: string = '',
  cutIdxInWin: number = -1
): [boolean, string, ReadResult | null] {
  // Use cache for repeated reads
  // Quality participates in usability. Omitting it made equal sequences with
  // different FASTQ qualities reuse whichever result happened to run first.
  const cacheKey = `${seq}|${qualityCacheKey(qual)}|${refWindow}|${phredThreshold}|${sgrnaSeq}|${cutIdxInWin}`;
  const cached = usabilityCache.get(cacheKey);
  if (cached) return [cached[0], cached[1], cached[2] ? { ...cached[2] } : null];

  const result = isReadUsableCached(seq, qual, refWindow, phredThreshold, sgrnaSeq, cutIdxInWin);
  // Only cache if the map isn't too large (prevents memory issues)
  // A cache entry retains the full read sequence. A large limit turns a
  // multi-file run into an accidental in-memory FASTQ copy, especially on
  // Safari. This still captures common duplicate reads without risking a tab
  // reload.
  if (usabilityCache.size >= MAX_USABILITY_CACHE_ENTRIES) {
    const oldest = usabilityCache.keys().next().value;
    if (oldest !== undefined) usabilityCache.delete(oldest);
  }
  usabilityCache.set(cacheKey, result);
  return [result[0], result[1], result[2] ? { ...result[2] } : null];
}

/** Illumina preprocessing needs each mate's own quality array, not the shared
 * Nanopore memoized result for another read with the same sequence. */
export function isReadUsableUncached(
  seq: string,
  qual: QualityScores | null,
  refWindow: string,
  phredThreshold: number,
  sgrnaSeq: string = '',
  cutIdxInWin: number = -1
): [boolean, string, ReadResult | null] {
  return isReadUsableCached(seq, qual, refWindow, phredThreshold, sgrnaSeq, cutIdxInWin);
}

export function clearClassifierCache(): void {
  usabilityCache.clear();
  alignmentScoreCache.clear();
  optimalIdentityCache.clear();
  referenceProjectionCache.clear();
}

// ─────────────────────────────────────────────────────────────────────────────
// Scoring (X-aware with optional Cut-Site Distance Weighting)
// ─────────────────────────────────────────────────────────────────────────────

// ─────────────────────────────────────────────────────────────────────────────
// Scoring (X-aware with optional Cut-Site Distance Weighting & Exclusion Flank)
// ─────────────────────────────────────────────────────────────────────────────

function calcSeagullWeight(d: number, maxDist: number, distanceWeight: number): number {
  if (distanceWeight <= 0.0) return 1.0;
  // Seagull / Sigmoidal S-curve centered around d = 10 bp
  // Rises steeply between d = 5 to d = 15 bp, then saturates
  const s0 = 1.0 / (1.0 + Math.exp(-0.35 * (-10)));
  const sMax = 1.0 / (1.0 + Math.exp(-0.35 * (maxDist - 10)));
  const sD = 1.0 / (1.0 + Math.exp(-0.35 * (d - 10)));
  const normSig = Math.max(0.0, Math.min(1.0, (sD - s0) / (sMax - s0)));
  return 1.0 + (distanceWeight * normSig);
}

function computeAlignmentScoreWithDynamicExclusion(
  strand: string,
  refUp: string,
  cutSitePos: number,
  maxDist: number,
  distanceWeight: number,
  exclusionFlank: number,
  verifiedOpcodes?: Opcode[]
): number {
  const cleanStrand = strand;
  const cleanRef = refUp.replace(/X/g, '');
  if (!cleanStrand || !cleanRef) return 0.0;

  const sm = new SequenceMatcher(null, cleanStrand, cleanRef, cleanRef.length < 200);
  const opcodes = verifiedOpcodes ?? sm.getOpcodes();
  const excludedRefIndices = new Set<number>();

  // 1. User-configured cut-site exclusion flank (disabled when set to 0)
  if (exclusionFlank > 0 && cutSitePos >= 0) {
    const sStart = Math.max(0, cutSitePos - exclusionFlank);
    const sEnd = Math.min(cleanRef.length - 1, cutSitePos + exclusionFlank);
    for (let r = sStart; r <= sEnd; r++) {
      excludedRefIndices.add(r);
    }
  }

  // 2. Expand the configured exclusion to the full span of an intersecting
  // indel or substitution. No mutation span is excluded at the default 0 bp.
  if (exclusionFlank > 0 && cutSitePos >= 0) {
    for (const [tag, , , j1, j2] of opcodes) {
      if (tag !== 'equal') {
        const nearCutSite = (j1 <= cutSitePos + exclusionFlank && j2 >= cutSitePos - exclusionFlank);
        if (nearCutSite) {
          for (let r = j1; r < j2; r++) {
            excludedRefIndices.add(r);
          }
        }
      }
    }
  }

  let matchedWeight = 0.0;
  let totalMaxWeight = 0.0;

  for (const [tag, i1, i2, j1, j2] of opcodes) {
    if (tag === 'equal') {
      const len = i2 - i1;
      for (let k = 0; k < len; k++) {
        const refIdx = j1 + k;
        if (excludedRefIndices.has(refIdx)) continue;
        const refChar = cleanRef[refIdx];
        if (refChar === 'N') continue;
        const d = Math.abs(refIdx - cutSitePos);
        const w = calcSeagullWeight(d, maxDist, distanceWeight);
        totalMaxWeight += w;
        matchedWeight += w;
      }
    // SequenceMatcher compares (read, reference): 'insert' consumes only
    // reference columns (a deletion in the read), while 'delete' consumes
    // only read columns (an insertion in the read). Reversing these tags
    // made both gap types add ZERO weight and could reward a worse match.
    } else if (tag === 'replace' || tag === 'insert') {
      const len = j2 - j1;
      for (let k = 0; k < len; k++) {
        const refIdx = j1 + k;
        if (excludedRefIndices.has(refIdx)) continue;
        const refChar = cleanRef[refIdx];
        if (refChar === 'N') continue;
        const d = Math.abs(refIdx - cutSitePos);
        const w = calcSeagullWeight(d, maxDist, distanceWeight);
        totalMaxWeight += w;
      }
    } else if (tag === 'delete') {
      const refIdx = j1;
      if (!excludedRefIndices.has(refIdx)) {
        const d = Math.abs(refIdx - cutSitePos);
        const w = calcSeagullWeight(d, maxDist, distanceWeight);
        totalMaxWeight += w * (i2 - i1);
      }
    }
  }

  return totalMaxWeight > 0.0 ? (matchedWeight / totalMaxWeight) : 0.0;
}

// Scores depend on sequence/settings, not FASTQ quality. Quality eligibility
// remains independently checked above; equal windows can safely reuse scores.
const alignmentScoreCache = new Map<string, number>();
const MAX_SCORE_CACHE_ENTRIES = 8192;
const optimalIdentityCache = new Map<string, number | null>();

/** Verify a fragmented greedy alignment, using the same match/column score.
 * Exact common ends do not need a matrix. Bound the remaining verification
 * work; this is NOT a biological window/indel-size filter. A large matrix
 * simply retains the existing classification instead of adding expensive DP.
 */
function optimalObservedIdentity(a: string, b: string, cut = -1, distanceWeight = 0, exclusionFlank = 0): number | null {
  const key = `${a}|${b}|${cut}|${distanceWeight}|${exclusionFlank}`;
  if (optimalIdentityCache.has(key)) return optimalIdentityCache.get(key)!;
  let prefix = 0, suffix = 0;
  while (prefix < Math.min(a.length,b.length) && a[prefix] === b[prefix]) prefix++;
  while (suffix < Math.min(a.length,b.length)-prefix &&
    a[a.length-1-suffix] === b[b.length-1-suffix]) suffix++;
  const x = a.slice(prefix,a.length-suffix), y = b.slice(prefix,b.length-suffix);
  let score: number | null = null;
  if ((x.length+1)*(y.length+1) <= 65536) {
    const directions = new Uint8Array((x.length+1)*(y.length+1));
    const width = y.length+1;
    let distances = new Uint32Array(y.length+1);
    let matches = new Uint32Array(y.length+1);
    for (let j=0;j<=y.length;j++) { distances[j]=j; directions[j]=2; }
    for (let i=1;i<=x.length;i++) {
      const nextDistances = new Uint32Array(y.length+1), nextMatches = new Uint32Array(y.length+1);
      nextDistances[0]=i;
      directions[i*width]=1;
      for (let j=1;j<=y.length;j++) {
        const equal = x[i-1]===y[j-1];
        let cost = distances[j-1]+(equal?0:1), count = matches[j-1]+(equal?1:0);
        let direction = 0;
        const deletion=distances[j]+1,insertion=nextDistances[j-1]+1;
        if (deletion<cost || (deletion===cost && matches[j]>count)) {cost=deletion;count=matches[j];direction=1;}
        if (insertion<cost || (insertion===cost && nextMatches[j-1]>count)) {cost=insertion;count=nextMatches[j-1];direction=2;}
        directions[i*width+j]=direction;
        nextDistances[j]=cost;nextMatches[j]=count;
      }
      distances=nextDistances;matches=nextMatches;
    }
    const reverse: Opcode[] = [];
    let i=x.length,j=y.length;
    while (i || j) {
      const direction=directions[i*width+j], i2=i,j2=j;
      if (direction===0) {i--;j--;}
      else if (direction===1) i--;
      else j--;
      const tag=direction===1?'delete':direction===2?'insert':x[i]===y[j]?'equal':'replace';
      reverse.push([tag,i+prefix,i2+prefix,j+prefix,j2+prefix]);
    }
    const ops: Opcode[] = [];
    const append=(op:Opcode)=>{
      const previous=ops.at(-1);
      if (previous && previous[0]===op[0] && previous[2]===op[1] && previous[4]===op[3]) {previous[2]=op[2];previous[4]=op[4];}
      else ops.push(op);
    };
    if (prefix) append(['equal',0,prefix,0,prefix]);
    for (const op of reverse.reverse()) append(op);
    if (suffix) append(['equal',a.length-suffix,a.length,b.length-suffix,b.length]);
    score=computeAlignmentScoreWithDynamicExclusion(a,b,cut,Math.max(cut,b.length-1-cut),distanceWeight,exclusionFlank,ops);
  }
  if (optimalIdentityCache.size >= MAX_SCORE_CACHE_ENTRIES) optimalIdentityCache.clear();
  optimalIdentityCache.set(key,score);
  return score;
}

/** Repair only demonstrated ranking contradictions on a common observed
 * homologous interval. Do not make sample-origin labels a classification rule.
 * Uncertain verification of an otherwise unchanged ranking keeps the legacy
 * result. Conflicting candidate crops with no verified margin remain ambiguous.
 */
function competitionScores(
  evidence: Map<ClassInfo,ReadResult>, distanceWeight: number, exclusionFlank: number
): Map<ClassInfo,number> {
  const entries = Array.from(evidence.entries());
  const scores = new Map(entries.map(([c,r]) => [c,scoreReadAgainstWindow(
    r.classification_window || r.read_window,c.ref_window,c.cut_index_in_window ?? -1,distanceWeight,exclusionFlank)]));
  if (entries.length<2) return scores;
  const first=entries[0][0], guide=first.sgrna_seq?.toUpperCase();
  if (!guide || entries.some(([c]) => c.sgrna_seq?.toUpperCase()!==guide ||
      c.ref_window.length!==first.ref_window.length || c.cut_index_in_window!==first.cut_index_in_window)) return scores;
  if (new Set(entries.map(([c])=>c.gene)).size!==entries.length) return scores;
  const physicalGuide=first.ref_window.includes(guide)?guide:reverseComplement(guide);
  if (entries.some(([c])=>c.ref_window.indexOf(physicalGuide)<0 ||
      c.ref_window.indexOf(physicalGuide)!==c.ref_window.lastIndexOf(physicalGuide))) return scores;
  const ranked=[...scores.entries()].sort((a,b)=>b[1]-a[1]);
  // Scores must not depend on the user's acceptance margin. Verification
  // repairs a strict ranking contradiction; the ordinary margin gate follows.
  if (ranked[0][1]===ranked[1][1]) return scores;
  const observations=[...new Set(entries.filter(([,r])=>!r.seed_recovered)
    .map(([,r])=>r.classification_window || r.read_window))];
  if (!observations.length || observations.some(o=>o.replace(/^X+|X+$/g,'').includes('X'))) return scores;
  const unknown=(o:string)=>(o.match(/^X+/)?.[0].length||0)+(o.match(/X+$/)?.[0].length||0);
  const minUnknown=Math.min(...observations.map(unknown));
  const selected=observations.filter(o=>unknown(o)===minUnknown);
  if (ranked[0][1]===1 && observations.length===1) return scores;
  let verified: Map<ClassInfo,number> | undefined;
  let verifiedWinner: ClassInfo | undefined;
  for (const observation of selected) {
    const left=observation.match(/^X+/)?.[0].length||0,right=observation.match(/X+$/)?.[0].length||0;
    const observed=observation.slice(left,observation.length-right);
    const comparison=new Map<ClassInfo,number>();
    for (const [c] of entries) {
      const value=optimalObservedIdentity(observed,c.ref_window.slice(left,c.ref_window.length-right),
        (c.cut_index_in_window ?? -1)-left,distanceWeight,exclusionFlank);
      if (value===null) return scores;
      comparison.set(c,value);
    }
    const optimal=[...comparison.entries()].sort((a,b)=>b[1]-a[1]);
    const clear=optimal[0][1]>optimal[1][1];
    if (!clear && observations.length===1) return scores;
    if (clear && optimal[0][0]===ranked[0][0] && observations.length===1) return scores;
    if (verifiedWinner && verifiedWinner!==optimal[0][0]) {
      return new Map(entries.map(([c])=>[c,ranked[0][1]]));
    }
    verifiedWinner=optimal[0][0];verified=comparison;
  }
  return verified || scores;
}

export function scoreReadAgainstWindow(
  read: string,
  refWindow: string,
  cutIndexInWindow: number = -1,
  distanceWeight: number = 0.0,
  exclusionFlank: number = 0
): number {
  const readUp = toStr(read).toUpperCase();
  // Explicit terminal X padding denotes unknown reference positions. Do not
  // penalize those positions as biological deletions when scoring gaps.
  const fullRef = toStr(refWindow).toUpperCase();
  const scoreKey = `${readUp}|${fullRef}|${cutIndexInWindow}|${distanceWeight}|${exclusionFlank}`;
  const cachedScore = alignmentScoreCache.get(scoreKey);
  if (cachedScore !== undefined) return cachedScore;
  const leftUnknown = Math.min(readUp.match(/^X+/)?.[0].length || 0, fullRef.length);
  const rightUnknown = Math.min(readUp.match(/X+$/)?.[0].length || 0, fullRef.length-leftUnknown);
  const refUp = fullRef.slice(leftUnknown, fullRef.length-rightUnknown);
  if (!readUp || !refUp) return 0.0;

  const cutSitePos = (cutIndexInWindow >= 0 ? cutIndexInWindow : Math.floor(fullRef.length / 2)) - leftUnknown;
  const maxDist = Math.max(cutSitePos, refUp.length - cutSitePos, 1);

  const segments = readUp.split(/X+/).filter(Boolean);
  if (!segments.length) return 0.0;

  const rcRefUp = reverseComplement(refUp);
  const rcCutSitePos = refUp.length - 1 - cutSitePos;
  const rcMaxDist = Math.max(rcCutSitePos, rcRefUp.length - rcCutSitePos, 1);
  let bestScore = 0.0;
  for (const segment of segments) {
    const fwScore = computeAlignmentScoreWithDynamicExclusion(segment, refUp, cutSitePos, maxDist, distanceWeight, exclusionFlank);
    const rcScore = computeAlignmentScoreWithDynamicExclusion(reverseComplement(segment), rcRefUp, rcCutSitePos, rcMaxDist, distanceWeight, exclusionFlank);
    bestScore = Math.max(bestScore, fwScore, rcScore);
  }

  const result = Math.min(1.0, bestScore);
  if (alignmentScoreCache.size >= MAX_SCORE_CACHE_ENTRIES) {
    const oldest = alignmentScoreCache.keys().next().value;
    if (oldest !== undefined) alignmentScoreCache.delete(oldest);
  }
  alignmentScoreCache.set(scoreKey, result);
  return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Classification
// ─────────────────────────────────────────────────────────────────────────────

export interface ClassInfo {
  gene: string;
  target: string;
  ref_window: string;
  sgrna_seq?: string;
  cut_index_in_window?: number;
}

export interface ClassificationResult {
  assigned: boolean;
  reason?: string;
  predicted_gene?: string;
  predicted_target?: string;
  top1_score?: number;
  top2_score?: number;
  gap?: number;
  debug?: any;
  validated_windows?: Record<string, ValidatedTargetWindow>;
}

const referenceProjectionCache = new Map<string, ReturnType<SequenceMatcher['getOpcodes']> | null>();

/** Project reference coordinates, not read lengths, between homologous windows.
 * Compute once per reference pair. The supplied guide must align as the same
 * exact locus in both references; never use unrelated matching fragments to
 * connect two targets. */
function projectObservedInterval(known: ClassInfo, candidate: ClassInfo, start: number, end: number): [number,number] | null {
  const guide = known.sgrna_seq?.toUpperCase();
  if (!guide || guide !== candidate.sgrna_seq?.toUpperCase()) return null;
  const a = known.ref_window.toUpperCase(), b = candidate.ref_window.toUpperCase();
  const physicalGuide = a.includes(guide) ? guide : reverseComplement(guide);
  const guideA = a.indexOf(physicalGuide), guideB = b.indexOf(physicalGuide);
  if (guideA < 0 || guideB < 0 || a.lastIndexOf(physicalGuide) !== guideA || b.lastIndexOf(physicalGuide) !== guideB) return null;
  const key = `${a}|${b}|${physicalGuide}`;
  let ops = referenceProjectionCache.get(key);
  if (ops === undefined) {
    const alignment = new SequenceMatcher(null,a,b,false).getOpcodes();
    const sharedGuide = alignment.some(([tag,i1,i2,j1]) => tag === 'equal' &&
      i1 <= guideA && i2 >= guideA + physicalGuide.length && j1 + guideA - i1 === guideB);
    ops = sharedGuide ? alignment : null;
    if (referenceProjectionCache.size >= 512) referenceProjectionCache.clear();
    referenceProjectionCache.set(key,ops);
  }
  if (!ops) return null;
  const mapBoundary = (position: number, rightBoundary: boolean): number => {
    for (const [tag,i1,i2,j1,j2] of ops!) {
      if (position < i1 || position > i2) continue;
      if (i1 === i2) { if (position === i1) return rightBoundary ? j1 : j2; continue; }
      if (tag === 'equal') return j1 + position - i1;
      if (position === i1) return j1;
      if (position === i2) return j2;
      // A boundary inside a divergent block has no defensible exact mapping.
      return -1;
    }
    return -1;
  };
  const left = mapBoundary(start,false), right = mapBoundary(end,true);
  const cut = candidate.cut_index_in_window ?? -1;
  if (left < 0 || right < left || left > cut - CUT_SITE_MIN_FLANK || right < cut + CUT_SITE_MIN_FLANK) return null;
  return [left,right];
}

/** An exact seed failure is not evidence that a competing homoeolog is absent.
 * Recheck only the already located, fully observed window, not the whole read.
 * Keep the ordinary eligibility path unchanged, including single-reference
 * analysis and reads for which no target window could be located at all. */
function classificationEvidence(
  readSeq: string, readQual: QualityScores | null, classes: ClassInfo[], phredThreshold: number
): Map<ClassInfo, ReadResult> {
  const evidence = new Map<ClassInfo, ReadResult>();
  const unseeded: ClassInfo[] = [];
  for (const c of classes) {
    const [usable, reason, res] = isReadUsable(readSeq, readQual, c.ref_window,
      phredThreshold, c.sgrna_seq || '', c.cut_index_in_window ?? -1);
    if (usable && res) evidence.set(c, res);
    else if (reason === 'no_alignment') unseeded.push(c);
  }
  // Snapshot: newly recovered competitors must not seed additional rescues.
  const located = Array.from(evidence.entries());
  for (const c of unseeded) {
    for (const [known, observation] of located) {
      if (!c.sgrna_seq || c.sgrna_seq.toUpperCase() !== known.sgrna_seq?.toUpperCase()) continue;
      const locatedWindow = observation.classification_window || observation.read_window;
      const leftX = locatedWindow.match(/^X+/)?.[0].length || 0;
      const rightX = locatedWindow.match(/X+$/)?.[0].length || 0;
      const observed = locatedWindow.slice(leftX,locatedWindow.length-rightX);
      if (!observed || observed.includes('X')) continue;
      let recoveredWindow = locatedWindow;
      if (c.ref_window.length !== known.ref_window.length ||
          (c.cut_index_in_window ?? -1) !== (known.cut_index_in_window ?? -1) || leftX || rightX) {
        const interval = projectObservedInterval(known,c,leftX,known.ref_window.length-rightX);
        if (!interval) continue;
        recoveredWindow = 'X'.repeat(interval[0]) + observed + 'X'.repeat(c.ref_window.length-interval[1]);
      }
      // The locus, observed flanks, and local quality were already validated.
      // Reuse exactly that observation instead of choosing reference-specific
      // crop boundaries or scanning unrelated parts of a long molecule.
      if (xawareAnchorCheck(recoveredWindow, c.ref_window)) {
        const preservePaired = observation.read_window !== locatedWindow;
        evidence.set(c, { ...observation,
          read_window: preservePaired ? observation.read_window : recoveredWindow,
          classification_window: recoveredWindow,
          left_x: preservePaired ? observation.left_x : (recoveredWindow.match(/^X+/)?.[0].length || 0),
          right_x: preservePaired ? observation.right_x : (recoveredWindow.match(/X+$/)?.[0].length || 0),
          seed_recovered: true });
        break;
      }
    }
  }
  return evidence;
}

export function applyClassification(
  readSeq: string,
  readQual: QualityScores | null,
  classes: ClassInfo[],
  phredThreshold: number,
  margin: number,
  cutSiteDistanceWeight: number = 0.0,
  cutSiteExclusionFlank: number = 0
): ClassificationResult {
  const eligibleClasses: Array<{ classInfo: ClassInfo; targetSeq: string }> = [];
  const evidence=classificationEvidence(readSeq, readQual, classes, phredThreshold);
  for (const [c, res] of evidence) {
    eligibleClasses.push({ classInfo: c, targetSeq: res.classification_window || res.read_window });
  }

  if (eligibleClasses.length === 0) {
    return { assigned: false, reason: 'filtered' };
  }

  const comparableScores=competitionScores(evidence,cutSiteDistanceWeight,cutSiteExclusionFlank);
  const scores: Array<[number, string, string]> = eligibleClasses.map(item => [
    comparableScores.get(item.classInfo)!,
    item.classInfo.gene,
    item.classInfo.target,
  ]);
  scores.sort((a, b) => b[0] - a[0]);

  const [top1Score, top1Gene, top1Target] = scores[0];

  if (classes.length === 1) {
    return {
      assigned: true,
      predicted_gene: top1Gene,
      predicted_target: top1Target,
      top1_score: Math.round(top1Score * 10000) / 10000,
    };
  }

  const top2Score = scores.length > 1 ? scores[1][0] : 0.0;
  const gap = top1Score - top2Score;

  if (gap >= margin) {
    return {
      assigned: true,
      predicted_gene: top1Gene,
      predicted_target: top1Target,
      top1_score: Math.round(top1Score * 10000) / 10000,
      top2_score: Math.round(top2Score * 10000) / 10000,
      gap: Math.round(gap * 10000) / 10000,
    };
  }

  return {
    assigned: false,
    reason: 'ambiguous',
    top1_score: Math.round(top1Score * 10000) / 10000,
  };
}

// ─────────────────────────────────────────────────────────────────────────────
// Gene-Level Classification
// ─────────────────────────────────────────────────────────────────────────────

export function applyGeneClassification(
  readSeq: string,
  readQual: QualityScores | null,
  geneClasses: Record<string, ClassInfo[]>,
  phredThreshold: number,
  margin: number,
  cutSiteDistanceWeight: number = 0.0,
  cutSiteExclusionFlank: number = 0
): ClassificationResult {
  const geneNames = Object.keys(geneClasses);

  // Unified gene-level classification: score each gene candidate by its best target window
  const geneScores: Array<[number, string]> = [];
  const geneDebug: Record<string, any> = {};
  let anyUsable = false;

  const evidence = classificationEvidence(readSeq, readQual, geneNames.flatMap(g => geneClasses[g]), phredThreshold);
  const comparableScores=competitionScores(evidence,cutSiteDistanceWeight,cutSiteExclusionFlank);

  for (const geneName of geneNames) {
    const targets = geneClasses[geneName];
    let bestScore = -1.0;
    let bestTarget: string | null = null;
    let geneUsable = false;
    let usableCount = 0;

    for (const t of targets) {
      const res = evidence.get(t);
      if (res) {
        geneUsable = true;
        anyUsable = true;
        usableCount++;
        const score = comparableScores.get(t)!;
        if (score > bestScore) {
          bestScore = score;
          bestTarget = t.target;
        }
      }
    }

    geneDebug[geneName] = {
      best_score: bestScore >= 0 ? Math.round(bestScore * 10000) / 10000 : null,
      best_target: bestTarget,
      usable: usableCount,
      total: targets.length,
    };

    if (geneUsable) geneScores.push([bestScore, geneName]);
  }

  if (!anyUsable) {
    return {
      assigned: false, reason: 'filtered',
      debug: {
        gene_scores: geneDebug,
        best_gene: null, second_gene: null, gap: null, margin,
        outcome: 'filtered_no_usable_window',
      },
    };
  }

  geneScores.sort((a, b) => b[0] - a[0]);
  const [top1Score, top1Gene] = geneScores[0];
  const recoveredWindows: Record<string, ValidatedTargetWindow> = {};
  for (const t of geneClasses[top1Gene]) {
    const res = evidence.get(t);
    if (res?.seed_recovered) recoveredWindows[t.target] = {
      ref_window: t.ref_window, cut_index_in_window: t.cut_index_in_window ?? -1, result: res
    };
  }
  const validatedWindows = Object.keys(recoveredWindows).length ? recoveredWindows : undefined;

  if (geneScores.length === 1) {
    return {
      assigned: true,
      predicted_gene: top1Gene,
      validated_windows: validatedWindows,
      top1_score: Math.round(top1Score * 10000) / 10000,
      debug: {
        gene_scores: geneDebug,
        best_gene: top1Gene, second_gene: null, gap: null, margin,
        outcome: 'assigned_only_one_gene_usable',
      },
    };
  }

  const [top2Score, top2Gene] = geneScores[1];
  const gap = top1Score - top2Score;

  const debugBlock: any = {
    gene_scores: geneDebug,
    best_gene: top1Gene, second_gene: top2Gene,
    gap: Math.round(gap * 10000) / 10000, margin,
  };

  if (gap >= margin) {
    debugBlock.outcome = 'assigned_margin_pass';
    return {
      assigned: true,
      predicted_gene: top1Gene,
      validated_windows: validatedWindows,
      top1_score: Math.round(top1Score * 10000) / 10000,
      top2_score: Math.round(top2Score * 10000) / 10000,
      gap: Math.round(gap * 10000) / 10000,
      debug: debugBlock,
    };
  }

  debugBlock.outcome = 'ambiguous_margin_fail';
  return {
    assigned: false, reason: 'ambiguous',
    top1_score: Math.round(top1Score * 10000) / 10000,
    top2_score: Math.round(top2Score * 10000) / 10000,
    gap: Math.round(gap * 10000) / 10000,
    debug: debugBlock,
  };
}
