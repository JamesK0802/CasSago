export interface AnnotationTarget {
  ref_sequence?: string;
  sgrna_seq?: string;
  display_sgrna_seq?: string;
  grna_start_index?: number;
  cut_site_index?: number;
  is_rc?: boolean;
  strand?: string;
}

const isNgg = (value: string): boolean => /^[ACGT]GG$/.test(value);
const isCcn = (value: string): boolean => /^CC[ACGT]$/.test(value);

/**
 * Resolve the marker used by the result annotation only.
 *
 * Analysis coordinates remain untouched. This handles both supported input
 * conventions: a spacer by itself, or a reference-oriented sequence that
 * already includes its three-base PAM.
 */
export function resolveAnnotationCutSite(target: AnnotationTarget | null | undefined): number | null {
  if (!target) return null;

  const reference = (target.ref_sequence || '').replace(/\s+/g, '').toUpperCase();
  const displayedGuide = (target.display_sgrna_seq || target.sgrna_seq || '').replace(/\s+/g, '').toUpperCase();
  const guideStart = Number(target.grna_start_index);
  const savedCut = target.cut_site_index;
  const fallback = typeof savedCut === 'number' && Number.isFinite(savedCut) ? savedCut : null;

  if (!reference || !displayedGuide || !Number.isInteger(guideStart) || guideStart < 0) return fallback;

  const guideEnd = guideStart + displayedGuide.length;
  if (guideEnd > reference.length) return fallback;

  const rightPam = reference.slice(guideEnd, guideEnd + 3);
  const leftPam = guideStart >= 3 ? reference.slice(guideStart - 3, guideStart) : '';
  const guideEndsWithPam = displayedGuide.length > 3 && isNgg(displayedGuide.slice(-3));
  const guideStartsWithPam = displayedGuide.length > 3 && isCcn(displayedGuide.slice(0, 3));

  // Prefer the strand selected during analysis when both orientations happen
  // to be possible. This changes only the annotation marker.
  const reverse = target.is_rc === true || target.strand === 'reverse';

  const forwardCut = isNgg(rightPam)
    ? guideEnd - 3
    : guideEndsWithPam
      ? guideEnd - 6
      : null;

  const reverseCut = isCcn(leftPam)
    ? guideStart + 3
    : guideStartsWithPam
      ? guideStart + 6
      : null;

  const resolved = reverse
    ? (reverseCut ?? forwardCut)
    : (forwardCut ?? reverseCut);

  return resolved !== null && resolved >= 0 && resolved <= reference.length
    ? resolved
    : fallback;
}
