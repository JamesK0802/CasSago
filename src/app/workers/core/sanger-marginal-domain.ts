/** Physical hypothesis domain shared by the current one-read group workflow.
 * These are existing cut-relative assumptions, not numerical time budgets,
 * genotype-count limits, GT-derived settings or synthetic-input selectors.
 * Deletions may span the cut with this additional origin slack; insertion
 * origins use the same slack. Substitution support is cut-relative per phase.
 * Changing either value changes the model and requires scientific validation.
 */
export const MARGINAL_INDEL_ORIGIN_RADIUS = 3;
export const MARGINAL_SUBSTITUTION_RADIUS = 45;

/** Cut positions are inter-base coordinates. Exactly R reference bases on each
 * side form the half-open genomic window [cut-R,cut+R); this is invariant to
 * reverse-complementing the product and must not be shifted by allele length.
 * Inserted bases have no reference donor and are already free insertion bases.
 */
export function marginalReferenceSubstitutionAllowed(donor:number|undefined,cut:number){
 return donor!==undefined&&donor>=cut-MARGINAL_SUBSTITUTION_RADIUS&&donor<cut+MARGINAL_SUBSTITUTION_RADIUS;
}
/** Project the fixed genomic support to ALL allowed read phases. The nominal
 * extent is retained separately so a coordinate-registration correction does
 * not silently recalibrate the existing noise allowance or representative prior.
 */
export function marginalTargetReadSupport(cut:number,incomingMinimum:number,incomingMaximum:number,maxDeletion:number,maxInsertion:number){
 return{first:cut-Math.max(MARGINAL_SUBSTITUTION_RADIUS,maxDeletion+MARGINAL_INDEL_ORIGIN_RADIUS)+incomingMinimum,
  last:cut+Math.max(MARGINAL_SUBSTITUTION_RADIUS-1+maxInsertion,maxInsertion+MARGINAL_INDEL_ORIGIN_RADIUS)+incomingMaximum,
  nominalLast:cut+Math.max(MARGINAL_SUBSTITUTION_RADIUS,maxInsertion+MARGINAL_INDEL_ORIGIN_RADIUS)+incomingMaximum};
}

/** The reporting/enumeration domain is exactly the fitted donor-coordinate
 * domain. Read positions after an indel cannot use an unshifted cut radius. */
export function marginalOriginReadBaseFree(q:number,origin:{start:number;incoming:number;outgoing:number;size:number},cut:number){
 const inserted=q>=origin.start+origin.incoming&&q<origin.start+origin.incoming+Math.max(0,origin.size);
 const donor=inserted?undefined:q<origin.start+origin.incoming?q-origin.incoming:q-origin.outgoing;
 return inserted||marginalReferenceSubstitutionAllowed(donor,cut);
}
