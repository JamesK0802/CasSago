/** A substituted molecule consumes at least one unit of the existing departure
 * budget, regardless of its incoming phase or how many substitutions it has.
 * These are within-reference composition bounds, not new fit constraints.
 */
export function zeroLengthCompositionBounds(mass:readonly [number,number],departureBudget:number,share:number){
 if(!(share>0)||!Number.isFinite(share)||!Number.isFinite(departureBudget)||departureBudget<0)throw Error('Composition bounds require a finite nonnegative budget and positive reference share.');
 const lo=Math.max(0,mass[0]),hi=Math.min(share,mass[1]),budget=departureBudget+1e-7;
 return{unmodifiedRange:[Math.max(0,lo-budget)/share,hi/share]as[number,number],substitutionRange:[0,Math.min(hi,budget)/share]as[number,number]};
}
