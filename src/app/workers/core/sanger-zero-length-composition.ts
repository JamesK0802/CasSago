/** A substituted molecule consumes at least one unit of the existing departure
 * budget, regardless of its incoming phase or how many substitutions it has.
 * These are within-reference composition bounds, not new fit constraints.
 */
export function zeroLengthCompositionBounds(mass:readonly [number,number],departureBudget:number,share:number,substitutionMass?:readonly[number,number]){
 if(!(share>0)||!Number.isFinite(share)||!Number.isFinite(departureBudget)||departureBudget<0)throw Error('Composition bounds require a finite nonnegative budget and positive reference share.');
 const lo=Math.max(0,mass[0]),hi=Math.min(share,mass[1]),budget=departureBudget+1e-7;
 // Sub carrier union and WT partition the SAME local zero-length mass.
 // Bounds remain conservative when mass and carrier amount covary; they are
 // not independent allocations or a molecular linkage claim.
 const subLo=Math.max(0,substitutionMass?substitutionMass[0]-1e-7:0),subHi=Math.min(hi,budget,substitutionMass?substitutionMass[1]+1e-7:Infinity);
 return{unmodifiedRange:[Math.max(0,lo-subHi)/share,Math.max(0,hi-subLo)/share]as[number,number],substitutionRange:[subLo/share,subHi/share]as[number,number]};
}
