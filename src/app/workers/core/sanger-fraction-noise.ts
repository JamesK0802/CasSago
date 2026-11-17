/** Delta-method mean dye-fraction variance for independent ADC dye errors.
 * Residual epsilon_b - fraction_b*sum(epsilon) carries covariance induced by
 * normalization. Inputs are photon-area variances and the SAME total units. */
export function fractionNoiseVariance(fractions:number[],areaVariances:number[],total:number){
 const variance=areaVariances.reduce((a,b)=>a+b,0);
 return fractions.reduce((n,y,b)=>n+(1-y)**2*areaVariances[b]+y*y*(variance-areaVariances[b]),0)/(4*Math.max(1e-30,total*total));
}
