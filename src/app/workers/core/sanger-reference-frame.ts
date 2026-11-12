/** Express the fitted signal coefficients in the paired-control reference
 * frame. This is an exact change of units: inverse emission scaling preserves
 * every predicted photon and every within-reference genotype proportion.
 * Fitted edited-signal coefficients are retained separately; they are not
 * silently substituted for the measured control reference contribution. */
export function controlReferenceFrame(lengths:{size:number;fraction:number}[][][],controlShares:number[],fittedShares:number[]){
 if(controlShares.some((v,h)=>!(v>0)||!(fittedShares[h]>0)))throw Error('A submitted reference lacks independently supported control or edited signal.');
 return{responseScale:controlShares.map((v,h)=>fittedShares[h]/v),lengths:lengths.map(w=>w.map((g,h)=>g.map(c=>({...c,fraction:c.fraction*controlShares[h]/fittedShares[h]}))))};
}
