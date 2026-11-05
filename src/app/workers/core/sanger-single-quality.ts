/** Net-size estimates can tolerate isolated discordant cycles; a single artifact
 * must not veto hundreds of concordant positions. All cycles remain in the fit
 * and raw R². Dense or sustained disagreement still prevents certification.
 * This gate certifies net-size fractions, never complete allele identity. */
export function singleResidualSupport(positionRmse:number[]):{passes:boolean;discordantFraction:number;maxBlockRmse:number}{
 if(!positionRmse.length||positionRmse.some(x=>!Number.isFinite(x)||x<0))return {passes:false,discordantFraction:1,maxBlockRmse:Infinity};
 const discordantFraction=positionRmse.filter(x=>x>.20).length/positionRmse.length;
 const blocks:number[]=[];for(let i=0;i+20<=positionRmse.length;i+=10)blocks.push(Math.sqrt(positionRmse.slice(i,i+20).reduce((s,x)=>s+x*x,0)/20));
 const maxBlockRmse=Math.max(0,...blocks);return {passes:discordantFraction<=.05&&maxBlockRmse<=.15,discordantFraction,maxBlockRmse};
}

/** Independent WT-template uncertainty. The edited read's Phred is never used:
 * overlapping real alleles can legitimately lower its base-call confidence.
 * All-zero/absent quality is unavailable metadata, not fabricated certainty.
 * A WT donor and its local emission context use adjacent calls, so their
 * template variance combines the four neighbouring Phred errors. Relative
 * inverse-variance weights preserve uniform-quality reads and keep every cycle. */
export function controlTemplateReliability(quality:Uint8Array|null,indices:number[],sigma:number):number[]{
 if(!quality?.some(q=>q>0))return indices.map(()=>1);
 const baseErrors=Array.from(quality,q=>Math.pow(10,-q/10)),errors=baseErrors.map((_,i)=>{let intact=1;for(let j=Math.max(0,i-2);j<=Math.min(quality.length-1,i+1);j++)intact*=1-baseErrors[j];return 1-intact;}),sorted=indices.flatMap(i=>i>=0&&i<errors.length?[errors[i]]:[]).sort((a,b)=>a-b),median=sorted.length?sorted[Math.floor(sorted.length/2)]:0,variance=Math.max(1e-10,sigma*sigma);
 return indices.map(i=>i>=0&&i<errors.length?Math.min(1,(variance+median)/(variance+errors[i])):1);
}
