/** A simultaneous residual guard for a family of measured dye channels.
 * A pointwise 3× noise tolerance need not cover hundreds of rows together.
 * The sub-Gaussian union-bound scale grows with measured rows, never with
 * genotype count or synthetic metadata. It is a model safety allowance, not
 * a claim that Sanger residuals provide an exact statistical confidence set. */
export function simultaneousResidualLimit(noiseRmse:number,channelRows:number,alpha=.05){
 if(!Number.isFinite(noiseRmse)||noiseRmse<0||!Number.isFinite(channelRows)||channelRows<1||alpha<=0||alpha>=1)throw Error('Invalid simultaneous residual guard inputs.');
 return Math.max(1e-5,noiseRmse*Math.sqrt(2*Math.log(2*channelRows/alpha)));
}
