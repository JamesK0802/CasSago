import type {PathEdit} from './sanger-path-common';
/** Noise-scaled structural regularization. Inserted bases add sequence
 * parameters; deletion length and distance from a requested cut add coding
 * cost. Reporting thresholds never enter this cost. */
export function alleleComplexity(events:PathEdit[],cuts:number[]){
 return events.reduce((cost,e)=>cost+1+2*e.inserted.length+Math.log2(1+e.end-e.start)+Math.min(...cuts.map(c=>Math.min(Math.abs(e.start-c),Math.abs(e.end-c)))),0);
}
