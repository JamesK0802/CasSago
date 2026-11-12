import type {LinearTerm} from './sanger-marginal-lp';
/** Replace a baseline by its algebraically identical state-flow expression.
 * The observed model still includes every departure and origin. Outside an
 * editing junction, conservation permits this exact sparse representation. */
export function rebasePhotonForms(original:LinearTerm[][],baseline:LinearTerm[][],states:LinearTerm[][]){
 return original.map((form,b)=>{const coefficients=new Map<number,number>();for(const[rows,sign]of [[form,1],[baseline[b],-1],[states[b],1]]as const)for(const[i,w]of rows)coefficients.set(i,(coefficients.get(i)||0)+sign*w);return[...coefficients].filter(([,w])=>Math.abs(w)>1e-14)as LinearTerm[];});
}
