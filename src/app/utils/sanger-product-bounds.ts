import {PooledAmplicon} from '../models/sanger-analysis';
import {reverseComplement} from './biology.utils';
/** Product orientation is forward-primer 5'→3'. Reverse primer is entered 5'→3'.
 * This verifies submitted reference boundaries, not experimental PCR specificity. */
export function validateSangerProductBounds(product:PooledAmplicon):boolean{
 const forward=(product.forwardPrimer||'').replace(/\s/g,'').toUpperCase(),reverse=(product.reversePrimer||'').replace(/\s/g,'').toUpperCase();
 if(!forward&&!reverse)return false;
 if(!/^[ACGT]{12,60}$/.test(forward)||!/^[ACGT]{12,60}$/.test(reverse))throw Error(product.name+': provide both PCR primers, 12–60 A/C/G/T bases each, in 5′→3′ orientation.');
 if(!product.sequence.startsWith(forward)||!product.sequence.endsWith(reverseComplement(reverse)))throw Error(product.name+': the submitted amplicon must start with the forward PCR primer and end with the reverse complement of the reverse PCR primer.');
 if(forward.length+reverse.length>=product.sequence.length)throw Error(product.name+': PCR primer-binding ends overlap.');return true;
}
