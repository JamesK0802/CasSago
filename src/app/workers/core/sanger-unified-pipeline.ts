import type {SignalPathInput} from './sanger-path-common';
import {analyzeMarginalSanger} from './sanger-marginal-engine';
import type {SangerLinearOptimizer} from './sanger-discrete-phase';
/** The production sequence-group workflow is identical for every input:
 * control-qualified observations, measured-control calibration, flank length
 * proportions, continuous local base-flow bounds, independent flank validation,
 * and reporting. Selected-group sequence queries reuse that same joint model. Synthetic metadata is never an inference selector. */
export function analyzeUnifiedSequenceGroups(input:SignalPathInput,progress:(stage:string)=>void,optimizer?:SangerLinearOptimizer){
 if(!optimizer)throw Error('Initialize the common sequence optimizer before group analysis.');
 return analyzeMarginalSanger(input,progress,optimizer);
}
