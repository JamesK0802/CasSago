import {it,expect} from 'vitest';
import {comparePointFamilies} from '../probes/compare-point-families';
it('preserves the unchanged point threshold, reference/window scope and extra-family criterion',()=>{
 const baseline=[{target:1,size:0,fraction:.6},{target:1,size:-1,fraction:.4}],current=[{reference:'a',target:2,size:0,fraction:.62},{reference:'a',target:2,size:-1,fraction:.38},{reference:'b',target:2,size:-1,fraction:.9}],targets=[{current:2,baseline:1}];expect(comparePointFamilies(current,'a',baseline,targets,.52,.5).pass).toBe(true);expect(comparePointFamilies([...current,{reference:'a',target:2,size:3,fraction:.02001}],'a',baseline,targets,.5,.5).pass).toBe(false);expect(comparePointFamilies(current,'a',baseline,targets,.521,.5).pass).toBe(false);expect(comparePointFamilies(current,'a',baseline,[{current:1,baseline:1}],.5,.5).pass).toBe(false);
});
