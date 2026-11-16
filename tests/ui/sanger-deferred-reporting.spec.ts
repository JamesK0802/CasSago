import '@angular/compiler';
import {expect,it} from 'vitest';
import {SangerAnnotationComponent} from '../../src/app/components/sanger-annotation/sanger-annotation.component';
it('does not treat an uncomputed report envelope as evidence of an abundant group',()=>{
 const ui=new SangerAnnotationComponent();ui.threshold=2;
 const deferred={edited:true,fraction:.001,fractionRange:[0,1]as[number,number],marginal:{boundMethod:'deferred-report-envelope'}};
 expect(ui.showGroup(deferred)).toBe(false);
 expect(ui.showGroup({...deferred,fraction:.02})).toBe(true);
 expect(ui.showGroup({...deferred,marginal:{boundMethod:'certified-outer-envelope'},fractionRange:[0,.03]})).toBe(true);
 ui.views=[{rows:[deferred]}]as any;
 expect(ui.smallGroups).toHaveLength(1);expect(ui.smallMass).toBe(.001);
 expect(deferred.fractionRange).toEqual([0,1]);
});
