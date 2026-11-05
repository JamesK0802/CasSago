import '@angular/compiler';
import {it,expect,vi,afterEach} from 'vitest';
import {SangerAnalysisPanelComponent} from '../../src/app/pages/sanger-analysis-panel/sanger-analysis-panel.component';
afterEach(()=>vi.unstubAllGlobals());
it('starts post-inference validation when benchmark is enabled later and keeps a single active worker',()=>{
 const created:any[]=[];vi.stubGlobal('Worker',class{postMessage=vi.fn();terminate=vi.fn();onmessage:any;onerror:any;constructor(){created.push(this);}});
 const p=new SangerAnalysisPanelComponent({markForCheck:vi.fn()}as any);p.loadDemo();p.run();
 expect(created[0].postMessage.mock.calls[0][0]).not.toHaveProperty('truth');
 const result:any={};created[0].onmessage({data:{type:'result',result}});
 expect(created[0].terminate).toHaveBeenCalledOnce();expect(p.busy).toBe(false);
 p.benchmarkTruth={truth:[],recipe:null}as any;p.toggleBenchmark(true);
 expect(created).toHaveLength(2);expect(created[1].postMessage.mock.calls[0][0]).toMatchObject({action:'benchmark',result});
 expect(created[1].postMessage.mock.calls[0][0].controlBytes).toBeInstanceOf(ArrayBuffer);
 p.chooseExactRepresentative({reference:'unknown',option:0});p.toggleBenchmark(true);
 expect(created).toHaveLength(2);expect(p.result).toBe(result);
 created[1].onmessage({data:{type:'benchmark',assessment:{certificates:[]}}});
 expect(created[1].terminate).toHaveBeenCalledOnce();expect(p.busy).toBe(false);p.ngOnDestroy();
});
