import '@angular/compiler';
import {it,expect,vi} from 'vitest';
import {readFileSync} from 'node:fs';
import {SangerDesignPageComponent} from '../../src/app/pages/sanger-design-page/sanger-design-page.component';
import {SangerPrimerMapComponent} from '../../src/app/components/sanger-primer-map/sanger-primer-map.component';
import {designSangerPrimers} from '../../src/app/workers/core/sanger-primer-designer';
import {resolveSangerReference} from '../../src/app/utils/sanger-design.utils';
import {synthReferences,generateSangerSynth} from '../../src/app/utils/sanger-synthetic-benchmark';
import {initSynthProfiles} from '../../src/app/utils/sanger-synth-profiles';
import {chooseOppositePrimers} from '../../src/app/utils/sanger-opposite-selection';
import {selectedAmplicons} from '../../src/app/utils/sanger-design-export';
const page=()=>new SangerDesignPageComponent({markForCheck(){}} as any);
function example(){const p=page();p.loadExample();p.settings={...p.settings,tmMin:52,tmMax:68,gcMin:30,gcMax:70};const refs=p.rows.map(r=>resolveSangerReference(r.name,r.sequence,r.targets,'guides'));return {p,refs,out:designSangerPrimers(refs,{...p.settings,sharedSearchStage:0})};}
it('restores default length 20-20 and retains multiple opposite choices under one sequencing primer',()=>{
 const {p,out}=example();expect([p.settings.primerMin,p.settings.primerMax]).toEqual([20,20]);expect(p.showMap).toBe(true);expect(out.candidates).toHaveLength(1);expect(out.candidates[0].commonOppositeChoices!.length).toBeGreaterThanOrEqual(3);expect(out.searchedStage).toBe(0);
},120000);
it('generates multi-target edited synth for every retained example opposite choice',()=>{
 const {out}=example(),c=out.candidates[0];
 for(const option of c.commonOppositeChoices!){const candidate=chooseOppositePrimers(out,c,option),refs=synthReferences(out,candidate);refs.forEach(ref=>{initSynthProfiles(ref);ref.profiles!.forEach(profile=>{profile.alleles=Array.from({length:5},(_,i)=>({id:i+1,percent:10,events:[{kind:'deletion' as const,target:profile.target,offset:-2,length:i+1,sequence:''}]}));});});const data=generateSangerSynth(out,candidate,refs,921);expect(data.cycles).toBeGreaterThanOrEqual(Math.max(...candidate.sites.flatMap(h=>h.cutOffsets))+Math.ceil(out.settings.window/2));}
},120000);
it('updates products and the selected primer set when another common opposite is chosen',()=>{
 const {p,out}=example();p.result=out;p.selected=out.candidates[0];const before=p.analysisAmplicons.map(a=>a.sequence);p.chooseCommonOpposite(1);expect(p.error).toBe('');expect(p.commonChoiceIndex(p.selected!)).toBe(1);expect(p.analysisAmplicons.map(a=>a.sequence)).not.toEqual(before);expect(p.result.candidates[0]).toBe(p.selected);expect(p.showMap).toBe(true);
},120000);
it('finds exactly the requested stage, preserving 3-prime ten-base protection and per-reference opposite choices',()=>{
 const data=JSON.parse(readFileSync(new URL('../fixtures/primer-cpc-try.json',import.meta.url),'utf8')),p=page();const refs=Object.entries(data).map(([name,s])=>resolveSangerReference(name,String(s),'AATATCTCTCTATCTCCTC','guides'));
 const none=designSangerPrimers(refs,{...p.settings,sharedSearchStage:0});expect(none.candidates).toHaveLength(0);expect(none.searchedStage).toBe(0);
 const out=designSangerPrimers(refs,{...p.settings,tmMin:52,tmMax:68,gcMin:30,gcMax:70,sharedSearchStage:2}),c=out.candidates[0];expect(out.searchedStage).toBe(2);expect(c.oppositeChoices).toHaveLength(2);expect(c.oppositeChoices!.every(choices=>choices.length>1)).toBe(true);
 const map=new SangerPrimerMapComponent();map.result=out;map.candidates=[c];map.selected=c;map.ngOnChanges();const g=map.groups[0],track=g.tracks[0],cells=map.primerCells(g,track,'primary');expect(cells.flatMap((x,i)=>x.mismatch?[i]:[])).toEqual([5]);expect(cells.slice(-10).every(x=>!x.mismatch)).toBe(true);map.hoverPrimer(g,track,'primary',0);expect(map.hover?.sequence).toBe(c.sequence);expect(map.hover?.cells.filter(x=>x.mismatch)).toHaveLength(1);expect(map.hover?.cells.map(x=>x.base).join('')).toBe(out.references[0].sequence.slice(track.primary.start,track.primary.end));expect(map.hover?.cells[5].base).not.toBe(c.sequence[5]);expect(map.hover?.top).toBeLessThan(map.y(0));
 const chosen=c.oppositeChoices![0].slice(1).map(q=>{try{return chooseOppositePrimers(out,c,[q,c.oppositeChoices![1].find(v=>v.primer.sequence===c.partner!.perReference![1].sequence)!]);}catch{return null;}}).find(Boolean);expect(chosen).toBeTruthy();expect(selectedAmplicons(out,chosen!)[0].reversePrimer).not.toBe(selectedAmplicons(out,c)[0].reversePrimer);
},120000);
it('requests the next stage only when Find more is clicked and preserves existing selection',()=>{
 const workers:any[]=[];vi.stubGlobal('Worker',class{onmessage:any;postMessage=vi.fn();terminate=vi.fn();constructor(){workers.push(this);}});
 try{const p=page();p.loadExample();p.search();expect(workers[0].postMessage.mock.calls[0][0].settings.sharedSearchStage).toBe(0);const result={candidates:[],complete:true,stats:{commonUnique:0},references:[],settings:p.settings,searchedStage:0};workers[0].onmessage({data:{type:'result',result}});expect(workers).toHaveLength(1);expect(p.canFindMore).toBe(true);p.findMore();expect(workers[1].postMessage.mock.calls[0][0].settings.sharedSearchStage).toBe(1);p.cancel();}finally{vi.unstubAllGlobals();}
});
it('opens the map by default and permits closing and reopening the same primer map',()=>{const p=page(),c={id:'one',sequence:'ACGT',direction:'forward'} as any;p.selectCandidate(c);expect(p.showMap).toBe(true);p.toggleCandidateMap(c);expect(p.showMap).toBe(false);p.toggleCandidateMap(c);expect(p.showMap).toBe(true);});
it('requires explicit opposite-set selection and starts with only three return options',()=>{
 const {p,out}=example();p.result=out;p.selected=out.candidates[0];expect(p.oppositeConfirmed).toBe(false);expect(p.oppositeVisibleCount).toBe(3);expect(p.oppositeTickets.slice(0,p.oppositeVisibleCount)).toHaveLength(3);
 const map=new SangerPrimerMapComponent();map.showOpposite=false;map.result=out;map.candidates=[p.selected];map.selected=p.selected;map.ngOnChanges();expect(map.groups[0].tracks.every(t=>!t.partner)).toBe(true);
 p.showOppositeChoices=true;p.chooseOppositeTicket(1);expect(p.error).toBe('');expect(p.oppositeConfirmed).toBe(true);expect(p.showOppositeChoices).toBe(false);expect(p.selectedOppositeRows.length).toBe(p.selected!.partner!.perReference?2:1);expect(p.ticketSelected(p.oppositeTickets[1])).toBe(true);map.showOpposite=true;map.candidates=[p.selected!];map.selected=p.selected;map.ngOnChanges();expect(map.groups[0].tracks.every(t=>!!t.partner)).toBe(true);
},120000);
it('offers complete validated individual-reference sets instead of independent unvalidated combinations',()=>{
 const input=JSON.parse(readFileSync(new URL('../fixtures/primer-cpc-try.json',import.meta.url),'utf8')),p=page(),refs=Object.entries(input).map(([name,s])=>resolveSangerReference(name,String(s),'AATATCTCTCTATCTCCTC','guides'));
 const out=designSangerPrimers(refs,{...p.settings,tmMin:52,sharedSearchStage:2});p.result=out;p.selected=out.candidates[0];expect(p.oppositeTickets.length).toBeGreaterThan(1);expect(p.oppositeTickets.every(t=>!t.common&&t.partner.perReference?.length===2)).toBe(true);
 p.showOppositeChoices=true;p.chooseOppositeTicket(1);expect(p.error).toBe('');expect(p.oppositeConfirmed).toBe(true);expect(p.showOppositeChoices).toBe(false);expect(p.selectedOppositeRows.length).toBe(p.selected!.partner!.perReference?2:1);expect(()=>generateSangerSynth(out,p.selected!,synthReferences(out,p.selected!),23)).not.toThrow();
},120000);

it('keeps CPC/TRY single-mismatch candidates subject to the current Tm bounds',()=>{const p=page();p.loadCpcTrySample();expect([p.settings.tmMin,p.settings.tmMax,p.settings.gcMin,p.settings.gcMax]).toEqual([55,65,40,60]);const refs=p.rows.map(r=>resolveSangerReference(r.name,r.sequence,r.targets,'guides'));const strict=designSangerPrimers(refs,{...p.settings,sharedSearchStage:2});expect(strict.candidates).toHaveLength(0);const relaxed=designSangerPrimers(refs,{...p.settings,tmMin:52,sharedSearchStage:2});expect(relaxed.candidates.some(c=>c.sequence==='GTGAGGAGTATCGAATGGGA'&&c.sites.every(h=>h.mismatchPositions?.length===1))).toBe(true);},120000);

it('requires two mismatches on every CPC/TRY reference in the final stage',()=>{const p=page();p.loadCpcTrySample();const refs=p.rows.map(r=>resolveSangerReference(r.name,r.sequence,r.targets,'guides'));const out=designSangerPrimers(refs,{...p.settings,sharedSearchStage:3});expect(out.candidates.length).toBeGreaterThan(0);expect(out.candidates.every(c=>c.sites.every(h=>h.mismatchPositions?.length===2&&h.mismatchPositions.every(n=>n<=c.sequence.length-10)))).toBe(true);},120000);
it('restores each sequencing primer’s confirmed opposite set after switching candidates',()=>{
 const p=page();p.loadCpcTrySample();const refs=p.rows.map(r=>resolveSangerReference(r.name,r.sequence,r.targets,'guides')),out=designSangerPrimers(refs,{...p.settings,sharedSearchStage:3});p.result=out;const [first,second,third]=out.candidates;expect(third).toBeTruthy();
 p.selectCandidate(first);p.chooseOppositeTicket(1);expect(p.error).toBe('');const firstChosen=p.selected!,firstProducts=[...firstChosen.partner!.products];
 p.selectCandidate(second);expect(p.oppositeConfirmed).toBe(false);p.chooseOppositeTicket(0);expect(p.error).toBe('');const secondChosen=p.selected!;
 p.selectCandidate(third);expect(p.oppositeConfirmed).toBe(false);
 p.selectCandidate(first);expect(p.selected).toBe(firstChosen);expect(p.selected!.partner!.products).toEqual(firstProducts);expect(p.oppositeConfirmed).toBe(true);expect(p.showOppositeChoices).toBe(false);
 p.selectCandidate(second);expect(p.selected).toBe(secondChosen);expect(p.oppositeConfirmed).toBe(true);
 p.changed();p.selectCandidate(first);expect(p.oppositeConfirmed).toBe(false);
},120000);
