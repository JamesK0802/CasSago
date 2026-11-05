import {GuideReference,GuideSettings,GuideSite,GuideCandidate,GuideResult} from '../../models/grna-design';
import {reverseComplement} from '../../utils/biology.utils';
/** Coordinates are zero-based half-open. Cut N is the boundary after input base N. */
export function scanSpCas9(sequence:string,reference:number):GuideSite[]{
 const sites:GuideSite[]=[];const n=sequence.length;
 for(const strand of ['+','-'] as const){const dna=strand==='+'?sequence:reverseComplement(sequence);for(let p=0;p+23<=n;p++){
  if(dna[p+21]!=='G'||dna[p+22]!=='G')continue;
  const spacer=dna.slice(p,p+20),pam=dna.slice(p+20,p+23);if(!/^[ACGT]{20}$/.test(spacer)||!/^[ACGT]GG$/.test(pam))continue;
  sites.push({reference,start:strand==='+'?p:n-p-20,end:strand==='+'?p+20:n-p,pamStart:strand==='+'?p+20:n-p-23,pamEnd:strand==='+'?p+23:n-p-20,pam,strand,cut:strand==='+'?p+17:n-p-17,sequence:spacer});
 }}return sites;
}
export function validateGuideInput(refs:GuideReference[],backgrounds:{name:string;sequence:string}[],s:GuideSettings){
 if(!refs.length||refs.length>8)throw new Error('Use 1–8 target references.');
 const all=[...refs,...backgrounds];if(all.some(r=>!r.name.trim())||new Set(all.map(r=>r.name)).size!==all.length)throw new Error('Give every target and background a unique name.');
 if(all.some(r=>!r.sequence||/[^ACGTN]/.test(r.sequence)))throw new Error('Use DNA references containing A/C/G/T/N.');
 if(all.some(r=>r.sequence.length>500000)||all.reduce((sum,r)=>sum+r.sequence.length,0)>2000000)throw new Error('Use ≤500 kb per reference and ≤2 Mb total, including backgrounds.');
 if(!['single','shared','specific'].includes(s.mode)||!Number.isInteger(s.mismatches)||s.mismatches<0||s.mismatches>3)throw new Error('Choose a valid target mode and 0–3 mismatches.');
 if(!Number.isFinite(s.gcMin)||!Number.isFinite(s.gcMax)||s.gcMin<0||s.gcMax>100||s.gcMin>s.gcMax)throw new Error('Use a valid GC interval within 0–100%.');
 if(!Number.isInteger(s.maxCandidates)||s.maxCandidates<1||s.maxCandidates>2000)throw new Error('Screen 1–2,000 ranked candidates.');
 const selected=refs.filter(r=>r.selected);if(s.mode==='single'?selected.length!==1:selected.length<2)throw new Error(s.mode==='single'?'Select exactly one target reference.':'Select at least two target references.');
 refs.forEach(r=>{if((r.regionStart!==null&&(!Number.isInteger(r.regionStart)||r.regionStart<1||r.regionStart>r.sequence.length))||(r.regionEnd!==null&&(!Number.isInteger(r.regionEnd)||r.regionEnd<1||r.regionEnd>r.sequence.length))||(r.regionStart!==null&&r.regionEnd!==null&&r.regionStart>r.regionEnd))throw new Error(`${r.name}: check the 1-based cut-region bounds.`);});
}
export function designGuides(refs:GuideReference[],backgrounds:{name:string;sequence:string}[],settings:GuideSettings,progress:(stage:string)=>void=()=>{}):GuideResult{
 validateGuideInput(refs,backgrounds,settings);const all=[...refs,...backgrounds];progress('Scanning both strands for NGG PAMs');const scanned=all.flatMap((r,i)=>scanSpCas9(r.sequence,i)),exact=new Map<string,GuideSite[]>();for(const site of scanned){if(!exact.has(site.sequence))exact.set(site.sequence,[]);exact.get(site.sequence)!.push(site);}
 const targets=refs.flatMap((r,i)=>r.selected?[i]:[]),pool:GuideCandidate[]=[];
 const inRegion=(h:GuideSite)=>{const r=refs[h.reference];return h.cut>=(r.regionStart??1)&&h.cut<=(r.regionEnd??r.sequence.length);};
 for(const [sequence,allSites] of exact){const gc=100*[...sequence].filter(b=>'GC'.includes(b)).length/20,polyT=sequence.includes('TTTT');if(gc<settings.gcMin||gc>settings.gcMax||settings.excludePolyT&&polyT)continue;
  const groups=targets.map(i=>allSites.filter(h=>h.reference===i));const make=(sites:GuideSite[],suffix:string)=>pool.push({id:sequence+suffix,sequence,gc,polyT,sites,otherHits:[],nearestOtherMismatch:null,exampleOtherSites:[]});
  if(settings.mode==='shared'){if(groups.every(g=>g.length&&g.some(inRegion)&&(!settings.uniqueTargets||g.length===1)))make(groups.flatMap(g=>g.filter(inRegion)),'-shared');}
  else for(let k=0;k<targets.length;k++){const g=groups[k];if(!g.length||settings.uniqueTargets&&g.length!==1)continue;const eligible=g.filter(inRegion);if(eligible.length)make(eligible,'-'+targets[k]);}
 }
 const distance=(c:GuideCandidate)=>c.sites.reduce((sum,h)=>{const r=refs[h.reference];return sum+(r.regionStart!==null||r.regionEnd!==null?Math.abs(h.cut-((r.regionStart??1)+(r.regionEnd??r.sequence.length))/2):0);},0)/c.sites.length;
 pool.sort((a,b)=>distance(a)-distance(b)||Math.abs(a.gc-50)-Math.abs(b.gc-50)||Number(a.polyT)-Number(b.polyT)||a.sequence.localeCompare(b.sequence)||a.id.localeCompare(b.id));
 // Pigeonhole retrieval: <=m substitutions guarantees one unchanged disjoint seed.
 const parts=settings.mismatches+1,index=Array.from({length:parts},()=>new Map<string,number[]>()),bounds=Array.from({length:parts},(_,i)=>[Math.floor(i*20/parts),Math.floor((i+1)*20/parts)]);
 scanned.forEach((h,i)=>bounds.forEach(([lo,hi],k)=>{const key=h.sequence.slice(lo,hi);if(!index[k].has(key))index[k].set(key,[]);index[k].get(key)!.push(i);}));
 const checked=pool.slice(0,settings.maxCandidates),candidates:GuideCandidate[]=[];
 checked.forEach((c,i)=>{if(i%50===0)progress(`Screening input-reference near sites ${i+1}/${checked.length}`);const ids=new Set<number>();bounds.forEach(([lo,hi],k)=>index[k].get(c.sequence.slice(lo,hi))?.forEach(id=>ids.add(id)));const intended=new Set(c.sites.map(h=>`${h.reference}:${h.start}:${h.strand}`));c.otherHits=Array(parts).fill(0);
  for(const id of ids){const h=scanned[id];if(intended.has(`${h.reference}:${h.start}:${h.strand}`))continue;let mm=0;for(let j=0;j<20;j++)if(h.sequence[j]!==c.sequence[j]&&++mm>settings.mismatches)break;if(mm>settings.mismatches)continue;c.otherHits[mm]++;c.nearestOtherMismatch=Math.min(c.nearestOtherMismatch??Infinity,mm);if(c.exampleOtherSites.length<20)c.exampleOtherSites.push({...h,mismatches:mm});}
  if(settings.mode!=='specific'||c.otherHits.every(n=>n===0))candidates.push(c);
 });
 candidates.sort((a,b)=>a.otherHits.reduce((s,n)=>s+n,0)-b.otherHits.reduce((s,n)=>s+n,0)||distance(a)-distance(b)||Math.abs(a.gc-50)-Math.abs(b.gc-50)||a.id.localeCompare(b.id));
 return {references:refs,backgrounds,settings:{...settings},candidates,eligible:pool.length,screened:checked.length,complete:checked.length===pool.length,sitesScanned:scanned.length};
}
