import {CommonPartner,OppositePrimerChoice,PrimerBindingScreen,PrimerProperties,PrimerSite,SangerCandidate,SangerDesignResult,SangerDesignSettings,SangerReference} from '../../models/sanger-design';
import {reverseComplement} from '../../utils/biology.utils';
import {complementarity,coversRead,minOffset,properties,site,traceDistance,windows} from './sanger-primer-designer';
import {screenPrimerBinding,usablePrimerRead} from './sanger-primer-safety';
import {screenSharedSangerWT} from './sanger-primer-separation';

type Direction='forward'|'reverse';
type Stage=NonNullable<SangerCandidate['sharedDesign']>['stage'];
const STAGES:Stage[]=['both-exact','one-exact','single-aligned','single-unrestricted','double'];
const LABELS=['Both primers exact and shared','Exact shared sequencing primer','One mismatch per reference, positions within 3 bp','One mismatch per reference, positions over 3 bp apart','Last resort: two mismatches per reference'];
const TERMINAL=10,MAX_QUERIES=30000,MAX_PRIMARIES=600,MAX_PAIR_WINDOWS=20000;

/** Position rules are independent of oligo quality, off-target and read guards. */
export function sharedMismatchStage(positions:number[][],length:number):Stage|null {
 if(!positions.length||positions.some(p=>p.length>2||p.some(n=>!Number.isInteger(n)||n<1||n>length-TERMINAL)||new Set(p).size!==p.length))return null;
 if(positions.every(p=>p.length===0))return 'one-exact';
 if(positions.every(p=>p.length===1))return Math.max(...positions.map(p=>p[0]))-Math.min(...positions.map(p=>p[0]))<=3?'single-aligned':'single-unrestricted';
 return positions.every(p=>p.length===2)?'double':null;
}
export function compareSharedCandidates(a:SangerCandidate,b:SangerCandidate){
 const am=a.sharedDesign!,bm=b.sharedDesign!;
 return STAGES.indexOf(am.stage)-STAGES.indexOf(bm.stage)||am.positionSpread-bm.positionSpread||am.worstFivePrimePosition-bm.worstFivePrimePosition||Math.abs(a.tm-60)-Math.abs(b.tm-60)||Math.max(...a.sites.flatMap(s=>s.cutOffsets))-Math.max(...b.sites.flatMap(s=>s.cutOffsets))||a.sequence.localeCompare(b.sequence);
}
/** Enumerate the intersection of two ungapped Hamming balls; no degenerate bases.
 * Ten terminal bases are immutable. Additional references are checked afterwards. */
function* variants(a:string,b:string,max:number,exactOne:boolean):Generator<string>{
 const limit=a.length-TERMINAL;if(limit<0||a.slice(limit)!==b.slice(limit))return;
 const diff=[...a.slice(0,limit)].flatMap((base,i)=>base!==b[i]?[i]:[]);if(diff.length>2*max)return;
 const chars=[...a];
 function* changedShared(at:number,left:number,da:number,db:number):Generator<string>{
  if(!exactOne||da===1&&db===1)yield chars.join('');
  if(!left)return;
  for(let j=at;j<limit;j++)if(a[j]===b[j])for(const base of 'ACGT')if(base!==a[j]){chars[j]=base;yield*changedShared(j+1,left-1,da+1,db+1);chars[j]=a[j];}
 }
 function* differing(at:number,da:number,db:number):Generator<string>{
  if(da>max||db>max)return;
  if(at===diff.length){yield*changedShared(0,Math.min(max-da,max-db),da,db);return;}
  const j=diff[at];for(const base of 'ACGT'){chars[j]=base;yield*differing(at+1,da+Number(base!==a[j]),db+Number(base!==b[j]));}chars[j]=a[j];
 }
 yield*differing(0,0,0);
}

export function designProgressiveSharedPrimers(refs:SangerReference[],settings:SangerDesignSettings,progress:(fraction:number,stage:string)=>void):SangerDesignResult {
 const s=settings,quality=new Map<string,PrimerProperties|null>(),scanCache=new Map<string,PrimerBindingScreen>();
 const stats={enumerated:0,qualityPassed:0,commonUnique:0,geometryPassed:0,offsetPassed:0,partnerChecked:0,separabilityRejected:0};
 const excludedCandidates:{id:string;reason:string}[]=[];
 let complete=true;
 const dna:Record<Direction,string[]>= {forward:refs.map(r=>r.sequence),reverse:refs.map(r=>reverseComplement(r.sequence))};
 const props=(seq:string)=>{if(!quality.has(seq))quality.set(seq,properties(seq,s));return quality.get(seq)!;};
 const scan=(seq:string,hits:PrimerSite[])=>{const key=seq+'|'+hits.map(h=>h.reference+':'+h.start+':'+h.end+':'+h.direction).join('|');if(!scanCache.has(key))scanCache.set(key,screenPrimerBinding(seq,refs,hits));return scanCache.get(key)!;};
 const safe=(v:PrimerBindingScreen)=>v.minimumExact&&v.minimumMismatches>=5&&!v.references.some(r=>r.hasUnknownBases);
 const geometry=(h:PrimerSite,i:number)=>coversRead([h],[refs[i]],s);
 const sharedReadReady=(sites:PrimerSite[])=>Math.max(...sites.flatMap(h=>h.cutOffsets))+Math.ceil(s.window/2)+s.indelRange<=usablePrimerRead(s,Math.min(...sites.map(h=>h.extension.length))).end;
 // Only target-local sites are indexed. The safety scan subsequently checks the
 // entire supplied references on BOTH strands, including other copies of a site.
 const locals=new Map<string,PrimerSite[]>(),indexes=new Map<string,Map<string,PrimerSite[]>>();
 const primarySites=(direction:Direction,k:number,i:number)=>{
  const key=direction+':'+k+':'+i;if(locals.has(key))return locals.get(key)!;
  const r=refs[i],oriented=dna[direction][i],cuts=r.cuts.map(c=>direction==='forward'?c:oriented.length-c),lo=Math.max(0,Math.max(...cuts)-s.readEnd-k),hi=Math.min(oriented.length-k,Math.min(...cuts)-k-122),out:PrimerSite[]=[];
  for(let at=lo;at<=hi;at++){const h=site(r,at,k,direction,s);if(geometry(h,i)&&!/[^ACGT]/.test(oriented.slice(at,at+k)))out.push(h);}
  locals.set(key,out);return out;
 };
 const queryAt=(h:PrimerSite,i:number)=>h.direction==='forward'?refs[i].sequence.slice(h.start,h.end):reverseComplement(refs[i].sequence.slice(h.start,h.end));
 const indexFor=(direction:Direction,k:number,i:number)=>{
  const key=direction+':'+k+':'+i;if(!indexes.has(key)){const map=new Map<string,PrimerSite[]>();for(const h of primarySites(direction,k,i)){const seed=queryAt(h,i).slice(-TERMINAL);map.set(seed,[...(map.get(seed)||[]),h]);}indexes.set(key,map);}return indexes.get(key)!;
 };
 const compatible=(seq:string,direction:Direction,i:number,max:number)=>{
  const found:PrimerSite[]=[];for(const h of indexFor(direction,seq.length,i).get(seq.slice(-TERMINAL))||[]){const target=queryAt(h,i),positions:number[]=[];for(let j=0;j<seq.length;j++)if(seq[j]!==target[j])positions.push(j+1);if(positions.length<=max)found.push({...h,mismatchPositions:positions});if(found.length>1)break;}return found;
 };
 const partnerPools=new Map<string,{props:PrimerProperties;hit:PrimerSite}[]>();
 const partnerPool=(direction:Direction,i:number)=>{
  const key=direction+':'+i;if(partnerPools.has(key))return partnerPools.get(key)!;
  const out:{props:PrimerProperties;hit:PrimerSite}[]=[],r=refs[i],bounds=windows(r,s),oriented=dna[direction][i];
  for(let k=s.primerMin;k<=s.primerMax;k++){
   const genomicLo=direction==='reverse'?Math.max(...bounds.map(w=>w[1]))+s.indelRange:Math.max(0,Math.max(...r.cuts)-s.productMax);
   const genomicHi=direction==='reverse'?Math.min(r.sequence.length-k,Math.min(...r.cuts)+s.productMax):Math.min(...bounds.map(w=>w[0]))-s.indelRange-k;
   const lo=direction==='forward'?genomicLo:Math.max(0,oriented.length-genomicHi-k),hi=direction==='forward'?genomicHi:Math.min(oriented.length-k,oriented.length-genomicLo-k);
   for(let at=lo;at<=hi;at++){const seq=oriented.slice(at,at+k),p=props(seq);if(!p)continue;out.push({props:p,hit:site(r,at,k,direction,s)});}
  }
  partnerPools.set(key,out);return out;
 };
 const finish=(primary:SangerCandidate,items:{props:PrimerProperties;hit:PrimerSite}[],common:boolean):SangerCandidate|null=>{
  const products=items.map((v,i)=>primary.direction==='forward'?v.hit.end-primary.sites[i].start:primary.sites[i].end-v.hit.start);
  const partner:CommonPartner={...items[0].props,...(!common?{perReference:items.map(v=>v.props)}:{}),sites:items.map((v,i)=>({...v.hit,extension:v.hit.extension.slice(0,products[i]-v.props.sequence.length)})),products,pairEnd:Math.max(...items.map(v=>complementarity(primary.sequence,v.props.sequence).end))};
  const candidate:SangerCandidate={...primary,partner,sites:primary.sites.map((h,i)=>({...h,extension:h.extension.slice(0,products[i]-primary.sequence.length)}))};
  if(!coversRead(candidate.sites,refs,s)||!sharedReadReady(candidate.sites))return null;
  const range=usablePrimerRead(s,Math.min(...candidate.sites.map(h=>h.extension.length))),wt=screenSharedSangerWT(candidate.sites,{readStart:range.start,readEnd:range.end});
  if(!wt.passed){stats.separabilityRejected++;if(excludedCandidates.length<100)excludedCandidates.push({id:candidate.id,reason:wt.reason});return null;}
  candidate.wtSeparation={minimumDifferences:wt.minimumDifferences,comparedBases:wt.comparedBases,condition:wt.condition,warnings:wt.warnings};
  [candidate.minTraceDifferences,candidate.traceComparedBases]=traceDistance(candidate.sites,s);
  const screens=common?[scan(partner.sequence,partner.sites)]:items.map((v,i)=>scan(v.props.sequence,[partner.sites[i]]));
  if(screens.some(v=>!safe(v)))return null;
  // Show the least separated partner screen, while retaining per-primer scans.
  candidate.primerSafety={primary:scan(candidate.sequence,candidate.sites),partner:screens.reduce((a,b)=>a.minimumMismatches<b.minimumMismatches?a:b)};
  candidate.partnerSafety=screens;
  candidate.id+='|'+items.map(v=>v.props.sequence+':'+v.hit.start).join('|');
  return candidate;
 };
 const pairPrimary=(c:SangerCandidate,commonOnly:boolean,individualOnly=false)=>{
  const direction:Direction=c.direction==='forward'?'reverse':'forward';
  const pools=refs.map((_,i)=>partnerPool(direction,i).filter(v=>{
   const product=c.direction==='forward'?v.hit.end-c.sites[i].start:c.sites[i].end-v.hit.start;
   if(product<s.productMin||product>s.productMax||Math.abs(v.props.tm-c.tm)>s.maxTmDifference)return false;
   const pair=complementarity(c.sequence,v.props.sequence);return pair.end<5&&pair.any<9&&geometry({...c.sites[i],extension:c.sites[i].extension.slice(0,product-c.sequence.length)},i)&&Math.max(...c.sites.flatMap(h=>h.cutOffsets))+Math.ceil(s.window/2)+s.indelRange<=usablePrimerRead(s,product-c.sequence.length).end;
  }));
  const maps=pools.map(pool=>{const map=new Map<string,typeof pool>();for(const v of pool)map.set(v.props.sequence,[...(map.get(v.props.sequence)||[]),v]);return map;});
  // Common, exact opposite primer always gets first preference.
  const shared=pools[0].filter(v=>maps.every(m=>m.get(v.props.sequence)?.length===1)).sort((a,b)=>Math.abs(a.props.tm-c.tm)-Math.abs(b.props.tm-c.tm));
  const commonCandidates:SangerCandidate[]=[];
  for(const v of individualOnly?[]:shared){const items=maps.map(m=>m.get(v.props.sequence)![0]),candidate=finish(c,items,true);if(candidate)commonCandidates.push(candidate);if(commonCandidates.length>=30){if(shared.length>30)complete=false;break;}}
  if(commonCandidates.length){const candidate=commonCandidates[0];candidate.commonOppositeChoices=commonCandidates.map(v=>v.partner!);return candidate;}
  if(commonOnly)return null;
  const options=pools.map(pool=>pool.filter(v=>safe(scan(v.props.sequence,[v.hit]))).sort((a,b)=>Math.abs(a.props.tm-c.tm)-Math.abs(b.props.tm-c.tm)||a.props.sequence.localeCompare(b.props.sequence)));
  if(options.some(p=>!p.length))return null;
  // Best-first full sets vary choices across every reference, rather than
  // exhausting the final reference before offering another choice for the first.
  const sets:SangerCandidate[]=[],queue:number[][]=[refs.map(()=>0)],seen=new Set<string>([refs.map(()=>0).join(',')]);
  let attempts=0;
  while(queue.length&&attempts++<300&&sets.length<30){
   queue.sort((a,b)=>a.reduce((n,v)=>n+v,0)-b.reduce((n,v)=>n+v,0)||a.join(',').localeCompare(b.join(',')));
   const indices=queue.shift()!,items=indices.map((at,i)=>options[i][at]);
   let compatible=true;for(let i=0;i<items.length;i++)for(let j=0;j<i;j++){const pair=complementarity(items[i].props.sequence,items[j].props.sequence);if(pair.end>=5||pair.any>=9)compatible=false;}
   if(compatible){const candidate=finish(c,items,false);if(candidate)sets.push(candidate);}
   for(let i=0;i<indices.length;i++){const next=[...indices];next[i]++;const key=next.join(',');if(next[i]<Math.min(40,options[i].length)&&!seen.has(key)){seen.add(key);queue.push(next);}}
  }
  if(attempts>=300&&queue.length)complete=false;
  const candidate=sets[0]||null;
  if(candidate)candidate.oppositeSetChoices=sets.map(c=>({partner:c.partner!,screens:c.partnerSafety!}));
  if(candidate)candidate.oppositeChoices=options.map((list,i)=>{const chosen=list.find(v=>v.props.sequence===candidate.partner!.perReference![i].sequence&&v.hit.start===candidate.partner!.sites[i].start)!,short=list.slice(0,30);if(!short.includes(chosen))short.push(chosen);return short.map(v=>({primer:v.props,site:v.hit,product:c.direction==='forward'?v.hit.end-c.sites[i].start:c.sites[i].end-v.hit.start,safety:scan(v.props.sequence,[v.hit])}));});
  return candidate;
 };
 const exactPrimaries:SangerCandidate[]=[];
 const build=(seq:string,direction:Direction,max:number,wanted:Stage|null)=>{
  stats.enumerated++;const p=props(seq);if(!p)return null;stats.qualityPassed++;
  const hits=refs.map((_,i)=>compatible(seq,direction,i,max));if(hits.some(h=>h.length!==1))return null;stats.commonUnique++;
  const sites=hits.map(h=>h[0]),positions=sites.map(h=>h.mismatchPositions!),stage=sharedMismatchStage(positions,seq.length);
  if(!stage||wanted&&stage!==wanted)return null;stats.geometryPassed++;
  const offset=minOffset(sites);if(offset<s.minOffset)return null;stats.offsetPassed++;
  const safety=scan(seq,sites);if(!safe(safety))return null;
  const flat=positions.flat(),spread=positions.every(p=>p.length===1)?Math.max(...flat)-Math.min(...flat):0;
  const candidate:SangerCandidate={...p,id:direction+'|'+seq,direction,sites,minOffset:offset,minTraceDifferences:0,traceComparedBases:0,indelOffsetOverlap:offset<=2*s.indelRange,partner:null,sharedDesign:{stage,positionSpread:spread,worstFivePrimePosition:Math.max(0,...flat),terminalExact:10},primerSafety:{primary:safety,partner:null}};
  return candidate;
 };
 let singlePrimaries:SangerCandidate[]|null=null;
 const mismatchPrimaries=(max:1|2)=>{
  const out:SangerCandidate[]=[],seen=new Set<string>();let queries=0,pairs=0,stopped=false;
  outer:for(let k=s.primerMin;k<=s.primerMax;k++)for(const direction of ['forward','reverse'] as const){
   if(s.direction!=='both'&&s.direction!==direction)continue;
   for(const a of primarySites(direction,k,0))for(const b of indexFor(direction,k,1).get(queryAt(a,0).slice(-TERMINAL))||[]){
    if(++pairs>MAX_PAIR_WINDOWS){stopped=true;break outer;}
    for(const seq of variants(queryAt(a,0),queryAt(b,1),max,max===1)){
     const key=direction+'|'+seq;if(seen.has(key))continue;seen.add(key);
     if(++queries>MAX_QUERIES){stopped=true;break outer;}
     const c=build(seq,direction,max,max===2?'double':null);if(c)out.push(c);
    }
   }
  }
  if(stopped)complete=false;return out.sort(compareSharedCandidates);
 };
 const requested=s.sharedSearchStage;
 const requestedPhases=requested===undefined?null:requested===0?[0]:requested===1?[1]:requested===2?[2,3]:[4];
 const collected:SangerCandidate[]=[];
 for(let phase=0;phase<STAGES.length;phase++){
  if(requestedPhases&&phase!==0&&!requestedPhases.includes(phase))continue;
  progress(phase/STAGES.length,LABELS[phase]);let primaries:SangerCandidate[];
  if(phase===0){
   const seen=new Set<string>();for(let k=s.primerMin;k<=s.primerMax;k++)for(const direction of ['forward','reverse'] as const){if(s.direction!=='both'&&s.direction!==direction)continue;for(const h of primarySites(direction,k,0)){const seq=queryAt(h,0),key=direction+seq;if(seen.has(key))continue;seen.add(key);const c=build(seq,direction,0,'one-exact');if(c)exactPrimaries.push(c);}}
   primaries=exactPrimaries.sort(compareSharedCandidates);
  }else if(phase===1)primaries=exactPrimaries;
  else if(phase<4){if(singlePrimaries===null)singlePrimaries=mismatchPrimaries(1);primaries=singlePrimaries.filter(c=>c.sharedDesign!.stage===STAGES[phase]);}
  else primaries=mismatchPrimaries(2);
  if(requestedPhases&&!requestedPhases.includes(phase))continue;
  const result:SangerCandidate[]=[];
  if(primaries.length>MAX_PRIMARIES)complete=false;
  for(let i=0;i<Math.min(primaries.length,MAX_PRIMARIES);i++){
   if(i%10===0)progress((phase+.9*i/Math.max(1,Math.min(primaries.length,MAX_PRIMARIES)))/STAGES.length,LABELS[phase]+' · '+(i+1)+' / '+Math.min(primaries.length,MAX_PRIMARIES));
   stats.partnerChecked++;const c=pairPrimary(primaries[i],phase===0,phase===1);if(c){c.sharedDesign={...c.sharedDesign!,stage:STAGES[phase]};result.push(c);}
   if(result.length>=100){if(i+1<primaries.length)complete=false;break;}
  }
  if(requested!==undefined)collected.push(...result);
  else if(result.length){progress(1,'Complete');return {candidates:result.sort(compareSharedCandidates),stats,excludedCandidates,complete,references:refs,settings:{...s}};}
 }
 progress(1,'Complete');return {candidates:collected.sort(compareSharedCandidates),stats,excludedCandidates,complete,references:refs,settings:{...s},...(requested!==undefined?{searchedStage:requested}:{})};
}
