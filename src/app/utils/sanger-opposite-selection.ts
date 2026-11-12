import {CommonPartner,OppositePrimerChoice,SangerCandidate,SangerDesignResult} from '../models/sanger-design';
import {coversRead,complementarity,site,traceDistance} from '../workers/core/sanger-primer-designer';
import {usablePrimerRead,screenPrimerBinding} from '../workers/core/sanger-primer-safety';
import {screenSharedSangerWT} from '../workers/core/sanger-primer-separation';

/** Change only the opposite PCR choice, preserving the sequencing-primer search. */
export function chooseOppositePrimers(result:SangerDesignResult,c:SangerCandidate,selection:CommonPartner|OppositePrimerChoice[]):SangerCandidate {
 const {settings:s}=result,refs=c.sites.map((_,i)=>result.references[c.referenceIndex??i]);
 const common=!Array.isArray(selection);
 const items=Array.isArray(selection)?selection:selection.sites.map((h,i)=>({primer:selection,site:h,product:selection.products[i],safety:screenPrimerBinding(selection.sequence,result.references,selection.sites)}));
 if(items.length!==c.sites.length)throw Error('Choose one opposite primer for each reference.');
 for(let i=0;i<items.length;i++){
  const q=items[i];if(!q||q.site.reference!==c.sites[i].reference||q.site.direction===c.direction)throw Error('Opposite primer does not match this reference and direction.');
  const actual=refs[i].sequence.slice(q.site.start,q.site.end),expected=q.site.direction==='forward'?actual:reverseComplement(actual);
  if(expected!==q.primer.sequence)throw Error('The opposite primer must match the reference exactly.');
  const physicalProduct=c.direction==='forward'?q.site.end-c.sites[i].start:c.sites[i].end-q.site.start;
  if(q.product!==physicalProduct)throw Error('Opposite primer product coordinates are inconsistent.');
  if(Math.abs(q.primer.tm-c.tm)>s.maxTmDifference||q.product<s.productMin||q.product>s.productMax)throw Error('Opposite primer does not meet the product and Tm limits.');
  if(c.sharedDesign&&(!q.safety.minimumExact||q.safety.minimumMismatches<5||q.safety.references.some(r=>r.hasUnknownBases)))throw Error('Opposite primer does not pass the five-base non-target difference screen.');
  for(const seq of [c.sequence,...items.slice(0,i).map(p=>p.primer.sequence)]){const pair=complementarity(seq,q.primer.sequence);if(pair.end>=5||pair.any>=9)throw Error('These PCR primers have excessive complementarity. Choose another opposite primer.');}
 }
 const sites=c.sites.map((h,i)=>{const at=c.direction==='forward'?h.start:refs[i].sequence.length-h.end,full=site(refs[i],at,c.sequence.length,c.direction,s);return {...full,mismatchPositions:h.mismatchPositions,extension:full.extension.slice(0,items[i].product-c.sequence.length)};});
 if(!coversRead(sites,refs,s))throw Error('Opposite primer ends before a required target window.');
 const range=usablePrimerRead(s,Math.min(...sites.map(h=>h.extension.length)));
 if(s.workflow==='shared-sanger'&&Math.max(...sites.flatMap(h=>h.cutOffsets))+Math.ceil(s.window/2)+s.indelRange>range.end)throw Error('Choose longer products to cover every target in the shared Sanger read.');
 const wt=s.workflow==='shared-sanger'?screenSharedSangerWT(sites,{readStart:range.start,readEnd:range.end}):null;
 if(wt&&!wt.passed)throw Error(wt.reason);
 const partner:CommonPartner={...items[0].primer,...(!common?{perReference:items.map(q=>q.primer)}:{}),sites:items.map(q=>q.site),products:items.map(q=>q.product),pairEnd:Math.max(...items.map(q=>complementarity(c.sequence,q.primer.sequence).end))};
 const out:SangerCandidate={...c,id:c.direction+'|'+c.sequence+'|'+items.map(q=>q.primer.sequence+':'+q.site.start).join('|'),partner,sites,partnerSafety:common?[items[0].safety]:items.map(q=>q.safety)};
 const screens=items.map(q=>q.safety),worst=screens.reduce((a,b)=>a.minimumMismatches<b.minimumMismatches?a:b);
 if(c.primerSafety)out.primerSafety={primary:c.primerSafety.primary,partner:worst};
 [out.minTraceDifferences,out.traceComparedBases]=traceDistance(sites,s);
 if(wt)out.wtSeparation={minimumDifferences:wt.minimumDifferences,comparedBases:wt.comparedBases,condition:wt.condition,warnings:wt.warnings};
 return out;
}
import {reverseComplement} from './biology.utils';
