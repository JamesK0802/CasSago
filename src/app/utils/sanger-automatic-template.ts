import {resolveSangerReference} from './sanger-design.utils';
/** Read design metadata without importing obsolete analysis-engine or read-range
 * choices. Measured coverage is established later from the uploaded traces. */
export function parseAutomaticSangerTemplate(value:unknown){
 const v=value as any;if(v?.format!=='casmango-sanger-analysis-v1'||!Array.isArray(v.amplicons)||!v.amplicons.length||v.amplicons.length>4)throw Error('Use a Sanger primer template with 1–4 references.');
 const references=v.amplicons.map((a:any)=>{if(!Array.isArray(a.cuts)||a.cuts.length<1||a.cuts.length>2)throw Error('Use one or two targets per reference.');const r=resolveSangerReference(a.name,a.sequence,a.cuts.join(','),'cuts');if(Array.isArray(a.guides)&&a.guides.length===r.cuts.length&&a.guides.every((g:unknown)=>typeof g==='string'&&/^[ACGT]{17,30}$/i.test(g))){const g=resolveSangerReference(a.name,a.sequence,a.guides.join(','),'guides');if(g.cuts.some((cut,i)=>cut!==r.cuts[i]))throw Error('Template gRNA and cut positions disagree.');return {name:r.name,sequence:r.sequence,targets:a.guides.join(','),mode:'guides' as 'guides'|'cuts'};}return {name:r.name,sequence:r.sequence,targets:r.cuts.join(','),mode:'cuts' as 'guides'|'cuts'};});
 if(new Set(v.amplicons.map((a:any)=>a.cuts.length)).size!==1)throw Error('Provide corresponding targets in each reference.');
 if(new Set(references.map((r:any)=>r.name)).size!==references.length)throw Error('Reference names must be unique.');
 const primer=typeof v.primer==='string'?v.primer.replace(/\s/g,'').toUpperCase():'';
 if((references.length>1||references.some((r:any)=>r.targets.includes(',')))&&!/^[ACGT]{12,40}$/.test(primer))throw Error('A common sequencing primer is required for multiple references.');
 const s=v.settings||{},automatic=s.analysisPolicy==='automatic',maxDeletion=automatic?(s.maxDeletion??s.maxIndel??10):s.largeDeletions?(s.maxDeletion??200):(s.maxIndel??10),maxInsertion=automatic?(s.maxInsertion??s.maxIndel??10):s.largeDeletions?(s.maxInsertion??s.maxIndel??10):(s.maxIndel??10);
 if(!Number.isInteger(maxDeletion)||maxDeletion<1||maxDeletion>500||!Number.isInteger(maxInsertion)||maxInsertion<1||maxInsertion>30)throw Error('Use deletion limit 1–500 and insertion limit 1–30 bp.');
 const oppositeCandidates=v.amplicons.map((a:any)=>[a.forwardPrimer,a.reversePrimer].find((p:unknown)=>typeof p==='string'&&p.replace(/\s/g,'').toUpperCase()!==primer));
 const oppositePrimer=typeof v.oppositePrimer==='string'?v.oppositePrimer.replace(/\s/g,'').toUpperCase():oppositeCandidates.every((p:unknown)=>typeof p==='string')&&new Set(oppositeCandidates).size===1?String(oppositeCandidates[0]).replace(/\s/g,'').toUpperCase():'';
 if(oppositePrimer&&!/^[ACGT]{12,40}$/.test(oppositePrimer))throw Error('The opposite sequencing primer must contain 12–40 A/C/G/T bases.');
 return {references,primer,oppositePrimer,maxDeletion,maxInsertion,window:Number.isInteger(s.window)&&s.window>=60&&s.window<=500?s.window:200};
}
