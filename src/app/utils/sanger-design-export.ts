import {validateSangerProductBounds} from './sanger-product-bounds';
import {SangerCandidate,SangerDesignResult} from '../models/sanger-design';
import {SANGER_ANALYSIS_DEFAULTS,SangerAnalysisSettings} from '../models/sanger-analysis';
import {cleanReference} from './sanger-design.utils';
import {exactBindingSites} from '../workers/core/sanger-primer-designer';
export interface DesignedAmplicon {forwardPrimer?:string;reversePrimer?:string;name:string;sequence:string;cuts:number[];targetNames?:string[];guides?:string[];sourceReference:string;start:number;end:number}
export interface SangerAnalysisTemplate {format:'casmango-sanger-analysis-v1';primer:string;oppositePrimer?:string;amplicons:DesignedAmplicon[];settings:SangerAnalysisSettings;scopeNote:string}
export function selectedAmplicons(result:SangerDesignResult,c:SangerCandidate):DesignedAmplicon[]{
 if(!c.partner)throw new Error('Select a candidate with an opposite PCR primer to export complete amplicons.');
 return c.sites.map((hit,i)=>{const ref=result.references[c.referenceIndex??i],other=c.partner!.sites[i],start=Math.min(hit.start,other.start),end=Math.max(hit.end,other.end);return {name:ref.name.replace(/[^a-zA-Z0-9_-]/g,'_')+'_'+(start+1)+'-'+end,sequence:ref.sequence.slice(start,end),forwardPrimer:c.direction==='forward'?c.sequence:c.partner!.sequence,reversePrimer:c.direction==='reverse'?c.sequence:c.partner!.sequence,cuts:ref.cuts.map(cut=>cut-start),sourceReference:ref.name,start,end,...(ref.targetNames?{targetNames:[...ref.targetNames]}:{}),...(ref.guides?{guides:[...ref.guides]}:{})};});
}
export function makeAnalysisTemplate(result:SangerDesignResult,c:SangerCandidate):SangerAnalysisTemplate{
 if(result.settings.platform==='ngs')throw new Error('This is an NGS design. Download Amplicon FASTA and Targets CSV for the main analysis.');
 const amplicons=selectedAmplicons(result,c),settings={...SANGER_ANALYSIS_DEFAULTS,analysisPolicy:'automatic' as const,maxDeletion:result.settings.indelRange,maxInsertion:Math.min(30,result.settings.indelRange),scope:amplicons.length===1?'single' as const:'pooled' as const,window:result.settings.window,maxIndel:Math.min(30,result.settings.indelRange),readStart:result.settings.readStart,readEnd:result.settings.readEnd,direction:amplicons.length===1?c.direction:'auto' as const};
 return validateAnalysisTemplate({format:'casmango-sanger-analysis-v1',primer:c.sequence,oppositePrimer:c.partner!.sequence,amplicons,settings,scopeNote:'Predicted PCR products. Confirm actual submitted amplicons after sequencing; no trace data or experimental validation included.'});
}
export function validateAnalysisTemplate(value:unknown):SangerAnalysisTemplate{
 const v=value as SangerAnalysisTemplate;
 if(v?.format!=='casmango-sanger-analysis-v1'||!Array.isArray(v.amplicons)||!v.amplicons.length||v.amplicons.length>4)throw new Error('Sanger templates support 1–4 submitted amplicons.');
 const a=v.amplicons.map(r=>{if(typeof r.name!=='string'||!r.name.trim()||typeof r.sequence!=='string'||!Array.isArray(r.cuts)||(r.cuts.length<1||r.cuts.length>2))throw new Error('Sanger analysis supports one or two cut sites per amplicon.');const sequence=cleanReference(r.sequence);if(sequence.length>20000||new Set(r.cuts).size!==r.cuts.length||r.cuts.some(c=>!Number.isInteger(c)||c<=0||c>=sequence.length))throw new Error('Invalid submitted amplicon or cut position.');const product={...r,sequence};validateSangerProductBounds({...product,cut:r.cuts[0]});return product;});
 if(new Set(a.map(r=>r.name)).size!==a.length)throw new Error('Template amplicon names must be unique.');
 const primer=typeof v.primer==='string'?v.primer.replace(/\s/g,'').toUpperCase():'';if(!/^[ACGT]{16,30}$/.test(primer))throw new Error('Template requires the sequencing primer (16–30 nt).');
 const s={...SANGER_ANALYSIS_DEFAULTS,...v.settings,scope:a.length===1?'single' as const:'pooled' as const,method:'mixture' as const};
 const automatic=s.analysisPolicy==='automatic';if(automatic){s.maxDeletion=v.settings?.maxDeletion??s.maxIndel;s.maxInsertion=v.settings?.maxInsertion??Math.min(30,s.maxIndel);}
 const reserve=automatic?Math.max(s.maxDeletion??s.maxIndel,s.maxInsertion??s.maxIndel):s.maxIndel;
 if(!Number.isInteger(s.window)||s.window<60||s.window>500||!Number.isInteger(s.maxIndel)||s.maxIndel<1||s.maxIndel>(automatic?30:a.length>1?15:30))throw new Error('Sanger analysis requires a 60–500 bp window and indel range ±1–15 bp (single: ±1–30). Adjust the design before making a template.');
 if(!Number.isInteger(s.readStart)||!Number.isInteger(s.readEnd)||s.readStart<0||s.readEnd>1000||s.readEnd<=s.readStart)throw new Error('Analysis templates require a usable Sanger range within 0–1,000 bp.');
 if(s.joinNearby!==undefined&&(!Number.isInteger(s.joinNearby)||s.joinNearby<0||s.joinNearby>10))throw new Error('Join nearby changes must be 0–10 matched bases.');
 if(!['auto','forward','reverse'].includes(s.direction)||!Number.isInteger(s.controlPhred)||s.controlPhred<0||s.controlPhred>60||!['classic','calibrated'].includes(s.signalModel||''))throw new Error('Invalid analysis settings.');
 a.forEach(r=>{const hits=exactBindingSites(primer,{...r,cuts:r.cuts});if(hits.length!==1)throw new Error(r.name+': sequencing primer must bind once in the submitted product.');const h=hits[0];for(const inputCut of r.cuts){const cut=h.direction==='forward'?inputCut-h.end:h.start-inputCut;const extensionLength=h.direction==='forward'?r.sequence.length-h.end:h.start;if(cut+Math.ceil(s.window/2)+reserve>extensionLength)throw new Error(r.name+': PCR product ends before the complete analysis region.');if((automatic||a.length>1||r.cuts.length>1||s.largeDeletions)&&cut<s.readStart+72)throw new Error(r.name+': mixed Sanger needs at least 72 upstream bases after the usable read start.');if(cut-Math.floor(s.window/2)-reserve<s.readStart||cut+Math.ceil(s.window/2)+reserve>s.readEnd)throw new Error(r.name+': target window exceeds the usable Sanger range.');}});
 if(s.largeDeletions!==undefined&&typeof s.largeDeletions!=='boolean')throw new Error('Invalid large-deletion setting.');
 if(s.maxDeletion!==undefined&&(!Number.isInteger(s.maxDeletion)||s.maxDeletion<(automatic||s.jointEngine&&s.jointEngine!=='staged'?1:s.maxIndel)||s.maxDeletion>500))throw new Error('Deletion search limit must be 1–500 bp for signal paths, or at least the small-indel range for staged models.');
 if(s.jointEngine!==undefined&&!['staged','signal-local','signal-path','signal-hybrid'].includes(s.jointEngine))throw new Error('Invalid joint inference method.');
 if(s.editStartRadius!==undefined&&(!Number.isInteger(s.editStartRadius)||s.editStartRadius<0||s.editStartRadius>100))throw new Error('Indel origin radius must be 0–100 bp.');
 if(s.maxInsertion!==undefined&&(!Number.isInteger(s.maxInsertion)||s.maxInsertion<1||s.maxInsertion>(automatic?30:60)))throw new Error('Insertion search limit must be 1–60 bp.');
 return {...v,primer,amplicons:a,settings:s};
}
