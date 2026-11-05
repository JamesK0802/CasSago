import {PrimerSite,SangerDesignSettings} from '../../models/sanger-design';import {gramCondition} from './sanger-signal';
/** Prospective WT screen using only submitted references and primer/product
 * geometry. No edited trace, outcome, truth label or future allele is used.
 * This is scoped to normalized mixed-Sanger discrimination, not PCR success. */
export function screenSharedSangerWT(sites:PrimerSite[],settings:Pick<SangerDesignSettings,'readStart'|'readEnd'>){
 if(sites.length<2)return {passed:true,reason:'Single product',minimumDifferences:0,comparedBases:0,condition:1,warnings:[]as string[]};
 const end=Math.min(settings.readEnd,...sites.map(s=>s.extension.length)),positions=Array.from({length:Math.max(0,end-settings.readStart)},(_,i)=>i+settings.readStart).filter(p=>sites.every(s=>/^[ACGT]$/.test(s.extension[p]||'')));
 if(!positions.length)return {passed:false,reason:'No shared known WT sequence in the usable Sanger range',minimumDifferences:0,comparedBases:0,condition:null,warnings:[]as string[]};
 let differences=Infinity;for(let a=0;a<sites.length;a++)for(let b=a+1;b<sites.length;b++)differences=Math.min(differences,positions.filter(p=>sites[a].extension[p]!==sites[b].extension[p]).length);
 if(differences===0)return {passed:false,reason:'WT products have identical known sequence throughout the shared usable read range',minimumDifferences:0,comparedBases:positions.length,condition:null,warnings:[]as string[]};
 const columns=sites.map(s=>positions.flatMap(p=>'ACGT'.split('').map(b=>Number(s.extension[p]===b)))),geometry=gramCondition(columns);
 if(geometry.rank<sites.length)return {passed:false,reason:'WT mixture contributions are linearly indistinguishable in the shared usable range',minimumDifferences:differences,comparedBases:positions.length,condition:geometry.condition,warnings:[]as string[]};
 const warnings:string[]=[];if(differences<3)warnings.push('Only one or two informative WT positions');if(geometry.condition===null||geometry.condition>100)warnings.push('Poorly conditioned WT mixture');
 return {passed:true,reason:'WT geometry screen passed; future edited mixtures are not guaranteed identifiable',minimumDifferences:differences,comparedBases:positions.length,condition:geometry.condition,warnings};
}
