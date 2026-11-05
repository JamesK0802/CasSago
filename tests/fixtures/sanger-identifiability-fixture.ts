import {profileFixture} from './sanger-profile-fixture';
import {synthReferences,SynthAllele,SynthEvent} from '../../src/app/utils/sanger-synthetic-benchmark';
import {initSynthProfiles} from '../../src/app/utils/sanger-synth-profiles';
/** Two distinct five-variant-per-target recipes with the same fluorescence.
 * All 30 variant percentages are distinct within either complete recipe. */
export function fiveGroupIdentifiabilityFixture(){
 const f=profileFixture(3),refs=synthReferences(f.result,f.candidate);
 refs.forEach((r,h)=>{
  initSynthProfiles(r);const cut=r.cuts[0],before=r.sequence.slice(cut,cut+2),after=before.replace(/[ACGT]/g,b=>b==='A'?'C':'A'),weights=[[13.01,12.03,7.07,11.09,17.13],[13.19,12.29,7.41,11.63,17.89],[13.37,12.59,7.73,11.97,18.17]][h],sub=(offset:number,length:number,percent:number,id:number):SynthAllele=>({id,percent,events:[{kind:'substitution',target:0,offset,length,selectionStart:cut+offset,before:before.slice(offset,offset+length),sequence:after.slice(offset,offset+length)}]});
  r.profiles![0].alleles=[{id:1,percent:weights[0],events:[{kind:'insertion',target:0,offset:3,length:3,sequence:'TGA'}]},{id:2,percent:weights[1],events:[{kind:'deletion',target:0,offset:-3,length:2,sequence:''}]},sub(0,1,weights[2],3),sub(1,1,weights[3],4),sub(0,2,weights[4],5)];
  const percentages=[4.01+.11*h,5.03+.13*h,6.07+.17*h,8.09+.19*h,9.11+.23*h].map(v=>Math.round(v*100)/100),events:SynthEvent[]=[{kind:'insertion',target:1,offset:-3,length:2,sequence:'AC'},{kind:'deletion',target:1,offset:-2,length:3,sequence:''},{kind:'insertion',target:1,offset:1,length:3,sequence:'GTA'},{kind:'deletion',target:1,offset:3,length:2,sequence:''},{kind:'substitution',target:1,offset:0,length:2,selectionStart:r.cuts[1],before:r.sequence.slice(r.cuts[1],r.cuts[1]+2),sequence:r.sequence.slice(r.cuts[1],r.cuts[1]+2).replace(/[ACGT]/g,b=>b==='A'?'C':'A')}];
  r.profiles![1].alleles=events.map((event,i)=>({id:101+i,percent:percentages[i],events:[event]}));
 });
 const alternate=structuredClone(refs);alternate.forEach((r,h)=>{const delta=[3.17,2.23,1.57][h];for(const i of [2,3])r.profiles![0].alleles[i].percent=Math.round((r.profiles![0].alleles[i].percent-delta)*100)/100;r.profiles![0].alleles[4].percent=Math.round((r.profiles![0].alleles[4].percent+delta)*100)/100;});
 return {...f,refs,alternate};
}
/** Stronger witness: disjoint insertion-sequence sets, no substitutions in
 * the changed profile. Five indel groups, distinct percentages, same trace. */
export function fiveIndelSequenceIdentifiabilityFixture(){
 const f=fiveGroupIdentifiabilityFixture(),refs=structuredClone(f.refs),alternate=structuredClone(f.refs);
 for(const [recipe,strings,weights]of [[refs,['AAA','ACC','CAC','CCA'],[11.5,12.5,9.5,6.5]],[alternate,['AAC','ACA','CAA','CCC'],[13.5,10.5,7.5,8.5]]]as const){
  recipe[0].profiles![0].alleles=[{id:1,percent:9,events:[{kind:'deletion',target:0,offset:-3,length:2,sequence:''}]},...strings.map((sequence,i)=>({id:i+2,percent:weights[i],events:[{kind:'insertion' as const,target:0,offset:1,length:3,sequence}]}))];
 }
 return {...f,refs,alternate};
}
