import {mixPhase,phaseNeighbor} from './sanger-signal';
/** Empirical fluorescence for a known WT (previous,current) base pair. This
 * models dye spill and sequence-dependent phase using CONTROL observations,
 * not a truth allele or a sample edit. Every third recurrence is held out. */
export interface TransitionResponse {markerValidation?:{oldMse:number;newMse:number;trainingPositions:number;validationPositions:number;accepted:boolean};used:boolean;table:number[][];pairYields:number[];trainingCounts:number[];validationCounts:number[];rawValidationMse:number;validationMse:number}
const bases='ACGT',median=(v:number[])=>{const s=[...v].sort((a,b)=>a-b);return s.length?s[Math.floor(s.length/2)]:0;};
export function fitTransitionResponse(raw:number[][],sequences:string[],positions:number[],lag:number,baseStrengths:number[]=[1,1,1,1]):TransitionResponse {
 const groups=Array.from({length:16},()=>[] as number[]);
 for(const p of positions){const pair=sequences[0].slice(p-1,p+1),l=bases.indexOf(pair[0]),b=bases.indexOf(pair[1]);if(pair.length!==2||l<0||b<0||!raw[p]?.some(v=>v>0)||!sequences.every(s=>s.slice(p-1,p+1)===pair))continue;groups[4*l+b].push(p);}
 const train=groups.map(g=>g.filter((_,i)=>i%3!==0)),valid=groups.map(g=>g.filter((_,i)=>i%3===0)),basic=groups.map((_,k)=>Array.from(bases,(b,i)=>(i===k%4?1-Math.abs(lag):0)+(lag<0?Math.abs(lag)/4:i===Math.floor(k/4)?lag:0)));
 const table=train.map((g,k)=>{if(g.length<4)return basic[k];const vector=bases.split('').map((_,b)=>median(g.map(p=>raw[p][b]))),total=vector.reduce((s,v)=>s+v,0);return vector.map(v=>total?v/total:0);});
 const score=(vectors:number[][])=>valid.reduce((sum,g,k)=>sum+g.reduce((s,p)=>s+raw[p].reduce((e,v,b)=>e+(v-vectors[k][b])**2,0),0),0)/Math.max(1,4*valid.flat().length);
 const before=score(basic),after=score(table),used=valid.flat().length>=30&&after<before*.95&&after+1e-7<before;
 const chosen=used?table:basic,pairYields=chosen.map((v,k)=>{
  const l=Math.floor(k/4),b=k%4;if(l===b)return baseStrengths[b];
  const current=chosen[4*b+b],previous=chosen[4*l+l],difference=previous.map((x,d)=>x-current[d]),norm=difference.reduce((s,x)=>s+x*x,0);
  const fraction=Math.max(0,Math.min(1,norm?difference.reduce((s,x,d)=>s+x*(v[d]-current[d]),0)/norm:0));
  // A normalized colour mixture loses its absolute photon yield. Recover that
  // yield from the two pure-base responses before mixing molecular fractions.
  return 1/((1-fraction)/baseStrengths[b]+fraction/baseStrengths[l]);
 });
 return {used,table:chosen,pairYields,trainingCounts:train.map(g=>g.length),validationCounts:valid.map(g=>g.length),rawValidationMse:before,validationMse:used?after:before};
}
export function transitionProbabilityEmission(left:number[],current:number[],right:number[],theta:number[],table:number[][],baseYields:number[]):number[]{
 const previous=left.some(v=>v>0)?left:[.25,.25,.25,.25],out=[0,0,0,0];
 for(let l=0;l<4;l++)if(previous[l])for(let b=0;b<4;b++)if(current[b]){const context=right.reduce((sum,v,r)=>v?sum+v*Math.exp((theta[4+r]||0)+(theta[8+16*l+4*b+r]||0)):sum,0)||1;const amount=previous[l]*current[b]*context*Math.exp(theta[l]||0)*(baseYields[4*l+b]??1),response=table[4*l+b];for(let dye=0;dye<4;dye++)out[dye]+=amount*response[dye];}
 return out;
}
export function transitionEmission(sequence:string,p:number,lag:number,theta:number[],table:number[][],baseYields:number[]):number[]{
 const prob=(at:number)=>Array.from(bases,b=>Number(sequence[at]===b)),one=(at:number)=>transitionProbabilityEmission(prob(at-1),prob(at),prob(at+1),theta,table,baseYields),a=one(p),b=one(phaseNeighbor(p,lag));return mixPhase(a,b,lag);
}
