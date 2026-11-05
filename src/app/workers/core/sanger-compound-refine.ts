import type {PathEdit} from './sanger-path-common';
/** Propose bases on the measured READ axis of a compound molecule. All edit
 * coordinates remain on the original reference; earlier indels move later
 * observed cycles. Residual evidence proposes candidates, never truth labels. */
export function compoundBaseProposals(events:PathEdit[],column:number[],residual:number[],weight:number,positionIndex:Int32Array):PathEdit[][]{
 if(events.length<2||weight<=0)return [];const proposals:PathEdit[][]=[];
 for(let e=0;e<events.length;e++){
  const event=events[e];if(!event.inserted.length)continue;
  const shift=events.filter((_,j)=>j!==e).filter(p=>p.end<=event.start).reduce((s,p)=>s+p.inserted.length-(p.end-p.start),0),scores:number[][]=[];let inserted='';
  for(let j=0;j<event.inserted.length;j++){const k=positionIndex[event.start+shift+j]??-1;if(k<0){inserted='';break;}const values=Array.from({length:4},(_,b)=>residual[k*4+b]+weight*column[k*4+b]);scores.push(values);inserted+='ACGT'[values.indexOf(Math.max(...values))];}
  if(inserted.length!==event.inserted.length)continue;
  const replace=(sequence:string)=>events.map((p,j)=>({...p,inserted:j===e?sequence:p.inserted}));if(inserted!==event.inserted)proposals.push(replace(inserted));
  for(let j=0;j<scores.length;j++)for(let b=0;b<4;b++)if('ACGT'[b]!==event.inserted[j]&&scores[j][b]>.1*weight)proposals.push(replace(event.inserted.slice(0,j)+'ACGT'[b]+event.inserted.slice(j+1)));
 }
 return proposals;
}
