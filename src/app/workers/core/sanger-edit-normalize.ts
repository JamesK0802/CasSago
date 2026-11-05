export interface NormalizableEdit {start:number;end:number;inserted:string}
/** Remove matched context and cancelling edit histories while preserving the
 * final molecule sequence exactly. Does not invent a unique physical repair
 * history or collapse alternative placements through microhomology. */
export function normalizeMoleculeEdits(reference:string,input:NormalizableEdit[]):NormalizableEdit[]{
 let events=input.map(e=>({...e})).sort((a,b)=>a.start-b.start);
 const trim=(e:NormalizableEdit)=>{let {start,end,inserted}=e;while(start<end&&inserted.length&&reference[start]===inserted[0]){start++;inserted=inserted.slice(1);}while(start<end&&inserted.length&&reference[end-1]===inserted.at(-1)){end--;inserted=inserted.slice(0,-1);}return {start,end,inserted};};
 events=events.map(trim).filter(e=>e.end!==e.start||e.inserted.length);
 for(let round=0;round<8;round++){let changed=false;
  outer:for(let i=0;i<events.length;i++)for(let j=events.length-1;j>i;j--){let at=events[i].start,sequence='',valid=true;for(let k=i;k<=j;k++){const e=events[k];if(e.start<at){valid=false;break;}sequence+=reference.slice(at,e.start)+e.inserted;at=e.end;}if(!valid)continue;const end=events[j].end;if(sequence===reference.slice(events[i].start,end)){events.splice(i,j-i+1);changed=true;break outer;}}
  if(!changed)break;
 }
 return events;
}
