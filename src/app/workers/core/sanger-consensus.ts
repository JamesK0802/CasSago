import {SangerVariant,TraceAlignment,TraceAlignmentColumn} from '../../models/sanger-analysis';
import {reverseComplement} from '../../utils/biology.utils';
/** Local affine-gap alignment for consensus variant representation only.
 * Registration/calibration retain their independent alignment policy. */
export function alignConsensusSequence(reference:string,query:string):{score:number;identity:number;columns:TraceAlignmentColumn[]}{
 const n=reference.length,m=query.length;if(n*m>24000000)throw new Error('Alignment region is too large.');
 const width=m+1,trace=[new Uint8Array((n+1)*width),new Uint8Array((n+1)*width),new Uint8Array((n+1)*width)];
 let previous=[new Int32Array(width),new Int32Array(width),new Int32Array(width)],current=[new Int32Array(width),new Int32Array(width),new Int32Array(width)],best=0,bi=0,bj=0,bs=0;
 for(let i=1;i<=n;i++){
  current.forEach(row=>row.fill(0));
  for(let j=1;j<=m;j++){
   const equal=reference[i-1]===query[j-1]&&/[ACGT]/.test(query[j-1]),sub=equal?3:-4,at=i*width+j;
   const choices=[[previous[0][j-1]+sub,previous[1][j-1]+sub,previous[2][j-1]+sub],[previous[0][j]-7,previous[1][j]-1,previous[2][j]-7],[current[0][j-1]-7,current[1][j-1]-7,current[2][j-1]-1]];
   for(let state=0;state<3;state++){let value=0,from=0;for(let k=0;k<3;k++)if(choices[state][k]>value){value=choices[state][k];from=k+1;}current[state][j]=value;trace[state][at]=from;if(value>best){best=value;bi=i;bj=j;bs=state;}}
  }
  [previous,current]=[current,previous];
 }
 const columns:TraceAlignmentColumn[]=[];let i=bi,j=bj,state=bs;
 while(i>0&&j>0){const from=trace[state][i*width+j];if(!from)break;if(state===0){columns.push({ref:i-1,query:j-1,a:reference[i-1],b:query[j-1]});i--;j--;}else if(state===1){columns.push({ref:i-1,query:null,a:reference[i-1],b:'-'});i--;}else{columns.push({ref:null,query:j-1,a:'-',b:query[j-1]});j--;}state=from-1;}
 columns.reverse();return {score:best,identity:columns.length?columns.filter(c=>c.a===c.b&&/[ACGT]/.test(c.a)).length/columns.length:0,columns};
}
/** Join nearby differences as one representation, not a claim of one repair event.
 * Keep primitive parts; only exact matched bridges inside cut±20 may be joined. */
export function consensusVariants(alignment:TraceAlignment,referenceLength:number,window:number,joinGap=3):SangerVariant[]{
 if(!Number.isInteger(joinGap)||joinGap<0||joinGap>10)throw new Error('Join nearby changes must be 0–10 matched bases.');
 type Block={columns:TraceAlignmentColumn[];start:number;end:number;parts:Block[]};
 const blocks:Block[]=[];let run:TraceAlignmentColumn[]=[],cursor=alignment.referenceStart;
 const flush=()=>{if(!run.length)return;const coords=run.flatMap(c=>c.ref===null?[]:[c.ref]);const start=coords.length?Math.min(...coords):cursor,end=coords.length?Math.max(...coords)+1:start;blocks.push({columns:run,start,end,parts:[]});run=[];};
 for(const c of alignment.columns){if(c.a===c.b&&/[ACGT]/.test(c.a)){flush();cursor=c.ref!+1;}else{if(!run.length&&c.ref!==null)cursor=c.ref;run.push(c);if(c.ref!==null)cursor=c.ref+1;}}flush();
 const groups:Block[]=[];
 for(const block of blocks){const last=groups.at(-1),gap=last?block.start-last.end:Infinity;
  if(last&&joinGap>0&&gap>=0&&gap<=joinGap&&last.start>=alignment.cut-20&&block.end<=alignment.cut+20&&last.columns.every(c=>!/[NRYWSKMBDHV]/.test(c.a+c.b))&&block.columns.every(c=>!/[NRYWSKMBDHV]/.test(c.a+c.b))){
   // Recover the exact matched bridge from the alignment, preserving all bases.
   const bridge=alignment.columns.filter(c=>c.ref!==null&&c.ref>=last.end&&c.ref<block.start&&c.a===c.b&&/[ACGT]/.test(c.a));
   if(bridge.length===gap){last.columns.push(...bridge,...block.columns);last.end=block.end;last.parts.push(block);continue;}
  }
  groups.push({...block,columns:[...block.columns],parts:[block]});
 }
 const represent=(block:Block):SangerVariant=>{let ref=block.columns.map(c=>c.a).join('').replace(/-/g,''),alt=block.columns.map(c=>c.b).join('').replace(/-/g,'');if(alignment.direction==='reverse'){ref=reverseComplement(ref);alt=reverseComplement(alt);}const netSize=alt.length-ref.length,hasInsertion=block.columns.some(c=>c.ref===null),hasDeletion=block.columns.some(c=>c.query===null),hasSubstitution=block.columns.some(c=>c.ref!==null&&c.query!==null&&c.a!==c.b),ambiguous=/[^ACGT]/.test(ref+alt);
  const type:SangerVariant['type']=ambiguous?'Ambiguous call':(hasInsertion&&hasDeletion)||((hasInsertion||hasDeletion)&&hasSubstitution)?'Complex replacement':netSize===0?'Substitution':netSize%3===0?'In-frame indel':'Out-of-frame indel';
  const position=alignment.direction==='forward'?block.start+1:referenceLength-block.end+1;
  return {position,endPosition:position+Math.max(1,ref.length)-1,reference:ref||'-',sample:alt||'-',netSize,type};
 };
 return groups.filter(b=>b.end>=alignment.cut-window/2&&b.start<=alignment.cut+window/2).map(g=>{const v=represent(g);if(g.parts.length>1){v.grouped=true;v.parts=g.parts.map(represent).sort((a,b)=>a.position-b.position);v.matchedBridgeBases=g.columns.filter(c=>c.a===c.b&&/[ACGT]/.test(c.a)).length;}return v;}).sort((a,b)=>a.position-b.position);
}
