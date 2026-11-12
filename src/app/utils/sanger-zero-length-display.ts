import type {SangerAnnotation,SangerMarginalGroup,SangerJointEvent} from '../models/sanger-analysis';

/** Display decomposition only. Preserve the joint zero-length mass and its
 * bounds; do not treat a zero net indel as WT or invent SNP linkage. */
export function zeroLengthDisplay(a:SangerAnnotation,g:SangerMarginalGroup){
 if(!g.bases.length||g.size!==0||!g.unmodifiedRange||!g.substitutionRange)return null;
 const first=a.direction==='reverse'?(a.cuts?.length||1):1,offset=g.baselineReadOffset??(g.target===first?0:null);if(offset===null)return null;
 const baseline=(p:number)=>'ACGT'.indexOf(a.reference[p+offset]);if(g.bases.some(p=>baseline(p.position)<0))return null;
 const f=g.fraction,flo=g.fractionRange[0],fhi=g.fractionRange[1],subRange:[number,number]=[Math.max(0,g.substitutionRange[0]),Math.min(fhi,g.substitutionRange[1])],wtRange:[number,number]=[Math.max(0,g.unmodifiedRange[0]),Math.min(fhi,g.unmodifiedRange[1])];
 const subPoint=Math.min(subRange[1],Math.max(subRange[0],...g.bases.map(p=>f*(1-p.estimate[baseline(p.position)])))),wtPoint=Math.max(0,f-subPoint),quantified=wtRange[1]-wtRange[0]<=.020001&&subRange[1]-subRange[0]<=.020001;
 if(subRange[0]<=1e-5||subPoint<=1e-7)return{wtPoint,subPoint,wtRange,subRange,quantified,split:false as const};
 const wtBases=g.bases.map(p=>{const b=baseline(p.position),one=[0,1,2,3].map(i=>Number(i===b));return{position:p.position,estimate:one,lower:one.slice(),upper:one.slice()};});
 const subBases=g.bases.map(p=>{const b=baseline(p.position),estimate=p.estimate.map((v,i)=>Math.max(0,Math.min(1,(f*v-(i===b?wtPoint:0))/subPoint))),lower=p.lower.map((v,i)=>Math.max(0,Math.min(1,i===b?1-fhi*(1-v)/Math.max(1e-30,subRange[0]):flo*v/Math.max(1e-30,subRange[1])))),upper=p.upper.map((v,i)=>Math.max(0,Math.min(1,i===b?1-flo*(1-v)/Math.max(1e-30,subRange[1]):fhi*v/Math.max(1e-30,subRange[0]))));return{position:p.position,estimate,lower,upper};});
 const events:SangerJointEvent[]=[],changed=g.bases.filter(p=>1-p.estimate[baseline(p.position)]>1e-5);
 for(const p of subBases){const b=p.estimate.indexOf(Math.max(...p.estimate)),at=p.position+offset;if(b!==baseline(p.position))events.push({target:g.target,start:at,end:at+1,inserted:1,insertedSequence:'ACGT'[b]});}
 events.sort((a,b)=>a.start-b.start);
 return{wtPoint,subPoint,wtRange,subRange,quantified,split:true as const,wtBases,subBases,events,mixture:changed.length!==1||subBases.some(p=>p.estimate.filter(v=>v>1e-5).length>1)};
}
