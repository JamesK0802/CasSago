import {PooledAmplicon,SangerAnalysisResult,SangerJointEvent} from '../../models/sanger-analysis';
/** A repeat can place the same final molecule on either side of a cut.
 * Enumerate sequence-equivalent junctions, not just paths found by the search.
 * Size/fraction evidence is retained; non-identifiable target attribution is withheld. */
export function guardEquivalentJunctions(result:SangerAnalysisResult,products:PooledAmplicon[],originRadius?:number):SangerAnalysisResult{
 if(!result.joint)return result;
 const affected=new Set<string>();
 for(const group of result.joint.groups){
  const explanations=group.alternatives.length?group.alternatives:group.reference?[{reference:group.reference,events:group.events}]:[];
  const expanded:{reference:string;events:SangerJointEvent[]}[]=[];
  for(const explanation of explanations){
   const product=products.find(p=>p.name===explanation.reference);if(!product)continue;
   const reference=product.sequence,cuts=[...(product.cuts??[product.cut])].sort((a,b)=>a-b);
   const target=(start:number,end:number)=>{if(cuts.length===2&&start<=cuts[0]&&end>=cuts[1])return 0;const distance=cuts.map(c=>Math.min(Math.abs(start-c),Math.abs(end-c)));const nearest=distance.indexOf(Math.min(...distance));return Math.min(...distance)>60?-1:nearest+1;};
   const semantic=(events:SangerJointEvent[])=>[...new Set(events.flatMap(e=>e.target===0?cuts.map((_,i)=>i+1):[e.target]))].sort().join('/');
   const original=semantic(explanation.events);expanded.push(explanation);
   // A junction inside both cut-origin zones cannot support a single-target
   // claim just because its fitted end falls one base short of the exact cut.
   if(originRadius!==undefined)for(const event of explanation.events){
    const compatible=cuts.flatMap((cut,i)=>event.start<=cut+originRadius&&event.end>=cut-originRadius?[i+1]:[]);
    if(event.end>event.start&&event.target>0&&compatible.length>1)affected.add(explanation.reference);
   }
   for(let i=0;i<explanation.events.length;i++){
    const event=explanation.events[i],insert=event.insertedSequence??'N'.repeat(event.inserted),deleted=event.end-event.start;
    if(deleted===insert.length||insert.includes('N'))continue;
    const molecule=reference.slice(0,event.start)+insert+reference.slice(event.end);let prefix=0,suffix=0;
    while(prefix<Math.min(reference.length,molecule.length)&&reference[prefix]===molecule[prefix])prefix++;
    while(suffix<Math.min(reference.length,molecule.length)&&reference[reference.length-1-suffix]===molecule[molecule.length-1-suffix])suffix++;
    const lo=Math.max(0,reference.length-suffix-deleted,...explanation.events.filter((_,k)=>k!==i).filter(e=>e.end<=event.start).map(e=>e.end));
    const hi=Math.min(prefix,reference.length-deleted,molecule.length-insert.length,...explanation.events.filter((_,k)=>k!==i).filter(e=>e.start>=event.end).map(e=>e.start-deleted));
    const represented=new Set<string>([original]);
    for(let start=lo;start<=hi;start++){
     if(originRadius!==undefined&&!cuts.some(c=>start<=c+originRadius&&start+deleted>=c-originRadius))continue;
     const shifted={...event,start,end:start+deleted,insertedSequence:molecule.slice(start,start+insert.length),target:target(start,start+deleted)},events=explanation.events.map((e,k)=>k===i?shifted:e),key=semantic(events);
     if(represented.has(key))continue;represented.add(key);expanded.push({reference:explanation.reference,events});affected.add(explanation.reference);
    }
   }
  }
  if(expanded.length>explanations.length){group.alternatives=expanded;group.events=[];}
 }
 if(affected.size){
  result.mixture!.passesChecks=false;
  for(const t of result.joint.targets)if(affected.has(t.reference))t.edited=null;
  for(const annotation of result.annotations??[])if(affected.has(annotation.name))for(const component of annotation.components)component.unresolved=true;
  result.warnings.push('Repeat-equivalent junctions or overlapping cut-origin zones leave target attribution unresolved for '+[...affected].join(', ')+'. Target editing estimates are withheld; fitted net sizes and fractions are retained.');
 }
 return result;
}
