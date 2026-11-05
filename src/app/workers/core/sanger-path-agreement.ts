import {PooledAmplicon,SangerAnalysisResult} from '../../models/sanger-analysis';
/** A fallback is not an independent replicate. If a similarly good first fit
 * explains the same reference shares and net-size groups but implies different
 * target marginals, retain group evidence and withhold the target claim.
 * Agreement thresholds are prototype resolution policies, not p-values. */
export function guardCompetingTargetPaths(primary:SangerAnalysisResult,alternative:SangerAnalysisResult,products:PooledAmplicon[]):SangerAnalysisResult{
 if(!primary.joint||!alternative.joint||!primary.mixture||!alternative.mixture||!primary.mixture.converged)return alternative;
 if(primary.mixture.r2<.995)return alternative;
 const affected=new Set<string>();
 for(const product of products){
  const a=primary.joint.groups.filter(g=>g.reference===product.name),b=alternative.joint.groups.filter(g=>g.reference===product.name),sa=a.reduce((s,g)=>s+g.fraction,0),sb=b.reduce((s,g)=>s+g.fraction,0);
  if(Math.max(sa,sb)<.02)continue;
  if(Math.min(sa,sb)<.02||Math.abs(sa-sb)>.03||a.some(g=>g.netSize===null)||b.some(g=>g.netSize===null)){affected.add(product.name);continue;}
  const sizes=new Set(a.concat(b).map(g=>g.netSize));let tv=0;
  for(const size of sizes)tv+=Math.abs(a.filter(g=>g.netSize===size).reduce((s,g)=>s+g.fraction,0)/sa-b.filter(g=>g.netSize===size).reduce((s,g)=>s+g.fraction,0)/sb)/2;
  if(tv>.05){affected.add(product.name);continue;}
  for(const target of alternative.joint.targets.filter(t=>t.reference===product.name&&t.edited!==null)){
   let lower=0,upper=0;
   for(const group of a){const explanations=group.alternatives.length?group.alternatives:[{reference:product.name,events:group.events}],present=explanations.map(e=>e.events.some(event=>event.target===0||event.target===target.target));if(present.every(Boolean))lower+=group.fraction/sa;if(present.some(Boolean))upper+=group.fraction/sa;}
   if(Math.max(Math.abs(target.edited!-lower),Math.abs(target.edited!-upper))>.05)affected.add(product.name);
  }
 }
 if(affected.size){
  alternative.mixture.passesChecks=false;
  for(const target of alternative.joint.targets)if(affected.has(target.reference))target.edited=null;
  for(const ref of alternative.pooled?.amplicons??[])if(affected.has(ref.name))for(const key of ['wt','indel','inFrame','outOfFrame'] as const)ref[key]=null;
  for(const annotation of alternative.annotations??[])if(affected.has(annotation.name)){annotation.estimable=false;for(const component of annotation.components)component.unresolved=true;}
  alternative.warnings.push('Competing high-fit paths disagree on reference, net-size or target attribution for '+[...affected].join(', ')+'. Net-size groups are retained; target percentages are withheld.');
 }
 return alternative;
}
