import type {LinearTerm} from './sanger-marginal-lp';
/** Only an objective certificate. The caller must ALSO freshly check every
 * original row, cone and variable box before reusing the proposed point. */
export function attainedBoxObjective(upper:number[],point:number[],terms:LinearTerm[],maximize=false){
 if(point.length!==upper.length)return;
 const weights=new Map<number,number>();for(const[i,w]of terms){if(!Number.isInteger(i)||i<0||i>=upper.length||!Number.isFinite(w))return;weights.set(i,(weights.get(i)||0)+w);}
 let value=0,valueCorrection=0,lower=0,lowerCorrection=0,magnitude=0;const sign=maximize?-1:1;
 for(const[i,w]of weights){if(!w)continue;const u=upper[i],x=point[i],c=sign*w;if(!Number.isFinite(c)||!Number.isFinite(u)||u<0||!Number.isFinite(x))return;
  const v=c*x,term=Math.min(0,c)*u;if(!Number.isFinite(v)||!Number.isFinite(term))return;
  const adjustedValue=v-valueCorrection,totalValue=value+adjustedValue;valueCorrection=(totalValue-value)-adjustedValue;value=totalValue;
  const adjustedLower=term-lowerCorrection,totalLower=lower+adjustedLower;lowerCorrection=(totalLower-lower)-adjustedLower;lower=totalLower;magnitude+=Math.abs(v)+Math.abs(term);
 }
 if(!Number.isFinite(magnitude))return;lower-=128*Number.EPSILON*(1+magnitude);const gap=value-lower;
 if(!Number.isFinite(gap)||gap< -1e-7||gap>1e-7)return;
 return{value:maximize?-value:value,lower};
}
