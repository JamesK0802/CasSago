/** Exact nearest subset with native numeric sorting. All distinct IEEE values
 * and their first enumeration indices are retained, including ties. */
export function frequencySubsetIndex(sums:number[]){
 const first=new Map<number,number>();for(let i=0;i<sums.length;i++)if(!first.has(sums[i]))first.set(sums[i],i);
 const values=Float64Array.from(first.keys()).sort();
 const mask=(value:number)=>{if(!values.length)return-1;let lo=0,hi=values.length;while(lo<hi){const mid=(lo+hi)>>>1;if(values[mid]<value)lo=mid+1;else hi=mid;}if(!lo)return first.get(values[0])!;if(lo===values.length)return first.get(values[lo-1])!;const a=values[lo-1],b=values[lo],da=Math.abs(a-value),db=Math.abs(b-value),ia=first.get(a)!,ib=first.get(b)!;return da<db?ia:db<da?ib:Math.min(ia,ib);};
 return{mask,distance:(value:number)=>{const i=mask(value);return i<0?Infinity:Math.abs(sums[i]-value);}};
}
