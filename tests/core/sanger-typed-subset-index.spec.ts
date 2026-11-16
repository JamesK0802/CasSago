import{it,expect}from'vitest';import{frequencySubsetIndex}from'../probes/experiments/sanger-typed-subset-index';
it('preserves exhaustive nearest-mask decisions including duplicate and equidistant sums',()=>{
 let state=190738;const random=()=>((state=(Math.imul(state,1664525)+1013904223)>>>0)/4294967296);
 for(let n=1;n<=13;n++){const sums=[0];for(let j=0;j<n;j++){const atom=Math.floor(random()*9)/100;for(let i=0,len=sums.length;i<len;i++)sums.push(sums[i]+atom);}const lookup=frequencySubsetIndex(sums);for(const value of[0,...sums.slice(0,20),...Array.from({length:100},()=>random())]){let mask=0;for(let i=1;i<sums.length;i++)if(Math.abs(sums[i]-value)<Math.abs(sums[mask]-value))mask=i;expect(lookup.mask(value)).toBe(mask);expect(lookup.distance(value)).toBe(Math.abs(sums[mask]-value));}}
 const lookup=frequencySubsetIndex([0,.2,.2,.4]);expect(lookup.mask(.3)).toBe(1);expect(lookup.mask(.2)).toBe(1);
});
