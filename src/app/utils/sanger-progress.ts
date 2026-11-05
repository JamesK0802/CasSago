/** Estimated stage budget, not a measurement of wall-clock completion. */
export function sangerProgress(stage:string):{percent:number;label:string}{
 if(stage==='Checking ground truth after inference')return{percent:99,label:'Ground-truth comparison'};
 if(stage.startsWith('Conditional prefix · ')){const inner=sangerProgress(stage.slice(21));return{percent:5+inner.percent*.2,label:'1 / 5 · First-target prefix'};}
 const refine=stage.match(/Refining first-target allocations · reference (\d+) \/ (\d+) · (\d+) \/ (\d+)/);
 if(refine){const fraction=(Number(refine[1])-1+Number(refine[3])/Number(refine[4]))/Number(refine[2]);return{percent:60+30*fraction,label:'3 / 4 · Joint base assignment · reference '+refine[1]+' / '+refine[2]+' · '+refine[3]+' / '+refine[4]};}
 const discrete=stage.match(/Joint base assignment · pass (\d+) · reference (\d+) \/ (\d+)/);
 if(discrete)return{percent:82+Math.min(12,2*(Number(discrete[1])-1+Number(discrete[2])/Number(discrete[3]))),label:'4 / 5 · Joint base assignment · reference '+discrete[2]+' / '+discrete[3]+' · pass '+discrete[1]};
 if(stage.includes('first-target measured prefix'))return{percent:5,label:'1 / 5 · First-target prefix'};
 if(stage.includes('downstream total-length'))return{percent:30,label:'2 / 5 · Downstream length fractions'};
 if(stage.includes('Proposing second-target')){const m=stage.match(/reference (\d+) \/ (\d+)/);return{percent:40+(m?(Number(m[1])-1)/Number(m[2])*25:0),label:'3 / 5 · Second-target candidates'+(m?' · reference '+m[1]+' / '+m[2]:'')};}
 const allocations=stage.match(/^(First target reference|Reference) (\d+) \/ (\d+) · Comparing length-class allocations · (\d+) \/ (\d+)/);
 if(allocations){const fraction=(Number(allocations[2])-1+Number(allocations[4])/Number(allocations[5]))/Number(allocations[3]);return{percent:allocations[1]==='Reference'?40+25*fraction:40+40*fraction,label:(allocations[1]==='Reference'?'3 / 5 · Second-target candidates':'3 / 4 · Sequence groups')+' · reference '+allocations[2]+' / '+allocations[3]+' · '+allocations[4]+' / '+allocations[5]};}
 const joint=stage.match(/Comparing joint target hypotheses · reference (\d+) \/ (\d+) · (\d+) \/ (\d+)/);
 if(joint){const fraction=(Number(joint[1])-1+Number(joint[3])/Number(joint[4]))/Number(joint[2]);return{percent:65+15*fraction,label:'4 / 5 · Joint signal fit · reference '+joint[1]+' / '+joint[2]+' · '+joint[3]+' / '+joint[4]};}
 const phase=stage.match(/Sequence phase · pass (\d+) · reference (\d+) \/ (\d+) · (\d+) \/ (\d+)/);
 if(phase){const fraction=(Number(phase[1])-1+(Number(phase[2])-1+Number(phase[4])/Number(phase[5]))/Number(phase[3]))/3;return{percent:82+12*fraction,label:'4 / 5 · Sequence groups · reference '+phase[2]+' / '+phase[3]+' · '+phase[4]+' / '+phase[5]};}
 if(stage.includes('Resolving joint sequence phase'))return{percent:82,label:'4 / 5 · Sequence groups'};
 if(stage.includes('Joint sequence phasing')){const pass=Number(stage.match(/pass (\d+)/)?.[1]||1);return{percent:82+Math.min(12,4*pass),label:'4 / 5 · Sequence groups · pass '+pass};}
 if(stage.includes('Jointly fitting free target pairings')){const pass=Number(stage.match(/pass (\d+)/)?.[1]||1);return{percent:65+Math.min(15,5*(pass-1)),label:'4 / 5 · Joint signal fit · pass '+pass};}
 if(stage.includes('Validating fitted conditional'))return{percent:97,label:'5 / 5 · Profile validation'};
 if(stage.includes('Registering offset traces'))return {percent:5,label:'1 / 4 · Trace registration'};
 if(stage.includes('Estimating indel length fractions'))return {percent:25,label:'2 / 4 · Length fractions'};
 if(stage.includes('Decoding offset')){const m=stage.match(/reference (\d+) \/ (\d+)/),done=m?Math.max(0,(Number(m[1])-1)/Number(m[2])):0;return {percent:40+40*done,label:'3 / 4 · Sequence groups'+(m?' · reference '+m[1]+' / '+m[2]:'')};}
 if(stage.includes('Validating offset profiles'))return {percent:95,label:'4 / 4 · Profile validation'};
 if(stage.includes('Validating fitted'))return {percent:97,label:'7 / 7 · Profile validation'};
 const round=Number(stage.match(/round (\d+)/)?.[1]||0);
 if(stage.includes('Rechecking sequence'))return {percent:92,label:'6 / 7 · Candidate recheck'};
 if(stage.includes('Jointly refitting'))return {percent:88,label:'5 / 7 · Joint signal fit'};
 if(stage.includes('substitution combinations'))return {percent:70,label:'4 / 7 · Substitution combinations'};
 if(stage.includes('sequence groups'))return {percent:20+Math.min(42,round/80*42),label:'3 / 7 · Sequence search'+(round?' · round '+round+' / 80':'')};
 if(/Calibrat|WT product|unedited pooled/.test(stage))return {percent:15,label:'2 / 7 · Control calibration'};
 if(/Fitting|separability/.test(stage))return {percent:35,label:'3 / 7 · Signal fitting'};
 return {percent:5,label:'1 / 7 · Trace registration'};
}
