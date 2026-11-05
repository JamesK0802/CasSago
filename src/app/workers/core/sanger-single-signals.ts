import {SangerTrace,TRACE_BASES} from '../../models/sanger-analysis';
import {contrastPeaks,integratedPeaks,calibrateSignal,correctGain} from './sanger-signal';

/** Select a matched fluorescence representation using intact upstream bases only.
 * Contrast removes background but assumes comparable dye peak shapes. If that
 * representation fails the existing upstream checks, integrated area preserves
 * unequal dye shapes. Never choose a representation from an edited fit or truth.
 */
export function matchedSingleSignals(control:SangerTrace,sample:SangerTrace,positions:number[],offset:number){
  const evaluate=(method:'contrast'|'area')=>{
    const extract=method==='contrast'?contrastPeaks:integratedPeaks;
    const controlSignal=extract(control),sampleSignal=extract(sample),wt=controlSignal.vectors,raw=sampleSignal.vectors;
    const calibration=calibrateSignal(wt.map((_,i)=>raw[i+offset]||[0,0,0,0]),positions,[wt]);
    const edited=raw.map(v=>correctGain(v,calibration.gains));
    const upstreamError=positions.reduce((sum,i)=>sum+wt[i].reduce((s,v,b)=>s+Math.abs(v-edited[i+offset][b]),0)/2,0)/positions.length;
    return {method,controlSignal,sampleSignal,wt,raw,edited,calibration,upstreamError};
  };
  const contrast=evaluate('contrast');
  const passes=(r:ReturnType<typeof evaluate>)=>r.calibration.validationR2>=.8&&r.upstreamError<=.1;
  const area=passes(contrast)?undefined:evaluate('area');
  const chosen=area&&passes(area)?area:contrast;
  // An area-based WT template may have genuine spectral overlap. Assess loss of
  // information relative to that control's own upstream dye dominance.
  const expected=TRACE_BASES.map((base,b)=>{
    const values=positions.filter(i=>control.sequence[i]===base).map(i=>chosen.wt[i][b]).sort((a,b)=>a-b);
    return chosen.method==='area'&&values.length>=8?values[Math.floor(values.length/2)]:1;
  });
  const templateUsable=chosen.wt.map((_,i)=>{
    const lo=Math.max(0,i-12),hi=Math.min(chosen.wt.length,i+13);
    return chosen.wt.slice(lo,hi).reduce((sum,v,k)=>{const b=TRACE_BASES.indexOf(control.sequence[lo+k] as typeof TRACE_BASES[number]);
      return sum+(b<0?0:Math.min(1,v[b]/Math.max(.01,expected[b])));},0)/(hi-lo)>=.9;
  });
  return {...chosen,templateUsable};
}
