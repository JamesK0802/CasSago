import {it,expect}from'vitest';
import {automaticSangerSettings}from'../../src/app/workers/core/sanger-automatic';
const base={references:[{name:'ref',sequence:'ACGT',cut:2}],maxDeletion:10,maxInsertion:10,sample:null,control:null}as any;
it('defaults inference to2 percent and preserves explicit zero and fractional thresholds',()=>{expect(automaticSangerSettings(base).indelThreshold).toBe(2);for(const indelThreshold of [0,.5,2,10,100])expect(automaticSangerSettings({...base,indelThreshold}).indelThreshold).toBe(indelThreshold);});
it('rejects invalid inference thresholds before discovery',()=>{for(const indelThreshold of [-1,100.1,NaN,Infinity])expect(()=>automaticSangerSettings({...base,indelThreshold})).toThrow('threshold');});
