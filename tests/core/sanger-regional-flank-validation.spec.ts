import {it,expect} from 'vitest';
import {regionalLengthInitialization,RegionalLengthWindow} from '../../src/app/workers/core/sanger-regional-length-fit';
const window=(positions:number[]):RegionalLengthWindow=>({classes:[{reference:0,size:0,column:positions.flatMap(()=>[1,0,0,0])}],observed:positions.flatMap(()=>[1,0,0,0])});
it('validates actual independent blocks rather than rejecting a short total flank',()=>{
 const positions=Array.from({length:9},(_,i)=>294+i),r=regionalLengthInitialization([window(positions)],[positions],[1],.005);
 expect(r.lengths[0][0]).toEqual([{size:0,fraction:1}]);expect(r.trainingPredictions[0]).toEqual(window(positions).observed);
});
it('rejects missing training or fewer than the existing five held-out positions',()=>{
 for(const positions of [[295,296,297,298,299],[294,295,296,297,298],[300,301,302,303,304,305,306,307,308,309]])expect(()=>regionalLengthInitialization([window(positions)],[positions],[1],.005)).toThrow('Insufficient independent flank validation.');
});
