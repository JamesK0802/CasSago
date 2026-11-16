import {it,expect,vi} from 'vitest';
const state=vi.hoisted(()=>({sessions:new WeakMap<object,object>(),analyze:vi.fn(),truth:vi.fn()}));
vi.mock('../../src/app/workers/core/sanger-marginal-engine',()=>({marginalSessions:state.sessions}));
vi.mock('../../src/app/workers/core/sanger-automatic',()=>({analyzeAutomaticSanger:state.analyze}));
vi.mock('../../src/app/workers/core/sanger-marginal-truth-retrieval',()=>({retrieveMarginalTruth:state.truth}));
vi.mock('../../src/app/workers/core/sanger-marginal-alternatives',()=>({marginalAlternativePage:vi.fn()}));
vi.mock('../../src/app/utils/sanger-native-signal-certificate',()=>({certifyNativeBenchmark:vi.fn()}));
vi.mock('../../src/app/utils/sanger-benchmark-assessment',()=>({assessSangerBenchmark:vi.fn()}));
it('releases the previous physical model before replacement and rejects stale queries after failure',async()=>{
 let receive:any;const post=vi.fn(),first={},second={};
 vi.stubGlobal('addEventListener',(_event:string,callback:any)=>receive=callback);vi.stubGlobal('postMessage',post);
 try{
  await import('../../src/app/workers/sanger-analysis.worker');
  state.analyze.mockImplementationOnce(()=>{state.sessions.set(first,{});return first;});
  await receive({data:{action:'automatic'}});expect(state.sessions.has(first)).toBe(true);
  state.analyze.mockImplementationOnce(()=>{expect(state.sessions.has(first)).toBe(false);state.sessions.set(second,{});return second;});
  await receive({data:{action:'automatic'}});expect(state.sessions.has(second)).toBe(true);
  state.analyze.mockImplementationOnce(()=>{expect(state.sessions.has(second)).toBe(false);throw Error('Replacement failed');});
  await receive({data:{action:'automatic'}});await receive({data:{action:'marginal-truth',requestId:7}});
  expect(state.truth).not.toHaveBeenCalled();expect(post).toHaveBeenLastCalledWith(expect.objectContaining({type:'marginal-truth-error',requestId:7}));
 }finally{vi.unstubAllGlobals();}
});
