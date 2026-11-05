/// <reference lib="webworker" />
import { designSangerPrimers } from './core/sanger-primer-designer';
addEventListener('message',({data}) => {
  try {
    const result = designSangerPrimers(data.references,data.settings,(fraction,stage) => postMessage({type:'progress',fraction,stage}));
    postMessage({type:'result',result});
  } catch (error) { postMessage({type:'error',message:error instanceof Error ? error.message : String(error)}); }
});
