/// <reference lib="webworker" />
import {designGuides} from './core/grna-designer';
addEventListener('message',({data})=>{try{const result=designGuides(data.references,data.backgrounds,data.settings,stage=>postMessage({type:'progress',stage}));postMessage({type:'result',result});}catch(e){postMessage({type:'error',message:e instanceof Error?e.message:String(e)});}});
