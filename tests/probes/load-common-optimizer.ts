import {readFileSync} from 'node:fs';
import loadHighs from '../../public/vendor/highs-1.15.4/highs.mjs';
import initialize,* as conic from '../../public/vendor/clarabel-faer-0.11.1/clarabel.mjs';
import initializeFast,* as fast from '../../public/vendor/clarabel-qdldl-0.11.1/clarabel.mjs';
export default async function loadCommonOptimizer(){await initialize({module_or_path:readFileSync('public/vendor/clarabel-faer-0.11.1/clarabel_wasm_bg.wasm')});await initializeFast({module_or_path:readFileSync('public/vendor/clarabel-qdldl-0.11.1/clarabel_wasm_bg.wasm')});const engine=await loadHighs();engine.conic={solve:fast.solve,fallback:conic};return engine;}
