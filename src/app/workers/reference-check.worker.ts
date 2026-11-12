/// <reference lib="webworker" />
import { compareReferences, comparisonKey } from './core/reference-comparison';
addEventListener('message', ({ data }) => {
  try {
    const results = data.pairs.map(([a, b]: [string, string]) => [comparisonKey(a, b), compareReferences(a, b)]);
    postMessage({ id: data.id, results });
  } catch (error) { postMessage({ id: data.id, error: String(error) }); }
});
