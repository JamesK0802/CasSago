// Global reference diagnostics only; never used to classify sequencing reads.
export type Opcode = ['equal' | 'replace' | 'delete' | 'insert', number, number, number, number];
export interface ReferenceComparison { opcodes: Opcode[]; similarity: number; }
const cache = new Map<string, ReferenceComparison>();
let cacheLimit = 256;
export function setComparisonCacheLimit(limit: number): void { cacheLimit = Math.max(256, limit); }
export const comparisonKey = (a: string, b: string) => JSON.stringify([a, b]);
export function rememberComparison(key: string, result: ReferenceComparison): void {
  if (cache.size >= cacheLimit && !cache.has(key)) cache.delete(cache.keys().next().value!);
  cache.set(key, result);
}
export function compareReferences(a: string, b: string): ReferenceComparison {
  const key = comparisonKey(a, b), saved = cache.get(key);
  if (saved) return saved;
  if (a > b) {
    const reverse = compareReferences(b, a);
    const result: ReferenceComparison = { similarity: reverse.similarity, opcodes: reverse.opcodes.map(
      ([tag, i, ii, j, jj]) => [tag === 'insert' ? 'delete' : tag === 'delete' ? 'insert' : tag, j, jj, i, ii]) };
    rememberComparison(key, result); return result;
  }
  if (a === b) {
    const result: ReferenceComparison = { opcodes: a.length ? [['equal', 0, a.length, 0, b.length]] : [], similarity: a.length ? 100 : 0 };
    rememberComparison(key, result); return result;
  }
  const width = b.length + 1;
  if ((a.length + 1) * width > 64_000_000) throw new Error('Reference comparison exceeds 64 million alignment cells. Use shorter reference sequences for this diagnostic view.');
  const trace = new Uint8Array((a.length + 1) * width);
  let previous = new Uint32Array(width);
  for (let j = 1; j < width; j++) { previous[j] = j; trace[j] = 2; }
  for (let i = 1; i <= a.length; i++) {
    const current = new Uint32Array(width); current[0] = i; trace[i * width] = 1;
    for (let j = 1; j < width; j++) {
      const diagonal = previous[j - 1] + Number(a[i - 1] !== b[j - 1]);
      const deletion = previous[j] + 1, insertion = current[j - 1] + 1;
      const best = Math.min(diagonal, deletion, insertion);
      current[j] = best;
      trace[i * width + j] = best === diagonal ? 0 : best === deletion ? 1 : 2;
    }
    previous = current;
  }
  let i = a.length, j = b.length;
  const steps: Opcode[] = [];
  while (i || j) {
    const direction = trace[i * width + j];
    if (direction === 1) { steps.push(['delete', i - 1, i, j, j]); i--; }
    else if (direction === 2) { steps.push(['insert', i, i, j - 1, j]); j--; }
    else { steps.push([a[i - 1] === b[j - 1] ? 'equal' : 'replace', i - 1, i, j - 1, j]); i--; j--; }
  }
  const opcodes: Opcode[] = [];
  let matches = 0;
  for (const step of steps.reverse()) {
    if (step[0] === 'equal') matches++;
    const last = opcodes[opcodes.length - 1];
    if (last && last[0] === step[0]) { last[2] = step[2]; last[4] = step[4]; }
    else opcodes.push([...step]);
  }
  const result = { opcodes, similarity: steps.length ? Math.round(matches / steps.length * 1000) / 10 : 0 };
  rememberComparison(key, result); return result;
}
export class ReferenceMatcher {
  constructor(_unused: unknown, private a: string, private b: string, _autojunk = false) {}
  getOpcodes(): Opcode[] { return compareReferences(this.a, this.b).opcodes; }
  ratio(): number { return compareReferences(this.a, this.b).similarity / 100; }
}
