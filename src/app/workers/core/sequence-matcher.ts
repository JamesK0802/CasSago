/**
 * sequence-matcher.ts — TypeScript port of Python's difflib.SequenceMatcher.
 *
 * Faithfully implements the Ratcliffe/Obershelp algorithm used by CPython's
 * difflib module.  Only the subset required by the CRISPR analysis pipeline
 * is ported: get_opcodes(), get_matching_blocks(), and ratio().
 */

export type Opcode = [tag: string, i1: number, i2: number, j1: number, j2: number];

export class SequenceMatcher {
  private a: string;
  private b: string;
  private b2j: Map<string, number[]> = new Map();
  private matchingBlocksCache: [number, number, number][] | null = null;
  private opcodesCache: Opcode[] | null = null;
  private fullBCount: Map<string, number> | null = null;
  private bjunk: Set<string> = new Set();
  private bpopular: Set<string> = new Set();
  private autojunk: boolean;

  constructor(isjunk: null, a: string = '', b: string = '', autojunk: boolean = true) {
    this.a = a;
    this.b = b;
    this.autojunk = autojunk;
    this.chainB();
  }

  setSeqs(a: string, b: string): void {
    this.setSeq1(a);
    this.setSeq2(b);
  }

  setSeq1(a: string): void {
    if (a === this.a) return;
    this.a = a;
    this.matchingBlocksCache = null;
    this.opcodesCache = null;
  }

  setSeq2(b: string): void {
    if (b === this.b) return;
    this.b = b;
    this.matchingBlocksCache = null;
    this.opcodesCache = null;
    this.chainB();
  }

  private chainB(): void {
    const b = this.b;
    const n = b.length;
    this.b2j = new Map();
    const b2j = this.b2j;

    for (let i = 0; i < n; i++) {
      const c = b[i];
      const indices = b2j.get(c);
      if (indices) {
        indices.push(i);
      } else {
        b2j.set(c, [i]);
      }
    }

    // Purge junk elements (none in our usage, isjunk is always null)
    this.bjunk = new Set();

    // Purge popular elements (autojunk heuristic)
    this.bpopular = new Set();
    if (this.autojunk && n >= 200) {
      const ntest = Math.floor(n / 100) + 1;
      for (const [elt, idxs] of b2j) {
        if (idxs.length > ntest) {
          this.bpopular.add(elt);
        }
      }
      for (const elt of this.bpopular) {
        b2j.delete(elt);
      }
    }
  }

  /**
   * Find longest matching block in a[alo:ahi] and b[blo:bhi].
   * Returns [i, j, k] such that a[i:i+k] == b[j:j+k].
   */
  findLongestMatch(alo: number, ahi: number, blo: number, bhi: number): [number, number, number] {
    const a = this.a;
    const b2j = this.b2j;

    let besti = alo;
    let bestj = blo;
    let bestsize = 0;

    // j2len[j] = length of longest match ending with a[i-1] and b[j]
    let j2len: Map<number, number> = new Map();

    for (let i = alo; i < ahi; i++) {
      const newj2len: Map<number, number> = new Map();
      const indices = b2j.get(a[i]);
      if (indices) {
        for (const j of indices) {
          if (j < blo) continue;
          if (j >= bhi) break;
          const k = (j2len.get(j - 1) || 0) + 1;
          newj2len.set(j, k);
          if (k > bestsize) {
            besti = i - k + 1;
            bestj = j - k + 1;
            bestsize = k;
          }
        }
      }
      j2len = newj2len;
    }

    // Extend the match as far as possible (handles junk)
    while (besti > alo && bestj > blo && !this.bjunk.has(this.b[bestj - 1]) && a[besti - 1] === this.b[bestj - 1]) {
      besti--;
      bestj--;
      bestsize++;
    }
    while (besti + bestsize < ahi && bestj + bestsize < bhi && !this.bjunk.has(this.b[bestj + bestsize]) && a[besti + bestsize] === this.b[bestj + bestsize]) {
      bestsize++;
    }

    // Extend with junk on both sides
    while (besti > alo && bestj > blo && this.bjunk.has(this.b[bestj - 1]) && a[besti - 1] === this.b[bestj - 1]) {
      besti--;
      bestj--;
      bestsize++;
    }
    while (besti + bestsize < ahi && bestj + bestsize < bhi && this.bjunk.has(this.b[bestj + bestsize]) && a[besti + bestsize] === this.b[bestj + bestsize]) {
      bestsize++;
    }

    return [besti, bestj, bestsize];
  }

  getMatchingBlocks(): [number, number, number][] {
    if (this.matchingBlocksCache !== null) return this.matchingBlocksCache;

    const la = this.a.length;
    const lb = this.b.length;

    const queue: [number, number, number, number][] = [[0, la, 0, lb]];
    const matching: [number, number, number][] = [];

    while (queue.length > 0) {
      const [alo, ahi, blo, bhi] = queue.pop()!;
      const [i, j, k] = this.findLongestMatch(alo, ahi, blo, bhi);
      if (k > 0) {
        matching.push([i, j, k]);
        if (alo < i && blo < j) queue.push([alo, i, blo, j]);
        if (i + k < ahi && j + k < bhi) queue.push([i + k, ahi, j + k, bhi]);
      }
    }

    matching.sort((a, b) => a[0] - b[0] || a[1] - b[1] || a[2] - b[2]);

    // Collapse adjacent equal blocks
    let i1 = 0, j1 = 0, k1 = 0;
    const nonAdjacent: [number, number, number][] = [];
    for (const [i2, j2, k2] of matching) {
      if (i1 + k1 === i2 && j1 + k1 === j2) {
        k1 += k2;
      } else {
        if (k1) nonAdjacent.push([i1, j1, k1]);
        i1 = i2;
        j1 = j2;
        k1 = k2;
      }
    }
    if (k1) nonAdjacent.push([i1, j1, k1]);

    nonAdjacent.push([la, lb, 0]); // sentinel
    this.matchingBlocksCache = nonAdjacent;
    return nonAdjacent;
  }

  getOpcodes(): Opcode[] {
    if (this.opcodesCache !== null) return this.opcodesCache;

    let i = 0, j = 0;
    const answer: Opcode[] = [];

    for (const [ai, bj, size] of this.getMatchingBlocks()) {
      let tag = '';
      if (i < ai && j < bj) tag = 'replace';
      else if (i < ai) tag = 'delete';
      else if (j < bj) tag = 'insert';

      if (tag) answer.push([tag, i, ai, j, bj]);
      i = ai + size;
      j = bj + size;
      if (size) answer.push(['equal', ai, ai + size, bj, bj + size]);
    }

    this.opcodesCache = answer;
    return answer;
  }

  ratio(): number {
    const matches = this.getMatchingBlocks().reduce((sum, [, , k]) => sum + k, 0);
    const length = this.a.length + this.b.length;
    return length ? (2.0 * matches) / length : 1.0;
  }
}
