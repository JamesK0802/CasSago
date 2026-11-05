/**
 * fastq-parser.ts — Streaming FASTQ parser for browser-side analysis.
 *
 * Reads FASTQ files using the File API without loading the entire file
 * into memory. Yields (sequence, quality_scores) tuples matching the
 * output of BioPython's SeqIO.parse(handle, "fastq").
 *
 * FASTQ format (4 lines per record):
 *   Line 1: @identifier
 *   Line 2: sequence
 *   Line 3: + (separator)
 *   Line 4: quality string (ASCII-encoded Phred scores)
 */

import { Gunzip, gunzipSync, strFromU8 } from 'fflate';

/** FASTQ Phred scores are bytes (0–93), not general-purpose numbers. */
export type QualityScores = number[] | Uint8Array;

export interface FastqRead {
  id?: string;
  seq: string;
  qual: QualityScores;
}

export function normalizeFastqReadId(identifier: string): string {
  const firstToken = identifier.trim().replace(/^@/, '').split(/\s+/)[0] || '';
  return firstToken.replace(/\/[12]$/, '');
}

/**
 * Parse an entire FASTQ file from a string.
 * Used in Web Workers where the file has already been read to text.
 */
export function parseFastqString(text: string): FastqRead[] {
  const results: FastqRead[] = [];
  const lines = text.split('\n');

  let i = 0;
  while (i < lines.length) {
    // Skip empty lines
    if (lines[i].trim() === '') {
      i++;
      continue;
    }

    // Line 1: identifier (starts with @)
    if (!lines[i].startsWith('@')) {
      i++;
      continue;
    }

    // Ensure we have all 4 lines
    if (i + 3 >= lines.length) break;

    // Line 2: sequence
    const seq = lines[i + 1].trim();

    // Line 3: separator (starts with +)
    // Line 4: quality string
    const qualStr = lines[i + 3].trim();

    if (seq.length > 0 && qualStr.length > 0) {
      // Convert ASCII quality to Phred scores (ASCII - 33)
      const qual = new Uint8Array(qualStr.length);
      for (let q = 0; q < qualStr.length; q++) {
        qual[q] = qualStr.charCodeAt(q) - 33;
      }
      results.push({ id: lines[i].trim().replace(/^@/, ''), seq, qual });
    }

    i += 4;
  }

  return results;
}

/**
 * Legacy text helper used only by the explicit “export filtered FASTQ” action.
 * Analysis itself uses `parseFastqFile()` above and never takes this path.
 */
export async function readFileAsText(file: File): Promise<string> {
  const buffer = await file.arrayBuffer();
  return /\.gz$/i.test(file.name)
    ? strFromU8(gunzipSync(new Uint8Array(buffer)))
    : new TextDecoder().decode(buffer);
}

/**
 * Read files in small chunks. The former whole-file `arrayBuffer()` path made
 * Safari keep the compressed file, decompressed text, line array, and parsed
 * reads alive together — enough to reload a tab during multi-file analysis.
 */
export async function parseFastqFileInBatches(
  file: File,
  onBatch: (reads: FastqRead[]) => Promise<void> | void,
  batchSize = 10000,
): Promise<number> {
  const decoder = new TextDecoder();
  const chunkSize = 1024 * 1024;
  let pending = '';
  let batch: FastqRead[] = [];
  let parsedCount = 0;
  const readyBatches: FastqRead[][] = [];

  const queueBatch = () => {
    if (!batch.length) return;
    readyBatches.push(batch);
    batch = [];
  };
  const flushReadyBatches = async () => {
    while (readyBatches.length) await onBatch(readyBatches.shift()!);
  };

  const consume = (text: string) => {
    pending += text;
    const lines = pending.split(/\r?\n/);
    pending = lines.pop() || '';
    let index = 0;
    while (index + 3 < lines.length) {
      if (!lines[index].trim() || !lines[index].startsWith('@')) {
        index++;
        continue;
      }
      const seq = lines[index + 1].trim();
      const quality = lines[index + 3].trim();
      if (seq && quality) {
        const qual = new Uint8Array(quality.length);
        for (let q = 0; q < quality.length; q++) qual[q] = quality.charCodeAt(q) - 33;
        batch.push({ id: lines[index].trim().replace(/^@/, ''), seq, qual });
        parsedCount++;
        if (batch.length >= batchSize) queueBatch();
      }
      index += 4;
    }
    // At most one incomplete four-line FASTQ record is kept for the next chunk.
    pending = lines.slice(index).join('\n') + (pending ? `\n${pending}` : '');
  };

  const emit = (chunk: Uint8Array, final = false) => {
    const text = decoder.decode(chunk, { stream: !final });
    if (text) consume(text);
  };
  const isGz = /\.gz$/i.test(file.name);
  const gunzip = isGz ? new Gunzip((chunk, final) => emit(chunk, final)) : null;

  for (let offset = 0; offset < file.size; offset += chunkSize) {
    const end = Math.min(offset + chunkSize, file.size);
    const chunk = new Uint8Array(await file.slice(offset, end).arrayBuffer());
    if (gunzip) gunzip.push(chunk, end === file.size);
    else emit(chunk, end === file.size);
    await flushReadyBatches();
  }
  if (pending) {
    const trailing = parseFastqString(pending);
    batch.push(...trailing);
    parsedCount += trailing.length;
  }
  queueBatch();
  await flushReadyBatches();
  return parsedCount;
}

/** Parse a FASTQ file when a complete in-memory read list is genuinely needed. */
export async function parseFastqFile(file: File): Promise<FastqRead[]> {
  const reads: FastqRead[] = [];
  await parseFastqFileInBatches(file, batch => { reads.push(...batch); });
  return reads;
}
