import {Component,Input,Output,EventEmitter,OnChanges,OnDestroy,ChangeDetectorRef} from '@angular/core';
import {CommonModule} from '@angular/common';import {FormsModule} from '@angular/forms';
import {findGrnaCutSite,extractWindow,cutIndexInWindow} from '../../workers/core/classifier';
import {findHomoeologGuideSite,homoeologSimilarity,projectHomoeologPosition,HomoeologGuideMatch} from '../../workers/core/homoeolog';
import {ReferenceMatcher as SequenceMatcher,ReferenceComparison,compareReferences,rememberComparison,setComparisonCacheLimit} from '../../workers/core/reference-comparison';
export interface ExtractedWindowItem {
  /** Stable FormArray position used internally; names may be edited or duplicated. */
  geneIndex: number;
  targetIndex: number;
  geneName: string;
  targetId: string;
  sequence: string;
  cutSiteIndex: number;
  cutSiteInWindow: number;
  pam: string;
  strand: string;

  // Visual Reference Map properties
  refLength?: number;
  grnaStart?: number;
  grnaLength?: number;
  grnaEnd?: number;
  winStart?: number;
  winEnd?: number;
  winLeftPercent?: number;
  winWidthPercent?: number;
  targetLeftPercent?: number;
  targetWidthPercent?: number;
  guideMatch?: HomoeologGuideMatch;
}

export interface AlignmentCharToken {
  char: string;
  isMatch: boolean;
  isGap: boolean;
  isCutSite?: boolean;
}

export interface PairwiseComparisonData {
  target1: ExtractedWindowItem;
  target2: ExtractedWindowItem;
  similarity: number;
  tokens1: AlignmentCharToken[];
  tokens2: AlignmentCharToken[];
  matchBar: string[];
  matchCount: number;
  mismatchCount: number;
  totalLen: number;
}

export interface HomoeologOverviewMarker {
  targetId: string;
  rowIndex: number;
  leftPercent?: number;
  windowLeftPercent?: number;
  windowWidthPercent?: number;
  guideLeftPercent?: number;
  guideWidthPercent?: number;
  cutSite?: number;
  exact: boolean;
  mismatches: number;
  identity: number;
  windowDifferenceBp?: number;
  error?: string;
}

export interface HomoeologDifferenceSegment {
  leftPercent: number;
  widthPercent: number;
  label: string;
}

export interface HomoeologOverviewLane {
  geneIndex: number;
  geneName: string;
  isAnchor: boolean;
  similarity: number;
  aligned: boolean;
  markers: HomoeologOverviewMarker[];
  differences: HomoeologDifferenceSegment[];
}

export interface HomoeologOverview {
  groupId: string;
  anchorName: string;
  targetIds: string[];
  lanes: HomoeologOverviewLane[];
}


@Component({selector:'app-sanger-window-check',standalone:true,imports:[CommonModule,FormsModule],templateUrl:'./sanger-window-check.component.html',styles:[':host{display:block;min-width:0}'],styleUrls:['../../pages/analysis-page/analysis-page.component.css']})
export class SangerWindowCheckComponent implements OnChanges,OnDestroy {

 @Input() references:{name:string;sequence:string;targets:string;mode:'cuts'|'guides';uiTargets?:{name:string;value:string}[]}[]=[];
 @Input() customWindowEnabled=false;@Input() customWindowLeft=45;@Input() customWindowRight=45;@Input() windowCheckSize=90;
 @Output() close=new EventEmitter<void>();showWindowCheck=true;homoeologMode=true;
 isCalculatingWindowCheck=false;windowCheckError='';extractedWindows:ExtractedWindowItem[]=[];similarityMatrix:number[][]=[];homoeologOverviews:HomoeologOverview[]=[];
 selectedPairRowIndex:number|null=null;selectedPairColIndex:number|null=null;selectedPairComparison:PairwiseComparisonData|null=null;
 private homoeologAnchorByGroup=new Map<string,number>();private windowCheckTimer?:ReturnType<typeof setTimeout>;private windowCheckRevision=0;
 private referenceCheckWorker?:Worker;private referenceRequestId=0;private referenceRequests=new Map<number,{resolve:()=>void;reject:(error:Error)=>void}>();
 isDraggingScroll=false;startX=0;scrollLeft=0;
 constructor(private cdr:ChangeDetectorRef){}
 private get genes(){return this.references.map((row,i)=>({gene_name:row.name,gene_reference:row.sequence,homoeolog_group:'Homoeolog 1',geneTargets:(row.uiTargets||row.targets.split(/[,;\s]+/).filter(Boolean).map((value,t)=>({name:'T'+(t+1),value}))).map(t=>({target_id:t.name,gRNA:t.value,...(row.mode==='cuts'?{importedCut:Number(t.value)}:{})}))}));}
 private state={analysisForm:{get:(name:string):{value:any}=>({value:name==='genes'?this.genes:name==='customWindowEnabled'?this.customWindowEnabled:name==='customWindowLeft'?this.customWindowLeft:name==='customWindowRight'?this.customWindowRight:this.windowCheckSize})}};
 ngOnChanges(){this.recalculateWindowCheck();}
 selectHomoeologAnchor(groupId:string,geneIndex:number){this.homoeologAnchorByGroup.set(groupId,geneIndex);this.recalculateWindowCheck();}
 ngOnDestroy(){if(this.windowCheckTimer)clearTimeout(this.windowCheckTimer);this.windowCheckRevision++;this.referenceCheckWorker?.terminate();this.referenceRequests.clear();}
  onWindowCheckSizeChange(newSize: any) {
    const parsed = Number(newSize);
    this.windowCheckSize = isNaN(parsed) || parsed < 1 ? 90 : Math.min(1000, Math.round(parsed));
    if (this.showWindowCheck) {
      this.recalculateWindowCheck();
    }
  }

  selectPairComparison(rIdx: number, cIdx: number) {
    if (rIdx < 0 || rIdx >= this.extractedWindows.length || cIdx < 0 || cIdx >= this.extractedWindows.length) {
      return;
    }
    this.selectedPairRowIndex = rIdx;
    this.selectedPairColIndex = cIdx;
    this.updatePairComparison();
  }

  clearPairComparison() {
    this.selectedPairRowIndex = null;
    this.selectedPairColIndex = null;
    this.selectedPairComparison = null;
  }

  // ── Drag to scroll handler for sequence alignment & heatmap table ───────
  startDragScroll(e: MouseEvent, element: HTMLElement) {
    this.isDraggingScroll = true;
    element.classList.add('dragging-scroll');
    this.startX = e.pageX - element.offsetLeft;
    this.scrollLeft = element.scrollLeft;
  }

  stopDragScroll(element: HTMLElement) {
    this.isDraggingScroll = false;
    element.classList.remove('dragging-scroll');
  }

  onDragScroll(e: MouseEvent, element: HTMLElement) {
    if (!this.isDraggingScroll) return;
    e.preventDefault();
    const x = e.pageX - element.offsetLeft;
    const walk = (x - this.startX) * 1.5;
    element.scrollLeft = this.scrollLeft - walk;
  }

  updatePairComparison() {
    if (this.selectedPairRowIndex === null || this.selectedPairColIndex === null) {
      this.selectedPairComparison = null;
      return;
    }
    const rIdx = this.selectedPairRowIndex;
    const cIdx = this.selectedPairColIndex;

    if (rIdx >= this.extractedWindows.length || cIdx >= this.extractedWindows.length) {
      this.selectedPairComparison = null;
      return;
    }

    const t1 = this.extractedWindows[rIdx];
    const t2 = this.extractedWindows[cIdx];
    const seq1 = t1.sequence;
    const seq2 = t2.sequence;
    const similarity = this.similarityMatrix[rIdx]?.[cIdx] ?? 0;

    const matcher = new SequenceMatcher(null, seq1, seq2, seq2.length < 200);
    const opcodes = matcher.getOpcodes();

    const tokens1: AlignmentCharToken[] = [];
    const tokens2: AlignmentCharToken[] = [];
    const matchBar: string[] = [];
    let matchCount = 0;
    let mismatchCount = 0;

    for (const [tag, i1, i2, j1, j2] of opcodes) {
      if (tag === 'equal') {
        for (let k = 0; k < (i2 - i1); k++) {
          const idx1 = i1 + k;
          const idx2 = j1 + k;
          const c1 = seq1[idx1];
          const c2 = seq2[idx2];
          const isCut1 = idx1 === t1.cutSiteInWindow;
          const isCut2 = idx2 === t2.cutSiteInWindow;
          tokens1.push({ char: c1, isMatch: true, isGap: false, isCutSite: isCut1 });
          tokens2.push({ char: c2, isMatch: true, isGap: false, isCutSite: isCut2 });
          matchBar.push((isCut1 || isCut2) ? '✂' : '|');
          matchCount++;
        }
      } else if (tag === 'replace') {
        const len1 = i2 - i1;
        const len2 = j2 - j1;
        const maxLen = Math.max(len1, len2);
        for (let k = 0; k < maxLen; k++) {
          const idx1 = i1 + k;
          const idx2 = j1 + k;
          const c1 = k < len1 ? seq1[idx1] : '-';
          const c2 = k < len2 ? seq2[idx2] : '-';
          const isM = c1 === c2;
          const isCut1 = k < len1 && idx1 === t1.cutSiteInWindow;
          const isCut2 = k < len2 && idx2 === t2.cutSiteInWindow;
          tokens1.push({ char: c1, isMatch: isM, isGap: c1 === '-', isCutSite: isCut1 });
          tokens2.push({ char: c2, isMatch: isM, isGap: c2 === '-', isCutSite: isCut2 });
          matchBar.push((isCut1 || isCut2) ? '✂' : (isM ? '|' : '•'));
          if (isM) matchCount++; else mismatchCount++;
        }
      } else if (tag === 'delete') {
        for (let k = 0; k < (i2 - i1); k++) {
          const idx1 = i1 + k;
          const c1 = seq1[idx1];
          const isCut1 = idx1 === t1.cutSiteInWindow;
          tokens1.push({ char: c1, isMatch: false, isGap: false, isCutSite: isCut1 });
          tokens2.push({ char: '-', isMatch: false, isGap: true });
          matchBar.push(isCut1 ? '✂' : '•');
          mismatchCount++;
        }
      } else if (tag === 'insert') {
        for (let k = 0; k < (j2 - j1); k++) {
          const idx2 = j1 + k;
          const c2 = seq2[idx2];
          const isCut2 = idx2 === t2.cutSiteInWindow;
          tokens1.push({ char: '-', isMatch: false, isGap: true });
          tokens2.push({ char: c2, isMatch: false, isGap: false, isCutSite: isCut2 });
          matchBar.push(isCut2 ? '✂' : '•');
          mismatchCount++;
        }
      }
    }

    this.selectedPairComparison = {
      target1: t1,
      target2: t2,
      similarity,
      tokens1,
      tokens2,
      matchBar,
      matchCount,
      mismatchCount,
      totalLen: tokens1.length
    };
  }

  recalculateWindowCheck() {
    if (this.windowCheckTimer) clearTimeout(this.windowCheckTimer);
    const revision = ++this.windowCheckRevision;
    this.isCalculatingWindowCheck = true;
    this.windowCheckError = '';
    this.cdr.detectChanges();

    this.windowCheckTimer = setTimeout(async () => {
      this.windowCheckTimer = undefined;
      try {
        const genesFormVal = this.state.analysisForm.get('genes')?.value || [];
        const customWindowEnabled = Boolean(this.state.analysisForm.get('customWindowEnabled')?.value);
        const customWindowLeft = Math.max(0, Number(this.state.analysisForm.get('customWindowLeft')?.value) || 0);
        const customWindowRight = Math.max(0, Number(this.state.analysisForm.get('customWindowRight')?.value) || 0);
        const displayWindowSize = customWindowEnabled ? customWindowLeft + customWindowRight : this.windowCheckSize;
        const extracted: ExtractedWindowItem[] = [];

        genesFormVal.forEach((g: any, gi: number) => {
          const geneName = g.gene_name?.trim() || `Gene ${gi + 1}`;
          const refSeq = (g.gene_reference || '').replace(/\s+/g, '').toUpperCase();
          if (!refSeq) return;

          (g.geneTargets || []).forEach((t: any, ti: number) => {
            const targetId = t.target_id?.trim() || `T${ti + 1}`;
            const grna = (t.gRNA || '').trim().toUpperCase();
            if (!grna) return;

            const importedCut=typeof t.importedCut==='number'?t.importedCut:undefined;
            const guideMatch = importedCut!==undefined?{matched:true,exact:true,mismatches:0,identity:100,strand:'forward' as const,grnaStart:Math.max(0,importedCut-17),grnaEnd:Math.min(refSeq.length,importedCut+3),cutSite:importedCut,pam:'',pamFound:false,referenceSpacer:'',analysisGuide:''}:findHomoeologGuideSite(refSeq, grna);
            const exactCutInfo = guideMatch ? null : findGrnaCutSite(refSeq, grna);
            const grnaStart = guideMatch ? guideMatch.grnaStart : exactCutInfo!.grna_start;
            const strand = guideMatch ? guideMatch.strand : exactCutInfo!.strand;
            const pam = guideMatch ? guideMatch.pam : exactCutInfo!.pam;
            let cutSite = guideMatch ? guideMatch.cutSite : exactCutInfo!.cut_site;
            if (grnaStart < 0 || cutSite < 0 || cutSite >= refSeq.length) {
              cutSite = Math.floor(refSeq.length / 2);
            }

            const winSeq = extractWindow(refSeq, cutSite, displayWindowSize,
              customWindowEnabled ? customWindowLeft : undefined,
              customWindowEnabled ? customWindowRight : undefined);
            const cutWinIdx = grnaStart !== -1 ? cutIndexInWindow(refSeq, cutSite, displayWindowSize,
              customWindowEnabled ? customWindowLeft : undefined,
              customWindowEnabled ? customWindowRight : undefined) : -1;

            const refLength = refSeq.length;
            const grnaLength = grna.length;
            const grnaEnd = grnaStart >= 0 ? grnaStart + grnaLength : -1;

            let winStart = 0;
            let winEnd = 0;
            if (cutSite >= 0) {
              const left = customWindowEnabled ? customWindowLeft : Math.floor(displayWindowSize / 2);
              const right = customWindowEnabled ? customWindowRight : Math.floor(displayWindowSize / 2);
              winStart = Math.max(0, cutSite - left);
              winEnd = Math.min(refLength, cutSite + right);
            }

            const winLeftPercent = refLength > 0 ? (winStart / refLength) * 100 : 0;
            const winWidthPercent = refLength > 0 ? ((winEnd - winStart) / refLength) * 100 : 0;
            const targetLeftPercent = (grnaStart >= 0 && refLength > 0) ? (grnaStart / refLength) * 100 : 0;
            const targetWidthPercent = (grnaStart >= 0 && refLength > 0) ? (grnaLength / refLength) * 100 : 0;

            extracted.push({
              geneIndex: gi,
              targetIndex: ti,
              geneName,
              targetId,
              sequence: winSeq,
              cutSiteIndex: grnaStart !== -1 ? cutSite : -1,
              cutSiteInWindow: cutWinIdx,
              pam,
              strand,
              refLength,
              grnaStart: grnaStart >= 0 ? grnaStart : undefined,
              grnaLength,
              grnaEnd: grnaEnd >= 0 ? grnaEnd : undefined,
              winStart,
              winEnd,
              winLeftPercent,
              winWidthPercent,
              targetLeftPercent: grnaStart >= 0 ? targetLeftPercent : undefined,
              targetWidthPercent: grnaStart >= 0 ? targetWidthPercent : undefined,
              guideMatch
            });
          });
        });

        // A delete/reorder can trigger a newer form update while a previous
        // calculation is queued. Only let the current FormArray snapshot
        // update the cards, matrix, and homoeolog lanes together.
        if (revision !== this.windowCheckRevision) return;
        const sequenceSets = [extracted.map(item => item.sequence),
          this.homoeologMode ? genesFormVal.map((gene: any) => String(gene.gene_reference || '').replace(/\s+/g, '').toUpperCase()) : []];
        await this.prepareReferenceComparisons(sequenceSets);
        if (revision !== this.windowCheckRevision) return;
        this.extractedWindows = extracted;
        const n = extracted.length;
        const matrix: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));

        for (let i = 0; i < n; i++) {
          matrix[i][i] = 100.0;
          for (let j = i + 1; j < n; j++) {
            const sim = this.calculateSymmetricSimilarity(extracted[i].sequence, extracted[j].sequence);
            matrix[i][j] = sim;
            matrix[j][i] = sim;
          }
        }
        this.similarityMatrix = matrix;
        this.homoeologOverviews = this.homoeologMode
          ? this.buildHomoeologOverviews(genesFormVal, extracted)
          : [];

        if (this.selectedPairRowIndex !== null && this.selectedPairColIndex !== null) {
          this.updatePairComparison();
        }
      } catch (err) {
        console.error('Error calculating window check matrix:', err);
        if (revision === this.windowCheckRevision) this.windowCheckError = err instanceof Error ? err.message : String(err);
      } finally {
        if (revision === this.windowCheckRevision) this.isCalculatingWindowCheck = false;
        this.cdr.detectChanges();
      }
    }, 150);
  }

  private prepareReferenceComparisons(sequenceSets: string[][]): Promise<void> {
    if (!this.referenceCheckWorker) {
      this.referenceCheckWorker = new Worker(new URL('../../workers/reference-check.worker', import.meta.url), { type: 'module' });
      this.referenceCheckWorker.onmessage = ({ data }) => {
        const request = this.referenceRequests.get(data.id);
        if (!request) return;
        this.referenceRequests.delete(data.id);
        if (data.error) request.reject(new Error(data.error));
        else {
          for (const [key, result] of data.results as [string, ReferenceComparison][]) rememberComparison(key, result);
          request.resolve();
        }
      };
      this.referenceCheckWorker.onerror = () => {
        this.referenceRequests.forEach(request => request.reject(new Error('Reference comparison worker failed')));
        this.referenceRequests.clear();
        this.referenceCheckWorker?.terminate(); this.referenceCheckWorker = undefined;
      };
    }
    const pairs: [string, string][] = [];
    for (const values of sequenceSets) {
      const sequences = [...new Set(values)].filter(Boolean);
      for (const a of sequences) for (const b of sequences) pairs.push([a, b]);
    }
    const id = ++this.referenceRequestId;
    setComparisonCacheLimit(pairs.length * 2);
    return new Promise((resolve, reject) => {
      this.referenceRequests.set(id, { resolve, reject });
      this.referenceCheckWorker!.postMessage({ id, pairs });
    });
  }

  calculateSymmetricSimilarity(seq1: string, seq2: string): number {
    if (!seq1 || !seq2) return 0;
    if (seq1 === seq2) return 100.0;

    return compareReferences(seq1, seq2).similarity;
  }

  private buildHomoeologOverviews(genes: any[], windows: ExtractedWindowItem[]): HomoeologOverview[] {
    const grouped = new Map<string, any[]>();
    for (const [geneIndex, gene] of genes.entries()) {
      const sequence = String(gene.gene_reference || '').replace(/\s+/g, '').toUpperCase();
      if (!sequence) continue;
      const groupId = String(gene.homoeolog_group || '').trim() || `Homoeolog ${geneIndex + 1}`;
      const list = grouped.get(groupId) || [];
      list.push({ ...gene, geneIndex, sequence, geneName: String(gene.gene_name || `Ref${list.length + 1}`).trim() });
      grouped.set(groupId, list);
    }

    return [...grouped.entries()].map(([groupId, refs]) => {
      const fallbackAnchor = [...refs].sort((a, b) => b.sequence.length - a.sequence.length)[0];
      const selectedAnchorIndex = this.homoeologAnchorByGroup.get(groupId);
      const anchor = refs.find(ref => ref.geneIndex === selectedAnchorIndex) || fallbackAnchor;
      this.homoeologAnchorByGroup.set(groupId, anchor.geneIndex);
      const targetIds = [...new Set(refs.flatMap(ref => (ref.geneTargets || []).map((t: any, i: number) => String(t.target_id || `T${i + 1}`).trim())))];
      const lanes: HomoeologOverviewLane[] = refs.map(ref => {
        const similarity = homoeologSimilarity(ref.sequence, anchor.sequence);
        const aligned = ref === anchor || similarity >= 45;
        const differences = aligned && ref !== anchor
          ? new SequenceMatcher(null, ref.sequence, anchor.sequence, false).getOpcodes()
            .filter(([tag]) => tag !== 'equal')
            .map(([tag, i1, i2, j1, j2]) => ({
              leftPercent: Math.max(0, Math.min(100, (j1 / Math.max(1, anchor.sequence.length)) * 100)),
              widthPercent: Math.max(0.55, ((j2 - j1) / Math.max(1, anchor.sequence.length)) * 100),
              label: `${tag} · ref ${i1 + 1}-${i2} · alignment ${j1 + 1}-${j2}`
            }))
          : [];
        const markers = targetIds.map((targetId, rowIndex) => {
          const item = windows.find(window => window.geneIndex === ref.geneIndex && window.targetId === targetId);
          const anchorItem = windows.find(window => window.geneIndex === anchor.geneIndex && window.targetId === targetId);
          if (!item || item.cutSiteIndex < 0 || !aligned) {
            return {
              targetId, rowIndex, exact: false, mismatches: item?.guideMatch?.mismatches ?? 0,
              identity: item?.guideMatch?.identity ?? 0,
              error: !aligned ? `Reference alignment ${similarity.toFixed(1)}%` : (item?.guideMatch?.error || 'Guide not matched')
            };
          }
          const guideStart = item.grnaStart ?? item.cutSiteIndex;
          const guideEnd = guideStart + (item.grnaLength || 1);
          const projectedGuideStart = projectHomoeologPosition(ref.sequence, anchor.sequence, guideStart);
          const projectedGuideEnd = projectHomoeologPosition(ref.sequence, anchor.sequence, guideEnd);
          const projectedWindowStart = projectHomoeologPosition(ref.sequence, anchor.sequence, item.winStart ?? 0);
          const projectedWindowEnd = projectHomoeologPosition(ref.sequence, anchor.sequence, item.winEnd ?? ref.sequence.length);
          const projected = projectedGuideStart === null || projectedGuideEnd === null
            ? null
            : (projectedGuideStart + projectedGuideEnd) / 2;
          const percentageRange = (start: number | null, end: number | null) => {
            if (start === null || end === null) return {};
            const left = Math.min(start, end);
            const right = Math.max(start, end);
            return {
              left: Math.max(0, Math.min(100, (left / Math.max(1, anchor.sequence.length)) * 100)),
              width: Math.max(0.55, ((right - left) / Math.max(1, anchor.sequence.length)) * 100)
            };
          };
          const guideRange = percentageRange(projectedGuideStart, projectedGuideEnd);
          const windowRange = percentageRange(projectedWindowStart, projectedWindowEnd);
          return {
            targetId,
            rowIndex,
            leftPercent: projected === null ? undefined : Math.max(0, Math.min(100, (projected / Math.max(1, anchor.sequence.length)) * 100)),
            guideLeftPercent: guideRange.left,
            guideWidthPercent: guideRange.width,
            windowLeftPercent: windowRange.left,
            windowWidthPercent: windowRange.width,
            cutSite: item.cutSiteIndex,
            exact: item.guideMatch?.exact ?? true,
            mismatches: item.guideMatch?.mismatches ?? 0,
            identity: item.guideMatch?.identity ?? 100,
            windowDifferenceBp: anchorItem ? this.sequenceDifferenceBp(item.sequence, anchorItem.sequence) : undefined,
            error: projected === null ? 'Alignment position unavailable' : undefined
          };
        });
        return {
          geneIndex: ref.geneIndex,
          geneName: ref.geneName,
          isAnchor: ref.geneIndex === anchor.geneIndex,
          similarity,
          aligned,
          markers,
          differences
        };
      });
      return { groupId, anchorName: anchor.geneName, targetIds, lanes };
    });
  }

  private sequenceDifferenceBp(left: string, right: string): number {
    if (!left || !right) return 0;
    if (left === right) return 0;
    return new SequenceMatcher(null, left, right, false).getOpcodes()
      .filter(([tag]) => tag !== 'equal')
      .reduce((total, [, i1, i2, j1, j2]) => total + Math.max(i2 - i1, j2 - j1), 0);
  }

  targetColor(index: number): string {
    return ['#7c3aed', '#0284c7', '#ea580c', '#16a34a', '#db2777', '#ca8a04', '#4f46e5'][index % 7];
  }

  itemLabel(w: ExtractedWindowItem): string {
    return `${w.geneName} - ${w.targetId}`;
  }

  getHeatmapColor(val: number, isDiagonal: boolean): string {
    if (isDiagonal) {
      return 'rgba(46, 204, 113, 0.25)'; // Soft green highlight for diagonal 100%
    }
    if (val >= 90) return `rgba(46, 204, 113, ${0.15 + (val - 90) * 0.015})`;
    if (val >= 75) return `rgba(54, 162, 235, ${0.12 + (val - 75) * 0.01})`;
    if (val >= 50) return `rgba(255, 206, 86, ${0.12 + (val - 50) * 0.008})`;
    return `rgba(240, 242, 245, 0.8)`;
  }

}
