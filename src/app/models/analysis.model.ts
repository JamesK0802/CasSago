export interface AnalysisBreakdown {
  no_indel: number;
  substitution: number;
  in_frame: number;
  out_of_frame: number;
  ambiguous?: number;
  failed?: number;
  failure_reasons?: FailureReasons;
}

export interface FailureReasons {
  fail_no_anchor?: number;
  fail_quality?: number;
  fail_no_coverage?: number;
  fail_no_alignment?: number;
  low_frequency_filtered?: number;
}

export interface AnalysisSummary {
  total_reads: number;
  matched_reads: number;
  aligned_reads: number;
  unmodified: number;
  modified: number;
  editing_efficiency: number;
  
  // Percentages (denominated by matched_reads)
  out_of_frame_pct: number;
  in_frame_pct: number;
  no_indel_pct: number;
  substitution_pct: number;
  substitution_reads?: number;
  indel_editing_efficiency?: number;
  substitution_policy?: string;
  read_ambiguous?: number;
  low_frequency_filtered?: number;
  failed_reads?: number;
  failure_reasons?: FailureReasons;
}

export interface ReadDetail {
  read_index: number;
  target_found: boolean;
  classification: 'out_of_frame' | 'in_frame' | 'substitution' | 'no_indel' | null;
}

export interface AlignmentToken {
  type: 'equal' | 'substitute' | 'delete' | 'insert' | 'unobserved';
  val: string;
}

export interface MutationGroup {
  group_rank: number;
  read_inner: string;
  read_count: number;
  read_pct: number;
  classification: string;
  net_indel: number;
  tokens: AlignmentToken[];
}

export interface TargetResult {
  target_id: string;
  summary: AnalysisSummary;
  breakdown: AnalysisBreakdown;
  read_details: ReadDetail[];
  ref_sequence?: string;
  cut_site_index?: number;
  sgrna_seq?: string;
  display_sgrna_seq?: string;
  grna_start_index?: number;
  is_rc?: boolean;
  strand?: string;
  top_groups?: MutationGroup[];
  /** Groups before the final frequency filter; no alignment needs to be rerun. */
  unfiltered_top_groups?: MutationGroup[];
  threshold_data_complete?: boolean;
  applied_indel_threshold?: number;
}

export interface FileResult {
  fastq_file: string;
  sample_name?: string;
  target_results: TargetResult[];
}

export interface AnalysisResponse {
  metadata: {
    data_type: string;
    phred_threshold: number;
    indel_threshold: number;
    is_multi_reference?: boolean;
    margin_threshold?: number;
  };
  results: FileResult[];
}

// ── Phase 3: Multi-Reference types ───────────────────────────────────────────
export interface GeneResult {
  gene: string;
  homoeolog_group?: string;
  reference_sequence?: string;
  assigned_read_count: number;
  ambiguous_excluded: boolean;
  is_ambiguous_derived?: boolean;
  is_rescued_derived?: boolean;
  analysis_result: {
    targets: TargetResult[];
    multi_target_summary?: MultiTargetSummary;
  };
}

export interface MultiReferenceResponse {
  genes: GeneResult[];
  ambiguous_read_count: number;
  debug?: {
    total_reads_parsed: number;
    phred_passed_count?: number;
    usable_for_assignment_count?: number;
    assignment_filtered_count?: number;
    anchor_matched_count?: number;
    assignment_margin_threshold_used: number;
    genes: Array<{
      gene: string;
      reference_length: number;
      assigned_reads_analyzed: number;
      number_of_targets_analyzed: number;
    }>;
  };
}

// ── Multi-Target Read-Level Types (One-Read-Many-Targets) ────────────────────

export interface TargetEvaluation {
  target_id: string;
  evaluated: boolean;
  eval_fail_reason: string;
  classification: string | null;
  has_substitution: boolean;
  net_indel: number;
  similarity: number;
  edited: boolean;
}

export interface ReadLevelResult {
  read_index: number;
  target_results: TargetEvaluation[];
}

export interface CoEditingEntry {
  target_a: string;
  target_b: string;
  co_edited_reads: number;
}

export interface MultiTargetSummary {
  total_reads: number;
  reads_evaluated: number;
  reads_with_0_edits: number;
  reads_with_1_edit: number;
  reads_with_2_edits: number;
  reads_with_3_plus_edits: number;
  edit_distribution: Record<number, number>;
  co_editing_matrix: CoEditingEntry[];
}

// A file result in multi-ref mode wraps the above
export interface MultiRefFileResult {
  fastq_file: string;
  multi_reference_result: MultiReferenceResponse;
}

// ── Benchmark Tab Types ───────────────────────────────────────────────────────

export interface BenchmarkRow {
  file: File | null;
  r1File?: File | null;
  r2File?: File | null;
  geneName: string;
  targetName: string;
  referenceSequence: string;
  grnaSequence: string;
}

export interface SplitPreviewRow {
  gene: string;
  target: string;
  total: number;
  train_count: number;
  test_count: number;
}

export interface SplitPreview {
  rows: SplitPreviewRow[];
  total: number;
  train_count: number;
  test_count: number;
}

export interface CutSiteInfo {
  gene: string;
  target: string;
  strand: string;
  cut_site: number;
  grna_start: number;
  grna_end: number;
  pam: string;
  pam_found: boolean;
}

export interface BenchmarkClassResult {
  gene: string;
  target: string;
  total: number;
  correct: number;
  wrong: number;
  ambiguous: number;
  correct_rate: number;
  filtered: number;
}

export interface BenchmarkResult {
  subset: string;
  platform?: 'nanopore' | 'illumina';
  preprocessing?: {
    input_molecules: number;
    normalized_molecules: number;
    filtered_molecules: number;
    consensus_molecules: number;
    padded_molecules: number;
  };
  // Raw counts
  total_reads: number;
  filtered_out: number;
  fail_no_anchor: number;
  fail_quality: number;
  fail_similarity: number;
  fail_no_coverage: number;
  fail_no_alignment: number;
  usable_reads: number;
  usable_rate: number;
  // Classification counts (denominator = usable_reads)
  correct_count: number;
  wrong_count: number;
  ambiguous_count: number;
  // Rates over usable reads
  correct_rate: number;
  wrong_rate: number;
  ambiguous_rate: number;
  // Rates over total subset reads
  correct_rate_total: number;
  wrong_rate_total: number;
  ambiguous_rate_total: number;
  split_info: SplitPreviewRow[];
  cut_sites: CutSiteInfo[];
  per_class: BenchmarkClassResult[];
}
