# CasSago: fixed engine validation criteria

These instructions record the user's agreed acceptance criteria (2026-11-15).
They apply to engine development, validation, reports, and handoffs in this tree.
Do not replace them with a different definition of PASS. A change to these
scientific acceptance criteria requires explicit agreement from the user.

## Accuracy and ambiguity: single and multiplex

- Accuracy requirements must not be weaker than those of relevant external tools.
  Compare matched inputs, settings, supported tasks, and metrics. Do not claim
  external-tool parity without the corresponding evidence.
- Recover distinguishable signal accurately: group structure, indel proportions,
  and position-wise base proportions must not be excused by calling them ambiguous.
- Chromatograms may admit different underlying sequence mixtures with the same
  observable signal. Do not force an unsupported unique decomposition.
- Display a representative sequence/group solution, with the supported base
  proportions/ranges and alternatives or ambiguity. The representative need not
  be the literal GT genotype when the data cannot distinguish equivalent mixtures.
- Exact GT sequence/group recovery is an important diagnostic, but is NOT by
  itself the mandatory overall acceptance criterion. Conversely, fitting the trace
  or finding a feasible alternative is NOT by itself acceptance either.
- GT coverage is necessary but not sufficient: audit the reported base and indel
  proportions/ranges against GT, AND their width, informativeness, and calibration.
  Vacuous 0–100% ranges, unjustified wide ranges, and blanket ambiguity do not pass.
- Check joint compatibility as well as marginal coverage. Independently attained
  position-wise extrema must not be presented as a jointly realizable genotype.
- Preserve the distinction between identifiable aggregate quantities and unresolved
  molecular linkage/phase. Do not label an entire task unrecoverable because one
  linkage is ambiguous, or claim full molecular recovery from accurate marginals.
- Apply the SAME principles to single-reference/single-target and multi-reference/
  multi-target inputs. Multiplexing is not permission to weaken acceptance.

## Matched individual-versus-pooled validation: a stricter comparison

The user specifically requires stronger evidence of separation when the SAME
component data are analyzed individually and as a pool.

- Freeze and document component traces, references, targets, mixing weights,
  construction/alignment steps, and primer-design cut-site geometry before testing.
  Maintain the design-required spacing across ALL ref/target cut-site combinations,
  including different genes or different gRNAs. Do not silently add new spacing rules.
- Compare every component/window in the pool with its matched standalone result.
  Normalize consistently: within-reference/window editing fractions are distinct
  from pooled reference abundance; account for known pooling weights explicitly.
- Require agreement in distinguishable group structure, base/indel proportions,
  and supported alternatives. Audit cross-reference/target leakage, lost groups,
  spurious groups, and any widening or loss of information relative to standalone.
- Overlapping broad intervals or each result separately containing GT is NOT
  sufficient evidence that the pool has been separated successfully.
- Where standalone analysis already has genuine ambiguity, compare the supported
  observable quantities and compatible solution sets; do not demand the same
  arbitrary representative or invent unique molecular linkage.
- Individual results are validation comparators, not privileged inputs to the blind
  pooled inference. Keep GT-conditioned feasibility checks separate from recovery.
- If a digitally constructed measured pool fails, investigate construction artifacts
  as well as inference. Do not tune the engine to invalid/chimeric data.

## Evidence and reporting rules

- Freeze acceptance metrics and numeric tolerances BEFORE using a test to claim
  PASS. Recover previously agreed thresholds from existing evidence where possible.
  If a threshold is unspecified, report measurements and mark acceptance undecided;
  do not invent a relaxed cutoff, retroactively widen it, or claim overall PASS.
- Retain original physical constraints, noise assumptions, and certified original
  primal/global optimality tolerance <= 1e-7. Do not widen noise or ranges, prune
  legitimate hypotheses, impose GT-derived allele counts, or use GT warmstarts to
  manufacture acceptance.
- Separate blind inference from posthoc GT-conditioned diagnostics, marginal
  coverage from joint feasibility, and exact representative identity from accuracy
  of identifiable quantities. Report all relevant failures, not just one score.
- Report matched external-tool comparisons separately from self-consistency tests.
  Neither substitutes for the other.
- Report native/multicore development timings, checkpoint/replay timings, and cold
  production-browser timings separately. Numeric completion/QC is not scientific PASS.
- Existing reports using only exact-group recovery or only GT interval coverage must
  be reclassified under these criteria; do not silently declare them all PASS/FAIL.
- Preserve this policy in future handoffs. Do not drift criteria between iterations.

## Development resources (user instruction, 2026-11-15)

- Perform development, diagnosis, optimization, and iterative validation using
  native/multicore resources. Do not use Chrome as the development compute loop.
- Reserve the actual production browser for FINAL verification of the selected
  implementation. Saved browser outputs may be audited offline without rerunning it.
- Use available CPU resources with bounded memory/concurrency; use GPU acceleration
  only where an implemented, numerically validated backend supports the workload.
  Do not claim unused GPU resources were used or equate native timings to web time.
- Iterate on representative cases and the hardest supported case, rather than
  rerunning every validation panel after each exploratory change.

## Keep percentage quantities distinct

- Distinguish (1) reference abundance, (2) total fraction for an indel length
  class, (3) fractions of different molecular sequences with the SAME length,
  and (4) a nucleotide fraction conditional on a group and coordinate. Always
  identify which quantity a number or range describes.
- Model-feasible extrema are not automatically statistical confidence intervals.
  Compare external-tool indel point estimates with CasSago indel point estimates;
  compare uncertainty measures only when their meaning and assumptions match.
  This does not waive accuracy, informative uncertainty, or joint-consistency gates.
- A reduction from a vacuous range to a still broad range is a diagnostic
  improvement, not sufficient evidence of useful quantification or overall PASS.

## Latest user execution order and priorities (2026-11-16)

- Use the same main scientific workflow for single/multi-reference, one/two targets,
  and any number of indel species; no synthetic-specific model or count cutoff.
- Prioritize correct, informative results over speed. Discard unnecessary numerical
  machinery only after proving it preserves results; optimize time after quality.
- Validate in order: single synthetic, single measured, multi synthetic, multi measured.
- Multiplex acceptance specifically requires the SAME component/window profiles
  as matched standalone analyses, including distinguishable group structure,
  proportions, supported alternatives and information content. Fitting or interval
  overlap alone is insufficient. Standalone results must never seed blind pools.
- Use bounded native/multicore development resources, and supported validated GPU
  only if implemented. Actual browser verification comes last for result equivalence
  and final web timing, not as an iterative development compute loop.

## Documented coordinate bug corrections

Preservation of the original domain forbids arbitrary pruning or noise changes
made to pass a test. It does not forbid correcting a demonstrated coordinate bug.
Keep pre-fix failures and input hashes, document the physical correction, freeze
fresh post-fix domains, and certify those actual domains anew at the unchanged
1e-7 tolerance. Never reuse an older domain's proof as a new-model certificate.

## Result-first validation clarification (user, 2026-11-16)

- Signal-prediction QC depends on input quality, available validation geometry,
  calibration and model adequacy. A signal-only regression failure is NOT an
  observed failure of recovered group/base/indel profiles. Report these separately.
- Do not block result-profile comparisons solely because a signal QC or historical
  regression limit failed. Preserve and report the QC failure; do not loosen noise,
  numeric certificates or existing comparison tolerances to erase it.
- Primary outcome evidence: synthetic output versus GT; measured output versus
  matched external results as corroboration (not assumed truth); and multiplex
  output versus matched standalone component/window profiles. Poor QC alone does
  not prove poor input quality or correct/incorrect recovered profiles.


## Latest user reset: point-first common workflow (2026-11-17)

- The current priority is to remove the severe single-sample time bottleneck by returning a minimal common scientific fit first, then improve representative/profile accuracy. Do not keep whole-model worst-case envelope enumeration in basic reporting.
- Default output must include representative DNA and fitted group/base percentages. Mark estimates as estimates; never equate omitted ranges with uniqueness or confidence.
- Range/alternative detail is postprocessing for the selected group only; alternatives remain 10 per page. No single/multi/synthetic/count-specific inference forks.
- Preserve all scientific candidates/physical constraints, postprocess-only reporting threshold, and original numerical certificates. Speed improvements must document which output work was deferred and demonstrate unchanged fitted results.
- Preserve posthoc GT/profile retrieval colouring, including genuinely unresolved alternatives. Independent broad interval inclusion or uncomputed [0,1] placeholders must never receive scientific PASS credit. Distinguishable wrong profiles must be rejected.
- The previous all-blue GT retrieval panels are not default exact sequence recovery or a proof that every difference is indistinguishable. Continue single representative selection and uncertainty calibration work before multiplex.

## Phase-2 evidence and interpretation (2026-11-17)

- Current private v117 intentionally changes representative point selection by a common normalized fluorescence objective and noise-scaled public-base stabilization. Phase-1 point-equality audits apply to deferring output work only; do not falsely require phase-2 points to equal phase-1 or claim they do. The unchanged physical model/candidates and strict 1e-7 original certificates remain required.
- Current source/evidence: validation-runs/single-reestablishment-20261116/photon-proximal-phase-2-summary.json. Aggregate fluorescence/base improvement is much larger than the conditional group-DNA improvement; report both. Five representative-DNA regressions remain failures of literal point recovery. Compatible whole-GT fits do not prove equal-score ambiguity or calibrated uncertainty.
- Independent holdout inference was frozen before generation. Later SNP-marker reporting changes must be disclosed; do not quote pre-change literal-DNA counts as final UI scores. No inference tuning from holdout GT.
- SNP marker partitions preserve every fitted marginal and minor signal. Their cross-site ranges use conservative Frechet envelopes, not demonstrated molecular linkage or attainable simultaneous endpoints. Point100% with shared potential channels must not be called uniquely determined.
- SINGLE remains unfinished; do not resume multiplex based on old all-blue panels or historical running controller identifiers. Native measured smoke completion is not external result-profile validation. No deployment/push authorized by a scientific PASS that has not occurred.

- Exact posthoc raw-ABI alias proofs now exist under validation-runs/single-reestablishment-20261116/photon-proximal-explicit-identical-chromatograms: selected 5-indel +8G and 12-indel +2T/+4A, input173,100% versus0%, identical all original DATA/PBAS/PLOC, no new Sub, class masses preserved. These prove those three ambiguities without a statistical/solver assumption, not all discrepancies. Preserve GT separation: these use GT only as a known forward population AFTER blind output. Do not feed alternatives into inference or count these as new blind success cases.

## Explicit checkpoint authorization (2026-11-17)
User requested Git push before continuing. This authorizes the WIP checkpoint branch, without asserting scientific PASS. Single-first validation remains active.

## User clarification: Sub representative blocks (2026-11-17)
- WT being the largest group does NOT mean maximizing WT by overlapping otherwise separate Sub patterns. Do not revive the rejected maximal-WT display convention.
- User example:adjacent -3..0 at10% and1..3 at20% must be separate representative10%/20% blocks, not a long10% block plus another10% remainder. Same-fraction contiguous sites can share a block; spatially separated or differing-fraction blocks stay separate in the representative where a disjoint partition is feasible. If their total exceeds the family, retain the unresolved mixture rather than discard signal.
- This is a representative-display rule, not measured phasing or an inference constraint. Preserve every parent fitted proportion, all original hypotheses/constraints and unresolved overlap ranges/alternatives. Selected alternatives remain10/page.
- Selected WT is the all-WT subset; a WT marker at one position is not equivalent to no Sub elsewhere. Keep global all-length Sub separate from selected class0 Sub.
