import {defineConfig} from 'vitest/config';
// Focused checkpoint regressions; not a replacement for scientific acceptance
// or an assertion that the historical release suite is green.
export default defineConfig({test:{maxWorkers:2,include:[
 'tests/core/sanger-representative-selection.spec.ts',
 'tests/core/sanger-representative-origin.spec.ts',
 'tests/core/sanger-point-marker-partition.spec.ts',
 'tests/core/sanger-marker-group-bounds.spec.ts',
 'tests/core/sanger-marginal-representative-coordinate.spec.ts',
 'tests/core/sanger-marginal-pagination.spec.ts',
 'tests/core/sanger-marginal-truth-retrieval.spec.ts',
 'tests/core/sanger-initial-marginal-fit.spec.ts',
 'tests/core/sanger-conic-objective-units.spec.ts',
 'tests/core/sanger-marginal-domain.spec.ts',
 'tests/core/sanger-raw-signal-identifiability.spec.ts',
 'tests/ui/sanger-point-first-display.spec.ts',
 'tests/ui/sanger-global-substitution-display.spec.ts',
 'tests/ui/sanger-ambiguous-base-table.spec.ts',
 'tests/ui/sanger-local-alternatives.spec.ts',
 'tests/ui/sanger-truth-result.spec.ts',
 'tests/ui/sanger-zero-length-display.spec.ts',
 'tests/ui/sanger-proportion-colors.spec.ts',
]}});
