import {defineConfig} from 'vitest/config';

// Publish the usable primer workflow while the analysis engine is developed
// and validated separately. Keep the broader research test configuration intact.
export default defineConfig({test:{
  maxWorkers:2,
  testTimeout:120000,
  include:[
    'tests/core/sanger-primer-designer.spec.ts',
    'tests/core/sanger-reference-orientation.spec.ts',
    'tests/core/sanger-primer-safety.spec.ts',
    'tests/core/sanger-primer-separation.spec.ts',
    'tests/core/sanger-progressive-primer.spec.ts',
    'tests/core/sanger-design-export.spec.ts',
    'tests/ui/sanger-primer-selection.spec.ts',
    'tests/ui/sanger-design-page.spec.ts',
    'tests/ui/cassago-navigation.spec.ts',
    'tests/ui/reference-library.spec.ts',
  ],
}});
