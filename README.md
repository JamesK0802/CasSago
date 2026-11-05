# CasSago

Sanger chromatogram analysis, sequencing-primer design and gRNA design in a standalone browser application.

- Analyze matched unedited-control and edited-sample chromatograms with submitted amplicon references and the actual sequencing primer.
- Inspect reference-specific mutation groups, target profiles and sequence annotations.
- Design shared or reference-specific primer pairs, and pass selected amplicons into the analysis workspace.
- Find gRNAs on supplied sequences or import current reference sequences from supported libraries.

Sequencing files are processed locally in the browser. Reference-library requests use the selected public reference service. Analysis, Results Viewer and Primer Design open directly. gRNA Design retains a client-side development access screen while under construction.

## Run

Use Node.js 22.12+.

```sh
npm ci
npm start
npm test
npm run build
```

`npm test` runs the release suite in `vitest.release.config.ts`. Historical and exploratory tests can be run with `npm run test:all`; they are not all release acceptance criteria.

## GitHub Pages

Pushes to `main` run the release tests, build with `/CasSago/` as the base path and deploy `dist/cassago/browser` through GitHub Actions. Set Settings → Pages → Source to **GitHub Actions**.

Site: https://jamesk0802.github.io/CasSago/

## Decoder checkpoint

The v42 conditional decoder passed five frozen, newly generated synthetic datasets on the same three-reference / two-target panel: 121 of 150 edited groups matched directly, 29 received strict signal-equivalence certificates, and none remained unexplained. These certificates are evaluated after blind inference. This is a fixed-panel synthetic checkpoint, not a claim of general performance on all experimental chromatograms. See `docs/validation-v42.txt` for scope and reproduction.

## Interface

The layout follows CasMango. Analysis sessions retain independent in-memory inputs and results, and navigation is locked during worker activity. The default display window is 90 bp; Advanced provides separate left/right bounds and Cover full ref. Results Viewer reopens CasSago result JSON. Create Curated View permits reversible group exclusion and selection among certified identical-signal representatives; retained fractions are normalized per reference and target. Original data remains available via Back to Original and is embedded in curated JSON exports. PDF and CSV exports follow the active selection. Figure and Guide are placeholders.
