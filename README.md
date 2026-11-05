# CasSago

Sanger chromatogram analysis, sequencing-primer design and gRNA design in a standalone browser application.

- Analyze matched unedited-control and edited-sample chromatograms with submitted amplicon references and the actual sequencing primer.
- Inspect reference-specific mutation groups, target profiles and sequence annotations.
- Design shared or reference-specific primer pairs, and pass selected amplicons into the analysis workspace.
- Find gRNAs on supplied sequences or import current reference sequences from supported libraries.

Sequencing files are processed locally in the browser. Reference-library requests use the selected public reference service. The development access screen is retained; it is a client-side development gate, not server-side authentication.

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
