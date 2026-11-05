# Deployment

This repository builds a standalone static application with no CasMango runtime or source dependency.

GitHub Pages uses the same build-and-deploy Actions pattern as the existing CasMango site, with `/CasSago/` as the base URL and `dist/cassago/browser` as the artifact path. Enable GitHub Actions as the Pages publishing source before the first deployment.

Local production build: `npm run build`.
Pages production build: `npm run build:pages`.

Hash routes: `#/sanger-design`, `#/grna-design`, `#/contact`.

Maintain the bundled HiGHS loader and WASM together under `public/vendor/highs-1.15.4`. The analysis web worker loads these assets relative to its own deployed URL.
