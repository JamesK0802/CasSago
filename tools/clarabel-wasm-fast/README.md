# Verified fast browser queries

These bindings use the official cvxjs binding source recorded in the vendor provenance, Clarabel 0.11.1, and its QDLDL linear solver. Build with Rust 1.100.0 and wasm-bindgen-cli 0.2.106:

```sh
cargo build --lib --target wasm32-unknown-unknown --release --locked
wasm-bindgen --target web --out-name clarabel_wasm --out-dir pkg target/wasm32-unknown-unknown/release/clarabel_wasm.wasm
cp pkg/clarabel_wasm.js ../../public/vendor/clarabel-qdldl-0.11.1/clarabel.mjs
cp pkg/clarabel_wasm_bg.wasm ../../public/vendor/clarabel-qdldl-0.11.1/
```

Cargo.lock pins the build dependencies. Unlike the established FAER build, this module does not need faer-sparse or the faer-traits WASM32 patch. Both builds retain complete finite early-stop witnesses. The public vendor directory contains artifact hashes and dependency notices.

The application uses this optional module for explicit linear queries. Every result must pass the unchanged original primal and global dual checks at1e-7. Failure invokes the established FAER kernel; a timed-out primary is skipped for the remaining queries on that unchanged geometry. Default inference and quadratic representative objectives remain on FAER. If the optional fast module cannot initialize, the established kernel remains usable. No allele, variable, row, cone, or noise allowance is removed by this numerical selection.

This is single-thread browser WASM. Native development and multi-worker endpoint batches are separate experiments; their wall times are not browser latency.
