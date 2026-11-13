# Browser conic optimizer

Reproducible bindings derived from the official [cvxjs source](https://github.com/SteveDiamond/cvxjs/tree/13f0994f153b12c0ac78c9decf7b622627b0dea7/wasm), with Clarabel 0.11.1 and FAER sparse linear algebra. This is a single-thread WebAssembly runtime. It uses neither a server nor a native optimizer in the product.

The binding returns finite full-sized primal and dual witnesses even when Clarabel terminates early. A returned witness or solver status does not establish success: `VerifiedConicSolver` independently checks every original constraint and an original global dual certificate, both with tolerance 1e-7. Equivalent positive objective scaling is undone before these checks.

The vendored `faer-traits` 0.21.5 includes only the attached WASM32 index compatibility patch. It enables the existing i64/u64 index implementations on 32-bit targets, preserving signed sentinel conversion. Live sparse matrix indices must fit the target pointer size.

Build with Rust 1.100.0, the wasm32-unknown-unknown target, and wasm-bindgen-cli 0.2.106:

```sh
cargo build --target wasm32-unknown-unknown --release --locked
wasm-bindgen --target web --out-name clarabel_wasm --out-dir pkg target/wasm32-unknown-unknown/release/clarabel_wasm.wasm
cp pkg/clarabel_wasm.js ../../public/vendor/clarabel-faer-0.11.1/clarabel.mjs
cp pkg/clarabel_wasm_bg.wasm ../../public/vendor/clarabel-faer-0.11.1/
```

Cargo.lock pins all dependencies. The public vendor directory records artifact SHA-256 hashes and dependency licenses. Preserve the original cvxjs Apache license and dependency notices when rebuilding. The older vendor directory remains for archived regression tests.

Node or native development benchmarks are not browser UI timings. Performance and numerical certification must be checked independently in the browser.
