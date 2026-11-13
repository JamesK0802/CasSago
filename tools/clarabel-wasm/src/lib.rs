use wasm_bindgen::prelude::*;
use clarabel::algebra::*;
use clarabel::solver::*;
use serde::{Deserialize, Serialize};
use std::panic;

/// Result of solving an optimization problem
#[derive(Serialize, Deserialize)]
pub struct SolveResult {
    /// Status: "optimal", "infeasible", "unbounded", "max_iterations", "unknown"
    pub status: String,
    /// Optimal objective value (if solved)
    pub obj_val: Option<f64>,
    /// Primal solution vector
    pub x: Option<Vec<f64>>,
    /// Dual solution vector (shadow prices for constraints)
    pub z: Option<Vec<f64>>,
    /// Solve time in seconds
    pub solve_time: f64,
    /// Number of iterations
    pub iterations: u32,
}

/// Cone specification from JavaScript
#[derive(Serialize, Deserialize)]
pub struct ConeSpec {
    /// Zero cone dimension (equality constraints)
    pub zero: usize,
    /// Nonnegative cone dimension (inequality constraints)
    pub nonneg: usize,
    /// Second-order cone dimensions
    pub soc: Vec<usize>,
    /// Number of exponential cones (each is 3D)
    pub exp: usize,
    /// Power cone alphas (each cone is 3D)
    pub power: Vec<f64>,
}

/// Solver settings from JavaScript
#[derive(Serialize, Deserialize, Default)]
pub struct SolverSettings {
    #[serde(default)]
    pub verbose: bool,
    #[serde(default = "default_max_iter")]
    pub max_iter: u32,
    #[serde(default = "default_time_limit")]
    pub time_limit: f64,
    #[serde(default = "default_tol")]
    pub tol_gap_abs: f64,
    #[serde(default = "default_tol")]
    pub tol_gap_rel: f64,
}

fn default_max_iter() -> u32 { 100 }
fn default_time_limit() -> f64 { f64::INFINITY }
fn default_tol() -> f64 { 1e-8 }

/// Build Clarabel cone specification from our ConeSpec
fn build_cones(spec: &ConeSpec) -> Vec<SupportedConeT<f64>> {
    let mut cones = Vec::new();

    if spec.zero > 0 {
        cones.push(ZeroConeT(spec.zero));
    }

    if spec.nonneg > 0 {
        cones.push(NonnegativeConeT(spec.nonneg));
    }

    for &dim in &spec.soc {
        cones.push(SecondOrderConeT(dim));
    }

    for _ in 0..spec.exp {
        cones.push(ExponentialConeT());
    }

    for &alpha in &spec.power {
        cones.push(PowerConeT(alpha));
    }

    cones
}

/// Convert solver status to string
fn status_to_string(status: SolverStatus) -> String {
    match status {
        SolverStatus::Solved => "optimal".to_string(),
        SolverStatus::PrimalInfeasible => "infeasible".to_string(),
        SolverStatus::DualInfeasible => "unbounded".to_string(),
        SolverStatus::MaxIterations => "max_iterations".to_string(),
        SolverStatus::MaxTime => "max_iterations".to_string(),
        _ => "unknown".to_string(),
    }
}

/// Solve a conic optimization problem.
///
/// Problem form:
///   minimize    (1/2) x' P x + q' x
///   subject to  A x + s = b
///               s in K
///
/// # Arguments
/// * `p_col_ptr`, `p_row_idx`, `p_values` - P matrix (upper triangular CSC)
/// * `q` - Linear cost vector
/// * `a_col_ptr`, `a_row_idx`, `a_values` - A matrix (CSC)
/// * `b` - Constraint vector
/// * `n` - Number of variables
/// * `m` - Number of constraints
/// * `cone_spec_json` - Cone specification as JSON string
/// * `settings_json` - Solver settings as JSON string
#[wasm_bindgen]
pub fn solve(
    p_col_ptr: &[u32],
    p_row_idx: &[u32],
    p_values: &[f64],
    q: &[f64],
    a_col_ptr: &[u32],
    a_row_idx: &[u32],
    a_values: &[f64],
    b: &[f64],
    n: u32,
    m: u32,
    cone_spec_json: &str,
    settings_json: &str,
) -> JsValue {
    let n = n as usize;
    let m = m as usize;

    // Convert u32 arrays to usize for Clarabel
    let p_col_ptr: Vec<usize> = p_col_ptr.iter().map(|&x| x as usize).collect();
    let p_row_idx: Vec<usize> = p_row_idx.iter().map(|&x| x as usize).collect();
    let a_col_ptr: Vec<usize> = a_col_ptr.iter().map(|&x| x as usize).collect();
    let a_row_idx: Vec<usize> = a_row_idx.iter().map(|&x| x as usize).collect();

    // Parse cone specification
    let cone_spec: ConeSpec = match serde_json::from_str(cone_spec_json) {
        Ok(spec) => spec,
        Err(e) => {
            let result = SolveResult {
                status: format!("error: invalid cone spec: {}", e),
                obj_val: None,
                x: None,
                z: None,
                solve_time: 0.0,
                iterations: 0,
            };
            return serde_wasm_bindgen::to_value(&result).unwrap();
        }
    };

    // Parse settings
    let settings: SolverSettings = serde_json::from_str(settings_json).unwrap_or_default();

    // Build P matrix (upper triangular)
    let p = CscMatrix::new(
        n,
        n,
        p_col_ptr.clone(),
        p_row_idx.clone(),
        p_values.to_vec(),
    );

    // Build A matrix
    let a = CscMatrix::new(
        m,
        n,
        a_col_ptr.clone(),
        a_row_idx.clone(),
        a_values.to_vec(),
    );

    // Build cones
    let cones = build_cones(&cone_spec);

    // Use a large but finite time_limit instead of infinity (WASM compatibility)
    let time_limit = if settings.time_limit.is_infinite() {
        1e10 // ~317 years, effectively infinite
    } else {
        settings.time_limit
    };

    let mut solver_settings = DefaultSettingsBuilder::default()
        .verbose(settings.verbose)
        .max_iter(settings.max_iter)
        .time_limit(time_limit)
        .tol_gap_abs(settings.tol_gap_abs)
        .tol_gap_rel(settings.tol_gap_rel)
        .build()
        .unwrap();
    solver_settings.max_threads = 1;
    solver_settings.direct_solve_method = "faer".to_string();
    solver_settings.tol_feas = settings.tol_gap_abs.min(1e-8);
    solver_settings.reduced_tol_gap_abs = settings.tol_gap_abs.min(1e-10);
    solver_settings.reduced_tol_gap_rel = settings.tol_gap_rel.min(1e-10);
    solver_settings.reduced_tol_feas = settings.tol_gap_abs.min(1e-10);

    // Use panic::catch_unwind to handle panics gracefully
    let solve_result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
        // Create and solve
        let mut solver = DefaultSolver::new(&p, q, &a, b, &cones, solver_settings).expect("Valid conic dimensions and settings");
        solver.solve();

        // Extract solution
        let solution = &solver.solution;
        let info = &solver.info;

        SolveResult {
            status: status_to_string(solution.status),
            obj_val: if solution.status == SolverStatus::Solved {
                Some(solution.obj_val)
            } else {
                None
            },
            x: if solution.x.len() == n && solution.x.iter().all(|v| v.is_finite()) {
                Some(solution.x.clone())
            } else {
                None
            },
            z: if solution.z.len() == m && solution.z.iter().all(|v| v.is_finite()) {
                Some(solution.z.clone())
            } else {
                None
            },
            solve_time: info.solve_time,
            iterations: info.iterations,
        }
    }));

    let result = match solve_result {
        Ok(result) => result,
        Err(e) => {
            let msg = if let Some(s) = e.downcast_ref::<&str>() {
                s.to_string()
            } else if let Some(s) = e.downcast_ref::<String>() {
                s.clone()
            } else {
                "Unknown panic".to_string()
            };
            SolveResult {
                status: format!("error: solver panic: {}", msg),
                obj_val: None,
                x: None,
                z: None,
                solve_time: 0.0,
                iterations: 0,
            }
        }
    };

    serde_wasm_bindgen::to_value(&result).unwrap()
}

/// Simple test function to verify WASM is working
#[wasm_bindgen]
pub fn test_wasm() -> String {
    "Clarabel WASM is working!".to_string()
}

/// Get version info
#[wasm_bindgen]
pub fn version() -> String {
    "clarabel-wasm 0.1.0".to_string()
}
