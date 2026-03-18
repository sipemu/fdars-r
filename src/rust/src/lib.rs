// Allow indexed loops - this is numerical code where indexed loops are clearer
#![allow(clippy::needless_range_loop)]

// Note: Ridge regression is temporarily disabled due to CRAN Rust version constraints
// anofox-regression requires faer 0.23+ which needs Rust 1.84+, but CRAN has Rust 1.81
// Once CRAN updates Rust, re-enable by uncommenting and adding dependencies:
// use anofox_regression::solvers::RidgeRegressor;
// use anofox_regression::{FittedRegressor, Regressor};

use extendr_api::prelude::*;
use nalgebra::{DMatrix, DVector, SVD};
use rand::prelude::*;
use rand_distr::StandardNormal;
use rayon::prelude::*;
use rustfft::{num_complex::Complex, FftPlanner};
use std::f64::consts::PI;

// Import shared utilities from fdars-core
use fdars_core::{hilbert_transform, FdMatrix};
use fdars_core::alignment;
use fdars_core::landmark;
use fdars_core::metric as core_metric;
use fdars_core::tolerance;
use fdars_core::streaming_depth::{self, StreamingDepth, SortedReferenceState, FullReferenceState};

// New modules from fdars-core 0.7.2
use fdars_core::scalar_on_function;
use fdars_core::function_on_scalar;
use fdars_core::classification;
use fdars_core::gmm;
use fdars_core::famm;
use fdars_core::depth as core_depth;

// New modules from fdars-core 0.8.1
use fdars_core::explain;
use fdars_core::explain_generic::FpcPredictor;
use fdars_core::elastic_regression;
use fdars_core::elastic_fpca;
use fdars_core::elastic_changepoint;
use fdars_core::elastic_explain;
use fdars_core::smooth_basis;

// =============================================================================
// Helper functions
// =============================================================================

/// Compute Simpson's rule integration weights for non-uniform grid
fn simpsons_weights(argvals: &[f64]) -> Vec<f64> {
    let n = argvals.len();
    if n < 2 {
        return vec![1.0; n];
    }

    let mut weights = vec![0.0; n];

    if n == 2 {
        // Trapezoidal rule
        let h = argvals[1] - argvals[0];
        weights[0] = h / 2.0;
        weights[1] = h / 2.0;
        return weights;
    }

    // For non-uniform spacing, use composite trapezoidal rule
    for i in 0..n {
        if i == 0 {
            weights[i] = (argvals[1] - argvals[0]) / 2.0;
        } else if i == n - 1 {
            weights[i] = (argvals[n - 1] - argvals[n - 2]) / 2.0;
        } else {
            weights[i] = (argvals[i + 1] - argvals[i - 1]) / 2.0;
        }
    }

    weights
}

/// Compute 2D integration weights using tensor product of 1D weights
/// Returns a flattened vector of weights for an m1 x m2 grid
fn simpsons_weights_2d(argvals_s: &[f64], argvals_t: &[f64]) -> Vec<f64> {
    let weights_s = simpsons_weights(argvals_s);
    let weights_t = simpsons_weights(argvals_t);
    let m1 = argvals_s.len();
    let m2 = argvals_t.len();

    let mut weights = vec![0.0; m1 * m2];
    for i in 0..m1 {
        for j in 0..m2 {
            weights[i * m2 + j] = weights_s[i] * weights_t[j];
        }
    }
    weights
}

// =============================================================================
// Classification probability helper (for conformal prediction)
// =============================================================================

/// Compute class probabilities from a ClassifMethod and score matrix.
/// Replaces the now-private `classif_predict_probs` from fdars-core.
fn classif_probs_from_method(
    method: &classification::ClassifMethod,
    scores: &FdMatrix,
) -> Vec<Vec<f64>> {
    let n = scores.nrows();
    let d = scores.ncols();

    fn softmax(v: &[f64]) -> Vec<f64> {
        let max = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let exps: Vec<f64> = v.iter().map(|&x| (x - max).exp()).collect();
        let sum: f64 = exps.iter().sum();
        exps.iter().map(|&e| e / sum).collect()
    }

    fn mahalanobis_sq(x: &[f64], mu: &[f64], chol: &[f64], d: usize) -> f64 {
        // Solve L * z = (x - mu)
        let diff: Vec<f64> = x.iter().zip(mu).map(|(a, b)| a - b).collect();
        let mut z = diff.clone();
        for i in 0..d {
            for j in 0..i {
                z[i] -= chol[i * d + j] * z[j];
            }
            z[i] /= chol[i * d + i];
        }
        z.iter().map(|v| v * v).sum()
    }

    match method {
        classification::ClassifMethod::Lda { class_means, cov_chol, priors, n_classes } => {
            let g = *n_classes;
            (0..n).map(|i| {
                let x: Vec<f64> = (0..d).map(|j| scores[(i, j)]).collect();
                let disc: Vec<f64> = (0..g).map(|c| {
                    priors[c].max(1e-15).ln() - 0.5 * mahalanobis_sq(&x, &class_means[c], cov_chol, d)
                }).collect();
                softmax(&disc)
            }).collect()
        }
        classification::ClassifMethod::Qda { class_means, class_chols, class_log_dets, priors, n_classes } => {
            let g = *n_classes;
            (0..n).map(|i| {
                let x: Vec<f64> = (0..d).map(|j| scores[(i, j)]).collect();
                let disc: Vec<f64> = (0..g).map(|c| {
                    priors[c].max(1e-15).ln() - 0.5 * (class_log_dets[c] + mahalanobis_sq(&x, &class_means[c], &class_chols[c], d))
                }).collect();
                softmax(&disc)
            }).collect()
        }
        classification::ClassifMethod::Knn { training_scores, training_labels, k, n_classes } => {
            let g = *n_classes;
            let n_train = training_scores.nrows();
            (0..n).map(|i| {
                let x: Vec<f64> = (0..d).map(|j| scores[(i, j)]).collect();
                let mut dists: Vec<(f64, usize)> = (0..n_train).map(|t| {
                    let dist: f64 = (0..d).map(|j| {
                        let diff = x[j] - training_scores[(t, j)];
                        diff * diff
                    }).sum();
                    (dist, training_labels[t])
                }).collect();
                dists.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                let mut counts = vec![0.0; g];
                for &(_, label) in dists.iter().take(*k) {
                    if label < g { counts[label] += 1.0; }
                }
                let sum: f64 = counts.iter().sum();
                if sum > 0.0 { counts.iter().map(|&c| c / sum).collect() }
                else { vec![1.0 / g as f64; g] }
            }).collect()
        }
        _ => {
            // Fallback for future ClassifMethod variants
            (0..n).map(|_| vec![1.0]).collect()
        }
    }
}

// =============================================================================
// fdata functions
// =============================================================================

/// Compute the mean function across all samples (1D)
#[extendr]
fn fdata_mean_1d(data: RMatrix<f64>) -> Robj {
    let nrow = data.nrows();
    let ncol = data.ncols();

    if nrow == 0 || ncol == 0 {
        return Robj::from(Vec::<f64>::new());
    }

    let data_slice = data.as_real_slice().unwrap();

    let mean: Vec<f64> = (0..ncol)
        .into_par_iter()
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..nrow {
                sum += data_slice[i + j * nrow];
            }
            sum / nrow as f64
        })
        .collect();

    Robj::from(mean)
}

/// Compute the mean function across all samples (2D surfaces)
/// Data is stored as n x (m1*m2) matrix where each row is a flattened surface
#[extendr]
fn fdata_mean_2d(data: RMatrix<f64>) -> Robj {
    let nrow = data.nrows();
    let ncol = data.ncols();

    if nrow == 0 || ncol == 0 {
        return Robj::from(Vec::<f64>::new());
    }

    let data_slice = data.as_real_slice().unwrap();

    // Same computation as 1D - just compute pointwise mean
    let mean: Vec<f64> = (0..ncol)
        .into_par_iter()
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..nrow {
                sum += data_slice[i + j * nrow];
            }
            sum / nrow as f64
        })
        .collect();

    Robj::from(mean)
}

/// Center functional data by subtracting the mean function
#[extendr]
fn fdata_center_1d(data: RMatrix<f64>) -> Robj {
    let nrow = data.nrows();
    let ncol = data.ncols();

    if nrow == 0 || ncol == 0 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data_slice = data.as_real_slice().unwrap();

    // First compute the mean for each column (parallelized)
    let means: Vec<f64> = (0..ncol)
        .into_par_iter()
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..nrow {
                sum += data_slice[i + j * nrow];
            }
            sum / nrow as f64
        })
        .collect();

    // Create centered data (parallelized by column)
    let centered: Vec<f64> = (0..ncol)
        .into_par_iter()
        .flat_map(|j| {
            (0..nrow)
                .map(|i| data_slice[i + j * nrow] - means[j])
                .collect::<Vec<_>>()
        })
        .collect();

    // Note: flat_map produces column-major order since we iterate columns first
    // But we need to reorder since flat_map gives [col0_row0, col0_row1, ..., col1_row0, ...]
    let result = RMatrix::new_matrix(nrow, ncol, |i, j| centered[j * nrow + i]);
    r!(result)
}

/// Compute Lp norm for each sample
#[extendr]
fn fdata_norm_lp_1d(data: RMatrix<f64>, argvals: Vec<f64>, p: f64) -> Robj {
    let nrow = data.nrows();
    let ncol = data.ncols();

    if nrow == 0 || ncol == 0 || argvals.len() != ncol {
        return Robj::from(Vec::<f64>::new());
    }

    let data_slice = data.as_real_slice().unwrap();
    let weights = simpsons_weights(&argvals);

    let norms: Vec<f64> = (0..nrow)
        .into_par_iter()
        .map(|i| {
            let mut integral = 0.0;
            for j in 0..ncol {
                let val = data_slice[i + j * nrow].abs().powf(p);
                integral += val * weights[j];
            }
            integral.powf(1.0 / p)
        })
        .collect();

    Robj::from(norms)
}

/// Compute numerical derivative of functional data (parallelized over rows)
#[extendr]
fn fdata_deriv_1d(data: RMatrix<f64>, argvals: Vec<f64>, nderiv: i32) -> Robj {
    let nrow = data.nrows();
    let ncol = data.ncols();

    if nrow == 0 || ncol == 0 || argvals.len() != ncol || nderiv < 1 {
        return r!(RMatrix::new_matrix(nrow, ncol, |_, _| 0.0));
    }

    let data_slice = data.as_real_slice().unwrap();
    let mut current: Vec<f64> = data_slice.to_vec();

    // Pre-compute step sizes for central differences
    let h0 = argvals[1] - argvals[0];
    let hn = argvals[ncol - 1] - argvals[ncol - 2];
    let h_central: Vec<f64> = (1..(ncol - 1))
        .map(|j| argvals[j + 1] - argvals[j - 1])
        .collect();

    for _ in 0..nderiv {
        // Compute derivative for each row in parallel
        let deriv: Vec<f64> = (0..nrow)
            .into_par_iter()
            .flat_map(|i| {
                let mut row_deriv = vec![0.0; ncol];

                // Forward difference at left boundary
                row_deriv[0] = (current[i + nrow] - current[i]) / h0;

                // Central differences for interior points
                for j in 1..(ncol - 1) {
                    row_deriv[j] = (current[i + (j + 1) * nrow] - current[i + (j - 1) * nrow])
                        / h_central[j - 1];
                }

                // Backward difference at right boundary
                row_deriv[ncol - 1] =
                    (current[i + (ncol - 1) * nrow] - current[i + (ncol - 2) * nrow]) / hn;

                row_deriv
            })
            .collect();

        // Reorder from row-major to column-major order
        current = vec![0.0; nrow * ncol];
        for i in 0..nrow {
            for j in 0..ncol {
                current[i + j * nrow] = deriv[i * ncol + j];
            }
        }
    }

    let result = RMatrix::new_matrix(nrow, ncol, |i, j| current[i + j * nrow]);
    r!(result)
}

/// Compute 2D partial derivatives for surface data
///
/// For a surface f(s,t), computes:
/// - ds: partial derivative with respect to s (∂f/∂s)
/// - dt: partial derivative with respect to t (∂f/∂t)
/// - dsdt: mixed partial derivative (∂²f/∂s∂t)
///
/// Data layout: n surfaces, each stored as m1*m2 values in row-major order (s varies fastest)
/// Return: list with three matrices (ds, dt, dsdt), each n x (m1*m2)
#[extendr]
fn fdata_deriv_2d(
    data: RMatrix<f64>,
    argvals_s: Vec<f64>,
    argvals_t: Vec<f64>,
    m1: i32,
    m2: i32,
) -> Robj {
    let n = data.nrows();
    let ncol = data.ncols();
    let m1 = m1 as usize;
    let m2 = m2 as usize;

    if n == 0 || ncol != m1 * m2 || argvals_s.len() != m1 || argvals_t.len() != m2 {
        // Return empty list on invalid input
        return list!(
            ds = r!(RMatrix::<f64>::new_matrix(0, 0, |_, _| 0.0)),
            dt = r!(RMatrix::<f64>::new_matrix(0, 0, |_, _| 0.0)),
            dsdt = r!(RMatrix::<f64>::new_matrix(0, 0, |_, _| 0.0))
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();

    // Helper to access data: surface i, position (si, ti)
    // Data is stored as: data[i + (si + ti * m1) * n]
    let get_val = |i: usize, si: usize, ti: usize| -> f64 { data_slice[i + (si + ti * m1) * n] };

    // Pre-compute step sizes for s direction
    let hs: Vec<f64> = (0..m1)
        .map(|j| {
            if j == 0 {
                argvals_s[1] - argvals_s[0]
            } else if j == m1 - 1 {
                argvals_s[m1 - 1] - argvals_s[m1 - 2]
            } else {
                argvals_s[j + 1] - argvals_s[j - 1]
            }
        })
        .collect();

    // Pre-compute step sizes for t direction
    let ht: Vec<f64> = (0..m2)
        .map(|j| {
            if j == 0 {
                argvals_t[1] - argvals_t[0]
            } else if j == m2 - 1 {
                argvals_t[m2 - 1] - argvals_t[m2 - 2]
            } else {
                argvals_t[j + 1] - argvals_t[j - 1]
            }
        })
        .collect();

    // Compute all derivatives in parallel over surfaces
    let results: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut ds = vec![0.0; m1 * m2];
            let mut dt = vec![0.0; m1 * m2];
            let mut dsdt = vec![0.0; m1 * m2];

            for ti in 0..m2 {
                for si in 0..m1 {
                    let idx = si + ti * m1;

                    // ∂f/∂s using finite differences
                    if si == 0 {
                        // Forward difference
                        ds[idx] = (get_val(i, 1, ti) - get_val(i, 0, ti)) / hs[si];
                    } else if si == m1 - 1 {
                        // Backward difference
                        ds[idx] = (get_val(i, m1 - 1, ti) - get_val(i, m1 - 2, ti)) / hs[si];
                    } else {
                        // Central difference
                        ds[idx] = (get_val(i, si + 1, ti) - get_val(i, si - 1, ti)) / hs[si];
                    }

                    // ∂f/∂t using finite differences
                    if ti == 0 {
                        // Forward difference
                        dt[idx] = (get_val(i, si, 1) - get_val(i, si, 0)) / ht[ti];
                    } else if ti == m2 - 1 {
                        // Backward difference
                        dt[idx] = (get_val(i, si, m2 - 1) - get_val(i, si, m2 - 2)) / ht[ti];
                    } else {
                        // Central difference
                        dt[idx] = (get_val(i, si, ti + 1) - get_val(i, si, ti - 1)) / ht[ti];
                    }

                    // ∂²f/∂s∂t (mixed partial) using finite differences
                    // Use central differences where possible, boundary adjustments otherwise
                    let denom = hs[si] * ht[ti];

                    if si == 0 && ti == 0 {
                        // Forward-forward
                        dsdt[idx] = (get_val(i, 1, 1) - get_val(i, 0, 1) - get_val(i, 1, 0)
                            + get_val(i, 0, 0))
                            / denom;
                    } else if si == m1 - 1 && ti == 0 {
                        // Backward-forward
                        dsdt[idx] =
                            (get_val(i, m1 - 1, 1) - get_val(i, m1 - 2, 1) - get_val(i, m1 - 1, 0)
                                + get_val(i, m1 - 2, 0))
                                / denom;
                    } else if si == 0 && ti == m2 - 1 {
                        // Forward-backward
                        dsdt[idx] =
                            (get_val(i, 1, m2 - 1) - get_val(i, 0, m2 - 1) - get_val(i, 1, m2 - 2)
                                + get_val(i, 0, m2 - 2))
                                / denom;
                    } else if si == m1 - 1 && ti == m2 - 1 {
                        // Backward-backward
                        dsdt[idx] = (get_val(i, m1 - 1, m2 - 1)
                            - get_val(i, m1 - 2, m2 - 1)
                            - get_val(i, m1 - 1, m2 - 2)
                            + get_val(i, m1 - 2, m2 - 2))
                            / denom;
                    } else if si == 0 {
                        // Forward-central
                        dsdt[idx] =
                            (get_val(i, 1, ti + 1) - get_val(i, 0, ti + 1) - get_val(i, 1, ti - 1)
                                + get_val(i, 0, ti - 1))
                                / denom;
                    } else if si == m1 - 1 {
                        // Backward-central
                        dsdt[idx] = (get_val(i, m1 - 1, ti + 1)
                            - get_val(i, m1 - 2, ti + 1)
                            - get_val(i, m1 - 1, ti - 1)
                            + get_val(i, m1 - 2, ti - 1))
                            / denom;
                    } else if ti == 0 {
                        // Central-forward
                        dsdt[idx] =
                            (get_val(i, si + 1, 1) - get_val(i, si - 1, 1) - get_val(i, si + 1, 0)
                                + get_val(i, si - 1, 0))
                                / denom;
                    } else if ti == m2 - 1 {
                        // Central-backward
                        dsdt[idx] = (get_val(i, si + 1, m2 - 1)
                            - get_val(i, si - 1, m2 - 1)
                            - get_val(i, si + 1, m2 - 2)
                            + get_val(i, si - 1, m2 - 2))
                            / denom;
                    } else {
                        // Central-central
                        dsdt[idx] = (get_val(i, si + 1, ti + 1)
                            - get_val(i, si - 1, ti + 1)
                            - get_val(i, si + 1, ti - 1)
                            + get_val(i, si - 1, ti - 1))
                            / denom;
                    }
                }
            }

            (ds, dt, dsdt)
        })
        .collect();

    // Convert to matrices in column-major format
    let ds_mat = RMatrix::new_matrix(n, m1 * m2, |i, j| results[i].0[j]);
    let dt_mat = RMatrix::new_matrix(n, m1 * m2, |i, j| results[i].1[j]);
    let dsdt_mat = RMatrix::new_matrix(n, m1 * m2, |i, j| results[i].2[j]);

    list!(ds = r!(ds_mat), dt = r!(dt_mat), dsdt = r!(dsdt_mat)).into()
}

/// Compute the geometric median (L1 median) of functional data using Weiszfeld's algorithm
/// The geometric median minimizes sum of L2 distances to all curves
#[extendr]
fn geometric_median_1d(data: RMatrix<f64>, argvals: Vec<f64>, max_iter: i32, tol: f64) -> Robj {
    let nrow = data.nrows();
    let ncol = data.ncols();

    if nrow == 0 || ncol == 0 || argvals.len() != ncol {
        return Robj::from(Vec::<f64>::new());
    }

    let data_slice = data.as_real_slice().unwrap();
    let weights = simpsons_weights(&argvals);

    // Initialize with the mean
    let mut median: Vec<f64> = (0..ncol)
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..nrow {
                sum += data_slice[i + j * nrow];
            }
            sum / nrow as f64
        })
        .collect();

    let max_iter = max_iter as usize;

    for _ in 0..max_iter {
        // Compute distances from current median to all curves
        let distances: Vec<f64> = (0..nrow)
            .map(|i| {
                let mut dist_sq = 0.0;
                for j in 0..ncol {
                    let diff = data_slice[i + j * nrow] - median[j];
                    dist_sq += diff * diff * weights[j];
                }
                dist_sq.sqrt()
            })
            .collect();

        // Compute weights (1/distance), handling zero distances
        let eps = 1e-10;
        let inv_distances: Vec<f64> = distances
            .iter()
            .map(|d| if *d > eps { 1.0 / d } else { 1.0 / eps })
            .collect();

        let sum_inv_dist: f64 = inv_distances.iter().sum();

        // Update median using Weiszfeld iteration
        let new_median: Vec<f64> = (0..ncol)
            .map(|j| {
                let mut weighted_sum = 0.0;
                for i in 0..nrow {
                    weighted_sum += data_slice[i + j * nrow] * inv_distances[i];
                }
                weighted_sum / sum_inv_dist
            })
            .collect();

        // Check convergence
        let diff: f64 = median
            .iter()
            .zip(new_median.iter())
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
            / ncol as f64;

        median = new_median;

        if diff < tol {
            break;
        }
    }

    Robj::from(median)
}

/// Compute the geometric median (L1 median) of 2D functional data using Weiszfeld's algorithm
/// Data is stored as n x (m1*m2) matrix where each row is a flattened surface
#[extendr]
fn geometric_median_2d(
    data: RMatrix<f64>,
    argvals_s: Vec<f64>,
    argvals_t: Vec<f64>,
    max_iter: i32,
    tol: f64,
) -> Robj {
    let nrow = data.nrows();
    let ncol = data.ncols();
    let expected_cols = argvals_s.len() * argvals_t.len();

    if nrow == 0 || ncol == 0 || ncol != expected_cols {
        return Robj::from(Vec::<f64>::new());
    }

    let data_slice = data.as_real_slice().unwrap();
    let weights = simpsons_weights_2d(&argvals_s, &argvals_t);

    // Initialize with the mean
    let mut median: Vec<f64> = (0..ncol)
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..nrow {
                sum += data_slice[i + j * nrow];
            }
            sum / nrow as f64
        })
        .collect();

    let max_iter = max_iter as usize;

    for _ in 0..max_iter {
        // Compute distances from current median to all surfaces
        let distances: Vec<f64> = (0..nrow)
            .map(|i| {
                let mut dist_sq = 0.0;
                for j in 0..ncol {
                    let diff = data_slice[i + j * nrow] - median[j];
                    dist_sq += diff * diff * weights[j];
                }
                dist_sq.sqrt()
            })
            .collect();

        // Compute weights (1/distance), handling zero distances
        let eps = 1e-10;
        let inv_distances: Vec<f64> = distances
            .iter()
            .map(|d| if *d > eps { 1.0 / d } else { 1.0 / eps })
            .collect();

        let sum_inv_dist: f64 = inv_distances.iter().sum();

        // Update median using Weiszfeld iteration
        let new_median: Vec<f64> = (0..ncol)
            .map(|j| {
                let mut weighted_sum = 0.0;
                for i in 0..nrow {
                    weighted_sum += data_slice[i + j * nrow] * inv_distances[i];
                }
                weighted_sum / sum_inv_dist
            })
            .collect();

        // Check convergence
        let diff: f64 = median
            .iter()
            .zip(new_median.iter())
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
            / ncol as f64;

        median = new_median;

        if diff < tol {
            break;
        }
    }

    Robj::from(median)
}

// =============================================================================
// depth functions
// =============================================================================

/// Compute Fraiman-Muniz depth
///
/// Uses the FM1 formula from fda.usc: d = 1 - |0.5 - Fn(x)|
/// With scale=TRUE (default): d = (1 - |0.5 - Fn(x)| - 0.5) * 2 = 2 * min(Fn(x), 1-Fn(x))
#[extendr]
fn depth_fm_1d(fdataobj: RMatrix<f64>, fdataori: RMatrix<f64>, _trim: f64, scale: bool) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();

    if ncol_obj != ncol_ori || nobj == 0 || nori == 0 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    // Scale factor: 2 if scale=TRUE to match fda.usc's FM1 with scale=TRUE
    let scale_factor = if scale { 2.0 } else { 1.0 };

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut depth_sum = 0.0;

            for t in 0..n_points {
                let x_t = data_obj[i + t * nobj];
                let mut le_count = 0; // count of y <= x (for Fn(x))

                for j in 0..nori {
                    let y_t = data_ori[j + t * nori];
                    if y_t <= x_t {
                        le_count += 1;
                    }
                }

                // Empirical CDF: Fn(x) = P(Y <= x) = count(Y <= x) / n
                // This matches R's ecdf behavior exactly
                let fn_x = le_count as f64 / nori as f64;
                // FM1 formula: 1 - |0.5 - Fn(x)| which equals min(Fn(x), 1-Fn(x)) + 0.5
                // With scale: 2 * min(Fn(x), 1-Fn(x))
                let univariate_depth = fn_x.min(1.0 - fn_x) * scale_factor;
                depth_sum += univariate_depth;
            }

            depth_sum / n_points as f64
        })
        .collect();

    // Note: trim parameter is used for trimmed mean computation in fda.usc,
    // not for zeroing out depth values. The depth values are always returned as-is.
    // The trim parameter is kept for API compatibility but not used here.
    Robj::from(depths)
}

/// Compute modal depth
#[extendr]
fn depth_mode_1d(fdataobj: RMatrix<f64>, fdataori: RMatrix<f64>, h: f64) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();

    if ncol_obj != ncol_ori || nobj == 0 || nori == 0 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut depth = 0.0;

            for j in 0..nori {
                let mut dist_sq = 0.0;
                for t in 0..n_points {
                    let diff = data_obj[i + t * nobj] - data_ori[j + t * nori];
                    dist_sq += diff * diff;
                }
                let dist = (dist_sq / n_points as f64).sqrt();
                let kernel_val = (-0.5 * (dist / h).powi(2)).exp();
                depth += kernel_val;
            }

            depth / nori as f64
        })
        .collect();

    Robj::from(depths)
}

/// Compute random projection depth
#[extendr]
fn depth_rp_1d(fdataobj: RMatrix<f64>, fdataori: RMatrix<f64>, nproj: i32) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();

    if ncol_obj != ncol_ori || nobj == 0 || nori == 0 || nproj < 1 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let nproj = nproj as usize;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    let mut rng = rand::thread_rng();
    let projections: Vec<Vec<f64>> = (0..nproj)
        .map(|_| {
            let mut proj: Vec<f64> = (0..n_points).map(|_| rng.sample(StandardNormal)).collect();
            let norm: f64 = proj.iter().map(|x| x * x).sum::<f64>().sqrt();
            proj.iter_mut().for_each(|x| *x /= norm);
            proj
        })
        .collect();

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut total_depth = 0.0;

            for proj in &projections {
                let mut proj_i = 0.0;
                for t in 0..n_points {
                    proj_i += data_obj[i + t * nobj] * proj[t];
                }

                let mut proj_ori: Vec<f64> = (0..nori)
                    .map(|j| {
                        let mut p = 0.0;
                        for t in 0..n_points {
                            p += data_ori[j + t * nori] * proj[t];
                        }
                        p
                    })
                    .collect();

                proj_ori.sort_by(|a, b| a.partial_cmp(b).unwrap());

                let below = proj_ori.iter().filter(|&&x| x < proj_i).count();
                let above = proj_ori.iter().filter(|&&x| x > proj_i).count();
                let depth = (below.min(above) as f64 + 1.0) / (nori as f64 + 1.0);

                total_depth += depth;
            }

            total_depth / nproj as f64
        })
        .collect();

    Robj::from(depths)
}

/// Compute random Tukey depth
#[extendr]
fn depth_rt_1d(fdataobj: RMatrix<f64>, fdataori: RMatrix<f64>, nproj: i32) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();

    if ncol_obj != ncol_ori || nobj == 0 || nori == 0 || nproj < 1 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let nproj = nproj as usize;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    let mut rng = rand::thread_rng();
    let projections: Vec<Vec<f64>> = (0..nproj)
        .map(|_| {
            let mut proj: Vec<f64> = (0..n_points).map(|_| rng.sample(StandardNormal)).collect();
            let norm: f64 = proj.iter().map(|x| x * x).sum::<f64>().sqrt();
            proj.iter_mut().for_each(|x| *x /= norm);
            proj
        })
        .collect();

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut min_depth = f64::INFINITY;

            for proj in &projections {
                let mut proj_i = 0.0;
                for t in 0..n_points {
                    proj_i += data_obj[i + t * nobj] * proj[t];
                }

                let proj_ori: Vec<f64> = (0..nori)
                    .map(|j| {
                        let mut p = 0.0;
                        for t in 0..n_points {
                            p += data_ori[j + t * nori] * proj[t];
                        }
                        p
                    })
                    .collect();

                let below = proj_ori.iter().filter(|&&x| x < proj_i).count();
                let above = proj_ori.iter().filter(|&&x| x > proj_i).count();
                let depth = (below.min(above) as f64 + 1.0) / (nori as f64 + 1.0);

                min_depth = min_depth.min(depth);
            }

            min_depth
        })
        .collect();

    Robj::from(depths)
}

/// Compute Functional Spatial Depth
#[extendr]
fn depth_fsd_1d(fdataobj: RMatrix<f64>, fdataori: RMatrix<f64>) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();

    if ncol_obj != ncol_ori || nobj == 0 || nori == 0 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut sum_unit = vec![0.0; n_points];

            for j in 0..nori {
                let mut direction = vec![0.0; n_points];
                let mut norm_sq = 0.0;

                for t in 0..n_points {
                    direction[t] = data_ori[j + t * nori] - data_obj[i + t * nobj];
                    norm_sq += direction[t] * direction[t];
                }

                let norm = norm_sq.sqrt();
                if norm > 1e-10 {
                    for t in 0..n_points {
                        sum_unit[t] += direction[t] / norm;
                    }
                }
            }

            let mut avg_norm_sq = 0.0;
            for t in 0..n_points {
                let avg = sum_unit[t] / nori as f64;
                avg_norm_sq += avg * avg;
            }

            1.0 - avg_norm_sq.sqrt()
        })
        .collect();

    Robj::from(depths)
}

/// Kernel Functional Spatial Depth (KFSD) for 1D functional data
/// Implements the RKHS-based formulation matching fda.usc
/// h is treated as the actual bandwidth, matching how fda.usc uses hq2
/// argvals is used for trapezoidal integration to compute L2 norms
#[extendr]
fn depth_kfsd_1d(fdataobj: RMatrix<f64>, fdataori: RMatrix<f64>, argvals: Robj, h: f64) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();

    if ncol_obj != ncol_ori || nobj == 0 || nori == 0 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    // Get argvals for trapezoidal integration
    let argvals_vec: Vec<f64> = argvals
        .as_real_vector()
        .unwrap_or_else(|| (0..n_points).map(|i| i as f64).collect());

    // Compute trapezoidal integration weights
    let weights = simpsons_weights(&argvals_vec);

    let h_sq = h * h;

    // Helper function to compute integrated L2 norm squared between two curves
    // using trapezoidal rule: ||f - g||^2 = integral((f(t) - g(t))^2 dt)
    let norm_sq_integrated = |data1: &[f64],
                              row1: usize,
                              nrow1: usize,
                              data2: &[f64],
                              row2: usize,
                              nrow2: usize|
     -> f64 {
        let mut sum = 0.0;
        for t in 0..n_points {
            let diff = data1[row1 + t * nrow1] - data2[row2 + t * nrow2];
            sum += weights[t] * diff * diff;
        }
        sum
    };

    // Gaussian kernel: K(x,y) = exp(-||x-y||^2 / h^2)
    let kern = |dist_sq: f64| -> f64 { (-dist_sq / h_sq).exp() };

    // Pre-compute M1[j,k] = K(X_j, X_k) for all pairs in reference data
    // This is symmetric, so we only compute upper triangle (parallelized)
    let m1_upper: Vec<(usize, usize, f64)> = (0..nori)
        .into_par_iter()
        .flat_map(|j| {
            ((j + 1)..nori)
                .map(|k| {
                    let mut sum = 0.0;
                    for t in 0..n_points {
                        let diff = data_ori[j + t * nori] - data_ori[k + t * nori];
                        sum += weights[t] * diff * diff;
                    }
                    let kval = (-sum / h_sq).exp();
                    (j, k, kval)
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let mut m1 = vec![vec![0.0; nori]; nori];
    for j in 0..nori {
        m1[j][j] = 1.0; // K(X_j, X_j) = 1
    }
    for (j, k, kval) in m1_upper {
        m1[j][k] = kval;
        m1[k][j] = kval;
    }

    let nori_f64 = nori as f64;

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            // K02[i] = K(x_i, x_i) = 1 for Gaussian kernel
            let k02_i = 1.0;

            // Compute M2[i,j] = K(x_i, X_j) for all j
            let m2: Vec<f64> = (0..nori)
                .map(|j| {
                    let d_sq = norm_sq_integrated(data_obj, i, nobj, data_ori, j, nori);
                    kern(d_sq)
                })
                .collect();

            // Compute sum of M[i,j,k] over all j,k
            // M[i,j,k] = (K02_i + M1[j,k] - M2[i,j] - M2[i,k]) /
            //            (sqrt(K02_i + K01_j - 2*M2[i,j]) * sqrt(K02_i + K01_k - 2*M2[i,k]))
            // where K01_j = K(X_j, X_j) = 1
            let mut total_sum = 0.0;
            let mut valid_count = 0;

            for j in 0..nori {
                let k01_j = 1.0;
                let denom_j_sq = k02_i + k01_j - 2.0 * m2[j];

                // Skip if denominator would be zero or negative (identical curves)
                if denom_j_sq < 1e-20 {
                    continue;
                }
                let denom_j = denom_j_sq.sqrt();

                for k in 0..nori {
                    let k01_k = 1.0;
                    let denom_k_sq = k02_i + k01_k - 2.0 * m2[k];

                    if denom_k_sq < 1e-20 {
                        continue;
                    }
                    let denom_k = denom_k_sq.sqrt();

                    let numerator = k02_i + m1[j][k] - m2[j] - m2[k];
                    let denom = denom_j * denom_k;

                    if denom > 1e-20 {
                        let m_ijk = numerator / denom;
                        if m_ijk.is_finite() {
                            total_sum += m_ijk;
                            valid_count += 1;
                        }
                    }
                }
            }

            // depth = 1 - sqrt(sum(M)) / n
            // Note: fda.usc doesn't account for NA values in the same way
            if valid_count > 0 && total_sum >= 0.0 {
                1.0 - total_sum.sqrt() / nori_f64
            } else if total_sum < 0.0 {
                // Handle numerical issues - clamp to valid range
                1.0
            } else {
                0.0
            }
        })
        .collect();

    Robj::from(depths)
}

/// Kernel Functional Spatial Depth (KFSD) for 2D functional data
/// Implements the RKHS-based formulation matching fda.usc
#[extendr]
fn depth_kfsd_2d(
    fdataobj: RMatrix<f64>,
    fdataori: RMatrix<f64>,
    _m1: i32,
    _m2: i32,
    h: f64,
) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();

    if ncol_obj != ncol_ori || nobj == 0 || nori == 0 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    let h_sq = h * h;

    // Helper function to compute L2 norm squared between two surfaces
    let norm_sq = |data1: &[f64],
                   row1: usize,
                   nrow1: usize,
                   data2: &[f64],
                   row2: usize,
                   nrow2: usize|
     -> f64 {
        let mut sum = 0.0;
        for t in 0..n_points {
            let diff = data1[row1 + t * nrow1] - data2[row2 + t * nrow2];
            sum += diff * diff;
        }
        sum
    };

    // Gaussian kernel: K(x,y) = exp(-||x-y||^2 / h^2)
    let kern = |dist_sq: f64| -> f64 { (-dist_sq / h_sq).exp() };

    // Pre-compute M1[j,k] = K(X_j, X_k) for all pairs in reference data (parallelized)
    let m1_upper: Vec<(usize, usize, f64)> = (0..nori)
        .into_par_iter()
        .flat_map(|j| {
            ((j + 1)..nori)
                .map(|k| {
                    let d_sq = norm_sq(data_ori, j, nori, data_ori, k, nori);
                    let kval = kern(d_sq);
                    (j, k, kval)
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let mut m1_mat = vec![vec![0.0; nori]; nori];
    for j in 0..nori {
        m1_mat[j][j] = 1.0;
    }
    for (j, k, kval) in m1_upper {
        m1_mat[j][k] = kval;
        m1_mat[k][j] = kval;
    }

    let nori_f64 = nori as f64;

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let k02_i = 1.0;

            // Compute M2[i,j] = K(x_i, X_j) for all j
            let m2: Vec<f64> = (0..nori)
                .map(|j| {
                    let d_sq = norm_sq(data_obj, i, nobj, data_ori, j, nori);
                    kern(d_sq)
                })
                .collect();

            let mut total_sum = 0.0;
            let mut valid_count = 0;

            for j in 0..nori {
                let k01_j = 1.0;
                let denom_j_sq = k02_i + k01_j - 2.0 * m2[j];

                if denom_j_sq < 1e-20 {
                    continue;
                }
                let denom_j = denom_j_sq.sqrt();

                for k in 0..nori {
                    let k01_k = 1.0;
                    let denom_k_sq = k02_i + k01_k - 2.0 * m2[k];

                    if denom_k_sq < 1e-20 {
                        continue;
                    }
                    let denom_k = denom_k_sq.sqrt();

                    let numerator = k02_i + m1_mat[j][k] - m2[j] - m2[k];
                    let denom = denom_j * denom_k;

                    if denom > 1e-20 {
                        let m_ijk = numerator / denom;
                        if m_ijk.is_finite() {
                            total_sum += m_ijk;
                            valid_count += 1;
                        }
                    }
                }
            }

            if valid_count > 0 && total_sum >= 0.0 {
                1.0 - total_sum.sqrt() / nori_f64
            } else if total_sum < 0.0 {
                1.0
            } else {
                0.0
            }
        })
        .collect();

    Robj::from(depths)
}

/// Band Depth (BD) for 1D functional data
/// BD(x) = proportion of pairs (i,j) where x lies within the band formed by curves i and j
/// A curve lies in the band if at every time point t, min(X_i(t), X_j(t)) <= x(t) <= max(X_i(t), X_j(t))
#[extendr]
fn depth_bd_1d(fdataobj: RMatrix<f64>, fdataori: RMatrix<f64>) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();

    if ncol_obj != ncol_ori || nobj == 0 || nori < 2 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    // Number of pairs C(n, 2) = n * (n-1) / 2
    let n_pairs = (nori * (nori - 1)) / 2;

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut count_in_band = 0usize;

            // For each pair (j, k) with j < k
            for j in 0..nori {
                for k in (j + 1)..nori {
                    // Check if curve i is completely inside the band formed by j and k
                    let mut inside_band = true;

                    for t in 0..n_points {
                        let x_t = data_obj[i + t * nobj];
                        let y_j_t = data_ori[j + t * nori];
                        let y_k_t = data_ori[k + t * nori];

                        let band_min = y_j_t.min(y_k_t);
                        let band_max = y_j_t.max(y_k_t);

                        if x_t < band_min || x_t > band_max {
                            inside_band = false;
                            break;
                        }
                    }

                    if inside_band {
                        count_in_band += 1;
                    }
                }
            }

            count_in_band as f64 / n_pairs as f64
        })
        .collect();

    Robj::from(depths)
}

/// Modified Band Depth (MBD) for 1D functional data
/// MBD(x) = average over pairs (i,j) of the proportion of the domain where x is inside the band
/// This is more robust than BD as it doesn't require complete containment
#[extendr]
fn depth_mbd_1d(fdataobj: RMatrix<f64>, fdataori: RMatrix<f64>) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();

    if ncol_obj != ncol_ori || nobj == 0 || nori < 2 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    // Number of pairs C(n, 2) = n * (n-1) / 2
    let n_pairs = (nori * (nori - 1)) / 2;

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut total_proportion = 0.0;

            // For each pair (j, k) with j < k
            for j in 0..nori {
                for k in (j + 1)..nori {
                    // Count the proportion of time points where curve i is inside the band
                    let mut count_inside = 0usize;

                    for t in 0..n_points {
                        let x_t = data_obj[i + t * nobj];
                        let y_j_t = data_ori[j + t * nori];
                        let y_k_t = data_ori[k + t * nori];

                        let band_min = y_j_t.min(y_k_t);
                        let band_max = y_j_t.max(y_k_t);

                        if x_t >= band_min && x_t <= band_max {
                            count_inside += 1;
                        }
                    }

                    total_proportion += count_inside as f64 / n_points as f64;
                }
            }

            total_proportion / n_pairs as f64
        })
        .collect();

    Robj::from(depths)
}

/// Modified Epigraph Index (MEI) for 1D functional data
/// MEI measures the proportion of time a curve is below other curves
/// MEI(x_i) = (1/n) * sum_j (1/m) * sum_t I(x_i(t) < x_j(t)) + 0.5*I(x_i(t) = x_j(t))
#[extendr]
fn depth_mei_1d(fdataobj: RMatrix<f64>, fdataori: RMatrix<f64>) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();

    if ncol_obj != ncol_ori || nobj == 0 || nori == 0 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    let mei: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut total = 0.0;

            for j in 0..nori {
                let mut count = 0.0;

                for t in 0..n_points {
                    let xi = data_obj[i + t * nobj];
                    let xj = data_ori[j + t * nori];

                    if xi < xj {
                        count += 1.0;
                    } else if (xi - xj).abs() < 1e-12 {
                        // Handle equality with small tolerance
                        count += 0.5;
                    }
                }

                total += count / n_points as f64;
            }

            total / nori as f64
        })
        .collect();

    Robj::from(mei)
}

// =============================================================================
// 2D depth functions (for surfaces)
// =============================================================================

/// Fraiman-Muniz depth for 2D functional data (surfaces)
/// Integrates univariate depth over (s,t) grid
#[extendr]
fn depth_fm_2d(
    fdataobj: RMatrix<f64>,
    fdataori: RMatrix<f64>,
    m1: i32,
    m2: i32,
    scale: bool,
) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();
    let m1 = m1 as usize;
    let m2 = m2 as usize;

    if ncol_obj != ncol_ori || ncol_obj != m1 * m2 || nobj == 0 || nori == 0 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    let scale_factor = if scale { 2.0 } else { 1.0 };

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut depth_sum = 0.0;

            // Iterate over all (s,t) grid points
            for idx in 0..n_points {
                let x_st = data_obj[i + idx * nobj];
                let mut le_count = 0;

                for j in 0..nori {
                    let y_st = data_ori[j + idx * nori];
                    if y_st <= x_st {
                        le_count += 1;
                    }
                }

                let fn_x = le_count as f64 / nori as f64;
                let univariate_depth = fn_x.min(1.0 - fn_x) * scale_factor;
                depth_sum += univariate_depth;
            }

            depth_sum / n_points as f64
        })
        .collect();

    Robj::from(depths)
}

/// Modal depth for 2D functional data (surfaces)
/// Uses L2 distance in the flattened surface space
#[extendr]
fn depth_mode_2d(fdataobj: RMatrix<f64>, fdataori: RMatrix<f64>, m1: i32, m2: i32, h: f64) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();
    let m1 = m1 as usize;
    let m2 = m2 as usize;

    if ncol_obj != ncol_ori || ncol_obj != m1 * m2 || nobj == 0 || nori == 0 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut depth = 0.0;

            for j in 0..nori {
                let mut dist_sq = 0.0;
                for idx in 0..n_points {
                    let diff = data_obj[i + idx * nobj] - data_ori[j + idx * nori];
                    dist_sq += diff * diff;
                }
                // Normalize by number of points like in 1D
                let dist = (dist_sq / n_points as f64).sqrt();
                let kernel_val = (-0.5 * (dist / h).powi(2)).exp();
                depth += kernel_val;
            }

            depth / nori as f64
        })
        .collect();

    Robj::from(depths)
}

/// Random projection depth for 2D functional data (surfaces)
/// Projects surfaces to scalars using random projections
#[extendr]
fn depth_rp_2d(
    fdataobj: RMatrix<f64>,
    fdataori: RMatrix<f64>,
    m1: i32,
    m2: i32,
    nproj: i32,
) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();
    let m1 = m1 as usize;
    let m2 = m2 as usize;

    if ncol_obj != ncol_ori || ncol_obj != m1 * m2 || nobj == 0 || nori == 0 || nproj < 1 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let nproj = nproj as usize;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    let mut rng = rand::thread_rng();
    let projections: Vec<Vec<f64>> = (0..nproj)
        .map(|_| {
            let mut proj: Vec<f64> = (0..n_points).map(|_| rng.sample(StandardNormal)).collect();
            let norm: f64 = proj.iter().map(|x| x * x).sum::<f64>().sqrt();
            proj.iter_mut().for_each(|x| *x /= norm);
            proj
        })
        .collect();

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut total_depth = 0.0;

            for proj in &projections {
                let mut proj_i = 0.0;
                for idx in 0..n_points {
                    proj_i += data_obj[i + idx * nobj] * proj[idx];
                }

                let mut proj_ori: Vec<f64> = (0..nori)
                    .map(|j| {
                        let mut p = 0.0;
                        for idx in 0..n_points {
                            p += data_ori[j + idx * nori] * proj[idx];
                        }
                        p
                    })
                    .collect();

                proj_ori.sort_by(|a, b| a.partial_cmp(b).unwrap());

                let mut le_count = 0;
                for &val in &proj_ori {
                    if val <= proj_i {
                        le_count += 1;
                    }
                }

                let fn_x = le_count as f64 / nori as f64;
                let depth = fn_x.min(1.0 - fn_x) * 2.0;
                total_depth += depth;
            }

            total_depth / nproj as f64
        })
        .collect();

    Robj::from(depths)
}

/// Random Tukey depth for 2D functional data (surfaces)
#[extendr]
fn depth_rt_2d(
    fdataobj: RMatrix<f64>,
    fdataori: RMatrix<f64>,
    m1: i32,
    m2: i32,
    nproj: i32,
) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();
    let m1 = m1 as usize;
    let m2 = m2 as usize;

    if ncol_obj != ncol_ori || ncol_obj != m1 * m2 || nobj == 0 || nori == 0 || nproj < 1 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let nproj = nproj as usize;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    let mut rng = rand::thread_rng();
    let projections: Vec<Vec<f64>> = (0..nproj)
        .map(|_| {
            let mut proj: Vec<f64> = (0..n_points).map(|_| rng.sample(StandardNormal)).collect();
            let norm: f64 = proj.iter().map(|x| x * x).sum::<f64>().sqrt();
            proj.iter_mut().for_each(|x| *x /= norm);
            proj
        })
        .collect();

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut min_depth = f64::INFINITY;

            for proj in &projections {
                let mut proj_i = 0.0;
                for idx in 0..n_points {
                    proj_i += data_obj[i + idx * nobj] * proj[idx];
                }

                let proj_ori: Vec<f64> = (0..nori)
                    .map(|j| {
                        let mut p = 0.0;
                        for idx in 0..n_points {
                            p += data_ori[j + idx * nori] * proj[idx];
                        }
                        p
                    })
                    .collect();

                let below = proj_ori.iter().filter(|&&x| x < proj_i).count();
                let above = proj_ori.iter().filter(|&&x| x > proj_i).count();
                let depth = (below.min(above) as f64 + 1.0) / (nori as f64 + 1.0);

                min_depth = min_depth.min(depth);
            }

            min_depth
        })
        .collect();

    Robj::from(depths)
}

/// Functional Spatial Depth for 2D functional data
#[extendr]
fn depth_fsd_2d(fdataobj: RMatrix<f64>, fdataori: RMatrix<f64>, m1: i32, m2: i32) -> Robj {
    let nobj = fdataobj.nrows();
    let nori = fdataori.nrows();
    let ncol_obj = fdataobj.ncols();
    let ncol_ori = fdataori.ncols();
    let m1 = m1 as usize;
    let m2 = m2 as usize;

    if ncol_obj != ncol_ori || ncol_obj != m1 * m2 || nobj == 0 || nori == 0 {
        return Robj::from(Vec::<f64>::new());
    }

    let n_points = ncol_obj;
    let data_obj = fdataobj.as_real_slice().unwrap();
    let data_ori = fdataori.as_real_slice().unwrap();

    let depths: Vec<f64> = (0..nobj)
        .into_par_iter()
        .map(|i| {
            let mut sum_unit = vec![0.0; n_points];

            for j in 0..nori {
                let mut direction = vec![0.0; n_points];
                let mut norm_sq = 0.0;

                for idx in 0..n_points {
                    direction[idx] = data_ori[j + idx * nori] - data_obj[i + idx * nobj];
                    norm_sq += direction[idx] * direction[idx];
                }

                let norm = norm_sq.sqrt();
                if norm > 1e-10 {
                    for idx in 0..n_points {
                        sum_unit[idx] += direction[idx] / norm;
                    }
                }
            }

            let mut avg_norm_sq = 0.0;
            for idx in 0..n_points {
                let avg = sum_unit[idx] / nori as f64;
                avg_norm_sq += avg * avg;
            }

            1.0 - avg_norm_sq.sqrt()
        })
        .collect();

    Robj::from(depths)
}

// =============================================================================
// metric functions
// =============================================================================

/// Compute Lp distance matrix between two sets of functional data
#[extendr]
fn metric_lp_1d(
    fdata1: RMatrix<f64>,
    fdata2: RMatrix<f64>,
    argvals: Vec<f64>,
    p: f64,
    w: Vec<f64>,
) -> Robj {
    let n1 = fdata1.nrows();
    let n2 = fdata2.nrows();
    let ncol1 = fdata1.ncols();
    let ncol2 = fdata2.ncols();

    if ncol1 != ncol2 || ncol1 != argvals.len() {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let n_points = ncol1;
    let data1 = fdata1.as_real_slice().unwrap();
    let data2 = fdata2.as_real_slice().unwrap();

    let base_weights = simpsons_weights(&argvals);
    let weights: Vec<f64> = if w.len() == n_points {
        base_weights
            .iter()
            .zip(w.iter())
            .map(|(b, u)| b * u)
            .collect()
    } else {
        base_weights
    };

    let distances: Vec<f64> = (0..n2)
        .into_par_iter()
        .flat_map(|j| {
            (0..n1)
                .map(|i| {
                    let mut integral = 0.0;
                    for k in 0..n_points {
                        let diff = (data1[i + k * n1] - data2[j + k * n2]).abs();
                        integral += diff.powf(p) * weights[k];
                    }
                    integral.powf(1.0 / p)
                })
                .collect::<Vec<f64>>()
        })
        .collect();

    let result = RMatrix::new_matrix(n1, n2, |i, j| distances[i + j * n1]);
    r!(result)
}

/// Compute Lp distance matrix for self-distances (symmetric)
#[extendr]
fn metric_lp_self_1d(fdata: RMatrix<f64>, argvals: Vec<f64>, p: f64, w: Vec<f64>) -> Robj {
    let n = fdata.nrows();
    let ncol = fdata.ncols();

    if ncol != argvals.len() {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let n_points = ncol;
    let data = fdata.as_real_slice().unwrap();

    let base_weights = simpsons_weights(&argvals);
    let weights: Vec<f64> = if w.len() == n_points {
        base_weights
            .iter()
            .zip(w.iter())
            .map(|(b, u)| b * u)
            .collect()
    } else {
        base_weights
    };

    let mut distances = vec![0.0; n * n];

    let upper_triangle: Vec<(usize, usize, f64)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            ((i + 1)..n)
                .map(|j| {
                    let mut integral = 0.0;
                    for k in 0..n_points {
                        let diff = (data[i + k * n] - data[j + k * n]).abs();
                        integral += diff.powf(p) * weights[k];
                    }
                    (i, j, integral.powf(1.0 / p))
                })
                .collect::<Vec<_>>()
        })
        .collect();

    for (i, j, dist) in upper_triangle {
        distances[i + j * n] = dist;
        distances[j + i * n] = dist;
    }

    let result = RMatrix::new_matrix(n, n, |i, j| distances[i + j * n]);
    r!(result)
}

/// Compute Hausdorff distance matrix for self-distances (symmetric)
///
/// The Hausdorff distance treats curves as sets of points (t, f(t)) in 2D space.
/// For each pair of curves, computes: max(max_s min_t d(s,t), max_t min_s d(s,t))
/// where d(s,t) = sqrt((x(s) - y(t))^2 + (s - t)^2)
#[extendr]
fn metric_hausdorff_1d(fdata: RMatrix<f64>, argvals: Vec<f64>) -> Robj {
    let n = fdata.nrows();
    let m = fdata.ncols();

    if m != argvals.len() || n == 0 || m == 0 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data = fdata.as_real_slice().unwrap();

    // Precompute squared time differences: Mtt[s][t] = (argvals[s] - argvals[t])^2
    let mtt: Vec<f64> = {
        let mut result = Vec::with_capacity(m * m);
        for s in 0..m {
            for t in 0..m {
                let diff = argvals[s] - argvals[t];
                result.push(diff * diff);
            }
        }
        result
    };

    // Compute upper triangle of distance matrix in parallel
    let upper_triangle: Vec<(usize, usize, f64)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            ((i + 1)..n)
                .map(|j| {
                    // Compute Hausdorff distance between curve i and curve j
                    // Mxy[s][t] = sqrt((x[s] - y[t])^2 + (s - t)^2)

                    // For each row s, find min over columns t
                    let max_row_min = (0..m)
                        .map(|s| {
                            let x_s = data[i + s * n];
                            (0..m)
                                .map(|t| {
                                    let y_t = data[j + t * n];
                                    let val_diff = x_s - y_t;
                                    (val_diff * val_diff + mtt[s * m + t]).sqrt()
                                })
                                .fold(f64::INFINITY, |a, b| a.min(b))
                        })
                        .fold(f64::NEG_INFINITY, |a, b| a.max(b));

                    // For each column t, find min over rows s
                    let max_col_min = (0..m)
                        .map(|t| {
                            let y_t = data[j + t * n];
                            (0..m)
                                .map(|s| {
                                    let x_s = data[i + s * n];
                                    let val_diff = x_s - y_t;
                                    (val_diff * val_diff + mtt[s * m + t]).sqrt()
                                })
                                .fold(f64::INFINITY, |a, b| a.min(b))
                        })
                        .fold(f64::NEG_INFINITY, |a, b| a.max(b));

                    let hausdorff_dist = max_row_min.max(max_col_min);
                    (i, j, hausdorff_dist)
                })
                .collect::<Vec<_>>()
        })
        .collect();

    // Build symmetric distance matrix
    let mut distances = vec![0.0; n * n];
    for (i, j, dist) in upper_triangle {
        distances[i + j * n] = dist;
        distances[j + i * n] = dist;
    }

    let result = RMatrix::new_matrix(n, n, |i, j| distances[i + j * n]);
    r!(result)
}

/// Compute Hausdorff distance matrix for cross-distances (n1 x n2)
///
/// The Hausdorff distance treats curves as sets of points (t, f(t)) in 2D space.
#[extendr]
fn metric_hausdorff_cross_1d(
    fdata1: RMatrix<f64>,
    fdata2: RMatrix<f64>,
    argvals: Vec<f64>,
) -> Robj {
    let n1 = fdata1.nrows();
    let n2 = fdata2.nrows();
    let m = fdata1.ncols();

    if m != argvals.len() || m != fdata2.ncols() || n1 == 0 || n2 == 0 || m == 0 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data1 = fdata1.as_real_slice().unwrap();
    let data2 = fdata2.as_real_slice().unwrap();

    // Precompute squared time differences
    let mtt: Vec<f64> = {
        let mut result = Vec::with_capacity(m * m);
        for s in 0..m {
            for t in 0..m {
                let diff = argvals[s] - argvals[t];
                result.push(diff * diff);
            }
        }
        result
    };

    // Compute all n1 x n2 distances in parallel
    let distances: Vec<f64> = (0..n1)
        .into_par_iter()
        .flat_map(|i| {
            (0..n2)
                .map(|j| {
                    // For each row s, find min over columns t
                    let max_row_min = (0..m)
                        .map(|s| {
                            let x_s = data1[i + s * n1];
                            (0..m)
                                .map(|t| {
                                    let y_t = data2[j + t * n2];
                                    let val_diff = x_s - y_t;
                                    (val_diff * val_diff + mtt[s * m + t]).sqrt()
                                })
                                .fold(f64::INFINITY, |a, b| a.min(b))
                        })
                        .fold(f64::NEG_INFINITY, |a, b| a.max(b));

                    // For each column t, find min over rows s
                    let max_col_min = (0..m)
                        .map(|t| {
                            let y_t = data2[j + t * n2];
                            (0..m)
                                .map(|s| {
                                    let x_s = data1[i + s * n1];
                                    let val_diff = x_s - y_t;
                                    (val_diff * val_diff + mtt[s * m + t]).sqrt()
                                })
                                .fold(f64::INFINITY, |a, b| a.min(b))
                        })
                        .fold(f64::NEG_INFINITY, |a, b| a.max(b));

                    max_row_min.max(max_col_min)
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let result = RMatrix::new_matrix(n1, n2, |i, j| distances[i * n2 + j]);
    r!(result)
}

/// Compute DTW distance between two time series
///
/// Uses Lp distance with optional Sakoe-Chiba window constraint.
/// This is a direct translation of fda.usc's cumgamma function.
///
/// R code being translated:
/// ```R
/// DTW = matrix(0, nrow = n + 1, ncol = m + 1)
/// DTW[1, ] = Inf
/// DTW[, 1] = Inf
/// DTW[1, 1] = 0
/// for (i in 2:(n + 1)) {
///     for (j in max(2, i - w):min(m + 1, i + w)) {
///         DTW[i, j] = 0
///     }
/// }
/// for (i in 2:(n + 1)) {
///     for (j in max(2, i - w):min(m + 1, i + w)) {
///         DTW[i, j] = D[i - 1, j - 1] + min(DTW[i - 1, j], DTW[i, j - 1], DTW[i - 1, j - 1])
///     }
/// }
/// ```
fn dtw_distance(x: &[f64], y: &[f64], p: f64, w: usize) -> f64 {
    let n = x.len();
    let m = y.len();

    // DTW matrix is (n+1) x (m+1) with 1-based indexing logic
    // Index 0 in Rust = Index 1 in R
    let mut dtw = vec![vec![0.0_f64; m + 1]; n + 1];

    // DTW[1, ] = Inf (R) -> DTW[0][*] = Inf (Rust)
    for j in 0..=m {
        dtw[0][j] = f64::INFINITY;
    }
    // DTW[, 1] = Inf (R) -> DTW[*][0] = Inf (Rust)
    for i in 0..=n {
        dtw[i][0] = f64::INFINITY;
    }
    // DTW[1, 1] = 0 (R) -> DTW[0][0] = 0 (Rust)
    dtw[0][0] = 0.0;

    // First pass: initialize cells within the window band to 0
    // R: for (i in 2:(n + 1))
    // Rust: for i in 1..=n (since our index 1 = R's index 2)
    for i in 1..=n {
        // R: max(2, i - w) where R's i goes from 2 to n+1
        // Our i goes from 1 to n, so R's i = our i + 1
        // R: max(2, (our_i+1) - w) -> Rust: max(1, i + 1 - w) but clamped to >=1
        let r_i = i + 1; // R's i value
        let j_start_r = 2.max(r_i as isize - w as isize) as usize; // R's j_start
        let j_end_r = (m + 1).min(r_i + w); // R's j_end

        // Convert R indices to Rust: R's j -> Rust's j-1
        for j_r in j_start_r..=j_end_r {
            let j = j_r - 1; // Convert to Rust index
            if j <= m {
                dtw[i][j] = 0.0;
            }
        }
    }

    // Second pass: fill the DTW matrix with actual values
    for i in 1..=n {
        let r_i = i + 1;
        let j_start_r = 2.max(r_i as isize - w as isize) as usize;
        let j_end_r = (m + 1).min(r_i + w);

        for j_r in j_start_r..=j_end_r {
            let j = j_r - 1;
            if j <= m && j >= 1 {
                // D[i - 1, j - 1] in R is cost at data indices (i-2, j-2) in R
                // which is (i-1, j-1) in 0-indexed data
                // Our dtw[i][j] corresponds to R's DTW[i+1][j+1]
                // So data indices are (i-1, j-1) for Rust
                let cost = (x[i - 1] - y[j - 1]).abs().powf(p);
                dtw[i][j] = cost + dtw[i - 1][j].min(dtw[i][j - 1]).min(dtw[i - 1][j - 1]);
            }
        }
    }

    dtw[n][m]
}

/// Compute DTW distance matrix for self-distances (symmetric)
#[extendr]
fn metric_dtw_self_1d(fdata: RMatrix<f64>, p: f64, w: i32) -> Robj {
    let n = fdata.nrows();
    let m = fdata.ncols();

    if n == 0 || m == 0 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data = fdata.as_real_slice().unwrap();
    let w = if w < 0 { m } else { w as usize };

    // Extract curves as vectors for easier access
    let curves: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..m).map(|j| data[i + j * n]).collect())
        .collect();

    // Compute upper triangle in parallel
    let upper_triangle: Vec<(usize, usize, f64)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            ((i + 1)..n)
                .map(|j| {
                    let dist = dtw_distance(&curves[i], &curves[j], p, w);
                    (i, j, dist)
                })
                .collect::<Vec<_>>()
        })
        .collect();

    // Build symmetric distance matrix
    let mut distances = vec![0.0; n * n];
    for (i, j, dist) in upper_triangle {
        distances[i + j * n] = dist;
        distances[j + i * n] = dist;
    }

    let result = RMatrix::new_matrix(n, n, |i, j| distances[i + j * n]);
    r!(result)
}

/// Compute DTW distance matrix for cross-distances (n1 x n2)
#[extendr]
fn metric_dtw_cross_1d(fdata1: RMatrix<f64>, fdata2: RMatrix<f64>, p: f64, w: i32) -> Robj {
    let n1 = fdata1.nrows();
    let n2 = fdata2.nrows();
    let m1 = fdata1.ncols();
    let m2 = fdata2.ncols();

    if n1 == 0 || n2 == 0 || m1 == 0 || m2 == 0 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data1 = fdata1.as_real_slice().unwrap();
    let data2 = fdata2.as_real_slice().unwrap();
    let w = if w < 0 { m1.min(m2) } else { w as usize };

    // Extract curves as vectors
    let curves1: Vec<Vec<f64>> = (0..n1)
        .map(|i| (0..m1).map(|j| data1[i + j * n1]).collect())
        .collect();
    let curves2: Vec<Vec<f64>> = (0..n2)
        .map(|i| (0..m2).map(|j| data2[i + j * n2]).collect())
        .collect();

    // Compute all n1 x n2 distances in parallel
    let distances: Vec<f64> = (0..n1)
        .into_par_iter()
        .flat_map(|i| {
            (0..n2)
                .map(|j| dtw_distance(&curves1[i], &curves2[j], p, w))
                .collect::<Vec<_>>()
        })
        .collect();

    let result = RMatrix::new_matrix(n1, n2, |i, j| distances[i * n2 + j]);
    r!(result)
}

// =============================================================================
// 2D metric functions (surfaces)
// =============================================================================

/// Compute Lp distance between two 2D functional data objects (surfaces)
///
/// For 2D functional data, the Lp distance is computed as:
/// d(f, g) = (∫∫ |f(s,t) - g(s,t)|^p ds dt)^(1/p)
///
/// Data is stored as flattened matrices: n x (m1 * m2) where m1 = len(argvals_s), m2 = len(argvals_t)
/// Data is stored in row-major order for the surface: f(s_i, t_j) is at index i * m2 + j
#[extendr]
fn metric_lp_2d(
    fdata1: RMatrix<f64>,
    fdata2: RMatrix<f64>,
    argvals_s: Vec<f64>,
    argvals_t: Vec<f64>,
    p: f64,
    w: Vec<f64>,
) -> Robj {
    let n1 = fdata1.nrows();
    let n2 = fdata2.nrows();
    let ncol1 = fdata1.ncols();
    let ncol2 = fdata2.ncols();
    let m1 = argvals_s.len();
    let m2 = argvals_t.len();
    let n_points = m1 * m2;

    if ncol1 != n_points || ncol2 != n_points || ncol1 != ncol2 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data1 = fdata1.as_real_slice().unwrap();
    let data2 = fdata2.as_real_slice().unwrap();

    // Compute 2D integration weights
    let base_weights = simpsons_weights_2d(&argvals_s, &argvals_t);

    // Combine with user weights
    let weights: Vec<f64> = if w.len() == n_points {
        base_weights
            .iter()
            .zip(w.iter())
            .map(|(bw, uw)| bw * uw)
            .collect()
    } else {
        base_weights
    };

    // Compute pairwise distances
    let distances: Vec<f64> = (0..n1)
        .into_par_iter()
        .flat_map(|i| {
            (0..n2)
                .map(|j| {
                    let mut sum = 0.0;
                    for k in 0..n_points {
                        let diff = (data1[i + k * n1] - data2[j + k * n2]).abs();
                        sum += weights[k] * diff.powf(p);
                    }
                    sum.powf(1.0 / p)
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let result = RMatrix::new_matrix(n1, n2, |i, j| distances[i * n2 + j]);
    r!(result)
}

/// Compute Lp self-distance matrix for 2D functional data (symmetric)
#[extendr]
fn metric_lp_self_2d(
    fdata: RMatrix<f64>,
    argvals_s: Vec<f64>,
    argvals_t: Vec<f64>,
    p: f64,
    w: Vec<f64>,
) -> Robj {
    let n = fdata.nrows();
    let ncol = fdata.ncols();
    let m1 = argvals_s.len();
    let m2 = argvals_t.len();
    let n_points = m1 * m2;

    if ncol != n_points {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data = fdata.as_real_slice().unwrap();

    // Compute 2D integration weights
    let base_weights = simpsons_weights_2d(&argvals_s, &argvals_t);

    // Combine with user weights
    let weights: Vec<f64> = if w.len() == n_points {
        base_weights
            .iter()
            .zip(w.iter())
            .map(|(bw, uw)| bw * uw)
            .collect()
    } else {
        base_weights
    };

    // Compute only upper triangle (symmetric matrix)
    let upper_triangle: Vec<(usize, usize, f64)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            ((i + 1)..n)
                .map(|j| {
                    let mut sum = 0.0;
                    for k in 0..n_points {
                        let diff = (data[i + k * n] - data[j + k * n]).abs();
                        sum += weights[k] * diff.powf(p);
                    }
                    (i, j, sum.powf(1.0 / p))
                })
                .collect::<Vec<_>>()
        })
        .collect();

    // Build symmetric distance matrix
    let mut distances = vec![0.0; n * n];
    for (i, j, dist) in upper_triangle {
        distances[i + j * n] = dist;
        distances[j + i * n] = dist;
    }

    let result = RMatrix::new_matrix(n, n, |i, j| distances[i + j * n]);
    r!(result)
}

/// Compute Hausdorff distance for 2D functional data (surfaces)
///
/// For surfaces, each sample is treated as a point cloud in 3D space:
/// {(s_i, t_j, f(s_i, t_j)) : for all grid points}
///
/// The Hausdorff distance measures how far apart two such point clouds are.
#[extendr]
fn metric_hausdorff_2d(fdata: RMatrix<f64>, argvals_s: Vec<f64>, argvals_t: Vec<f64>) -> Robj {
    let n = fdata.nrows();
    let ncol = fdata.ncols();
    let m1 = argvals_s.len();
    let m2 = argvals_t.len();
    let n_points = m1 * m2;

    if n == 0 || ncol != n_points {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data = fdata.as_real_slice().unwrap();

    // Build 3D point representation for each surface
    // Point (i, j) has coordinates (s[i], t[j], f(s[i], t[j]))
    let surfaces: Vec<Vec<(f64, f64, f64)>> = (0..n)
        .map(|curve| {
            let mut points = Vec::with_capacity(n_points);
            for i in 0..m1 {
                for j in 0..m2 {
                    let k = i * m2 + j;
                    points.push((argvals_s[i], argvals_t[j], data[curve + k * n]));
                }
            }
            points
        })
        .collect();

    // Compute upper triangle of symmetric distance matrix
    let upper_triangle: Vec<(usize, usize, f64)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            ((i + 1)..n)
                .map(|j| {
                    let dist = hausdorff_3d(&surfaces[i], &surfaces[j]);
                    (i, j, dist)
                })
                .collect::<Vec<_>>()
        })
        .collect();

    // Build symmetric distance matrix
    let mut distances = vec![0.0; n * n];
    for (i, j, dist) in upper_triangle {
        distances[i + j * n] = dist;
        distances[j + i * n] = dist;
    }

    let result = RMatrix::new_matrix(n, n, |i, j| distances[i + j * n]);
    r!(result)
}

/// Compute Hausdorff distance between two 3D point clouds
fn hausdorff_3d(points1: &[(f64, f64, f64)], points2: &[(f64, f64, f64)]) -> f64 {
    // Compute directed Hausdorff distance from points1 to points2
    let h12 = points1
        .iter()
        .map(|p1| {
            points2
                .iter()
                .map(|p2| {
                    let ds = p1.0 - p2.0;
                    let dt = p1.1 - p2.1;
                    let df = p1.2 - p2.2;
                    (ds * ds + dt * dt + df * df).sqrt()
                })
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0.0_f64, f64::max);

    // Compute directed Hausdorff distance from points2 to points1
    let h21 = points2
        .iter()
        .map(|p2| {
            points1
                .iter()
                .map(|p1| {
                    let ds = p1.0 - p2.0;
                    let dt = p1.1 - p2.1;
                    let df = p1.2 - p2.2;
                    (ds * ds + dt * dt + df * df).sqrt()
                })
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0.0_f64, f64::max);

    // Hausdorff distance is the maximum of the two directed distances
    h12.max(h21)
}

/// Compute Hausdorff cross-distances for 2D functional data
#[extendr]
fn metric_hausdorff_cross_2d(
    fdata1: RMatrix<f64>,
    fdata2: RMatrix<f64>,
    argvals_s: Vec<f64>,
    argvals_t: Vec<f64>,
) -> Robj {
    let n1 = fdata1.nrows();
    let n2 = fdata2.nrows();
    let m1 = argvals_s.len();
    let m2 = argvals_t.len();
    let n_points = m1 * m2;

    if n1 == 0 || n2 == 0 || fdata1.ncols() != n_points || fdata2.ncols() != n_points {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data1 = fdata1.as_real_slice().unwrap();
    let data2 = fdata2.as_real_slice().unwrap();

    // Build 3D point representation for each surface
    let surfaces1: Vec<Vec<(f64, f64, f64)>> = (0..n1)
        .map(|curve| {
            let mut points = Vec::with_capacity(n_points);
            for i in 0..m1 {
                for j in 0..m2 {
                    let k = i * m2 + j;
                    points.push((argvals_s[i], argvals_t[j], data1[curve + k * n1]));
                }
            }
            points
        })
        .collect();

    let surfaces2: Vec<Vec<(f64, f64, f64)>> = (0..n2)
        .map(|curve| {
            let mut points = Vec::with_capacity(n_points);
            for i in 0..m1 {
                for j in 0..m2 {
                    let k = i * m2 + j;
                    points.push((argvals_s[i], argvals_t[j], data2[curve + k * n2]));
                }
            }
            points
        })
        .collect();

    // Compute n1 x n2 distances
    let distances: Vec<f64> = (0..n1)
        .into_par_iter()
        .flat_map(|i| {
            (0..n2)
                .map(|j| hausdorff_3d(&surfaces1[i], &surfaces2[j]))
                .collect::<Vec<_>>()
        })
        .collect();

    let result = RMatrix::new_matrix(n1, n2, |i, j| distances[i * n2 + j]);
    r!(result)
}

// =============================================================================
// Semimetric functions (FFT-based)
// =============================================================================

/// Compute Fourier coefficients for a curve using FFT
/// Returns magnitude of first nfreq Fourier coefficients
fn fft_coefficients(data: &[f64], nfreq: usize) -> Vec<f64> {
    let n = data.len();
    let nfreq = nfreq.min(n / 2); // Can't have more than n/2 meaningful frequencies

    // Create a planner and get the FFT
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(n);

    // Create complex input buffer
    let mut buffer: Vec<Complex<f64>> = data.iter().map(|&x| Complex::new(x, 0.0)).collect();

    // Perform FFT in place
    fft.process(&mut buffer);

    // Return magnitudes of first nfreq frequencies (including DC)
    // Note: FFT output has DC at index 0, then positive frequencies 1 to n/2
    buffer
        .iter()
        .take(nfreq + 1) // DC + nfreq frequencies
        .map(|c| c.norm() / n as f64) // Normalize by n
        .collect()
}

/// Compute semimetric based on Fourier coefficients for self-distances (symmetric)
/// Uses FFT to compute Fourier coefficients and then L2 distance on coefficients
#[extendr]
fn semimetric_fourier_self_1d(fdata: RMatrix<f64>, nfreq: i32) -> Robj {
    let n = fdata.nrows();
    let m = fdata.ncols();
    let nfreq = nfreq as usize;

    if n == 0 || m == 0 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data = fdata.as_real_slice().unwrap();

    // Extract curves as vectors
    let curves: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..m).map(|j| data[i + j * n]).collect())
        .collect();

    // Compute FFT coefficients for all curves in parallel
    let coeffs: Vec<Vec<f64>> = curves
        .into_par_iter()
        .map(|curve| fft_coefficients(&curve, nfreq))
        .collect();

    // Compute upper triangle of distance matrix in parallel
    let upper_triangle: Vec<(usize, usize, f64)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            ((i + 1)..n)
                .map(|j| {
                    // L2 distance between Fourier coefficients
                    let dist_sq: f64 = coeffs[i]
                        .iter()
                        .zip(coeffs[j].iter())
                        .map(|(a, b)| (a - b).powi(2))
                        .sum();
                    (i, j, dist_sq.sqrt())
                })
                .collect::<Vec<_>>()
        })
        .collect();

    // Build symmetric distance matrix
    let mut distances = vec![0.0; n * n];
    for (i, j, dist) in upper_triangle {
        distances[i + j * n] = dist;
        distances[j + i * n] = dist;
    }

    let result = RMatrix::new_matrix(n, n, |i, j| distances[i + j * n]);
    r!(result)
}

/// Compute semimetric based on Fourier coefficients for cross-distances
#[extendr]
fn semimetric_fourier_cross_1d(fdata1: RMatrix<f64>, fdata2: RMatrix<f64>, nfreq: i32) -> Robj {
    let n1 = fdata1.nrows();
    let n2 = fdata2.nrows();
    let m1 = fdata1.ncols();
    let m2 = fdata2.ncols();
    let nfreq = nfreq as usize;

    if n1 == 0 || n2 == 0 || m1 == 0 || m2 == 0 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    // If curves have different lengths, FFT coefficients might not be directly comparable
    // For simplicity, we proceed assuming same length or truncate
    let m = m1.min(m2);

    let data1 = fdata1.as_real_slice().unwrap();
    let data2 = fdata2.as_real_slice().unwrap();

    // Extract curves as vectors
    let curves1: Vec<Vec<f64>> = (0..n1)
        .map(|i| (0..m).map(|j| data1[i + j * n1]).collect())
        .collect();
    let curves2: Vec<Vec<f64>> = (0..n2)
        .map(|i| (0..m).map(|j| data2[i + j * n2]).collect())
        .collect();

    // Compute FFT coefficients for all curves in parallel
    let coeffs1: Vec<Vec<f64>> = curves1
        .into_par_iter()
        .map(|curve| fft_coefficients(&curve, nfreq))
        .collect();
    let coeffs2: Vec<Vec<f64>> = curves2
        .into_par_iter()
        .map(|curve| fft_coefficients(&curve, nfreq))
        .collect();

    // Compute all n1 x n2 distances in parallel
    let distances: Vec<f64> = (0..n1)
        .into_par_iter()
        .flat_map(|i| {
            (0..n2)
                .map(|j| {
                    let dist_sq: f64 = coeffs1[i]
                        .iter()
                        .zip(coeffs2[j].iter())
                        .map(|(a, b)| (a - b).powi(2))
                        .sum();
                    dist_sq.sqrt()
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let result = RMatrix::new_matrix(n1, n2, |i, j| distances[i * n2 + j]);
    r!(result)
}

/// Compute minimum L2 distance after horizontal shift between two curves
/// Returns the minimum distance over all shifts in [-max_shift, max_shift]
fn hshift_distance(x: &[f64], y: &[f64], weights: &[f64], max_shift: usize) -> f64 {
    let n = x.len();
    if n == 0 || y.len() != n || weights.len() != n {
        return f64::INFINITY;
    }

    let mut min_dist = f64::INFINITY;

    // Try all shifts from -max_shift to +max_shift
    for shift in -(max_shift as i32)..=(max_shift as i32) {
        let mut sum = 0.0;
        let mut valid_points = 0;

        for i in 0..n {
            let j = i as i32 + shift;
            if j >= 0 && (j as usize) < n {
                let diff = x[i] - y[j as usize];
                sum += weights[i] * diff * diff;
                valid_points += 1;
            }
        }

        // Only consider shifts with at least half the points overlapping
        if valid_points >= n / 2 {
            let dist = sum.sqrt();
            if dist < min_dist {
                min_dist = dist;
            }
        }
    }

    min_dist
}

/// Compute semimetric based on horizontal shift for self-distances (symmetric)
/// This finds the minimum L2 distance after optimally shifting one curve horizontally
#[extendr]
fn semimetric_hshift_self_1d(fdata: RMatrix<f64>, argvals: Vec<f64>, max_shift: i32) -> Robj {
    let n = fdata.nrows();
    let m = fdata.ncols();

    if n == 0 || m == 0 || argvals.len() != m {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data = fdata.as_real_slice().unwrap();
    let weights = simpsons_weights(&argvals);
    let max_shift = if max_shift < 0 {
        m / 4
    } else {
        max_shift as usize
    };

    // Extract curves as vectors
    let curves: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..m).map(|j| data[i + j * n]).collect())
        .collect();

    // Compute upper triangle of distance matrix in parallel
    let upper_triangle: Vec<(usize, usize, f64)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            ((i + 1)..n)
                .map(|j| {
                    let dist = hshift_distance(&curves[i], &curves[j], &weights, max_shift);
                    (i, j, dist)
                })
                .collect::<Vec<_>>()
        })
        .collect();

    // Build symmetric distance matrix
    let mut distances = vec![0.0; n * n];
    for (i, j, dist) in upper_triangle {
        distances[i + j * n] = dist;
        distances[j + i * n] = dist;
    }

    let result = RMatrix::new_matrix(n, n, |i, j| distances[i + j * n]);
    r!(result)
}

/// Compute semimetric based on horizontal shift for cross-distances
#[extendr]
fn semimetric_hshift_cross_1d(
    fdata1: RMatrix<f64>,
    fdata2: RMatrix<f64>,
    argvals: Vec<f64>,
    max_shift: i32,
) -> Robj {
    let n1 = fdata1.nrows();
    let n2 = fdata2.nrows();
    let m1 = fdata1.ncols();
    let m2 = fdata2.ncols();

    if n1 == 0 || n2 == 0 || m1 == 0 || m2 == 0 || m1 != m2 || argvals.len() != m1 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let m = m1;
    let data1 = fdata1.as_real_slice().unwrap();
    let data2 = fdata2.as_real_slice().unwrap();
    let weights = simpsons_weights(&argvals);
    let max_shift = if max_shift < 0 {
        m / 4
    } else {
        max_shift as usize
    };

    // Extract curves as vectors
    let curves1: Vec<Vec<f64>> = (0..n1)
        .map(|i| (0..m).map(|j| data1[i + j * n1]).collect())
        .collect();
    let curves2: Vec<Vec<f64>> = (0..n2)
        .map(|i| (0..m).map(|j| data2[i + j * n2]).collect())
        .collect();

    // Compute all n1 x n2 distances in parallel
    let distances: Vec<f64> = (0..n1)
        .into_par_iter()
        .flat_map(|i| {
            (0..n2)
                .map(|j| hshift_distance(&curves1[i], &curves2[j], &weights, max_shift))
                .collect::<Vec<_>>()
        })
        .collect();

    let result = RMatrix::new_matrix(n1, n2, |i, j| distances[i * n2 + j]);
    r!(result)
}

// =============================================================================
// KL divergence metric
// =============================================================================

/// Helper function to compute KL divergence between two normalized curves
/// KL(p||q) = sum_i p_i * log(p_i / q_i)
/// Uses symmetric KL: (KL(p||q) + KL(q||p)) / 2
fn kl_divergence(p: &[f64], q: &[f64], weights: &[f64], eps: f64) -> f64 {
    let mut kl_pq = 0.0;
    let mut kl_qp = 0.0;

    for i in 0..p.len() {
        let p_i = (p[i] + eps).max(eps);
        let q_i = (q[i] + eps).max(eps);
        kl_pq += weights[i] * p_i * (p_i / q_i).ln();
        kl_qp += weights[i] * q_i * (q_i / p_i).ln();
    }

    // Return symmetric KL divergence
    (kl_pq + kl_qp) / 2.0
}

/// Compute symmetric KL divergence matrix for self-distances (1D)
/// Curves are first normalized to be valid probability distributions
#[extendr]
fn metric_kl_self_1d(fdata: RMatrix<f64>, argvals: Vec<f64>, eps: f64, normalize: bool) -> Robj {
    let n = fdata.nrows();
    let m = fdata.ncols();

    if n == 0 || m == 0 || argvals.len() != m {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data = fdata.as_real_slice().unwrap();
    let weights = simpsons_weights(&argvals);
    let eps = if eps <= 0.0 { 1e-10 } else { eps };

    // Extract and optionally normalize curves
    let curves: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let curve: Vec<f64> = (0..m).map(|j| data[i + j * n]).collect();
            if normalize {
                // Normalize to make it a probability distribution
                let min_val = curve.iter().cloned().fold(f64::INFINITY, f64::min);
                let shifted: Vec<f64> = curve.iter().map(|&x| x - min_val + eps).collect();
                let sum: f64 = shifted.iter().zip(weights.iter()).map(|(v, w)| v * w).sum();
                shifted.iter().map(|&x| x / sum).collect()
            } else {
                curve
            }
        })
        .collect();

    // Compute upper triangle of symmetric distance matrix
    let pairs: Vec<(usize, usize)> = (0..n).flat_map(|i| (i..n).map(move |j| (i, j))).collect();

    let distances: Vec<(usize, usize, f64)> = pairs
        .into_par_iter()
        .map(|(i, j)| {
            let d = if i == j {
                0.0
            } else {
                kl_divergence(&curves[i], &curves[j], &weights, eps)
            };
            (i, j, d)
        })
        .collect();

    // Build symmetric matrix
    let mut dist_matrix = vec![0.0; n * n];
    for (i, j, d) in distances {
        dist_matrix[i * n + j] = d;
        dist_matrix[j * n + i] = d;
    }

    let result = RMatrix::new_matrix(n, n, |i, j| dist_matrix[i * n + j]);
    r!(result)
}

/// Compute symmetric KL divergence matrix for cross-distances (1D)
#[extendr]
fn metric_kl_cross_1d(
    fdata1: RMatrix<f64>,
    fdata2: RMatrix<f64>,
    argvals: Vec<f64>,
    eps: f64,
    normalize: bool,
) -> Robj {
    let n1 = fdata1.nrows();
    let n2 = fdata2.nrows();
    let m1 = fdata1.ncols();
    let m2 = fdata2.ncols();

    if n1 == 0 || n2 == 0 || m1 == 0 || m2 == 0 || m1 != m2 || argvals.len() != m1 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let m = m1;
    let data1 = fdata1.as_real_slice().unwrap();
    let data2 = fdata2.as_real_slice().unwrap();
    let weights = simpsons_weights(&argvals);
    let eps = if eps <= 0.0 { 1e-10 } else { eps };

    // Extract and optionally normalize curves
    let normalize_curve = |curve: Vec<f64>| -> Vec<f64> {
        if normalize {
            let min_val = curve.iter().cloned().fold(f64::INFINITY, f64::min);
            let shifted: Vec<f64> = curve.iter().map(|&x| x - min_val + eps).collect();
            let sum: f64 = shifted.iter().zip(weights.iter()).map(|(v, w)| v * w).sum();
            shifted.iter().map(|&x| x / sum).collect()
        } else {
            curve
        }
    };

    let curves1: Vec<Vec<f64>> = (0..n1)
        .map(|i| {
            let curve: Vec<f64> = (0..m).map(|j| data1[i + j * n1]).collect();
            normalize_curve(curve)
        })
        .collect();

    let curves2: Vec<Vec<f64>> = (0..n2)
        .map(|i| {
            let curve: Vec<f64> = (0..m).map(|j| data2[i + j * n2]).collect();
            normalize_curve(curve)
        })
        .collect();

    // Compute all n1 x n2 distances in parallel
    let distances: Vec<f64> = (0..n1)
        .into_par_iter()
        .flat_map(|i| {
            (0..n2)
                .map(|j| kl_divergence(&curves1[i], &curves2[j], &weights, eps))
                .collect::<Vec<_>>()
        })
        .collect();

    let result = RMatrix::new_matrix(n1, n2, |i, j| distances[i * n2 + j]);
    r!(result)
}

// =============================================================================
// stats functions
// =============================================================================

/// Compute the Adot matrix (parallelized)
#[extendr]
fn compute_adot(n: i32, inprod: Vec<f64>) -> Robj {
    let n = n as usize;

    if n == 0 {
        return Robj::from(Vec::<f64>::new());
    }

    let expected_len = (n * n + n) / 2;
    if inprod.len() != expected_len {
        return Robj::from(Vec::<f64>::new());
    }

    let out_len = (n * n - n + 2) / 2;
    let mut adot_vec = vec![0.0; out_len];

    adot_vec[0] = PI * (n + 1) as f64;

    // Collect all (i, j) pairs for parallel processing
    let pairs: Vec<(usize, usize)> = (2..=n).flat_map(|i| (1..i).map(move |j| (i, j))).collect();

    // Compute adot values in parallel
    let results: Vec<(usize, f64)> = pairs
        .into_par_iter()
        .map(|(i, j)| {
            let mut sumr = 0.0;

            for r in 1..=n {
                if i == r || j == r {
                    sumr += PI;
                } else {
                    let auxi = i * (i - 1) / 2;
                    let auxj = j * (j - 1) / 2;
                    let auxr = r * (r - 1) / 2;

                    let ij = auxi + j - 1;
                    let ii = auxi + i - 1;
                    let jj = auxj + j - 1;
                    let rr = auxr + r - 1;

                    let ir = if i > r { auxi + r - 1 } else { auxr + i - 1 };
                    let rj = if r > j { auxr + j - 1 } else { auxj + r - 1 };
                    let jr = rj;

                    let num = inprod[ij] - inprod[ir] - inprod[rj] + inprod[rr];
                    let aux1 = (inprod[ii] - 2.0 * inprod[ir] + inprod[rr]).sqrt();
                    let aux2 = (inprod[jj] - 2.0 * inprod[jr] + inprod[rr]).sqrt();
                    let den = aux1 * aux2;

                    let mut quo = if den.abs() > 1e-10 { num / den } else { 0.0 };
                    quo = quo.clamp(-1.0, 1.0);

                    sumr += (PI - quo.acos()).abs();
                }
            }

            let idx = 1 + ((i - 1) * (i - 2) / 2) + j - 1;
            (idx, sumr)
        })
        .collect();

    // Fill in the results
    for (idx, val) in results {
        adot_vec[idx] = val;
    }

    Robj::from(adot_vec)
}

/// Compute the PCvM statistic
#[extendr]
fn pcvm_statistic(adot_vec: Vec<f64>, residuals: Vec<f64>) -> f64 {
    let n = residuals.len();

    if n == 0 || adot_vec.is_empty() {
        return 0.0;
    }

    let mut sums = 0.0;
    for i in 2..=n {
        for j in 1..i {
            let idx = 1 + ((i - 1) * (i - 2) / 2) + j - 1;
            if idx < adot_vec.len() {
                sums += residuals[i - 1] * adot_vec[idx] * residuals[j - 1];
            }
        }
    }

    let diag_sum: f64 = residuals.iter().map(|r| r * r).sum();
    adot_vec[0] * diag_sum + 2.0 * sums
}

/// Compute random projection statistics (parallelized over projections)
#[extendr]
fn rp_stat(proj_x_ord: Vec<i32>, residuals: Vec<f64>, n_proj: i32) -> Robj {
    let n = residuals.len();
    let n_proj = n_proj as usize;

    if n == 0 || n_proj == 0 {
        return list!(cvm = Vec::<f64>::new(), ks = Vec::<f64>::new()).into();
    }

    // Process projections in parallel
    let stats: Vec<(f64, f64)> = (0..n_proj)
        .into_par_iter()
        .map(|p| {
            let mut y = vec![0.0; n];
            let mut cumsum = 0.0;

            for i in 0..n {
                let idx = proj_x_ord[p * n + i] as usize - 1;
                if idx < n {
                    cumsum += residuals[idx];
                }
                y[i] = cumsum;
            }

            let sum_y_sq: f64 = y.iter().map(|yi| yi * yi).sum();
            let cvm = sum_y_sq / (n * n) as f64;

            let max_abs_y = y.iter().map(|yi| yi.abs()).fold(0.0, f64::max);
            let ks = max_abs_y / (n as f64).sqrt();

            (cvm, ks)
        })
        .collect();

    let cvm_stats: Vec<f64> = stats.iter().map(|(cvm, _)| *cvm).collect();
    let ks_stats: Vec<f64> = stats.iter().map(|(_, ks)| *ks).collect();

    list!(cvm = cvm_stats, ks = ks_stats).into()
}

// =============================================================================
// PCA / PLS / Basis functions
// =============================================================================

/// Perform functional PCA via SVD on centered data
/// Returns: singular values, rotation matrix (loadings), scores, mean
#[extendr]
fn fdata2pc_1d(data: RMatrix<f64>, ncomp: i32, _lambda: f64) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m == 0 || ncomp < 1 {
        return list!(
            d = Vec::<f64>::new(),
            rotation = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            scores = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            mean = Vec::<f64>::new(),
            centered = RMatrix::new_matrix(0, 0, |_, _| 0.0)
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let ncomp = (ncomp as usize).min(n).min(m);

    // Compute column means
    let means: Vec<f64> = (0..m)
        .into_par_iter()
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..n {
                sum += data_slice[i + j * n];
            }
            sum / n as f64
        })
        .collect();

    // Center the data and convert to nalgebra matrix (column-major)
    let centered_data: Vec<f64> = (0..(n * m))
        .map(|idx| {
            let i = idx % n;
            let j = idx / n;
            data_slice[i + j * n] - means[j]
        })
        .collect();

    // Create nalgebra DMatrix (column-major, nrows x ncols)
    let matrix = DMatrix::from_column_slice(n, m, &centered_data);

    // Compute SVD: X = U * S * V^T
    let svd = SVD::new(matrix, true, true);

    // Extract singular values
    let singular_values: Vec<f64> = svd.singular_values.iter().take(ncomp).cloned().collect();

    // Extract V (right singular vectors, the loadings/rotation)
    // V is m x min(n, m), we want first ncomp columns
    let v_t = svd.v_t.as_ref().unwrap();
    let rotation_data: Vec<f64> = (0..ncomp)
        .flat_map(|k| (0..m).map(move |j| v_t[(k, j)]))
        .collect();

    // Compute scores: X_centered * V = U * S
    // scores[i, k] = sum_j centered_data[i, j] * rotation[j, k]
    let u = svd.u.as_ref().unwrap();
    let mut scores_data: Vec<f64> = Vec::with_capacity(n * ncomp);
    for k in 0..ncomp {
        let sv_k = singular_values[k];
        for i in 0..n {
            scores_data.push(u[(i, k)] * sv_k);
        }
    }

    // Return results
    let rotation = RMatrix::new_matrix(m, ncomp, |j, k| rotation_data[k * m + j]);
    let scores = RMatrix::new_matrix(n, ncomp, |i, k| scores_data[k * n + i]);
    let centered = RMatrix::new_matrix(n, m, |i, j| centered_data[i + j * n]);

    list!(
        d = singular_values,
        rotation = rotation,
        scores = scores,
        mean = means,
        centered = centered
    )
    .into()
}

/// Perform PLS via NIPALS algorithm
/// Returns: weights, scores, loadings
#[extendr]
fn fdata2pls_1d(data: RMatrix<f64>, y: Vec<f64>, ncomp: i32, _lambda: f64) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m == 0 || y.len() != n || ncomp < 1 {
        return list!(
            weights = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            scores = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            loadings = RMatrix::new_matrix(0, 0, |_, _| 0.0)
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let ncomp = (ncomp as usize).min(n).min(m);

    // Center X and y
    let x_means: Vec<f64> = (0..m)
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..n {
                sum += data_slice[i + j * n];
            }
            sum / n as f64
        })
        .collect();

    let y_mean: f64 = y.iter().sum::<f64>() / n as f64;

    // Centered data
    let mut x_cen: Vec<f64> = (0..(n * m))
        .map(|idx| {
            let i = idx % n;
            let j = idx / n;
            data_slice[i + j * n] - x_means[j]
        })
        .collect();

    let mut y_cen: Vec<f64> = y.iter().map(|&yi| yi - y_mean).collect();

    let mut weights = vec![0.0; m * ncomp];
    let mut scores = vec![0.0; n * ncomp];
    let mut loadings = vec![0.0; m * ncomp];

    // NIPALS algorithm
    for k in 0..ncomp {
        // w = X'y / ||X'y||
        let mut w: Vec<f64> = (0..m)
            .map(|j| {
                let mut sum = 0.0;
                for i in 0..n {
                    sum += x_cen[i + j * n] * y_cen[i];
                }
                sum
            })
            .collect();

        let w_norm: f64 = w.iter().map(|&wi| wi * wi).sum::<f64>().sqrt();
        if w_norm > 1e-10 {
            for wi in &mut w {
                *wi /= w_norm;
            }
        }

        // t = Xw
        let t: Vec<f64> = (0..n)
            .map(|i| {
                let mut sum = 0.0;
                for j in 0..m {
                    sum += x_cen[i + j * n] * w[j];
                }
                sum
            })
            .collect();

        let t_norm_sq: f64 = t.iter().map(|&ti| ti * ti).sum();

        // p = X't / (t't)
        let p: Vec<f64> = (0..m)
            .map(|j| {
                let mut sum = 0.0;
                for i in 0..n {
                    sum += x_cen[i + j * n] * t[i];
                }
                sum / t_norm_sq.max(1e-10)
            })
            .collect();

        // Store results
        for j in 0..m {
            weights[j + k * m] = w[j];
            loadings[j + k * m] = p[j];
        }
        for i in 0..n {
            scores[i + k * n] = t[i];
        }

        // Deflate X: X = X - t * p'
        for j in 0..m {
            for i in 0..n {
                x_cen[i + j * n] -= t[i] * p[j];
            }
        }

        // Deflate y: y = y - t * (t'y / t't)
        let t_y: f64 = t.iter().zip(y_cen.iter()).map(|(&ti, &yi)| ti * yi).sum();
        let q = t_y / t_norm_sq.max(1e-10);
        for i in 0..n {
            y_cen[i] -= t[i] * q;
        }
    }

    let weights_mat = RMatrix::new_matrix(m, ncomp, |j, k| weights[j + k * m]);
    let scores_mat = RMatrix::new_matrix(n, ncomp, |i, k| scores[i + k * n]);
    let loadings_mat = RMatrix::new_matrix(m, ncomp, |j, k| loadings[j + k * m]);

    list!(
        weights = weights_mat,
        scores = scores_mat,
        loadings = loadings_mat
    )
    .into()
}

/// Compute B-spline basis matrix for given knots and grid points
fn bspline_basis(t: &[f64], nknots: usize, order: usize) -> Vec<f64> {
    let n = t.len();
    // Total knots: order (left) + nknots (interior) + order (right) = 2*order + nknots
    // Number of B-spline basis functions: total_knots - order = nknots + order
    let nbasis = nknots + order;

    // Create knot vector with appropriate padding
    let t_min = t.iter().cloned().fold(f64::INFINITY, f64::min);
    let t_max = t.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let dt = (t_max - t_min) / (nknots - 1) as f64;

    let mut knots = Vec::with_capacity(nknots + 2 * order);
    // Pad at beginning (order knots)
    for i in 0..order {
        knots.push(t_min - (order - i) as f64 * dt);
    }
    // Interior knots
    for i in 0..nknots {
        knots.push(t_min + i as f64 * dt);
    }
    // Pad at end (order knots)
    for i in 1..=order {
        knots.push(t_max + i as f64 * dt);
    }

    // Index of t_max in knot vector: it's the last interior knot
    let t_max_knot_idx = order + nknots - 1;

    // Cox-de Boor recursion
    let mut basis = vec![0.0; n * nbasis];

    for (ti, &t_val) in t.iter().enumerate() {
        // Initialize order 1 (constant) splines
        let mut b0 = vec![0.0; knots.len() - 1];
        // Find which interval t_val belongs to
        // Use half-open intervals [knots[j], knots[j+1]) except at t_max
        // where we use the closed interval [knots[j], knots[j+1]]
        for j in 0..(knots.len() - 1) {
            let in_interval = if j == t_max_knot_idx - 1 {
                // Last interior interval: use closed [t_max - dt, t_max]
                t_val >= knots[j] && t_val <= knots[j + 1]
            } else {
                // Normal half-open interval [knots[j], knots[j+1])
                t_val >= knots[j] && t_val < knots[j + 1]
            };

            if in_interval {
                b0[j] = 1.0;
                break;
            }
        }

        // Recursively compute higher orders
        let mut b = b0;
        for k in 2..=order {
            let mut b_new = vec![0.0; knots.len() - k];
            for j in 0..(knots.len() - k) {
                let d1 = knots[j + k - 1] - knots[j];
                let d2 = knots[j + k] - knots[j + 1];

                let left = if d1.abs() > 1e-10 {
                    (t_val - knots[j]) / d1 * b[j]
                } else {
                    0.0
                };
                let right = if d2.abs() > 1e-10 {
                    (knots[j + k] - t_val) / d2 * b[j + 1]
                } else {
                    0.0
                };
                b_new[j] = left + right;
            }
            b = b_new;
        }

        // Store the nbasis values
        for j in 0..nbasis {
            basis[ti + j * n] = b[j];
        }
    }

    basis
}

/// Compute Fourier basis matrix
fn fourier_basis(t: &[f64], nbasis: usize) -> Vec<f64> {
    let n = t.len();
    let t_min = t.iter().cloned().fold(f64::INFINITY, f64::min);
    let t_max = t.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let period = t_max - t_min;

    let mut basis = vec![0.0; n * nbasis];

    for (i, &ti) in t.iter().enumerate() {
        let x = 2.0 * PI * (ti - t_min) / period;

        // First basis function is constant
        basis[i] = 1.0;

        // Remaining basis functions are sin/cos pairs
        let mut k = 1;
        let mut freq = 1;
        while k < nbasis {
            if k < nbasis {
                basis[i + k * n] = (freq as f64 * x).sin();
                k += 1;
            }
            if k < nbasis {
                basis[i + k * n] = (freq as f64 * x).cos();
                k += 1;
            }
            freq += 1;
        }
    }

    basis
}

/// Convert functional data to basis coefficients
/// type: 0 = bspline, 1 = fourier
#[extendr]
fn fdata2basis_1d(data: RMatrix<f64>, argvals: Vec<f64>, nbasis: i32, basis_type: i32) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let nbasis = nbasis as usize;

    if n == 0 || m == 0 || argvals.len() != m || nbasis < 2 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data_slice = data.as_real_slice().unwrap();

    // Compute basis matrix (m x nbasis)
    let basis = if basis_type == 1 {
        fourier_basis(&argvals, nbasis)
    } else {
        // For order 4 B-splines: nbasis = nknots + order, so nknots = nbasis - 4
        // Minimum nknots is 2, so minimum nbasis is 6 for cubic B-splines
        bspline_basis(&argvals, (nbasis.saturating_sub(4)).max(2), 4)
    };

    let actual_nbasis = basis.len() / m;

    // Create nalgebra matrix for basis
    let b_mat = DMatrix::from_column_slice(m, actual_nbasis, &basis);

    // Compute (B'B)^-1 B' for least squares: coefs = (B'B)^-1 B' X'
    let btb = &b_mat.transpose() * &b_mat;

    // Use SVD for pseudo-inverse (more stable than direct inversion)
    let btb_svd = SVD::new(btb, true, true);

    // Check for near-singular matrix
    let max_sv = btb_svd.singular_values.iter().cloned().fold(0.0, f64::max);
    let eps = 1e-10 * max_sv;

    // Pseudo-inverse of singular values
    let s_inv: Vec<f64> = btb_svd
        .singular_values
        .iter()
        .map(|&s| if s > eps { 1.0 / s } else { 0.0 })
        .collect();

    // Compute pseudo-inverse: V * S_inv * U^T
    let v = btb_svd.v_t.as_ref().unwrap().transpose();
    let u_t = btb_svd.u.as_ref().unwrap().transpose();

    // btb_inv = V * diag(s_inv) * U^T
    let mut btb_inv = DMatrix::zeros(actual_nbasis, actual_nbasis);
    for i in 0..actual_nbasis {
        for j in 0..actual_nbasis {
            let mut sum = 0.0;
            for k in 0..actual_nbasis.min(s_inv.len()) {
                sum += v[(i, k)] * s_inv[k] * u_t[(k, j)];
            }
            btb_inv[(i, j)] = sum;
        }
    }

    // proj = btb_inv * B^T
    let proj = btb_inv * b_mat.transpose();

    // Compute coefficients for each curve (parallelized)
    let coefs: Vec<f64> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            // Extract curve i
            let curve: Vec<f64> = (0..m).map(|j| data_slice[i + j * n]).collect();

            // coefs_i = proj * curve
            (0..actual_nbasis)
                .map(|k| {
                    let mut sum = 0.0;
                    for j in 0..m {
                        sum += proj[(k, j)] * curve[j];
                    }
                    sum
                })
                .collect::<Vec<_>>()
        })
        .collect();

    // Return n x nbasis coefficient matrix
    let result = RMatrix::new_matrix(n, actual_nbasis, |i, k| coefs[i * actual_nbasis + k]);
    r!(result)
}

// =============================================================================
// Basis representation: reconstruction, GCV, AIC, BIC, P-splines
// =============================================================================

/// Reconstruct functional data from basis coefficients
/// Returns data matrix (n x m)
#[extendr]
fn basis2fdata_1d(coefs: RMatrix<f64>, argvals: Vec<f64>, nbasis: i32, basis_type: i32) -> Robj {
    let n = coefs.nrows();
    let nbasis_in = coefs.ncols();
    let m = argvals.len();
    let nbasis = nbasis as usize;

    if n == 0 || m == 0 || nbasis < 2 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let coefs_slice = coefs.as_real_slice().unwrap();

    // Compute basis matrix (m x nbasis)
    let basis = if basis_type == 1 {
        fourier_basis(&argvals, nbasis)
    } else {
        // For order 4 B-splines: nbasis = nknots + order, so nknots = nbasis - 4
        bspline_basis(&argvals, (nbasis.saturating_sub(4)).max(2), 4)
    };

    let actual_nbasis = basis.len() / m;

    // Reconstruct: data = coefs * B^T
    // coefs: [n x nbasis], B: [m x nbasis] => data: [n x m]
    let data: Vec<f64> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            (0..m)
                .map(|j| {
                    let mut sum = 0.0;
                    for k in 0..actual_nbasis.min(nbasis_in) {
                        sum += coefs_slice[i + k * n] * basis[j + k * m];
                    }
                    sum
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let result = RMatrix::new_matrix(n, m, |i, j| data[i * m + j]);
    r!(result)
}

/// Compute difference matrix for P-spline penalty
/// order 1: first differences, order 2: second differences, etc.
fn difference_matrix(n: usize, order: usize) -> DMatrix<f64> {
    if order == 0 {
        return DMatrix::identity(n, n);
    }

    // First-order difference matrix [n-1 x n]
    let mut d = DMatrix::zeros(n - 1, n);
    for i in 0..(n - 1) {
        d[(i, i)] = -1.0;
        d[(i, i + 1)] = 1.0;
    }

    // Apply recursively for higher orders
    let mut result = d;
    for _ in 1..order {
        if result.nrows() <= 1 {
            break;
        }
        let rows = result.nrows() - 1;
        let cols = result.ncols();
        let mut d_next = DMatrix::zeros(rows, cols);
        for i in 0..rows {
            for j in 0..cols {
                d_next[(i, j)] = -result[(i, j)] + result[(i + 1, j)];
            }
        }
        result = d_next;
    }

    result
}

/// Compute GCV score for basis fit
/// GCV = RSS/n / (1 - edf/n)^2
/// When pooled=true: compute single GCV across all curves
/// When pooled=false: compute per-curve GCV and return mean
#[extendr]
fn basis_gcv_1d(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    nbasis: i32,
    basis_type: i32,
    lambda: f64,
    pooled: bool,
) -> f64 {
    let n = data.nrows();
    let m = data.ncols();
    let nbasis = nbasis as usize;

    if n == 0 || m == 0 || nbasis < 2 {
        return f64::INFINITY;
    }

    let data_slice = data.as_real_slice().unwrap();

    // Compute basis matrix
    let basis = if basis_type == 1 {
        fourier_basis(&argvals, nbasis)
    } else {
        // For order 4 B-splines: nbasis = nknots + order, so nknots = nbasis - 4
        bspline_basis(&argvals, (nbasis.saturating_sub(4)).max(2), 4)
    };

    let actual_nbasis = basis.len() / m;
    let b_mat = DMatrix::from_column_slice(m, actual_nbasis, &basis);

    // Compute B'B
    let btb = &b_mat.transpose() * &b_mat;

    // Add penalty if lambda > 0
    let btb_penalized = if lambda > 0.0 {
        let d = difference_matrix(actual_nbasis, 2);
        let penalty = &d.transpose() * &d;
        btb + lambda * penalty
    } else {
        btb
    };

    // Solve for coefficients using SVD pseudoinverse (handles singular matrices)
    let svd = SVD::new(btb_penalized.clone(), true, true);
    let eps = 1e-10 * svd.singular_values.iter().cloned().fold(0.0_f64, f64::max);

    let u = match svd.u.as_ref() {
        Some(u) => u,
        None => return f64::INFINITY,
    };
    let v_t = match svd.v_t.as_ref() {
        Some(v_t) => v_t,
        None => return f64::INFINITY,
    };

    // Compute pseudoinverse: V * S^-1 * U^T
    let s_inv: Vec<f64> = svd
        .singular_values
        .iter()
        .map(|&s| if s > eps { 1.0 / s } else { 0.0 })
        .collect();

    let mut btb_inv = DMatrix::zeros(actual_nbasis, actual_nbasis);
    for i in 0..actual_nbasis {
        for j in 0..actual_nbasis {
            let mut sum = 0.0;
            for k in 0..actual_nbasis.min(s_inv.len()) {
                sum += v_t[(k, i)] * s_inv[k] * u[(j, k)];
            }
            btb_inv[(i, j)] = sum;
        }
    }

    // Hat matrix H = B * (B'B + lambda*P)^-1 * B'
    // edf = trace(H) = trace(B * btb_inv * B')
    let proj = &btb_inv * b_mat.transpose();
    let h_mat = &b_mat * &proj;
    let edf: f64 = (0..m).map(|i| h_mat[(i, i)]).sum();

    // GCV denominator (same for all curves since hat matrix is shared)
    let gcv_denom = 1.0 - edf / m as f64;
    if gcv_denom.abs() < 1e-10 {
        return f64::INFINITY;
    }
    let gcv_denom_sq = gcv_denom * gcv_denom;

    if pooled {
        // Pooled mode: compute single GCV across all curves
        let mut total_rss = 0.0;
        let total_points = n * m;

        for i in 0..n {
            let curve: Vec<f64> = (0..m).map(|j| data_slice[i + j * n]).collect();
            let curve_vec = DVector::from_vec(curve.clone());
            let bt_y = b_mat.transpose() * &curve_vec;
            let coefs = &btb_inv * bt_y;
            let fitted = &b_mat * &coefs;

            for j in 0..m {
                let resid = curve[j] - fitted[j];
                total_rss += resid * resid;
            }
        }

        // GCV = (RSS/n_total) / (1 - edf/m)^2
        (total_rss / total_points as f64) / gcv_denom_sq
    } else {
        // Per-curve mode: compute GCV for each curve, return mean
        let mut gcv_sum = 0.0;

        for i in 0..n {
            let curve: Vec<f64> = (0..m).map(|j| data_slice[i + j * n]).collect();
            let curve_vec = DVector::from_vec(curve.clone());
            let bt_y = b_mat.transpose() * &curve_vec;
            let coefs = &btb_inv * bt_y;
            let fitted = &b_mat * &coefs;

            let mut curve_rss = 0.0;
            for j in 0..m {
                let resid = curve[j] - fitted[j];
                curve_rss += resid * resid;
            }

            // GCV for this curve
            gcv_sum += (curve_rss / m as f64) / gcv_denom_sq;
        }

        gcv_sum / n as f64
    }
}

/// Compute AIC for basis fit
/// AIC = n * log(RSS/n) + 2 * total_edf
/// Where total_edf = n_curves * edf (each curve has edf parameters)
/// When pooled=true: compute single AIC across all curves
/// When pooled=false: compute per-curve AIC and return mean
#[extendr]
fn basis_aic_1d(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    nbasis: i32,
    basis_type: i32,
    lambda: f64,
    pooled: bool,
) -> f64 {
    let n = data.nrows();
    let m = data.ncols();
    let nbasis = nbasis as usize;

    if n == 0 || m == 0 || nbasis < 2 {
        return f64::INFINITY;
    }

    let data_slice = data.as_real_slice().unwrap();

    // Compute basis matrix
    let basis = if basis_type == 1 {
        fourier_basis(&argvals, nbasis)
    } else {
        // For order 4 B-splines: nbasis = nknots + order, so nknots = nbasis - 4
        bspline_basis(&argvals, (nbasis.saturating_sub(4)).max(2), 4)
    };

    let actual_nbasis = basis.len() / m;
    let b_mat = DMatrix::from_column_slice(m, actual_nbasis, &basis);

    let btb = &b_mat.transpose() * &b_mat;
    let btb_penalized = if lambda > 0.0 {
        let d = difference_matrix(actual_nbasis, 2);
        let penalty = &d.transpose() * &d;
        btb + lambda * penalty
    } else {
        btb
    };

    // Solve for coefficients using SVD pseudoinverse (handles singular matrices)
    let svd = SVD::new(btb_penalized.clone(), true, true);
    let eps = 1e-10 * svd.singular_values.iter().cloned().fold(0.0_f64, f64::max);

    let u = match svd.u.as_ref() {
        Some(u) => u,
        None => return f64::INFINITY,
    };
    let v_t = match svd.v_t.as_ref() {
        Some(v_t) => v_t,
        None => return f64::INFINITY,
    };

    let s_inv: Vec<f64> = svd
        .singular_values
        .iter()
        .map(|&s| if s > eps { 1.0 / s } else { 0.0 })
        .collect();

    let mut btb_inv = DMatrix::zeros(actual_nbasis, actual_nbasis);
    for i in 0..actual_nbasis {
        for j in 0..actual_nbasis {
            let mut sum = 0.0;
            for k in 0..actual_nbasis.min(s_inv.len()) {
                sum += v_t[(k, i)] * s_inv[k] * u[(j, k)];
            }
            btb_inv[(i, j)] = sum;
        }
    }

    // EDF (effective degrees of freedom per curve)
    let proj = &btb_inv * b_mat.transpose();
    let h_mat = &b_mat * &proj;
    let edf: f64 = (0..m).map(|i| h_mat[(i, i)]).sum();

    if pooled {
        // Pooled mode: compute single AIC across all curves
        let mut total_rss = 0.0;
        let total_points = n * m;

        for i in 0..n {
            let curve: Vec<f64> = (0..m).map(|j| data_slice[i + j * n]).collect();
            let curve_vec = DVector::from_vec(curve.clone());
            let bt_y = b_mat.transpose() * &curve_vec;
            let coefs = &btb_inv * bt_y;
            let fitted = &b_mat * &coefs;

            for j in 0..m {
                let resid = curve[j] - fitted[j];
                total_rss += resid * resid;
            }
        }

        // AIC = n_total * log(RSS/n_total) + 2 * total_edf
        // total_edf = n_curves * edf (each curve has edf effective parameters)
        let mse = total_rss / total_points as f64;
        let total_edf = (n as f64) * edf;
        (total_points as f64) * mse.ln() + 2.0 * total_edf
    } else {
        // Per-curve mode: compute AIC for each curve, return mean
        let mut aic_sum = 0.0;

        for i in 0..n {
            let curve: Vec<f64> = (0..m).map(|j| data_slice[i + j * n]).collect();
            let curve_vec = DVector::from_vec(curve.clone());
            let bt_y = b_mat.transpose() * &curve_vec;
            let coefs = &btb_inv * bt_y;
            let fitted = &b_mat * &coefs;

            let mut curve_rss = 0.0;
            for j in 0..m {
                let resid = curve[j] - fitted[j];
                curve_rss += resid * resid;
            }

            // AIC for this curve: m * log(RSS/m) + 2 * edf
            let curve_mse = curve_rss / m as f64;
            aic_sum += (m as f64) * curve_mse.ln() + 2.0 * edf;
        }

        aic_sum / n as f64
    }
}

/// Compute BIC for basis fit
/// BIC = n * log(RSS/n) + log(n) * total_edf
/// Where total_edf = n_curves * edf (each curve has edf parameters)
/// When pooled=true: compute single BIC across all curves
/// When pooled=false: compute per-curve BIC and return mean
#[extendr]
fn basis_bic_1d(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    nbasis: i32,
    basis_type: i32,
    lambda: f64,
    pooled: bool,
) -> f64 {
    let n = data.nrows();
    let m = data.ncols();
    let nbasis = nbasis as usize;

    if n == 0 || m == 0 || nbasis < 2 {
        return f64::INFINITY;
    }

    let data_slice = data.as_real_slice().unwrap();

    let basis = if basis_type == 1 {
        fourier_basis(&argvals, nbasis)
    } else {
        // For order 4 B-splines: nbasis = nknots + order, so nknots = nbasis - 4
        bspline_basis(&argvals, (nbasis.saturating_sub(4)).max(2), 4)
    };

    let actual_nbasis = basis.len() / m;
    let b_mat = DMatrix::from_column_slice(m, actual_nbasis, &basis);

    let btb = &b_mat.transpose() * &b_mat;
    let btb_penalized = if lambda > 0.0 {
        let d = difference_matrix(actual_nbasis, 2);
        let penalty = &d.transpose() * &d;
        btb + lambda * penalty
    } else {
        btb
    };

    // Solve for coefficients using SVD pseudoinverse (handles singular matrices)
    let svd = SVD::new(btb_penalized.clone(), true, true);
    let eps = 1e-10 * svd.singular_values.iter().cloned().fold(0.0_f64, f64::max);

    let u = match svd.u.as_ref() {
        Some(u) => u,
        None => return f64::INFINITY,
    };
    let v_t = match svd.v_t.as_ref() {
        Some(v_t) => v_t,
        None => return f64::INFINITY,
    };

    let s_inv: Vec<f64> = svd
        .singular_values
        .iter()
        .map(|&s| if s > eps { 1.0 / s } else { 0.0 })
        .collect();

    let mut btb_inv = DMatrix::zeros(actual_nbasis, actual_nbasis);
    for i in 0..actual_nbasis {
        for j in 0..actual_nbasis {
            let mut sum = 0.0;
            for k in 0..actual_nbasis.min(s_inv.len()) {
                sum += v_t[(k, i)] * s_inv[k] * u[(j, k)];
            }
            btb_inv[(i, j)] = sum;
        }
    }

    // EDF (effective degrees of freedom per curve)
    let proj = &btb_inv * b_mat.transpose();
    let h_mat = &b_mat * &proj;
    let edf: f64 = (0..m).map(|i| h_mat[(i, i)]).sum();

    if pooled {
        // Pooled mode: compute single BIC across all curves
        let mut total_rss = 0.0;
        let total_points = n * m;

        for i in 0..n {
            let curve: Vec<f64> = (0..m).map(|j| data_slice[i + j * n]).collect();
            let curve_vec = DVector::from_vec(curve.clone());
            let bt_y = b_mat.transpose() * &curve_vec;
            let coefs = &btb_inv * bt_y;
            let fitted = &b_mat * &coefs;

            for j in 0..m {
                let resid = curve[j] - fitted[j];
                total_rss += resid * resid;
            }
        }

        // BIC = n_total * log(RSS/n_total) + log(n_total) * total_edf
        // total_edf = n_curves * edf (each curve has edf effective parameters)
        let mse = total_rss / total_points as f64;
        let total_edf = (n as f64) * edf;
        (total_points as f64) * mse.ln() + (total_points as f64).ln() * total_edf
    } else {
        // Per-curve mode: compute BIC for each curve, return mean
        let mut bic_sum = 0.0;

        for i in 0..n {
            let curve: Vec<f64> = (0..m).map(|j| data_slice[i + j * n]).collect();
            let curve_vec = DVector::from_vec(curve.clone());
            let bt_y = b_mat.transpose() * &curve_vec;
            let coefs = &btb_inv * bt_y;
            let fitted = &b_mat * &coefs;

            let mut curve_rss = 0.0;
            for j in 0..m {
                let resid = curve[j] - fitted[j];
                curve_rss += resid * resid;
            }

            // BIC for this curve: m * log(RSS/m) + log(m) * edf
            let curve_mse = curve_rss / m as f64;
            bic_sum += (m as f64) * curve_mse.ln() + (m as f64).ln() * edf;
        }

        bic_sum / n as f64
    }
}

/// P-spline fitting: returns coefficients, fitted values, and diagnostics
#[extendr]
fn pspline_fit_1d(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    nbasis: i32,
    lambda: f64,
    order: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let nbasis = nbasis as usize;
    let order = order as usize;

    if n == 0 || m == 0 || nbasis < 2 {
        return list!(
            coefs = r!(RMatrix::new_matrix(0, 0, |_, _| 0.0)),
            fitted = r!(RMatrix::new_matrix(0, 0, |_, _| 0.0)),
            edf = 0.0,
            rss = 0.0,
            gcv = f64::INFINITY,
            aic = f64::INFINITY,
            bic = f64::INFINITY
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();

    // B-spline basis for P-splines
    // For order 4 B-splines: nbasis = nknots + order, so nknots = nbasis - 4
    let basis = bspline_basis(&argvals, (nbasis.saturating_sub(4)).max(2), 4);
    let actual_nbasis = basis.len() / m;
    let b_mat = DMatrix::from_column_slice(m, actual_nbasis, &basis);

    // Penalty matrix
    let d = difference_matrix(actual_nbasis, order);
    let penalty = &d.transpose() * &d;

    // B'B + lambda * D'D
    let btb = &b_mat.transpose() * &b_mat;
    let btb_penalized = &btb + lambda * &penalty;

    // Use SVD pseudoinverse for robustness with singular matrices
    let svd = SVD::new(btb_penalized.clone(), true, true);
    let eps = 1e-10 * svd.singular_values.iter().cloned().fold(0.0_f64, f64::max);

    let (u, v_t) = match (svd.u.as_ref(), svd.v_t.as_ref()) {
        (Some(u), Some(v_t)) => (u, v_t),
        _ => {
            return list!(
                coefs = r!(RMatrix::new_matrix(0, 0, |_, _| 0.0)),
                fitted = r!(RMatrix::new_matrix(0, 0, |_, _| 0.0)),
                edf = 0.0,
                rss = 0.0,
                gcv = f64::INFINITY,
                aic = f64::INFINITY,
                bic = f64::INFINITY
            )
            .into();
        }
    };

    let s_inv: Vec<f64> = svd
        .singular_values
        .iter()
        .map(|&s| if s > eps { 1.0 / s } else { 0.0 })
        .collect();

    let mut btb_inv = DMatrix::zeros(actual_nbasis, actual_nbasis);
    for i in 0..actual_nbasis {
        for j in 0..actual_nbasis {
            let mut sum = 0.0;
            for k in 0..actual_nbasis.min(s_inv.len()) {
                sum += v_t[(k, i)] * s_inv[k] * u[(j, k)];
            }
            btb_inv[(i, j)] = sum;
        }
    }

    // EDF = trace(H) where H = B * (B'B + lambda*P)^-1 * B'
    let proj = &btb_inv * b_mat.transpose();
    let h_mat = &b_mat * &proj;
    let edf: f64 = (0..m).map(|i| h_mat[(i, i)]).sum();

    // Fit each curve
    let mut all_coefs = vec![0.0; n * actual_nbasis];
    let mut all_fitted = vec![0.0; n * m];
    let mut total_rss = 0.0;

    for i in 0..n {
        let curve: Vec<f64> = (0..m).map(|j| data_slice[i + j * n]).collect();
        let curve_vec = DVector::from_vec(curve.clone());

        // Coefficients
        let bt_y = b_mat.transpose() * &curve_vec;
        let coefs = &btb_inv * bt_y;

        // Store coefficients
        for k in 0..actual_nbasis {
            all_coefs[i + k * n] = coefs[k];
        }

        // Fitted values
        let fitted = &b_mat * &coefs;
        for j in 0..m {
            all_fitted[i + j * n] = fitted[j];
            let resid = curve[j] - fitted[j];
            total_rss += resid * resid;
        }
    }

    let total_points = (n * m) as f64;

    // GCV
    let gcv_denom = 1.0 - edf / m as f64;
    let gcv = if gcv_denom.abs() > 1e-10 {
        (total_rss / total_points) / (gcv_denom * gcv_denom)
    } else {
        f64::INFINITY
    };

    // AIC and BIC
    let mse = total_rss / total_points;
    let aic = total_points * mse.ln() + 2.0 * edf;
    let bic = total_points * mse.ln() + total_points.ln() * edf;

    let coefs_mat = RMatrix::new_matrix(n, actual_nbasis, |i, k| all_coefs[i + k * n]);
    let fitted_mat = RMatrix::new_matrix(n, m, |i, j| all_fitted[i + j * n]);

    list!(
        coefs = r!(coefs_mat),
        fitted = r!(fitted_mat),
        edf = edf,
        rss = total_rss,
        gcv = gcv,
        aic = aic,
        bic = bic
    )
    .into()
}

// =============================================================================
// 2D Tensor Product Basis Functions
// =============================================================================

/// Compute tensor product basis matrix for 2D functional data
fn tensor_product_basis(
    argvals_s: &[f64],
    argvals_t: &[f64],
    nbasis_s: usize,
    nbasis_t: usize,
    basis_type: i32,
) -> Vec<f64> {
    let m1 = argvals_s.len();
    let m2 = argvals_t.len();

    // Compute 1D bases
    let b_s = if basis_type == 1 {
        fourier_basis(argvals_s, nbasis_s)
    } else {
        bspline_basis(argvals_s, nbasis_s + 2, 4)
    };

    let b_t = if basis_type == 1 {
        fourier_basis(argvals_t, nbasis_t)
    } else {
        bspline_basis(argvals_t, nbasis_t + 2, 4)
    };

    let actual_nbasis_s = b_s.len() / m1;
    let actual_nbasis_t = b_t.len() / m2;
    let total_basis = actual_nbasis_s * actual_nbasis_t;

    // Tensor product: B_2d[(i*m2+j), (k*nbasis_t+l)] = B_s[i,k] * B_t[j,l]
    let mut b_2d = vec![0.0; m1 * m2 * total_basis];

    for i in 0..m1 {
        for j in 0..m2 {
            let row_idx = i * m2 + j;
            for k in 0..actual_nbasis_s {
                for l in 0..actual_nbasis_t {
                    let col_idx = k * actual_nbasis_t + l;
                    b_2d[row_idx + col_idx * (m1 * m2)] = b_s[i + k * m1] * b_t[j + l * m2];
                }
            }
        }
    }

    b_2d
}

/// Project 2D functional data to tensor product basis coefficients
#[extendr]
fn fdata2basis_2d(
    data: RMatrix<f64>,
    argvals_s: Vec<f64>,
    argvals_t: Vec<f64>,
    nbasis_s: i32,
    nbasis_t: i32,
    basis_type: i32,
) -> Robj {
    let n = data.nrows();
    let m_total = data.ncols();
    let m1 = argvals_s.len();
    let m2 = argvals_t.len();

    if n == 0 || m_total != m1 * m2 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let nbasis_s = nbasis_s as usize;
    let nbasis_t = nbasis_t as usize;
    let data_slice = data.as_real_slice().unwrap();

    // Compute tensor product basis
    let basis = tensor_product_basis(&argvals_s, &argvals_t, nbasis_s, nbasis_t, basis_type);

    // Determine actual basis sizes
    let b_s_test = if basis_type == 1 {
        fourier_basis(&argvals_s, nbasis_s)
    } else {
        bspline_basis(&argvals_s, nbasis_s + 2, 4)
    };
    let actual_nbasis_s = b_s_test.len() / m1;

    let b_t_test = if basis_type == 1 {
        fourier_basis(&argvals_t, nbasis_t)
    } else {
        bspline_basis(&argvals_t, nbasis_t + 2, 4)
    };
    let actual_nbasis_t = b_t_test.len() / m2;

    let total_basis = actual_nbasis_s * actual_nbasis_t;
    let total_points = m1 * m2;

    let b_mat = DMatrix::from_column_slice(total_points, total_basis, &basis);

    // Least squares projection
    let btb = &b_mat.transpose() * &b_mat;
    let btb_svd = SVD::new(btb, true, true);

    let max_sv = btb_svd.singular_values.iter().cloned().fold(0.0, f64::max);
    let eps = 1e-10 * max_sv;

    let s_inv: Vec<f64> = btb_svd
        .singular_values
        .iter()
        .map(|&s| if s > eps { 1.0 / s } else { 0.0 })
        .collect();

    let v = btb_svd.v_t.as_ref().unwrap().transpose();
    let u_t = btb_svd.u.as_ref().unwrap().transpose();

    let mut btb_inv = DMatrix::zeros(total_basis, total_basis);
    for i in 0..total_basis {
        for j in 0..total_basis {
            let mut sum = 0.0;
            for k in 0..total_basis.min(s_inv.len()) {
                sum += v[(i, k)] * s_inv[k] * u_t[(k, j)];
            }
            btb_inv[(i, j)] = sum;
        }
    }

    let proj = &btb_inv * b_mat.transpose();

    // Compute coefficients for each surface
    let coefs: Vec<f64> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            let surface: Vec<f64> = (0..total_points).map(|j| data_slice[i + j * n]).collect();
            (0..total_basis)
                .map(|k| {
                    let mut sum = 0.0;
                    for j in 0..total_points {
                        sum += proj[(k, j)] * surface[j];
                    }
                    sum
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let result = RMatrix::new_matrix(n, total_basis, |i, k| coefs[i * total_basis + k]);
    r!(result)
}

/// Reconstruct 2D functional data from tensor product basis coefficients
#[extendr]
fn basis2fdata_2d(
    coefs: RMatrix<f64>,
    argvals_s: Vec<f64>,
    argvals_t: Vec<f64>,
    nbasis_s: i32,
    nbasis_t: i32,
    basis_type: i32,
) -> Robj {
    let n = coefs.nrows();
    let m1 = argvals_s.len();
    let m2 = argvals_t.len();
    let total_points = m1 * m2;

    let nbasis_s = nbasis_s as usize;
    let nbasis_t = nbasis_t as usize;
    let coefs_slice = coefs.as_real_slice().unwrap();
    let coefs_ncols = coefs.ncols();

    // Compute tensor product basis
    let basis = tensor_product_basis(&argvals_s, &argvals_t, nbasis_s, nbasis_t, basis_type);
    let total_basis = basis.len() / total_points;
    let k_max = total_basis.min(coefs_ncols);

    // Reconstruct: data = coefs * B^T
    let data: Vec<f64> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            (0..total_points)
                .map(|j| {
                    let mut sum = 0.0;
                    for k in 0..k_max {
                        sum += coefs_slice[i + k * n] * basis[j + k * total_points];
                    }
                    sum
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let result = RMatrix::new_matrix(n, total_points, |i, j| data[i * total_points + j]);
    r!(result)
}

/// 2D P-spline fitting with anisotropic penalties
#[extendr]
#[allow(clippy::too_many_arguments)]
fn pspline_fit_2d(
    data: RMatrix<f64>,
    argvals_s: Vec<f64>,
    argvals_t: Vec<f64>,
    nbasis_s: i32,
    nbasis_t: i32,
    lambda_s: f64,
    lambda_t: f64,
    order: i32,
) -> Robj {
    let n = data.nrows();
    let m1 = argvals_s.len();
    let m2 = argvals_t.len();
    let total_points = m1 * m2;

    if n == 0 || data.ncols() != total_points {
        return list!(
            coefs = r!(RMatrix::new_matrix(0, 0, |_, _| 0.0)),
            fitted = r!(RMatrix::new_matrix(0, 0, |_, _| 0.0)),
            edf = 0.0,
            gcv = f64::INFINITY,
            aic = f64::INFINITY,
            bic = f64::INFINITY
        )
        .into();
    }

    let nbasis_s = nbasis_s as usize;
    let nbasis_t = nbasis_t as usize;
    let order = order as usize;
    let data_slice = data.as_real_slice().unwrap();

    // Tensor product B-spline basis
    let basis = tensor_product_basis(&argvals_s, &argvals_t, nbasis_s, nbasis_t, 0);

    // Get actual basis dimensions
    let b_s_test = bspline_basis(&argvals_s, nbasis_s + 2, 4);
    let actual_nbasis_s = b_s_test.len() / m1;

    let b_t_test = bspline_basis(&argvals_t, nbasis_t + 2, 4);
    let actual_nbasis_t = b_t_test.len() / m2;

    let total_basis = actual_nbasis_s * actual_nbasis_t;
    let b_mat = DMatrix::from_column_slice(total_points, total_basis, &basis);

    // Create 1D penalty matrices
    let d_s = difference_matrix(actual_nbasis_s, order);
    let p_s = &d_s.transpose() * &d_s;

    let d_t = difference_matrix(actual_nbasis_t, order);
    let p_t = &d_t.transpose() * &d_t;

    // Kronecker product penalty: P = lambda_s*(I_t ⊗ P_s) + lambda_t*(P_t ⊗ I_s)
    // I_t ⊗ P_s
    let mut penalty_s = DMatrix::zeros(total_basis, total_basis);
    for it in 0..actual_nbasis_t {
        for jt in 0..actual_nbasis_t {
            if it == jt {
                for is in 0..actual_nbasis_s {
                    for js in 0..actual_nbasis_s {
                        let row = is * actual_nbasis_t + it;
                        let col = js * actual_nbasis_t + jt;
                        penalty_s[(row, col)] = p_s[(is, js)];
                    }
                }
            }
        }
    }

    // P_t ⊗ I_s
    let mut penalty_t = DMatrix::zeros(total_basis, total_basis);
    for is in 0..actual_nbasis_s {
        for js in 0..actual_nbasis_s {
            if is == js {
                for it in 0..actual_nbasis_t {
                    for jt in 0..actual_nbasis_t {
                        let row = is * actual_nbasis_t + it;
                        let col = js * actual_nbasis_t + jt;
                        penalty_t[(row, col)] = p_t[(it, jt)];
                    }
                }
            }
        }
    }

    let penalty = lambda_s * penalty_s + lambda_t * penalty_t;

    let btb = &b_mat.transpose() * &b_mat;
    let btb_penalized = &btb + &penalty;

    let btb_inv = match btb_penalized.clone().try_inverse() {
        Some(inv) => inv,
        None => {
            return list!(
                coefs = r!(RMatrix::new_matrix(0, 0, |_, _| 0.0)),
                fitted = r!(RMatrix::new_matrix(0, 0, |_, _| 0.0)),
                edf = 0.0,
                gcv = f64::INFINITY,
                aic = f64::INFINITY,
                bic = f64::INFINITY
            )
            .into();
        }
    };

    // EDF
    let proj = &btb_inv * b_mat.transpose();
    let h_mat = &b_mat * &proj;
    let edf: f64 = (0..total_points).map(|i| h_mat[(i, i)]).sum();

    // Fit each surface
    let mut all_coefs = vec![0.0; n * total_basis];
    let mut all_fitted = vec![0.0; n * total_points];
    let mut total_rss = 0.0;

    for i in 0..n {
        let surface: Vec<f64> = (0..total_points).map(|j| data_slice[i + j * n]).collect();
        let surface_vec = DVector::from_vec(surface.clone());

        let bt_y = b_mat.transpose() * &surface_vec;
        let coefs = &btb_inv * bt_y;

        for k in 0..total_basis {
            all_coefs[i + k * n] = coefs[k];
        }

        let fitted = &b_mat * &coefs;
        for j in 0..total_points {
            all_fitted[i + j * n] = fitted[j];
            let resid = surface[j] - fitted[j];
            total_rss += resid * resid;
        }
    }

    let total_obs = (n * total_points) as f64;

    let gcv_denom = 1.0 - edf / total_points as f64;
    let gcv = if gcv_denom.abs() > 1e-10 {
        (total_rss / total_obs) / (gcv_denom * gcv_denom)
    } else {
        f64::INFINITY
    };

    let mse = total_rss / total_obs;
    let aic = total_obs * mse.ln() + 2.0 * edf;
    let bic = total_obs * mse.ln() + total_obs.ln() * edf;

    let coefs_mat = RMatrix::new_matrix(n, total_basis, |i, k| all_coefs[i + k * n]);
    let fitted_mat = RMatrix::new_matrix(n, total_points, |i, j| all_fitted[i + j * n]);

    list!(
        coefs = r!(coefs_mat),
        fitted = r!(fitted_mat),
        edf = edf,
        rss = total_rss,
        gcv = gcv,
        aic = aic,
        bic = bic
    )
    .into()
}

/// Automatic basis selection for each curve individually.
///
/// This function compares Fourier and P-spline bases for each curve,
/// selecting the optimal basis type and number of basis functions using
/// model selection criteria (GCV, AIC, or BIC).
///
/// # Arguments
/// * `data` - Functional data matrix (n x m)
/// * `argvals` - Evaluation points
/// * `criterion` - Model selection criterion: 0=GCV (default), 1=AIC, 2=BIC
/// * `nbasis_min` - Minimum number of basis functions (0 for auto)
/// * `nbasis_max` - Maximum number of basis functions (0 for auto)
/// * `lambda_pspline` - Smoothing parameter for P-spline (negative for auto-select)
/// * `use_seasonal_hint` - Whether to use FFT to detect seasonality
///
/// # Returns
/// Named list with:
/// - basis_type: Integer vector (0=pspline, 1=fourier)
/// - nbasis: Integer vector of selected nbasis per curve
/// - score: Numeric vector of criterion scores
/// - coefficients: List of coefficient vectors
/// - fitted: Fitted values matrix (n x m)
/// - edf: Numeric vector of effective degrees of freedom
/// - seasonal_detected: Logical vector
/// - lambda: Numeric vector of lambda values (NA for Fourier)
/// - criterion: Character (criterion name used)
#[extendr]
fn select_basis_auto(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    criterion: i32,
    nbasis_min: i32,
    nbasis_max: i32,
    lambda_pspline: f64,
    use_seasonal_hint: bool,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m == 0 || argvals.len() != m {
        return list!(
            basis_type = r!(Vec::<i32>::new()),
            nbasis = r!(Vec::<i32>::new()),
            score = r!(Vec::<f64>::new()),
            coefficients = r!(extendr_api::List::new(0)),
            fitted = r!(RMatrix::new_matrix(0, 0, |_, _| 0.0)),
            edf = r!(Vec::<f64>::new()),
            seasonal_detected = r!(Vec::<bool>::new()),
            lambda = r!(Vec::<f64>::new()),
            criterion_name = "GCV"
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let result = fdars_core::basis::select_basis_auto_1d(
        &fd,
        &argvals,
        criterion,
        nbasis_min as usize,
        nbasis_max as usize,
        lambda_pspline,
        use_seasonal_hint,
    );

    // Extract results into R vectors
    let basis_types: Vec<i32> = result.selections.iter().map(|s| s.basis_type.to_i32()).collect();
    let nbasis_vec: Vec<i32> = result.selections.iter().map(|s| s.nbasis as i32).collect();
    let scores: Vec<f64> = result.selections.iter().map(|s| s.score).collect();
    let edfs: Vec<f64> = result.selections.iter().map(|s| s.edf).collect();
    let seasonal_detected: Vec<bool> = result
        .selections
        .iter()
        .map(|s| s.seasonal_detected)
        .collect();
    let lambdas: Vec<f64> = result.selections.iter().map(|s| s.lambda).collect();

    // Build coefficients list
    let coef_list: Vec<Robj> = result
        .selections
        .iter()
        .map(|s| r!(s.coefficients.clone()))
        .collect();

    // Build fitted matrix
    let mut fitted_data = vec![0.0; n * m];
    for (i, sel) in result.selections.iter().enumerate() {
        for (j, &val) in sel.fitted.iter().enumerate() {
            if j < m {
                fitted_data[i + j * n] = val;
            }
        }
    }
    let fitted_mat = RMatrix::new_matrix(n, m, |i, j| fitted_data[i + j * n]);

    let criterion_name = match result.criterion {
        0 => "GCV",
        1 => "AIC",
        _ => "BIC",
    };

    list!(
        basis_type = r!(basis_types),
        nbasis = r!(nbasis_vec),
        score = r!(scores),
        coefficients = extendr_api::List::from_values(coef_list),
        fitted = r!(fitted_mat),
        edf = r!(edfs),
        seasonal_detected = r!(seasonal_detected),
        lambda = r!(lambdas),
        criterion_name = criterion_name
    )
    .into()
}

// =============================================================================
// Outlier detection (LRT-based)
// =============================================================================

/// Compute bootstrap threshold for LRT outlier detection
/// Highly parallelized across bootstrap iterations
#[extendr]
fn outliers_thres_lrt(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    nb: i32,
    smo: f64,
    trim: f64,
    seed: u64,
    percentile: f64,
) -> f64 {
    let n = data.nrows();
    let m = data.ncols();

    if n < 3 || m == 0 || argvals.len() != m {
        return 0.0;
    }

    let data_slice = data.as_real_slice().unwrap();
    let nb = nb as usize;
    let n_keep = ((1.0 - trim) * n as f64).ceil() as usize;

    // First, compute FM depth on the original data and identify the trimmed
    // (deep) set. Bootstrap only from these curves to avoid the masking effect
    // where outliers inflate the bootstrap distribution threshold.
    let orig_depths = compute_fm_depth_internal(data_slice, n, m);
    let mut orig_depth_idx: Vec<(usize, f64)> = orig_depths
        .iter()
        .enumerate()
        .map(|(i, &d)| (i, d))
        .collect();
    orig_depth_idx
        .sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let clean_idx: Vec<usize> = orig_depth_idx
        .iter()
        .take(n_keep)
        .map(|(i, _)| *i)
        .collect();
    let n_clean = clean_idx.len();

    // Compute column standard deviations from the clean set only
    let col_vars: Vec<f64> = (0..m)
        .map(|j| {
            let mut sum = 0.0;
            let mut sum_sq = 0.0;
            for &i in &clean_idx {
                let val = data_slice[i + j * n];
                sum += val;
                sum_sq += val * val;
            }
            let mean = sum / n_clean as f64;
            let var = sum_sq / n_clean as f64 - mean * mean;
            var.max(0.0).sqrt()
        })
        .collect();

    // Run bootstrap iterations in parallel
    // Resample from the clean (trimmed) set to compute the null distribution
    let max_dists: Vec<f64> = (0..nb)
        .into_par_iter()
        .map(|b| {
            let mut rng = StdRng::seed_from_u64(seed + b as u64);

            // Resample with replacement from the clean set only
            let indices: Vec<usize> =
                (0..n_clean).map(|_| clean_idx[rng.gen_range(0..n_clean)]).collect();

            // Build resampled data with smoothing noise
            let mut boot_data = vec![0.0; n_clean * m];
            for (new_i, &old_i) in indices.iter().enumerate() {
                for j in 0..m {
                    let noise: f64 = rng.sample::<f64, _>(StandardNormal) * smo * col_vars[j];
                    boot_data[new_i + j * n_clean] = data_slice[old_i + j * n] + noise;
                }
            }

            // Compute FM depth for bootstrap sample
            let depths = compute_fm_depth_internal(&boot_data, n_clean, m);

            // Get indices of top n_keep_boot curves by depth
            let n_keep_boot = ((1.0 - trim) * n_clean as f64).ceil() as usize;
            let mut depth_idx: Vec<(usize, f64)> =
                depths.iter().enumerate().map(|(i, &d)| (i, d)).collect();
            depth_idx.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            let keep_idx: Vec<usize> = depth_idx.iter().take(n_keep_boot).map(|(i, _)| *i).collect();

            // Compute trimmed mean and variance of kept curves
            let mut trimmed_mean = vec![0.0; m];
            for j in 0..m {
                for &i in &keep_idx {
                    trimmed_mean[j] += boot_data[i + j * n_clean];
                }
                trimmed_mean[j] /= n_keep_boot as f64;
            }

            let mut trimmed_var = vec![0.0; m];
            for j in 0..m {
                for &i in &keep_idx {
                    let diff = boot_data[i + j * n_clean] - trimmed_mean[j];
                    trimmed_var[j] += diff * diff;
                }
                trimmed_var[j] /= n_keep_boot as f64;
                trimmed_var[j] = trimmed_var[j].max(1e-10); // Avoid division by zero
            }

            // Compute max normalized distance to trimmed mean
            let mut max_dist = 0.0;
            for i in 0..n_clean {
                let mut dist = 0.0;
                for j in 0..m {
                    let diff = boot_data[i + j * n_clean] - trimmed_mean[j];
                    dist += diff * diff / trimmed_var[j];
                }
                dist = (dist / m as f64).sqrt();
                if dist > max_dist {
                    max_dist = dist;
                }
            }

            max_dist
        })
        .collect();

    // Return threshold at specified percentile
    let mut sorted_dists = max_dists;
    sorted_dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let idx = ((nb as f64 * percentile) as usize).min(nb.saturating_sub(1));
    sorted_dists.get(idx).copied().unwrap_or(0.0)
}

/// Internal helper: compute threshold AND bootstrap distribution
fn outliers_thres_lrt_with_dist(
    data_slice: &[f64],
    n: usize,
    m: usize,
    nb: usize,
    smo: f64,
    trim: f64,
    seed: u64,
    percentile: f64,
) -> (f64, Vec<f64>) {
    if n < 3 || m == 0 {
        return (0.0, vec![]);
    }
    let n_keep = ((1.0 - trim) * n as f64).ceil().max(1.0) as usize;
    let n_keep = n_keep.min(n);

    let orig_depths = compute_fm_depth_internal(data_slice, n, m);
    let mut orig_depth_idx: Vec<(usize, f64)> = orig_depths
        .iter()
        .enumerate()
        .map(|(i, &d)| (i, d))
        .collect();
    orig_depth_idx
        .sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let clean_idx: Vec<usize> = orig_depth_idx
        .iter()
        .take(n_keep)
        .map(|(i, _)| *i)
        .collect();
    let n_clean = clean_idx.len();

    let col_vars: Vec<f64> = (0..m)
        .map(|j| {
            let mut sum = 0.0;
            let mut sum_sq = 0.0;
            for &i in &clean_idx {
                let val = data_slice[i + j * n];
                sum += val;
                sum_sq += val * val;
            }
            let mean = sum / n_clean as f64;
            let var = sum_sq / n_clean as f64 - mean * mean;
            var.max(0.0).sqrt()
        })
        .collect();

    let max_dists: Vec<f64> = (0..nb)
        .into_par_iter()
        .map(|b| {
            let mut rng = StdRng::seed_from_u64(seed + b as u64);
            let indices: Vec<usize> =
                (0..n_clean).map(|_| clean_idx[rng.gen_range(0..n_clean)]).collect();
            let mut boot_data = vec![0.0; n_clean * m];
            for (new_i, &old_i) in indices.iter().enumerate() {
                for j in 0..m {
                    let noise: f64 = rng.sample::<f64, _>(StandardNormal) * smo * col_vars[j];
                    boot_data[new_i + j * n_clean] = data_slice[old_i + j * n] + noise;
                }
            }
            let depths = compute_fm_depth_internal(&boot_data, n_clean, m);
            let n_keep_boot = ((1.0 - trim) * n_clean as f64).ceil().max(1.0) as usize;
            let n_keep_boot = n_keep_boot.min(n_clean);
            let mut depth_idx: Vec<(usize, f64)> =
                depths.iter().enumerate().map(|(i, &d)| (i, d)).collect();
            depth_idx.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            let keep_idx: Vec<usize> = depth_idx.iter().take(n_keep_boot).map(|(i, _)| *i).collect();
            let mut trimmed_mean = vec![0.0; m];
            for j in 0..m {
                for &i in &keep_idx {
                    trimmed_mean[j] += boot_data[i + j * n_clean];
                }
                trimmed_mean[j] /= n_keep_boot as f64;
            }
            let mut trimmed_var = vec![0.0; m];
            for j in 0..m {
                for &i in &keep_idx {
                    let diff = boot_data[i + j * n_clean] - trimmed_mean[j];
                    trimmed_var[j] += diff * diff;
                }
                trimmed_var[j] /= n_keep_boot as f64;
                trimmed_var[j] = trimmed_var[j].max(1e-10);
            }
            let mut max_dist = 0.0;
            for i in 0..n_clean {
                let mut dist = 0.0;
                for j in 0..m {
                    let diff = boot_data[i + j * n_clean] - trimmed_mean[j];
                    dist += diff * diff / trimmed_var[j];
                }
                dist = (dist / m as f64).sqrt();
                if dist > max_dist {
                    max_dist = dist;
                }
            }
            max_dist
        })
        .collect();

    let mut sorted_dists = max_dists;
    sorted_dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let idx = ((nb as f64 * percentile) as usize).min(nb.saturating_sub(1));
    let threshold = sorted_dists.get(idx).copied().unwrap_or(0.0);
    (threshold, sorted_dists)
}

/// Helper to compute FM depth internally
fn compute_fm_depth_internal(data: &[f64], n: usize, m: usize) -> Vec<f64> {
    // Compute FM depth for each curve
    (0..n)
        .into_par_iter()
        .map(|i| {
            let mut depth_sum = 0.0;
            for j in 0..m {
                let val = data[i + j * n];
                let mut n_below = 0;
                let mut n_equal = 0;
                for k in 0..n {
                    let other_val = data[k + j * n];
                    if other_val < val {
                        n_below += 1;
                    } else if (other_val - val).abs() < 1e-10 {
                        n_equal += 1;
                    }
                }
                // Univariate depth at this point
                let prop_below = (n_below as f64 + 0.5 * n_equal as f64) / n as f64;
                let univ_depth = prop_below.min(1.0 - prop_below);
                depth_sum += univ_depth;
            }
            // Average over all evaluation points, scale by 2 for FM1 formula
            2.0 * depth_sum / m as f64
        })
        .collect()
}

/// LRT-based outlier detection
/// Returns indices of detected outliers
#[extendr]
fn outliers_lrt(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    nb: i32,
    smo: f64,
    trim: f64,
    seed: u64,
    percentile: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n < 3 || m == 0 || argvals.len() != m {
        return list!(
            outliers = Vec::<i32>::new(),
            distances = Vec::<f64>::new(),
            threshold = 0.0
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap().to_vec();
    let _n_keep = ((1.0 - trim) * n as f64).ceil() as usize;

    // Compute threshold and bootstrap null distribution
    let (threshold, boot_dist) = outliers_thres_lrt_with_dist(&data_slice, n, m, nb as usize, smo, trim, seed, percentile);

    // Iterative outlier detection (up to 5 iterations)
    let mut current_mask = vec![true; n];
    let mut detected_outliers = Vec::new();

    for _iter in 0..5 {
        let current_n: usize = current_mask.iter().filter(|&&x| x).count();
        if current_n < 3 {
            break;
        }

        // Get current active indices
        let active_idx: Vec<usize> = current_mask
            .iter()
            .enumerate()
            .filter(|(_, &active)| active)
            .map(|(i, _)| i)
            .collect();

        // Compute FM depth on current active set
        let active_n = active_idx.len();
        let mut active_data: Vec<f64> = Vec::with_capacity(active_n * m);
        for j in 0..m {
            for &i in &active_idx {
                active_data.push(data_slice[i + j * n]);
            }
        }

        let depths = compute_fm_depth_internal(&active_data, active_n, m);

        // Get indices of top n_keep_cur curves by depth
        let n_keep_cur = ((1.0 - trim) * active_n as f64).ceil() as usize;
        let mut depth_idx: Vec<(usize, f64)> =
            depths.iter().enumerate().map(|(i, &d)| (i, d)).collect();
        depth_idx.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let keep_idx: Vec<usize> = depth_idx.iter().take(n_keep_cur).map(|(i, _)| *i).collect();

        // Compute trimmed mean and variance
        let mut trimmed_mean = vec![0.0; m];
        for j in 0..m {
            for &local_i in &keep_idx {
                trimmed_mean[j] += active_data[local_i + j * active_n];
            }
            trimmed_mean[j] /= n_keep_cur as f64;
        }

        let mut trimmed_var = vec![0.0; m];
        for j in 0..m {
            for &local_i in &keep_idx {
                let diff = active_data[local_i + j * active_n] - trimmed_mean[j];
                trimmed_var[j] += diff * diff;
            }
            trimmed_var[j] /= n_keep_cur as f64;
            trimmed_var[j] = trimmed_var[j].max(1e-10);
        }

        // Check each curve against threshold
        let mut new_outliers = Vec::new();
        for (local_i, &global_i) in active_idx.iter().enumerate() {
            let mut dist = 0.0;
            for j in 0..m {
                let diff = active_data[local_i + j * active_n] - trimmed_mean[j];
                dist += diff * diff / trimmed_var[j];
            }
            dist = (dist / m as f64).sqrt();

            if dist > threshold {
                new_outliers.push(global_i);
                current_mask[global_i] = false;
            }
        }

        if new_outliers.is_empty() {
            break;
        }

        detected_outliers.extend(new_outliers);
    }

    // Compute final distances for all curves
    let final_active_idx: Vec<usize> = current_mask
        .iter()
        .enumerate()
        .filter(|(_, &active)| active)
        .map(|(i, _)| i)
        .collect();

    let final_n = final_active_idx.len();
    if final_n < 3 {
        return list!(
            outliers = detected_outliers
                .iter()
                .map(|&i| (i + 1) as i32)
                .collect::<Vec<_>>(),
            distances = vec![0.0; n],
            threshold = threshold,
            boot_dist = boot_dist
        )
        .into();
    }

    let mut final_data: Vec<f64> = Vec::with_capacity(final_n * m);
    for j in 0..m {
        for &i in &final_active_idx {
            final_data.push(data_slice[i + j * n]);
        }
    }

    let final_depths = compute_fm_depth_internal(&final_data, final_n, m);

    let n_keep_final = ((1.0 - trim) * final_n as f64).ceil() as usize;
    let mut depth_idx: Vec<(usize, f64)> = final_depths
        .iter()
        .enumerate()
        .map(|(i, &d)| (i, d))
        .collect();
    depth_idx.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let keep_idx: Vec<usize> = depth_idx
        .iter()
        .take(n_keep_final)
        .map(|(i, _)| *i)
        .collect();

    let mut final_mean = vec![0.0; m];
    let mut final_var = vec![0.0; m];
    for j in 0..m {
        for &local_i in &keep_idx {
            final_mean[j] += final_data[local_i + j * final_n];
        }
        final_mean[j] /= n_keep_final as f64;
    }
    for j in 0..m {
        for &local_i in &keep_idx {
            let diff = final_data[local_i + j * final_n] - final_mean[j];
            final_var[j] += diff * diff;
        }
        final_var[j] /= n_keep_final as f64;
        final_var[j] = final_var[j].max(1e-10);
    }

    // Compute distances for all original curves
    let distances: Vec<f64> = (0..n)
        .map(|i| {
            let mut dist = 0.0;
            for j in 0..m {
                let diff = data_slice[i + j * n] - final_mean[j];
                dist += diff * diff / final_var[j];
            }
            (dist / m as f64).sqrt()
        })
        .collect();

    // Return 1-indexed outlier indices (for R compatibility)
    let outlier_indices: Vec<i32> = detected_outliers.iter().map(|&i| (i + 1) as i32).collect();

    list!(
        outliers = outlier_indices,
        distances = distances,
        threshold = threshold,
        boot_dist = boot_dist
    )
    .into()
}

// =============================================================================
// Smoothing Functions
// =============================================================================

/// Kernel function enum for smoothing
fn kernel_value(u: f64, kernel_type: &str) -> f64 {
    match kernel_type {
        "norm" | "normal" | "gaussian" => (-0.5 * u * u).exp() / (2.0 * PI).sqrt(),
        "epa" | "epanechnikov" => {
            if u.abs() <= 1.0 {
                0.75 * (1.0 - u * u)
            } else {
                0.0
            }
        }
        "tri" | "triweight" => {
            if u.abs() <= 1.0 {
                (35.0 / 32.0) * (1.0 - u * u).powi(3)
            } else {
                0.0
            }
        }
        "quar" | "quartic" | "biweight" => {
            if u.abs() <= 1.0 {
                (15.0 / 16.0) * (1.0 - u * u).powi(2)
            } else {
                0.0
            }
        }
        "cos" | "cosine" => {
            if u.abs() <= 1.0 {
                (PI / 4.0) * (PI * u / 2.0).cos()
            } else {
                0.0
            }
        }
        "unif" | "uniform" => {
            if u.abs() <= 1.0 {
                0.5
            } else {
                0.0
            }
        }
        _ => (-0.5 * u * u).exp() / (2.0 * PI).sqrt(), // default to normal
    }
}

/// Nadaraya-Watson smoother matrix
/// S_ij = K((t_i - t_j)/h) * w_j / sum_k(K((t_i - t_k)/h) * w_k)
#[extendr]
fn s_nw(argvals: Vec<f64>, h: f64, kernel_type: &str, weights: Vec<f64>, cv: bool) -> Robj {
    let n = argvals.len();
    if n == 0 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let w = if weights.len() == n {
        weights
    } else {
        vec![1.0; n]
    };

    // Compute smoother matrix (parallelized over rows)
    let s_mat: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut row = vec![0.0; n];
            let mut sum = 0.0;

            for j in 0..n {
                if cv && i == j {
                    // Leave-one-out: skip diagonal
                    row[j] = 0.0;
                } else {
                    let u = (argvals[i] - argvals[j]) / h;
                    let k = kernel_value(u, kernel_type) * w[j];
                    row[j] = k;
                    sum += k;
                }
            }

            // Normalize row
            if sum > 0.0 {
                for val in &mut row {
                    *val /= sum;
                }
            }

            row
        })
        .collect();

    // Flatten to column-major for R
    let result = RMatrix::new_matrix(n, n, |i, j| s_mat[i][j]);
    r!(result)
}

/// Local Linear Regression smoother matrix
/// Uses weighted least squares with degree-1 polynomial
#[extendr]
fn s_llr(argvals: Vec<f64>, h: f64, kernel_type: &str, weights: Vec<f64>, cv: bool) -> Robj {
    let n = argvals.len();
    if n < 2 {
        return r!(RMatrix::new_matrix(n, n, |i, j| if i == j {
            1.0
        } else {
            0.0
        }));
    }

    let w = if weights.len() == n {
        weights
    } else {
        vec![1.0; n]
    };

    // Compute smoother matrix (parallelized over rows)
    let s_mat: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let t0 = argvals[i];
            let mut row = vec![0.0; n];

            // Compute kernel weights and moments
            let mut s0 = 0.0;
            let mut s1 = 0.0;
            let mut s2 = 0.0;

            for j in 0..n {
                if cv && i == j {
                    continue;
                }
                let u = (argvals[j] - t0) / h;
                let k = kernel_value(u, kernel_type) * w[j];
                let delta = argvals[j] - t0;
                s0 += k;
                s1 += k * delta;
                s2 += k * delta * delta;
            }

            let denom = s0 * s2 - s1 * s1;

            if denom.abs() < 1e-15 {
                // Fall back to NW weights
                for j in 0..n {
                    if cv && i == j {
                        row[j] = 0.0;
                    } else {
                        let u = (argvals[j] - t0) / h;
                        row[j] = kernel_value(u, kernel_type) * w[j];
                    }
                }
                let sum: f64 = row.iter().sum();
                if sum > 0.0 {
                    for val in &mut row {
                        *val /= sum;
                    }
                }
            } else {
                // Local linear weights
                for j in 0..n {
                    if cv && i == j {
                        row[j] = 0.0;
                    } else {
                        let u = (argvals[j] - t0) / h;
                        let k = kernel_value(u, kernel_type) * w[j];
                        let delta = argvals[j] - t0;
                        row[j] = k * (s2 - s1 * delta) / denom;
                    }
                }
            }

            row
        })
        .collect();

    let result = RMatrix::new_matrix(n, n, |i, j| s_mat[i][j]);
    r!(result)
}

/// Local Polynomial Regression smoother matrix
/// Solves (p+1)×(p+1) weighted least squares system for each point
#[extendr]
fn s_lpr(
    argvals: Vec<f64>,
    h: f64,
    p: i32,
    kernel_type: &str,
    weights: Vec<f64>,
    cv: bool,
) -> Robj {
    let n = argvals.len();
    let degree = p as usize;

    if n <= degree {
        return r!(RMatrix::new_matrix(n, n, |i, j| if i == j {
            1.0
        } else {
            0.0
        }));
    }

    let w = if weights.len() == n {
        weights
    } else {
        vec![1.0; n]
    };

    // Compute smoother matrix (parallelized over rows)
    let s_mat: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let t0 = argvals[i];
            let mut row = vec![0.0; n];

            // Build design matrix and weight vector for local polynomial
            // X_j = [1, (t_j - t0), (t_j - t0)^2, ..., (t_j - t0)^p]
            // We want the first row of (X'WX)^{-1} X'W

            let p1 = degree + 1;

            // Compute X'WX matrix
            let mut xtw = vec![0.0; p1 * n]; // (p+1) x n
            let mut xtwx = vec![0.0; p1 * p1]; // (p+1) x (p+1)

            for j in 0..n {
                let skip = cv && i == j;
                if skip {
                    continue;
                }

                let u = (argvals[j] - t0) / h;
                let k = kernel_value(u, kernel_type) * w[j];
                let delta = argvals[j] - t0;

                // Compute powers of delta
                let mut powers = vec![1.0; p1];
                for d in 1..p1 {
                    powers[d] = powers[d - 1] * delta;
                }

                // Accumulate X'W and X'WX
                for d1 in 0..p1 {
                    xtw[d1 * n + j] = powers[d1] * k;
                    for d2 in 0..p1 {
                        xtwx[d1 * p1 + d2] += powers[d1] * powers[d2] * k;
                    }
                }
            }

            // Solve for first row of (X'WX)^{-1} X'W using nalgebra
            let xtwx_mat = DMatrix::from_row_slice(p1, p1, &xtwx);

            // We only need first row of inverse, i.e., e1' * (X'WX)^{-1}
            // Then multiply by X'W to get smoother weights
            match xtwx_mat.try_inverse() {
                Some(inv) => {
                    // First row of inverse times X'W
                    for j in 0..n {
                        let mut sum = 0.0;
                        for d in 0..p1 {
                            sum += inv[(0, d)] * xtw[d * n + j];
                        }
                        row[j] = sum;
                    }
                }
                None => {
                    // Fall back to NW weights
                    for j in 0..n {
                        if cv && i == j {
                            row[j] = 0.0;
                        } else {
                            let u = (argvals[j] - t0) / h;
                            row[j] = kernel_value(u, kernel_type) * w[j];
                        }
                    }
                    let sum: f64 = row.iter().sum();
                    if sum > 0.0 {
                        for val in &mut row {
                            *val /= sum;
                        }
                    }
                }
            }

            row
        })
        .collect();

    let result = RMatrix::new_matrix(n, n, |i, j| s_mat[i][j]);
    r!(result)
}

/// K-Nearest Neighbors smoother matrix
#[extendr]
fn s_knn(argvals: Vec<f64>, k: i32, kernel_type: &str, weights: Vec<f64>, cv: bool) -> Robj {
    let n = argvals.len();
    let knn = (k as usize).min(n);

    if n == 0 || knn == 0 {
        return r!(RMatrix::new_matrix(n, n, |i, j| if i == j {
            1.0
        } else {
            0.0
        }));
    }

    let w = if weights.len() == n {
        weights
    } else {
        vec![1.0; n]
    };

    // Compute smoother matrix (parallelized over rows)
    let s_mat: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let t0 = argvals[i];
            let mut row = vec![0.0; n];

            // Compute distances and find k-th nearest
            let mut distances: Vec<(usize, f64)> = argvals
                .iter()
                .enumerate()
                .filter(|(j, _)| !cv || *j != i)
                .map(|(j, &t)| (j, (t - t0).abs()))
                .collect();

            distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

            // Adaptive bandwidth: distance to k-th neighbor
            let h = if distances.len() >= knn {
                distances[knn - 1].1.max(1e-10)
            } else if !distances.is_empty() {
                distances.last().unwrap().1.max(1e-10)
            } else {
                1.0
            };

            // Apply kernel with adaptive bandwidth
            let mut sum = 0.0;
            for &(j, dist) in &distances {
                let u = dist / h;
                let k_val = kernel_value(u, kernel_type) * w[j];
                row[j] = k_val;
                sum += k_val;
            }

            // Normalize
            if sum > 0.0 {
                for val in &mut row {
                    *val /= sum;
                }
            }

            row
        })
        .collect();

    let result = RMatrix::new_matrix(n, n, |i, j| s_mat[i][j]);
    r!(result)
}

// =============================================================================
// Clustering Functions
// =============================================================================

/// Helper: compute distance matrix using specified metric
fn compute_distance_matrix(
    data: &[f64],
    n: usize,
    m: usize,
    argvals: &[f64],
    metric: &str,
) -> Vec<Vec<f64>> {
    let weights = simpsons_weights(argvals);

    (0..n)
        .into_par_iter()
        .map(|i| {
            let mut row = vec![0.0; n];
            for j in 0..n {
                if i == j {
                    row[j] = 0.0;
                } else {
                    let mut dist = 0.0;
                    match metric {
                        "L2" | "euclidean" => {
                            for k in 0..m {
                                let diff = data[i + k * n] - data[j + k * n];
                                dist += diff * diff * weights[k];
                            }
                            dist = dist.sqrt();
                        }
                        "L1" | "manhattan" => {
                            for k in 0..m {
                                let diff = (data[i + k * n] - data[j + k * n]).abs();
                                dist += diff * weights[k];
                            }
                        }
                        "Linf" | "supremum" => {
                            for k in 0..m {
                                let diff = (data[i + k * n] - data[j + k * n]).abs();
                                if diff > dist {
                                    dist = diff;
                                }
                            }
                        }
                        _ => {
                            // Default to L2
                            for k in 0..m {
                                let diff = data[i + k * n] - data[j + k * n];
                                dist += diff * diff * weights[k];
                            }
                            dist = dist.sqrt();
                        }
                    }
                    row[j] = dist;
                }
            }
            row
        })
        .collect()
}

/// Functional k-means clustering
#[extendr]
fn kmeans_fd(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    nclusters: i32,
    max_iter: i32,
    nstart: i32,
    metric: &str,
    seed: Nullable<i32>,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let k = nclusters as usize;

    if n == 0 || m == 0 || k == 0 || k > n {
        return r!(list!(
            cluster = Vec::<i32>::new(),
            centers = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            withinss = Vec::<f64>::new(),
            tot_withinss = 0.0,
            size = Vec::<i32>::new()
        ));
    }

    let data_slice = data.as_real_slice().unwrap();

    // Initialize RNG
    let mut rng = match seed {
        Nullable::NotNull(s) => rand::rngs::StdRng::seed_from_u64(s as u64),
        Nullable::Null => rand::rngs::StdRng::from_entropy(),
    };

    // Compute full distance matrix
    let dist_matrix = compute_distance_matrix(data_slice, n, m, &argvals, metric);

    // Run k-means nstart times and keep best result
    let mut best_cluster: Vec<i32> = vec![0; n];
    let mut best_centers: Vec<f64> = vec![0.0; k * m];
    let mut best_withinss = f64::MAX;
    let mut best_sizes: Vec<i32> = vec![0; k];

    for _ in 0..nstart.max(1) {
        // K-means++ initialization
        let mut centers_idx: Vec<usize> = Vec::with_capacity(k);

        // First center: random
        centers_idx.push(rng.gen_range(0..n));

        // Remaining centers: weighted by D^2
        for _ in 1..k {
            let mut min_dists: Vec<f64> = vec![f64::MAX; n];
            for i in 0..n {
                for &c in &centers_idx {
                    let d = dist_matrix[i][c];
                    if d < min_dists[i] {
                        min_dists[i] = d;
                    }
                }
            }

            // Square distances
            let d2: Vec<f64> = min_dists.iter().map(|d| d * d).collect();
            let total: f64 = d2.iter().sum();

            if total <= 0.0 {
                // All points are centers, pick random
                centers_idx.push(rng.gen_range(0..n));
            } else {
                // Weighted random selection
                let mut r = rng.gen::<f64>() * total;
                let mut chosen = 0;
                for (i, &d) in d2.iter().enumerate() {
                    r -= d;
                    if r <= 0.0 {
                        chosen = i;
                        break;
                    }
                }
                centers_idx.push(chosen);
            }
        }

        // Initialize centers from selected indices
        let mut centers: Vec<f64> = vec![0.0; k * m];
        for (c, &idx) in centers_idx.iter().enumerate() {
            for j in 0..m {
                centers[c + j * k] = data_slice[idx + j * n];
            }
        }

        let mut cluster: Vec<usize> = vec![0; n];
        let mut sizes: Vec<usize> = vec![0; k];

        // Iterate until convergence or max_iter
        for _ in 0..max_iter.max(1) {
            // Assignment step: assign each curve to nearest center
            let new_cluster: Vec<usize> = (0..n)
                .into_par_iter()
                .map(|i| {
                    let mut best_c = 0;
                    let mut best_d = f64::MAX;

                    for c in 0..k {
                        let mut d = 0.0;
                        match metric {
                            "L2" | "euclidean" => {
                                let weights = simpsons_weights(&argvals);
                                for j in 0..m {
                                    let diff = data_slice[i + j * n] - centers[c + j * k];
                                    d += diff * diff * weights[j];
                                }
                                d = d.sqrt();
                            }
                            "L1" | "manhattan" => {
                                let weights = simpsons_weights(&argvals);
                                for j in 0..m {
                                    let diff = (data_slice[i + j * n] - centers[c + j * k]).abs();
                                    d += diff * weights[j];
                                }
                            }
                            _ => {
                                let weights = simpsons_weights(&argvals);
                                for j in 0..m {
                                    let diff = data_slice[i + j * n] - centers[c + j * k];
                                    d += diff * diff * weights[j];
                                }
                                d = d.sqrt();
                            }
                        }

                        if d < best_d {
                            best_d = d;
                            best_c = c;
                        }
                    }
                    best_c
                })
                .collect();

            // Check for convergence
            let converged = new_cluster == cluster;
            cluster = new_cluster;

            // Update step: compute new centers as mean of assigned curves
            sizes = vec![0; k];
            for &c in &cluster {
                sizes[c] += 1;
            }

            let mut new_centers = vec![0.0; k * m];
            for (i, &c) in cluster.iter().enumerate() {
                for j in 0..m {
                    new_centers[c + j * k] += data_slice[i + j * n];
                }
            }
            for c in 0..k {
                if sizes[c] > 0 {
                    for j in 0..m {
                        new_centers[c + j * k] /= sizes[c] as f64;
                    }
                }
            }
            centers = new_centers;

            if converged {
                break;
            }
        }

        // Compute within-cluster sum of squares
        let mut withinss = vec![0.0; k];
        let weights = simpsons_weights(&argvals);
        for (i, &c) in cluster.iter().enumerate() {
            let mut d2 = 0.0;
            for j in 0..m {
                let diff = data_slice[i + j * n] - centers[c + j * k];
                d2 += diff * diff * weights[j];
            }
            withinss[c] += d2;
        }
        let tot_withinss: f64 = withinss.iter().sum();

        // Keep best result
        if tot_withinss < best_withinss {
            best_withinss = tot_withinss;
            best_cluster = cluster.iter().map(|&c| (c + 1) as i32).collect(); // 1-indexed
            best_centers = centers;
            best_sizes = sizes.iter().map(|&s| s as i32).collect();
        }
    }

    // Create centers matrix (k x m)
    let centers_mat = RMatrix::new_matrix(k, m, |i, j| best_centers[i + j * k]);

    // Create withinss vector
    let weights = simpsons_weights(&argvals);
    let mut withinss = vec![0.0; k];
    for (i, &c) in best_cluster.iter().enumerate() {
        let ci = (c - 1) as usize; // Convert back to 0-indexed
        let mut d2 = 0.0;
        for j in 0..m {
            let diff = data_slice[i + j * n] - best_centers[ci + j * k];
            d2 += diff * diff * weights[j];
        }
        withinss[ci] += d2;
    }

    list!(
        cluster = best_cluster,
        centers = centers_mat,
        withinss = withinss,
        tot_withinss = best_withinss,
        size = best_sizes
    )
    .into()
}

/// Fuzzy C-Means clustering for functional data
/// m_fuzz is the fuzziness parameter (typically 2)
#[extendr]
fn fuzzycmeans_fd(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    nclusters: i32,
    m_fuzz: f64,
    max_iter: i32,
    tol: f64,
    seed: Nullable<i32>,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let k = nclusters as usize;

    if n == 0 || m == 0 || k == 0 || k > n || m_fuzz <= 1.0 {
        return r!(list!(
            membership = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            centers = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            cluster = Vec::<i32>::new(),
            objective = 0.0
        ));
    }

    let data_slice = data.as_real_slice().unwrap();
    let weights = simpsons_weights(&argvals);

    // Initialize RNG
    let mut rng = match seed {
        Nullable::NotNull(s) => rand::rngs::StdRng::seed_from_u64(s as u64),
        Nullable::Null => rand::rngs::StdRng::from_entropy(),
    };

    // Initialize membership matrix randomly and normalize rows to sum to 1
    let mut membership: Vec<f64> = (0..(n * k)).map(|_| rng.gen::<f64>()).collect();

    // Normalize each row
    for i in 0..n {
        let row_sum: f64 = (0..k).map(|c| membership[i * k + c]).sum();
        for c in 0..k {
            membership[i * k + c] /= row_sum;
        }
    }

    // Helper to compute squared L2 distance
    let dist_sq = |i: usize, center_j: usize, centers: &[f64]| -> f64 {
        let mut d = 0.0;
        for t in 0..m {
            let diff = data_slice[i + t * n] - centers[center_j + t * k];
            d += diff * diff * weights[t];
        }
        d
    };

    let mut centers = vec![0.0; k * m];
    let max_iter = max_iter as usize;
    let exponent = 2.0 / (m_fuzz - 1.0);

    for _ in 0..max_iter {
        // Update centers
        // c_j = sum_i (u_ij^m * x_i) / sum_i (u_ij^m)
        for c in 0..k {
            let mut weighted_sum = vec![0.0; m];
            let mut weight_total = 0.0;

            for i in 0..n {
                let u_pow = membership[i * k + c].powf(m_fuzz);
                weight_total += u_pow;
                for t in 0..m {
                    weighted_sum[t] += u_pow * data_slice[i + t * n];
                }
            }

            if weight_total > 1e-10 {
                for t in 0..m {
                    centers[c + t * k] = weighted_sum[t] / weight_total;
                }
            }
        }

        // Update membership matrix
        let old_membership = membership.clone();

        for i in 0..n {
            // Compute distances to all centers
            let dists: Vec<f64> = (0..k).map(|c| dist_sq(i, c, &centers)).collect();

            // Check for zero distances
            let has_zero = dists.iter().any(|&d| d < 1e-20);

            if has_zero {
                // Assign full membership to closest center
                for c in 0..k {
                    membership[i * k + c] = if dists[c] < 1e-20 { 1.0 } else { 0.0 };
                }
            } else {
                // Standard FCM update
                for c in 0..k {
                    let mut sum = 0.0;
                    for c2 in 0..k {
                        sum += (dists[c] / dists[c2]).powf(exponent);
                    }
                    membership[i * k + c] = 1.0 / sum;
                }
            }
        }

        // Check convergence
        let diff: f64 = membership
            .iter()
            .zip(old_membership.iter())
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
            / (n * k) as f64;

        if diff < tol {
            break;
        }
    }

    // Compute final objective function
    let mut objective = 0.0;
    for i in 0..n {
        for c in 0..k {
            let u_pow = membership[i * k + c].powf(m_fuzz);
            objective += u_pow * dist_sq(i, c, &centers);
        }
    }

    // Compute hard cluster assignments (argmax of membership)
    let cluster: Vec<i32> = (0..n)
        .map(|i| {
            let mut best_c = 0;
            let mut best_u = membership[i * k];
            for c in 1..k {
                if membership[i * k + c] > best_u {
                    best_u = membership[i * k + c];
                    best_c = c;
                }
            }
            (best_c + 1) as i32 // 1-indexed
        })
        .collect();

    // Create output matrices
    let membership_mat = RMatrix::new_matrix(n, k, |i, c| membership[i * k + c]);
    let centers_mat = RMatrix::new_matrix(k, m, |c, t| centers[c + t * k]);

    list!(
        membership = membership_mat,
        centers = centers_mat,
        cluster = cluster,
        objective = objective
    )
    .into()
}

/// Shift registration: find optimal horizontal shift for each curve
/// to align with a target (usually the mean)
#[extendr]
fn register_shift_1d(
    data: RMatrix<f64>,
    target: Vec<f64>,
    argvals: Vec<f64>,
    max_shift: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m == 0 || argvals.len() != m || target.len() != m {
        return r!(list!(
            shifts = Vec::<f64>::new(),
            registered = RMatrix::new_matrix(0, 0, |_, _| 0.0)
        ));
    }

    let data_slice = data.as_real_slice().unwrap();
    let dt = if m > 1 { argvals[1] - argvals[0] } else { 1.0 };

    // Maximum shift in indices
    let max_shift_idx = ((max_shift / dt).abs().ceil() as usize).min(m / 4);

    // Find optimal shift for each curve using cross-correlation
    let results: Vec<(f64, Vec<f64>)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let curve: Vec<f64> = (0..m).map(|j| data_slice[i + j * n]).collect();

            let mut best_shift = 0i64;
            let mut best_corr = f64::NEG_INFINITY;

            // Try different shifts
            for shift in -(max_shift_idx as i64)..=(max_shift_idx as i64) {
                let mut corr = 0.0;
                let mut count = 0;

                for j in 0..m {
                    let shifted_j = (j as i64 + shift) as usize;
                    if shifted_j < m {
                        corr += curve[shifted_j] * target[j];
                        count += 1;
                    }
                }

                if count > 0 {
                    corr /= count as f64;
                    if corr > best_corr {
                        best_corr = corr;
                        best_shift = shift;
                    }
                }
            }

            // Apply the shift (using linear interpolation for sub-index shifts)
            let shift_amount = best_shift as f64 * dt;
            let registered: Vec<f64> = (0..m)
                .map(|j| {
                    let src_idx = j as f64 - best_shift as f64;
                    if src_idx < 0.0 || src_idx >= (m - 1) as f64 {
                        // Extrapolate with boundary value
                        if src_idx < 0.0 {
                            curve[0]
                        } else {
                            curve[m - 1]
                        }
                    } else {
                        // Linear interpolation
                        let lo = src_idx.floor() as usize;
                        let hi = (lo + 1).min(m - 1);
                        let frac = src_idx - lo as f64;
                        curve[lo] * (1.0 - frac) + curve[hi] * frac
                    }
                })
                .collect();

            (shift_amount, registered)
        })
        .collect();

    // Collect results
    let shifts: Vec<f64> = results.iter().map(|(s, _)| *s).collect();
    let registered_mat = RMatrix::new_matrix(n, m, |i, j| results[i].1[j]);

    list!(shifts = shifts, registered = registered_mat).into()
}

// =============================================================================
// Utility Functions
// =============================================================================

/// Simpson's rule integration for functional data
/// Integrates each curve over the domain
#[extendr]
fn int_simpson(data: RMatrix<f64>, argvals: Vec<f64>) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m == 0 || argvals.len() != m {
        return Robj::from(Vec::<f64>::new());
    }

    let data_slice = data.as_real_slice().unwrap();
    let weights = simpsons_weights(&argvals);

    let integrals: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut sum = 0.0;
            for j in 0..m {
                sum += data_slice[i + j * n] * weights[j];
            }
            sum
        })
        .collect();

    Robj::from(integrals)
}

/// Inner product of two functional data objects
/// <f, g> = integral(f(t) * g(t) dt)
#[extendr]
fn inprod_fdata(data1: RMatrix<f64>, data2: RMatrix<f64>, argvals: Vec<f64>) -> Robj {
    let n1 = data1.nrows();
    let n2 = data2.nrows();
    let m = data1.ncols();

    if m != data2.ncols() || argvals.len() != m {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    let data1_slice = data1.as_real_slice().unwrap();
    let data2_slice = data2.as_real_slice().unwrap();
    let weights = simpsons_weights(&argvals);

    // Compute n1 x n2 matrix of inner products
    let result: Vec<Vec<f64>> = (0..n1)
        .into_par_iter()
        .map(|i| {
            (0..n2)
                .map(|j| {
                    let mut sum = 0.0;
                    for k in 0..m {
                        sum += data1_slice[i + k * n1] * data2_slice[j + k * n2] * weights[k];
                    }
                    sum
                })
                .collect()
        })
        .collect();

    let mat = RMatrix::new_matrix(n1, n2, |i, j| result[i][j]);
    r!(mat)
}

// =============================================================================
// kNN Regression functions (global and local cross-validation)
// =============================================================================

/// Kernel prediction with fixed bandwidth for prediction on new data
#[extendr]
fn knn_predict(
    dist_matrix: RMatrix<f64>,
    response: Vec<f64>,
    k: i32,
    local_k: Nullable<Vec<i32>>,
) -> Robj {
    let n = dist_matrix.nrows();
    let m = dist_matrix.ncols();
    let dist_slice = dist_matrix.as_real_slice().unwrap();

    // Extract distance matrix in column-major format
    let get_dist = |i: usize, j: usize| dist_slice[i + j * n];

    let predictions: Vec<f64> = (0..m)
        .into_par_iter()
        .map(|j| {
            // Get distances for column j
            let mut dists: Vec<(usize, f64)> = (0..n).map(|i| (i, get_dist(i, j))).collect();
            dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

            // Determine bandwidth based on k-th neighbor
            let k_val = match &local_k {
                Nullable::NotNull(k_vec) => {
                    // Local: find the nearest training point and use its k
                    let nearest_train_idx = dists[0].0;
                    k_vec[nearest_train_idx] as usize
                }
                Nullable::Null => k as usize,
            };

            // Bandwidth is midpoint between k-th and (k+1)-th neighbor
            let h = if k_val + 1 < n {
                0.5 * (dists[k_val].1 + dists[k_val + 1].1)
            } else {
                dists[k_val].1 * 1.1
            };

            // Epanechnikov kernel weights
            let mut sum_ky = 0.0;
            let mut sum_k = 0.0;
            for i in 0..n {
                let u = get_dist(i, j) / h;
                if u <= 1.0 {
                    let k_val = 1.0 - u * u;
                    sum_ky += k_val * response[i];
                    sum_k += k_val;
                }
            }

            if sum_k > 0.0 {
                sum_ky / sum_k
            } else {
                0.0
            }
        })
        .collect();

    Robj::from(predictions)
}

/// k-NN with Global Cross-Validation
/// Finds a single optimal k for all observations
#[extendr]
fn knn_gcv(dist_matrix: RMatrix<f64>, response: Vec<f64>, max_k: i32) -> Robj {
    let n = dist_matrix.nrows();
    let dist_slice = dist_matrix.as_real_slice().unwrap();

    // Extract distance matrix
    let get_dist = |i: usize, j: usize| dist_slice[i + j * n];

    // Sort distances for each observation (leave-one-out style)
    let mut sorted_indices: Vec<Vec<usize>> = Vec::with_capacity(n);
    let mut sorted_dists: Vec<Vec<f64>> = Vec::with_capacity(n);

    for j in 0..n {
        let mut dists: Vec<(usize, f64)> = (0..n).map(|i| (i, get_dist(i, j))).collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        sorted_indices.push(dists.iter().map(|x| x.0).collect());
        sorted_dists.push(dists.iter().map(|x| x.1).collect());
    }

    let max_k = (max_k as usize).min(n - 2);
    let mut mse_vec: Vec<f64> = Vec::with_capacity(max_k);

    // Try each k value
    for k in 1..=max_k {
        let mse: f64 = (0..n)
            .into_par_iter()
            .map(|j| {
                // Leave-one-out: skip self (index 0 in sorted is self with distance 0)
                // Bandwidth from k-th and (k+1)-th neighbors (excluding self)
                let h = 0.5 * (sorted_dists[j][k] + sorted_dists[j][k + 1]);

                let mut sum_ky = 0.0;
                let mut sum_k = 0.0;

                for i in 0..n {
                    if i == j {
                        continue;
                    } // Leave-one-out
                    let u = get_dist(i, j) / h;
                    if u <= 1.0 {
                        let kernel_val = 1.0 - u * u;
                        sum_ky += kernel_val * response[i];
                        sum_k += kernel_val;
                    }
                }

                let pred = if sum_k > 0.0 {
                    sum_ky / sum_k
                } else {
                    response[j]
                };
                (pred - response[j]).powi(2)
            })
            .sum();

        mse_vec.push(mse / n as f64);
    }

    // Find optimal k
    let k_opt = mse_vec
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i + 1)
        .unwrap_or(1);

    // Compute final predictions with optimal k
    let yhat: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|j| {
            let h = 0.5 * (sorted_dists[j][k_opt] + sorted_dists[j][k_opt + 1]);

            let mut sum_ky = 0.0;
            let mut sum_k = 0.0;

            for i in 0..n {
                let u = get_dist(i, j) / h;
                if u <= 1.0 {
                    let kernel_val = 1.0 - u * u;
                    sum_ky += kernel_val * response[i];
                    sum_k += kernel_val;
                }
            }

            if sum_k > 0.0 {
                sum_ky / sum_k
            } else {
                response[j]
            }
        })
        .collect();

    list!(k_opt = k_opt as i32, mse = mse_vec, yhat = yhat).into()
}

/// k-NN with Local Cross-Validation
/// Finds an optimal k for each observation
#[extendr]
fn knn_lcv(dist_matrix: RMatrix<f64>, response: Vec<f64>, max_k: i32) -> Robj {
    let n = dist_matrix.nrows();
    let dist_slice = dist_matrix.as_real_slice().unwrap();

    let get_dist = |i: usize, j: usize| dist_slice[i + j * n];

    let max_k = (max_k as usize).min(n - 2);

    // For each observation, find optimal local k
    let results: Vec<(i32, f64)> = (0..n)
        .into_par_iter()
        .map(|j| {
            // Sort neighbors by distance
            let mut dists: Vec<(usize, f64)> = (0..n)
                .filter(|&i| i != j)
                .map(|i| (i, get_dist(i, j)))
                .collect();
            dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

            let mut best_k = 1usize;
            let mut best_error = f64::MAX;
            let mut best_pred = response[j];

            // Try each k
            for k in 1..=max_k.min(dists.len() - 1) {
                let h = 0.5 * (dists[k - 1].1 + dists[k].1);

                let mut sum_ky = 0.0;
                let mut sum_k = 0.0;

                for (idx, d) in &dists {
                    let u = d / h;
                    if u <= 1.0 {
                        let kernel_val = 1.0 - u * u;
                        sum_ky += kernel_val * response[*idx];
                        sum_k += kernel_val;
                    }
                }

                let pred = if sum_k > 0.0 {
                    sum_ky / sum_k
                } else {
                    response[j]
                };
                let error = (pred - response[j]).abs();

                if error < best_error {
                    best_error = error;
                    best_k = k;
                    best_pred = pred;
                }
            }

            (best_k as i32, best_pred)
        })
        .collect();

    let k_opt: Vec<i32> = results.iter().map(|x| x.0).collect();
    let yhat: Vec<f64> = results.iter().map(|x| x.1).collect();

    // Compute MSE
    let mse: f64 = yhat
        .iter()
        .zip(response.iter())
        .map(|(pred, actual)| (pred - actual).powi(2))
        .sum::<f64>()
        / n as f64;

    list!(k_opt = k_opt, mse = mse, yhat = yhat).into()
}

// =============================================================================
// Cluster validation functions
// =============================================================================

/// Compute silhouette score for clustering
/// Returns the mean silhouette coefficient across all samples
#[extendr]
fn silhouette_score(dist_matrix: RMatrix<f64>, clusters: Vec<i32>) -> f64 {
    let n = dist_matrix.nrows();
    let dist_slice = dist_matrix.as_real_slice().unwrap();

    let get_dist = |i: usize, j: usize| dist_slice[i + j * n];

    // Find unique clusters
    let mut unique_clusters: Vec<i32> = clusters.clone();
    unique_clusters.sort();
    unique_clusters.dedup();
    let k = unique_clusters.len();

    if k <= 1 {
        return 0.0; // Silhouette undefined for single cluster
    }

    // Compute silhouette for each point
    let silhouettes: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let ci = clusters[i];

            // a(i) = mean distance to points in same cluster
            let same_cluster: Vec<usize> =
                (0..n).filter(|&j| j != i && clusters[j] == ci).collect();

            if same_cluster.is_empty() {
                return 0.0; // Single point in cluster
            }

            let a_i: f64 = same_cluster.iter().map(|&j| get_dist(i, j)).sum::<f64>()
                / same_cluster.len() as f64;

            // b(i) = min mean distance to points in other clusters
            let mut b_i = f64::MAX;
            for &ck in &unique_clusters {
                if ck == ci {
                    continue;
                }

                let other_cluster: Vec<usize> = (0..n).filter(|&j| clusters[j] == ck).collect();

                if other_cluster.is_empty() {
                    continue;
                }

                let mean_dist: f64 = other_cluster.iter().map(|&j| get_dist(i, j)).sum::<f64>()
                    / other_cluster.len() as f64;

                if mean_dist < b_i {
                    b_i = mean_dist;
                }
            }

            // s(i) = (b(i) - a(i)) / max(a(i), b(i))
            let max_ab = a_i.max(b_i);
            if max_ab == 0.0 {
                0.0
            } else {
                (b_i - a_i) / max_ab
            }
        })
        .collect();

    silhouettes.iter().sum::<f64>() / n as f64
}

/// Compute Calinski-Harabasz index (variance ratio criterion)
/// Higher values indicate better defined clusters
#[extendr]
fn calinski_harabasz(data: RMatrix<f64>, clusters: Vec<i32>) -> f64 {
    let n = data.nrows();
    let m = data.ncols();
    let data_slice = data.as_real_slice().unwrap();

    let get_val = |i: usize, j: usize| data_slice[i + j * n];

    // Find unique clusters
    let mut unique_clusters: Vec<i32> = clusters.clone();
    unique_clusters.sort();
    unique_clusters.dedup();
    let k = unique_clusters.len();

    if k <= 1 || k >= n {
        return 0.0;
    }

    // Compute overall centroid
    let mut global_centroid = vec![0.0; m];
    for j in 0..m {
        for i in 0..n {
            global_centroid[j] += get_val(i, j);
        }
        global_centroid[j] /= n as f64;
    }

    // Compute cluster centroids and sizes
    let mut cluster_centroids: Vec<Vec<f64>> = vec![vec![0.0; m]; k];
    let mut cluster_sizes: Vec<usize> = vec![0; k];

    for i in 0..n {
        let ci = clusters[i];
        let cluster_idx = unique_clusters.iter().position(|&x| x == ci).unwrap();
        cluster_sizes[cluster_idx] += 1;
        for j in 0..m {
            cluster_centroids[cluster_idx][j] += get_val(i, j);
        }
    }

    for (idx, centroid) in cluster_centroids.iter_mut().enumerate() {
        if cluster_sizes[idx] > 0 {
            for j in 0..m {
                centroid[j] /= cluster_sizes[idx] as f64;
            }
        }
    }

    // Between-cluster dispersion (SSB)
    let mut ssb = 0.0;
    for (idx, centroid) in cluster_centroids.iter().enumerate() {
        let mut dist_sq = 0.0;
        for j in 0..m {
            let diff = centroid[j] - global_centroid[j];
            dist_sq += diff * diff;
        }
        ssb += cluster_sizes[idx] as f64 * dist_sq;
    }

    // Within-cluster dispersion (SSW)
    let mut ssw = 0.0;
    for i in 0..n {
        let ci = clusters[i];
        let cluster_idx = unique_clusters.iter().position(|&x| x == ci).unwrap();
        let centroid = &cluster_centroids[cluster_idx];

        for j in 0..m {
            let diff = get_val(i, j) - centroid[j];
            ssw += diff * diff;
        }
    }

    // CH = (SSB / (k-1)) / (SSW / (n-k))
    if ssw == 0.0 {
        return f64::MAX;
    }

    (ssb / (k - 1) as f64) / (ssw / (n - k) as f64)
}

// =============================================================================
// Ridge Regression using anofox-regression
// TEMPORARILY DISABLED: Requires faer 0.23+ which needs Rust 1.84+
// CRAN Windows currently has Rust 1.81. Re-enable when CRAN updates.
// =============================================================================

// /// Fit ridge regression model
// ///
// /// Uses the anofox-regression crate for efficient L2-regularized regression.
// ///
// /// # Arguments
// /// * `x` - Design matrix (n x m), column-major
// /// * `y` - Response vector (n)
// /// * `lambda` - Regularization parameter (L2 penalty)
// /// * `with_intercept` - Whether to include an intercept term
// ///
// /// # Returns
// /// A list with coefficients, intercept, fitted values, residuals, and R-squared
// #[extendr]
// fn ridge_regression_fit(x: RMatrix<f64>, y: Vec<f64>, lambda: f64, with_intercept: bool) -> Robj {
//     ... (function body commented out for CRAN compatibility)
// }

// =============================================================================
// Seasonal Analysis Functions
// =============================================================================

/// Estimate period using FFT periodogram
#[extendr]
fn seasonal_estimate_period_fft(data: RMatrix<f64>, argvals: Vec<f64>) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 || argvals.len() != m {
        return list!(
            period = f64::NAN,
            frequency = f64::NAN,
            power = 0.0,
            confidence = 0.0
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();

    // Compute mean curve first
    let mean_curve: Vec<f64> = (0..m)
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..n {
                sum += data_slice[i + j * n];
            }
            sum / n as f64
        })
        .collect();

    // Compute periodogram
    let dt = (argvals[m - 1] - argvals[0]) / (m - 1) as f64;
    let fs = 1.0 / dt;

    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(m);

    let mut buffer: Vec<Complex<f64>> = mean_curve.iter().map(|&x| Complex::new(x, 0.0)).collect();
    fft.process(&mut buffer);

    let n_freq = m / 2 + 1;
    let mut max_power = 0.0;
    let mut max_idx = 1;
    let mut total_power = 0.0;

    for k in 1..n_freq {
        let p = buffer[k].norm_sqr() / (m as f64 * m as f64);
        let p = if k < m / 2 { 2.0 * p } else { p };
        total_power += p;
        if p > max_power {
            max_power = p;
            max_idx = k;
        }
    }

    let dominant_freq = max_idx as f64 * fs / m as f64;
    let period = if dominant_freq > 1e-15 {
        1.0 / dominant_freq
    } else {
        f64::INFINITY
    };

    let mean_power = total_power / (n_freq - 1) as f64;
    let confidence = if mean_power > 1e-15 {
        max_power / mean_power
    } else {
        0.0
    };

    list!(
        period = period,
        frequency = dominant_freq,
        power = max_power,
        confidence = confidence
    )
    .into()
}

/// Estimate period using autocorrelation
#[extendr]
fn seasonal_estimate_period_acf(data: RMatrix<f64>, argvals: Vec<f64>, max_lag: i32) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let max_lag = max_lag as usize;

    if n == 0 || m < 4 || argvals.len() != m {
        return list!(
            period = f64::NAN,
            frequency = f64::NAN,
            power = 0.0,
            confidence = 0.0
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();

    // Compute mean curve
    let mean_curve: Vec<f64> = (0..m)
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..n {
                sum += data_slice[i + j * n];
            }
            sum / n as f64
        })
        .collect();

    // Compute ACF
    let mean: f64 = mean_curve.iter().sum::<f64>() / m as f64;
    let var: f64 = mean_curve.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / m as f64;

    if var < 1e-15 {
        return list!(
            period = f64::NAN,
            frequency = f64::NAN,
            power = 0.0,
            confidence = 0.0
        )
        .into();
    }

    let max_lag = max_lag.min(m - 1);
    let mut acf = Vec::with_capacity(max_lag + 1);
    for lag in 0..=max_lag {
        let mut sum = 0.0;
        for i in 0..(m - lag) {
            sum += (mean_curve[i] - mean) * (mean_curve[i + lag] - mean);
        }
        acf.push(sum / (m as f64 * var));
    }

    // Find first peak after lag 0
    let min_lag = 2;
    let mut peak_lag = 0;
    for i in (min_lag + 1)..(acf.len() - 1) {
        if acf[i] > acf[i - 1] && acf[i] > acf[i + 1] {
            peak_lag = i;
            break;
        }
    }

    if peak_lag == 0 {
        return list!(
            period = f64::NAN,
            frequency = f64::NAN,
            power = 0.0,
            confidence = 0.0
        )
        .into();
    }

    let dt = (argvals[m - 1] - argvals[0]) / (m - 1) as f64;
    let period = peak_lag as f64 * dt;
    let frequency = if period > 1e-15 { 1.0 / period } else { 0.0 };

    list!(
        period = period,
        frequency = frequency,
        power = acf[peak_lag],
        confidence = acf[peak_lag].abs()
    )
    .into()
}

/// Detect multiple concurrent periodicities using iterative residual subtraction
#[extendr]
fn seasonal_detect_multiple_periods(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    max_periods: i32,
    min_confidence: f64,
    min_strength: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 || argvals.len() != m || max_periods <= 0 {
        return list!(
            periods = Vec::<f64>::new(),
            confidence = Vec::<f64>::new(),
            strength = Vec::<f64>::new(),
            amplitude = Vec::<f64>::new(),
            phase = Vec::<f64>::new(),
            iteration = Vec::<i32>::new()
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();

    let detected = fdars_core::seasonal::detect_multiple_periods(
        data_slice,
        n,
        m,
        &argvals,
        max_periods as usize,
        min_confidence,
        min_strength,
    );

    // Convert to R-friendly format
    let periods: Vec<f64> = detected.iter().map(|d| d.period).collect();
    let confidences: Vec<f64> = detected.iter().map(|d| d.confidence).collect();
    let strengths: Vec<f64> = detected.iter().map(|d| d.strength).collect();
    let amplitudes: Vec<f64> = detected.iter().map(|d| d.amplitude).collect();
    let phases: Vec<f64> = detected.iter().map(|d| d.phase).collect();
    let iterations: Vec<i32> = detected.iter().map(|d| d.iteration as i32).collect();

    list!(
        period = periods,
        confidence = confidences,
        strength = strengths,
        amplitude = amplitudes,
        phase = phases,
        iteration = iterations
    )
    .into()
}

/// Detect peaks in functional data using Fourier basis smoothing
#[extendr]
fn seasonal_detect_peaks(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    min_distance: Robj,
    min_prominence: Robj,
    smooth_first: bool,
    smooth_nbasis: Robj,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 3 || argvals.len() != m {
        return list!(
            peaks = list!(),
            inter_peak_distances = list!(),
            mean_period = f64::NAN
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let min_dist: Option<f64> = if min_distance.is_null() {
        None
    } else {
        min_distance.as_real()
    };

    let min_prom: Option<f64> = if min_prominence.is_null() {
        None
    } else {
        min_prominence.as_real()
    };

    // Parse smooth_nbasis: NULL = auto GCV selection, otherwise use provided value
    let nbasis: Option<usize> = if smooth_nbasis.is_null() {
        None // Triggers auto GCV selection of Fourier basis in Rust core
    } else {
        smooth_nbasis.as_real().map(|v| v as usize)
    };

    // Use the Rust core detect_peaks function with Fourier basis smoothing
    let result = fdars_core::seasonal::detect_peaks(
        &fd,
        &argvals,
        min_dist,
        min_prom,
        smooth_first,
        nbasis,
    );

    // Convert result to R format
    let mut all_peaks: Vec<Robj> = Vec::with_capacity(result.peaks.len());
    let mut all_distances: Vec<Robj> = Vec::with_capacity(result.inter_peak_distances.len());

    for (peaks, distances) in result.peaks.iter().zip(result.inter_peak_distances.iter()) {
        let times: Vec<f64> = peaks.iter().map(|p| p.time).collect();
        let values: Vec<f64> = peaks.iter().map(|p| p.value).collect();
        let proms: Vec<f64> = peaks.iter().map(|p| p.prominence).collect();

        all_peaks.push(list!(time = times, value = values, prominence = proms).into());
        all_distances.push(Robj::from(distances.clone()));
    }

    list!(
        peaks = List::from_values(all_peaks),
        inter_peak_distances = List::from_values(all_distances),
        mean_period = result.mean_period
    )
    .into()
}

/// Measure seasonal strength using variance decomposition
#[extendr]
fn seasonal_strength_variance(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    period: f64,
    n_harmonics: i32,
) -> f64 {
    let n = data.nrows();
    let m = data.ncols();
    let n_harmonics = n_harmonics as usize;

    if n == 0 || m < 4 || argvals.len() != m || period <= 0.0 {
        return f64::NAN;
    }

    let data_slice = data.as_real_slice().unwrap();

    // Compute mean curve
    let mean_curve: Vec<f64> = (0..m)
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..n {
                sum += data_slice[i + j * n];
            }
            sum / n as f64
        })
        .collect();

    // Total variance
    let global_mean: f64 = mean_curve.iter().sum::<f64>() / m as f64;
    let total_var: f64 = mean_curve
        .iter()
        .map(|&x| (x - global_mean).powi(2))
        .sum::<f64>()
        / m as f64;

    if total_var < 1e-15 {
        return 0.0;
    }

    // Fit Fourier basis with explicit period
    let nbasis = 1 + 2 * n_harmonics;
    let t_min = argvals.iter().cloned().fold(f64::INFINITY, f64::min);

    let mut basis = vec![0.0; m * nbasis];
    for (i, &ti) in argvals.iter().enumerate() {
        let x = 2.0 * PI * (ti - t_min) / period;
        basis[i] = 1.0;

        let mut k = 1;
        let mut freq = 1;
        while k < nbasis {
            if k < nbasis {
                basis[i + k * m] = (freq as f64 * x).sin();
                k += 1;
            }
            if k < nbasis {
                basis[i + k * m] = (freq as f64 * x).cos();
                k += 1;
            }
            freq += 1;
        }
    }

    // Project to seasonal (skip DC)
    let mut seasonal = vec![0.0; m];
    for k in 1..nbasis {
        let b_sum: f64 = (0..m).map(|j| basis[j + k * m].powi(2)).sum();
        if b_sum > 1e-15 {
            let coef: f64 = (0..m)
                .map(|j| mean_curve[j] * basis[j + k * m])
                .sum::<f64>()
                / b_sum;
            for j in 0..m {
                seasonal[j] += coef * basis[j + k * m];
            }
        }
    }

    let seasonal_mean: f64 = seasonal.iter().sum::<f64>() / m as f64;
    let seasonal_var: f64 = seasonal
        .iter()
        .map(|&x| (x - seasonal_mean).powi(2))
        .sum::<f64>()
        / m as f64;

    (seasonal_var / total_var).min(1.0)
}

/// Measure seasonal strength using spectral method
#[extendr]
fn seasonal_strength_spectral(data: RMatrix<f64>, argvals: Vec<f64>, period: f64) -> f64 {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 || argvals.len() != m || period <= 0.0 {
        return f64::NAN;
    }

    let data_slice = data.as_real_slice().unwrap();

    // Compute mean curve
    let mean_curve: Vec<f64> = (0..m)
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..n {
                sum += data_slice[i + j * n];
            }
            sum / n as f64
        })
        .collect();

    // Compute periodogram
    let dt = (argvals[m - 1] - argvals[0]) / (m - 1) as f64;
    let fs = 1.0 / dt;

    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(m);

    let mut buffer: Vec<Complex<f64>> = mean_curve.iter().map(|&x| Complex::new(x, 0.0)).collect();
    fft.process(&mut buffer);

    let n_freq = m / 2 + 1;
    let fundamental_freq = 1.0 / period;
    let mut seasonal_power = 0.0;
    let mut total_power = 0.0;

    for k in 1..n_freq {
        let freq = k as f64 * fs / m as f64;
        let p = buffer[k].norm_sqr() / (m as f64 * m as f64);
        let p = if k < m / 2 { 2.0 * p } else { p };

        total_power += p;

        // Check if near harmonic of fundamental
        let ratio = freq / fundamental_freq;
        let nearest = ratio.round();
        if (ratio - nearest).abs() < 0.1 && nearest >= 1.0 {
            seasonal_power += p;
        }
    }

    if total_power < 1e-15 {
        return 0.0;
    }

    (seasonal_power / total_power).min(1.0)
}

/// Measure seasonal strength using wavelet (Morlet) method
///
/// Uses Continuous Wavelet Transform with Morlet wavelet to measure
/// power at the specified seasonal period.
///
/// @param data Matrix of functional data (n x m)
/// @param argvals Vector of evaluation points (length m)
/// @param period Seasonal period in argvals units
/// @return Seasonal strength as ratio of wavelet power to total variance (0 to 1)
#[extendr]
fn seasonal_strength_wavelet(data: RMatrix<f64>, argvals: Vec<f64>, period: f64) -> f64 {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 || argvals.len() != m || period <= 0.0 {
        return f64::NAN;
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    fdars_core::seasonal::seasonal_strength_wavelet(&fd, &argvals, period)
}

/// Time-varying seasonal strength using sliding windows
#[extendr]
fn seasonal_strength_windowed(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    period: f64,
    window_size: f64,
    method: &str,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 || argvals.len() != m || period <= 0.0 || window_size <= 0.0 {
        return Robj::from(Vec::<f64>::new());
    }

    let data_slice = data.as_real_slice().unwrap();
    let dt = (argvals[m - 1] - argvals[0]) / (m - 1) as f64;
    let half_window_points = ((window_size / 2.0) / dt).round() as usize;

    // Compute mean curve
    let mean_curve: Vec<f64> = (0..m)
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..n {
                sum += data_slice[i + j * n];
            }
            sum / n as f64
        })
        .collect();

    let use_spectral = method == "spectral";

    let strength: Vec<f64> = (0..m)
        .map(|center| {
            let start = center.saturating_sub(half_window_points);
            let end = (center + half_window_points + 1).min(m);
            let window_m = end - start;

            if window_m < 4 {
                return f64::NAN;
            }

            let window_data: Vec<f64> = mean_curve[start..end].to_vec();
            let window_argvals: Vec<f64> = argvals[start..end].to_vec();

            if use_spectral {
                // Spectral method inline
                let w_dt =
                    (window_argvals[window_m - 1] - window_argvals[0]) / (window_m - 1) as f64;
                let fs = 1.0 / w_dt;

                let mut planner = FftPlanner::<f64>::new();
                let fft = planner.plan_fft_forward(window_m);
                let mut buffer: Vec<Complex<f64>> =
                    window_data.iter().map(|&x| Complex::new(x, 0.0)).collect();
                fft.process(&mut buffer);

                let n_freq = window_m / 2 + 1;
                let fundamental = 1.0 / period;
                let mut seasonal_p = 0.0;
                let mut total_p = 0.0;

                for k in 1..n_freq {
                    let freq = k as f64 * fs / window_m as f64;
                    let p = buffer[k].norm_sqr() / (window_m as f64 * window_m as f64);
                    let p = if k < window_m / 2 { 2.0 * p } else { p };
                    total_p += p;
                    let ratio = freq / fundamental;
                    let nearest = ratio.round();
                    if (ratio - nearest).abs() < 0.1 && nearest >= 1.0 {
                        seasonal_p += p;
                    }
                }
                if total_p < 1e-15 {
                    0.0
                } else {
                    (seasonal_p / total_p).min(1.0)
                }
            } else {
                // Variance method inline
                let global_mean: f64 = window_data.iter().sum::<f64>() / window_m as f64;
                let total_var: f64 = window_data
                    .iter()
                    .map(|&x| (x - global_mean).powi(2))
                    .sum::<f64>()
                    / window_m as f64;

                if total_var < 1e-15 {
                    return 0.0;
                }

                let nbasis = 7; // 1 + 2*3 harmonics
                let t_min = window_argvals[0];
                let mut basis = vec![0.0; window_m * nbasis];
                for (i, &ti) in window_argvals.iter().enumerate() {
                    let x = 2.0 * PI * (ti - t_min) / period;
                    basis[i] = 1.0;
                    let mut k = 1;
                    let mut freq = 1;
                    while k < nbasis {
                        if k < nbasis {
                            basis[i + k * window_m] = (freq as f64 * x).sin();
                            k += 1;
                        }
                        if k < nbasis {
                            basis[i + k * window_m] = (freq as f64 * x).cos();
                            k += 1;
                        }
                        freq += 1;
                    }
                }

                let mut seasonal = vec![0.0; window_m];
                for k in 1..nbasis {
                    let b_sum: f64 = (0..window_m).map(|j| basis[j + k * window_m].powi(2)).sum();
                    if b_sum > 1e-15 {
                        let coef: f64 = (0..window_m)
                            .map(|j| window_data[j] * basis[j + k * window_m])
                            .sum::<f64>()
                            / b_sum;
                        for j in 0..window_m {
                            seasonal[j] += coef * basis[j + k * window_m];
                        }
                    }
                }

                let s_mean: f64 = seasonal.iter().sum::<f64>() / window_m as f64;
                let s_var: f64 =
                    seasonal.iter().map(|&x| (x - s_mean).powi(2)).sum::<f64>() / window_m as f64;

                (s_var / total_var).min(1.0)
            }
        })
        .collect();

    Robj::from(strength)
}

/// Detect seasonality changes (onset/cessation)
#[extendr]
fn seasonal_detect_changes(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    period: f64,
    threshold: f64,
    window_size: f64,
    min_duration: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 || argvals.len() != m {
        return list!(
            change_times = Vec::<f64>::new(),
            change_types = Vec::<String>::new(),
            strength_before = Vec::<f64>::new(),
            strength_after = Vec::<f64>::new(),
            strength_curve = Vec::<f64>::new()
        )
        .into();
    }

    // Compute windowed strength
    let data_slice = data.as_real_slice().unwrap();
    let dt = (argvals[m - 1] - argvals[0]) / (m - 1) as f64;
    let half_window = ((window_size / 2.0) / dt).round() as usize;
    let min_dur_points = (min_duration / dt).round() as usize;

    let mean_curve: Vec<f64> = (0..m)
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..n {
                sum += data_slice[i + j * n];
            }
            sum / n as f64
        })
        .collect();

    // Compute strength curve (variance method)
    let strength_curve: Vec<f64> = (0..m)
        .map(|center| {
            let start = center.saturating_sub(half_window);
            let end = (center + half_window + 1).min(m);
            let wm = end - start;
            if wm < 4 {
                return f64::NAN;
            }

            let wd: Vec<f64> = mean_curve[start..end].to_vec();
            let wa: Vec<f64> = argvals[start..end].to_vec();

            let gm: f64 = wd.iter().sum::<f64>() / wm as f64;
            let tv: f64 = wd.iter().map(|&x| (x - gm).powi(2)).sum::<f64>() / wm as f64;
            if tv < 1e-15 {
                return 0.0;
            }

            let nbasis = 7;
            let t0 = wa[0];
            let mut basis = vec![0.0; wm * nbasis];
            for (i, &ti) in wa.iter().enumerate() {
                let x = 2.0 * PI * (ti - t0) / period;
                basis[i] = 1.0;
                let mut k = 1;
                let mut f = 1;
                while k < nbasis {
                    if k < nbasis {
                        basis[i + k * wm] = (f as f64 * x).sin();
                        k += 1;
                    }
                    if k < nbasis {
                        basis[i + k * wm] = (f as f64 * x).cos();
                        k += 1;
                    }
                    f += 1;
                }
            }

            let mut seasonal = vec![0.0; wm];
            for k in 1..nbasis {
                let bs: f64 = (0..wm).map(|j| basis[j + k * wm].powi(2)).sum();
                if bs > 1e-15 {
                    let c: f64 = (0..wm).map(|j| wd[j] * basis[j + k * wm]).sum::<f64>() / bs;
                    for j in 0..wm {
                        seasonal[j] += c * basis[j + k * wm];
                    }
                }
            }

            let sm: f64 = seasonal.iter().sum::<f64>() / wm as f64;
            let sv: f64 = seasonal.iter().map(|&x| (x - sm).powi(2)).sum::<f64>() / wm as f64;
            (sv / tv).min(1.0)
        })
        .collect();

    // Detect threshold crossings
    let mut change_times: Vec<f64> = Vec::new();
    let mut change_types: Vec<String> = Vec::new();
    let mut strength_before: Vec<f64> = Vec::new();
    let mut strength_after: Vec<f64> = Vec::new();

    let first_valid = strength_curve
        .iter()
        .position(|&x| !x.is_nan())
        .unwrap_or(0);
    let mut in_seasonal = strength_curve[first_valid] > threshold;
    let mut last_change: Option<usize> = None;

    for (i, &ss) in strength_curve.iter().enumerate().skip(first_valid + 1) {
        if ss.is_nan() {
            continue;
        }

        let now_seasonal = ss > threshold;
        if now_seasonal != in_seasonal {
            if let Some(last) = last_change {
                if i - last < min_dur_points {
                    continue;
                }
            }

            let sb = if i > 0 && !strength_curve[i - 1].is_nan() {
                strength_curve[i - 1]
            } else {
                ss
            };

            change_times.push(argvals[i]);
            change_types.push(if now_seasonal {
                "onset".to_string()
            } else {
                "cessation".to_string()
            });
            strength_before.push(sb);
            strength_after.push(ss);

            in_seasonal = now_seasonal;
            last_change = Some(i);
        }
    }

    list!(
        change_times = change_times,
        change_types = change_types,
        strength_before = strength_before,
        strength_after = strength_after,
        strength_curve = strength_curve
    )
    .into()
}

/// Estimate instantaneous period using Hilbert transform
#[extendr]
fn seasonal_instantaneous_period(data: RMatrix<f64>, argvals: Vec<f64>) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 || argvals.len() != m {
        return list!(
            period = Vec::<f64>::new(),
            frequency = Vec::<f64>::new(),
            amplitude = Vec::<f64>::new()
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();

    // Compute mean curve
    let mean_curve: Vec<f64> = (0..m)
        .map(|j| {
            let mut sum = 0.0;
            for i in 0..n {
                sum += data_slice[i + j * n];
            }
            sum / n as f64
        })
        .collect();

    // Detrend
    let dc: f64 = mean_curve.iter().sum::<f64>() / m as f64;
    let detrended: Vec<f64> = mean_curve.iter().map(|&x| x - dc).collect();

    // Hilbert transform using shared implementation
    let analytic = hilbert_transform(&detrended);

    // Extract amplitude and phase
    let amplitude: Vec<f64> = analytic.iter().map(|c| c.norm()).collect();
    let phase: Vec<f64> = analytic.iter().map(|c| c.im.atan2(c.re)).collect();

    // Unwrap phase
    let mut unwrapped = vec![phase[0]];
    let mut cumulative = 0.0;
    for i in 1..phase.len() {
        let diff = phase[i] - phase[i - 1];
        if diff > PI {
            cumulative -= 2.0 * PI;
        } else if diff < -PI {
            cumulative += 2.0 * PI;
        }
        unwrapped.push(phase[i] + cumulative);
    }

    // Instantaneous frequency
    let dt = (argvals[m - 1] - argvals[0]) / (m - 1) as f64;
    let mut inst_freq = vec![0.0; m];

    if m > 1 {
        inst_freq[0] = (unwrapped[1] - unwrapped[0]) / dt / (2.0 * PI);
    }
    for j in 1..(m - 1) {
        inst_freq[j] = (unwrapped[j + 1] - unwrapped[j - 1]) / (2.0 * dt) / (2.0 * PI);
    }
    if m > 1 {
        inst_freq[m - 1] = (unwrapped[m - 1] - unwrapped[m - 2]) / dt / (2.0 * PI);
    }

    // Period
    let period: Vec<f64> = inst_freq
        .iter()
        .map(|&f| {
            if f.abs() > 1e-10 {
                (1.0 / f).abs()
            } else {
                f64::INFINITY
            }
        })
        .collect();

    list!(
        period = period,
        frequency = inst_freq,
        amplitude = amplitude
    )
    .into()
}

/// Analyze peak timing variability across cycles (uses Fourier smoothing)
#[extendr]
fn seasonal_analyze_peak_timing(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    period: f64,
    smooth_nbasis: Robj,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 3 || argvals.len() != m || period <= 0.0 {
        return list!(
            peak_times = Vec::<f64>::new(),
            peak_values = Vec::<f64>::new(),
            normalized_timing = Vec::<f64>::new(),
            mean_timing = f64::NAN,
            std_timing = f64::NAN,
            range_timing = f64::NAN,
            variability_score = f64::NAN,
            timing_trend = f64::NAN,
            cycle_indices = Vec::<i32>::new()
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let nbasis: Option<usize> = if smooth_nbasis.is_null() {
        None
    } else {
        smooth_nbasis.as_real().map(|v| v as usize)
    };

    let result =
        fdars_core::seasonal::analyze_peak_timing(&fd, &argvals, period, nbasis);

    let cycle_indices: Vec<i32> = result.cycle_indices.iter().map(|&i| i as i32).collect();

    list!(
        peak_times = result.peak_times,
        peak_values = result.peak_values,
        normalized_timing = result.normalized_timing,
        mean_timing = result.mean_timing,
        std_timing = result.std_timing,
        range_timing = result.range_timing,
        variability_score = result.variability_score,
        timing_trend = result.timing_trend,
        cycle_indices = cycle_indices
    )
    .into()
}

/// Classify seasonality type
#[extendr]
fn seasonal_classify_seasonality(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    period: f64,
    strength_threshold: Robj,
    timing_threshold: Robj,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 || argvals.len() != m || period <= 0.0 {
        return list!(
            is_seasonal = false,
            has_stable_timing = false,
            timing_variability = f64::NAN,
            seasonal_strength = f64::NAN,
            cycle_strengths = Vec::<f64>::new(),
            weak_seasons = Vec::<i32>::new(),
            classification = "NonSeasonal",
            peak_timing = Robj::from(())
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let str_thresh: Option<f64> = if strength_threshold.is_null() {
        None
    } else {
        strength_threshold.as_real()
    };

    let tim_thresh: Option<f64> = if timing_threshold.is_null() {
        None
    } else {
        timing_threshold.as_real()
    };

    let result = fdars_core::seasonal::classify_seasonality(
        &fd, &argvals, period, str_thresh, tim_thresh,
    );

    let weak_seasons: Vec<i32> = result.weak_seasons.iter().map(|&i| i as i32).collect();

    let classification_str = match result.classification {
        fdars_core::seasonal::SeasonalType::StableSeasonal => "StableSeasonal",
        fdars_core::seasonal::SeasonalType::VariableTiming => "VariableTiming",
        fdars_core::seasonal::SeasonalType::IntermittentSeasonal => "IntermittentSeasonal",
        fdars_core::seasonal::SeasonalType::NonSeasonal => "NonSeasonal",
        _ => "Unknown",
    };

    let peak_timing_robj = if let Some(pt) = result.peak_timing {
        let cycle_indices: Vec<i32> = pt.cycle_indices.iter().map(|&i| i as i32).collect();

        list!(
            peak_times = pt.peak_times,
            peak_values = pt.peak_values,
            normalized_timing = pt.normalized_timing,
            mean_timing = pt.mean_timing,
            std_timing = pt.std_timing,
            range_timing = pt.range_timing,
            variability_score = pt.variability_score,
            timing_trend = pt.timing_trend,
            cycle_indices = cycle_indices
        )
        .into()
    } else {
        Robj::from(())
    };

    list!(
        is_seasonal = result.is_seasonal,
        has_stable_timing = result.has_stable_timing,
        timing_variability = result.timing_variability,
        seasonal_strength = result.seasonal_strength,
        cycle_strengths = result.cycle_strengths,
        weak_seasons = weak_seasons,
        classification = classification_str,
        peak_timing = peak_timing_robj
    )
    .into()
}

/// Detect seasonality changes with automatic threshold
#[extendr]
fn seasonal_detect_changes_auto(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    period: f64,
    threshold_method: &str,
    threshold_value: Robj,
    window_size: f64,
    min_duration: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 || argvals.len() != m {
        return list!(
            change_times = Vec::<f64>::new(),
            change_types = Vec::<String>::new(),
            strength_before = Vec::<f64>::new(),
            strength_after = Vec::<f64>::new(),
            strength_curve = Vec::<f64>::new(),
            computed_threshold = f64::NAN
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let method = match threshold_method {
        "fixed" => {
            let val = threshold_value.as_real().unwrap_or(0.3);
            fdars_core::seasonal::ThresholdMethod::Fixed(val)
        }
        "percentile" => {
            let val = threshold_value.as_real().unwrap_or(20.0);
            fdars_core::seasonal::ThresholdMethod::Percentile(val)
        }
        "otsu" | _ => fdars_core::seasonal::ThresholdMethod::Otsu,
    };

    let result = fdars_core::seasonal::detect_seasonality_changes_auto(
        &fd,
        &argvals,
        period,
        method,
        window_size,
        min_duration,
    );

    let change_times: Vec<f64> = result.change_points.iter().map(|cp| cp.time).collect();
    let change_types: Vec<String> = result
        .change_points
        .iter()
        .map(|cp| match cp.change_type {
            fdars_core::seasonal::ChangeType::Onset => "onset".to_string(),
            fdars_core::seasonal::ChangeType::Cessation => "cessation".to_string(),
            _ => "unknown".to_string(),
        })
        .collect();
    let strength_before: Vec<f64> = result
        .change_points
        .iter()
        .map(|cp| cp.strength_before)
        .collect();
    let strength_after: Vec<f64> = result
        .change_points
        .iter()
        .map(|cp| cp.strength_after)
        .collect();

    // Compute the threshold that was used
    let computed_threshold = match method {
        fdars_core::seasonal::ThresholdMethod::Fixed(t) => t,
        _ => {
            // For auto methods, we need to compute the threshold from the strength curve
            let valid: Vec<f64> = result
                .strength_curve
                .iter()
                .copied()
                .filter(|x| x.is_finite())
                .collect();
            if valid.is_empty() {
                0.5
            } else if matches!(method, fdars_core::seasonal::ThresholdMethod::Percentile(_)) {
                let p = if let fdars_core::seasonal::ThresholdMethod::Percentile(p) = method {
                    p
                } else {
                    20.0
                };
                let mut sorted = valid.clone();
                sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let idx = ((p / 100.0) * sorted.len() as f64) as usize;
                sorted[idx.min(sorted.len() - 1)]
            } else {
                // Otsu - simplified computation
                otsu_threshold_inline(&valid)
            }
        }
    };

    list!(
        change_times = change_times,
        change_types = change_types,
        strength_before = strength_before,
        strength_after = strength_after,
        strength_curve = result.strength_curve,
        computed_threshold = computed_threshold
    )
    .into()
}

/// Helper function for Otsu threshold computation
fn otsu_threshold_inline(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.5;
    }

    let min_val = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_val = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    if (max_val - min_val).abs() < 1e-10 {
        return (min_val + max_val) / 2.0;
    }

    let n_bins = 256;
    let bin_width = (max_val - min_val) / n_bins as f64;
    let mut histogram = vec![0usize; n_bins];

    for &v in values {
        let bin = ((v - min_val) / bin_width).min(n_bins as f64 - 1.0) as usize;
        histogram[bin] += 1;
    }

    let total = values.len() as f64;
    let mut sum_total = 0.0;
    for (i, &count) in histogram.iter().enumerate() {
        sum_total += i as f64 * count as f64;
    }

    let mut best_threshold = min_val;
    let mut best_variance = 0.0;
    let mut sum_b = 0.0;
    let mut weight_b = 0.0;

    for t in 0..n_bins {
        weight_b += histogram[t] as f64;
        if weight_b == 0.0 {
            continue;
        }

        let weight_f = total - weight_b;
        if weight_f == 0.0 {
            break;
        }

        sum_b += t as f64 * histogram[t] as f64;
        let mean_b = sum_b / weight_b;
        let mean_f = (sum_total - sum_b) / weight_f;
        let variance = weight_b * weight_f * (mean_b - mean_f).powi(2);

        if variance > best_variance {
            best_variance = variance;
            best_threshold = min_val + (t as f64 + 0.5) * bin_width;
        }
    }

    best_threshold
}

/// Detect amplitude modulation in seasonal time series using Hilbert transform
///
/// @param data Matrix of functional data (n x m)
/// @param argvals Vector of evaluation points (length m)
/// @param period Seasonal period in argvals units
/// @param modulation_threshold CV threshold for detecting modulation (default: 0.15)
/// @param seasonality_threshold Strength threshold for seasonality (default: 0.3)
/// @return List with detection results
#[extendr]
fn seasonal_detect_amplitude_modulation(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    period: f64,
    modulation_threshold: f64,
    seasonality_threshold: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 || argvals.len() != m || period <= 0.0 {
        return list!(
            is_seasonal = false,
            seasonal_strength = f64::NAN,
            has_modulation = false,
            modulation_type = "unknown",
            modulation_score = f64::NAN,
            amplitude_trend = f64::NAN,
            strength_curve = Vec::<f64>::new(),
            time_points = Vec::<f64>::new(),
            min_strength = f64::NAN,
            max_strength = f64::NAN
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let result = fdars_core::seasonal::detect_amplitude_modulation(
        &fd,
        &argvals,
        period,
        modulation_threshold,
        seasonality_threshold,
    );

    let modulation_type = match result.modulation_type {
        fdars_core::seasonal::ModulationType::Stable => "stable",
        fdars_core::seasonal::ModulationType::Emerging => "emerging",
        fdars_core::seasonal::ModulationType::Fading => "fading",
        fdars_core::seasonal::ModulationType::Oscillating => "oscillating",
        fdars_core::seasonal::ModulationType::NonSeasonal => "non_seasonal",
        _ => "unknown",
    };

    list!(
        is_seasonal = result.is_seasonal,
        seasonal_strength = result.seasonal_strength,
        has_modulation = result.has_modulation,
        modulation_type = modulation_type,
        modulation_score = result.modulation_score,
        amplitude_trend = result.amplitude_trend,
        strength_curve = result.strength_curve,
        time_points = result.time_points,
        min_strength = result.min_strength,
        max_strength = result.max_strength
    )
    .into()
}

/// Detect amplitude modulation using wavelet transform (Morlet wavelet)
///
/// @param data Matrix of functional data (n x m)
/// @param argvals Vector of evaluation points (length m)
/// @param period Seasonal period in argvals units
/// @param modulation_threshold CV threshold for detecting modulation (default: 0.15)
/// @param seasonality_threshold Strength threshold for seasonality (default: 0.3)
/// @return List with detection results including wavelet amplitude
#[extendr]
fn seasonal_detect_amplitude_modulation_wavelet(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    period: f64,
    modulation_threshold: f64,
    seasonality_threshold: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 || argvals.len() != m || period <= 0.0 {
        return list!(
            is_seasonal = false,
            seasonal_strength = f64::NAN,
            has_modulation = false,
            modulation_type = "unknown",
            modulation_score = f64::NAN,
            amplitude_trend = f64::NAN,
            wavelet_amplitude = Vec::<f64>::new(),
            time_points = Vec::<f64>::new(),
            scale = f64::NAN
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let result = fdars_core::seasonal::detect_amplitude_modulation_wavelet(
        &fd,
        &argvals,
        period,
        modulation_threshold,
        seasonality_threshold,
    );

    let modulation_type = match result.modulation_type {
        fdars_core::seasonal::ModulationType::Stable => "stable",
        fdars_core::seasonal::ModulationType::Emerging => "emerging",
        fdars_core::seasonal::ModulationType::Fading => "fading",
        fdars_core::seasonal::ModulationType::Oscillating => "oscillating",
        fdars_core::seasonal::ModulationType::NonSeasonal => "non_seasonal",
        _ => "unknown",
    };

    list!(
        is_seasonal = result.is_seasonal,
        seasonal_strength = result.seasonal_strength,
        has_modulation = result.has_modulation,
        modulation_type = modulation_type,
        modulation_score = result.modulation_score,
        amplitude_trend = result.amplitude_trend,
        wavelet_amplitude = result.wavelet_amplitude,
        time_points = result.time_points,
        scale = result.scale
    )
    .into()
}

/// CFDAutoperiod: Clustered Filtered Detrended Autoperiod
/// Uses differencing for detrending and clustering for robust period detection
#[extendr]
fn seasonal_cfd_autoperiod(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    cluster_tolerance: f64,
    min_cluster_size: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 8 || argvals.len() != m {
        return list!(
            period = f64::NAN,
            confidence = 0.0,
            acf_validation = 0.0,
            n_periods = 0i32,
            periods = Vec::<f64>::new(),
            confidences = Vec::<f64>::new()
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let tol = if cluster_tolerance > 0.0 {
        Some(cluster_tolerance)
    } else {
        None
    };
    let min_size = if min_cluster_size > 0 {
        Some(min_cluster_size as usize)
    } else {
        None
    };

    let result = fdars_core::cfd_autoperiod_fdata(&fd, &argvals, tol, min_size);

    list!(
        period = result.period,
        confidence = result.confidence,
        acf_validation = result.acf_validation,
        n_periods = result.periods.len() as i32,
        periods = result.periods,
        confidences = result.confidences
    )
    .into()
}

/// Autoperiod: Hybrid FFT + ACF period detection with gradient ascent refinement
/// Returns period, confidence, FFT power, ACF validation score, and candidates
#[extendr]
fn seasonal_autoperiod(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    n_candidates: i32,
    gradient_steps: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 8 || argvals.len() != m {
        return list!(
            period = f64::NAN,
            confidence = 0.0,
            fft_power = 0.0,
            acf_validation = 0.0,
            n_candidates = 0i32,
            candidates = list!()
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let n_cand = if n_candidates > 0 {
        Some(n_candidates as usize)
    } else {
        None
    };
    let steps = if gradient_steps > 0 {
        Some(gradient_steps as usize)
    } else {
        None
    };

    let result = fdars_core::autoperiod_fdata(&fd, &argvals, n_cand, steps);

    // Convert candidates to R list
    let candidate_periods: Vec<f64> = result.candidates.iter().map(|c| c.period).collect();
    let candidate_fft_powers: Vec<f64> = result.candidates.iter().map(|c| c.fft_power).collect();
    let candidate_acf_scores: Vec<f64> = result.candidates.iter().map(|c| c.acf_score).collect();
    let candidate_combined_scores: Vec<f64> =
        result.candidates.iter().map(|c| c.combined_score).collect();

    list!(
        period = result.period,
        confidence = result.confidence,
        fft_power = result.fft_power,
        acf_validation = result.acf_validation,
        n_candidates = result.candidates.len() as i32,
        candidates = list!(
            period = candidate_periods,
            fft_power = candidate_fft_powers,
            acf_score = candidate_acf_scores,
            combined_score = candidate_combined_scores
        )
    )
    .into()
}

/// SAZED: Spectral-ACF Zero-crossing Ensemble Detection
/// A parameter-free ensemble method for robust period detection
/// Returns period, confidence, component periods, and agreeing component count
#[extendr]
fn seasonal_sazed(data: RMatrix<f64>, argvals: Vec<f64>, tolerance: f64) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 8 || argvals.len() != m {
        return list!(
            period = f64::NAN,
            confidence = 0.0,
            agreeing_components = 0i32,
            components = list!(
                spectral = f64::NAN,
                acf_peak = f64::NAN,
                acf_average = f64::NAN,
                zero_crossing = f64::NAN,
                spectral_diff = f64::NAN
            )
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let tol = if tolerance > 0.0 {
        Some(tolerance)
    } else {
        None
    };

    let result = fdars_core::sazed_fdata(&fd, &argvals, tol);

    list!(
        period = result.period,
        confidence = result.confidence,
        agreeing_components = result.agreeing_components as i32,
        components = list!(
            spectral = result.component_periods.spectral,
            acf_peak = result.component_periods.acf_peak,
            acf_average = result.component_periods.acf_average,
            zero_crossing = result.component_periods.zero_crossing,
            spectral_diff = result.component_periods.spectral_diff
        )
    )
    .into()
}

// =============================================================================
// Lomb-Scargle Periodogram
// =============================================================================

/// Lomb-Scargle periodogram for irregularly sampled data
/// Computes the power spectrum and significance for period detection
#[extendr]
fn seasonal_lomb_scargle(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    oversampling: f64,
    nyquist_factor: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 3 || argvals.len() != m {
        return list!(
            frequencies = Vec::<f64>::new(),
            periods = Vec::<f64>::new(),
            power = Vec::<f64>::new(),
            peak_period = f64::NAN,
            peak_frequency = f64::NAN,
            peak_power = f64::NAN,
            false_alarm_probability = 1.0,
            significance = 0.0
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let oversample = if oversampling > 0.0 {
        Some(oversampling)
    } else {
        None
    };
    let nyquist = if nyquist_factor > 0.0 {
        Some(nyquist_factor)
    } else {
        None
    };

    let result =
        fdars_core::seasonal::lomb_scargle_fdata(&fd, &argvals, oversample, nyquist);

    list!(
        frequencies = result.frequencies,
        periods = result.periods,
        power = result.power,
        peak_period = result.peak_period,
        peak_frequency = result.peak_frequency,
        peak_power = result.peak_power,
        false_alarm_probability = result.false_alarm_probability,
        significance = result.significance
    )
    .into()
}

// =============================================================================
// Matrix Profile (STOMP Algorithm)
// =============================================================================

/// Matrix Profile for motif discovery and period detection
/// Uses STOMP algorithm for efficient computation
#[extendr]
fn seasonal_matrix_profile(
    data: RMatrix<f64>,
    subsequence_length: i32,
    exclusion_zone: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 8 {
        return list!(
            profile = Vec::<f64>::new(),
            profile_index = Vec::<i32>::new(),
            subsequence_length = 0i32,
            detected_periods = Vec::<f64>::new(),
            arc_counts = Vec::<i32>::new(),
            primary_period = f64::NAN,
            confidence = 0.0
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let subseq_len = if subsequence_length > 0 {
        Some(subsequence_length as usize)
    } else {
        None
    };

    let exc_zone = if exclusion_zone > 0.0 {
        Some(exclusion_zone)
    } else {
        None
    };

    let result = fdars_core::seasonal::matrix_profile_fdata(&fd, subseq_len, exc_zone);

    // Convert profile_index to i32 (R uses 1-based indexing)
    let profile_index_r: Vec<i32> = result
        .profile_index
        .iter()
        .map(|&i| (i + 1) as i32) // Convert to 1-based
        .collect();

    // Convert arc_counts to i32
    let arc_counts_r: Vec<i32> = result.arc_counts.iter().map(|&c| c as i32).collect();

    list!(
        profile = result.profile,
        profile_index = profile_index_r,
        subsequence_length = result.subsequence_length as i32,
        detected_periods = result.detected_periods,
        arc_counts = arc_counts_r,
        primary_period = result.primary_period,
        confidence = result.confidence
    )
    .into()
}

// =============================================================================
// Singular Spectrum Analysis (SSA)
// =============================================================================

/// Singular Spectrum Analysis for time series decomposition
/// Extracts trend, seasonal, and noise components via SVD
#[extendr]
fn seasonal_ssa(data: RMatrix<f64>, window_length: i32, n_components: i32) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 {
        return list!(
            trend = Vec::<f64>::new(),
            seasonal = Vec::<f64>::new(),
            noise = Vec::<f64>::new(),
            singular_values = Vec::<f64>::new(),
            contributions = Vec::<f64>::new(),
            window_length = 0i32,
            n_components = 0i32,
            detected_period = f64::NAN,
            confidence = 0.0
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let win_len = if window_length > 0 {
        Some(window_length as usize)
    } else {
        None
    };
    let n_comp = if n_components > 0 {
        Some(n_components as usize)
    } else {
        None
    };

    let result = fdars_core::seasonal::ssa_fdata(&fd, win_len, n_comp);

    // Convert trend, seasonal, noise to R matrices
    let trend_mat = RMatrix::new_matrix(n, m, |_r, c| result.trend.get(c).copied().unwrap_or(0.0));
    let seasonal_mat =
        RMatrix::new_matrix(n, m, |_r, c| result.seasonal.get(c).copied().unwrap_or(0.0));
    let noise_mat = RMatrix::new_matrix(n, m, |_r, c| result.noise.get(c).copied().unwrap_or(0.0));

    list!(
        trend = trend_mat,
        seasonal = seasonal_mat,
        noise = noise_mat,
        singular_values = result.singular_values,
        contributions = result.contributions,
        window_length = result.window_length as i32,
        n_components = result.n_components as i32,
        detected_period = result.detected_period,
        confidence = result.confidence
    )
    .into()
}

// =============================================================================
// STL Decomposition
// =============================================================================

/// STL (Seasonal and Trend decomposition using LOESS)
/// Implements Cleveland et al. 1990 algorithm
#[extendr]
fn seasonal_stl(
    data: RMatrix<f64>,
    period: i32,
    s_window: i32,
    t_window: i32,
    robust: bool,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 4 || period < 2 {
        return list!(
            trend = RMatrix::<f64>::new_matrix(0, 0, |_, _| 0.0),
            seasonal = RMatrix::<f64>::new_matrix(0, 0, |_, _| 0.0),
            remainder = RMatrix::<f64>::new_matrix(0, 0, |_, _| 0.0),
            weights = RMatrix::<f64>::new_matrix(0, 0, |_, _| 1.0),
            period = 0i32,
            s_window = 0i32,
            t_window = 0i32,
            inner_iterations = 0i32,
            outer_iterations = 0i32
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let s_win = if s_window > 0 {
        Some(s_window as usize)
    } else {
        None
    };
    let t_win = if t_window > 0 {
        Some(t_window as usize)
    } else {
        None
    };

    let result = fdars_core::detrend::stl_decompose(
        &fd,
        period as usize,
        s_win,
        t_win,
        None,
        robust,
        None,
        None,
    );

    // Convert to R matrices
    let trend_s = result.trend.as_slice();
    let seasonal_s = result.seasonal.as_slice();
    let remainder_s = result.remainder.as_slice();
    let weights_s = result.weights.as_slice();
    let trend_mat = RMatrix::new_matrix(n, m, |r, c| trend_s[r + c * n]);
    let seasonal_mat = RMatrix::new_matrix(n, m, |r, c| seasonal_s[r + c * n]);
    let remainder_mat = RMatrix::new_matrix(n, m, |r, c| remainder_s[r + c * n]);
    let weights_mat = RMatrix::new_matrix(n, m, |r, c| weights_s[r + c * n]);

    list!(
        trend = trend_mat,
        seasonal = seasonal_mat,
        remainder = remainder_mat,
        weights = weights_mat,
        period = result.period as i32,
        s_window = result.s_window as i32,
        t_window = result.t_window as i32,
        inner_iterations = result.inner_iterations as i32,
        outer_iterations = result.outer_iterations as i32
    )
    .into()
}

// =============================================================================
// Detrending functions
// =============================================================================

/// Detrend functional data using specified method
/// Returns trend, detrended data, method used, RSS per curve, and number of parameters
#[extendr]
fn seasonal_detrend(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    method: &str,
    degree: i32,
    bandwidth: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 3 || argvals.len() != m {
        return list!(
            trend = RMatrix::<f64>::new(0, 0),
            detrended = RMatrix::<f64>::new(0, 0),
            method = "none",
            rss = Vec::<f64>::new(),
            n_params = 0i32
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let result = match method {
        "linear" => fdars_core::detrend::detrend_linear(&fd, &argvals),
        "polynomial" => {
            fdars_core::detrend::detrend_polynomial(&fd, &argvals, degree as usize)
        }
        "diff1" => fdars_core::detrend::detrend_diff(&fd, 1),
        "diff2" => fdars_core::detrend::detrend_diff(&fd, 2),
        "loess" => fdars_core::detrend::detrend_loess(&fd, &argvals, bandwidth, 1),
        "auto" => fdars_core::detrend::auto_detrend(&fd, &argvals),
        _ => fdars_core::detrend::detrend_linear(&fd, &argvals),
    };

    // Determine output dimensions based on method
    // v0.4.0 returns full n x m matrices for diff methods (zero-padded),
    // but R wrapper expects truncated n x (m - order) matrices
    let out_m = if method == "diff1" {
        m - 1
    } else if method == "diff2" {
        m - 2
    } else {
        m
    };

    // Convert to R matrices (only first out_m columns for diff methods)
    let trend_s = result.trend.as_slice();
    let detrended_s = result.detrended.as_slice();
    let trend_mat = RMatrix::new_matrix(n, out_m, |r, c| trend_s[r + c * n]);
    let detrended_mat = RMatrix::new_matrix(n, out_m, |r, c| detrended_s[r + c * n]);

    list!(
        trend = trend_mat,
        detrended = detrended_mat,
        method = result.method.as_ref(),
        rss = result.rss,
        n_params = result.n_params as i32
    )
    .into()
}

/// Decompose functional data into trend, seasonal, and remainder components
#[extendr]
fn seasonal_decompose(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    period: f64,
    method: &str,
    trend_method: &str,
    bandwidth: f64,
    n_harmonics: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();

    if n == 0 || m < 3 || argvals.len() != m || period <= 0.0 {
        return list!(
            trend = RMatrix::<f64>::new(0, 0),
            seasonal = RMatrix::<f64>::new(0, 0),
            remainder = RMatrix::<f64>::new(0, 0),
            period = f64::NAN,
            method = "none"
        )
        .into();
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let result = match method {
        "multiplicative" => fdars_core::detrend::decompose_multiplicative(
            &fd,
            &argvals,
            period,
            trend_method,
            bandwidth,
            n_harmonics as usize,
        ),
        "additive" | _ => fdars_core::detrend::decompose_additive(
            &fd,
            &argvals,
            period,
            trend_method,
            bandwidth,
            n_harmonics as usize,
        ),
    };

    // Convert to R matrices
    let trend_s = result.trend.as_slice();
    let seasonal_s = result.seasonal.as_slice();
    let remainder_s = result.remainder.as_slice();
    let trend_mat = RMatrix::new_matrix(n, m, |r, c| trend_s[r + c * n]);
    let seasonal_mat = RMatrix::new_matrix(n, m, |r, c| seasonal_s[r + c * n]);
    let remainder_mat = RMatrix::new_matrix(n, m, |r, c| remainder_s[r + c * n]);

    list!(
        trend = trend_mat,
        seasonal = seasonal_mat,
        remainder = remainder_mat,
        period = result.period,
        method = result.method.as_ref()
    )
    .into()
}

// =============================================================================
// Simulation functions
// =============================================================================

/// Compute eigenfunction basis values
/// efun_type: 0 = Fourier, 1 = Poly, 2 = PolyHigh, 3 = Wiener
#[extendr]
fn eigenfunctions_1d(argvals: Vec<f64>, m: i32, efun_type: i32) -> Robj {
    let m = m as usize;
    let efun = match fdars_core::simulation::EFunType::from_i32(efun_type) {
        Ok(t) => t,
        Err(_) => {
            return Robj::from(RMatrix::<f64>::new_matrix(0, 0, |_, _| 0.0));
        }
    };

    let phi = fdars_core::simulation::eigenfunctions(&argvals, m, efun);
    let n = phi.nrows();
    let mc = phi.ncols();
    let s = phi.as_slice();

    // Return as n x m matrix (column-major)
    let result = RMatrix::new_matrix(n, mc, |i, j| s[i + j * n]);
    Robj::from(result)
}

/// Generate eigenvalue sequence
/// eval_type: 0 = linear, 1 = exponential, 2 = wiener
#[extendr]
fn eigenvalues_1d(m: i32, eval_type: i32) -> Robj {
    let m = m as usize;
    let eval = match fdars_core::simulation::EValType::from_i32(eval_type) {
        Ok(t) => t,
        Err(_) => {
            return Robj::from(Vec::<f64>::new());
        }
    };

    let lambda = fdars_core::simulation::eigenvalues(m, eval);
    Robj::from(lambda)
}

/// Simulate functional data via Karhunen-Loève expansion
#[extendr]
fn sim_kl_1d(n: i32, phi: RMatrix<f64>, lambda: Vec<f64>, seed: Robj) -> Robj {
    let n_samples = n as usize;
    let m = phi.nrows();
    let big_m = phi.ncols();
    let phi_slice = phi.as_real_slice().unwrap();

    let seed_val = if seed.is_null() {
        None
    } else {
        seed.as_integer().map(|s| s as u64)
    };

    let phi_mat = FdMatrix::from_slice(phi_slice, m, big_m).expect("Invalid phi matrix dimensions");
    let data = fdars_core::simulation::sim_kl(n_samples, &phi_mat, big_m, &lambda, seed_val);
    let s = data.as_slice();

    // Return as n x m matrix (column-major)
    let result = RMatrix::new_matrix(data.nrows(), data.ncols(), |i, j| s[i + j * data.nrows()]);
    Robj::from(result)
}

/// Add pointwise Gaussian noise to functional data
#[extendr]
fn add_error_pointwise_1d(data: RMatrix<f64>, sd: f64, seed: Robj) -> Robj {
    let nrow = data.nrows();
    let ncol = data.ncols();
    let data_slice = data.as_real_slice().unwrap();

    let seed_val = if seed.is_null() {
        None
    } else {
        seed.as_integer().map(|s| s as u64)
    };

    let fd = FdMatrix::from_slice(data_slice, nrow, ncol).expect("Invalid matrix dimensions");
    let noisy = fdars_core::simulation::add_error_pointwise(&fd, sd, seed_val);
    let s = noisy.as_slice();

    let result = RMatrix::new_matrix(noisy.nrows(), noisy.ncols(), |i, j| s[i + j * noisy.nrows()]);
    Robj::from(result)
}

/// Add curve-level Gaussian noise to functional data
#[extendr]
fn add_error_curve_1d(data: RMatrix<f64>, sd: f64, seed: Robj) -> Robj {
    let nrow = data.nrows();
    let ncol = data.ncols();
    let data_slice = data.as_real_slice().unwrap();

    let seed_val = if seed.is_null() {
        None
    } else {
        seed.as_integer().map(|s| s as u64)
    };

    let fd = FdMatrix::from_slice(data_slice, nrow, ncol).expect("Invalid matrix dimensions");
    let noisy = fdars_core::simulation::add_error_curve(&fd, sd, seed_val);
    let s = noisy.as_slice();

    let result = RMatrix::new_matrix(noisy.nrows(), noisy.ncols(), |i, j| s[i + j * noisy.nrows()]);
    Robj::from(result)
}

// =============================================================================
// Irregular functional data functions
// =============================================================================

/// Helper to construct IrregFdata from R inputs
fn make_irreg_fdata(offsets: &[i32], argvals: Vec<f64>, values: Vec<f64>) -> Option<fdars_core::irreg_fdata::IrregFdata> {
    let offsets_usize: Vec<usize> = offsets.iter().map(|&x| x as usize).collect();
    // Compute rangeval from argvals
    let range_min = argvals.iter().cloned().fold(f64::INFINITY, f64::min);
    let range_max = argvals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    fdars_core::irreg_fdata::IrregFdata::from_flat(offsets_usize, argvals, values, [range_min, range_max]).ok()
}

/// Compute integral for each curve in irregular functional data
#[extendr]
fn irreg_integrate(offsets: Vec<i32>, argvals: Vec<f64>, values: Vec<f64>) -> Robj {
    let ifd = match make_irreg_fdata(&offsets, argvals, values) {
        Some(ifd) => ifd,
        None => return Robj::from(Vec::<f64>::new()),
    };

    let integrals = fdars_core::irreg_fdata::integrate_irreg(&ifd);
    Robj::from(integrals)
}

/// Compute Lp norm for each curve in irregular functional data
#[extendr]
fn irreg_norm_lp(offsets: Vec<i32>, argvals: Vec<f64>, values: Vec<f64>, p: f64) -> Robj {
    let ifd = match make_irreg_fdata(&offsets, argvals, values) {
        Some(ifd) => ifd,
        None => return Robj::from(Vec::<f64>::new()),
    };

    let norms = fdars_core::irreg_fdata::norm_lp_irreg(&ifd, p);
    Robj::from(norms)
}

/// Estimate mean function for irregular data using kernel smoothing
#[extendr]
fn irreg_mean_kernel(
    offsets: Vec<i32>,
    argvals: Vec<f64>,
    values: Vec<f64>,
    target_argvals: Vec<f64>,
    bandwidth: f64,
    kernel_type: i32,
) -> Robj {
    let ifd = match make_irreg_fdata(&offsets, argvals, values) {
        Some(ifd) => ifd,
        None => return Robj::from(Vec::<f64>::new()),
    };

    let kt = match kernel_type {
        1 => fdars_core::irreg_fdata::KernelType::Gaussian,
        _ => fdars_core::irreg_fdata::KernelType::Epanechnikov,
    };

    let mean = fdars_core::irreg_fdata::mean_irreg(&ifd, &target_argvals, bandwidth, kt);
    Robj::from(mean)
}

/// Compute pairwise Lp distances for irregular functional data
#[extendr]
fn irreg_metric_lp(offsets: Vec<i32>, argvals: Vec<f64>, values: Vec<f64>, p: f64) -> Robj {
    let n = offsets.len() - 1;
    let ifd = match make_irreg_fdata(&offsets, argvals, values) {
        Some(ifd) => ifd,
        None => return Robj::from(RMatrix::new_matrix(0, 0, |_, _| 0.0)),
    };

    let dist = fdars_core::irreg_fdata::metric_lp_irreg(&ifd, p);
    let ds = dist.as_slice();

    let result = RMatrix::new_matrix(n, n, |i, j| ds[i + j * n]);
    Robj::from(result)
}

/// Convert irregular data to regular grid via interpolation
#[extendr]
fn irreg_to_regular(
    offsets: Vec<i32>,
    argvals: Vec<f64>,
    values: Vec<f64>,
    target_grid: Vec<f64>,
) -> Robj {
    let n = offsets.len() - 1;
    let m = target_grid.len();
    let ifd = match make_irreg_fdata(&offsets, argvals, values) {
        Some(ifd) => ifd,
        None => return Robj::from(RMatrix::new_matrix(0, 0, |_, _| 0.0)),
    };

    let data = fdars_core::irreg_fdata::to_regular_grid(&ifd, &target_grid);
    let ds = data.as_slice();

    let result = RMatrix::new_matrix(n, m, |i, j| ds[i + j * n]);
    Robj::from(result)
}

/// Fit basis functions to irregular functional data
/// Each curve is individually fitted via least squares at its own observation points
/// basis_type: 0 = bspline, 1 = fourier
#[extendr]
fn irreg_fdata2basis(
    offsets: Vec<i32>,
    argvals: Vec<f64>,
    values: Vec<f64>,
    nbasis: i32,
    basis_type: i32,
) -> Robj {
    let offsets: Vec<usize> = offsets.iter().map(|&x| x as usize).collect();
    let n = offsets.len() - 1;
    let nbasis = nbasis as usize;

    if n == 0 || nbasis < 2 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }

    // Get overall domain range from all observations
    let t_min = argvals.iter().cloned().fold(f64::INFINITY, f64::min);
    let t_max = argvals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Result matrix: n curves x nbasis coefficients
    let mut coefs = vec![0.0; n * nbasis];

    // Process each curve individually
    for i in 0..n {
        let start = offsets[i];
        let end = offsets[i + 1];
        let m_i = end - start;

        if m_i < 2 {
            // Not enough points to fit basis, leave coefficients as zeros
            continue;
        }

        let curve_argvals = &argvals[start..end];
        let curve_values = &values[start..end];

        // Compute basis matrix for this curve's observation points
        // Scale argvals to domain [t_min, t_max] for consistent basis
        let basis = if basis_type == 1 {
            fourier_basis_with_range(curve_argvals, nbasis, t_min, t_max)
        } else {
            bspline_basis_with_range(curve_argvals, nbasis, t_min, t_max)
        };

        let actual_nbasis = basis.len() / m_i;
        if actual_nbasis == 0 {
            continue;
        }

        // Create nalgebra matrices for least squares
        let b_mat = DMatrix::from_column_slice(m_i, actual_nbasis, &basis);
        let y_vec = DVector::from_column_slice(curve_values);

        // Solve least squares: (B'B)^-1 B' y
        let btb = &b_mat.transpose() * &b_mat;
        let bty = &b_mat.transpose() * &y_vec;

        // Use SVD for pseudo-inverse (stable for potentially ill-conditioned systems)
        let btb_svd = SVD::new(btb.clone(), true, true);

        let max_sv = btb_svd.singular_values.iter().cloned().fold(0.0, f64::max);
        let eps = 1e-10 * max_sv.max(1e-10);

        // Pseudo-inverse of singular values
        let s_inv: Vec<f64> = btb_svd
            .singular_values
            .iter()
            .map(|&s| if s > eps { 1.0 / s } else { 0.0 })
            .collect();

        // Compute pseudo-inverse: V * S_inv * U^T
        if let (Some(ref v_t), Some(ref u)) = (&btb_svd.v_t, &btb_svd.u) {
            let v = v_t.transpose();
            let u_t = u.transpose();

            // btb_inv = V * diag(s_inv) * U^T
            let mut btb_inv = DMatrix::zeros(actual_nbasis, actual_nbasis);
            for ii in 0..actual_nbasis {
                for jj in 0..actual_nbasis {
                    let mut sum = 0.0;
                    for k in 0..actual_nbasis {
                        sum += v[(ii, k)] * s_inv[k] * u_t[(k, jj)];
                    }
                    btb_inv[(ii, jj)] = sum;
                }
            }

            // Compute coefficients: btb_inv * bty
            let curve_coefs = btb_inv * bty;

            // Store coefficients (column-major for R)
            for k in 0..actual_nbasis.min(nbasis) {
                coefs[i + k * n] = curve_coefs[k];
            }
        }
    }

    let result = RMatrix::new_matrix(n, nbasis, |i, j| coefs[i + j * n]);
    Robj::from(result)
}

/// B-spline basis with specified range (for irregular data with common domain)
fn bspline_basis_with_range(t: &[f64], nbasis: usize, t_min: f64, t_max: f64) -> Vec<f64> {
    let n = t.len();
    let order = 4; // Cubic B-splines
    let nknots = (nbasis.saturating_sub(order)).max(2);
    let actual_nbasis = nknots + order;

    let dt = (t_max - t_min) / (nknots - 1).max(1) as f64;

    let mut knots = Vec::with_capacity(nknots + 2 * order);
    // Pad at beginning
    for i in 0..order {
        knots.push(t_min - (order - i) as f64 * dt);
    }
    // Interior knots
    for i in 0..nknots {
        knots.push(t_min + i as f64 * dt);
    }
    // Pad at end
    for i in 1..=order {
        knots.push(t_max + i as f64 * dt);
    }

    let t_max_knot_idx = order + nknots - 1;

    // Cox-de Boor recursion
    let mut basis = vec![0.0; n * actual_nbasis];

    for (ti, &t_val) in t.iter().enumerate() {
        let mut b0 = vec![0.0; knots.len() - 1];
        for j in 0..(knots.len() - 1) {
            let in_interval = if j == t_max_knot_idx - 1 {
                t_val >= knots[j] && t_val <= knots[j + 1]
            } else {
                t_val >= knots[j] && t_val < knots[j + 1]
            };
            if in_interval {
                b0[j] = 1.0;
            }
        }

        let mut b_prev = b0;
        for k in 1..order {
            let mut b_next = vec![0.0; knots.len() - 1 - k];
            for j in 0..b_next.len() {
                let left_denom = knots[j + k] - knots[j];
                let left = if left_denom.abs() > 1e-15 {
                    (t_val - knots[j]) / left_denom * b_prev[j]
                } else {
                    0.0
                };

                let right_denom = knots[j + k + 1] - knots[j + 1];
                let right = if right_denom.abs() > 1e-15 {
                    (knots[j + k + 1] - t_val) / right_denom * b_prev[j + 1]
                } else {
                    0.0
                };

                b_next[j] = left + right;
            }
            b_prev = b_next;
        }

        for (j, &val) in b_prev.iter().enumerate().take(actual_nbasis) {
            basis[ti + j * n] = val;
        }
    }

    basis
}

/// Fourier basis with specified range (for irregular data with common domain)
fn fourier_basis_with_range(t: &[f64], nbasis: usize, t_min: f64, t_max: f64) -> Vec<f64> {
    let n = t.len();
    let period = t_max - t_min;

    let mut basis = vec![0.0; n * nbasis];

    for (i, &ti) in t.iter().enumerate() {
        let x = 2.0 * PI * (ti - t_min) / period;

        // First basis function is constant
        basis[i] = 1.0;

        // Remaining basis functions are sin/cos pairs
        let mut k = 1;
        let mut freq = 1;
        while k < nbasis {
            if k < nbasis {
                basis[i + k * n] = (freq as f64 * x).sin();
                k += 1;
            }
            if k < nbasis {
                basis[i + k * n] = (freq as f64 * x).cos();
                k += 1;
            }
            freq += 1;
        }
    }

    basis
}

// =============================================================================
// Alignment functions
// =============================================================================

/// SRSF transform of functional data
#[extendr]
fn alignment_srsf_transform(data: RMatrix<f64>, argvals: Vec<f64>) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let result = alignment::srsf_transform(&fd, &argvals);
    let s = result.as_slice();
    let mat = RMatrix::new_matrix(result.nrows(), result.ncols(), |i, j| s[i + j * result.nrows()]);
    Robj::from(mat)
}

/// Inverse SRSF: reconstruct curve from SRSF representation
#[extendr]
fn alignment_srsf_inverse(q: Vec<f64>, argvals: Vec<f64>, f0: f64) -> Robj {
    if q.len() != argvals.len() || q.is_empty() {
        return Robj::from(Vec::<f64>::new());
    }
    Robj::from(alignment::srsf_inverse(&q, &argvals, f0))
}

/// Apply warping function to reparameterize a curve
#[extendr]
fn alignment_reparameterize(f: Vec<f64>, argvals: Vec<f64>, gamma: Vec<f64>) -> Robj {
    if f.len() != argvals.len() || f.len() != gamma.len() || f.is_empty() {
        return Robj::from(Vec::<f64>::new());
    }
    Robj::from(alignment::reparameterize_curve(&f, &argvals, &gamma))
}

/// Compose two warping functions
#[extendr]
fn alignment_compose_warps(gamma1: Vec<f64>, gamma2: Vec<f64>, argvals: Vec<f64>) -> Robj {
    if gamma1.len() != argvals.len() || gamma2.len() != argvals.len() || argvals.is_empty() {
        return Robj::from(Vec::<f64>::new());
    }
    Robj::from(alignment::compose_warps(&gamma1, &gamma2, &argvals))
}

/// Elastic alignment of one curve to another
#[extendr]
fn alignment_elastic_pair(f1: Vec<f64>, f2: Vec<f64>, argvals: Vec<f64>) -> Robj {
    if f1.len() != argvals.len() || f2.len() != argvals.len() || argvals.len() < 2 {
        return list!(gamma = Vec::<f64>::new(), f_aligned = Vec::<f64>::new(), distance = f64::NAN).into();
    }
    let res = alignment::elastic_align_pair(&f1, &f2, &argvals, 0.0);
    list!(gamma = res.gamma, f_aligned = res.f_aligned, distance = res.distance).into()
}

/// Elastic (Fisher-Rao) distance between two curves
#[extendr]
fn alignment_elastic_distance(f1: Vec<f64>, f2: Vec<f64>, argvals: Vec<f64>) -> Robj {
    if f1.len() != argvals.len() || f2.len() != argvals.len() || argvals.len() < 2 {
        return Robj::from(f64::NAN);
    }
    Robj::from(alignment::elastic_distance(&f1, &f2, &argvals, 0.0))
}

/// Align all curves to a target curve
#[extendr]
fn alignment_align_to_target(data: RMatrix<f64>, target: Vec<f64>, argvals: Vec<f64>) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m || target.len() != m {
        return list!(
            gammas = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            aligned_data = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            distances = Vec::<f64>::new()
        ).into();
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let res = alignment::align_to_target(&fd, &target, &argvals, 0.0);
    let gs = res.gammas.as_slice();
    let as_ = res.aligned_data.as_slice();
    let g_mat = RMatrix::new_matrix(res.gammas.nrows(), res.gammas.ncols(), |i, j| gs[i + j * res.gammas.nrows()]);
    let a_mat = RMatrix::new_matrix(res.aligned_data.nrows(), res.aligned_data.ncols(), |i, j| as_[i + j * res.aligned_data.nrows()]);
    list!(gammas = g_mat, aligned_data = a_mat, distances = res.distances).into()
}

/// Elastic self-distance matrix
#[extendr]
fn alignment_self_dist(data: RMatrix<f64>, argvals: Vec<f64>) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let result = alignment::elastic_self_distance_matrix(&fd, &argvals, 0.0);
    let s = result.as_slice();
    let mat = RMatrix::new_matrix(result.nrows(), result.ncols(), |i, j| s[i + j * result.nrows()]);
    Robj::from(mat)
}

/// Elastic cross-distance matrix
#[extendr]
fn alignment_cross_dist(data1: RMatrix<f64>, data2: RMatrix<f64>, argvals: Vec<f64>) -> Robj {
    let n1 = data1.nrows();
    let m1 = data1.ncols();
    let n2 = data2.nrows();
    let m2 = data2.ncols();
    if n1 == 0 || n2 == 0 || m1 < 2 || m1 != m2 || argvals.len() != m1 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }
    let d1 = data1.as_real_slice().unwrap();
    let d2 = data2.as_real_slice().unwrap();
    let fd1 = FdMatrix::from_slice(d1, n1, m1).expect("Invalid matrix dimensions");
    let fd2 = FdMatrix::from_slice(d2, n2, m2).expect("Invalid matrix dimensions");
    let result = alignment::elastic_cross_distance_matrix(&fd1, &fd2, &argvals, 0.0);
    let s = result.as_slice();
    let mat = RMatrix::new_matrix(result.nrows(), result.ncols(), |i, j| s[i + j * result.nrows()]);
    Robj::from(mat)
}

/// Karcher (Fréchet) mean in elastic metric
#[extendr]
fn alignment_karcher_mean(data: RMatrix<f64>, argvals: Vec<f64>, max_iter: i32, tol: f64) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return list!(
            mean = Vec::<f64>::new(),
            mean_srsf = Vec::<f64>::new(),
            gammas = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            aligned_data = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            n_iter = 0i32,
            converged = false
        ).into();
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let res = alignment::karcher_mean(&fd, &argvals, max_iter as usize, tol, 0.0);
    let gs = res.gammas.as_slice();
    let as_ = res.aligned_data.as_slice();
    let g_mat = RMatrix::new_matrix(res.gammas.nrows(), res.gammas.ncols(), |i, j| gs[i + j * res.gammas.nrows()]);
    let a_mat = RMatrix::new_matrix(res.aligned_data.nrows(), res.aligned_data.ncols(), |i, j| as_[i + j * res.aligned_data.nrows()]);
    list!(
        mean = res.mean,
        mean_srsf = res.mean_srsf,
        gammas = g_mat,
        aligned_data = a_mat,
        n_iter = res.n_iter as i32,
        converged = res.converged
    ).into()
}

// =============================================================================
// Soft-DTW metric functions
// =============================================================================

/// Soft-DTW self-distance matrix
#[extendr]
fn metric_soft_dtw_self_1d(fdata: RMatrix<f64>, gamma: f64) -> Robj {
    let n = fdata.nrows();
    let m = fdata.ncols();
    if n == 0 || m == 0 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }
    let data_slice = fdata.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let result = core_metric::soft_dtw_self_1d(&fd, gamma);
    let s = result.as_slice();
    let mat = RMatrix::new_matrix(result.nrows(), result.ncols(), |i, j| s[i + j * result.nrows()]);
    Robj::from(mat)
}

/// Soft-DTW cross-distance matrix
#[extendr]
fn metric_soft_dtw_cross_1d(fdata1: RMatrix<f64>, fdata2: RMatrix<f64>, gamma: f64) -> Robj {
    let n1 = fdata1.nrows();
    let m1 = fdata1.ncols();
    let n2 = fdata2.nrows();
    let m2 = fdata2.ncols();
    if n1 == 0 || n2 == 0 || m1 == 0 || m2 == 0 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }
    let d1 = fdata1.as_real_slice().unwrap();
    let d2 = fdata2.as_real_slice().unwrap();
    let fd1 = FdMatrix::from_slice(d1, n1, m1).expect("Invalid matrix dimensions");
    let fd2 = FdMatrix::from_slice(d2, n2, m2).expect("Invalid matrix dimensions");
    let result = core_metric::soft_dtw_cross_1d(&fd1, &fd2, gamma);
    let s = result.as_slice();
    let mat = RMatrix::new_matrix(result.nrows(), result.ncols(), |i, j| s[i + j * result.nrows()]);
    Robj::from(mat)
}

/// Soft-DTW divergence self-distance matrix
#[extendr]
fn metric_soft_dtw_div_self_1d(fdata: RMatrix<f64>, gamma: f64) -> Robj {
    let n = fdata.nrows();
    let m = fdata.ncols();
    if n == 0 || m == 0 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }
    let data_slice = fdata.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let result = core_metric::soft_dtw_div_self_1d(&fd, gamma);
    let s = result.as_slice();
    let mat = RMatrix::new_matrix(result.nrows(), result.ncols(), |i, j| s[i + j * result.nrows()]);
    Robj::from(mat)
}

/// Soft-DTW divergence cross-distance matrix
#[extendr]
fn metric_soft_dtw_div_cross_1d(fdata1: RMatrix<f64>, fdata2: RMatrix<f64>, gamma: f64) -> Robj {
    let n1 = fdata1.nrows();
    let m1 = fdata1.ncols();
    let n2 = fdata2.nrows();
    let m2 = fdata2.ncols();
    if n1 == 0 || n2 == 0 || m1 == 0 || m2 == 0 {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }
    let d1 = fdata1.as_real_slice().unwrap();
    let d2 = fdata2.as_real_slice().unwrap();
    let fd1 = FdMatrix::from_slice(d1, n1, m1).expect("Invalid matrix dimensions");
    let fd2 = FdMatrix::from_slice(d2, n2, m2).expect("Invalid matrix dimensions");
    let result = core_metric::soft_dtw_div_cross_1d(&fd1, &fd2, gamma);
    let s = result.as_slice();
    let mat = RMatrix::new_matrix(result.nrows(), result.ncols(), |i, j| s[i + j * result.nrows()]);
    Robj::from(mat)
}

/// Soft-DTW barycenter computation
#[extendr]
fn metric_soft_dtw_barycenter(fdata: RMatrix<f64>, gamma: f64, max_iter: i32, tol: f64) -> Robj {
    let n = fdata.nrows();
    let m = fdata.ncols();
    if n == 0 || m == 0 {
        return list!(barycenter = Vec::<f64>::new(), n_iter = 0i32, converged = true).into();
    }
    let data_slice = fdata.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let result = core_metric::soft_dtw_barycenter(&fd, gamma, max_iter as usize, tol);
    list!(barycenter = result.barycenter, n_iter = result.n_iter as i32, converged = result.converged).into()
}

// =============================================================================
// Landmark registration functions
// =============================================================================

/// Detect landmarks in a single curve
#[extendr]
fn landmark_detect(curve: Vec<f64>, argvals: Vec<f64>, kind: &str, min_prominence: f64) -> Robj {
    if curve.len() != argvals.len() || curve.len() < 3 {
        return list!(position = Vec::<f64>::new(), kind = Vec::<String>::new(), value = Vec::<f64>::new(), prominence = Vec::<f64>::new()).into();
    }
    let lk = match kind {
        "peak" => landmark::LandmarkKind::Peak,
        "valley" => landmark::LandmarkKind::Valley,
        "zero" | "zero_crossing" => landmark::LandmarkKind::ZeroCrossing,
        "inflection" => landmark::LandmarkKind::Inflection,
        _ => landmark::LandmarkKind::Peak,
    };
    let lms = landmark::detect_landmarks(&curve, &argvals, lk, min_prominence);
    let positions: Vec<f64> = lms.iter().map(|l| l.position).collect();
    let kinds: Vec<String> = lms.iter().map(|l| match l.kind {
        landmark::LandmarkKind::Peak => "peak".to_string(),
        landmark::LandmarkKind::Valley => "valley".to_string(),
        landmark::LandmarkKind::ZeroCrossing => "zero_crossing".to_string(),
        landmark::LandmarkKind::Inflection => "inflection".to_string(),
        landmark::LandmarkKind::Custom => "custom".to_string(),
        _ => "unknown".to_string(),
    }).collect();
    let values: Vec<f64> = lms.iter().map(|l| l.value).collect();
    let prominences: Vec<f64> = lms.iter().map(|l| l.prominence).collect();
    list!(position = positions, kind = kinds, value = values, prominence = prominences).into()
}

/// Detect landmarks and register curves
#[extendr]
fn landmark_register_curves(data: RMatrix<f64>, argvals: Vec<f64>, kind: &str, min_prominence: f64, expected_count: i32) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 3 || argvals.len() != m {
        return list!(
            registered = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            gammas = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            target_landmarks = Vec::<f64>::new(),
            n_landmarks = Vec::<i32>::new()
        ).into();
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let lk = match kind {
        "peak" => landmark::LandmarkKind::Peak,
        "valley" => landmark::LandmarkKind::Valley,
        "zero" | "zero_crossing" => landmark::LandmarkKind::ZeroCrossing,
        "inflection" => landmark::LandmarkKind::Inflection,
        _ => landmark::LandmarkKind::Peak,
    };
    let result = landmark::detect_and_register(&fd, &argvals, lk, min_prominence, expected_count as usize);
    let rs = result.registered.as_slice();
    let gs = result.gammas.as_slice();
    let r_mat = RMatrix::new_matrix(result.registered.nrows(), result.registered.ncols(), |i, j| rs[i + j * result.registered.nrows()]);
    let g_mat = RMatrix::new_matrix(result.gammas.nrows(), result.gammas.ncols(), |i, j| gs[i + j * result.gammas.nrows()]);
    let n_landmarks: Vec<i32> = result.landmarks.iter().map(|lms| lms.len() as i32).collect();
    list!(registered = r_mat, gammas = g_mat, target_landmarks = result.target_landmarks, n_landmarks = n_landmarks).into()
}

// =============================================================================
// TSRVF functions
// =============================================================================

/// Full TSRVF transform: compute Karcher mean + transport to tangent space
#[extendr]
fn alignment_tsrvf_transform(data: RMatrix<f64>, argvals: Vec<f64>, max_iter: i32, tol: f64, lambda: f64) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return list!(
            tangent_vectors = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            mean = Vec::<f64>::new(),
            mean_srsf = Vec::<f64>::new(),
            gammas = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            converged = false
        ).into();
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let res = alignment::tsrvf_transform(&fd, &argvals, max_iter as usize, tol, lambda);
    let tv_s = res.tangent_vectors.as_slice();
    let g_s = res.gammas.as_slice();
    let tv_mat = RMatrix::new_matrix(res.tangent_vectors.nrows(), res.tangent_vectors.ncols(), |i, j| tv_s[i + j * res.tangent_vectors.nrows()]);
    let g_mat = RMatrix::new_matrix(res.gammas.nrows(), res.gammas.ncols(), |i, j| g_s[i + j * res.gammas.nrows()]);
    list!(
        tangent_vectors = tv_mat,
        mean = res.mean,
        mean_srsf = res.mean_srsf,
        mean_srsf_norm = res.mean_srsf_norm,
        srsf_norms = res.srsf_norms,
        gammas = g_mat,
        converged = res.converged
    ).into()
}

/// Compute TSRVF from a pre-computed Karcher mean
#[extendr]
fn alignment_tsrvf_from_karcher(
    mean: Vec<f64>, mean_srsf: Vec<f64>,
    gammas: RMatrix<f64>, aligned_data: RMatrix<f64>,
    argvals: Vec<f64>, n_iter: i32, converged: bool
) -> Robj {
    let n = aligned_data.nrows();
    let m = aligned_data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return list!(
            tangent_vectors = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            mean = Vec::<f64>::new(),
            mean_srsf = Vec::<f64>::new(),
            gammas = RMatrix::new_matrix(0, 0, |_, _| 0.0),
            converged = false
        ).into();
    }
    // Reconstruct KarcherMeanResult
    let g_slice = gammas.as_real_slice().unwrap();
    let a_slice = aligned_data.as_real_slice().unwrap();
    let g_fd = FdMatrix::from_slice(g_slice, n, m).expect("Invalid gammas dimensions");
    let a_fd = FdMatrix::from_slice(a_slice, n, m).expect("Invalid aligned_data dimensions");
    let karcher = alignment::KarcherMeanResult::new(
        mean.clone(), mean_srsf.clone(), g_fd, a_fd,
        n_iter as usize, converged, None,
    );
    let res = alignment::tsrvf_from_alignment(&karcher, &argvals);
    let tv_s = res.tangent_vectors.as_slice();
    let tv_mat = RMatrix::new_matrix(res.tangent_vectors.nrows(), res.tangent_vectors.ncols(), |i, j| tv_s[i + j * res.tangent_vectors.nrows()]);
    list!(
        tangent_vectors = tv_mat,
        mean = res.mean,
        mean_srsf = res.mean_srsf,
        mean_srsf_norm = res.mean_srsf_norm,
        srsf_norms = res.srsf_norms,
        gammas = gammas,
        converged = res.converged
    ).into()
}

/// Inverse TSRVF: reconstruct curves from tangent vectors
#[extendr]
fn alignment_tsrvf_inverse(
    tangent_vectors: RMatrix<f64>, _mean: Vec<f64>,
    mean_srsf: Vec<f64>, mean_srsf_norm: f64,
    srsf_norms: Vec<f64>, _gammas: RMatrix<f64>,
    argvals: Vec<f64>, _converged: bool
) -> Robj {
    let n = tangent_vectors.nrows();
    let m = tangent_vectors.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }
    let tv_slice = tangent_vectors.as_real_slice().unwrap();
    let tv_fd = FdMatrix::from_slice(tv_slice, n, m).expect("Invalid tangent_vectors dimensions");

    // Inline tsrvf_inverse logic since TsrvfResult is non-exhaustive
    let time: Vec<f64> = (0..m).map(|i| i as f64 / (m - 1) as f64).collect();
    let mu_unit: Vec<f64> = if mean_srsf_norm > 1e-10 {
        mean_srsf.iter().map(|&q| q / mean_srsf_norm).collect()
    } else {
        vec![0.0; m]
    };
    let initial_values = vec![0.0; n];
    let curves: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let vi = tv_fd.row(i);
            let qi_unit = fdars_core::warping::exp_map_sphere(&mu_unit, &vi, &time);
            let qi: Vec<f64> = qi_unit.iter().map(|&q| q * srsf_norms[i]).collect();
            alignment::srsf_inverse(&qi, &argvals, initial_values[i])
        })
        .collect();
    let mut result = FdMatrix::zeros(n, m);
    for i in 0..n {
        for j in 0..m {
            result[(i, j)] = curves[i][j];
        }
    }
    let s = result.as_slice();
    let mat = RMatrix::new_matrix(result.nrows(), result.ncols(), |i, j| s[i + j * result.nrows()]);
    Robj::from(mat)
}

// =============================================================================
// Alignment quality functions
// =============================================================================

/// Compute alignment quality metrics
#[extendr]
fn alignment_quality_compute(
    data: RMatrix<f64>, mean: Vec<f64>, mean_srsf: Vec<f64>,
    gammas: RMatrix<f64>, aligned_data: RMatrix<f64>,
    argvals: Vec<f64>, n_iter: i32, converged: bool
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return list!(
            warp_complexity = Vec::<f64>::new(),
            mean_warp_complexity = 0.0,
            warp_smoothness = Vec::<f64>::new(),
            mean_warp_smoothness = 0.0,
            total_variance = 0.0,
            amplitude_variance = 0.0,
            phase_variance = 0.0,
            phase_amplitude_ratio = 0.0,
            pointwise_variance_ratio = Vec::<f64>::new(),
            mean_variance_reduction = 0.0
        ).into();
    }
    let data_slice = data.as_real_slice().unwrap();
    let g_slice = gammas.as_real_slice().unwrap();
    let a_slice = aligned_data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid data dimensions");
    let g_fd = FdMatrix::from_slice(g_slice, n, m).expect("Invalid gammas dimensions");
    let a_fd = FdMatrix::from_slice(a_slice, n, m).expect("Invalid aligned_data dimensions");
    let karcher = alignment::KarcherMeanResult::new(
        mean, mean_srsf, g_fd, a_fd,
        n_iter as usize, converged, None,
    );
    let q = alignment::alignment_quality(&fd, &karcher, &argvals);
    list!(
        warp_complexity = q.warp_complexity,
        mean_warp_complexity = q.mean_warp_complexity,
        warp_smoothness = q.warp_smoothness,
        mean_warp_smoothness = q.mean_warp_smoothness,
        total_variance = q.total_variance,
        amplitude_variance = q.amplitude_variance,
        phase_variance = q.phase_variance,
        phase_amplitude_ratio = q.phase_amplitude_ratio,
        pointwise_variance_ratio = q.pointwise_variance_ratio,
        mean_variance_reduction = q.mean_variance_reduction
    ).into()
}

/// Compute warp complexity (geodesic distance from identity)
#[extendr]
fn alignment_warp_complexity(gamma: Vec<f64>, argvals: Vec<f64>) -> f64 {
    if gamma.len() != argvals.len() || gamma.len() < 2 {
        return 0.0;
    }
    alignment::warp_complexity(&gamma, &argvals)
}

/// Compute warp smoothness (bending energy)
#[extendr]
fn alignment_warp_smoothness(gamma: Vec<f64>, argvals: Vec<f64>) -> f64 {
    if gamma.len() != argvals.len() || gamma.len() < 3 {
        return 0.0;
    }
    alignment::warp_smoothness(&gamma, &argvals)
}

// =============================================================================
// Decomposition and distance functions
// =============================================================================

/// Elastic phase-amplitude decomposition
#[extendr]
fn alignment_decomposition(f1: Vec<f64>, f2: Vec<f64>, argvals: Vec<f64>, lambda: f64) -> Robj {
    if f1.len() != argvals.len() || f2.len() != argvals.len() || argvals.len() < 2 {
        return list!(gamma = Vec::<f64>::new(), f_aligned = Vec::<f64>::new(), distance = f64::NAN, d_amplitude = f64::NAN, d_phase = f64::NAN).into();
    }
    let res = alignment::elastic_decomposition(&f1, &f2, &argvals, lambda);
    list!(
        gamma = res.alignment.gamma,
        f_aligned = res.alignment.f_aligned,
        distance = res.alignment.distance,
        d_amplitude = res.d_amplitude,
        d_phase = res.d_phase
    ).into()
}

/// Amplitude self-distance matrix
#[extendr]
fn alignment_amplitude_dist(data: RMatrix<f64>, argvals: Vec<f64>, lambda: f64) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let result = alignment::amplitude_self_distance_matrix(&fd, &argvals, lambda);
    let s = result.as_slice();
    let mat = RMatrix::new_matrix(result.nrows(), result.ncols(), |i, j| s[i + j * result.nrows()]);
    Robj::from(mat)
}

/// Phase self-distance matrix
#[extendr]
fn alignment_phase_dist(data: RMatrix<f64>, argvals: Vec<f64>, lambda: f64) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let result = alignment::phase_self_distance_matrix(&fd, &argvals, lambda);
    let s = result.as_slice();
    let mat = RMatrix::new_matrix(result.nrows(), result.ncols(), |i, j| s[i + j * result.nrows()]);
    Robj::from(mat)
}

/// Pairwise alignment consistency
#[extendr]
fn alignment_pairwise_consistency(data: RMatrix<f64>, argvals: Vec<f64>, lambda: f64, max_triplets: i32) -> f64 {
    let n = data.nrows();
    let m = data.ncols();
    if n < 3 || m < 2 || argvals.len() != m {
        return 0.0;
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    alignment::pairwise_consistency(&fd, &argvals, lambda, max_triplets as usize)
}

// =============================================================================
// Constrained alignment functions
// =============================================================================

/// Elastic alignment with explicit landmark constraints
#[extendr]
fn alignment_constrained(f1: Vec<f64>, f2: Vec<f64>, argvals: Vec<f64>, landmark_targets: Vec<f64>, landmark_sources: Vec<f64>, lambda: f64) -> Robj {
    if f1.len() != argvals.len() || f2.len() != argvals.len() || argvals.len() < 2 {
        return list!(gamma = Vec::<f64>::new(), f_aligned = Vec::<f64>::new(), distance = f64::NAN, enforced_landmarks_target = Vec::<f64>::new(), enforced_landmarks_source = Vec::<f64>::new()).into();
    }
    let pairs: Vec<(f64, f64)> = landmark_targets.iter().zip(landmark_sources.iter()).map(|(&t, &s)| (t, s)).collect();
    let res = alignment::elastic_align_pair_constrained(&f1, &f2, &argvals, &pairs, lambda);
    let et: Vec<f64> = res.enforced_landmarks.iter().map(|&(t, _)| t).collect();
    let es: Vec<f64> = res.enforced_landmarks.iter().map(|&(_, s)| s).collect();
    list!(gamma = res.gamma, f_aligned = res.f_aligned, distance = res.distance, enforced_landmarks_target = et, enforced_landmarks_source = es).into()
}

/// Elastic alignment with automatic landmark detection
#[extendr]
fn alignment_with_landmarks(f1: Vec<f64>, f2: Vec<f64>, argvals: Vec<f64>, kind: &str, min_prominence: f64, expected_count: i32, lambda: f64) -> Robj {
    if f1.len() != argvals.len() || f2.len() != argvals.len() || argvals.len() < 2 {
        return list!(gamma = Vec::<f64>::new(), f_aligned = Vec::<f64>::new(), distance = f64::NAN, enforced_landmarks_target = Vec::<f64>::new(), enforced_landmarks_source = Vec::<f64>::new()).into();
    }
    let lk = match kind {
        "peak" => landmark::LandmarkKind::Peak,
        "valley" => landmark::LandmarkKind::Valley,
        "zero" | "zero_crossing" => landmark::LandmarkKind::ZeroCrossing,
        "inflection" => landmark::LandmarkKind::Inflection,
        _ => landmark::LandmarkKind::Peak,
    };
    let res = alignment::elastic_align_pair_with_landmarks(&f1, &f2, &argvals, lk, min_prominence, expected_count as usize, lambda);
    let et: Vec<f64> = res.enforced_landmarks.iter().map(|&(t, _)| t).collect();
    let es: Vec<f64> = res.enforced_landmarks.iter().map(|&(_, s)| s).collect();
    list!(gamma = res.gamma, f_aligned = res.f_aligned, distance = res.distance, enforced_landmarks_target = et, enforced_landmarks_source = es).into()
}

// =============================================================================
// Tolerance band functions
// =============================================================================

/// FPCA-based tolerance band
#[extendr]
fn tolerance_fpca(data: RMatrix<f64>, ncomp: i32, nb: i32, coverage: f64, band_type: &str, seed: Robj) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let seed_val = if seed.is_null() { 42u64 } else { seed.as_integer().unwrap_or(42) as u64 };
    let bt = match band_type {
        "simultaneous" => tolerance::BandType::Simultaneous,
        _ => tolerance::BandType::Pointwise,
    };
    if n == 0 || m == 0 {
        return list!(lower = Vec::<f64>::new(), upper = Vec::<f64>::new(), center = Vec::<f64>::new(), half_width = Vec::<f64>::new(), success = false).into();
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    match tolerance::fpca_tolerance_band(&fd, ncomp as usize, nb as usize, coverage, bt, seed_val) {
        Ok(band) => list!(lower = band.lower, upper = band.upper, center = band.center, half_width = band.half_width, success = true).into(),
        Err(_) => list!(lower = Vec::<f64>::new(), upper = Vec::<f64>::new(), center = Vec::<f64>::new(), half_width = Vec::<f64>::new(), success = false).into(),
    }
}

/// Conformal prediction band
#[extendr]
fn tolerance_conformal(data: RMatrix<f64>, cal_fraction: f64, coverage: f64, score_type: &str, seed: Robj) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let seed_val = if seed.is_null() { 42u64 } else { seed.as_integer().unwrap_or(42) as u64 };
    let st = match score_type {
        "l2" | "L2" => tolerance::NonConformityScore::L2,
        _ => tolerance::NonConformityScore::SupNorm,
    };
    if n == 0 || m == 0 {
        return list!(lower = Vec::<f64>::new(), upper = Vec::<f64>::new(), center = Vec::<f64>::new(), half_width = Vec::<f64>::new(), success = false).into();
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    match tolerance::conformal_prediction_band(&fd, cal_fraction, coverage, st, seed_val) {
        Some(band) => list!(lower = band.lower, upper = band.upper, center = band.center, half_width = band.half_width, success = true).into(),
        None => list!(lower = Vec::<f64>::new(), upper = Vec::<f64>::new(), center = Vec::<f64>::new(), half_width = Vec::<f64>::new(), success = false).into(),
    }
}

/// SCB mean confidence band (Degras method)
#[extendr]
fn tolerance_scb_degras(data: RMatrix<f64>, argvals: Vec<f64>, bandwidth: f64, nb: i32, confidence: f64, multiplier: &str) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let mult = match multiplier {
        "rademacher" | "Rademacher" => tolerance::MultiplierDistribution::Rademacher,
        _ => tolerance::MultiplierDistribution::Gaussian,
    };
    if n == 0 || m == 0 || argvals.len() != m {
        return list!(lower = Vec::<f64>::new(), upper = Vec::<f64>::new(), center = Vec::<f64>::new(), half_width = Vec::<f64>::new(), success = false).into();
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    match tolerance::scb_mean_degras(&fd, &argvals, bandwidth, nb as usize, confidence, mult) {
        Ok(band) => list!(lower = band.lower, upper = band.upper, center = band.center, half_width = band.half_width, success = true).into(),
        Err(_) => list!(lower = Vec::<f64>::new(), upper = Vec::<f64>::new(), center = Vec::<f64>::new(), half_width = Vec::<f64>::new(), success = false).into(),
    }
}

/// Exponential family tolerance band
#[extendr]
fn tolerance_exponential(data: RMatrix<f64>, family: &str, ncomp: i32, nb: i32, coverage: f64, seed: Robj) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let seed_val = if seed.is_null() { 42u64 } else { seed.as_integer().unwrap_or(42) as u64 };
    let fam = match family {
        "binomial" => tolerance::ExponentialFamily::Binomial,
        "poisson" => tolerance::ExponentialFamily::Poisson,
        _ => tolerance::ExponentialFamily::Gaussian,
    };
    if n == 0 || m == 0 {
        return list!(lower = Vec::<f64>::new(), upper = Vec::<f64>::new(), center = Vec::<f64>::new(), half_width = Vec::<f64>::new(), success = false).into();
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    match tolerance::exponential_family_tolerance_band(&fd, fam, ncomp as usize, nb as usize, coverage, seed_val) {
        Ok(band) => list!(lower = band.lower, upper = band.upper, center = band.center, half_width = band.half_width, success = true).into(),
        Err(_) => list!(lower = Vec::<f64>::new(), upper = Vec::<f64>::new(), center = Vec::<f64>::new(), half_width = Vec::<f64>::new(), success = false).into(),
    }
}

/// Elastic tolerance band (alignment + FPCA)
#[extendr]
fn tolerance_elastic(data: RMatrix<f64>, argvals: Vec<f64>, ncomp: i32, nb: i32, coverage: f64, band_type: &str, max_iter: i32, seed: Robj) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let seed_val = if seed.is_null() { 42u64 } else { seed.as_integer().unwrap_or(42) as u64 };
    let bt = match band_type {
        "simultaneous" => tolerance::BandType::Simultaneous,
        _ => tolerance::BandType::Pointwise,
    };
    if n == 0 || m == 0 || argvals.len() != m {
        return list!(lower = Vec::<f64>::new(), upper = Vec::<f64>::new(), center = Vec::<f64>::new(), half_width = Vec::<f64>::new(), success = false).into();
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    match tolerance::elastic_tolerance_band(&fd, &argvals, ncomp as usize, nb as usize, coverage, bt, max_iter as usize, seed_val) {
        Ok(band) => list!(lower = band.lower, upper = band.upper, center = band.center, half_width = band.half_width, success = true).into(),
        Err(_) => list!(lower = Vec::<f64>::new(), upper = Vec::<f64>::new(), center = Vec::<f64>::new(), half_width = Vec::<f64>::new(), success = false).into(),
    }
}

// =============================================================================
// Phase & Elastic Config Tolerance Bands
// =============================================================================

/// Phase tolerance band (warping variation)
#[extendr]
fn tolerance_phase_rust(data: RMatrix<f64>, argvals: Vec<f64>, ncomp: i32, nb: i32, coverage: f64, band_type: &str, max_iter: i32, seed: Robj) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let seed_val = if seed.is_null() { 42u64 } else { seed.as_integer().unwrap_or(42) as u64 };
    let bt = match band_type {
        "simultaneous" => tolerance::BandType::Simultaneous,
        _ => tolerance::BandType::Pointwise,
    };
    if n == 0 || m == 0 || argvals.len() != m {
        return list!(success = false).into();
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    match tolerance::phase_tolerance_band(&fd, &argvals, ncomp as usize, nb as usize, coverage, bt, max_iter as usize, seed_val) {
        Ok(r) => list!(
            gamma_lower = r.gamma_lower,
            gamma_upper = r.gamma_upper,
            gamma_center = r.gamma_center,
            tangent_lower = r.tangent_band.lower,
            tangent_upper = r.tangent_band.upper,
            tangent_center = r.tangent_band.center,
            tangent_half_width = r.tangent_band.half_width,
            success = true
        ).into(),
        Err(_) => list!(success = false).into(),
    }
}

/// Elastic tolerance band with full config (amplitude + phase)
#[extendr]
fn tolerance_elastic_config_rust(data: RMatrix<f64>, argvals: Vec<f64>, ncomp_amplitude: i32, ncomp_phase: i32, nb: i32, coverage: f64, band_type: &str, max_iter: i32, tol: f64, seed: Robj) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let seed_val = if seed.is_null() { 42u64 } else { seed.as_integer().unwrap_or(42) as u64 };
    let bt = match band_type {
        "simultaneous" => tolerance::BandType::Simultaneous,
        _ => tolerance::BandType::Pointwise,
    };
    if n == 0 || m == 0 || argvals.len() != m {
        return list!(success = false).into();
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let mut config = tolerance::ElasticToleranceConfig::default();
    config.ncomp_amplitude = ncomp_amplitude as usize;
    config.ncomp_phase = ncomp_phase as usize;
    config.nb = nb as usize;
    config.coverage = coverage;
    config.band_type = bt;
    config.max_iter = max_iter as usize;
    config.tol = tol;
    config.seed = seed_val;

    match tolerance::elastic_tolerance_band_with_config(&fd, &argvals, &config) {
        Ok(r) => list!(
            amp_lower = r.amplitude.lower,
            amp_upper = r.amplitude.upper,
            amp_center = r.amplitude.center,
            amp_half_width = r.amplitude.half_width,
            phase_gamma_lower = r.phase.gamma_lower,
            phase_gamma_upper = r.phase.gamma_upper,
            phase_gamma_center = r.phase.gamma_center,
            success = true
        ).into(),
        Err(_) => list!(success = false).into(),
    }
}

// =============================================================================
// Streaming depth functions
// =============================================================================

/// Streaming depth: batch self-depth computation
#[extendr]
fn streaming_depth_batch(data: RMatrix<f64>, method: &str) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m == 0 {
        return Robj::from(Vec::<f64>::new());
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");

    let depths = match method {
        "MBD" => {
            let state = SortedReferenceState::from_reference(&fd);
            streaming_depth::StreamingMbd::new(state).depth_batch(&fd)
        }
        "BD" => {
            let state = FullReferenceState::from_reference(&fd);
            streaming_depth::StreamingBd::new(state).depth_batch(&fd)
        }
        _ => {
            // FM (default)
            let state = SortedReferenceState::from_reference(&fd);
            streaming_depth::StreamingFraimanMuniz::new(state, true).depth_batch(&fd)
        }
    };
    Robj::from(depths)
}

/// Streaming depth: new data against reference
#[extendr]
fn streaming_depth_vs_ref(ref_data: RMatrix<f64>, new_data: RMatrix<f64>, method: &str) -> Robj {
    let nr = ref_data.nrows();
    let mr = ref_data.ncols();
    let nn = new_data.nrows();
    let mn = new_data.ncols();
    if nr == 0 || nn == 0 || mr == 0 || mr != mn {
        return Robj::from(Vec::<f64>::new());
    }
    let ref_slice = ref_data.as_real_slice().unwrap();
    let new_slice = new_data.as_real_slice().unwrap();
    let ref_fd = FdMatrix::from_slice(ref_slice, nr, mr).expect("Invalid matrix dimensions");
    let new_fd = FdMatrix::from_slice(new_slice, nn, mn).expect("Invalid matrix dimensions");

    let depths = match method {
        "MBD" => {
            let state = SortedReferenceState::from_reference(&ref_fd);
            streaming_depth::StreamingMbd::new(state).depth_batch(&new_fd)
        }
        "BD" => {
            let state = FullReferenceState::from_reference(&ref_fd);
            streaming_depth::StreamingBd::new(state).depth_batch(&new_fd)
        }
        _ => {
            let state = SortedReferenceState::from_reference(&ref_fd);
            streaming_depth::StreamingFraimanMuniz::new(state, true).depth_batch(&new_fd)
        }
    };
    Robj::from(depths)
}

/// Streaming depth: single curve against reference
#[extendr]
fn streaming_depth_one(ref_data: RMatrix<f64>, curve: Vec<f64>, method: &str) -> Robj {
    let nr = ref_data.nrows();
    let mr = ref_data.ncols();
    if nr == 0 || mr == 0 || curve.len() != mr {
        return Robj::from(f64::NAN);
    }
    let ref_slice = ref_data.as_real_slice().unwrap();
    let ref_fd = FdMatrix::from_slice(ref_slice, nr, mr).expect("Invalid matrix dimensions");

    let depth = match method {
        "MBD" => {
            let state = SortedReferenceState::from_reference(&ref_fd);
            streaming_depth::StreamingMbd::new(state).depth_one(&curve)
        }
        "BD" => {
            let state = FullReferenceState::from_reference(&ref_fd);
            streaming_depth::StreamingBd::new(state).depth_one(&curve)
        }
        _ => {
            let state = SortedReferenceState::from_reference(&ref_fd);
            streaming_depth::StreamingFraimanMuniz::new(state, true).depth_one(&curve)
        }
    };
    Robj::from(depth)
}

// =============================================================================
// Scalar-on-function regression
// =============================================================================

/// Functional linear model (FPC-based)
#[extendr]
fn fregre_lm_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    scalar_covariates: Nullable<RMatrix<f64>>,
    ncomp: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid data dimensions");

    let scov = match scalar_covariates {
        Nullable::NotNull(sc) => {
            let sc_s = sc.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(sc_s, sc.nrows(), sc.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    let result = scalar_on_function::fregre_lm(&fd, &y, scov.as_ref(), ncomp as usize);
    match result {
        Ok(r) => {
            let beta_t = Robj::from(r.beta_t.clone());
            let gamma = Robj::from(r.gamma.clone());
            let fitted = Robj::from(r.fitted_values.clone());
            let resid = Robj::from(r.residuals.clone());
            let std_err = Robj::from(r.std_errors.clone());
            let coefs = Robj::from(r.coefficients.clone());
            let beta_se = Robj::from(r.beta_se.clone());
            r!(list!(
                intercept = r.intercept,
                beta_t = beta_t,
                beta_se = beta_se,
                gamma = gamma,
                fitted_values = fitted,
                residuals = resid,
                r_squared = r.r_squared,
                r_squared_adj = r.r_squared_adj,
                std_errors = std_err,
                ncomp = r.ncomp as i32,
                coefficients = coefs,
                residual_se = r.residual_se,
                gcv = r.gcv,
                aic = r.aic,
                bic = r.bic,
                // FPCA fields for prediction and explainability
                fpca_mean = Robj::from(r.fpca.mean.clone()),
                fpca_rotation_data = Robj::from(r.fpca.rotation.as_slice().to_vec()),
                fpca_rotation_nrow = r.fpca.rotation.nrows() as i32,
                fpca_rotation_ncol = r.fpca.rotation.ncols() as i32,
                fpca_scores_data = Robj::from(r.fpca.scores.as_slice().to_vec()),
                fpca_scores_nrow = r.fpca.scores.nrows() as i32,
                fpca_scores_ncol = r.fpca.scores.ncols() as i32
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Nonparametric functional regression with mixed predictors
#[extendr]
fn fregre_np_mixed_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    argvals: Vec<f64>,
    scalar_covariates: Nullable<RMatrix<f64>>,
    h_func: f64,
    h_scalar: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid data dimensions");

    let scov = match scalar_covariates {
        Nullable::NotNull(sc) => {
            let sc_s = sc.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(sc_s, sc.nrows(), sc.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    let result = scalar_on_function::fregre_np_mixed(&fd, &y, &argvals, scov.as_ref(), h_func, h_scalar);
    match result {
        Ok(r) => {
            r!(list!(
                fitted_values = Robj::from(r.fitted_values),
                residuals = Robj::from(r.residuals),
                r_squared = r.r_squared,
                h_func = r.h_func,
                h_scalar = r.h_scalar,
                cv_error = r.cv_error
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Functional logistic regression
#[extendr]
fn functional_logistic_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    scalar_covariates: Nullable<RMatrix<f64>>,
    ncomp: i32,
    max_iter: i32,
    tol: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid data dimensions");

    let scov = match scalar_covariates {
        Nullable::NotNull(sc) => {
            let sc_s = sc.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(sc_s, sc.nrows(), sc.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    let result = scalar_on_function::functional_logistic(
        &fd, &y, scov.as_ref(), ncomp as usize, max_iter as usize, tol,
    );
    match result {
        Ok(r) => {
            let pred_classes: Vec<i32> = r.predicted_classes.iter().map(|&c| c as i32).collect();
            r!(list!(
                intercept = r.intercept,
                beta_t = Robj::from(r.beta_t),
                gamma = Robj::from(r.gamma),
                probabilities = Robj::from(r.probabilities),
                predicted_classes = Robj::from(pred_classes),
                ncomp = r.ncomp as i32,
                accuracy = r.accuracy,
                coefficients = Robj::from(r.coefficients),
                log_likelihood = r.log_likelihood,
                aic = r.aic,
                bic = r.bic,
                iterations = r.iterations as i32,
                fpca_mean = Robj::from(r.fpca.mean.clone()),
                fpca_rotation_data = Robj::from(r.fpca.rotation.as_slice().to_vec()),
                fpca_rotation_nrow = r.fpca.rotation.nrows() as i32,
                fpca_rotation_ncol = r.fpca.rotation.ncols() as i32,
                fpca_scores_data = Robj::from(r.fpca.scores.as_slice().to_vec()),
                fpca_scores_nrow = r.fpca.scores.nrows() as i32,
                fpca_scores_ncol = r.fpca.scores.ncols() as i32,
                beta_se = Robj::from(r.beta_se),
                std_errors = Robj::from(r.std_errors)
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Cross-validation for FPC component selection
#[extendr]
fn fregre_cv_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    scalar_covariates: Nullable<RMatrix<f64>>,
    k_min: i32,
    k_max: i32,
    n_folds: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid data dimensions");

    let scov = match scalar_covariates {
        Nullable::NotNull(sc) => {
            let sc_s = sc.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(sc_s, sc.nrows(), sc.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    let result = scalar_on_function::fregre_cv(
        &fd, &y, scov.as_ref(), k_min as usize, k_max as usize, n_folds as usize,
    );
    match result {
        Ok(r) => {
            let k_values: Vec<i32> = r.k_values.iter().map(|&k| k as i32).collect();
            r!(list!(
                k_values = Robj::from(k_values),
                cv_errors = Robj::from(r.cv_errors),
                optimal_k = r.optimal_k as i32,
                min_cv_error = r.min_cv_error
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Bootstrap confidence intervals for β(t) in functional linear model
#[extendr]
fn bootstrap_ci_fregre_lm_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    scalar_covariates: Nullable<RMatrix<f64>>,
    ncomp: i32,
    n_boot: i32,
    alpha: f64,
    seed: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid data dimensions");

    let scov = match scalar_covariates {
        Nullable::NotNull(sc) => {
            let sc_s = sc.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(sc_s, sc.nrows(), sc.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    let result = scalar_on_function::bootstrap_ci_fregre_lm(
        &fd, &y, scov.as_ref(), ncomp as usize, n_boot as usize, alpha, seed as u64,
    );
    match result {
        Ok(r) => {
            r!(list!(
                lower = Robj::from(r.lower),
                upper = Robj::from(r.upper),
                center = Robj::from(r.center),
                sim_lower = Robj::from(r.sim_lower),
                sim_upper = Robj::from(r.sim_upper),
                n_boot_success = r.n_boot_success as i32
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Bootstrap confidence intervals for β(t) in functional logistic model
#[extendr]
fn bootstrap_ci_functional_logistic_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    scalar_covariates: Nullable<RMatrix<f64>>,
    ncomp: i32,
    n_boot: i32,
    alpha: f64,
    seed: i32,
    max_iter: i32,
    tol: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid data dimensions");

    let scov = match scalar_covariates {
        Nullable::NotNull(sc) => {
            let sc_s = sc.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(sc_s, sc.nrows(), sc.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    let result = scalar_on_function::bootstrap_ci_functional_logistic(
        &fd, &y, scov.as_ref(), ncomp as usize, n_boot as usize, alpha, seed as u64,
        max_iter as usize, tol,
    );
    match result {
        Ok(r) => {
            r!(list!(
                lower = Robj::from(r.lower),
                upper = Robj::from(r.upper),
                center = Robj::from(r.center),
                sim_lower = Robj::from(r.sim_lower),
                sim_upper = Robj::from(r.sim_upper),
                n_boot_success = r.n_boot_success as i32
            ))
        }
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Function-on-scalar regression
// =============================================================================

/// Helper to return FdMatrix as R list with data, nrow, ncol
fn fdmatrix_to_robj(mat: &FdMatrix) -> Robj {
    let s = mat.as_slice();
    RMatrix::new_matrix(mat.nrows(), mat.ncols(), |i, j| s[i + j * mat.nrows()]).into()
}

/// Function-on-scalar regression (penalized)
#[extendr]
fn fosr_rust(data: RMatrix<f64>, predictors: RMatrix<f64>, lambda: f64) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");

    let ps = predictors.as_real_slice().unwrap();
    let pred = FdMatrix::from_slice(ps, predictors.nrows(), predictors.ncols()).expect("Invalid predictors");

    let result = function_on_scalar::fosr(&fd, &pred, lambda);
    match result {
        Ok(r) => {
            r!(list!(
                intercept = Robj::from(r.intercept),
                beta = fdmatrix_to_robj(&r.beta),
                fitted = fdmatrix_to_robj(&r.fitted),
                residuals = fdmatrix_to_robj(&r.residuals),
                r_squared_t = Robj::from(r.r_squared_t),
                r_squared = r.r_squared,
                beta_se = fdmatrix_to_robj(&r.beta_se),
                lambda = r.lambda,
                gcv = r.gcv
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// FPC-based function-on-scalar regression
#[extendr]
fn fosr_fpc_rust(data: RMatrix<f64>, predictors: RMatrix<f64>, ncomp: i32) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");

    let ps = predictors.as_real_slice().unwrap();
    let pred = FdMatrix::from_slice(ps, predictors.nrows(), predictors.ncols()).expect("Invalid predictors");

    let result = function_on_scalar::fosr_fpc(&fd, &pred, ncomp as usize);
    match result {
        Ok(r) => {
            r!(list!(
                intercept = Robj::from(r.intercept),
                beta = fdmatrix_to_robj(&r.beta),
                fitted = fdmatrix_to_robj(&r.fitted),
                residuals = fdmatrix_to_robj(&r.residuals),
                r_squared_t = Robj::from(r.r_squared_t),
                r_squared = r.r_squared,
                ncomp = r.ncomp as i32
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Functional ANOVA
#[extendr]
fn fanova_rust(data: RMatrix<f64>, groups: Vec<i32>, n_perm: i32) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");

    let groups_usize: Vec<usize> = groups.iter().map(|&g| g as usize).collect();

    let result = function_on_scalar::fanova(&fd, &groups_usize, n_perm as usize);
    match result {
        Ok(r) => {
            let group_labels: Vec<i32> = r.group_labels.iter().map(|&g| g as i32).collect();
            r!(list!(
                group_means = fdmatrix_to_robj(&r.group_means),
                overall_mean = Robj::from(r.overall_mean),
                f_statistic_t = Robj::from(r.f_statistic_t),
                global_statistic = r.global_statistic,
                p_value = r.p_value,
                n_perm = r.n_perm as i32,
                n_groups = r.n_groups as i32,
                group_labels = Robj::from(group_labels)
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Predict from function-on-scalar regression
#[extendr]
fn predict_fosr_rust(
    intercept: Vec<f64>,
    beta_data: RMatrix<f64>,
    new_predictors: RMatrix<f64>,
    _lambda: f64,
) -> Robj {
    let bs = beta_data.as_real_slice().unwrap();
    let beta = FdMatrix::from_slice(bs, beta_data.nrows(), beta_data.ncols()).expect("Invalid beta");

    let m = intercept.len();
    let ps = new_predictors.as_real_slice().unwrap();
    let pred = FdMatrix::from_slice(ps, new_predictors.nrows(), new_predictors.ncols()).expect("Invalid predictors");

    // Inline predict_fosr logic since FosrResult is non-exhaustive
    let n_new = pred.nrows();
    let p = beta.nrows();
    let mut predicted = FdMatrix::zeros(n_new, m);
    for i in 0..n_new {
        for t in 0..m {
            let mut yhat = intercept[t];
            for j in 0..p {
                yhat += pred[(i, j)] * beta[(j, t)];
            }
            predicted[(i, t)] = yhat;
        }
    }
    fdmatrix_to_robj(&predicted)
}

// =============================================================================
// Classification
// =============================================================================

/// Helper: convert ClassifResult to R list
fn classif_result_to_robj(r: &classification::ClassifResult) -> Robj {
    let predicted: Vec<i32> = r.predicted.iter().map(|&c| c as i32).collect();
    let confusion: Vec<Vec<i32>> = r.confusion.iter()
        .map(|row| row.iter().map(|&c| c as i32).collect())
        .collect();
    let conf_flat: Vec<i32> = confusion.iter().flat_map(|row| row.iter().copied()).collect();
    let nc = r.n_classes;

    let probs = match &r.probabilities {
        Some(p) => fdmatrix_to_robj(p),
        None => r!(NULL),
    };

    let conf_mat = RMatrix::new_matrix(nc, nc, |i, j| conf_flat[i * nc + j] as f64);

    r!(list!(
        predicted = Robj::from(predicted),
        probabilities = probs,
        accuracy = r.accuracy,
        confusion = Robj::from(conf_mat),
        n_classes = r.n_classes as i32,
        ncomp = r.ncomp as i32
    ))
}

/// LDA classification
#[extendr]
fn fclassif_lda_rust(
    data: RMatrix<f64>,
    y: Vec<i32>,
    covariates: Nullable<RMatrix<f64>>,
    ncomp: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");
    let y_usize: Vec<usize> = y.iter().map(|&c| c as usize).collect();

    let cov = match covariates {
        Nullable::NotNull(c) => {
            let cs = c.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(cs, c.nrows(), c.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    match classification::fclassif_lda(&fd, &y_usize, cov.as_ref(), ncomp as usize) {
        Ok(r) => classif_result_to_robj(&r),
        Err(_) => r!(NULL),
    }
}

/// QDA classification
#[extendr]
fn fclassif_qda_rust(
    data: RMatrix<f64>,
    y: Vec<i32>,
    covariates: Nullable<RMatrix<f64>>,
    ncomp: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");
    let y_usize: Vec<usize> = y.iter().map(|&c| c as usize).collect();

    let cov = match covariates {
        Nullable::NotNull(c) => {
            let cs = c.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(cs, c.nrows(), c.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    match classification::fclassif_qda(&fd, &y_usize, cov.as_ref(), ncomp as usize) {
        Ok(r) => classif_result_to_robj(&r),
        Err(_) => r!(NULL),
    }
}

/// kNN classification
#[extendr]
fn fclassif_knn_rust(
    data: RMatrix<f64>,
    y: Vec<i32>,
    covariates: Nullable<RMatrix<f64>>,
    ncomp: i32,
    k_nn: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");
    let y_usize: Vec<usize> = y.iter().map(|&c| c as usize).collect();

    let cov = match covariates {
        Nullable::NotNull(c) => {
            let cs = c.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(cs, c.nrows(), c.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    match classification::fclassif_knn(&fd, &y_usize, cov.as_ref(), ncomp as usize, k_nn as usize) {
        Ok(r) => classif_result_to_robj(&r),
        Err(_) => r!(NULL),
    }
}

/// Kernel classification
#[extendr]
fn fclassif_kernel_rust(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    y: Vec<i32>,
    covariates: Nullable<RMatrix<f64>>,
    h_func: f64,
    h_scalar: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");
    let y_usize: Vec<usize> = y.iter().map(|&c| c as usize).collect();

    let cov = match covariates {
        Nullable::NotNull(c) => {
            let cs = c.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(cs, c.nrows(), c.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    match classification::fclassif_kernel(&fd, &y_usize, &argvals, cov.as_ref(), h_func, h_scalar) {
        Ok(r) => classif_result_to_robj(&r),
        Err(_) => r!(NULL),
    }
}

/// DD-plot classification
#[extendr]
fn fclassif_dd_rust(
    data: RMatrix<f64>,
    y: Vec<i32>,
    covariates: Nullable<RMatrix<f64>>,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");
    let y_usize: Vec<usize> = y.iter().map(|&c| c as usize).collect();

    let cov = match covariates {
        Nullable::NotNull(c) => {
            let cs = c.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(cs, c.nrows(), c.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    match classification::fclassif_dd(&fd, &y_usize, cov.as_ref()) {
        Ok(r) => classif_result_to_robj(&r),
        Err(_) => r!(NULL),
    }
}

/// Cross-validated classification
#[extendr]
fn fclassif_cv_rust(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    y: Vec<i32>,
    covariates: Nullable<RMatrix<f64>>,
    method: &str,
    ncomp: i32,
    nfold: i32,
    seed: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");
    let y_usize: Vec<usize> = y.iter().map(|&c| c as usize).collect();

    let cov = match covariates {
        Nullable::NotNull(c) => {
            let cs = c.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(cs, c.nrows(), c.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    match classification::fclassif_cv(
        &fd, &argvals, &y_usize, cov.as_ref(), method,
        ncomp as usize, nfold as usize, seed as u64,
    ) {
        Ok(r) => {
            r!(list!(
                error_rate = r.error_rate,
                fold_errors = Robj::from(r.fold_errors),
                best_ncomp = r.best_ncomp as i32
            ))
        }
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// GMM clustering
// =============================================================================

/// GMM clustering with automatic K selection
#[extendr]
fn gmm_cluster_rust(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    covariates: Nullable<RMatrix<f64>>,
    k_range: Vec<i32>,
    nbasis: i32,
    basis_type: i32,
    cov_type: &str,
    cov_weight: f64,
    max_iter: i32,
    tol: f64,
    n_init: i32,
    seed: i32,
    use_icl: bool,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");

    let cov = match covariates {
        Nullable::NotNull(c) => {
            let cs = c.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(cs, c.nrows(), c.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    let k_range_usize: Vec<usize> = k_range.iter().map(|&k| k as usize).collect();
    let ct = match cov_type {
        "diagonal" | "diag" => gmm::CovType::Diagonal,
        _ => gmm::CovType::Full,
    };

    let bt = fdars_core::basis::ProjectionBasisType::from_i32(basis_type);
    let result = gmm::gmm_cluster(
        &fd, &argvals, cov.as_ref(), &k_range_usize,
        nbasis as usize, bt, ct, cov_weight,
        max_iter as usize, tol, n_init as usize, seed as u64, use_icl,
    );

    match result {
        Ok(r) => {
            let best = &r.best;
            let cluster: Vec<i32> = best.cluster.iter().map(|&c| c as i32 + 1).collect();
            let membership = fdmatrix_to_robj(&best.membership);
            let weights = Robj::from(best.weights.clone());
            let bic_k: Vec<i32> = r.bic_values.iter().map(|&(k, _)| k as i32).collect();
            let bic_v: Vec<f64> = r.bic_values.iter().map(|&(_, v)| v).collect();
            let icl_v: Vec<f64> = r.icl_values.iter().map(|&(_, v)| v).collect();

            // Flatten means for R
            let means_flat: Vec<f64> = best.means.iter().flat_map(|m| m.iter().copied()).collect();
            let k = best.k;
            let d = best.d;
            let means_mat = RMatrix::new_matrix(k, d, |i, j| means_flat[i * d + j]);

            // Flatten covariances
            let cov_flat: Vec<f64> = best.covariances.iter().flat_map(|c| c.iter().copied()).collect();

            r!(list!(
                cluster = Robj::from(cluster),
                membership = membership,
                means = Robj::from(means_mat),
                weights = weights,
                covariances = Robj::from(cov_flat),
                log_likelihood = best.log_likelihood,
                bic = best.bic,
                icl = best.icl,
                iterations = best.iterations as i32,
                converged = best.converged,
                k = best.k as i32,
                d = best.d as i32,
                bic_k = Robj::from(bic_k),
                bic_values = Robj::from(bic_v),
                icl_values = Robj::from(icl_v),
                cov_type = cov_type,
                nbasis = nbasis,
                basis_type = basis_type,
                cov_weight = cov_weight
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Raw GMM EM on feature matrix
#[extendr]
fn gmm_em_rust(
    features: RMatrix<f64>,
    k: i32,
    cov_type: &str,
    max_iter: i32,
    tol: f64,
    seed: i32,
) -> Robj {
    let n = features.nrows();
    let d = features.ncols();
    let fs = features.as_real_slice().unwrap();

    // Convert R matrix (col-major) to Vec<Vec<f64>> (row-major)
    let feat_vecs: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..d).map(|j| fs[i + j * n]).collect())
        .collect();

    let ct = match cov_type {
        "diagonal" | "diag" => gmm::CovType::Diagonal,
        _ => gmm::CovType::Full,
    };

    let result = gmm::gmm_em(&feat_vecs, k as usize, ct, max_iter as usize, tol, seed as u64);
    match result {
        Ok(r) => {
            let cluster: Vec<i32> = r.cluster.iter().map(|&c| c as i32 + 1).collect();
            r!(list!(
                cluster = Robj::from(cluster),
                membership = fdmatrix_to_robj(&r.membership),
                weights = Robj::from(r.weights),
                log_likelihood = r.log_likelihood,
                bic = r.bic,
                icl = r.icl,
                converged = r.converged,
                k = r.k as i32,
                d = r.d as i32
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Predict from GMM
#[extendr]
fn predict_gmm_rust(
    new_data: RMatrix<f64>,
    argvals: Vec<f64>,
    new_covariates: Nullable<RMatrix<f64>>,
    means_data: RMatrix<f64>,
    covariances: Vec<f64>,
    weights: Vec<f64>,
    k: i32,
    d: i32,
    nbasis: i32,
    basis_type: i32,
    cov_weight: f64,
    cov_type: &str,
) -> Robj {
    let n = new_data.nrows();
    let m = new_data.ncols();
    let ds = new_data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");

    let ncov = match new_covariates {
        Nullable::NotNull(c) => {
            let cs = c.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(cs, c.nrows(), c.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    let ku = k as usize;
    let du = d as usize;
    let ct = match cov_type {
        "diagonal" | "diag" => gmm::CovType::Diagonal,
        _ => gmm::CovType::Full,
    };

    // Reconstruct means from R matrix
    let ms = means_data.as_real_slice().unwrap();
    let means_vec: Vec<Vec<f64>> = (0..ku)
        .map(|i| (0..du).map(|j| ms[i + j * ku]).collect())
        .collect();

    // Reconstruct covariances
    let cov_per = if ct == gmm::CovType::Full { du * du } else { du };
    let covs: Vec<Vec<f64>> = (0..ku)
        .map(|i| covariances[i * cov_per..(i + 1) * cov_per].to_vec())
        .collect();

    // Build a GmmResult for prediction
    let gmm_result = gmm::GmmResult::from_parts(
        vec![0; 1],
        FdMatrix::zeros(1, ku),
        means_vec,
        covs,
        weights,
        0.0,
        0.0,
        0.0,
        0,
        true,
        ku,
        du,
    );

    let bt = fdars_core::basis::ProjectionBasisType::from_i32(basis_type);
    let result = gmm::predict_gmm(&fd, &argvals, ncov.as_ref(), &gmm_result, nbasis as usize, bt, cov_weight, ct);
    match result {
        Ok((cluster, membership)) => {
            let cluster_r: Vec<i32> = cluster.iter().map(|&c| c as i32 + 1).collect();
            r!(list!(
                cluster = Robj::from(cluster_r),
                membership = fdmatrix_to_robj(&membership)
            ))
        }
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Functional mixed models (FAMM)
// =============================================================================

/// Functional mixed model
#[extendr]
fn fmm_rust(
    data: RMatrix<f64>,
    subject_ids: Vec<i32>,
    covariates: Nullable<RMatrix<f64>>,
    ncomp: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");
    let sids: Vec<usize> = subject_ids.iter().map(|&s| s as usize).collect();

    let cov = match covariates {
        Nullable::NotNull(c) => {
            let cs = c.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(cs, c.nrows(), c.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    let result = famm::fmm(&fd, &sids, cov.as_ref(), ncomp as usize);
    match result {
        Ok(r) => {
            r!(list!(
                mean_function = Robj::from(r.mean_function),
                beta_functions = fdmatrix_to_robj(&r.beta_functions),
                random_effects = fdmatrix_to_robj(&r.random_effects),
                fitted = fdmatrix_to_robj(&r.fitted),
                residuals = fdmatrix_to_robj(&r.residuals),
                random_variance = Robj::from(r.random_variance),
                sigma2_eps = r.sigma2_eps,
                sigma2_u = Robj::from(r.sigma2_u),
                ncomp = r.ncomp as i32,
                n_subjects = r.n_subjects as i32,
                eigenvalues = Robj::from(r.eigenvalues)
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Predict from functional mixed model
#[extendr]
fn fmm_predict_rust(
    mean_function: Vec<f64>,
    beta_data: RMatrix<f64>,
    new_covariates: Nullable<RMatrix<f64>>,
) -> Robj {
    let m = mean_function.len();
    let bs = beta_data.as_real_slice().unwrap();
    let beta = FdMatrix::from_slice(bs, beta_data.nrows(), beta_data.ncols()).expect("Invalid beta");

    let ncov = match new_covariates {
        Nullable::NotNull(c) => {
            let cs = c.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(cs, c.nrows(), c.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    // Construct a minimal FmmResult for prediction
    let fmm_result = famm::FmmResult {
        mean_function: mean_function.clone(),
        beta_functions: beta.clone(),
        random_effects: FdMatrix::zeros(1, m),
        fitted: FdMatrix::zeros(1, m),
        residuals: FdMatrix::zeros(1, m),
        random_variance: vec![0.0; m],
        sigma2_eps: 0.0,
        sigma2_u: vec![0.0],
        ncomp: 0,
        n_subjects: 0,
        eigenvalues: vec![0.0],
    };

    let result = famm::fmm_predict(&fmm_result, ncov.as_ref());
    fdmatrix_to_robj(&result)
}

/// Permutation test for fixed effects in FMM
#[extendr]
fn fmm_test_fixed_rust(
    data: RMatrix<f64>,
    subject_ids: Vec<i32>,
    covariates: RMatrix<f64>,
    ncomp: i32,
    n_perm: i32,
    seed: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");
    let sids: Vec<usize> = subject_ids.iter().map(|&s| s as usize).collect();

    let cs = covariates.as_real_slice().unwrap();
    let cov = FdMatrix::from_slice(cs, covariates.nrows(), covariates.ncols()).expect("Invalid covariates");

    let result = famm::fmm_test_fixed(
        &fd, &sids, &cov, ncomp as usize, n_perm as usize, seed as u64,
    );
    match result {
        Ok(r) => {
            r!(list!(
                f_statistics = Robj::from(r.f_statistics),
                p_values = Robj::from(r.p_values)
            ))
        }
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Seeded depth variants
// =============================================================================

/// Random projection depth with optional seed
#[extendr]
fn depth_rp_1d_seeded(
    data_obj: RMatrix<f64>,
    data_ori: RMatrix<f64>,
    nproj: i32,
    seed: Nullable<i32>,
) -> Robj {
    let n_obj = data_obj.nrows();
    let m = data_obj.ncols();
    let n_ori = data_ori.nrows();

    let ds_obj = data_obj.as_real_slice().unwrap();
    let ds_ori = data_ori.as_real_slice().unwrap();

    let fd_obj = FdMatrix::from_slice(ds_obj, n_obj, m).expect("Invalid data_obj");
    let fd_ori = FdMatrix::from_slice(ds_ori, n_ori, m).expect("Invalid data_ori");

    let seed_val = match seed {
        Nullable::NotNull(s) => Some(s as u64),
        Nullable::Null => None,
    };

    let depths = core_depth::random_projection_1d_seeded(&fd_obj, &fd_ori, nproj as usize, seed_val);
    Robj::from(depths)
}

/// Random Tukey depth with optional seed
#[extendr]
fn depth_rt_1d_seeded(
    data_obj: RMatrix<f64>,
    data_ori: RMatrix<f64>,
    nproj: i32,
    seed: Nullable<i32>,
) -> Robj {
    let n_obj = data_obj.nrows();
    let m = data_obj.ncols();
    let n_ori = data_ori.nrows();

    let ds_obj = data_obj.as_real_slice().unwrap();
    let ds_ori = data_ori.as_real_slice().unwrap();

    let fd_obj = FdMatrix::from_slice(ds_obj, n_obj, m).expect("Invalid data_obj");
    let fd_ori = FdMatrix::from_slice(ds_ori, n_ori, m).expect("Invalid data_ori");

    let seed_val = match seed {
        Nullable::NotNull(s) => Some(s as u64),
        Nullable::Null => None,
    };

    let depths = core_depth::random_tukey_1d_seeded(&fd_obj, &fd_ori, nproj as usize, seed_val);
    Robj::from(depths)
}

// =============================================================================
// RPD Depth
// =============================================================================

/// Random Projection Depth with Derivatives (seeded)
#[extendr]
fn depth_rpd_1d_rust(
    data_obj: RMatrix<f64>,
    data_ori: RMatrix<f64>,
    argvals: Vec<f64>,
    nproj: i32,
    nderiv: i32,
    seed: Nullable<i32>,
) -> Robj {
    let n_obj = data_obj.nrows();
    let m = data_obj.ncols();
    let n_ori = data_ori.nrows();

    let ds_obj = data_obj.as_real_slice().unwrap();
    let fd_obj = FdMatrix::from_slice(ds_obj, n_obj, m).expect("Invalid data_obj");
    let ds_ori = data_ori.as_real_slice().unwrap();
    let fd_ori = FdMatrix::from_slice(ds_ori, n_ori, m).expect("Invalid data_ori");

    let seed_val = match seed {
        Nullable::NotNull(s) => Some(s as u64),
        Nullable::Null => None,
    };

    let depths = core_depth::rpd_depth_1d_seeded(
        &fd_obj, &fd_ori, &argvals, nproj as usize, nderiv as usize, seed_val,
    );
    Robj::from(depths)
}

// =============================================================================
// Helper: Reconstruct FregreLmResult from flattened R parameters
// =============================================================================

/// Reconstruct a FregreLmResult from the 20 flattened parameters that .explain_call() sends.
fn reconstruct_fregre_lm(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
) -> fdars_core::scalar_on_function::FregreLmResult {
    let rotation = FdMatrix::from_slice(
        &fpca_rotation_data,
        fpca_rotation_nrow as usize,
        fpca_rotation_ncol as usize,
    )
    .unwrap_or_else(|_| FdMatrix::zeros(0, 0));

    let scores = FdMatrix::from_slice(
        &fpca_scores_data,
        fpca_scores_nrow as usize,
        fpca_scores_ncol as usize,
    )
    .unwrap_or_else(|_| FdMatrix::zeros(0, 0));

    let n = scores.nrows();
    let m = rotation.nrows();

    fdars_core::scalar_on_function::FregreLmResult {
        intercept,
        beta_t,
        beta_se,
        gamma,
        fitted_values,
        residuals,
        r_squared,
        r_squared_adj,
        std_errors,
        ncomp: ncomp as usize,
        coefficients,
        residual_se,
        gcv,
        aic: 0.0,
        bic: 0.0,
        fpca: fdars_core::regression::FpcaResult {
            singular_values: vec![],
            rotation,
            scores,
            mean: fpca_mean,
            centered: FdMatrix::zeros(n, m),
        },
    }
}

/// Reconstruct a FunctionalLogisticResult from the 19 flattened parameters.
fn reconstruct_functional_logistic(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    probabilities: Vec<f64>,
    predicted_classes: Vec<f64>,
    ncomp: i32,
    accuracy: f64,
    std_errors: Vec<f64>,
    coefficients: Vec<f64>,
    log_likelihood: f64,
    iterations: i32,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
) -> fdars_core::scalar_on_function::FunctionalLogisticResult {
    let rotation = FdMatrix::from_slice(
        &fpca_rotation_data,
        fpca_rotation_nrow as usize,
        fpca_rotation_ncol as usize,
    )
    .unwrap_or_else(|_| FdMatrix::zeros(0, 0));

    let scores = FdMatrix::from_slice(
        &fpca_scores_data,
        fpca_scores_nrow as usize,
        fpca_scores_ncol as usize,
    )
    .unwrap_or_else(|_| FdMatrix::zeros(0, 0));

    let n = scores.nrows();
    let m = rotation.nrows();

    fdars_core::scalar_on_function::FunctionalLogisticResult {
        intercept,
        beta_t,
        beta_se,
        gamma,
        probabilities,
        predicted_classes: predicted_classes.iter().map(|&x| x as usize).collect(),
        ncomp: ncomp as usize,
        accuracy,
        std_errors,
        coefficients,
        log_likelihood,
        iterations: iterations as usize,
        aic: 0.0,
        bic: 0.0,
        fpca: fdars_core::regression::FpcaResult {
            singular_values: vec![],
            rotation,
            scores,
            mean: fpca_mean,
            centered: FdMatrix::zeros(n, m),
        },
    }
}

/// Reconstruct a KarcherMeanResult from flattened R parameters.
fn reconstruct_karcher_mean(
    aligned_data: RMatrix<f64>,
    gammas_data: RMatrix<f64>,
    mean_curve: Vec<f64>,
    mean_srsf: Vec<f64>,
    _aligned_srsfs: Nullable<RMatrix<f64>>,
    _argvals: Vec<f64>,
) -> fdars_core::alignment::KarcherMeanResult {
    let n = aligned_data.nrows();
    let m = aligned_data.ncols();
    let ad_s = aligned_data.as_real_slice().unwrap();
    let fd_aligned = FdMatrix::from_slice(ad_s, n, m).expect("Invalid aligned_data");

    let gd_s = gammas_data.as_real_slice().unwrap();
    let fd_gammas = FdMatrix::from_slice(gd_s, gammas_data.nrows(), gammas_data.ncols())
        .expect("Invalid gammas");

    let aligned_srsfs = match _aligned_srsfs {
        Nullable::NotNull(s) => {
            let ss = s.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid aligned_srsfs"))
        }
        Nullable::Null => None,
    };

    fdars_core::alignment::KarcherMeanResult::new(
        mean_curve, mean_srsf, fd_gammas, fd_aligned,
        0, true, aligned_srsfs,
    )
}

// =============================================================================
// Explain wrappers (27 functions)
// =============================================================================

/// Partial dependence plot for a single FPC component
#[extendr]
fn explain_pdp_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
    component: i32,
    n_grid: i32,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::functional_pdp(&fit, &fd, None, component as usize, n_grid as usize) {
        Ok(r) => r!(list!(
            grid_values = Robj::from(r.grid_values),
            pdp_curve = Robj::from(r.pdp_curve),
            ice_curves = fdmatrix_to_robj(&r.ice_curves),
            component = r.component as i32
        )),
        Err(_) => r!(NULL),
    }
}

/// Beta decomposition by FPC components
#[extendr]
fn explain_beta_decomposition_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    match explain::beta_decomposition(&fit) {
        Ok(r) => {
            // Convert Vec<Vec<f64>> to an R matrix (ncomp x m)
            let nc = r.components.len();
            let m = if nc > 0 { r.components[0].len() } else { 0 };
            let mut flat = vec![0.0; nc * m];
            for i in 0..nc {
                for j in 0..m {
                    flat[j * nc + i] = r.components[i][j]; // column-major
                }
            }
            let comp_mat = RMatrix::new_matrix(nc, m, |i, j| flat[j * nc + i]);
            r!(list!(
                components = Robj::from(comp_mat),
                coefficients = Robj::from(r.coefficients),
                variance_proportion = Robj::from(r.variance_proportion)
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Pointwise importance of beta(t)
#[extendr]
fn explain_pointwise_importance_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    match explain::pointwise_importance(&fit) {
        Ok(r) => r!(list!(
            importance = Robj::from(r.importance),
            importance_normalized = Robj::from(r.importance_normalized),
            component_importance = fdmatrix_to_robj(&r.component_importance),
            score_variance = Robj::from(r.score_variance)
        )),
        Err(_) => r!(NULL),
    }
}

/// Significant regions from beta(t) and its standard errors
#[extendr]
fn explain_significant_regions_rust(
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    z_alpha: f64,
) -> Robj {
    match explain::significant_regions_from_se(&beta_t, &beta_se, z_alpha) {
        Ok(regions) => {
            let starts: Vec<i32> = regions.iter().map(|r| r.start_idx as i32).collect();
            let ends: Vec<i32> = regions.iter().map(|r| r.end_idx as i32).collect();
            let dirs: Vec<String> = regions
                .iter()
                .map(|r| match r.direction {
                    fdars_core::explain::SignificanceDirection::Positive => "positive".to_string(),
                    fdars_core::explain::SignificanceDirection::Negative => "negative".to_string(),
                    _ => "unknown".to_string(),
                })
                .collect();
            r!(list!(
                start_idx = Robj::from(starts),
                end_idx = Robj::from(ends),
                direction = Robj::from(dirs)
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// SHAP values for FPC scores
#[extendr]
fn explain_shap_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::fpc_shap_values(&fit, &fd, None) {
        Ok(r) => r!(list!(
            values = fdmatrix_to_robj(&r.values),
            base_value = r.base_value,
            mean_scores = Robj::from(r.mean_scores)
        )),
        Err(_) => r!(NULL),
    }
}

/// Influence diagnostics (leverage, Cook's distance)
#[extendr]
fn explain_influence_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::influence_diagnostics(&fit, &fd, None) {
        Ok(r) => r!(list!(
            leverage = Robj::from(r.leverage),
            cooks_distance = Robj::from(r.cooks_distance),
            p = r.p as i32,
            mse = r.mse
        )),
        Err(_) => r!(NULL),
    }
}

/// Variance inflation factors for FPC scores
#[extendr]
fn explain_vif_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::fpc_vif(&fit, &fd, None) {
        Ok(r) => r!(list!(
            vif = Robj::from(r.vif),
            labels = Robj::from(r.labels),
            mean_vif = r.mean_vif,
            n_moderate = r.n_moderate as i32,
            n_severe = r.n_severe as i32
        )),
        Err(_) => r!(NULL),
    }
}

/// DFBETAS and DFFITS diagnostics
#[extendr]
fn explain_dfbetas_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::dfbetas_dffits(&fit, &fd, None) {
        Ok(r) => r!(list!(
            dfbetas = fdmatrix_to_robj(&r.dfbetas),
            dffits = Robj::from(r.dffits),
            studentized_residuals = Robj::from(r.studentized_residuals),
            p = r.p as i32,
            dfbetas_cutoff = r.dfbetas_cutoff,
            dffits_cutoff = r.dffits_cutoff
        )),
        Err(_) => r!(NULL),
    }
}

/// Functional saliency map
#[extendr]
fn explain_saliency_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::functional_saliency(&fit, &fd, None) {
        Ok(r) => r!(list!(
            saliency_map = fdmatrix_to_robj(&r.saliency_map),
            mean_absolute_saliency = Robj::from(r.mean_absolute_saliency)
        )),
        Err(_) => r!(NULL),
    }
}

/// FPC permutation importance
#[extendr]
fn explain_importance_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
    y: Vec<f64>,
    n_perm: i32,
    seed: i32,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::fpc_permutation_importance(&fit, &fd, &y, n_perm as usize, seed as u64) {
        Ok(r) => r!(list!(
            importance = Robj::from(r.importance),
            baseline_metric = r.baseline_metric,
            permuted_metric = Robj::from(r.permuted_metric)
        )),
        Err(_) => r!(NULL),
    }
}

/// Leave-one-out cross-validation and PRESS
#[extendr]
fn explain_loo_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
    y: Vec<f64>,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::loo_cv_press(&fit, &fd, &y, None) {
        Ok(r) => r!(list!(
            loo_residuals = Robj::from(r.loo_residuals),
            press = r.press,
            loo_r_squared = r.loo_r_squared,
            leverage = Robj::from(r.leverage),
            tss = r.tss
        )),
        Err(_) => r!(NULL),
    }
}

/// Sobol sensitivity indices
#[extendr]
fn explain_sobol_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
    y: Vec<f64>,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::sobol_indices(&fit, &fd, &y, None) {
        Ok(r) => r!(list!(
            first_order = Robj::from(r.first_order),
            total_order = Robj::from(r.total_order),
            var_y = r.var_y,
            component_variance = Robj::from(r.component_variance)
        )),
        Err(_) => r!(NULL),
    }
}

/// Accumulated Local Effects
#[extendr]
fn explain_ale_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
    component: i32,
    n_bins: i32,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::fpc_ale(&fit, &fd, None, component as usize, n_bins as usize) {
        Ok(r) => {
            let bin_counts: Vec<i32> = r.bin_counts.iter().map(|&x| x as i32).collect();
            r!(list!(
                bin_midpoints = Robj::from(r.bin_midpoints),
                ale_values = Robj::from(r.ale_values),
                bin_edges = Robj::from(r.bin_edges),
                bin_counts = Robj::from(bin_counts),
                component = r.component as i32
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Domain selection (important intervals)
#[extendr]
fn explain_domain_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    window_width: i32,
    threshold: f64,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );

    match explain::domain_selection(&fit, window_width as usize, threshold) {
        Ok(r) => {
            let starts: Vec<i32> = r.intervals.iter().map(|iv| iv.start_idx as i32).collect();
            let ends: Vec<i32> = r.intervals.iter().map(|iv| iv.end_idx as i32).collect();
            let importances: Vec<f64> = r.intervals.iter().map(|iv| iv.importance).collect();
            r!(list!(
                pointwise_importance = Robj::from(r.pointwise_importance),
                interval_starts = Robj::from(starts),
                interval_ends = Robj::from(ends),
                interval_importances = Robj::from(importances),
                window_width = r.window_width as i32,
                threshold = r.threshold
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Prediction intervals for new observations
#[extendr]
fn explain_prediction_intervals_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    train_data: RMatrix<f64>,
    new_data: RMatrix<f64>,
    confidence: f64,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let td_s = train_data.as_real_slice().unwrap();
    let fd_train = FdMatrix::from_slice(td_s, train_data.nrows(), train_data.ncols())
        .expect("Invalid train_data");
    let nd_s = new_data.as_real_slice().unwrap();
    let fd_new = FdMatrix::from_slice(nd_s, new_data.nrows(), new_data.ncols())
        .expect("Invalid new_data");

    match explain::prediction_intervals(&fit, &fd_train, None, &fd_new, None, confidence) {
        Ok(r) => r!(list!(
            predictions = Robj::from(r.predictions),
            lower = Robj::from(r.lower),
            upper = Robj::from(r.upper),
            prediction_se = Robj::from(r.prediction_se),
            confidence_level = r.confidence_level,
            t_critical = r.t_critical,
            residual_se = r.residual_se
        )),
        Err(_) => r!(NULL),
    }
}

/// LIME local explanation
#[extendr]
fn explain_lime_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
    observation: i32,
    n_samples: i32,
    kernel_width: f64,
    seed: i32,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::lime_explanation(
        &fit, &fd, None, observation as usize, n_samples as usize, kernel_width, seed as u64,
    ) {
        Ok(r) => r!(list!(
            observation = r.observation as i32,
            attributions = Robj::from(r.attributions),
            local_intercept = r.local_intercept,
            local_r_squared = r.local_r_squared,
            kernel_width = r.kernel_width
        )),
        Err(_) => r!(NULL),
    }
}

/// Counterfactual explanation
#[extendr]
fn explain_counterfactual_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
    observation: i32,
    target_value: f64,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::counterfactual_regression(&fit, &fd, None, observation as usize, target_value) {
        Ok(r) => r!(list!(
            observation = r.observation as i32,
            original_scores = Robj::from(r.original_scores),
            counterfactual_scores = Robj::from(r.counterfactual_scores),
            delta_scores = Robj::from(r.delta_scores),
            delta_function = Robj::from(r.delta_function),
            distance = r.distance,
            original_prediction = r.original_prediction,
            counterfactual_prediction = r.counterfactual_prediction,
            found = r.found
        )),
        Err(_) => r!(NULL),
    }
}

/// Anchor explanation
#[extendr]
fn explain_anchor_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
    observation: i32,
    precision: f64,
    n_bins: i32,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::anchor_explanation(&fit, &fd, None, observation as usize, precision, n_bins as usize) {
        Ok(r) => {
            let components: Vec<i32> = r.rule.conditions.iter().map(|c| c.component as i32).collect();
            let lower_bounds: Vec<f64> = r.rule.conditions.iter().map(|c| c.lower_bound).collect();
            let upper_bounds: Vec<f64> = r.rule.conditions.iter().map(|c| c.upper_bound).collect();
            r!(list!(
                observation = r.observation as i32,
                predicted_value = r.predicted_value,
                precision = r.rule.precision,
                coverage = r.rule.coverage,
                n_matching = r.rule.n_matching as i32,
                condition_components = Robj::from(components),
                condition_lower = Robj::from(lower_bounds),
                condition_upper = Robj::from(upper_bounds)
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Conformal prediction intervals
#[extendr]
fn explain_conformal_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    train_data: RMatrix<f64>,
    train_y: Vec<f64>,
    test_data: RMatrix<f64>,
    cal_fraction: f64,
    alpha: f64,
    seed: i32,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let td_s = train_data.as_real_slice().unwrap();
    let fd_train = FdMatrix::from_slice(td_s, train_data.nrows(), train_data.ncols())
        .expect("Invalid train_data");
    let ts_s = test_data.as_real_slice().unwrap();
    let fd_test = FdMatrix::from_slice(ts_s, test_data.nrows(), test_data.ncols())
        .expect("Invalid test_data");

    match explain::conformal_prediction_residuals(
        &fit, &fd_train, &train_y, &fd_test, None, None, cal_fraction, alpha, seed as u64,
    ) {
        Ok(r) => r!(list!(
            predictions = Robj::from(r.predictions),
            lower = Robj::from(r.lower),
            upper = Robj::from(r.upper),
            residual_quantile = r.residual_quantile,
            coverage = r.coverage,
            calibration_scores = Robj::from(r.calibration_scores)
        )),
        Err(_) => r!(NULL),
    }
}

/// Explanation stability via bootstrap
#[extendr]
fn explain_stability_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    ncomp: i32,
    n_boot: i32,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::explanation_stability(&fd, &y, None, ncomp as usize, n_boot as usize, seed as u64) {
        Ok(r) => r!(list!(
            beta_t_std = Robj::from(r.beta_t_std),
            coefficient_std = Robj::from(r.coefficient_std),
            metric_std = r.metric_std,
            beta_t_cv = Robj::from(r.beta_t_cv),
            importance_stability = r.importance_stability,
            n_boot_success = r.n_boot_success as i32
        )),
        Err(_) => r!(NULL),
    }
}

/// Friedman H-statistic for interaction detection
#[extendr]
fn explain_friedman_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
    component_j: i32,
    component_k: i32,
    n_grid: i32,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::friedman_h_statistic(&fit, &fd, component_j as usize, component_k as usize, n_grid as usize) {
        Ok(r) => r!(list!(
            component_j = r.component_j as i32,
            component_k = r.component_k as i32,
            h_squared = r.h_squared,
            grid_j = Robj::from(r.grid_j),
            grid_k = Robj::from(r.grid_k),
            pdp_2d = fdmatrix_to_robj(&r.pdp_2d)
        )),
        Err(_) => r!(NULL),
    }
}

/// Regression depth
#[extendr]
fn explain_regression_depth_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
    y: Vec<f64>,
    n_boot: i32,
    seed: i32,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::regression_depth(
        &fit, &fd, &y, None, n_boot as usize,
        fdars_core::explain::DepthType::FraimanMuniz, seed as u64,
    ) {
        Ok(r) => r!(list!(
            beta_depth = r.beta_depth,
            score_depths = Robj::from(r.score_depths),
            mean_score_depth = r.mean_score_depth,
            n_boot_success = r.n_boot_success as i32
        )),
        Err(_) => r!(NULL),
    }
}

/// Prototype and criticism selection
#[extendr]
fn explain_prototype_rust(
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    ncomp: i32,
    n_prototypes: i32,
    n_criticisms: i32,
) -> Robj {
    let scores = FdMatrix::from_slice(
        &fpca_scores_data,
        fpca_scores_nrow as usize,
        fpca_scores_ncol as usize,
    )
    .unwrap_or_else(|_| FdMatrix::zeros(0, 0));

    // Build a minimal FpcaResult with just the scores
    let fpca = fdars_core::regression::FpcaResult {
        singular_values: vec![],
        rotation: FdMatrix::zeros(0, 0),
        scores,
        mean: vec![],
        centered: FdMatrix::zeros(0, 0),
    };

    match explain::prototype_criticism(&fpca, ncomp as usize, n_prototypes as usize, n_criticisms as usize) {
        Ok(r) => {
            let proto_idx: Vec<i32> = r.prototype_indices.iter().map(|&x| x as i32).collect();
            let crit_idx: Vec<i32> = r.criticism_indices.iter().map(|&x| x as i32).collect();
            r!(list!(
                prototype_indices = Robj::from(proto_idx),
                prototype_witness = Robj::from(r.prototype_witness),
                criticism_indices = Robj::from(crit_idx),
                criticism_witness = Robj::from(r.criticism_witness),
                bandwidth = r.bandwidth
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Calibration diagnostics for logistic model
#[extendr]
fn explain_calibration_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    probabilities: Vec<f64>,
    predicted_classes: Vec<f64>,
    ncomp: i32,
    accuracy: f64,
    std_errors: Vec<f64>,
    coefficients: Vec<f64>,
    log_likelihood: f64,
    iterations: i32,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    y: Vec<f64>,
    n_groups: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );

    match explain::calibration_diagnostics(&fit, &y, n_groups as usize) {
        Ok(r) => {
            let rel_predicted: Vec<f64> = r.reliability_bins.iter().map(|b| b.0).collect();
            let rel_observed: Vec<f64> = r.reliability_bins.iter().map(|b| b.1).collect();
            let bin_counts: Vec<i32> = r.bin_counts.iter().map(|&x| x as i32).collect();
            r!(list!(
                brier_score = r.brier_score,
                log_loss = r.log_loss,
                hosmer_lemeshow_chi2 = r.hosmer_lemeshow_chi2,
                hosmer_lemeshow_df = r.hosmer_lemeshow_df as i32,
                n_groups = r.n_groups as i32,
                reliability_predicted = Robj::from(rel_predicted),
                reliability_observed = Robj::from(rel_observed),
                bin_counts = Robj::from(bin_counts)
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Expected Calibration Error for logistic model
#[extendr]
fn explain_ece_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    probabilities: Vec<f64>,
    predicted_classes: Vec<f64>,
    ncomp: i32,
    accuracy: f64,
    std_errors: Vec<f64>,
    coefficients: Vec<f64>,
    log_likelihood: f64,
    iterations: i32,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    y: Vec<f64>,
    n_bins: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );

    match explain::expected_calibration_error(&fit, &y, n_bins as usize) {
        Ok(r) => r!(list!(
            ece = r.ece,
            mce = r.mce,
            ace = r.ace,
            n_bins = r.n_bins as i32,
            bin_ece_contributions = Robj::from(r.bin_ece_contributions)
        )),
        Err(_) => r!(NULL),
    }
}

/// Conditional permutation importance
#[extendr]
fn explain_conditional_importance_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    r_squared_adj: f64,
    std_errors: Vec<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    residual_se: f64,
    gcv: f64,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    data: RMatrix<f64>,
    y: Vec<f64>,
    n_bins: i32,
    n_perm: i32,
    seed: i32,
) -> Robj {
    let fit = reconstruct_fregre_lm(
        intercept, beta_t, beta_se, gamma, fitted_values, residuals,
        r_squared, r_squared_adj, std_errors, ncomp, coefficients, residual_se, gcv,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match explain::conditional_permutation_importance(
        &fit, &fd, &y, None, n_bins as usize, n_perm as usize, seed as u64,
    ) {
        Ok(r) => r!(list!(
            importance = Robj::from(r.importance),
            baseline_metric = r.baseline_metric,
            permuted_metric = Robj::from(r.permuted_metric),
            unconditional_importance = Robj::from(r.unconditional_importance)
        )),
        Err(_) => r!(NULL),
    }
}

/// Elastic PCR attribution (amplitude vs phase importance)
#[extendr]
fn elastic_pcr_attribution_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    argvals: Vec<f64>,
    ncomp: i32,
    pca_method: &str,
    lambda: f64,
    max_iter: i32,
    tol: f64,
    n_perm: i32,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    let method = match pca_method {
        "vertical" => fdars_core::elastic_regression::PcaMethod::Vertical,
        "horizontal" => fdars_core::elastic_regression::PcaMethod::Horizontal,
        _ => fdars_core::elastic_regression::PcaMethod::Joint,
    };

    // First fit elastic PCR, then do attribution
    let pcr_result = elastic_regression::elastic_pcr(
        &fd, &y, &argvals, ncomp as usize, method, lambda, max_iter as usize, tol,
    );

    match pcr_result {
        Ok(pcr) => {
            match elastic_explain::elastic_pcr_attribution(&pcr, &y, ncomp as usize, n_perm as usize, seed as u64) {
                Ok(r) => r!(list!(
                    amplitude_contribution = Robj::from(r.amplitude_contribution),
                    phase_contribution = Robj::from(r.phase_contribution),
                    amplitude_importance = r.amplitude_importance,
                    phase_importance = r.phase_importance
                )),
                Err(_) => r!(NULL),
            }
        }
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Logistic explain wrappers (17 functions)
// =============================================================================

/// PDP for logistic model
#[extendr]
fn explain_pdp_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    data: RMatrix<f64>, component: i32, n_grid: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::functional_pdp_logistic(&fit, &fd, None, component as usize, n_grid as usize) {
        Ok(r) => r!(list!(
            grid_values = Robj::from(r.grid_values),
            pdp_curve = Robj::from(r.pdp_curve),
            ice_curves = fdmatrix_to_robj(&r.ice_curves),
            component = r.component as i32
        )),
        Err(_) => r!(NULL),
    }
}

/// Beta decomposition for logistic model
#[extendr]
fn explain_beta_decomposition_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    match explain::beta_decomposition_logistic(&fit) {
        Ok(r) => {
            let nc = r.components.len();
            let m = if nc > 0 { r.components[0].len() } else { 0 };
            let mut flat = vec![0.0; nc * m];
            for i in 0..nc {
                for j in 0..m {
                    flat[j * nc + i] = r.components[i][j];
                }
            }
            let comp_mat = RMatrix::new_matrix(nc, m, |i, j| flat[j * nc + i]);
            r!(list!(
                components = Robj::from(comp_mat),
                coefficients = Robj::from(r.coefficients),
                variance_proportion = Robj::from(r.variance_proportion)
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Pointwise importance for logistic model
#[extendr]
fn explain_pointwise_importance_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    match explain::pointwise_importance_logistic(&fit) {
        Ok(r) => r!(list!(
            importance = Robj::from(r.importance),
            importance_normalized = Robj::from(r.importance_normalized),
            component_importance = fdmatrix_to_robj(&r.component_importance),
            score_variance = Robj::from(r.score_variance)
        )),
        Err(_) => r!(NULL),
    }
}

/// SHAP values for logistic model
#[extendr]
fn explain_shap_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    data: RMatrix<f64>, n_samples: i32, seed: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::fpc_shap_values_logistic(&fit, &fd, None, n_samples as usize, seed as u64) {
        Ok(r) => r!(list!(
            values = fdmatrix_to_robj(&r.values),
            base_value = r.base_value,
            mean_scores = Robj::from(r.mean_scores)
        )),
        Err(_) => r!(NULL),
    }
}

/// VIF for logistic model
#[extendr]
fn explain_vif_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    data: RMatrix<f64>,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::fpc_vif_logistic(&fit, &fd, None) {
        Ok(r) => r!(list!(
            vif = Robj::from(r.vif),
            labels = Robj::from(r.labels),
            mean_vif = r.mean_vif,
            n_moderate = r.n_moderate as i32,
            n_severe = r.n_severe as i32
        )),
        Err(_) => r!(NULL),
    }
}

/// Permutation importance for logistic model
#[extendr]
fn explain_importance_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    data: RMatrix<f64>, y: Vec<f64>, n_perm: i32, seed: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::fpc_permutation_importance_logistic(&fit, &fd, &y, n_perm as usize, seed as u64) {
        Ok(r) => r!(list!(
            importance = Robj::from(r.importance),
            baseline_metric = r.baseline_metric,
            permuted_metric = Robj::from(r.permuted_metric)
        )),
        Err(_) => r!(NULL),
    }
}

/// Sobol indices for logistic model (no y, uses n_samples/seed)
#[extendr]
fn explain_sobol_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    data: RMatrix<f64>, n_samples: i32, seed: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::sobol_indices_logistic(&fit, &fd, None, n_samples as usize, seed as u64) {
        Ok(r) => r!(list!(
            first_order = Robj::from(r.first_order),
            total_order = Robj::from(r.total_order),
            var_y = r.var_y,
            component_variance = Robj::from(r.component_variance)
        )),
        Err(_) => r!(NULL),
    }
}

/// ALE for logistic model
#[extendr]
fn explain_ale_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    data: RMatrix<f64>, component: i32, n_bins: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::fpc_ale_logistic(&fit, &fd, None, component as usize, n_bins as usize) {
        Ok(r) => {
            let bin_counts: Vec<i32> = r.bin_counts.iter().map(|&x| x as i32).collect();
            r!(list!(
                bin_midpoints = Robj::from(r.bin_midpoints),
                ale_values = Robj::from(r.ale_values),
                bin_edges = Robj::from(r.bin_edges),
                bin_counts = Robj::from(bin_counts),
                component = r.component as i32
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Domain selection for logistic model
#[extendr]
fn explain_domain_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    window_width: i32, threshold: f64,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    match explain::domain_selection_logistic(&fit, window_width as usize, threshold) {
        Ok(r) => {
            let starts: Vec<i32> = r.intervals.iter().map(|iv| iv.start_idx as i32).collect();
            let ends: Vec<i32> = r.intervals.iter().map(|iv| iv.end_idx as i32).collect();
            let importances: Vec<f64> = r.intervals.iter().map(|iv| iv.importance).collect();
            r!(list!(
                pointwise_importance = Robj::from(r.pointwise_importance),
                interval_starts = Robj::from(starts),
                interval_ends = Robj::from(ends),
                interval_importances = Robj::from(importances),
                window_width = r.window_width as i32,
                threshold = r.threshold
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Saliency map for logistic model (fit-only, no data)
#[extendr]
fn explain_saliency_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    match explain::functional_saliency_logistic(&fit) {
        Ok(r) => r!(list!(
            saliency_map = fdmatrix_to_robj(&r.saliency_map),
            mean_absolute_saliency = Robj::from(r.mean_absolute_saliency)
        )),
        Err(_) => r!(NULL),
    }
}

/// LIME for logistic model
#[extendr]
fn explain_lime_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    data: RMatrix<f64>, observation: i32, n_samples: i32, kernel_width: f64, seed: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::lime_explanation_logistic(
        &fit, &fd, None, observation as usize, n_samples as usize, kernel_width, seed as u64,
    ) {
        Ok(r) => r!(list!(
            observation = r.observation as i32,
            attributions = Robj::from(r.attributions),
            local_intercept = r.local_intercept,
            local_r_squared = r.local_r_squared,
            kernel_width = r.kernel_width
        )),
        Err(_) => r!(NULL),
    }
}

/// Counterfactual for logistic model (max_iter/step_size instead of target_value)
#[extendr]
fn explain_counterfactual_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    data: RMatrix<f64>, observation: i32, max_iter: i32, step_size: f64,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::counterfactual_logistic(
        &fit, &fd, None, observation as usize, max_iter as usize, step_size,
    ) {
        Ok(r) => r!(list!(
            observation = r.observation as i32,
            original_scores = Robj::from(r.original_scores),
            counterfactual_scores = Robj::from(r.counterfactual_scores),
            delta_scores = Robj::from(r.delta_scores),
            delta_function = Robj::from(r.delta_function),
            distance = r.distance,
            original_prediction = r.original_prediction,
            counterfactual_prediction = r.counterfactual_prediction,
            found = r.found
        )),
        Err(_) => r!(NULL),
    }
}

/// Anchor for logistic model
#[extendr]
fn explain_anchor_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    data: RMatrix<f64>, observation: i32, precision: f64, n_bins: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::anchor_explanation_logistic(
        &fit, &fd, None, observation as usize, precision, n_bins as usize,
    ) {
        Ok(r) => {
            let components: Vec<i32> = r.rule.conditions.iter().map(|c| c.component as i32).collect();
            let lower_bounds: Vec<f64> = r.rule.conditions.iter().map(|c| c.lower_bound).collect();
            let upper_bounds: Vec<f64> = r.rule.conditions.iter().map(|c| c.upper_bound).collect();
            r!(list!(
                observation = r.observation as i32,
                predicted_value = r.predicted_value,
                precision = r.rule.precision,
                coverage = r.rule.coverage,
                n_matching = r.rule.n_matching as i32,
                condition_components = Robj::from(components),
                condition_lower = Robj::from(lower_bounds),
                condition_upper = Robj::from(upper_bounds)
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Stability for logistic model (no fit param)
#[extendr]
fn explain_stability_logistic_rust(
    data: RMatrix<f64>, y: Vec<f64>,
    ncomp: i32, n_boot: i32, seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::explanation_stability_logistic(&fd, &y, None, ncomp as usize, n_boot as usize, seed as u64) {
        Ok(r) => r!(list!(
            beta_t_std = Robj::from(r.beta_t_std),
            coefficient_std = Robj::from(r.coefficient_std),
            metric_std = r.metric_std,
            beta_t_cv = Robj::from(r.beta_t_cv),
            importance_stability = r.importance_stability,
            n_boot_success = r.n_boot_success as i32
        )),
        Err(_) => r!(NULL),
    }
}

/// Friedman H-statistic for logistic model
#[extendr]
fn explain_friedman_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    data: RMatrix<f64>, component_j: i32, component_k: i32, n_grid: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::friedman_h_statistic_logistic(
        &fit, &fd, None, component_j as usize, component_k as usize, n_grid as usize,
    ) {
        Ok(r) => r!(list!(
            component_j = r.component_j as i32,
            component_k = r.component_k as i32,
            h_squared = r.h_squared,
            grid_j = Robj::from(r.grid_j),
            grid_k = Robj::from(r.grid_k),
            pdp_2d = fdmatrix_to_robj(&r.pdp_2d)
        )),
        Err(_) => r!(NULL),
    }
}

/// Regression depth for logistic model
#[extendr]
fn explain_regression_depth_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    data: RMatrix<f64>, y: Vec<f64>, n_boot: i32, seed: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::regression_depth_logistic(
        &fit, &fd, &y, None, n_boot as usize,
        fdars_core::explain::DepthType::FraimanMuniz, seed as u64,
    ) {
        Ok(r) => r!(list!(
            beta_depth = r.beta_depth,
            score_depths = Robj::from(r.score_depths),
            mean_score_depth = r.mean_score_depth,
            n_boot_success = r.n_boot_success as i32
        )),
        Err(_) => r!(NULL),
    }
}

/// Conditional permutation importance for logistic model
#[extendr]
fn explain_conditional_importance_logistic_rust(
    intercept: f64, beta_t: Vec<f64>, beta_se: Vec<f64>, gamma: Vec<f64>,
    probabilities: Vec<f64>, predicted_classes: Vec<f64>,
    ncomp: i32, accuracy: f64, std_errors: Vec<f64>,
    coefficients: Vec<f64>, log_likelihood: f64, iterations: i32,
    fpca_mean: Vec<f64>, fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32, fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>, fpca_scores_nrow: i32, fpca_scores_ncol: i32,
    data: RMatrix<f64>, y: Vec<f64>, n_bins: i32, n_perm: i32, seed: i32,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match explain::conditional_permutation_importance_logistic(
        &fit, &fd, &y, None, n_bins as usize, n_perm as usize, seed as u64,
    ) {
        Ok(r) => r!(list!(
            importance = Robj::from(r.importance),
            baseline_metric = r.baseline_metric,
            permuted_metric = Robj::from(r.permuted_metric),
            unconditional_importance = Robj::from(r.unconditional_importance)
        )),
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Elastic regression wrappers (3 functions)
// =============================================================================

/// Elastic regression (function-on-scalar with alignment)
#[extendr]
fn elastic_regression_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    argvals: Vec<f64>,
    ncomp_beta: i32,
    lambda: f64,
    max_iter: i32,
    tol: f64,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match elastic_regression::elastic_regression(&fd, &y, &argvals, ncomp_beta as usize, lambda, max_iter as usize, tol) {
        Ok(r) => r!(list!(
            alpha = r.alpha,
            beta = Robj::from(r.beta),
            fitted_values = Robj::from(r.fitted_values),
            residuals = Robj::from(r.residuals),
            sse = r.sse,
            r_squared = r.r_squared,
            gammas = fdmatrix_to_robj(&r.gammas),
            aligned_srsfs = fdmatrix_to_robj(&r.aligned_srsfs),
            n_iter = r.n_iter as i32
        )),
        Err(_) => r!(NULL),
    }
}

/// Elastic logistic regression
#[extendr]
fn elastic_logistic_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    argvals: Vec<f64>,
    ncomp_beta: i32,
    lambda: f64,
    max_iter: i32,
    tol: f64,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    let y_i8: Vec<i8> = y.iter().map(|&v| v as i8).collect();

    match elastic_regression::elastic_logistic(&fd, &y_i8, &argvals, ncomp_beta as usize, lambda, max_iter as usize, tol) {
        Ok(r) => {
            let pred_classes: Vec<i32> = r.predicted_classes.iter().map(|&x| x as i32).collect();
            r!(list!(
                alpha = r.alpha,
                beta = Robj::from(r.beta),
                probabilities = Robj::from(r.probabilities),
                predicted_classes = Robj::from(pred_classes),
                accuracy = r.accuracy,
                loss = r.loss,
                gammas = fdmatrix_to_robj(&r.gammas),
                aligned_srsfs = fdmatrix_to_robj(&r.aligned_srsfs),
                n_iter = r.n_iter as i32
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Elastic principal component regression
#[extendr]
fn elastic_pcr_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    argvals: Vec<f64>,
    ncomp: i32,
    pca_method: &str,
    lambda: f64,
    max_iter: i32,
    tol: f64,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    let method = match pca_method {
        "vertical" => fdars_core::elastic_regression::PcaMethod::Vertical,
        "horizontal" => fdars_core::elastic_regression::PcaMethod::Horizontal,
        _ => fdars_core::elastic_regression::PcaMethod::Joint,
    };

    match elastic_regression::elastic_pcr(&fd, &y, &argvals, ncomp as usize, method, lambda, max_iter as usize, tol) {
        Ok(r) => {
            // Return core results; FPCA sub-results are complex so return key summary fields
            r!(list!(
                alpha = r.alpha,
                coefficients = Robj::from(r.coefficients),
                fitted_values = Robj::from(r.fitted_values),
                sse = r.sse,
                r_squared = r.r_squared,
                karcher_mean = Robj::from(r.karcher.mean.clone()),
                karcher_mean_srsf = Robj::from(r.karcher.mean_srsf.clone()),
                karcher_gammas = fdmatrix_to_robj(&r.karcher.gammas),
                karcher_aligned = fdmatrix_to_robj(&r.karcher.aligned_data)
            ))
        }
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Elastic FPCA wrappers (3 functions)
// =============================================================================

/// Vertical FPCA on elastic-aligned data
#[extendr]
fn elastic_vert_fpca_rust(
    aligned_data: RMatrix<f64>,
    gammas_data: RMatrix<f64>,
    mean_curve: Vec<f64>,
    mean_srsf: Vec<f64>,
    aligned_srsfs: Nullable<RMatrix<f64>>,
    argvals: Vec<f64>,
    ncomp: i32,
) -> Robj {
    let karcher = reconstruct_karcher_mean(aligned_data, gammas_data, mean_curve, mean_srsf, aligned_srsfs, argvals.clone());

    match elastic_fpca::vert_fpca(&karcher, &argvals, ncomp as usize) {
        Ok(r) => r!(list!(
            scores = fdmatrix_to_robj(&r.scores),
            eigenfunctions_q = fdmatrix_to_robj(&r.eigenfunctions_q),
            eigenfunctions_f = fdmatrix_to_robj(&r.eigenfunctions_f),
            eigenvalues = Robj::from(r.eigenvalues),
            cumulative_variance = Robj::from(r.cumulative_variance),
            mean_q = Robj::from(r.mean_q)
        )),
        Err(_) => r!(NULL),
    }
}

/// Horizontal FPCA on warping functions
#[extendr]
fn elastic_horiz_fpca_rust(
    aligned_data: RMatrix<f64>,
    gammas_data: RMatrix<f64>,
    mean_curve: Vec<f64>,
    mean_srsf: Vec<f64>,
    aligned_srsfs: Nullable<RMatrix<f64>>,
    argvals: Vec<f64>,
    ncomp: i32,
) -> Robj {
    let karcher = reconstruct_karcher_mean(aligned_data, gammas_data, mean_curve, mean_srsf, aligned_srsfs, argvals.clone());

    match elastic_fpca::horiz_fpca(&karcher, &argvals, ncomp as usize) {
        Ok(r) => r!(list!(
            scores = fdmatrix_to_robj(&r.scores),
            eigenfunctions_psi = fdmatrix_to_robj(&r.eigenfunctions_psi),
            eigenfunctions_gam = fdmatrix_to_robj(&r.eigenfunctions_gam),
            eigenvalues = Robj::from(r.eigenvalues),
            cumulative_variance = Robj::from(r.cumulative_variance),
            mean_psi = Robj::from(r.mean_psi),
            shooting_vectors = fdmatrix_to_robj(&r.shooting_vectors)
        )),
        Err(_) => r!(NULL),
    }
}

/// Joint FPCA (amplitude + phase)
#[extendr]
fn elastic_joint_fpca_rust(
    aligned_data: RMatrix<f64>,
    gammas_data: RMatrix<f64>,
    mean_curve: Vec<f64>,
    mean_srsf: Vec<f64>,
    aligned_srsfs: Nullable<RMatrix<f64>>,
    argvals: Vec<f64>,
    ncomp: i32,
) -> Robj {
    let karcher = reconstruct_karcher_mean(aligned_data, gammas_data, mean_curve, mean_srsf, aligned_srsfs, argvals.clone());

    match elastic_fpca::joint_fpca(&karcher, &argvals, ncomp as usize, None) {
        Ok(r) => r!(list!(
            scores = fdmatrix_to_robj(&r.scores),
            eigenvalues = Robj::from(r.eigenvalues),
            cumulative_variance = Robj::from(r.cumulative_variance),
            balance_c = r.balance_c,
            vert_component = fdmatrix_to_robj(&r.vert_component),
            horiz_component = fdmatrix_to_robj(&r.horiz_component)
        )),
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Elastic changepoint wrappers (3 functions)
// =============================================================================

/// Amplitude changepoint detection
#[extendr]
fn elastic_amp_changepoint_rust(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    lambda: f64,
    max_iter: i32,
    n_mc: i32,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match elastic_changepoint::elastic_amp_changepoint(
        &fd, &argvals, lambda, max_iter as usize, n_mc as usize, seed as u64,
    ) {
        Ok(r) => r!(list!(
            changepoint = r.changepoint as i32,
            test_statistic = r.test_statistic,
            p_value = r.p_value,
            cusum_values = Robj::from(r.cusum_values)
        )),
        Err(_) => r!(NULL),
    }
}

/// Phase changepoint detection
#[extendr]
fn elastic_ph_changepoint_rust(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    lambda: f64,
    max_iter: i32,
    n_mc: i32,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    match elastic_changepoint::elastic_ph_changepoint(
        &fd, &argvals, lambda, max_iter as usize, n_mc as usize, seed as u64,
    ) {
        Ok(r) => r!(list!(
            changepoint = r.changepoint as i32,
            test_statistic = r.test_statistic,
            p_value = r.p_value,
            cusum_values = Robj::from(r.cusum_values)
        )),
        Err(_) => r!(NULL),
    }
}

/// FPCA-based changepoint detection
#[extendr]
fn elastic_fpca_changepoint_rust(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    pca_method: &str,
    ncomp: i32,
    lambda: f64,
    max_iter: i32,
    n_mc: i32,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    let method = match pca_method {
        "horizontal" => fdars_core::elastic_changepoint::FpcaChangepointMethod::Horizontal,
        "joint" => fdars_core::elastic_changepoint::FpcaChangepointMethod::Joint,
        _ => fdars_core::elastic_changepoint::FpcaChangepointMethod::Vertical,
    };

    match elastic_changepoint::elastic_fpca_changepoint(
        &fd, &argvals, method, ncomp as usize, lambda, max_iter as usize, n_mc as usize, seed as u64,
    ) {
        Ok(r) => r!(list!(
            changepoint = r.changepoint as i32,
            test_statistic = r.test_statistic,
            p_value = r.p_value,
            cusum_values = Robj::from(r.cusum_values)
        )),
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Smooth basis wrappers (2 functions)
// =============================================================================

/// Smooth functional data using B-spline or Fourier basis with fixed penalty
#[extendr]
fn smooth_basis_rust(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    basis_type: &str,
    nbasis: i32,
    lambda: f64,
    lfd_order: i32,
    period: f64,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    let bt = match basis_type {
        "fourier" => smooth_basis::BasisType::Fourier { period },
        _ => smooth_basis::BasisType::Bspline { order: 4 },
    };

    let pen = match &bt {
        smooth_basis::BasisType::Bspline { order } => {
            smooth_basis::bspline_penalty_matrix(&argvals, nbasis as usize, *order, lfd_order as usize)
        }
        smooth_basis::BasisType::Fourier { period } => {
            smooth_basis::fourier_penalty_matrix(nbasis as usize, *period, lfd_order as usize)
        }
    };

    let fdpar = smooth_basis::FdPar {
        basis_type: bt,
        nbasis: nbasis as usize,
        lambda,
        lfd_order: lfd_order as usize,
        penalty_matrix: pen,
    };

    match smooth_basis::smooth_basis(&fd, &argvals, &fdpar) {
        Ok(r) => r!(list!(
            coefficients = fdmatrix_to_robj(&r.coefficients),
            fitted = fdmatrix_to_robj(&r.fitted),
            edf = r.edf,
            gcv = r.gcv,
            aic = r.aic,
            bic = r.bic,
            nbasis = r.nbasis as i32
        )),
        Err(_) => r!(NULL),
    }
}

/// Smooth functional data with GCV-selected penalty
#[extendr]
fn smooth_basis_gcv_rust(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    basis_type: &str,
    nbasis: i32,
    lfd_order: i32,
    log_lambda_min: f64,
    log_lambda_max: f64,
    n_grid: i32,
    period: f64,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    let bt = match basis_type {
        "fourier" => smooth_basis::BasisType::Fourier { period },
        _ => smooth_basis::BasisType::Bspline { order: 4 },
    };

    match smooth_basis::smooth_basis_gcv(
        &fd, &argvals, &bt, nbasis as usize, lfd_order as usize,
        (log_lambda_min, log_lambda_max), n_grid as usize,
    ) {
        Some(r) => r!(list!(
            coefficients = fdmatrix_to_robj(&r.coefficients),
            fitted = fdmatrix_to_robj(&r.fitted),
            edf = r.edf,
            gcv = r.gcv,
            aic = r.aic,
            bic = r.bic,
            nbasis = r.nbasis as i32
        )),
        None => r!(NULL),
    }
}

// =============================================================================
// Phase B1: 2D Function-on-Scalar Regression
// =============================================================================

/// 2D Function-on-scalar regression with tensor-product penalties
#[extendr]
fn fosr_2d_rust(
    data: RMatrix<f64>,
    predictors: RMatrix<f64>,
    argvals_s: Vec<f64>,
    argvals_t: Vec<f64>,
    lambda_s: f64,
    lambda_t: f64,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    let ps = predictors.as_real_slice().unwrap();
    let pred = FdMatrix::from_slice(ps, predictors.nrows(), predictors.ncols()).expect("Invalid predictors");

    let grid = fdars_core::function_on_scalar_2d::Grid2d::new(argvals_s, argvals_t);

    match fdars_core::function_on_scalar_2d::fosr_2d(&fd, &pred, &grid, lambda_s, lambda_t) {
        Ok(r) => {
            let beta_se_robj = match &r.beta_se {
                Some(bse) => fdmatrix_to_robj(bse),
                None => r!(NULL),
            };
            r!(list!(
                intercept = Robj::from(r.intercept),
                beta = fdmatrix_to_robj(&r.beta),
                fitted = fdmatrix_to_robj(&r.fitted),
                residuals = fdmatrix_to_robj(&r.residuals),
                r_squared_pointwise = Robj::from(r.r_squared_pointwise),
                r_squared = r.r_squared,
                beta_se = beta_se_robj,
                lambda_s = r.lambda_s,
                lambda_t = r.lambda_t,
                gcv = r.gcv,
                grid_m1 = r.grid.m1() as i32,
                grid_m2 = r.grid.m2() as i32
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Predict from 2D function-on-scalar regression
#[extendr]
fn predict_fosr_2d_rust(
    intercept: Vec<f64>,
    beta_data: RMatrix<f64>,
    new_predictors: RMatrix<f64>,
    argvals_s: Vec<f64>,
    argvals_t: Vec<f64>,
    lambda_s: f64,
    lambda_t: f64,
) -> Robj {
    let bs = beta_data.as_real_slice().unwrap();
    let beta = FdMatrix::from_slice(bs, beta_data.nrows(), beta_data.ncols()).expect("Invalid beta");

    let ps = new_predictors.as_real_slice().unwrap();
    let pred = FdMatrix::from_slice(ps, new_predictors.nrows(), new_predictors.ncols()).expect("Invalid predictors");

    let grid = fdars_core::function_on_scalar_2d::Grid2d::new(argvals_s, argvals_t);
    let m = grid.total();

    let fosr_result = fdars_core::function_on_scalar_2d::FosrResult2d::from_parts(
        intercept.clone(),
        beta.clone(),
        FdMatrix::zeros(1, m),
        FdMatrix::zeros(1, m),
        vec![0.0; m],
        0.0,
        None,
        lambda_s,
        lambda_t,
        0.0,
        grid,
    );

    match fdars_core::function_on_scalar_2d::predict_fosr_2d(&fosr_result, &pred) {
        Ok(fitted) => fdmatrix_to_robj(&fitted),
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Phase B2: Predict Methods
// =============================================================================

/// Predict from functional logistic regression
#[extendr]
fn predict_functional_logistic_rust(
    intercept: f64,
    beta_t: Vec<f64>,
    beta_se: Vec<f64>,
    gamma: Vec<f64>,
    probabilities: Vec<f64>,
    predicted_classes: Vec<f64>,
    ncomp: i32,
    accuracy: f64,
    std_errors: Vec<f64>,
    coefficients: Vec<f64>,
    log_likelihood: f64,
    iterations: i32,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    fpca_scores_data: Vec<f64>,
    fpca_scores_nrow: i32,
    fpca_scores_ncol: i32,
    new_data: RMatrix<f64>,
    new_scalar: Nullable<RMatrix<f64>>,
) -> Robj {
    let fit = reconstruct_functional_logistic(
        intercept, beta_t, beta_se, gamma, probabilities, predicted_classes,
        ncomp, accuracy, std_errors, coefficients, log_likelihood, iterations,
        fpca_mean, fpca_rotation_data, fpca_rotation_nrow, fpca_rotation_ncol,
        fpca_scores_data, fpca_scores_nrow, fpca_scores_ncol,
    );

    let nds = new_data.as_real_slice().unwrap();
    let nfd = FdMatrix::from_slice(nds, new_data.nrows(), new_data.ncols()).expect("Invalid new_data");

    let ns = match new_scalar {
        Nullable::NotNull(s) => {
            let ss = s.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid new_scalar"))
        }
        Nullable::Null => None,
    };

    let probs = scalar_on_function::predict_functional_logistic(&fit, &nfd, ns.as_ref());
    Robj::from(probs)
}

/// Predict from elastic regression
#[extendr]
fn predict_elastic_regression_rust(
    alpha: f64,
    beta: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    sse: f64,
    r_squared: f64,
    gammas_data: RMatrix<f64>,
    aligned_srsfs_data: RMatrix<f64>,
    n_iter: i32,
    new_data: RMatrix<f64>,
    argvals: Vec<f64>,
) -> Robj {
    let gs = gammas_data.as_real_slice().unwrap();
    let gammas = FdMatrix::from_slice(gs, gammas_data.nrows(), gammas_data.ncols()).expect("Invalid gammas");

    let as_ = aligned_srsfs_data.as_real_slice().unwrap();
    let aligned_srsfs = FdMatrix::from_slice(as_, aligned_srsfs_data.nrows(), aligned_srsfs_data.ncols()).expect("Invalid aligned_srsfs");

    let fit = elastic_regression::ElasticRegressionResult::from_parts(
        alpha,
        beta,
        fitted_values,
        residuals,
        sse,
        r_squared,
        gammas,
        aligned_srsfs,
        n_iter as usize,
    );

    let nds = new_data.as_real_slice().unwrap();
    let nfd = FdMatrix::from_slice(nds, new_data.nrows(), new_data.ncols()).expect("Invalid new_data");

    let preds = elastic_regression::predict_elastic_regression(&fit, &nfd, &argvals);
    Robj::from(preds)
}

/// Predict from elastic logistic regression
#[extendr]
fn predict_elastic_logistic_rust(
    alpha: f64,
    beta: Vec<f64>,
    probabilities: Vec<f64>,
    predicted_classes: Vec<i32>,
    accuracy: f64,
    loss: f64,
    gammas_data: RMatrix<f64>,
    aligned_srsfs_data: RMatrix<f64>,
    n_iter: i32,
    new_data: RMatrix<f64>,
    argvals: Vec<f64>,
) -> Robj {
    let gs = gammas_data.as_real_slice().unwrap();
    let gammas = FdMatrix::from_slice(gs, gammas_data.nrows(), gammas_data.ncols()).expect("Invalid gammas");

    let as_ = aligned_srsfs_data.as_real_slice().unwrap();
    let aligned_srsfs = FdMatrix::from_slice(as_, aligned_srsfs_data.nrows(), aligned_srsfs_data.ncols()).expect("Invalid aligned_srsfs");

    let pred_usize: Vec<usize> = predicted_classes.iter().map(|&x| x as usize).collect();
    let fit = elastic_regression::ElasticLogisticResult::from_parts(
        alpha,
        beta,
        probabilities,
        pred_usize,
        accuracy,
        loss,
        gammas,
        aligned_srsfs,
        n_iter as usize,
    );

    let nds = new_data.as_real_slice().unwrap();
    let nfd = FdMatrix::from_slice(nds, new_data.nrows(), new_data.ncols()).expect("Invalid new_data");

    let probs = elastic_regression::predict_elastic_logistic(&fit, &nfd, &argvals);
    Robj::from(probs)
}

// =============================================================================
// Phase B3: Model Selection
// =============================================================================

/// Model selection for ncomp via AIC/BIC/GCV
#[extendr]
fn model_selection_ncomp_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    scalar_covariates: Nullable<RMatrix<f64>>,
    max_ncomp: i32,
    criterion: &str,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");

    let scov = match scalar_covariates {
        Nullable::NotNull(sc) => {
            let sc_s = sc.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(sc_s, sc.nrows(), sc.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    let crit = match criterion {
        "bic" => scalar_on_function::SelectionCriterion::Bic,
        "gcv" => scalar_on_function::SelectionCriterion::Gcv,
        _ => scalar_on_function::SelectionCriterion::Aic,
    };

    match scalar_on_function::model_selection_ncomp(&fd, &y, scov.as_ref(), max_ncomp as usize, crit) {
        Ok(r) => {
            let ncomp_vals: Vec<i32> = r.criteria.iter().map(|&(k, _, _, _)| k as i32).collect();
            let aic_vals: Vec<f64> = r.criteria.iter().map(|&(_, a, _, _)| a).collect();
            let bic_vals: Vec<f64> = r.criteria.iter().map(|&(_, _, b, _)| b).collect();
            let gcv_vals: Vec<f64> = r.criteria.iter().map(|&(_, _, _, g)| g).collect();
            r!(list!(
                best_ncomp = r.best_ncomp as i32,
                ncomp = Robj::from(ncomp_vals),
                aic = Robj::from(aic_vals),
                bic = Robj::from(bic_vals),
                gcv = Robj::from(gcv_vals)
            ))
        }
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Phase B4: Conformal Prediction
// =============================================================================

/// Conformal prediction for functional linear model
#[extendr]
fn conformal_fregre_lm_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    test_data: RMatrix<f64>,
    scalar_train: Nullable<RMatrix<f64>>,
    scalar_test: Nullable<RMatrix<f64>>,
    ncomp: i32,
    cal_fraction: f64,
    alpha: f64,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");

    let st = match scalar_train {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_train")) }
        Nullable::Null => None,
    };
    let ste = match scalar_test {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_test")) }
        Nullable::Null => None,
    };

    match fdars_core::conformal::conformal_fregre_lm(&fd, &y, &tfd, st.as_ref(), ste.as_ref(), ncomp as usize, cal_fraction, alpha, seed as u64) {
        Ok(r) => r!(list!(
            predictions = Robj::from(r.predictions),
            lower = Robj::from(r.lower),
            upper = Robj::from(r.upper),
            residual_quantile = r.residual_quantile,
            coverage = r.coverage,
            calibration_scores = Robj::from(r.calibration_scores)
        )),
        Err(_) => r!(NULL),
    }
}

/// Conformal prediction for nonparametric functional regression
#[extendr]
fn conformal_fregre_np_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    test_data: RMatrix<f64>,
    argvals: Vec<f64>,
    scalar_train: Nullable<RMatrix<f64>>,
    scalar_test: Nullable<RMatrix<f64>>,
    h_func: f64,
    h_scalar: f64,
    cal_fraction: f64,
    alpha: f64,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");

    let st = match scalar_train {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_train")) }
        Nullable::Null => None,
    };
    let ste = match scalar_test {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_test")) }
        Nullable::Null => None,
    };

    match fdars_core::conformal::conformal_fregre_np(&fd, &y, &tfd, &argvals, st.as_ref(), ste.as_ref(), h_func, h_scalar, cal_fraction, alpha, seed as u64) {
        Ok(r) => r!(list!(
            predictions = Robj::from(r.predictions),
            lower = Robj::from(r.lower),
            upper = Robj::from(r.upper),
            residual_quantile = r.residual_quantile,
            coverage = r.coverage,
            calibration_scores = Robj::from(r.calibration_scores)
        )),
        Err(_) => r!(NULL),
    }
}

/// Conformal prediction for elastic regression
#[extendr]
fn conformal_elastic_regression_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    test_data: RMatrix<f64>,
    argvals: Vec<f64>,
    ncomp_beta: i32,
    lambda: f64,
    cal_fraction: f64,
    alpha: f64,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");

    match fdars_core::conformal::conformal_elastic_regression(&fd, &y, &tfd, &argvals, ncomp_beta as usize, lambda, cal_fraction, alpha, seed as u64) {
        Ok(r) => r!(list!(
            predictions = Robj::from(r.predictions),
            lower = Robj::from(r.lower),
            upper = Robj::from(r.upper),
            residual_quantile = r.residual_quantile,
            coverage = r.coverage,
            calibration_scores = Robj::from(r.calibration_scores)
        )),
        Err(_) => r!(NULL),
    }
}

/// Conformal prediction for elastic PCR
#[extendr]
fn conformal_elastic_pcr_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    test_data: RMatrix<f64>,
    argvals: Vec<f64>,
    ncomp: i32,
    pca_method: &str,
    lambda: f64,
    cal_fraction: f64,
    alpha: f64,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");

    let pm = match pca_method {
        "horizontal" => elastic_regression::PcaMethod::Horizontal,
        "joint" => elastic_regression::PcaMethod::Joint,
        _ => elastic_regression::PcaMethod::Vertical,
    };

    match fdars_core::conformal::conformal_elastic_pcr(&fd, &y, &tfd, &argvals, ncomp as usize, pm, lambda, cal_fraction, alpha, seed as u64) {
        Ok(r) => r!(list!(
            predictions = Robj::from(r.predictions),
            lower = Robj::from(r.lower),
            upper = Robj::from(r.upper),
            residual_quantile = r.residual_quantile,
            coverage = r.coverage,
            calibration_scores = Robj::from(r.calibration_scores)
        )),
        Err(_) => r!(NULL),
    }
}

/// Conformal classification
#[extendr]
fn conformal_classif_rust(
    data: RMatrix<f64>,
    y: Vec<i32>,
    test_data: RMatrix<f64>,
    covariates_train: Nullable<RMatrix<f64>>,
    covariates_test: Nullable<RMatrix<f64>>,
    ncomp: i32,
    classifier: &str,
    k_nn: i32,
    score_type: &str,
    cal_fraction: f64,
    alpha: f64,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");
    let y_usize: Vec<usize> = y.iter().map(|&c| c as usize).collect();

    let ct = match covariates_train {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid covariates_train")) }
        Nullable::Null => None,
    };
    let cte = match covariates_test {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid covariates_test")) }
        Nullable::Null => None,
    };

    let st = match score_type {
        "aps" => fdars_core::conformal::ClassificationScore::Aps,
        _ => fdars_core::conformal::ClassificationScore::Lac,
    };

    match fdars_core::conformal::conformal_classif(&fd, &y_usize, &tfd, ct.as_ref(), cte.as_ref(), ncomp as usize, classifier, k_nn as usize, st, cal_fraction, alpha, seed as u64) {
        Ok(r) => {
            let pred: Vec<i32> = r.predicted_classes.iter().map(|&c| c as i32).collect();
            let set_sizes: Vec<i32> = r.set_sizes.iter().map(|&s| s as i32).collect();
            r!(list!(
                predicted_classes = Robj::from(pred),
                set_sizes = Robj::from(set_sizes),
                average_set_size = r.average_set_size,
                coverage = r.coverage,
                score_quantile = r.score_quantile
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Conformal prediction for logistic regression
#[extendr]
fn conformal_logistic_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    test_data: RMatrix<f64>,
    scalar_train: Nullable<RMatrix<f64>>,
    scalar_test: Nullable<RMatrix<f64>>,
    ncomp: i32,
    max_iter: i32,
    tol: f64,
    score_type: &str,
    cal_fraction: f64,
    alpha: f64,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");

    let st = match scalar_train {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_train")) }
        Nullable::Null => None,
    };
    let ste = match scalar_test {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_test")) }
        Nullable::Null => None,
    };

    let score = match score_type {
        "aps" => fdars_core::conformal::ClassificationScore::Aps,
        _ => fdars_core::conformal::ClassificationScore::Lac,
    };

    match fdars_core::conformal::conformal_logistic(&fd, &y, &tfd, st.as_ref(), ste.as_ref(), ncomp as usize, max_iter as usize, tol, score, cal_fraction, alpha, seed as u64) {
        Ok(r) => {
            let pred: Vec<i32> = r.predicted_classes.iter().map(|&c| c as i32).collect();
            let set_sizes: Vec<i32> = r.set_sizes.iter().map(|&s| s as i32).collect();
            r!(list!(
                predicted_classes = Robj::from(pred),
                set_sizes = Robj::from(set_sizes),
                average_set_size = r.average_set_size,
                coverage = r.coverage,
                score_quantile = r.score_quantile
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Conformal prediction for elastic logistic regression
#[extendr]
fn conformal_elastic_logistic_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    test_data: RMatrix<f64>,
    argvals: Vec<f64>,
    lambda: f64,
    score_type: &str,
    cal_fraction: f64,
    alpha: f64,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");

    let y_i8: Vec<i8> = y.iter().map(|&v| v as i8).collect();
    let score = match score_type {
        "aps" => fdars_core::conformal::ClassificationScore::Aps,
        _ => fdars_core::conformal::ClassificationScore::Lac,
    };

    match fdars_core::conformal::conformal_elastic_logistic(&fd, &y_i8, &tfd, &argvals, lambda, score, cal_fraction, alpha, seed as u64) {
        Ok(r) => {
            let pred: Vec<i32> = r.predicted_classes.iter().map(|&c| c as i32).collect();
            let set_sizes: Vec<i32> = r.set_sizes.iter().map(|&s| s as i32).collect();
            r!(list!(
                predicted_classes = Robj::from(pred),
                set_sizes = Robj::from(set_sizes),
                average_set_size = r.average_set_size,
                coverage = r.coverage,
                score_quantile = r.score_quantile
            ))
        }
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Phase C1: CV-Conformal & Jackknife+ (closure-based conformal methods)
// =============================================================================

/// CV-conformal regression using functional linear model
#[extendr]
fn cv_conformal_fregre_lm_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    test_data: RMatrix<f64>,
    scalar_train: Nullable<RMatrix<f64>>,
    scalar_test: Nullable<RMatrix<f64>>,
    ncomp: i32,
    n_folds: i32,
    alpha: f64,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");

    let st = match scalar_train {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_train")) }
        Nullable::Null => None,
    };
    let ste = match scalar_test {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_test")) }
        Nullable::Null => None,
    };

    let nc = ncomp as usize;
    let fit_predict = |train_d: &FdMatrix, train_y: &[f64], train_sc: Option<&FdMatrix>,
                       pred_d: &FdMatrix, pred_sc: Option<&FdMatrix>|
     -> Option<(Vec<f64>, Vec<f64>)> {
        let model = scalar_on_function::fregre_lm(train_d, train_y, train_sc, nc).ok()?;
        let preds = scalar_on_function::predict_fregre_lm(&model, pred_d, pred_sc);
        Some((preds, vec![]))
    };

    match fdars_core::conformal::cv_conformal_regression(
        &fd, &y, &tfd, st.as_ref(), ste.as_ref(),
        fit_predict, n_folds as usize, alpha, seed as u64,
    ) {
        Ok(r) => r!(list!(
            predictions = Robj::from(r.predictions),
            lower = Robj::from(r.lower),
            upper = Robj::from(r.upper),
            residual_quantile = r.residual_quantile,
            coverage = r.coverage,
            calibration_scores = Robj::from(r.calibration_scores)
        )),
        Err(_) => r!(NULL),
    }
}

/// CV-conformal regression using nonparametric kernel regression
#[extendr]
fn cv_conformal_fregre_np_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    test_data: RMatrix<f64>,
    argvals: Vec<f64>,
    scalar_train: Nullable<RMatrix<f64>>,
    scalar_test: Nullable<RMatrix<f64>>,
    h_func: f64,
    h_scalar: f64,
    n_folds: i32,
    alpha: f64,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");

    let st = match scalar_train {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_train")) }
        Nullable::Null => None,
    };
    let ste = match scalar_test {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_test")) }
        Nullable::Null => None,
    };

    let av = argvals.clone();
    let hf = h_func;
    let hs = h_scalar;
    let fit_predict = move |train_d: &FdMatrix, train_y: &[f64], train_sc: Option<&FdMatrix>,
                            pred_d: &FdMatrix, pred_sc: Option<&FdMatrix>|
     -> Option<(Vec<f64>, Vec<f64>)> {
        let _ = scalar_on_function::fregre_np_mixed(train_d, train_y, &av, train_sc, hf, hs).ok()?;
        let preds = scalar_on_function::predict_fregre_np(
            train_d, train_y, train_sc, pred_d, pred_sc, &av, hf, hs,
        );
        Some((preds, vec![]))
    };

    match fdars_core::conformal::cv_conformal_regression(
        &fd, &y, &tfd, st.as_ref(), ste.as_ref(),
        fit_predict, n_folds as usize, alpha, seed as u64,
    ) {
        Ok(r) => r!(list!(
            predictions = Robj::from(r.predictions),
            lower = Robj::from(r.lower),
            upper = Robj::from(r.upper),
            residual_quantile = r.residual_quantile,
            coverage = r.coverage,
            calibration_scores = Robj::from(r.calibration_scores)
        )),
        Err(_) => r!(NULL),
    }
}

/// CV-conformal classification
#[extendr]
fn cv_conformal_classif_rust(
    data: RMatrix<f64>,
    y: Vec<i32>,
    test_data: RMatrix<f64>,
    covariates_train: Nullable<RMatrix<f64>>,
    covariates_test: Nullable<RMatrix<f64>>,
    ncomp: i32,
    classifier: &str,
    k_nn: i32,
    score_type: &str,
    n_folds: i32,
    alpha: f64,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");
    let y_usize: Vec<usize> = y.iter().map(|&c| c as usize).collect();

    let ct = match covariates_train {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid covariates_train")) }
        Nullable::Null => None,
    };
    let cte = match covariates_test {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid covariates_test")) }
        Nullable::Null => None,
    };

    let st = match score_type {
        "aps" => fdars_core::conformal::ClassificationScore::Aps,
        _ => fdars_core::conformal::ClassificationScore::Lac,
    };

    let nc = ncomp as usize;
    let kn = k_nn as usize;
    let clf = classifier.to_string();
    let fit_predict_probs = move |train_d: &FdMatrix, train_y: &[usize], train_sc: Option<&FdMatrix>,
                                   pred_d: &FdMatrix, _pred_sc: Option<&FdMatrix>|
     -> Option<(Vec<Vec<f64>>, Vec<Vec<f64>>)> {
        let fit = match clf.as_str() {
            "lda" => classification::fclassif_lda_fit(train_d, train_y, train_sc, nc).ok()?,
            "qda" => classification::fclassif_qda_fit(train_d, train_y, train_sc, nc).ok()?,
            "knn" => classification::fclassif_knn_fit(train_d, train_y, train_sc, nc, kn).ok()?,
            _ => classification::fclassif_lda_fit(train_d, train_y, train_sc, nc).ok()?,
        };
        // Project: center and multiply by rotation
        let n_pred = pred_d.nrows();
        let m = pred_d.ncols();
        let nc_fit = fit.ncomp;
        let mut scores_data = vec![0.0; n_pred * nc_fit];
        for i in 0..n_pred {
            for j in 0..nc_fit {
                let mut s = 0.0;
                for k in 0..m {
                    s += (pred_d[(i, k)] - fit.fpca_mean[k]) * fit.fpca_rotation[(k, j)];
                }
                scores_data[i * nc_fit + j] = s;
            }
        }
        let scores_mat = FdMatrix::from_slice(&scores_data, n_pred, nc_fit).unwrap();
        // Compute probabilities from ClassifMethod
        let probs = classif_probs_from_method(&fit.method, &scores_mat);
        Some((probs, vec![]))
    };

    match fdars_core::conformal::cv_conformal_classification(
        &fd, &y_usize, &tfd, ct.as_ref(), cte.as_ref(),
        fit_predict_probs, n_folds as usize, st, alpha, seed as u64,
    ) {
        Ok(r) => {
            let pred: Vec<i32> = r.predicted_classes.iter().map(|&c| c as i32).collect();
            let set_sizes: Vec<i32> = r.set_sizes.iter().map(|&s| s as i32).collect();
            r!(list!(
                predicted_classes = Robj::from(pred),
                set_sizes = Robj::from(set_sizes),
                average_set_size = r.average_set_size,
                coverage = r.coverage,
                score_quantile = r.score_quantile
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Jackknife+ regression using functional linear model
#[extendr]
fn jackknife_plus_fregre_lm_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    test_data: RMatrix<f64>,
    scalar_train: Nullable<RMatrix<f64>>,
    scalar_test: Nullable<RMatrix<f64>>,
    ncomp: i32,
    alpha: f64,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");

    let st = match scalar_train {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_train")) }
        Nullable::Null => None,
    };
    let ste = match scalar_test {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_test")) }
        Nullable::Null => None,
    };

    let nc = ncomp as usize;
    let fit_predict = |train_d: &FdMatrix, train_y: &[f64], train_sc: Option<&FdMatrix>,
                       pred_d: &FdMatrix, pred_sc: Option<&FdMatrix>|
     -> Option<(Vec<f64>, Vec<f64>)> {
        let model = scalar_on_function::fregre_lm(train_d, train_y, train_sc, nc).ok()?;
        let preds = scalar_on_function::predict_fregre_lm(&model, pred_d, pred_sc);
        Some((preds, vec![]))
    };

    match fdars_core::conformal::jackknife_plus_regression(
        &fd, &y, &tfd, st.as_ref(), ste.as_ref(),
        fit_predict, alpha,
    ) {
        Ok(r) => r!(list!(
            predictions = Robj::from(r.predictions),
            lower = Robj::from(r.lower),
            upper = Robj::from(r.upper),
            residual_quantile = r.residual_quantile,
            coverage = r.coverage,
            calibration_scores = Robj::from(r.calibration_scores)
        )),
        Err(_) => r!(NULL),
    }
}

/// Jackknife+ regression using nonparametric kernel regression
#[extendr]
fn jackknife_plus_fregre_np_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    test_data: RMatrix<f64>,
    argvals: Vec<f64>,
    scalar_train: Nullable<RMatrix<f64>>,
    scalar_test: Nullable<RMatrix<f64>>,
    h_func: f64,
    h_scalar: f64,
    alpha: f64,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");

    let st = match scalar_train {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_train")) }
        Nullable::Null => None,
    };
    let ste = match scalar_test {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_test")) }
        Nullable::Null => None,
    };

    let av = argvals.clone();
    let hf = h_func;
    let hs = h_scalar;
    let fit_predict = move |train_d: &FdMatrix, train_y: &[f64], train_sc: Option<&FdMatrix>,
                            pred_d: &FdMatrix, pred_sc: Option<&FdMatrix>|
     -> Option<(Vec<f64>, Vec<f64>)> {
        let _ = scalar_on_function::fregre_np_mixed(train_d, train_y, &av, train_sc, hf, hs).ok()?;
        let preds = scalar_on_function::predict_fregre_np(
            train_d, train_y, train_sc, pred_d, pred_sc, &av, hf, hs,
        );
        Some((preds, vec![]))
    };

    match fdars_core::conformal::jackknife_plus_regression(
        &fd, &y, &tfd, st.as_ref(), ste.as_ref(),
        fit_predict, alpha,
    ) {
        Ok(r) => r!(list!(
            predictions = Robj::from(r.predictions),
            lower = Robj::from(r.lower),
            upper = Robj::from(r.upper),
            residual_quantile = r.residual_quantile,
            coverage = r.coverage,
            calibration_scores = Robj::from(r.calibration_scores)
        )),
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Phase C2: Conformal Generic (trait-based)
// =============================================================================

/// Helper struct implementing FpcPredictor for regression models reconstructed from R
struct RFregreLmPredictor {
    mean: Vec<f64>,
    rotation: FdMatrix,
    scores: FdMatrix,
    nc: usize,
    coefficients: Vec<f64>,
    gamma: Vec<f64>,
}

impl fdars_core::explain_generic::FpcPredictor for RFregreLmPredictor {
    fn fpca_mean(&self) -> &[f64] { &self.mean }
    fn fpca_rotation(&self) -> &FdMatrix { &self.rotation }
    fn ncomp(&self) -> usize { self.nc }
    fn training_scores(&self) -> &FdMatrix { &self.scores }
    fn task_type(&self) -> fdars_core::explain_generic::TaskType {
        fdars_core::explain_generic::TaskType::Regression
    }
    fn predict_from_scores(&self, scores: &[f64], scalar_covariates: Option<&[f64]>) -> f64 {
        let mut yhat = self.coefficients[0]; // intercept
        for k in 0..self.nc {
            yhat += self.coefficients[1 + k] * scores[k];
        }
        if let Some(sc) = scalar_covariates {
            for j in 0..self.gamma.len() {
                yhat += self.gamma[j] * sc[j];
            }
        }
        yhat
    }
}

/// Helper struct implementing FpcPredictor for logistic models reconstructed from R
struct RLogisticPredictor {
    mean: Vec<f64>,
    rotation: FdMatrix,
    scores: FdMatrix,
    nc: usize,
    intercept: f64,
    coefficients: Vec<f64>,
    gamma: Vec<f64>,
}

impl fdars_core::explain_generic::FpcPredictor for RLogisticPredictor {
    fn fpca_mean(&self) -> &[f64] { &self.mean }
    fn fpca_rotation(&self) -> &FdMatrix { &self.rotation }
    fn ncomp(&self) -> usize { self.nc }
    fn training_scores(&self) -> &FdMatrix { &self.scores }
    fn task_type(&self) -> fdars_core::explain_generic::TaskType {
        fdars_core::explain_generic::TaskType::BinaryClassification
    }
    fn predict_from_scores(&self, scores: &[f64], scalar_covariates: Option<&[f64]>) -> f64 {
        let mut eta = self.intercept;
        for k in 0..self.nc {
            eta += self.coefficients[1 + k] * scores[k];
        }
        if let Some(sc) = scalar_covariates {
            for j in 0..self.gamma.len() {
                eta += self.gamma[j] * sc[j];
            }
        }
        1.0 / (1.0 + (-eta).exp()) // sigmoid
    }
}

/// Generic conformal regression for a pre-fitted fregre.lm model
#[extendr]
fn conformal_generic_regression_lm_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    test_data: RMatrix<f64>,
    scalar_train: Nullable<RMatrix<f64>>,
    scalar_test: Nullable<RMatrix<f64>>,
    fpca_mean: Vec<f64>,
    fpca_rotation: RMatrix<f64>,
    fpca_scores: RMatrix<f64>,
    ncomp: i32,
    coefficients: Vec<f64>,
    gamma: Vec<f64>,
    calibration_indices: Nullable<Vec<i32>>,
    cal_fraction: f64,
    alpha: f64,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");

    let st = match scalar_train {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_train")) }
        Nullable::Null => None,
    };
    let ste = match scalar_test {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_test")) }
        Nullable::Null => None,
    };

    let rot_s = fpca_rotation.as_real_slice().unwrap();
    let rot = FdMatrix::from_slice(rot_s, fpca_rotation.nrows(), fpca_rotation.ncols()).expect("Invalid rotation");
    let sc_s = fpca_scores.as_real_slice().unwrap();
    let sc = FdMatrix::from_slice(sc_s, fpca_scores.nrows(), fpca_scores.ncols()).expect("Invalid scores");

    let predictor = RFregreLmPredictor {
        mean: fpca_mean,
        rotation: rot,
        scores: sc,
        nc: ncomp as usize,
        coefficients,
        gamma,
    };

    let cal_idx: Option<Vec<usize>> = match calibration_indices {
        Nullable::NotNull(v) => Some(v.iter().map(|&i| i as usize).collect()),
        Nullable::Null => None,
    };

    match fdars_core::conformal::conformal_generic_regression(
        &predictor, &fd, &y, &tfd, st.as_ref(), ste.as_ref(),
        cal_idx.as_deref(), cal_fraction, alpha, seed as u64,
    ) {
        Ok(r) => r!(list!(
            predictions = Robj::from(r.predictions),
            lower = Robj::from(r.lower),
            upper = Robj::from(r.upper),
            residual_quantile = r.residual_quantile,
            coverage = r.coverage,
            calibration_scores = Robj::from(r.calibration_scores)
        )),
        Err(_) => r!(NULL),
    }
}

/// Generic conformal classification for a pre-fitted logistic model
#[extendr]
fn conformal_generic_classification_logistic_rust(
    data: RMatrix<f64>,
    y: Vec<i32>,
    test_data: RMatrix<f64>,
    scalar_train: Nullable<RMatrix<f64>>,
    scalar_test: Nullable<RMatrix<f64>>,
    fpca_mean: Vec<f64>,
    fpca_rotation: RMatrix<f64>,
    fpca_scores: RMatrix<f64>,
    ncomp: i32,
    intercept: f64,
    coefficients: Vec<f64>,
    gamma: Vec<f64>,
    score_type: &str,
    calibration_indices: Nullable<Vec<i32>>,
    cal_fraction: f64,
    alpha: f64,
    seed: i32,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let tds = test_data.as_real_slice().unwrap();
    let tfd = FdMatrix::from_slice(tds, test_data.nrows(), test_data.ncols()).expect("Invalid test_data");
    let y_usize: Vec<usize> = y.iter().map(|&c| c as usize).collect();

    let st = match scalar_train {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_train")) }
        Nullable::Null => None,
    };
    let ste = match scalar_test {
        Nullable::NotNull(s) => { let ss = s.as_real_slice().unwrap(); Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid scalar_test")) }
        Nullable::Null => None,
    };

    let rot_s = fpca_rotation.as_real_slice().unwrap();
    let rot = FdMatrix::from_slice(rot_s, fpca_rotation.nrows(), fpca_rotation.ncols()).expect("Invalid rotation");
    let sc_s = fpca_scores.as_real_slice().unwrap();
    let sc = FdMatrix::from_slice(sc_s, fpca_scores.nrows(), fpca_scores.ncols()).expect("Invalid scores");

    let score = match score_type {
        "aps" => fdars_core::conformal::ClassificationScore::Aps,
        _ => fdars_core::conformal::ClassificationScore::Lac,
    };

    let predictor = RLogisticPredictor {
        mean: fpca_mean,
        rotation: rot,
        scores: sc,
        nc: ncomp as usize,
        intercept,
        coefficients,
        gamma,
    };

    let cal_idx: Option<Vec<usize>> = match calibration_indices {
        Nullable::NotNull(v) => Some(v.iter().map(|&i| i as usize).collect()),
        Nullable::Null => None,
    };

    match fdars_core::conformal::conformal_generic_classification(
        &predictor, &fd, &y_usize, &tfd, st.as_ref(), ste.as_ref(),
        score, cal_idx.as_deref(), cal_fraction, alpha, seed as u64,
    ) {
        Ok(r) => {
            let pred: Vec<i32> = r.predicted_classes.iter().map(|&c| c as i32).collect();
            let set_sizes: Vec<i32> = r.set_sizes.iter().map(|&s| s as i32).collect();
            r!(list!(
                predicted_classes = Robj::from(pred),
                set_sizes = Robj::from(set_sizes),
                average_set_size = r.average_set_size,
                coverage = r.coverage,
                score_quantile = r.score_quantile
            ))
        }
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Phase C3: Outliers threshold with distribution
// =============================================================================

/// LRT outlier threshold with full bootstrap distribution
#[extendr]
fn outliers_thres_lrt_with_dist_rust(
    data: RMatrix<f64>,
    _argvals: Vec<f64>,
    nb: i32,
    smo: f64,
    trim: f64,
    seed: f64,
    percentile: f64,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    let (threshold, boot_dist) = fdars_core::outliers::outliers_threshold_lrt_with_dist(
        &fd, nb as usize, smo, trim, seed as u64, percentile,
    );

    r!(list!(
        threshold = threshold,
        boot_distribution = Robj::from(boot_dist)
    ))
}

// =============================================================================
// Phase C4: Gradient functions
// =============================================================================

/// High-accuracy gradient using 5-point stencil (uniform) or 3-point Lagrange (non-uniform)
#[extendr]
fn fdata_gradient_rust(data: RMatrix<f64>, argvals: Vec<f64>) -> Robj {
    let nrow = data.nrows();
    let ncol = data.ncols();
    let ds = data.as_real_slice().unwrap();

    // Data arrives column-major from R; process each row
    let mut result = vec![0.0; nrow * ncol];
    for i in 0..nrow {
        let row: Vec<f64> = (0..ncol).map(|j| ds[j * nrow + i]).collect();
        let grad = fdars_core::helpers::gradient(&row, &argvals);
        for j in 0..ncol {
            result[j * nrow + i] = grad[j];
        }
    }

    let result_mat = RMatrix::new_matrix(nrow, ncol, |r, c| result[c * nrow + r]);
    r!(result_mat)
}

// =============================================================================
// Robust Regression
// =============================================================================

/// L1 (LAD) functional regression
#[extendr]
fn fregre_l1_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    scalar_covariates: Nullable<RMatrix<f64>>,
    ncomp: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data dimensions");

    let scov = match scalar_covariates {
        Nullable::NotNull(sc) => {
            let sc_s = sc.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(sc_s, sc.nrows(), sc.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    match scalar_on_function::fregre_l1(&fd, &y, scov.as_ref(), ncomp as usize) {
        Ok(r) => {
            let beta_t = Robj::from(r.beta_t.clone());
            let fitted = Robj::from(r.fitted_values.clone());
            let resid = Robj::from(r.residuals.clone());
            let coefs = Robj::from(r.coefficients.clone());
            let wts = Robj::from(r.weights.clone());
            r!(list!(
                intercept = r.intercept,
                beta_t = beta_t,
                fitted_values = fitted,
                residuals = resid,
                coefficients = coefs,
                ncomp = r.ncomp as i32,
                iterations = r.iterations as i32,
                converged = r.converged,
                weights = wts,
                r_squared = r.r_squared,
                fpca_mean = Robj::from(r.fpca.mean.clone()),
                fpca_rotation_data = Robj::from(r.fpca.rotation.as_slice().to_vec()),
                fpca_rotation_nrow = r.fpca.rotation.nrows() as i32,
                fpca_rotation_ncol = r.fpca.rotation.ncols() as i32,
                fpca_scores_data = Robj::from(r.fpca.scores.as_slice().to_vec()),
                fpca_scores_nrow = r.fpca.scores.nrows() as i32,
                fpca_scores_ncol = r.fpca.scores.ncols() as i32
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Huber M-estimation functional regression
#[extendr]
fn fregre_huber_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    scalar_covariates: Nullable<RMatrix<f64>>,
    ncomp: i32,
    k: f64,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data dimensions");

    let scov = match scalar_covariates {
        Nullable::NotNull(sc) => {
            let sc_s = sc.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(sc_s, sc.nrows(), sc.ncols()).expect("Invalid covariates"))
        }
        Nullable::Null => None,
    };

    match scalar_on_function::fregre_huber(&fd, &y, scov.as_ref(), ncomp as usize, k) {
        Ok(r) => {
            let beta_t = Robj::from(r.beta_t.clone());
            let fitted = Robj::from(r.fitted_values.clone());
            let resid = Robj::from(r.residuals.clone());
            let coefs = Robj::from(r.coefficients.clone());
            let wts = Robj::from(r.weights.clone());
            r!(list!(
                intercept = r.intercept,
                beta_t = beta_t,
                fitted_values = fitted,
                residuals = resid,
                coefficients = coefs,
                ncomp = r.ncomp as i32,
                iterations = r.iterations as i32,
                converged = r.converged,
                weights = wts,
                r_squared = r.r_squared,
                fpca_mean = Robj::from(r.fpca.mean.clone()),
                fpca_rotation_data = Robj::from(r.fpca.rotation.as_slice().to_vec()),
                fpca_rotation_nrow = r.fpca.rotation.nrows() as i32,
                fpca_rotation_ncol = r.fpca.rotation.ncols() as i32,
                fpca_scores_data = Robj::from(r.fpca.scores.as_slice().to_vec()),
                fpca_scores_nrow = r.fpca.scores.nrows() as i32,
                fpca_scores_ncol = r.fpca.scores.ncols() as i32
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Predict from robust (L1/Huber) regression
#[extendr]
fn predict_fregre_robust_rust(
    intercept: f64,
    coefficients: Vec<f64>,
    fpca_mean: Vec<f64>,
    fpca_rotation_data: Vec<f64>,
    fpca_rotation_nrow: i32,
    fpca_rotation_ncol: i32,
    ncomp: i32,
    new_data: RMatrix<f64>,
    new_scalar: Nullable<RMatrix<f64>>,
) -> Robj {
    let rot_nrow = fpca_rotation_nrow as usize;
    let rot_ncol = fpca_rotation_ncol as usize;
    let rot = FdMatrix::from_slice(&fpca_rotation_data, rot_nrow, rot_ncol)
        .expect("Invalid rotation dimensions");
    let nds = new_data.as_real_slice().unwrap();
    let n_new = new_data.nrows();
    let m = new_data.ncols();
    let nfd = FdMatrix::from_slice(nds, n_new, m).expect("Invalid new_data dimensions");
    let nc = ncomp as usize;

    let ns = match new_scalar {
        Nullable::NotNull(s) => {
            let ss = s.as_real_slice().unwrap();
            Some(FdMatrix::from_slice(ss, s.nrows(), s.ncols()).expect("Invalid new_scalar"))
        }
        Nullable::Null => None,
    };

    // Project new data onto FPCA space and predict
    let p_scalar = coefficients.len() - 1 - nc;
    let mut predictions = vec![0.0; n_new];
    for i in 0..n_new {
        let mut yhat = intercept;
        for k in 0..nc {
            let mut s = 0.0;
            for j in 0..m {
                s += (nfd[(i, j)] - fpca_mean[j]) * rot[(j, k)];
            }
            yhat += coefficients[1 + k] * s;
        }
        if let Some(ref sc) = ns {
            for j in 0..p_scalar {
                yhat += coefficients[1 + nc + j] * sc[(i, j)];
            }
        }
        predictions[i] = yhat;
    }

    Robj::from(predictions)
}

// =============================================================================
// Semimetric functions (PCA, deriv, basis)
// =============================================================================

/// PCA-based semimetric (self-distances)
#[extendr]
fn semimetric_pca_self_1d(data: RMatrix<f64>, ncomp: i32) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    match core_metric::pca_self_1d(&fd, ncomp as usize) {
        Ok(result) => fdmatrix_to_robj(&result),
        Err(_) => r!(NULL),
    }
}

/// PCA-based semimetric (cross-distances)
#[extendr]
fn semimetric_pca_cross_1d(data1: RMatrix<f64>, data2: RMatrix<f64>, ncomp: i32) -> Robj {
    let ds1 = data1.as_real_slice().unwrap();
    let fd1 = FdMatrix::from_slice(ds1, data1.nrows(), data1.ncols()).expect("Invalid data1");
    let ds2 = data2.as_real_slice().unwrap();
    let fd2 = FdMatrix::from_slice(ds2, data2.nrows(), data2.ncols()).expect("Invalid data2");
    match core_metric::pca_cross_1d(&fd1, &fd2, ncomp as usize) {
        Ok(result) => fdmatrix_to_robj(&result),
        Err(_) => r!(NULL),
    }
}

/// Derivative-based semimetric (self-distances)
#[extendr]
fn semimetric_deriv_self_1d(data: RMatrix<f64>, argvals: Vec<f64>, nderiv: i32, weights: Vec<f64>) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let result = core_metric::deriv_self_1d(&fd, &argvals, nderiv as usize, &weights);
    fdmatrix_to_robj(&result)
}

/// Derivative-based semimetric (cross-distances)
#[extendr]
fn semimetric_deriv_cross_1d(data1: RMatrix<f64>, data2: RMatrix<f64>, argvals: Vec<f64>, nderiv: i32, weights: Vec<f64>) -> Robj {
    let ds1 = data1.as_real_slice().unwrap();
    let fd1 = FdMatrix::from_slice(ds1, data1.nrows(), data1.ncols()).expect("Invalid data1");
    let ds2 = data2.as_real_slice().unwrap();
    let fd2 = FdMatrix::from_slice(ds2, data2.nrows(), data2.ncols()).expect("Invalid data2");
    let result = core_metric::deriv_cross_1d(&fd1, &fd2, &argvals, nderiv as usize, &weights);
    fdmatrix_to_robj(&result)
}

/// Basis coefficient semimetric (self-distances)
#[extendr]
fn semimetric_basis_self_1d(data: RMatrix<f64>, argvals: Vec<f64>, nbasis: i32, basis_type: i32) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");
    let bt = fdars_core::basis::ProjectionBasisType::from_i32(basis_type);
    let result = core_metric::basis_coef_self_1d(&fd, &argvals, nbasis as usize, bt);
    fdmatrix_to_robj(&result)
}

/// Basis coefficient semimetric (cross-distances)
#[extendr]
fn semimetric_basis_cross_1d(data1: RMatrix<f64>, data2: RMatrix<f64>, argvals: Vec<f64>, nbasis: i32, basis_type: i32) -> Robj {
    let ds1 = data1.as_real_slice().unwrap();
    let fd1 = FdMatrix::from_slice(ds1, data1.nrows(), data1.ncols()).expect("Invalid data1");
    let ds2 = data2.as_real_slice().unwrap();
    let fd2 = FdMatrix::from_slice(ds2, data2.nrows(), data2.ncols()).expect("Invalid data2");
    let bt = fdars_core::basis::ProjectionBasisType::from_i32(basis_type);
    let result = core_metric::basis_coef_cross_1d(&fd1, &fd2, &argvals, nbasis as usize, bt);
    fdmatrix_to_robj(&result)
}

// =============================================================================
// Scalar-on-Shape Regression
// =============================================================================

/// Fit scalar-on-shape regression model
#[extendr]
fn scalar_on_shape_rust(
    data: RMatrix<f64>,
    y: Vec<f64>,
    argvals: Vec<f64>,
    nbasis: i32,
    lambda: f64,
    index_method: &str,
    index_param: f64,
    g_degree: i32,
    dp_lambda: f64,
    max_iter_outer: i32,
    max_iter_inner: i32,
    tol: f64,
) -> Robj {
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, data.nrows(), data.ncols()).expect("Invalid data");

    let idx = match index_method {
        "polynomial" => fdars_core::elastic_regression::IndexMethod::Polynomial(index_param as usize),
        "nadaraya.watson" => fdars_core::elastic_regression::IndexMethod::NadarayaWatson(index_param),
        _ => fdars_core::elastic_regression::IndexMethod::Identity,
    };

    let config = fdars_core::elastic_regression::ScalarOnShapeConfig {
        nbasis: nbasis as usize,
        lambda,
        lfd_order: 2,
        index_method: idx,
        g_degree: g_degree as usize,
        max_iter_outer: max_iter_outer as usize,
        max_iter_inner: max_iter_inner as usize,
        tol,
        dp_lambda,
    };

    match fdars_core::elastic_regression::scalar_on_shape(&fd, &y, &argvals, &config) {
        Ok(r) => {
            let gammas = fdmatrix_to_robj(&r.gammas);
            let idx_str = match &r.index_method {
                fdars_core::elastic_regression::IndexMethod::Identity => "identity",
                fdars_core::elastic_regression::IndexMethod::Polynomial(_) => "polynomial",
                fdars_core::elastic_regression::IndexMethod::NadarayaWatson(_) => "nadaraya.watson",
                _ => "unknown",
            };
            let idx_param = match &r.index_method {
                fdars_core::elastic_regression::IndexMethod::Polynomial(d) => *d as f64,
                fdars_core::elastic_regression::IndexMethod::NadarayaWatson(bw) => *bw,
                _ => 0.0,
            };
            r!(list!(
                beta = Robj::from(r.beta),
                beta_coefficients = Robj::from(r.beta_coefficients),
                gammas = gammas,
                shape_scores = Robj::from(r.shape_scores),
                h_coefficients = Robj::from(r.h_coefficients),
                g_coefficients = Robj::from(r.g_coefficients),
                fitted_values = Robj::from(r.fitted_values),
                residuals = Robj::from(r.residuals),
                sse = r.sse,
                r_squared = r.r_squared,
                n_iter_outer = r.n_iter_outer as i32,
                n_iter_inner = r.n_iter_inner as i32,
                index_method = idx_str,
                index_param = idx_param
            ))
        }
        Err(_) => r!(NULL),
    }
}

/// Predict from a scalar-on-shape model
#[extendr]
fn predict_scalar_on_shape_rust(
    beta: Vec<f64>,
    beta_coefficients: Vec<f64>,
    gammas: RMatrix<f64>,
    shape_scores: Vec<f64>,
    h_coefficients: Vec<f64>,
    g_coefficients: Vec<f64>,
    fitted_values: Vec<f64>,
    residuals: Vec<f64>,
    sse: f64,
    r_squared: f64,
    n_iter_outer: i32,
    n_iter_inner: i32,
    index_method: &str,
    index_param: f64,
    new_data: RMatrix<f64>,
    argvals: Vec<f64>,
) -> Robj {
    let gs = gammas.as_real_slice().unwrap();
    let gm = FdMatrix::from_slice(gs, gammas.nrows(), gammas.ncols()).expect("Invalid gammas");

    let idx = match index_method {
        "polynomial" => fdars_core::elastic_regression::IndexMethod::Polynomial(index_param as usize),
        "nadaraya.watson" => fdars_core::elastic_regression::IndexMethod::NadarayaWatson(index_param),
        _ => fdars_core::elastic_regression::IndexMethod::Identity,
    };

    let fit = fdars_core::elastic_regression::ScalarOnShapeResult::from_parts(
        beta,
        beta_coefficients,
        gm,
        shape_scores,
        h_coefficients,
        g_coefficients,
        fitted_values,
        residuals,
        sse,
        r_squared,
        n_iter_outer as usize,
        n_iter_inner as usize,
        idx,
    );

    let nds = new_data.as_real_slice().unwrap();
    let nfd = FdMatrix::from_slice(nds, new_data.nrows(), new_data.ncols()).expect("Invalid new_data");

    match fit.predict(&nfd, &argvals) {
        Ok(preds) => r!(Robj::from(preds)),
        Err(_) => r!(NULL),
    }
}

// =============================================================================
// Statistical Process Monitoring (SPM) wrappers
// =============================================================================

/// SPM Phase I: Build a univariate control chart from in-control functional data.
#[extendr]
fn spm_phase1_rust(data: RMatrix<f64>, argvals: Vec<f64>,
                   ncomp: i32, alpha: f64, tuning_fraction: f64, seed: i32) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");

    let config = fdars_core::spm::SpmConfig {
        ncomp: ncomp as usize,
        alpha,
        tuning_fraction,
        seed: seed as u64,
    };

    match fdars_core::spm::spm_phase1(&fd, &argvals, &config) {
        Ok(chart) => {
            r!(list!(
                rotation = fdmatrix_to_robj(&chart.fpca.rotation),
                scores = fdmatrix_to_robj(&chart.fpca.scores),
                mean = Robj::from(chart.fpca.mean.clone()),
                singular_values = Robj::from(chart.fpca.singular_values.clone()),
                centered = fdmatrix_to_robj(&chart.fpca.centered),
                eigenvalues = Robj::from(chart.eigenvalues.clone()),
                t2_phase1 = Robj::from(chart.t2_phase1.clone()),
                spe_phase1 = Robj::from(chart.spe_phase1.clone()),
                t2_ucl = chart.t2_limit.ucl,
                t2_alpha = chart.t2_limit.alpha,
                t2_description = chart.t2_limit.description.clone(),
                spe_ucl = chart.spe_limit.ucl,
                spe_alpha = chart.spe_limit.alpha,
                spe_description = chart.spe_limit.description.clone(),
                ncomp = chart.eigenvalues.len() as i32,
                config_ncomp = config.ncomp as i32,
                config_alpha = config.alpha,
                config_tuning_fraction = config.tuning_fraction,
                config_seed = config.seed as i32
            ))
        }
        Err(e) => {
            rprintln!("spm_phase1 failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// SPM Phase II: Monitor new data against an established chart.
#[extendr]
fn spm_monitor_rust(
    new_data: RMatrix<f64>, argvals: Vec<f64>,
    rotation: RMatrix<f64>, mean: Vec<f64>, singular_values: Vec<f64>,
    centered: RMatrix<f64>,
    eigenvalues: Vec<f64>,
    t2_ucl: f64, t2_alpha: f64, t2_description: &str,
    spe_ucl: f64, spe_alpha: f64, spe_description: &str,
    ncomp: i32, config_alpha: f64, config_tuning_fraction: f64, config_seed: i32,
) -> Robj {
    use fdars_core::regression::FpcaResult;
    use fdars_core::spm::phase::SpmChart;
    use fdars_core::spm::ControlLimit;

    let n = new_data.nrows();
    let m = new_data.ncols();
    let nds = new_data.as_real_slice().unwrap();
    let nd = FdMatrix::from_slice(nds, n, m).expect("Invalid new_data");

    // Reconstruct FpcaResult
    let rs = rotation.as_real_slice().unwrap();
    let rot = FdMatrix::from_slice(rs, rotation.nrows(), rotation.ncols()).expect("Invalid rotation");

    let cs = centered.as_real_slice().unwrap();
    let cen = FdMatrix::from_slice(cs, centered.nrows(), centered.ncols()).expect("Invalid centered");

    // Build scores placeholder (not needed for projection, but required for struct)
    let n_tune = centered.nrows();
    let nc = ncomp as usize;
    let scores_placeholder = FdMatrix::zeros(n_tune, nc);

    let fpca = FpcaResult {
        singular_values,
        rotation: rot,
        scores: scores_placeholder,
        mean,
        centered: cen,
    };

    let t2_limit = ControlLimit::from_parts(t2_ucl, t2_alpha, t2_description.to_string());
    let spe_limit = ControlLimit::from_parts(spe_ucl, spe_alpha, spe_description.to_string());
    let config = fdars_core::spm::SpmConfig {
        ncomp: nc,
        alpha: config_alpha,
        tuning_fraction: config_tuning_fraction,
        seed: config_seed as u64,
    };

    let chart = SpmChart::from_parts(
        fpca, eigenvalues, vec![], vec![], t2_limit, spe_limit, config, true,
    );

    match fdars_core::spm::spm_monitor(&chart, &nd, &argvals) {
        Ok(result) => {
            let t2_alarm_int: Vec<i32> = result.t2_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            let spe_alarm_int: Vec<i32> = result.spe_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            r!(list!(
                t2 = Robj::from(result.t2),
                spe = Robj::from(result.spe),
                t2_alarm = Robj::from(t2_alarm_int),
                spe_alarm = Robj::from(spe_alarm_int),
                scores = fdmatrix_to_robj(&result.scores)
            ))
        }
        Err(e) => {
            rprintln!("spm_monitor failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Multivariate FPCA.
#[extendr]
fn mfpca_rust(data_list: List, ncomp: i32, weighted: bool) -> Robj {
    // data_list is a list of matrices from R
    let n_vars = data_list.len();
    if n_vars == 0 {
        rprintln!("mfpca: empty data_list");
        return r!(NULL);
    }

    let mut matrices: Vec<FdMatrix> = Vec::with_capacity(n_vars);
    for (_, item) in data_list.iter() {
        let rmat: RMatrix<f64> = item.try_into().expect("data_list elements must be matrices");
        let s = rmat.as_real_slice().unwrap();
        let fd = FdMatrix::from_slice(s, rmat.nrows(), rmat.ncols()).expect("Invalid matrix in data_list");
        matrices.push(fd);
    }

    let refs: Vec<&FdMatrix> = matrices.iter().collect();
    let config = fdars_core::spm::MfpcaConfig {
        ncomp: ncomp as usize,
        weighted,
    };

    match fdars_core::spm::mfpca(&refs, &config) {
        Ok(result) => {
            // Convert eigenfunctions to a list of matrices
            let ef_list = List::from_values(
                result.eigenfunctions.iter().map(|ef| fdmatrix_to_robj(ef))
            );

            // Convert means to a list of vectors
            let means_list = List::from_values(
                result.means.iter().map(|m| Robj::from(m.clone()))
            );

            let grid_sizes_i32: Vec<i32> = result.grid_sizes.iter().map(|&g| g as i32).collect();

            r!(list!(
                scores = fdmatrix_to_robj(&result.scores),
                eigenfunctions = Robj::from(ef_list),
                eigenvalues = Robj::from(result.eigenvalues.clone()),
                means = Robj::from(means_list),
                scales = Robj::from(result.scales.clone()),
                grid_sizes = Robj::from(grid_sizes_i32)
            ))
        }
        Err(e) => {
            rprintln!("mfpca failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// EWMA-based SPM monitoring.
#[extendr]
fn spm_ewma_rust(
    sequential_data: RMatrix<f64>, argvals: Vec<f64>,
    // Chart fields (same as spm_monitor_rust):
    rotation: RMatrix<f64>, mean: Vec<f64>, singular_values: Vec<f64>,
    centered: RMatrix<f64>,
    eigenvalues: Vec<f64>,
    t2_ucl: f64, t2_alpha: f64, t2_description: &str,
    spe_ucl: f64, spe_alpha: f64, spe_description: &str,
    chart_ncomp: i32, config_alpha: f64, config_tuning_fraction: f64, config_seed: i32,
    // EWMA config:
    lambda: f64, ewma_ncomp: i32, ewma_alpha: f64,
) -> Robj {
    use fdars_core::regression::FpcaResult;
    use fdars_core::spm::phase::SpmChart;
    use fdars_core::spm::ControlLimit;

    let n = sequential_data.nrows();
    let m = sequential_data.ncols();
    let sds = sequential_data.as_real_slice().unwrap();
    let sd = FdMatrix::from_slice(sds, n, m).expect("Invalid sequential_data");

    let rs = rotation.as_real_slice().unwrap();
    let rot = FdMatrix::from_slice(rs, rotation.nrows(), rotation.ncols()).expect("Invalid rotation");
    let cs = centered.as_real_slice().unwrap();
    let cen = FdMatrix::from_slice(cs, centered.nrows(), centered.ncols()).expect("Invalid centered");

    let n_tune = centered.nrows();
    let nc = chart_ncomp as usize;
    let scores_placeholder = FdMatrix::zeros(n_tune, nc);

    let fpca = FpcaResult {
        singular_values,
        rotation: rot,
        scores: scores_placeholder,
        mean,
        centered: cen,
    };

    let t2_limit = ControlLimit::from_parts(t2_ucl, t2_alpha, t2_description.to_string());
    let spe_limit = ControlLimit::from_parts(spe_ucl, spe_alpha, spe_description.to_string());
    let config = fdars_core::spm::SpmConfig {
        ncomp: nc,
        alpha: config_alpha,
        tuning_fraction: config_tuning_fraction,
        seed: config_seed as u64,
    };

    let chart = SpmChart::from_parts(
        fpca, eigenvalues, vec![], vec![], t2_limit, spe_limit, config, true,
    );

    let ewma_config = fdars_core::spm::EwmaConfig {
        lambda,
        ncomp: ewma_ncomp as usize,
        alpha: ewma_alpha,
        exact_covariance: false,
    };

    match fdars_core::spm::spm_ewma_monitor(&chart, &sd, &argvals, &ewma_config) {
        Ok(result) => {
            let t2_alarm_int: Vec<i32> = result.t2_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            let spe_alarm_int: Vec<i32> = result.spe_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            r!(list!(
                smoothed_scores = fdmatrix_to_robj(&result.smoothed_scores),
                t2 = Robj::from(result.t2),
                spe = Robj::from(result.spe),
                t2_limit = result.t2_limit,
                spe_limit = result.spe_limit,
                t2_alarm = Robj::from(t2_alarm_int),
                spe_alarm = Robj::from(spe_alarm_int)
            ))
        }
        Err(e) => {
            rprintln!("spm_ewma failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// T-squared contribution diagnostics.
#[extendr]
fn spm_t2_contrib_rust(scores: RMatrix<f64>, eigenvalues: Vec<f64>, grid_sizes: Vec<i32>) -> Robj {
    let n = scores.nrows();
    let nc = scores.ncols();
    let ss = scores.as_real_slice().unwrap();
    let sc = FdMatrix::from_slice(ss, n, nc).expect("Invalid scores");

    let gs: Vec<usize> = grid_sizes.iter().map(|&g| g as usize).collect();

    match fdars_core::spm::t2_contributions(&sc, &eigenvalues, &gs) {
        Ok(contrib) => fdmatrix_to_robj(&contrib),
        Err(e) => {
            rprintln!("t2_contributions failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// SPE contribution diagnostics.
#[extendr]
fn spm_spe_contrib_rust(
    std_list: List, recon_list: List, argvals_list: List,
) -> Robj {
    let p = std_list.len();
    if p == 0 {
        rprintln!("spe_contributions: empty input");
        return r!(NULL);
    }

    let mut std_mats: Vec<FdMatrix> = Vec::with_capacity(p);
    let mut recon_mats: Vec<FdMatrix> = Vec::with_capacity(p);
    let mut argvals_vecs: Vec<Vec<f64>> = Vec::with_capacity(p);

    for (_, item) in std_list.iter() {
        let rmat: RMatrix<f64> = item.try_into().expect("std_list elements must be matrices");
        let s = rmat.as_real_slice().unwrap();
        let fd = FdMatrix::from_slice(s, rmat.nrows(), rmat.ncols()).expect("Invalid matrix");
        std_mats.push(fd);
    }

    for (_, item) in recon_list.iter() {
        let rmat: RMatrix<f64> = item.try_into().expect("recon_list elements must be matrices");
        let s = rmat.as_real_slice().unwrap();
        let fd = FdMatrix::from_slice(s, rmat.nrows(), rmat.ncols()).expect("Invalid matrix");
        recon_mats.push(fd);
    }

    for (_, item) in argvals_list.iter() {
        let v: Vec<f64> = item.as_real_vector().expect("argvals_list elements must be numeric vectors");
        argvals_vecs.push(v);
    }

    let std_refs: Vec<&FdMatrix> = std_mats.iter().collect();
    let recon_refs: Vec<&FdMatrix> = recon_mats.iter().collect();
    let argvals_slices: Vec<&[f64]> = argvals_vecs.iter().map(|v| v.as_slice()).collect();

    match fdars_core::spm::spe_contributions(&std_refs, &recon_refs, &argvals_slices) {
        Ok(contrib) => fdmatrix_to_robj(&contrib),
        Err(e) => {
            rprintln!("spe_contributions failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// FRCC Phase I: Build a functional regression control chart.
#[extendr]
fn frcc_phase1_rust(
    y_curves: RMatrix<f64>, predictors: RMatrix<f64>, argvals: Vec<f64>,
    ncomp: i32, fosr_lambda: f64, alpha: f64, tuning_fraction: f64, seed: i32,
) -> Robj {
    let n = y_curves.nrows();
    let m = y_curves.ncols();
    let ys = y_curves.as_real_slice().unwrap();
    let y = FdMatrix::from_slice(ys, n, m).expect("Invalid y_curves");

    let ps = predictors.as_real_slice().unwrap();
    let pred = FdMatrix::from_slice(ps, predictors.nrows(), predictors.ncols()).expect("Invalid predictors");

    let config = fdars_core::spm::FrccConfig {
        ncomp: ncomp as usize,
        fosr_lambda,
        alpha,
        tuning_fraction,
        seed: seed as u64,
        min_r_squared: 0.0,
    };

    match fdars_core::spm::frcc_phase1(&y, &pred, &argvals, &config) {
        Ok(chart) => {
            r!(list!(
                // FOSR fields
                fosr_intercept = Robj::from(chart.fosr.intercept.clone()),
                fosr_beta = fdmatrix_to_robj(&chart.fosr.beta),
                fosr_fitted = fdmatrix_to_robj(&chart.fosr.fitted),
                fosr_residuals = fdmatrix_to_robj(&chart.fosr.residuals),
                fosr_r_squared_t = Robj::from(chart.fosr.r_squared_t.clone()),
                fosr_r_squared = chart.fosr.r_squared,
                fosr_beta_se = fdmatrix_to_robj(&chart.fosr.beta_se),
                fosr_lambda = chart.fosr.lambda,
                fosr_gcv = chart.fosr.gcv,
                // Residual FPCA fields
                resid_rotation = fdmatrix_to_robj(&chart.residual_fpca.rotation),
                resid_scores = fdmatrix_to_robj(&chart.residual_fpca.scores),
                resid_mean = Robj::from(chart.residual_fpca.mean.clone()),
                resid_singular_values = Robj::from(chart.residual_fpca.singular_values.clone()),
                resid_centered = fdmatrix_to_robj(&chart.residual_fpca.centered),
                // Eigenvalues
                eigenvalues = Robj::from(chart.eigenvalues.clone()),
                // Control limits
                t2_ucl = chart.t2_limit.ucl,
                t2_alpha = chart.t2_limit.alpha,
                t2_description = chart.t2_limit.description.clone(),
                spe_ucl = chart.spe_limit.ucl,
                spe_alpha = chart.spe_limit.alpha,
                spe_description = chart.spe_limit.description.clone(),
                // Config
                ncomp = chart.eigenvalues.len() as i32,
                config_ncomp = config.ncomp as i32,
                config_fosr_lambda = config.fosr_lambda,
                config_alpha = config.alpha,
                config_tuning_fraction = config.tuning_fraction,
                config_seed = config.seed as i32
            ))
        }
        Err(e) => {
            rprintln!("frcc_phase1 failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// FRCC Phase II: Monitor new data against a functional regression control chart.
#[extendr]
fn frcc_monitor_rust(
    new_y: RMatrix<f64>, new_predictors: RMatrix<f64>, argvals: Vec<f64>,
    // FOSR fields
    fosr_intercept: Vec<f64>, fosr_beta: RMatrix<f64>,
    fosr_lambda: f64,
    // Residual FPCA fields
    resid_rotation: RMatrix<f64>, resid_mean: Vec<f64>,
    resid_singular_values: Vec<f64>, resid_centered: RMatrix<f64>,
    // Eigenvalues & limits
    eigenvalues: Vec<f64>,
    t2_ucl: f64, t2_alpha: f64, t2_description: &str,
    spe_ucl: f64, spe_alpha: f64, spe_description: &str,
    ncomp: i32, config_fosr_lambda: f64, config_alpha: f64,
    config_tuning_fraction: f64, config_seed: i32,
) -> Robj {
    use fdars_core::regression::FpcaResult;
    use fdars_core::spm::ControlLimit;

    let n = new_y.nrows();
    let m = new_y.ncols();
    let nys = new_y.as_real_slice().unwrap();
    let ny = FdMatrix::from_slice(nys, n, m).expect("Invalid new_y");

    let nps = new_predictors.as_real_slice().unwrap();
    let np = FdMatrix::from_slice(nps, new_predictors.nrows(), new_predictors.ncols()).expect("Invalid new_predictors");

    // Reconstruct FosrResult (only intercept and beta needed for predict_fosr)
    let fbs = fosr_beta.as_real_slice().unwrap();
    let fb = FdMatrix::from_slice(fbs, fosr_beta.nrows(), fosr_beta.ncols()).expect("Invalid fosr_beta");

    let fosr = fdars_core::function_on_scalar::FosrResult::from_parts(
        fosr_intercept,
        fb,
        FdMatrix::zeros(1, 1),
        FdMatrix::zeros(1, 1),
        vec![],
        0.0,
        FdMatrix::zeros(1, 1),
        fosr_lambda,
        0.0,
    );

    // Reconstruct residual FPCA
    let rrs = resid_rotation.as_real_slice().unwrap();
    let rr = FdMatrix::from_slice(rrs, resid_rotation.nrows(), resid_rotation.ncols()).expect("Invalid resid_rotation");
    let rcs = resid_centered.as_real_slice().unwrap();
    let rc = FdMatrix::from_slice(rcs, resid_centered.nrows(), resid_centered.ncols()).expect("Invalid resid_centered");

    let nc = ncomp as usize;
    let n_cal = resid_centered.nrows();
    let resid_scores_placeholder = FdMatrix::zeros(n_cal, nc);

    let residual_fpca = FpcaResult {
        singular_values: resid_singular_values,
        rotation: rr,
        scores: resid_scores_placeholder,
        mean: resid_mean,
        centered: rc,
    };

    let t2_limit = ControlLimit::from_parts(t2_ucl, t2_alpha, t2_description.to_string());
    let spe_limit = ControlLimit::from_parts(spe_ucl, spe_alpha, spe_description.to_string());

    let frcc_config = fdars_core::spm::FrccConfig {
        ncomp: nc,
        fosr_lambda: config_fosr_lambda,
        alpha: config_alpha,
        tuning_fraction: config_tuning_fraction,
        seed: config_seed as u64,
        min_r_squared: 0.0,
    };

    let chart = fdars_core::spm::FrccChart::from_parts(
        fosr, residual_fpca, eigenvalues, t2_limit, spe_limit, 0.0, frcc_config,
    );

    match fdars_core::spm::frcc_monitor(&chart, &ny, &np, &argvals) {
        Ok(result) => {
            let t2_alarm_int: Vec<i32> = result.t2_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            let spe_alarm_int: Vec<i32> = result.spe_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            r!(list!(
                t2 = Robj::from(result.t2),
                spe = Robj::from(result.spe),
                t2_alarm = Robj::from(t2_alarm_int),
                spe_alarm = Robj::from(spe_alarm_int),
                residual_scores = fdmatrix_to_robj(&result.residual_scores)
            ))
        }
        Err(e) => {
            rprintln!("frcc_monitor failed: {:?}", e);
            r!(NULL)
        }
    }
}

// =============================================================================
// Phase 6-9: Alignment & Shape Analysis Wrappers
// =============================================================================

/// Lambda cross-validation for elastic alignment regularisation
#[extendr]
fn alignment_lambda_cv_rust(data: RMatrix<f64>, argvals: Vec<f64>, lambdas: Vec<f64>, n_folds: i32, max_iter: i32, tol: f64, seed: i32) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return r!(NULL);
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let config = alignment::LambdaCvConfig {
        lambdas,
        n_folds: n_folds as usize,
        max_iter: max_iter as usize,
        tol,
        seed: seed as u64,
    };
    match alignment::lambda_cv(&fd, &argvals, &config) {
        Ok(res) => {
            list!(
                best_lambda = res.best_lambda,
                cv_scores = res.cv_scores,
                lambdas = res.lambdas
            ).into()
        }
        Err(e) => {
            rprintln!("alignment_lambda_cv failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Warp statistics: mean, variance, confidence bands, Karcher mean warp
#[extendr]
fn alignment_warp_statistics_rust(gammas: RMatrix<f64>, argvals: Vec<f64>, confidence_level: f64) -> Robj {
    let n = gammas.nrows();
    let m = gammas.ncols();
    if n < 2 || m < 2 || argvals.len() != m {
        return r!(NULL);
    }
    let data_slice = gammas.as_real_slice().unwrap();
    let gammas_fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    match alignment::warp_statistics(&gammas_fd, &argvals, confidence_level) {
        Ok(res) => {
            list!(
                mean = res.mean,
                variance = res.variance,
                std_dev = res.std_dev,
                lower_band = res.lower_band,
                upper_band = res.upper_band,
                karcher_mean_warp = res.karcher_mean_warp,
                geodesic_distances = res.geodesic_distances
            ).into()
        }
        Err(e) => {
            rprintln!("alignment_warp_statistics failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Phase boxplot for warping functions
#[extendr]
fn alignment_phase_boxplot_rust(gammas: RMatrix<f64>, argvals: Vec<f64>, factor: f64) -> Robj {
    let n = gammas.nrows();
    let m = gammas.ncols();
    if n < 3 || m < 2 || argvals.len() != m {
        return r!(NULL);
    }
    let data_slice = gammas.as_real_slice().unwrap();
    let gammas_fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    match alignment::phase_boxplot(&gammas_fd, &argvals, factor) {
        Ok(res) => {
            let outlier_indices_r: Vec<i32> = res.outlier_indices.iter().map(|&i| (i + 1) as i32).collect();
            list!(
                median = res.median,
                median_index = (res.median_index + 1) as i32,
                central_lower = res.central_lower,
                central_upper = res.central_upper,
                whisker_lower = res.whisker_lower,
                whisker_upper = res.whisker_upper,
                depths = res.depths,
                outlier_indices = outlier_indices_r,
                factor = res.factor
            ).into()
        }
        Err(e) => {
            rprintln!("alignment_phase_boxplot failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Invert a warping function
#[extendr]
fn alignment_invert_warp_rust(gamma: Vec<f64>, argvals: Vec<f64>) -> Robj {
    if gamma.len() != argvals.len() || argvals.len() < 2 {
        return r!(NULL);
    }
    match alignment::invert_warp(&gamma, &argvals) {
        Ok(inv) => Robj::from(inv),
        Err(e) => {
            rprintln!("alignment_invert_warp failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Warp inverse error: max |gamma(gamma_inv(t)) - t|
#[extendr]
fn alignment_warp_inverse_error_rust(gamma: Vec<f64>, argvals: Vec<f64>) -> Robj {
    if gamma.len() != argvals.len() || argvals.len() < 2 {
        return Robj::from(f64::NAN);
    }
    match alignment::invert_warp(&gamma, &argvals) {
        Ok(gamma_inv) => {
            let err = alignment::warp_inverse_error(&gamma, &gamma_inv, &argvals);
            Robj::from(err)
        }
        Err(e) => {
            rprintln!("alignment_warp_inverse_error failed: {:?}", e);
            Robj::from(f64::NAN)
        }
    }
}

/// Penalized elastic alignment with configurable penalty type
#[extendr]
fn alignment_elastic_pair_penalized_rust(f1: Vec<f64>, f2: Vec<f64>, argvals: Vec<f64>, penalty_type: &str, lambda: f64) -> Robj {
    if f1.len() != argvals.len() || f2.len() != argvals.len() || argvals.len() < 2 {
        return list!(gamma = Vec::<f64>::new(), f_aligned = Vec::<f64>::new(), distance = f64::NAN).into();
    }
    let pt = match penalty_type {
        "first_order" => alignment::WarpPenaltyType::FirstOrder,
        "second_order" => alignment::WarpPenaltyType::SecondOrder,
        "combined" => alignment::WarpPenaltyType::Combined { second_order_weight: lambda.max(0.01) },
        _ => alignment::WarpPenaltyType::FirstOrder,
    };
    let res = alignment::elastic_align_pair_penalized(&f1, &f2, &argvals, lambda, pt);
    list!(gamma = res.gamma, f_aligned = res.f_aligned, distance = res.distance).into()
}

/// Diagnose alignment quality for all curves after Karcher mean computation
#[extendr]
fn alignment_diagnose_rust(
    data: RMatrix<f64>,
    karcher_mean: Vec<f64>,
    karcher_mean_srsf: Vec<f64>,
    gammas: RMatrix<f64>,
    aligned_data: RMatrix<f64>,
    n_iter: i32,
    converged: bool,
    argvals: Vec<f64>,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return r!(NULL);
    }

    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid data dimensions");

    let gs = gammas.as_real_slice().unwrap();
    let gammas_fd = FdMatrix::from_slice(gs, gammas.nrows(), gammas.ncols()).expect("Invalid gammas dimensions");

    let as_ = aligned_data.as_real_slice().unwrap();
    let aligned_fd = FdMatrix::from_slice(as_, aligned_data.nrows(), aligned_data.ncols()).expect("Invalid aligned_data dimensions");

    let karcher = alignment::KarcherMeanResult::new(
        karcher_mean,
        karcher_mean_srsf,
        gammas_fd,
        aligned_fd,
        n_iter as usize,
        converged,
        None,
    );

    let config = alignment::DiagnosticConfig::default();
    match alignment::diagnose_alignment(&fd, &karcher, &argvals, &config) {
        Ok(summary) => {
            // Build per-curve diagnostic lists
            let diag_list: Vec<Robj> = summary.diagnostics.iter().map(|d| {
                let issues_vec: Vec<String> = d.issues.clone();
                list!(
                    curve_index = (d.curve_index + 1) as i32,
                    warp_complexity = d.warp_complexity,
                    warp_smoothness = d.warp_smoothness,
                    is_under_aligned = if d.is_under_aligned { 1i32 } else { 0i32 },
                    is_over_aligned = if d.is_over_aligned { 1i32 } else { 0i32 },
                    has_non_monotone = if d.has_non_monotone { 1i32 } else { 0i32 },
                    residual = d.residual,
                    distance_ratio = d.distance_ratio,
                    flagged = if d.flagged { 1i32 } else { 0i32 },
                    issues = issues_vec
                ).into()
            }).collect();

            let flagged_r: Vec<i32> = summary.flagged_indices.iter().map(|&i| (i + 1) as i32).collect();

            list!(
                diagnostics = List::from_values(diag_list),
                flagged_indices = flagged_r,
                n_flagged = summary.n_flagged as i32,
                health_score = summary.health_score
            ).into()
        }
        Err(e) => {
            rprintln!("alignment_diagnose failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Diagnose a single pairwise alignment
#[extendr]
fn alignment_diagnose_pairwise_rust(
    f1: Vec<f64>,
    f2: Vec<f64>,
    gamma: Vec<f64>,
    f_aligned: Vec<f64>,
    distance: f64,
    argvals: Vec<f64>,
) -> Robj {
    if f1.len() != argvals.len() || f2.len() != argvals.len() || argvals.len() < 2 {
        return r!(NULL);
    }

    let alignment_result = alignment::AlignmentResult::from_parts(gamma, f_aligned, distance);

    let config = alignment::DiagnosticConfig::default();
    let d = alignment::diagnose_pairwise(&f1, &f2, &alignment_result, &argvals, &config);

    let issues_vec: Vec<String> = d.issues.clone();
    list!(
        curve_index = (d.curve_index + 1) as i32,
        warp_complexity = d.warp_complexity,
        warp_smoothness = d.warp_smoothness,
        is_under_aligned = if d.is_under_aligned { 1i32 } else { 0i32 },
        is_over_aligned = if d.is_over_aligned { 1i32 } else { 0i32 },
        has_non_monotone = if d.has_non_monotone { 1i32 } else { 0i32 },
        residual = d.residual,
        distance_ratio = d.distance_ratio,
        flagged = if d.flagged { 1i32 } else { 0i32 },
        issues = issues_vec
    ).into()
}

/// Compute the orbit representative of a curve in a quotient space
#[extendr]
fn alignment_orbit_representative_rust(f: Vec<f64>, argvals: Vec<f64>, quotient: &str) -> Robj {
    if f.len() != argvals.len() || argvals.len() < 2 {
        return r!(NULL);
    }
    let q = match quotient {
        "reparameterization" => alignment::ShapeQuotient::Reparameterization,
        "translation" => alignment::ShapeQuotient::ReparameterizationTranslation,
        "scale" => alignment::ShapeQuotient::ReparameterizationTranslationScale,
        _ => alignment::ShapeQuotient::Reparameterization,
    };
    match alignment::orbit_representative(&f, &argvals, q) {
        Ok(res) => {
            list!(
                representative = res.representative,
                representative_srsf = res.representative_srsf,
                gamma = res.gamma,
                translation = res.translation,
                scale = res.scale
            ).into()
        }
        Err(e) => {
            rprintln!("alignment_orbit_representative failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Elastic shape distance between two curves
#[extendr]
fn alignment_shape_distance_rust(f1: Vec<f64>, f2: Vec<f64>, argvals: Vec<f64>, quotient: &str, lambda: f64) -> Robj {
    if f1.len() != argvals.len() || f2.len() != argvals.len() || argvals.len() < 2 {
        return r!(NULL);
    }
    let q = match quotient {
        "reparameterization" => alignment::ShapeQuotient::Reparameterization,
        "translation" => alignment::ShapeQuotient::ReparameterizationTranslation,
        "scale" => alignment::ShapeQuotient::ReparameterizationTranslationScale,
        _ => alignment::ShapeQuotient::Reparameterization,
    };
    match alignment::shape_distance(&f1, &f2, &argvals, q, lambda) {
        Ok(res) => {
            list!(
                distance = res.distance,
                gamma = res.gamma,
                f2_aligned = res.f2_aligned
            ).into()
        }
        Err(e) => {
            rprintln!("alignment_shape_distance failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Shape mean (Karcher mean in quotient space)
#[extendr]
fn alignment_shape_mean_rust(data: RMatrix<f64>, argvals: Vec<f64>, quotient: &str, lambda: f64, max_iter: i32, tol: f64) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return r!(NULL);
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let q = match quotient {
        "reparameterization" => alignment::ShapeQuotient::Reparameterization,
        "translation" => alignment::ShapeQuotient::ReparameterizationTranslation,
        "scale" => alignment::ShapeQuotient::ReparameterizationTranslationScale,
        _ => alignment::ShapeQuotient::Reparameterization,
    };
    match alignment::shape_mean(&fd, &argvals, q, lambda, max_iter as usize, tol) {
        Ok(res) => {
            list!(
                mean = res.mean,
                mean_srsf = res.mean_srsf,
                gammas = fdmatrix_to_robj(&res.gammas),
                aligned_data = fdmatrix_to_robj(&res.aligned_data),
                n_iter = res.n_iter as i32,
                converged = res.converged
            ).into()
        }
        Err(e) => {
            rprintln!("alignment_shape_mean failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Shape self-distance matrix
#[extendr]
fn alignment_shape_self_distance_matrix_rust(data: RMatrix<f64>, argvals: Vec<f64>, quotient: &str, lambda: f64) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return r!(RMatrix::new_matrix(0, 0, |_, _| 0.0));
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let q = match quotient {
        "reparameterization" => alignment::ShapeQuotient::Reparameterization,
        "translation" => alignment::ShapeQuotient::ReparameterizationTranslation,
        "scale" => alignment::ShapeQuotient::ReparameterizationTranslationScale,
        _ => alignment::ShapeQuotient::Reparameterization,
    };
    match alignment::shape_self_distance_matrix(&fd, &argvals, q, lambda) {
        Ok(result) => fdmatrix_to_robj(&result),
        Err(e) => {
            rprintln!("alignment_shape_self_distance_matrix failed: {:?}", e);
            r!(RMatrix::new_matrix(0, 0, |_, _| 0.0))
        }
    }
}

/// Elastic k-means clustering
#[extendr]
fn elastic_kmeans_rust(
    data: RMatrix<f64>,
    argvals: Vec<f64>,
    k: i32,
    max_iter: i32,
    tol: f64,
    karcher_max_iter: i32,
    karcher_tol: f64,
    lambda: f64,
    seed: i32,
) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n == 0 || m < 2 || argvals.len() != m {
        return r!(NULL);
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let config = alignment::ElasticClusterConfig {
        k: k as usize,
        lambda,
        max_iter: max_iter as usize,
        tol,
        karcher_max_iter: karcher_max_iter as usize,
        karcher_tol,
        seed: seed as u64,
    };
    match alignment::elastic_kmeans(&fd, &argvals, &config) {
        Ok(res) => {
            // Convert 0-indexed labels to 1-indexed for R
            let labels_r: Vec<i32> = res.labels.iter().map(|&l| (l + 1) as i32).collect();

            // Convert centers (KarcherMeanResult) to R list
            let centers_list: Vec<Robj> = res.centers.iter().map(|c| {
                list!(
                    mean = c.mean.clone(),
                    mean_srsf = c.mean_srsf.clone(),
                    gammas = fdmatrix_to_robj(&c.gammas),
                    aligned_data = fdmatrix_to_robj(&c.aligned_data),
                    n_iter = c.n_iter as i32,
                    converged = c.converged
                ).into()
            }).collect();

            list!(
                labels = labels_r,
                centers = List::from_values(centers_list),
                within_distances = res.within_distances,
                total_within_distance = res.total_within_distance,
                n_iter = res.n_iter as i32,
                converged = res.converged
            ).into()
        }
        Err(e) => {
            rprintln!("elastic_kmeans failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Elastic hierarchical clustering
#[extendr]
fn elastic_hierarchical_rust(data: RMatrix<f64>, argvals: Vec<f64>, method: &str, lambda: f64) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    if n < 2 || m < 2 || argvals.len() != m {
        return r!(NULL);
    }
    let data_slice = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(data_slice, n, m).expect("Invalid matrix dimensions");
    let cluster_method = match method {
        "complete" => alignment::ElasticClusterMethod::HierarchicalComplete,
        "single" => alignment::ElasticClusterMethod::HierarchicalSingle,
        "average" => alignment::ElasticClusterMethod::HierarchicalAverage,
        _ => alignment::ElasticClusterMethod::HierarchicalSingle,
    };
    match alignment::elastic_hierarchical(&fd, &argvals, cluster_method, lambda) {
        Ok(dendro) => {
            // Convert merges to a 3-column matrix (i, j, distance) with 1-indexed
            let n_merges = dendro.merges.len();
            let merges_mat = RMatrix::new_matrix(n_merges, 3, |i, j| {
                match j {
                    0 => (dendro.merges[i].0 + 1) as f64,
                    1 => (dendro.merges[i].1 + 1) as f64,
                    2 => dendro.merges[i].2,
                    _ => 0.0,
                }
            });
            list!(
                merges = merges_mat,
                distance_matrix = fdmatrix_to_robj(&dendro.distance_matrix)
            ).into()
        }
        Err(e) => {
            rprintln!("elastic_hierarchical failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Cut a hierarchical dendrogram to get cluster labels
#[extendr]
fn elastic_cut_dendrogram_rust(merges_data: RMatrix<f64>, distance_matrix: RMatrix<f64>, k: i32) -> Robj {
    let n_merges = merges_data.nrows();
    let n = distance_matrix.nrows();
    let m_dist = distance_matrix.ncols();
    if n == 0 || n != m_dist || n_merges == 0 {
        return r!(NULL);
    }

    // Reconstruct merges from the 3-column matrix (1-indexed back to 0-indexed)
    let merges_slice = merges_data.as_real_slice().unwrap();
    let mut merges: Vec<(usize, usize, f64)> = Vec::with_capacity(n_merges);
    for i in 0..n_merges {
        let ci = merges_slice[i] as usize - 1;           // column 0
        let cj = merges_slice[n_merges + i] as usize - 1; // column 1
        let dist = merges_slice[2 * n_merges + i];         // column 2
        merges.push((ci, cj, dist));
    }

    // Reconstruct distance_matrix
    let dist_slice = distance_matrix.as_real_slice().unwrap();
    let dist_fd = FdMatrix::from_slice(dist_slice, n, n).expect("Invalid distance matrix dimensions");

    let dendro = alignment::ElasticDendrogram::from_parts(merges, dist_fd);

    match alignment::cut_dendrogram(&dendro, k as usize) {
        Ok(labels) => {
            // Convert 0-indexed labels to 1-indexed for R
            let labels_r: Vec<i32> = labels.iter().map(|&l| (l + 1) as i32).collect();
            Robj::from(labels_r)
        }
        Err(e) => {
            rprintln!("elastic_cut_dendrogram failed: {:?}", e);
            r!(NULL)
        }
    }
}

// =============================================================================
// SPM v0.9.0 wrappers: ncomp, contrib, rules, bootstrap, ARL, CUSUM,
// MEWMA, AMEWMA, iterative, profile, partial, elastic SPM
// =============================================================================

/// Helper: reconstruct an SpmChart from its serialized parts.
fn reconstruct_spm_chart(
    rotation: RMatrix<f64>, mean: Vec<f64>, singular_values: Vec<f64>,
    centered: RMatrix<f64>, eigenvalues: Vec<f64>,
    t2_ucl: f64, t2_alpha: f64, t2_description: &str,
    spe_ucl: f64, spe_alpha: f64, spe_description: &str,
    ncomp: i32, config_alpha: f64, config_tuning_fraction: f64, config_seed: i32,
) -> fdars_core::spm::phase::SpmChart {
    use fdars_core::regression::FpcaResult;
    use fdars_core::spm::phase::SpmChart;
    use fdars_core::spm::ControlLimit;

    let rs = rotation.as_real_slice().unwrap();
    let rot = FdMatrix::from_slice(rs, rotation.nrows(), rotation.ncols()).expect("Invalid rotation");
    let cs = centered.as_real_slice().unwrap();
    let cen = FdMatrix::from_slice(cs, centered.nrows(), centered.ncols()).expect("Invalid centered");

    let n_tune = centered.nrows();
    let nc = ncomp as usize;
    let scores_placeholder = FdMatrix::zeros(n_tune, nc);

    let fpca = FpcaResult {
        singular_values,
        rotation: rot,
        scores: scores_placeholder,
        mean,
        centered: cen,
    };

    let t2_limit = ControlLimit::from_parts(t2_ucl, t2_alpha, t2_description.to_string());
    let spe_limit = ControlLimit::from_parts(spe_ucl, spe_alpha, spe_description.to_string());
    let config = fdars_core::spm::SpmConfig {
        ncomp: nc,
        alpha: config_alpha,
        tuning_fraction: config_tuning_fraction,
        seed: config_seed as u64,
    };

    SpmChart::from_parts(
        fpca, eigenvalues, vec![], vec![], t2_limit, spe_limit, config, true,
    )
}

// --- Phase 1: ncomp, contrib, rules, bootstrap ---

/// Select the number of principal components.
#[extendr]
fn spm_select_ncomp_rust(eigenvalues: Vec<f64>, method: &str, threshold: f64) -> Robj {
    use fdars_core::spm::ncomp::{select_ncomp, NcompMethod};

    let method_enum = match method {
        "variance90" | "cumulative_variance" => NcompMethod::CumulativeVariance(threshold),
        "kaiser" => NcompMethod::Kaiser,
        "elbow" => NcompMethod::Elbow,
        "fixed" => NcompMethod::Fixed(threshold as usize),
        _ => NcompMethod::CumulativeVariance(threshold),
    };

    match select_ncomp(&eigenvalues, &method_enum) {
        Ok(k) => r!(k as i32),
        Err(e) => {
            rprintln!("select_ncomp failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Per-PC T-squared contributions for a single observation or batch.
#[extendr]
fn spm_t2_pc_contrib_rust(scores: RMatrix<f64>, eigenvalues: Vec<f64>) -> Robj {
    let n = scores.nrows();
    let nc = scores.ncols();
    let ss = scores.as_real_slice().unwrap();
    let sc = FdMatrix::from_slice(ss, n, nc).expect("Invalid scores");

    match fdars_core::spm::t2_pc_contributions(&sc, &eigenvalues) {
        Ok(contrib) => fdmatrix_to_robj(&contrib),
        Err(e) => {
            rprintln!("t2_pc_contributions failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Evaluate control chart rules (Western Electric, Nelson, or custom).
#[extendr]
fn spm_evaluate_rules_rust(values: Vec<f64>, center: f64, sigma: f64, rule_set: &str) -> Robj {
    use fdars_core::spm::rules::{western_electric_rules, nelson_rules, evaluate_rules, ChartRule};

    let violations = match rule_set {
        "western_electric" => western_electric_rules(&values, center, sigma),
        "nelson" => nelson_rules(&values, center, sigma),
        _ => evaluate_rules(&values, center, sigma, &[
            ChartRule::WE1, ChartRule::WE2, ChartRule::WE3, ChartRule::WE4,
        ]),
    };

    match violations {
        Ok(viols) => {
            let rule_names: Vec<String> = viols.iter().map(|v| {
                match &v.rule {
                    ChartRule::WE1 => "WE1".to_string(),
                    ChartRule::WE2 => "WE2".to_string(),
                    ChartRule::WE3 => "WE3".to_string(),
                    ChartRule::WE4 => "WE4".to_string(),
                    ChartRule::Nelson5 => "Nelson5".to_string(),
                    ChartRule::Nelson6 => "Nelson6".to_string(),
                    ChartRule::Nelson7 => "Nelson7".to_string(),
                    ChartRule::CustomRun { n_points, k_sigma } =>
                        format!("CustomRun({},{})", n_points, k_sigma),
                    _ => "Unknown".to_string(),
                }
            }).collect();

            let indices_list = List::from_values(
                viols.iter().map(|v| {
                    let idx: Vec<i32> = v.indices.iter().map(|&i| (i + 1) as i32).collect();
                    Robj::from(idx)
                })
            );

            r!(list!(
                rules = Robj::from(rule_names),
                indices = Robj::from(indices_list)
            ))
        }
        Err(e) => {
            rprintln!("evaluate_rules failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Robust T-squared control limit.
#[extendr]
fn spm_t2_limit_robust_rust(t2_values: Vec<f64>, ncomp: i32, alpha: f64,
                             method: &str, n_bootstrap: i32, seed: i32) -> Robj {
    use fdars_core::spm::bootstrap::{t2_limit_robust, ControlLimitMethod};

    let method_enum = match method {
        "parametric" => ControlLimitMethod::Parametric,
        "empirical" => ControlLimitMethod::Empirical,
        "bootstrap" => ControlLimitMethod::Bootstrap {
            n_bootstrap: n_bootstrap as usize,
            seed: seed as u64,
        },
        "kde" => ControlLimitMethod::KernelDensity { bandwidth: None },
        _ => ControlLimitMethod::Parametric,
    };

    match t2_limit_robust(&t2_values, ncomp as usize, alpha, &method_enum) {
        Ok(limit) => {
            r!(list!(
                ucl = limit.ucl,
                alpha = limit.alpha,
                description = limit.description.clone()
            ))
        }
        Err(e) => {
            rprintln!("t2_limit_robust failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Robust SPE control limit.
#[extendr]
fn spm_spe_limit_robust_rust(spe_values: Vec<f64>, alpha: f64,
                              method: &str, n_bootstrap: i32, seed: i32) -> Robj {
    use fdars_core::spm::bootstrap::{spe_limit_robust, ControlLimitMethod};

    let method_enum = match method {
        "parametric" => ControlLimitMethod::Parametric,
        "empirical" => ControlLimitMethod::Empirical,
        "bootstrap" => ControlLimitMethod::Bootstrap {
            n_bootstrap: n_bootstrap as usize,
            seed: seed as u64,
        },
        "kde" => ControlLimitMethod::KernelDensity { bandwidth: None },
        _ => ControlLimitMethod::Parametric,
    };

    match spe_limit_robust(&spe_values, alpha, &method_enum) {
        Ok(limit) => {
            r!(list!(
                ucl = limit.ucl,
                alpha = limit.alpha,
                description = limit.description.clone()
            ))
        }
        Err(e) => {
            rprintln!("spe_limit_robust failed: {:?}", e);
            r!(NULL)
        }
    }
}

// --- Phase 2: ARL and CUSUM ---

/// In-control ARL for T-squared chart.
#[extendr]
fn spm_arl0_t2_rust(eigenvalues: Vec<f64>, ucl: f64,
                     n_sim: i32, max_rl: i32, seed: i32) -> Robj {
    use fdars_core::spm::arl::{arl0_t2, ArlConfig};

    let config = ArlConfig {
        n_simulations: n_sim as usize,
        max_run_length: max_rl as usize,
        seed: seed as u64,
    };

    match arl0_t2(&eigenvalues, ucl, &config) {
        Ok(result) => {
            let rl_i32: Vec<i32> = result.run_lengths.iter().map(|&r| r as i32).collect();
            r!(list!(
                arl = result.arl,
                std_dev = result.std_dev,
                median_rl = result.median_rl,
                run_lengths = Robj::from(rl_i32)
            ))
        }
        Err(e) => {
            rprintln!("arl0_t2 failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Out-of-control ARL for T-squared chart.
#[extendr]
fn spm_arl1_t2_rust(eigenvalues: Vec<f64>, ucl: f64, shift: Vec<f64>,
                     n_sim: i32, max_rl: i32, seed: i32) -> Robj {
    use fdars_core::spm::arl::{arl1_t2, ArlConfig};

    let config = ArlConfig {
        n_simulations: n_sim as usize,
        max_run_length: max_rl as usize,
        seed: seed as u64,
    };

    match arl1_t2(&eigenvalues, ucl, &shift, &config) {
        Ok(result) => {
            let rl_i32: Vec<i32> = result.run_lengths.iter().map(|&r| r as i32).collect();
            r!(list!(
                arl = result.arl,
                std_dev = result.std_dev,
                median_rl = result.median_rl,
                run_lengths = Robj::from(rl_i32)
            ))
        }
        Err(e) => {
            rprintln!("arl1_t2 failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// In-control ARL for EWMA-T-squared chart.
#[extendr]
fn spm_arl0_ewma_t2_rust(eigenvalues: Vec<f64>, ucl: f64, lambda: f64,
                          n_sim: i32, max_rl: i32, seed: i32) -> Robj {
    use fdars_core::spm::arl::{arl0_ewma_t2, ArlConfig};

    let config = ArlConfig {
        n_simulations: n_sim as usize,
        max_run_length: max_rl as usize,
        seed: seed as u64,
    };

    match arl0_ewma_t2(&eigenvalues, ucl, lambda, &config) {
        Ok(result) => {
            let rl_i32: Vec<i32> = result.run_lengths.iter().map(|&r| r as i32).collect();
            r!(list!(
                arl = result.arl,
                std_dev = result.std_dev,
                median_rl = result.median_rl,
                run_lengths = Robj::from(rl_i32)
            ))
        }
        Err(e) => {
            rprintln!("arl0_ewma_t2 failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// In-control ARL for SPE chart.
#[extendr]
fn spm_arl0_spe_rust(spe_df: f64, spe_scale: f64, ucl: f64,
                      n_sim: i32, max_rl: i32, seed: i32) -> Robj {
    use fdars_core::spm::arl::{arl0_spe, ArlConfig};

    let config = ArlConfig {
        n_simulations: n_sim as usize,
        max_run_length: max_rl as usize,
        seed: seed as u64,
    };

    match arl0_spe(spe_df, spe_scale, ucl, &config) {
        Ok(result) => {
            let rl_i32: Vec<i32> = result.run_lengths.iter().map(|&r| r as i32).collect();
            r!(list!(
                arl = result.arl,
                std_dev = result.std_dev,
                median_rl = result.median_rl,
                run_lengths = Robj::from(rl_i32)
            ))
        }
        Err(e) => {
            rprintln!("arl0_spe failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// CUSUM monitoring on sequential functional data.
#[extendr]
fn spm_cusum_monitor_rust(
    sequential_data: RMatrix<f64>, argvals: Vec<f64>,
    rotation: RMatrix<f64>, mean: Vec<f64>, singular_values: Vec<f64>,
    centered: RMatrix<f64>, eigenvalues: Vec<f64>,
    t2_ucl: f64, t2_alpha: f64, t2_description: &str,
    spe_ucl: f64, spe_alpha: f64, spe_description: &str,
    chart_ncomp: i32, config_alpha: f64, config_tuning_fraction: f64, config_seed: i32,
    k: f64, h: f64, cusum_ncomp: i32, cusum_alpha: f64, multivariate: bool,
) -> Robj {
    let chart = reconstruct_spm_chart(
        rotation, mean, singular_values, centered, eigenvalues,
        t2_ucl, t2_alpha, t2_description,
        spe_ucl, spe_alpha, spe_description,
        chart_ncomp, config_alpha, config_tuning_fraction, config_seed,
    );

    let n = sequential_data.nrows();
    let m = sequential_data.ncols();
    let sds = sequential_data.as_real_slice().unwrap();
    let sd = FdMatrix::from_slice(sds, n, m).expect("Invalid sequential_data");

    let cusum_config = fdars_core::spm::CusumConfig {
        k,
        h,
        ncomp: cusum_ncomp as usize,
        alpha: cusum_alpha,
        multivariate,
        restart: false,
    };

    match fdars_core::spm::spm_cusum_monitor(&chart, &sd, &argvals, &cusum_config) {
        Ok(result) => {
            let alarm_int: Vec<i32> = result.alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            let spe_alarm_int: Vec<i32> = result.spe_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            r!(list!(
                cusum_statistic = Robj::from(result.cusum_statistic),
                ucl = result.ucl,
                alarm = Robj::from(alarm_int),
                scores = fdmatrix_to_robj(&result.scores),
                spe = Robj::from(result.spe),
                spe_limit = result.spe_limit,
                spe_alarm = Robj::from(spe_alarm_int)
            ))
        }
        Err(e) => {
            rprintln!("spm_cusum_monitor failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// CUSUM monitoring with restart after each alarm.
#[extendr]
fn spm_cusum_monitor_restart_rust(
    sequential_data: RMatrix<f64>, argvals: Vec<f64>,
    rotation: RMatrix<f64>, mean: Vec<f64>, singular_values: Vec<f64>,
    centered: RMatrix<f64>, eigenvalues: Vec<f64>,
    t2_ucl: f64, t2_alpha: f64, t2_description: &str,
    spe_ucl: f64, spe_alpha: f64, spe_description: &str,
    chart_ncomp: i32, config_alpha: f64, config_tuning_fraction: f64, config_seed: i32,
    k: f64, h: f64, cusum_ncomp: i32, cusum_alpha: f64, multivariate: bool,
) -> Robj {
    let chart = reconstruct_spm_chart(
        rotation, mean, singular_values, centered, eigenvalues,
        t2_ucl, t2_alpha, t2_description,
        spe_ucl, spe_alpha, spe_description,
        chart_ncomp, config_alpha, config_tuning_fraction, config_seed,
    );

    let n = sequential_data.nrows();
    let m = sequential_data.ncols();
    let sds = sequential_data.as_real_slice().unwrap();
    let sd = FdMatrix::from_slice(sds, n, m).expect("Invalid sequential_data");

    let cusum_config = fdars_core::spm::CusumConfig {
        k,
        h,
        ncomp: cusum_ncomp as usize,
        alpha: cusum_alpha,
        multivariate,
        restart: true,
    };

    match fdars_core::spm::spm_cusum_monitor_with_restart(&chart, &sd, &argvals, &cusum_config) {
        Ok(result) => {
            let alarm_int: Vec<i32> = result.alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            let spe_alarm_int: Vec<i32> = result.spe_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            r!(list!(
                cusum_statistic = Robj::from(result.cusum_statistic),
                ucl = result.ucl,
                alarm = Robj::from(alarm_int),
                scores = fdmatrix_to_robj(&result.scores),
                spe = Robj::from(result.spe),
                spe_limit = result.spe_limit,
                spe_alarm = Robj::from(spe_alarm_int)
            ))
        }
        Err(e) => {
            rprintln!("spm_cusum_monitor_with_restart failed: {:?}", e);
            r!(NULL)
        }
    }
}

// --- Phase 3: MEWMA, AMEWMA, Iterative ---

/// MEWMA monitoring on sequential functional data.
#[extendr]
fn spm_mewma_monitor_rust(
    sequential_data: RMatrix<f64>, argvals: Vec<f64>,
    rotation: RMatrix<f64>, mean: Vec<f64>, singular_values: Vec<f64>,
    centered: RMatrix<f64>, eigenvalues: Vec<f64>,
    t2_ucl: f64, t2_alpha: f64, t2_description: &str,
    spe_ucl: f64, spe_alpha: f64, spe_description: &str,
    chart_ncomp: i32, config_alpha: f64, config_tuning_fraction: f64, config_seed: i32,
    lambda: f64, mewma_ncomp: i32, mewma_alpha: f64, asymptotic: bool,
) -> Robj {
    let chart = reconstruct_spm_chart(
        rotation, mean, singular_values, centered, eigenvalues,
        t2_ucl, t2_alpha, t2_description,
        spe_ucl, spe_alpha, spe_description,
        chart_ncomp, config_alpha, config_tuning_fraction, config_seed,
    );

    let n = sequential_data.nrows();
    let m = sequential_data.ncols();
    let sds = sequential_data.as_real_slice().unwrap();
    let sd = FdMatrix::from_slice(sds, n, m).expect("Invalid sequential_data");

    let mewma_config = fdars_core::spm::MewmaConfig {
        lambda,
        ncomp: mewma_ncomp as usize,
        alpha: mewma_alpha,
        asymptotic,
    };

    match fdars_core::spm::spm_mewma_monitor(&chart, &sd, &argvals, &mewma_config) {
        Ok(result) => {
            let alarm_int: Vec<i32> = result.alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            let spe_alarm_int: Vec<i32> = result.spe_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            r!(list!(
                smoothed_scores = fdmatrix_to_robj(&result.smoothed_scores),
                mewma_statistic = Robj::from(result.mewma_statistic),
                ucl = result.ucl,
                alarm = Robj::from(alarm_int),
                spe = Robj::from(result.spe),
                spe_limit = result.spe_limit,
                spe_alarm = Robj::from(spe_alarm_int)
            ))
        }
        Err(e) => {
            rprintln!("spm_mewma_monitor failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Adaptive EWMA (AMEWMA) monitoring on sequential functional data.
#[extendr]
fn spm_amewma_monitor_rust(
    sequential_data: RMatrix<f64>, argvals: Vec<f64>,
    rotation: RMatrix<f64>, mean: Vec<f64>, singular_values: Vec<f64>,
    centered: RMatrix<f64>, eigenvalues: Vec<f64>,
    t2_ucl: f64, t2_alpha: f64, t2_description: &str,
    spe_ucl: f64, spe_alpha: f64, spe_description: &str,
    chart_ncomp: i32, config_alpha: f64, config_tuning_fraction: f64, config_seed: i32,
    lambda_min: f64, lambda_max: f64, lambda_init: f64, eta: f64,
    amewma_ncomp: i32, amewma_alpha: f64,
) -> Robj {
    let chart = reconstruct_spm_chart(
        rotation, mean, singular_values, centered, eigenvalues,
        t2_ucl, t2_alpha, t2_description,
        spe_ucl, spe_alpha, spe_description,
        chart_ncomp, config_alpha, config_tuning_fraction, config_seed,
    );

    let n = sequential_data.nrows();
    let m = sequential_data.ncols();
    let sds = sequential_data.as_real_slice().unwrap();
    let sd = FdMatrix::from_slice(sds, n, m).expect("Invalid sequential_data");

    let amewma_config = fdars_core::spm::AmewmaConfig {
        lambda_min,
        lambda_max,
        lambda_init,
        eta,
        ncomp: amewma_ncomp as usize,
        alpha: amewma_alpha,
        error_clip: None,
    };

    match fdars_core::spm::spm_amewma_monitor(&chart, &sd, &argvals, &amewma_config) {
        Ok(result) => {
            let alarm_int: Vec<i32> = result.alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            r!(list!(
                smoothed_scores = fdmatrix_to_robj(&result.smoothed_scores),
                t2_statistic = Robj::from(result.t2_statistic),
                lambda_t = Robj::from(result.lambda_t),
                ucl = result.ucl,
                alarm = Robj::from(alarm_int)
            ))
        }
        Err(e) => {
            rprintln!("spm_amewma_monitor failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Iterative Phase I chart construction.
#[extendr]
fn spm_phase1_iterative_rust(data: RMatrix<f64>, argvals: Vec<f64>,
                              ncomp: i32, alpha: f64, tuning_fraction: f64, seed: i32,
                              max_iterations: i32, max_removal_fraction: f64) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");

    let config = fdars_core::spm::IterativePhase1Config {
        spm: fdars_core::spm::SpmConfig {
            ncomp: ncomp as usize,
            alpha,
            tuning_fraction,
            seed: seed as u64,
        },
        max_iterations: max_iterations as usize,
        remove_t2_outliers: true,
        remove_spe_outliers: true,
        max_removal_fraction,
    };

    match fdars_core::spm::spm_phase1_iterative(&fd, &argvals, &config) {
        Ok(result) => {
            let chart = &result.chart;
            let removed_i32: Vec<i32> = result.removed_indices.iter().map(|&i| (i + 1) as i32).collect();

            // Convert removal_history to a list of integer vectors (1-indexed)
            let history_list = List::from_values(
                result.removal_history.iter().map(|v| {
                    let vi: Vec<i32> = v.iter().map(|&i| (i + 1) as i32).collect();
                    Robj::from(vi)
                })
            );

            r!(list!(
                rotation = fdmatrix_to_robj(&chart.fpca.rotation),
                scores = fdmatrix_to_robj(&chart.fpca.scores),
                mean = Robj::from(chart.fpca.mean.clone()),
                singular_values = Robj::from(chart.fpca.singular_values.clone()),
                centered = fdmatrix_to_robj(&chart.fpca.centered),
                eigenvalues = Robj::from(chart.eigenvalues.clone()),
                t2_phase1 = Robj::from(chart.t2_phase1.clone()),
                spe_phase1 = Robj::from(chart.spe_phase1.clone()),
                t2_ucl = chart.t2_limit.ucl,
                t2_alpha = chart.t2_limit.alpha,
                t2_description = chart.t2_limit.description.clone(),
                spe_ucl = chart.spe_limit.ucl,
                spe_alpha = chart.spe_limit.alpha,
                spe_description = chart.spe_limit.description.clone(),
                ncomp = chart.eigenvalues.len() as i32,
                config_ncomp = config.spm.ncomp as i32,
                config_alpha = config.spm.alpha,
                config_tuning_fraction = config.spm.tuning_fraction,
                config_seed = config.spm.seed as i32,
                n_iterations = result.n_iterations as i32,
                removed_indices = Robj::from(removed_i32),
                n_remaining = result.n_remaining as i32,
                removal_history = Robj::from(history_list),
                removal_rates = Robj::from(result.removal_rates)
            ))
        }
        Err(e) => {
            rprintln!("spm_phase1_iterative failed: {:?}", e);
            r!(NULL)
        }
    }
}

// --- Phase 4: Profile and Partial ---

/// Profile monitoring Phase I.
#[extendr]
fn spm_profile_phase1_rust(y_curves: RMatrix<f64>, predictors: RMatrix<f64>, argvals: Vec<f64>,
                            fosr_lambda: f64, ncomp: i32, alpha: f64,
                            window_size: i32, step_size: i32) -> Robj {
    let n = y_curves.nrows();
    let m = y_curves.ncols();
    let ys = y_curves.as_real_slice().unwrap();
    let y = FdMatrix::from_slice(ys, n, m).expect("Invalid y_curves");

    let ps = predictors.as_real_slice().unwrap();
    let pred = FdMatrix::from_slice(ps, predictors.nrows(), predictors.ncols()).expect("Invalid predictors");

    let config = fdars_core::spm::ProfileMonitorConfig {
        fosr_lambda,
        ncomp: ncomp as usize,
        alpha,
        window_size: window_size as usize,
        step_size: step_size as usize,
    };

    match fdars_core::spm::profile_phase1(&y, &pred, &argvals, &config) {
        Ok(chart) => {
            r!(list!(
                // Reference FOSR
                ref_fosr_intercept = Robj::from(chart.reference_fosr.intercept.clone()),
                ref_fosr_beta = fdmatrix_to_robj(&chart.reference_fosr.beta),
                ref_fosr_lambda = chart.reference_fosr.lambda,
                // Beta FPCA
                beta_rotation = fdmatrix_to_robj(&chart.beta_fpca.rotation),
                beta_scores = fdmatrix_to_robj(&chart.beta_fpca.scores),
                beta_mean = Robj::from(chart.beta_fpca.mean.clone()),
                beta_singular_values = Robj::from(chart.beta_fpca.singular_values.clone()),
                beta_centered = fdmatrix_to_robj(&chart.beta_fpca.centered),
                // Eigenvalues & limits
                eigenvalues = Robj::from(chart.eigenvalues.clone()),
                t2_ucl = chart.t2_limit.ucl,
                t2_alpha = chart.t2_limit.alpha,
                t2_description = chart.t2_limit.description.clone(),
                // Diagnostics
                lag1_autocorrelation = chart.lag1_autocorrelation,
                effective_n_windows = chart.effective_n_windows,
                // Config
                config_fosr_lambda = config.fosr_lambda,
                config_ncomp = config.ncomp as i32,
                config_alpha = config.alpha,
                config_window_size = config.window_size as i32,
                config_step_size = config.step_size as i32
            ))
        }
        Err(e) => {
            rprintln!("profile_phase1 failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Profile monitoring Phase II.
#[extendr]
fn spm_profile_monitor_rust(
    new_y: RMatrix<f64>, new_predictors: RMatrix<f64>, argvals: Vec<f64>,
    // Beta FPCA fields for chart reconstruction
    beta_rotation: RMatrix<f64>, beta_mean: Vec<f64>,
    beta_singular_values: Vec<f64>, beta_centered: RMatrix<f64>,
    // Reference FOSR fields
    ref_fosr_intercept: Vec<f64>, ref_fosr_beta: RMatrix<f64>, ref_fosr_lambda: f64,
    // Eigenvalues & limits
    eigenvalues: Vec<f64>,
    t2_ucl: f64, t2_alpha: f64, t2_description: &str,
    lag1_autocorrelation: f64, effective_n_windows: f64,
    // Config
    fosr_lambda: f64, ncomp: i32, alpha: f64,
    window_size: i32, step_size: i32,
) -> Robj {
    use fdars_core::regression::FpcaResult;
    use fdars_core::spm::ControlLimit;
    use fdars_core::spm::profile::{ProfileChart, ProfileMonitorConfig, profile_monitor};

    let n = new_y.nrows();
    let m = new_y.ncols();
    let nys = new_y.as_real_slice().unwrap();
    let ny = FdMatrix::from_slice(nys, n, m).expect("Invalid new_y");

    let nps = new_predictors.as_real_slice().unwrap();
    let np = FdMatrix::from_slice(nps, new_predictors.nrows(), new_predictors.ncols()).expect("Invalid new_predictors");

    // Reconstruct beta FPCA
    let brs = beta_rotation.as_real_slice().unwrap();
    let br = FdMatrix::from_slice(brs, beta_rotation.nrows(), beta_rotation.ncols()).expect("Invalid beta_rotation");
    let bcs = beta_centered.as_real_slice().unwrap();
    let bc = FdMatrix::from_slice(bcs, beta_centered.nrows(), beta_centered.ncols()).expect("Invalid beta_centered");

    let n_windows = beta_centered.nrows();
    let nc = eigenvalues.len().min(ncomp as usize);
    let beta_scores_placeholder = FdMatrix::zeros(n_windows, nc);

    let beta_fpca = FpcaResult {
        singular_values: beta_singular_values,
        rotation: br,
        scores: beta_scores_placeholder,
        mean: beta_mean,
        centered: bc,
    };

    // Reconstruct reference FOSR (minimal, just intercept + beta needed)
    let fbs = ref_fosr_beta.as_real_slice().unwrap();
    let fb = FdMatrix::from_slice(fbs, ref_fosr_beta.nrows(), ref_fosr_beta.ncols()).expect("Invalid ref_fosr_beta");

    let reference_fosr = fdars_core::function_on_scalar::FosrResult::from_parts(
        ref_fosr_intercept,
        fb,
        FdMatrix::zeros(1, 1),
        FdMatrix::zeros(1, 1),
        vec![],
        0.0,
        FdMatrix::zeros(1, 1),
        ref_fosr_lambda,
        0.0,
    );

    let t2_limit = ControlLimit::from_parts(t2_ucl, t2_alpha, t2_description.to_string());

    let monitor_config = ProfileMonitorConfig {
        fosr_lambda,
        ncomp: ncomp as usize,
        alpha,
        window_size: window_size as usize,
        step_size: step_size as usize,
    };

    let chart = ProfileChart::from_parts(
        reference_fosr, beta_fpca, eigenvalues, t2_limit,
        lag1_autocorrelation, effective_n_windows, monitor_config.clone(),
    );

    match profile_monitor(&chart, &ny, &np, &argvals, &monitor_config) {
        Ok(result) => {
            let t2_alarm_int: Vec<i32> = result.t2_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            r!(list!(
                betas = fdmatrix_to_robj(&result.betas),
                t2 = Robj::from(result.t2),
                t2_alarm = Robj::from(t2_alarm_int),
                beta_scores = fdmatrix_to_robj(&result.beta_scores)
            ))
        }
        Err(e) => {
            rprintln!("profile_monitor failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Monitor a single partially-observed curve.
#[extendr]
fn spm_monitor_partial_rust(
    partial_values: Vec<f64>, argvals: Vec<f64>, n_observed: i32,
    rotation: RMatrix<f64>, mean: Vec<f64>, singular_values: Vec<f64>,
    centered: RMatrix<f64>, eigenvalues: Vec<f64>,
    t2_ucl: f64, t2_alpha: f64, t2_description: &str,
    spe_ucl: f64, spe_alpha: f64, spe_description: &str,
    chart_ncomp: i32, config_alpha: f64, config_tuning_fraction: f64, config_seed: i32,
    partial_ncomp: i32, partial_alpha: f64, completion_method: &str,
) -> Robj {
    let chart = reconstruct_spm_chart(
        rotation, mean, singular_values, centered, eigenvalues,
        t2_ucl, t2_alpha, t2_description,
        spe_ucl, spe_alpha, spe_description,
        chart_ncomp, config_alpha, config_tuning_fraction, config_seed,
    );

    let completion = match completion_method {
        "partial_projection" => fdars_core::spm::DomainCompletion::PartialProjection,
        "conditional_expectation" | "blup" => fdars_core::spm::DomainCompletion::ConditionalExpectation,
        "zero_pad" => fdars_core::spm::DomainCompletion::ZeroPad,
        _ => fdars_core::spm::DomainCompletion::ConditionalExpectation,
    };

    let config = fdars_core::spm::PartialDomainConfig {
        ncomp: partial_ncomp as usize,
        alpha: partial_alpha,
        completion,
        regularization_eps: 1e-10,
    };

    match fdars_core::spm::spm_monitor_partial(&chart, &partial_values, &argvals, n_observed as usize, &config) {
        Ok(result) => {
            let completed = match result.completed_curve {
                Some(ref c) => Robj::from(c.clone()),
                None => r!(NULL),
            };
            r!(list!(
                scores = Robj::from(result.scores),
                t2 = result.t2,
                t2_alarm = if result.t2_alarm { 1 } else { 0 },
                domain_fraction = result.domain_fraction,
                completed_curve = completed
            ))
        }
        Err(e) => {
            rprintln!("spm_monitor_partial failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Monitor a batch of partially-observed curves.
#[extendr]
fn spm_monitor_partial_batch_rust(
    partial_data_list: List, n_observed_vec: Vec<i32>, argvals: Vec<f64>,
    rotation: RMatrix<f64>, mean: Vec<f64>, singular_values: Vec<f64>,
    centered: RMatrix<f64>, eigenvalues: Vec<f64>,
    t2_ucl: f64, t2_alpha: f64, t2_description: &str,
    spe_ucl: f64, spe_alpha: f64, spe_description: &str,
    chart_ncomp: i32, config_alpha: f64, config_tuning_fraction: f64, config_seed: i32,
    partial_ncomp: i32, partial_alpha: f64, completion_method: &str,
) -> Robj {
    let chart = reconstruct_spm_chart(
        rotation, mean, singular_values, centered, eigenvalues,
        t2_ucl, t2_alpha, t2_description,
        spe_ucl, spe_alpha, spe_description,
        chart_ncomp, config_alpha, config_tuning_fraction, config_seed,
    );

    let completion = match completion_method {
        "partial_projection" => fdars_core::spm::DomainCompletion::PartialProjection,
        "conditional_expectation" | "blup" => fdars_core::spm::DomainCompletion::ConditionalExpectation,
        "zero_pad" => fdars_core::spm::DomainCompletion::ZeroPad,
        _ => fdars_core::spm::DomainCompletion::ConditionalExpectation,
    };

    let config = fdars_core::spm::PartialDomainConfig {
        ncomp: partial_ncomp as usize,
        alpha: partial_alpha,
        completion,
        regularization_eps: 1e-10,
    };

    // Collect partial data from the R list
    let mut partial_vecs: Vec<Vec<f64>> = Vec::with_capacity(partial_data_list.len());
    for (_, item) in partial_data_list.iter() {
        let v: Vec<f64> = item.as_real_vector().expect("partial_data_list elements must be numeric vectors");
        partial_vecs.push(v);
    }

    // Build the slice-of-tuples that spm_monitor_partial_batch expects
    let partial_data: Vec<(&[f64], usize)> = partial_vecs.iter().zip(n_observed_vec.iter())
        .map(|(v, &n_obs)| (v.as_slice(), n_obs as usize))
        .collect();

    match fdars_core::spm::spm_monitor_partial_batch(&chart, &partial_data, &argvals, &config) {
        Ok(results) => {
            let result_list = List::from_values(
                results.iter().map(|result| {
                    let completed = match result.completed_curve {
                        Some(ref c) => Robj::from(c.clone()),
                        None => r!(NULL),
                    };
                    r!(list!(
                        scores = Robj::from(result.scores.clone()),
                        t2 = result.t2,
                        t2_alarm = if result.t2_alarm { 1 } else { 0 },
                        domain_fraction = result.domain_fraction,
                        completed_curve = completed
                    ))
                })
            );
            Robj::from(result_list)
        }
        Err(e) => {
            rprintln!("spm_monitor_partial_batch failed: {:?}", e);
            r!(NULL)
        }
    }
}

// --- Phase 5: Elastic SPM ---

/// Elastic SPM Phase I.
#[extendr]
fn elastic_spm_phase1_rust(data: RMatrix<f64>, argvals: Vec<f64>,
                            ncomp: i32, alpha: f64, tuning_fraction: f64, seed: i32,
                            align_lambda: f64, monitor_phase: bool, warp_ncomp: i32) -> Robj {
    let n = data.nrows();
    let m = data.ncols();
    let ds = data.as_real_slice().unwrap();
    let fd = FdMatrix::from_slice(ds, n, m).expect("Invalid data");

    let config = fdars_core::spm::ElasticSpmConfig {
        spm: fdars_core::spm::SpmConfig {
            ncomp: ncomp as usize,
            alpha,
            tuning_fraction,
            seed: seed as u64,
        },
        align_lambda,
        monitor_phase,
        warp_ncomp: warp_ncomp as usize,
        max_karcher_iterations: 20,
        karcher_tolerance: 1e-4,
    };

    match fdars_core::spm::elastic_spm_phase1(&fd, &argvals, &config) {
        Ok(chart) => {
            // Amplitude chart fields
            let amp = &chart.amplitude_chart;
            let mut result_list = list!(
                karcher_mean = Robj::from(chart.karcher_mean.clone()),
                mean_alignment_residual = chart.mean_alignment_residual,
                // Amplitude chart
                amp_rotation = fdmatrix_to_robj(&amp.fpca.rotation),
                amp_scores = fdmatrix_to_robj(&amp.fpca.scores),
                amp_mean = Robj::from(amp.fpca.mean.clone()),
                amp_singular_values = Robj::from(amp.fpca.singular_values.clone()),
                amp_centered = fdmatrix_to_robj(&amp.fpca.centered),
                amp_eigenvalues = Robj::from(amp.eigenvalues.clone()),
                amp_t2_phase1 = Robj::from(amp.t2_phase1.clone()),
                amp_spe_phase1 = Robj::from(amp.spe_phase1.clone()),
                amp_t2_ucl = amp.t2_limit.ucl,
                amp_t2_alpha = amp.t2_limit.alpha,
                amp_t2_description = amp.t2_limit.description.clone(),
                amp_spe_ucl = amp.spe_limit.ucl,
                amp_spe_alpha = amp.spe_limit.alpha,
                amp_spe_description = amp.spe_limit.description.clone(),
                amp_ncomp = amp.eigenvalues.len() as i32,
                // Config
                config_ncomp = config.spm.ncomp as i32,
                config_alpha = config.spm.alpha,
                config_tuning_fraction = config.spm.tuning_fraction,
                config_seed = config.spm.seed as i32,
                config_align_lambda = config.align_lambda,
                config_monitor_phase = if config.monitor_phase { 1 } else { 0 },
                config_warp_ncomp = config.warp_ncomp as i32
            );

            // Phase chart fields (if monitoring phase)
            if let Some(ref phase) = chart.phase_chart {
                let phase_list = list!(
                    rotation = fdmatrix_to_robj(&phase.fpca.rotation),
                    scores = fdmatrix_to_robj(&phase.fpca.scores),
                    mean = Robj::from(phase.fpca.mean.clone()),
                    singular_values = Robj::from(phase.fpca.singular_values.clone()),
                    centered = fdmatrix_to_robj(&phase.fpca.centered),
                    eigenvalues = Robj::from(phase.eigenvalues.clone()),
                    t2_phase1 = Robj::from(phase.t2_phase1.clone()),
                    spe_phase1 = Robj::from(phase.spe_phase1.clone()),
                    t2_ucl = phase.t2_limit.ucl,
                    t2_alpha = phase.t2_limit.alpha,
                    t2_description = phase.t2_limit.description.clone(),
                    spe_ucl = phase.spe_limit.ucl,
                    spe_alpha = phase.spe_limit.alpha,
                    spe_description = phase.spe_limit.description.clone(),
                    ncomp = phase.eigenvalues.len() as i32
                );
                result_list.set_attrib("phase_chart", Robj::from(phase_list)).ok();
            }

            r!(result_list)
        }
        Err(e) => {
            rprintln!("elastic_spm_phase1 failed: {:?}", e);
            r!(NULL)
        }
    }
}

/// Elastic SPM Phase II monitoring.
#[extendr]
fn elastic_spm_monitor_rust(
    new_data: RMatrix<f64>, argvals: Vec<f64>,
    karcher_mean: Vec<f64>, align_lambda: f64, monitor_phase: bool, warp_ncomp: i32,
    // Amplitude chart fields
    amp_rotation: RMatrix<f64>, amp_mean: Vec<f64>, amp_singular_values: Vec<f64>,
    amp_centered: RMatrix<f64>, amp_eigenvalues: Vec<f64>,
    amp_t2_ucl: f64, amp_t2_alpha: f64, amp_t2_description: &str,
    amp_spe_ucl: f64, amp_spe_alpha: f64, amp_spe_description: &str,
    amp_ncomp: i32, amp_config_alpha: f64, amp_config_tuning_fraction: f64, amp_config_seed: i32,
    // Phase chart fields (can be NULL if monitor_phase is false)
    phase_rotation: Robj, phase_mean: Robj, phase_singular_values: Robj,
    phase_centered: Robj, phase_eigenvalues: Robj,
    phase_t2_ucl: f64, phase_t2_alpha: f64, phase_t2_description: &str,
    phase_spe_ucl: f64, phase_spe_alpha: f64, phase_spe_description: &str,
    phase_ncomp: i32, phase_config_alpha: f64, phase_config_tuning_fraction: f64, phase_config_seed: i32,
) -> Robj {
    use fdars_core::spm::elastic_spm::{ElasticSpmChart, ElasticSpmConfig, elastic_spm_monitor};

    let n = new_data.nrows();
    let m = new_data.ncols();
    let nds = new_data.as_real_slice().unwrap();
    let nd = FdMatrix::from_slice(nds, n, m).expect("Invalid new_data");

    // Reconstruct amplitude chart
    let amplitude_chart = reconstruct_spm_chart(
        amp_rotation, amp_mean, amp_singular_values, amp_centered, amp_eigenvalues,
        amp_t2_ucl, amp_t2_alpha, amp_t2_description,
        amp_spe_ucl, amp_spe_alpha, amp_spe_description,
        amp_ncomp, amp_config_alpha, amp_config_tuning_fraction, amp_config_seed,
    );

    // Reconstruct phase chart if monitoring phase
    let phase_chart = if monitor_phase && !phase_rotation.is_null() {
        let pr: RMatrix<f64> = phase_rotation.try_into().expect("phase_rotation must be matrix");
        let pm: Vec<f64> = phase_mean.as_real_vector().expect("phase_mean must be numeric");
        let psv: Vec<f64> = phase_singular_values.as_real_vector().expect("phase_singular_values must be numeric");
        let pc: RMatrix<f64> = phase_centered.try_into().expect("phase_centered must be matrix");
        let pe: Vec<f64> = phase_eigenvalues.as_real_vector().expect("phase_eigenvalues must be numeric");

        Some(reconstruct_spm_chart(
            pr, pm, psv, pc, pe,
            phase_t2_ucl, phase_t2_alpha, phase_t2_description,
            phase_spe_ucl, phase_spe_alpha, phase_spe_description,
            phase_ncomp, phase_config_alpha, phase_config_tuning_fraction, phase_config_seed,
        ))
    } else {
        None
    };

    let spm_config = fdars_core::spm::SpmConfig {
        ncomp: amp_ncomp as usize,
        alpha: amp_config_alpha,
        tuning_fraction: amp_config_tuning_fraction,
        seed: amp_config_seed as u64,
    };

    let elastic_config = ElasticSpmConfig {
        spm: spm_config,
        align_lambda,
        monitor_phase,
        warp_ncomp: warp_ncomp as usize,
        max_karcher_iterations: 20,
        karcher_tolerance: 1e-4,
    };

    let elastic_chart = ElasticSpmChart::from_parts(
        karcher_mean, amplitude_chart, phase_chart, elastic_config, 0.0,
    );

    match elastic_spm_monitor(&elastic_chart, &nd, &argvals) {
        Ok(result) => {
            let amp_t2_alarm_int: Vec<i32> = result.amplitude.t2_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
            let amp_spe_alarm_int: Vec<i32> = result.amplitude.spe_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();

            let mut result_list = list!(
                amp_t2 = Robj::from(result.amplitude.t2),
                amp_spe = Robj::from(result.amplitude.spe),
                amp_t2_alarm = Robj::from(amp_t2_alarm_int),
                amp_spe_alarm = Robj::from(amp_spe_alarm_int),
                amp_scores = fdmatrix_to_robj(&result.amplitude.scores),
                aligned_data = fdmatrix_to_robj(&result.aligned_data),
                warping_functions = fdmatrix_to_robj(&result.warping_functions)
            );

            if let Some(ref phase_result) = result.phase {
                let ph_t2_alarm_int: Vec<i32> = phase_result.t2_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
                let ph_spe_alarm_int: Vec<i32> = phase_result.spe_alarm.iter().map(|&b| if b { 1 } else { 0 }).collect();
                let phase_list = list!(
                    t2 = Robj::from(phase_result.t2.clone()),
                    spe = Robj::from(phase_result.spe.clone()),
                    t2_alarm = Robj::from(ph_t2_alarm_int),
                    spe_alarm = Robj::from(ph_spe_alarm_int),
                    scores = fdmatrix_to_robj(&phase_result.scores)
                );
                result_list.set_attrib("phase", Robj::from(phase_list)).ok();
            }

            r!(result_list)
        }
        Err(e) => {
            rprintln!("elastic_spm_monitor failed: {:?}", e);
            r!(NULL)
        }
    }
}

// =============================================================================
// Module exports
// =============================================================================

extendr_module! {
    mod fdars;

    fn fdata_mean_1d;
    fn fdata_mean_2d;
    fn fdata_center_1d;
    fn fdata_norm_lp_1d;
    fn fdata_deriv_1d;
    fn fdata_deriv_2d;
    fn geometric_median_1d;
    fn geometric_median_2d;

    fn depth_fm_1d;
    fn depth_mode_1d;
    fn depth_rp_1d;
    fn depth_rt_1d;
    fn depth_fsd_1d;
    fn depth_kfsd_1d;
    fn depth_bd_1d;
    fn depth_mbd_1d;
    fn depth_mei_1d;

    fn depth_fm_2d;
    fn depth_mode_2d;
    fn depth_rp_2d;
    fn depth_rt_2d;
    fn depth_kfsd_2d;
    fn depth_fsd_2d;

    fn metric_lp_1d;
    fn metric_lp_self_1d;
    fn metric_hausdorff_1d;
    fn metric_hausdorff_cross_1d;
    fn metric_dtw_self_1d;
    fn metric_dtw_cross_1d;

    fn metric_lp_2d;
    fn metric_lp_self_2d;
    fn metric_hausdorff_2d;
    fn metric_hausdorff_cross_2d;

    fn compute_adot;
    fn pcvm_statistic;
    fn rp_stat;

    fn semimetric_fourier_self_1d;
    fn semimetric_fourier_cross_1d;
    fn semimetric_hshift_self_1d;
    fn semimetric_hshift_cross_1d;

    // Semimetric functions (PCA, deriv, basis)
    fn semimetric_pca_self_1d;
    fn semimetric_pca_cross_1d;
    fn semimetric_deriv_self_1d;
    fn semimetric_deriv_cross_1d;
    fn semimetric_basis_self_1d;
    fn semimetric_basis_cross_1d;

    fn metric_kl_self_1d;
    fn metric_kl_cross_1d;

    fn fdata2pc_1d;
    fn fdata2pls_1d;
    fn fdata2basis_1d;

    // Basis representation functions
    fn basis2fdata_1d;
    fn basis_gcv_1d;
    fn basis_aic_1d;
    fn basis_bic_1d;
    fn pspline_fit_1d;
    fn fdata2basis_2d;
    fn basis2fdata_2d;
    fn pspline_fit_2d;
    fn select_basis_auto;

    fn outliers_thres_lrt;
    fn outliers_lrt;

    // Smoothing functions
    fn s_nw;
    fn s_llr;
    fn s_lpr;
    fn s_knn;

    // Clustering functions
    fn kmeans_fd;
    fn fuzzycmeans_fd;

    // Registration functions
    fn register_shift_1d;

    // Utility functions
    fn int_simpson;
    fn inprod_fdata;

    // kNN regression functions
    fn knn_predict;
    fn knn_gcv;
    fn knn_lcv;

    // Cluster validation functions
    fn silhouette_score;
    fn calinski_harabasz;

    // Ridge regression (temporarily disabled - requires Rust 1.84+)
    // fn ridge_regression_fit;

    // Seasonal analysis functions
    fn seasonal_estimate_period_fft;
    fn seasonal_estimate_period_acf;
    fn seasonal_detect_multiple_periods;
    fn seasonal_detect_peaks;
    fn seasonal_strength_variance;
    fn seasonal_strength_spectral;
    fn seasonal_strength_wavelet;
    fn seasonal_strength_windowed;
    fn seasonal_detect_changes;
    fn seasonal_instantaneous_period;
    fn seasonal_analyze_peak_timing;
    fn seasonal_classify_seasonality;
    fn seasonal_detect_changes_auto;
    fn seasonal_detect_amplitude_modulation;
    fn seasonal_detect_amplitude_modulation_wavelet;
    fn seasonal_cfd_autoperiod;
    fn seasonal_autoperiod;
    fn seasonal_sazed;
    fn seasonal_lomb_scargle;
    fn seasonal_matrix_profile;
    fn seasonal_ssa;
    fn seasonal_stl;

    // Detrending functions
    fn seasonal_detrend;
    fn seasonal_decompose;

    // Simulation functions
    fn eigenfunctions_1d;
    fn eigenvalues_1d;
    fn sim_kl_1d;
    fn add_error_pointwise_1d;
    fn add_error_curve_1d;

    // Irregular functional data functions
    fn irreg_integrate;
    fn irreg_norm_lp;
    fn irreg_mean_kernel;
    fn irreg_metric_lp;
    fn irreg_to_regular;
    fn irreg_fdata2basis;

    // Alignment functions
    fn alignment_srsf_transform;
    fn alignment_srsf_inverse;
    fn alignment_reparameterize;
    fn alignment_compose_warps;
    fn alignment_elastic_pair;
    fn alignment_elastic_distance;
    fn alignment_align_to_target;
    fn alignment_self_dist;
    fn alignment_cross_dist;
    fn alignment_karcher_mean;

    // Soft-DTW functions
    fn metric_soft_dtw_self_1d;
    fn metric_soft_dtw_cross_1d;
    fn metric_soft_dtw_div_self_1d;
    fn metric_soft_dtw_div_cross_1d;
    fn metric_soft_dtw_barycenter;

    // Landmark registration functions
    fn landmark_detect;
    fn landmark_register_curves;

    // TSRVF functions
    fn alignment_tsrvf_transform;
    fn alignment_tsrvf_from_karcher;
    fn alignment_tsrvf_inverse;

    // Alignment quality functions
    fn alignment_quality_compute;
    fn alignment_warp_complexity;
    fn alignment_warp_smoothness;

    // Decomposition & distance functions
    fn alignment_decomposition;
    fn alignment_amplitude_dist;
    fn alignment_phase_dist;
    fn alignment_pairwise_consistency;

    // Constrained alignment functions
    fn alignment_constrained;
    fn alignment_with_landmarks;

    // Lambda CV, warp statistics, phase boxplots
    fn alignment_lambda_cv_rust;
    fn alignment_warp_statistics_rust;
    fn alignment_phase_boxplot_rust;

    // Warp inversion, penalized alignment, diagnostics
    fn alignment_invert_warp_rust;
    fn alignment_warp_inverse_error_rust;
    fn alignment_elastic_pair_penalized_rust;
    fn alignment_diagnose_rust;
    fn alignment_diagnose_pairwise_rust;

    // Shape analysis
    fn alignment_orbit_representative_rust;
    fn alignment_shape_distance_rust;
    fn alignment_shape_mean_rust;
    fn alignment_shape_self_distance_matrix_rust;

    // Elastic clustering
    fn elastic_kmeans_rust;
    fn elastic_hierarchical_rust;
    fn elastic_cut_dendrogram_rust;

    // Tolerance band functions
    fn tolerance_fpca;
    fn tolerance_conformal;
    fn tolerance_scb_degras;
    fn tolerance_exponential;
    fn tolerance_elastic;

    // Phase & elastic config tolerance bands
    fn tolerance_phase_rust;
    fn tolerance_elastic_config_rust;

    // Streaming depth functions
    fn streaming_depth_batch;
    fn streaming_depth_vs_ref;
    fn streaming_depth_one;

    // Scalar-on-function regression
    fn fregre_lm_rust;
    fn fregre_np_mixed_rust;
    fn functional_logistic_rust;
    fn fregre_cv_rust;
    fn bootstrap_ci_fregre_lm_rust;
    fn bootstrap_ci_functional_logistic_rust;

    // Robust regression
    fn fregre_l1_rust;
    fn fregre_huber_rust;
    fn predict_fregre_robust_rust;

    // Function-on-scalar regression
    fn fosr_rust;
    fn fosr_fpc_rust;
    fn fanova_rust;
    fn predict_fosr_rust;

    // Classification
    fn fclassif_lda_rust;
    fn fclassif_qda_rust;
    fn fclassif_knn_rust;
    fn fclassif_kernel_rust;
    fn fclassif_dd_rust;
    fn fclassif_cv_rust;

    // GMM clustering
    fn gmm_cluster_rust;
    fn gmm_em_rust;
    fn predict_gmm_rust;

    // Functional mixed models
    fn fmm_rust;
    fn fmm_predict_rust;
    fn fmm_test_fixed_rust;

    // Seeded depth variants
    fn depth_rp_1d_seeded;
    fn depth_rt_1d_seeded;

    // RPD depth
    fn depth_rpd_1d_rust;

    // Explainability wrappers
    fn explain_pdp_rust;
    fn explain_beta_decomposition_rust;
    fn explain_pointwise_importance_rust;
    fn explain_significant_regions_rust;
    fn explain_shap_rust;
    fn explain_influence_rust;
    fn explain_vif_rust;
    fn explain_dfbetas_rust;
    fn explain_saliency_rust;
    fn explain_importance_rust;
    fn explain_loo_rust;
    fn explain_sobol_rust;
    fn explain_ale_rust;
    fn explain_domain_rust;
    fn explain_prediction_intervals_rust;
    fn explain_lime_rust;
    fn explain_counterfactual_rust;
    fn explain_anchor_rust;
    fn explain_conformal_rust;
    fn explain_stability_rust;
    fn explain_friedman_rust;
    fn explain_regression_depth_rust;
    fn explain_prototype_rust;
    fn explain_calibration_rust;
    fn explain_ece_rust;
    fn explain_conditional_importance_rust;
    fn elastic_pcr_attribution_rust;

    // Logistic explainability wrappers
    fn explain_pdp_logistic_rust;
    fn explain_beta_decomposition_logistic_rust;
    fn explain_pointwise_importance_logistic_rust;
    fn explain_shap_logistic_rust;
    fn explain_vif_logistic_rust;
    fn explain_importance_logistic_rust;
    fn explain_sobol_logistic_rust;
    fn explain_ale_logistic_rust;
    fn explain_domain_logistic_rust;
    fn explain_saliency_logistic_rust;
    fn explain_lime_logistic_rust;
    fn explain_counterfactual_logistic_rust;
    fn explain_anchor_logistic_rust;
    fn explain_stability_logistic_rust;
    fn explain_friedman_logistic_rust;
    fn explain_regression_depth_logistic_rust;
    fn explain_conditional_importance_logistic_rust;

    // Elastic regression wrappers
    fn elastic_regression_rust;
    fn elastic_logistic_rust;
    fn elastic_pcr_rust;

    // Elastic FPCA wrappers
    fn elastic_vert_fpca_rust;
    fn elastic_horiz_fpca_rust;
    fn elastic_joint_fpca_rust;

    // Elastic changepoint wrappers
    fn elastic_amp_changepoint_rust;
    fn elastic_ph_changepoint_rust;
    fn elastic_fpca_changepoint_rust;

    // Smooth basis wrappers
    fn smooth_basis_rust;
    fn smooth_basis_gcv_rust;

    // Phase B1: 2D Function-on-Scalar
    fn fosr_2d_rust;
    fn predict_fosr_2d_rust;

    // Phase B2: Predict Methods
    fn predict_functional_logistic_rust;
    fn predict_elastic_regression_rust;
    fn predict_elastic_logistic_rust;

    // Phase B3: Model Selection
    fn model_selection_ncomp_rust;

    // Phase B4: Conformal Prediction
    fn conformal_fregre_lm_rust;
    fn conformal_fregre_np_rust;
    fn conformal_elastic_regression_rust;
    fn conformal_elastic_pcr_rust;
    fn conformal_classif_rust;
    fn conformal_logistic_rust;
    fn conformal_elastic_logistic_rust;

    // Phase C1: CV-Conformal & Jackknife+
    fn cv_conformal_fregre_lm_rust;
    fn cv_conformal_fregre_np_rust;
    fn cv_conformal_classif_rust;
    fn jackknife_plus_fregre_lm_rust;
    fn jackknife_plus_fregre_np_rust;

    // Phase C2: Conformal Generic
    fn conformal_generic_regression_lm_rust;
    fn conformal_generic_classification_logistic_rust;

    // Phase C3: Outliers with distribution
    fn outliers_thres_lrt_with_dist_rust;

    // Phase C4: Gradient
    fn fdata_gradient_rust;

    // Scalar-on-shape regression
    fn scalar_on_shape_rust;
    fn predict_scalar_on_shape_rust;

    // Statistical Process Monitoring (SPM)
    fn spm_phase1_rust;
    fn spm_monitor_rust;
    fn mfpca_rust;
    fn spm_ewma_rust;
    fn spm_t2_contrib_rust;
    fn spm_spe_contrib_rust;
    fn frcc_phase1_rust;
    fn frcc_monitor_rust;

    // SPM v0.9.0: ncomp, contrib, rules, bootstrap
    fn spm_select_ncomp_rust;
    fn spm_t2_pc_contrib_rust;
    fn spm_evaluate_rules_rust;
    fn spm_t2_limit_robust_rust;
    fn spm_spe_limit_robust_rust;

    // SPM v0.9.0: ARL
    fn spm_arl0_t2_rust;
    fn spm_arl1_t2_rust;
    fn spm_arl0_ewma_t2_rust;
    fn spm_arl0_spe_rust;

    // SPM v0.9.0: CUSUM
    fn spm_cusum_monitor_rust;
    fn spm_cusum_monitor_restart_rust;

    // SPM v0.9.0: MEWMA, AMEWMA, Iterative
    fn spm_mewma_monitor_rust;
    fn spm_amewma_monitor_rust;
    fn spm_phase1_iterative_rust;

    // SPM v0.9.0: Profile monitoring
    fn spm_profile_phase1_rust;
    fn spm_profile_monitor_rust;

    // SPM v0.9.0: Partial-domain monitoring
    fn spm_monitor_partial_rust;
    fn spm_monitor_partial_batch_rust;

    // SPM v0.9.0: Elastic SPM
    fn elastic_spm_phase1_rust;
    fn elastic_spm_monitor_rust;
}
