/**
 * @file sensor_calibration.c
 * @brief Sensor calibration algorithms — least-squares, polynomial fit, temperature comp
 *
 * Level: L5 — Algorithms/Methods
 */

#include "sensor_calibration.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ═══════════════════════════════════════════════════════════════════════
 * L5: Offset & Two-Point Calibration
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief One-point (offset-only) calibration.
 *
 * offset = mean(y_ref[i] - y_raw[i])
 *
 * This is valid when sensor sensitivity is known and stable,
 * or when operating near a single calibration point
 * (e.g., ice-point check for thermocouples).
 */
int calibration_offset(const double *x_raw_before, const double *x_ref,
                        size_t n_cal_points, double *offset)
{
    if (!x_raw_before || !x_ref || !offset || n_cal_points == 0) return -1;

    double sum = 0.0;
    for (size_t i = 0; i < n_cal_points; i++) {
        sum += x_ref[i] - x_raw_before[i];
    }
    *offset = sum / (double)n_cal_points;
    return 0;
}

int calibration_two_point(double x_raw_1, double x_ref_1,
                           double x_raw_2, double x_ref_2,
                           double *slope, double *intercept)
{
    if (!slope || !intercept) return -1;

    double dx_raw = x_raw_2 - x_raw_1;
    if (fabs(dx_raw) < 1e-30) return -1;

    *slope = (x_ref_2 - x_ref_1) / dx_raw;
    *intercept = x_ref_1 - *slope * x_raw_1;
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L5: Linear Least-Squares
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Linear least-squares with R² and RMS residual.
 *
 * Given n data pairs (x_i, y_i), the best-fit line y = m·x + b
 * minimizes Σ(y_i - ŷ_i)².
 *
 * Normal equations (Gauss, 1809):
 *   n·b + m·Σx = Σy
 *   b·Σx + m·Σx² = Σxy
 *
 * Solution:
 *   m = (n·Σxy - Σx·Σy) / (n·Σx² - (Σx)²)
 *   b = (Σy - m·Σx) / n
 *
 * R² (coefficient of determination):
 *   SS_tot = Σ(y_i - ȳ)²
 *   SS_res = Σ(y_i - ŷ_i)²
 *   R² = 1 - SS_res / SS_tot
 *
 * R² = 1 means perfect fit; R² = 0 means no better than horizontal line.
 *
 * Reference: Gauss, "Theoria motus corporum coelestium", 1809
 */
int calibration_linear_ls(const double *x_raw, const double *y_ref,
                           size_t n, double *slope, double *intercept,
                           double *r_squared, double *rms_residual)
{
    if (!x_raw || !y_ref || !slope || !intercept || n < 2) return -1;

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum_x  += x_raw[i];
        sum_y  += y_ref[i];
        sum_xy += x_raw[i] * y_ref[i];
        sum_x2 += x_raw[i] * x_raw[i];
    }

    double denom = n * sum_x2 - sum_x * sum_x;
    if (fabs(denom) < 1e-30) return -1;

    double m = (n * sum_xy - sum_x * sum_y) / denom;
    double b = (sum_y - m * sum_x) / n;

    *slope = m;
    *intercept = b;

    /* Compute R² and RMS residual */
    double y_mean = sum_y / n;
    double ss_tot = 0.0, ss_res = 0.0;
    for (size_t i = 0; i < n; i++) {
        double y_pred = m * x_raw[i] + b;
        double res = y_ref[i] - y_pred;
        ss_res += res * res;
        double dev = y_ref[i] - y_mean;
        ss_tot += dev * dev;
    }

    if (r_squared) {
        *r_squared = (ss_tot > 1e-30) ? (1.0 - ss_res / ss_tot) : 1.0;
    }
    if (rms_residual) {
        *rms_residual = sqrt(ss_res / n);
    }

    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L5: Polynomial Least-Squares
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Solve a small linear system A·x = b using Gaussian elimination
 *        with partial pivoting (Cholesky not strictly needed for small systems).
 *
 * This is a utility for polynomial least-squares.
 */
static int solve_linear_system(double *A, double *b, int n)
{
    /* Gaussian elimination with partial pivoting */
    for (int col = 0; col < n; col++) {
        /* Find pivot */
        int pivot_row = col;
        double pivot_val = fabs(A[col * n + col]);
        for (int row = col + 1; row < n; row++) {
            double val = fabs(A[row * n + col]);
            if (val > pivot_val) {
                pivot_val = val;
                pivot_row = row;
            }
        }
        if (pivot_val < 1e-30) return -1;

        /* Swap rows if needed */
        if (pivot_row != col) {
            for (int j = 0; j < n; j++) {
                double tmp = A[col * n + j];
                A[col * n + j] = A[pivot_row * n + j];
                A[pivot_row * n + j] = tmp;
            }
            double tmp = b[col];
            b[col] = b[pivot_row];
            b[pivot_row] = tmp;
        }

        /* Eliminate below */
        for (int row = col + 1; row < n; row++) {
            double factor = A[row * n + col] / A[col * n + col];
            for (int j = col; j < n; j++) {
                A[row * n + j] -= factor * A[col * n + j];
            }
            b[row] -= factor * b[col];
        }
    }

    /* Back substitution */
    for (int i = n - 1; i >= 0; i--) {
        double sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= A[i * n + j] * b[j];
        }
        b[i] = sum / A[i * n + i];
    }
    return 0;
}

/**
 * @brief Polynomial least-squares regression.
 *
 * Fits y = a₀ + a₁·x + ... + a_order·x^order.
 *
 * Method: Normal equations (X^T·X)·a = X^T·y
 *
 * Vandermonde matrix V[i][j] = x_i^j for j = 0..order
 *
 * (V^T·V)_{jk} = Σ_i x_i^{j+k}
 * (V^T·y)_j = Σ_i x_i^j · y_i
 *
 * Complexity: O(n·p² + p³) where p = order + 1
 *
 * Condition note: For high-order (>5) or poorly scaled x data,
 * the normal equations become ill-conditioned. QR decomposition
 * (Householder/Givens) is preferred for order > 7.
 */
int calibration_polynomial_ls(const double *x, const double *y,
                               size_t n, int order, double *coeff,
                               double *rms_residual)
{
    if (!x || !y || !coeff || n < (size_t)(order + 2) || order < 1 || order > 9)
        return -1;

    int p = order + 1;

    /* Build normal equations: A_{jk} = Σ x_i^{j+k},  b_j = Σ x_i^j · y_i */
    double *A = (double *)calloc(p * p, sizeof(double));
    double *b = (double *)calloc(p, sizeof(double));
    if (!A || !b) {
        free(A); free(b);
        return -1;
    }

    for (size_t i = 0; i < n; i++) {
        double x_pow = 1.0;
        for (int j = 0; j < p; j++) {
            b[j] += x_pow * y[i];
            double x_pow_k = 1.0;
            for (int k = 0; k < p; k++) {
                A[j * p + k] += x_pow * x_pow_k;
                x_pow_k *= x[i];
            }
            x_pow *= x[i];
        }
    }

    /* Solve */
    int ret = solve_linear_system(A, b, p);

    /* Copy coefficients: b now contains the solution */
    for (int j = 0; j < p; j++) {
        coeff[j] = b[j];
    }

    /* Compute RMS residual */
    if (rms_residual && ret == 0) {
        double ss_res = 0.0;
        for (size_t i = 0; i < n; i++) {
            double y_pred = 0.0;
            double x_pow = 1.0;
            for (int j = 0; j < p; j++) {
                y_pred += coeff[j] * x_pow;
                x_pow *= x[i];
            }
            double res = y[i] - y_pred;
            ss_res += res * res;
        }
        *rms_residual = sqrt(ss_res / n);
    }

    free(A);
    free(b);
    return ret;
}

/**
 * @brief Compute fit quality metrics for a polynomial calibration.
 *
 * R² = 1 - SS_res / SS_tot
 * Adjusted R² = 1 - (1-R²)·(n-1)/(n-p-1)
 *
 * The adjusted R² penalizes adding unnecessary polynomial terms,
 * preventing overfitting. For a good calibration, adj R² should be
 * close to R², indicating that the extra terms are justified.
 */
int calibration_fit_metrics(const double *x, const double *y,
                             size_t n, const double *coeff, int order,
                             double *r_squared, double *adj_r_squared,
                             double *rms_residual, double *max_residual)
{
    if (!x || !y || !coeff || n < 2 || order < 0) return -1;

    int p = order + 1;
    double y_mean = 0.0;
    for (size_t i = 0; i < n; i++) y_mean += y[i];
    y_mean /= n;

    double ss_tot = 0.0, ss_res = 0.0;
    double max_res = 0.0;

    for (size_t i = 0; i < n; i++) {
        double y_pred = 0.0;
        double x_pow = 1.0;
        for (int j = 0; j < p; j++) {
            y_pred += coeff[j] * x_pow;
            x_pow *= x[i];
        }
        double res = y[i] - y_pred;
        double ares = fabs(res);
        if (ares > max_res) max_res = ares;
        ss_res += res * res;
        double dev = y[i] - y_mean;
        ss_tot += dev * dev;
    }

    if (r_squared) {
        *r_squared = (ss_tot > 1e-30) ? (1.0 - ss_res / ss_tot) : 1.0;
    }
    if (adj_r_squared) {
        double R2 = (ss_tot > 1e-30) ? (1.0 - ss_res / ss_tot) : 1.0;
        double n_minus_p = (double)n - (double)p - 1.0;
        if (n_minus_p > 0.0) {
            *adj_r_squared = 1.0 - (1.0 - R2) * (n - 1.0) / n_minus_p;
        } else {
            *adj_r_squared = R2;
        }
    }
    if (rms_residual) {
        *rms_residual = sqrt(ss_res / n);
    }
    if (max_residual) {
        *max_residual = max_res;
    }

    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L5: Weighted Linear Least-Squares
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Weighted least-squares for linear calibration.
 *
 * Each data point has weight w_i (typically w_i = 1/σ_i²).
 *
 * This ensures that reference points with smaller uncertainty
 * have greater influence on the calibration coefficients.
 *
 * This is critical when combining calibration data from multiple
 * standards with different accuracies (e.g., NIST traceable at some
 * points, in-house reference at others).
 */
int calibration_weighted_linear_ls(const double *x, const double *y,
                                    const double *weights, size_t n,
                                    double *slope, double *intercept,
                                    double *rms_residual)
{
    if (!x || !y || !weights || !slope || !intercept || n < 2) return -1;

    double sum_w = 0.0, sum_wx = 0.0, sum_wy = 0.0;
    double sum_wxy = 0.0, sum_wx2 = 0.0;

    for (size_t i = 0; i < n; i++) {
        double wi = weights[i];
        if (wi < 0.0) return -1;
        sum_w   += wi;
        sum_wx  += wi * x[i];
        sum_wy  += wi * y[i];
        sum_wxy += wi * x[i] * y[i];
        sum_wx2 += wi * x[i] * x[i];
    }

    double denom = sum_w * sum_wx2 - sum_wx * sum_wx;
    if (fabs(denom) < 1e-30) return -1;

    double m = (sum_w * sum_wxy - sum_wx * sum_wy) / denom;
    double b = (sum_wy - m * sum_wx) / sum_w;

    *slope = m;
    *intercept = b;

    /* Weighted RMS residual */
    if (rms_residual) {
        double wss_res = 0.0;
        for (size_t i = 0; i < n; i++) {
            double y_pred = m * x[i] + b;
            double res = y[i] - y_pred;
            wss_res += weights[i] * res * res;
        }
        *rms_residual = sqrt(wss_res / sum_w);
    }

    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L5: Temperature Compensation
 * ═══════════════════════════════════════════════════════════════════════ */

int calibration_temp_compensation(const double *temps,
                                   const double *slopes,
                                   const double *offsets,
                                   size_t n_temps, double T0,
                                   double *S0_out, double *alpha_out,
                                   double *b0_out, double *beta_out)
{
    if (!temps || !slopes || !offsets || n_temps < 2) return -1;

    /* Fit: slope vs temperature → S(T) = S₀ + α·(T-T₀)
     *      offset vs temperature → b(T) = b₀ + β·(T-T₀)
     *
     * Use linear least-squares separately. */

    /* Prepare shifted temperature: t = T - T₀ */
    double sum_t = 0.0, sum_t2 = 0.0;
    double sum_s = 0.0, sum_st = 0.0;
    double sum_b = 0.0, sum_bt = 0.0;

    for (size_t i = 0; i < n_temps; i++) {
        double dt = temps[i] - T0;
        sum_t  += dt;
        sum_t2 += dt * dt;
        sum_s  += slopes[i];
        sum_st += dt * slopes[i];
        sum_b  += offsets[i];
        sum_bt += dt * offsets[i];
    }

    double denom = n_temps * sum_t2 - sum_t * sum_t;
    double alpha_val = 0.0, S0_val = 0.0, beta_val = 0.0, b0_val = 0.0;

    if (fabs(denom) > 1e-30) {
        /* α = (n·ΣtS - Σt·ΣS) / (n·Σt² - (Σt)²) */
        alpha_val = (n_temps * sum_st - sum_t * sum_s) / denom;
        S0_val = (sum_s - alpha_val * sum_t) / n_temps;

        beta_val = (n_temps * sum_bt - sum_t * sum_b) / denom;
        b0_val = (sum_b - beta_val * sum_t) / n_temps;
    } else {
        /* All temperatures are the same (degenerate), use mean */
        S0_val = sum_s / n_temps;
        b0_val = sum_b / n_temps;
        alpha_val = 0.0;
        beta_val = 0.0;
    }

    if (S0_out) *S0_out = S0_val;
    if (alpha_out) *alpha_out = alpha_val;
    if (b0_out) *b0_out = b0_val;
    if (beta_out) *beta_out = beta_val;

    return 0;
}

/**
 * @brief Apply temperature compensation to a raw reading.
 *
 * y_comp = (y_raw - b_eff) / S_eff
 *
 * where:
 *   b_eff = b₀ + β·(T_actual - T₀)  [effective offset at T_actual]
 *   S_eff = S₀ + α·(T_actual - T₀)  [effective sensitivity at T_actual]
 *
 * Example: Pressure sensor with α = 500 ppm/°C, T₀ = 25°C,
 *          T_actual = 50°C → S_eff = S₀ · (1 + 500e-6 · 25) = 1.0125·S₀
 *          → ~1.25% sensitivity change — significant for precision!
 */
double calibration_apply_temp_comp(double y_raw, double T_actual,
                                    double T0, double S0, double alpha,
                                    double b0, double beta)
{
    double dT = T_actual - T0;
    double S_eff = S0 + alpha * dT;
    double b_eff = b0 + beta * dT;

    if (fabs(S_eff) < 1e-30) return y_raw;
    return (y_raw - b_eff) / S_eff;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L5: Cross-Sensitivity Compensation
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Single-parameter cross-sensitivity compensation.
 *
 * When a sensor responds to both the measurand x and an interference
 * quantity z, the reading must be corrected.
 *
 * y = S_x·x + S_z·z + b
 * → x = (y - b - S_z·z) / S_x
 *
 * Common examples:
 *   - Humidity sensor (also sensitive to temperature):
 *     RH_comp = (RH_raw - b - S_T·T) / S_RH
 *   - O₂ sensor (also sensitive to pressure):
 *     O₂_comp = (O₂_raw - b - S_P·P) / S_O2
 *   - Strain gauge (also sensitive to temperature: apparent strain)
 */
double calibration_cross_sensitivity_comp(double y_raw, double z_measured,
                                           double S_x, double S_z, double b)
{
    if (fabs(S_x) < 1e-30) return y_raw;
    return (y_raw - b - S_z * z_measured) / S_x;
}

/**
 * @brief Multi-sensor cross-sensitivity matrix compensation.
 *
 * For N sensors measuring N quantities with cross-sensitivity:
 *   y = S·x + b
 *
 * Then: x = S⁻¹·(y - b)
 *
 * This is a small linear system solve (N ≤ 5 typically).
 * Uses Gaussian elimination for the inverse.
 *
 * Example: 3-axis accelerometer with cross-axis sensitivity:
 *   S = [[S_xx, S_xy, S_xz],    (actual: diagonal dominate, off-diagonal ~1-5%)
 *        [S_yx, S_yy, S_yz],
 *        [S_zx, S_zy, S_zz]]
 */
int calibration_cross_sensitivity_matrix(const double *y, const double *S,
                                          const double *b, int n,
                                          double *x_compensated)
{
    if (!y || !S || !b || !x_compensated || n < 1 || n > 8) return -1;

    /* Solve S·x = y - b using Gaussian elimination on the augmented matrix */
    double *aug = (double *)calloc(n * (n + 1), sizeof(double));
    if (!aug) return -1;

    /* Copy S and RHS into augmented matrix */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i * (n + 1) + j] = S[i * n + j];
        }
        aug[i * (n + 1) + n] = y[i] - b[i];
    }

    /* Forward elimination */
    for (int col = 0; col < n; col++) {
        /* Find pivot */
        int pivot = col;
        double max_val = fabs(aug[col * (n + 1) + col]);
        for (int row = col + 1; row < n; row++) {
            double val = fabs(aug[row * (n + 1) + col]);
            if (val > max_val) {
                max_val = val;
                pivot = row;
            }
        }
        if (max_val < 1e-30) { free(aug); return -1; }

        /* Swap pivot row to current row */
        if (pivot != col) {
            for (int j = 0; j <= n; j++) {
                double tmp = aug[col * (n + 1) + j];
                aug[col * (n + 1) + j] = aug[pivot * (n + 1) + j];
                aug[pivot * (n + 1) + j] = tmp;
            }
        }

        /* Eliminate below and above (Gauss-Jordan for simplicity) */
        double pivot_val = aug[col * (n + 1) + col];
        /* Normalize pivot row */
        for (int j = 0; j <= n; j++) {
            aug[col * (n + 1) + j] /= pivot_val;
        }

        /* Eliminate all other rows */
        for (int row = 0; row < n; row++) {
            if (row == col) continue;
            double factor = aug[row * (n + 1) + col];
            for (int j = 0; j <= n; j++) {
                aug[row * (n + 1) + j] -= factor * aug[col * (n + 1) + j];
            }
        }
    }

    /* Extract solution */
    for (int i = 0; i < n; i++) {
        x_compensated[i] = aug[i * (n + 1) + n];
    }

    free(aug);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L5: Auto-Zero & Piecewise Linearization
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Auto-zero (correlated double sampling) simulation.
 *
 * In each cycle:
 *   Phase 1: Measure zero reference → z[k]
 *   Phase 2: Measure signal → s[k]
 *   Corrected: y[k] = s[k] - z[k]
 *
 * This technique is fundamental to precision measurement:
 *   - Chopper-stabilized op-amps use this to cancel offset and 1/f noise
 *   - Correlated double sampling (CDS) in CCD/CMOS image sensors
 *   - AC bridge excitation with synchronous demodulation
 *
 * Effective noise reduction for 1/f noise: noise reduced by factor
 * proportional to (f_chopper / f_signal) for frequencies below f_chopper/2.
 */
int calibration_auto_zero(const double *signal_with_offset,
                           const double *zero_readings,
                           size_t n_cycles, double *corrected_signal)
{
    if (!signal_with_offset || !zero_readings || !corrected_signal || n_cycles == 0)
        return -1;

    for (size_t i = 0; i < n_cycles; i++) {
        corrected_signal[i] = signal_with_offset[i] - zero_readings[i];
    }
    return 0;
}

int calibration_piecewise_3point(double x0, double y0,
                                  double x1, double y1,
                                  double x2, double y2,
                                  double *m1, double *b1,
                                  double *m2, double *b2)
{
    if (!m1 || !b1 || !m2 || !b2) return -1;

    double dx1 = x1 - x0;
    double dx2 = x2 - x1;

    if (fabs(dx1) < 1e-30 || fabs(dx2) < 1e-30) return -1;

    *m1 = (y1 - y0) / dx1;
    *b1 = y0 - *m1 * x0;

    *m2 = (y2 - y1) / dx2;
    *b2 = y1 - *m2 * x1;

    return 0;
}
