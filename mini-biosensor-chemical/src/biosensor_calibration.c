/**
 * biosensor_calibration.c — Calibration methods and data reduction
 *
 * L5: Linear regression (OLS, WLS), 4PL/5PL Levenberg-Marquardt fitting,
 *     standard addition method, Grubbs & Dixon outlier tests,
 *     CV% computation, Mandel lack-of-fit test, Durbin-Watson test.
 * L6: Two-point recalibration, sensor lifetime estimation.
 * L8: Piecewise direct standardization for calibration transfer.
 */

#include "biosensor_calibration.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * L5: Linear regression calibration
 * ============================================================================ */

/**
 * Ordinary least squares: y = a + b·x
 *
 * b = S_xy / S_xx = Σ(x_i - x̄)(y_i - ȳ) / Σ(x_i - x̄)²
 * a = ȳ - b·x̄
 *
 * Standard errors:
 *   s_yx = √(Σ(y_i - ŷ_i)² / (N - 2))
 *   SE(b) = s_yx / √S_xx
 *   SE(a) = s_yx · √(1/N + x̄²/S_xx)
 */
double linear_calibration_fit(const CalibrationStandards *standards,
                              LinearCalibration *result)
{
    if (!standards || !result || standards->standard_count < 2) {
        if (result) linear_calibration_init(result);
        return 0.0;
    }

    int n = standards->standard_count;
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0, sum_yy = 0.0;

    for (int i = 0; i < n; i++) {
        double x = standards->concentrations[i];
        double y = standards->signals[i];
        sum_x  += x;
        sum_y  += y;
        sum_xy += x * y;
        sum_xx += x * x;
        sum_yy += y * y;
    }

    double N = (double)n;
    double denom = N * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1.0e-30) {
        linear_calibration_init(result);
        return 0.0;
    }

    double slope     = (N * sum_xy - sum_x * sum_y) / denom;
    double intercept = (sum_y - slope * sum_x) / N;

    /* Residuals and R² */
    double ss_res = 0.0, ss_tot = 0.0;
    double y_mean = sum_y / N;
    for (int i = 0; i < n; i++) {
        double pred = slope * standards->concentrations[i] + intercept;
        double res  = standards->signals[i] - pred;
        ss_res += res * res;
        ss_tot += (standards->signals[i] - y_mean) *
                  (standards->signals[i] - y_mean);
    }
    double r_sq = (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;

    /* Standard errors */
    double s_yx = sqrt(ss_res / (N - 2.0));
    double s_xx = sum_xx - sum_x * sum_x / N;
    double x_mean = sum_x / N;

    result->slope     = slope;
    result->intercept = intercept;
    result->r_squared = r_sq;
    result->standard_error_slope     = s_yx / sqrt(s_xx);
    result->standard_error_intercept = s_yx * sqrt(1.0 / N + x_mean * x_mean / s_xx);
    result->prediction_interval_95   = 2.0 * s_yx;  /* approximate 95% PI */

    return r_sq;
}

/**
 * Weighted least squares: minimize Σ w_i·(y_i - (a + b·x_i))²
 *
 * w_i = 1/σ_i²
 *
 * The solution is the normal equations with weights:
 *   b = (Σw · Σwxy - Σwx · Σwy) / (Σw · Σwxx - (Σwx)²)
 *   a = (Σwy - b·Σwx) / Σw
 */
double weighted_linear_calibration_fit(const CalibrationStandards *standards,
                                       LinearCalibration *result)
{
    if (!standards || !result || standards->standard_count < 2) {
        if (result) linear_calibration_init(result);
        return 0.0;
    }

    int n = standards->standard_count;

    /* Check for standard deviations; fall back to OLS if missing */
    bool has_stds = true;
    for (int i = 0; i < n; i++) {
        if (standards->signal_stds[i] <= 0.0) { has_stds = false; break; }
    }
    if (!has_stds) {
        return linear_calibration_fit(standards, result);
    }

    double sum_w = 0.0, sum_wx = 0.0, sum_wy = 0.0;
    double sum_wxx = 0.0, sum_wxy = 0.0, sum_wyy = 0.0;

    for (int i = 0; i < n; i++) {
        double w = 1.0 / (standards->signal_stds[i] *
                          standards->signal_stds[i]);
        double x = standards->concentrations[i];
        double y = standards->signals[i];

        sum_w   += w;
        sum_wx  += w * x;
        sum_wy  += w * y;
        sum_wxx += w * x * x;
        sum_wxy += w * x * y;
        sum_wyy += w * y * y;
    }

    double denom = sum_w * sum_wxx - sum_wx * sum_wx;
    if (fabs(denom) < 1.0e-30) {
        linear_calibration_init(result);
        return 0.0;
    }

    double slope     = (sum_w * sum_wxy - sum_wx * sum_wy) / denom;
    double intercept = (sum_wy - slope * sum_wx) / sum_w;

    /* Weighted R² */
    double ss_res = 0.0, ss_tot = 0.0;
    double y_wmean = sum_wy / sum_w;
    for (int i = 0; i < n; i++) {
        double w = 1.0 / (standards->signal_stds[i] *
                          standards->signal_stds[i]);
        double pred = slope * standards->concentrations[i] + intercept;
        double res  = standards->signals[i] - pred;
        ss_res += w * res * res;
        ss_tot += w * (standards->signals[i] - y_wmean) *
                     (standards->signals[i] - y_wmean);
    }
    double r_sq = (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;

    result->slope     = slope;
    result->intercept = intercept;
    result->r_squared = r_sq;
    result->standard_error_slope     = sqrt(sum_w / denom);
    result->standard_error_intercept = sqrt(sum_wxx / denom);
    result->prediction_interval_95   = 2.0 * sqrt(ss_res / (sum_w - 2.0));

    return r_sq;
}

/**
 * Predict concentration from linear calibration with prediction interval.
 *
 * x_pred = (y - a) / b
 *
 * Prediction interval accounts for:
 * 1. Uncertainty in a and b (calibration uncertainty)
 * 2. Measurement uncertainty of y (sample signal)
 */
void linear_calibration_predict(const LinearCalibration *lc,
                                double signal,
                                const CalibrationStandards *standards,
                                double *concentration,
                                double *lower_bound,
                                double *upper_bound)
{
    if (!lc || !concentration || !lower_bound || !upper_bound) return;

    if (fabs(lc->slope) < 1.0e-15) {
        *concentration = 0.0;
        *lower_bound   = 0.0;
        *upper_bound   = 0.0;
        return;
    }

    *concentration = (signal - lc->intercept) / lc->slope;
    double half_interval = lc->prediction_interval_95;

    if (standards && standards->standard_count > 2) {
        /* Compute leverage for this prediction */
        int n = standards->standard_count;
        double sum_x = 0.0, sum_xx = 0.0;
        for (int i = 0; i < n; i++) {
            sum_x  += standards->concentrations[i];
            sum_xx += standards->concentrations[i] *
                      standards->concentrations[i];
        }
        double x_mean = sum_x / n;
        double s_xx = sum_xx - sum_x * sum_x / n;
        if (s_xx > 0.0) {
            double leverage = 1.0 + 1.0 / n +
                (*concentration - x_mean) * (*concentration - x_mean) / s_xx;
            half_interval = lc->prediction_interval_95 * sqrt(leverage);
        }
    }

    *lower_bound = *concentration - half_interval;
    *upper_bound = *concentration + half_interval;
}

/* ============================================================================
 * L5: 4PL and 5PL Logistic Fitting
 * ============================================================================ */

/**
 * 4PL model: y = a + (d - a) / (1 + (x/c)^b)
 *
 * Jacobian for Levenberg-Marquardt:
 *   dy/da = 1 - 1/(1+(x/c)^b)
 *   dy/db = -(d-a)·(x/c)^b·ln(x/c) / (1+(x/c)^b)²
 *   dy/dc = (d-a)·b·(x/c)^b / (c·(1+(x/c)^b)²)
 *   dy/dd = 1/(1+(x/c)^b)
 */
int fourpl_fit(const CalibrationStandards *standards,
               FourPLLogisticFit *fit,
               int max_iter,
               double lambda_init)
{
    if (!standards || !fit || standards->standard_count < 4 ||
        max_iter < 1) return 0;

    int n = standards->standard_count;
    double a = fit->a, b = fit->b, c = fit->c, d = fit->d;
    double lambda = lambda_init;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Build Jacobian J and residual vector r */
        double jtj[4][4] = {{0}};
        double jtr[4] = {0};
        double ss_res = 0.0;

        for (int i = 0; i < n; i++) {
            double x = standards->concentrations[i];
            double y_meas = standards->signals[i];
            if (x <= 0.0) continue;

            double ratio = x / c;
            double ratio_b = pow(ratio, b);
            double denom = 1.0 + ratio_b;
            double denom_sq = denom * denom;

            double y_pred = a + (d - a) / denom;
            double residual = y_meas - y_pred;

            /* Partial derivatives */
            double da = 1.0 - 1.0 / denom;
            double db = -(d - a) * ratio_b * log(ratio) / denom_sq;
            double dc = (d - a) * b * ratio_b / (c * denom_sq);
            double dd = 1.0 / denom;

            double j_row[4] = {da, db, dc, dd};

            /* Accumulate J^T J and J^T r */
            for (int r = 0; r < 4; r++) {
                for (int c = 0; c < 4; c++) {
                    jtj[r][c] += j_row[r] * j_row[c];
                }
                jtr[r] += j_row[r] * residual;
            }
            ss_res += residual * residual;
        }

        /* Levenberg-Marquardt: (J^T·J + λ·diag(J^T·J)) · Δθ = J^T·r */
        double jtj_damped[4][4];
        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 4; c++) {
                jtj_damped[r][c] = jtj[r][c];
            }
            jtj_damped[r][r] += lambda * jtj[r][r];
            if (fabs(jtj_damped[r][r]) < 1.0e-15) {
                jtj_damped[r][r] = 1.0e-15;
            }
        }

        /* Solve 4×4 using Gaussian elimination */
        double delta[4];
        /* Copy augmented matrix */
        double aug[4][5];
        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 4; c++) {
                aug[r][c] = jtj_damped[r][c];
            }
            aug[r][4] = jtr[r];
        }

        /* Forward elimination */
        for (int p = 0; p < 4; p++) {
            /* Find pivot */
            int max_row = p;
            double max_val = fabs(aug[p][p]);
            for (int r = p + 1; r < 4; r++) {
                if (fabs(aug[r][p]) > max_val) {
                    max_val = fabs(aug[r][p]);
                    max_row = r;
                }
            }
            if (max_row != p) {
                for (int c = 0; c < 5; c++) {
                    double tmp = aug[p][c];
                    aug[p][c] = aug[max_row][c];
                    aug[max_row][c] = tmp;
                }
            }

            if (fabs(aug[p][p]) < 1.0e-15) continue;

            for (int r = p + 1; r < 4; r++) {
                double factor = aug[r][p] / aug[p][p];
                for (int c = p; c < 5; c++) {
                    aug[r][c] -= factor * aug[p][c];
                }
            }
        }

        /* Back substitution */
        for (int r = 3; r >= 0; r--) {
            double sum = aug[r][4];
            for (int c = r + 1; c < 4; c++) {
                sum -= aug[r][c] * delta[c];
            }
            delta[r] = (fabs(aug[r][r]) > 1.0e-15) ? sum / aug[r][r] : 0.0;
        }

        /* Update parameters */
        double new_a = a + delta[0];
        double new_b = b + delta[1];
        double new_c = c + delta[2];
        double new_d = d + delta[3];

        /* Compute new residuals */
        double new_ss_res = 0.0;
        for (int i = 0; i < n; i++) {
            double x = standards->concentrations[i];
            if (x <= 0.0) continue;
            double ratio = x / new_c;
            double ratio_b = pow(ratio, new_b);
            double y_pred = new_a + (new_d - new_a) / (1.0 + ratio_b);
            double res = standards->signals[i] - y_pred;
            new_ss_res += res * res;
        }

        /* Accept or adjust lambda */
        if (new_ss_res < ss_res && new_c > 0.0 && new_b > 0.0) {
            a = new_a; b = new_b; c = new_c; d = new_d;
            lambda *= 0.1;  /* reduce damping */
        } else {
            lambda *= 10.0; /* increase damping */
        }

        /* Convergence check */
        double rel_change = fabs(delta[0]/a) + fabs(delta[1]/b) +
                            fabs(delta[2]/c) + fabs(delta[3]/d);
        if (rel_change < 1.0e-6) {
            fit->a = a; fit->b = b; fit->c = c; fit->d = d;
            fit->r_squared = 1.0 - ss_res / (n * (n > 1 ? 1.0 : 1.0));
            fit->rmse = sqrt(ss_res / n);
            return iter + 1;
        }
    }

    fit->a = a; fit->b = b; fit->c = c; fit->d = d;
    fit->r_squared = 0.0; /* approximate */
    fit->rmse = 0.0;
    return max_iter;
}

/**
 * 5PL model: y = a + (d - a) / (1 + (x/c)^b)^g
 *
 * The additional parameter g models asymmetry.
 */
int fivepl_fit(const CalibrationStandards *standards,
               FivePLLogisticFit *fit,
               int max_iter,
               double lambda_init)
{
    if (!standards || !fit || standards->standard_count < 5 ||
        max_iter < 1) return 0;

    /* For brevity and numerical stability, first fit a 4PL,
     * then treat g=1.0 as initial and perform full 5D LM.
     * Here we implement a simplified approach:
     * fit 4PL, then optimize g via grid search. */

    /* Step 1: Fit 4PL to get initial a, b, c, d */
    FourPLLogisticFit f4;
    f4.a = fit->a; f4.b = fit->b; f4.c = fit->c; f4.d = fit->d;
    int fourpl_iters = fourpl_fit(standards, &f4, max_iter, lambda_init);
    if (fourpl_iters == 0) return 0;

    /* Step 2: Grid search for optimal g */
    int n = standards->standard_count;
    double best_g = 1.0;
    double best_ss = INFINITY;

    for (int gi = 0; gi < 100; gi++) {
        double g_test = 0.1 + 0.05 * gi;  /* 0.1 to 5.05 */
        double ss = 0.0;
        for (int i = 0; i < n; i++) {
            double x = standards->concentrations[i];
            if (x <= 0.0) continue;
            double ratio = x / f4.c;
            double y_pred = f4.a + (f4.d - f4.a) /
                pow(1.0 + pow(ratio, f4.b), g_test);
            double res = standards->signals[i] - y_pred;
            ss += res * res;
        }
        if (ss < best_ss) { best_ss = ss; best_g = g_test; }
    }

    fit->a = f4.a; fit->b = f4.b; fit->c = f4.c; fit->d = f4.d;
    fit->g = best_g;
    fit->r_squared = 1.0 - best_ss / (n * (n > 1 ? 1.0 : 1.0));
    fit->rmse = sqrt(best_ss / n);
    return fourpl_iters;
}

/**
 * Log-logistic linearization for initial 4PL parameter guesses.
 *
 * log((d - y)/(y - a)) = b · log(x) - b · log(c)
 *
 * Using known or estimated asymptotes a, d, a linear regression
 * of log((d-y)/(y-a)) vs log(x) yields b (slope) and log(c) (intercept).
 */
double fourpl_linearized_initial_guess(const CalibrationStandards *standards,
                                       FourPLLogisticFit *fit)
{
    if (!standards || !fit || standards->standard_count < 4) return 0.0;

    int n = standards->standard_count;

    /* Use extremes as asymptote estimates if not already set */
    if (fit->a <= 0.0 || fit->d <= fit->a) {
        fit->a = standards->signals[0];       /* lowest standard signal */
        fit->d = standards->signals[n - 1];   /* highest standard signal */
    }

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0, sum_yy = 0.0;
    int valid = 0;

    for (int i = 0; i < n; i++) {
        double signal = standards->signals[i];
        /* Avoid log of non-positive values */
        if (signal <= fit->a || signal >= fit->d) continue;

        double x = log(standards->concentrations[i]);
        double y = log((fit->d - signal) / (signal - fit->a));

        sum_x  += x;
        sum_y  += y;
        sum_xy += x * y;
        sum_xx += x * x;
        sum_yy += y * y;
        valid++;
    }

    if (valid < 2) return 0.0;

    double N = (double)valid;
    double denom = N * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1.0e-30) return 0.0;

    double slope     = (N * sum_xy - sum_x * sum_y) / denom;  /* b */
    double intercept = (sum_y - slope * sum_x) / N;           /* -b·log(c) */

    fit->b = slope;
    fit->c = exp(-intercept / slope);  /* c = exp(-intercept/b) */

    double y_mean = sum_y / N;
    double ss_res = 0.0, ss_tot = 0.0;
    for (int i = 0; i < n; i++) {
        double signal = standards->signals[i];
        if (signal <= fit->a || signal >= fit->d) continue;
        double y_pred = slope * log(standards->concentrations[i]) + intercept;
        double y_meas = log((fit->d - signal) / (signal - fit->a));
        ss_res += (y_meas - y_pred) * (y_meas - y_pred);
        ss_tot += (y_meas - y_mean) * (y_meas - y_mean);
    }
    return (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;
}

/* ============================================================================
 * L5: Standard addition method
 * ============================================================================ */

/**
 * Standard addition for complex matrices.
 *
 * In each aliquot, the sample is spiked with a known addition of analyte.
 * The total concentration in aliquot i is:
 *   C_i = (C_orig · V_sample + C_spike_i · V_spike_i) / V_total_i
 *
 * Linear regression of signal vs added concentration gives:
 *   x-intercept = -C_orig (when signal = 0, extrapolated)
 *
 * Alternatively: C_orig = |intercept / slope|
 */
bool standard_addition_method(const double *sample_volume_ml,
                              const double *spike_concentrations,
                              const double *spike_volumes_ml,
                              const double *total_volumes_ml,
                              const double *measured_signals,
                              int n_aliquots,
                              double *original_conc,
                              double *r_squared)
{
    if (!sample_volume_ml || !spike_concentrations ||
        !spike_volumes_ml || !total_volumes_ml ||
        !measured_signals || !original_conc || !r_squared ||
        n_aliquots < 3) return false;

    /* Compute added concentrations (the known x-axis values) */
    double *added_conc = (double *)malloc((size_t)n_aliquots * sizeof(double));
    if (!added_conc) return false;

    for (int i = 0; i < n_aliquots; i++) {
        if (sample_volume_ml[i] <= 0.0 || total_volumes_ml[i] <= 0.0) {
            free(added_conc);
            return false;
        }
        added_conc[i] = spike_concentrations[i] * spike_volumes_ml[i] /
                        total_volumes_ml[i];
    }

    /* Linear regression: signal = slope · added_conc + intercept */
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0, sum_yy = 0.0;
    for (int i = 0; i < n_aliquots; i++) {
        double x = added_conc[i];
        double y = measured_signals[i];
        sum_x  += x;
        sum_y  += y;
        sum_xy += x * y;
        sum_xx += x * x;
        sum_yy += y * y;
    }

    double N = (double)n_aliquots;
    double denom = N * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1.0e-30) { free(added_conc); return false; }

    double slope     = (N * sum_xy - sum_x * sum_y) / denom;
    double intercept = (sum_y - slope * sum_x) / N;

    /* Standard addition regression:
     *   signal = slope * (C_sample_diluted + C_added)
     *          = slope * C_added + slope * C_sample_diluted
     *   intercept = slope * C_sample_diluted
     *   C_sample_diluted = intercept / slope
     *
     * Then correct for dilution: C_orig = C_sample_diluted * V_total / V_sample */
    if (fabs(slope) < 1.0e-15) { free(added_conc); return false; }

    double c_in_aliquot = intercept / slope;
    *original_conc = c_in_aliquot * total_volumes_ml[0] / sample_volume_ml[0];

    if (*original_conc < 0.0) {
        /* Physical constraint: concentration cannot be negative */
        *original_conc = 0.0;
    }

    /* R² */
    double ss_res = 0.0, ss_tot = 0.0;
    double y_mean = sum_y / N;
    for (int i = 0; i < n_aliquots; i++) {
        double pred = slope * added_conc[i] + intercept;
        ss_res += (measured_signals[i] - pred) *
                  (measured_signals[i] - pred);
        ss_tot += (measured_signals[i] - y_mean) *
                  (measured_signals[i] - y_mean);
    }
    *r_squared = (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;

    free(added_conc);
    return (*r_squared > 0.95);
}

/* ============================================================================
 * L5: Outlier detection
 * ============================================================================ */

/**
 * Grubbs' test for a single outlier.
 *
 * G = |x_suspect - x̄| / s
 *
 * Compare against critical value:
 *   G_crit = (N-1)/√N · √(t²/(N-2+t²))
 *   where t = t(α/(2N), N-2) is the two-sided critical t-value.
 *
 * α = 0.05 (95% confidence) is standard.
 */
bool grubbs_outlier_test(const double *values,
                         int n_values,
                         double alpha,
                         int *outlier_index)
{
    if (!values || !outlier_index || n_values < 3) {
        if (outlier_index) *outlier_index = -1;
        return false;
    }

    *outlier_index = -1;

    /* Compute mean */
    double sum = 0.0;
    for (int i = 0; i < n_values; i++) sum += values[i];
    double mean = sum / n_values;

    /* Compute standard deviation */
    double sum_sq_diff = 0.0;
    for (int i = 0; i < n_values; i++) {
        double diff = values[i] - mean;
        sum_sq_diff += diff * diff;
    }
    double std = sqrt(sum_sq_diff / (n_values - 1));
    if (std < 1.0e-15) return false;  /* all values identical */

    /* Find the most extreme value */
    double max_g = 0.0;
    int max_idx = -1;
    for (int i = 0; i < n_values; i++) {
        double g = fabs(values[i] - mean) / std;
        if (g > max_g) { max_g = g; max_idx = i; }
    }

    /* Critical G values for α = 0.05 (two-sided):
     * N=3:1.155, N=4:1.481, N=5:1.715, N=6:1.887, N=7:2.020,
     * N=8:2.126, N=9:2.215, N=10:2.290, N=12:2.412, N=15:2.549,
     * N=20:2.708, N=30:2.908 */
    double g_crit;
    if (n_values <= 30) {
        /* Interpolate from tabulated values */
        double g_tab[] = {1.155, 1.481, 1.715, 1.887, 2.020, 2.126, 2.215,
                          2.290, 2.355, 2.412, 2.462, 2.507, 2.549,
                          2.585, 2.620, 2.708, 2.768, 2.841, 2.908};
        if (n_values <= 12) {
            g_crit = g_tab[n_values - 3];
        } else if (n_values <= 15) {
            g_crit = g_tab[10 + (n_values - 13) / 3];  /* approximate */
        } else {
            g_crit = g_tab[15 + (n_values - 20) / 10]; /* approximate */
        }
    } else {
        g_crit = 3.0;  /* asymptotic value for large N */
    }

    /* Adjust for alpha */
    if (alpha < 0.05) g_crit *= 1.05;
    if (alpha > 0.05) g_crit *= 0.95;

    if (max_g > g_crit) {
        *outlier_index = max_idx;
        return true;
    }
    return false;
}

/**
 * Dixon's Q-test for small datasets (N = 3-7).
 *
 * Testing high outlier: Q = (x_N - x_{N-1}) / (x_N - x_1)
 * Testing low outlier:  Q = (x_2 - x_1) / (x_N - x_1)
 *
 * Critical Q values at α = 0.05:
 *   N=3:0.941, N=4:0.765, N=5:0.642, N=6:0.560, N=7:0.507
 */
bool dixon_q_test(const double *values,
                  int n_values,
                  double alpha,
                  bool is_high_outlier)
{
    if (!values || n_values < 3 || n_values > 7) return false;

    double x1 = values[0];
    double xN = values[n_values - 1];
    double range = xN - x1;
    if (fabs(range) < 1.0e-15) return false;

    double q_stat;
    if (is_high_outlier) {
        q_stat = (xN - values[n_values - 2]) / range;
    } else {
        q_stat = (values[1] - x1) / range;
    }

    /* Critical values at α = 0.05 */
    double q_crit_tab[] = {0.941, 0.765, 0.642, 0.560, 0.507};
    double q_crit = q_crit_tab[n_values - 3];

    if (alpha > 0.05) q_crit *= 0.95;
    if (alpha < 0.05) q_crit *= 1.05;

    return (q_stat > q_crit);
}

/**
 * Coefficient of variation: CV% = (s / x̄) × 100%
 */
double coefficient_of_variation(const double *values, int n_values)
{
    if (!values || n_values < 2) return 0.0;

    double sum = 0.0, sum_sq = 0.0;
    for (int i = 0; i < n_values; i++) {
        sum    += values[i];
        sum_sq += values[i] * values[i];
    }

    double mean = sum / n_values;
    if (fabs(mean) < 1.0e-15) return INFINITY;

    double variance = (sum_sq - sum * sum / n_values) / (n_values - 1.0);
    if (variance < 0.0) variance = 0.0;

    return sqrt(variance) / mean * 100.0;
}

/* ============================================================================
 * L6: Sensor drift and aging
 * ============================================================================ */

/**
 * Two-point recalibration.
 *
 * slope_new = (sig_high - sig_low) / (conc_high - conc_low)
 * intercept_new = sig_low - slope_new · conc_low
 */
void two_point_calibration_adjust(double conc_low, double sig_low,
                                  double conc_high, double sig_high,
                                  const LinearCalibration *old_cal,
                                  LinearCalibration *new_cal)
{
    if (!new_cal) return;

    if (fabs(conc_high - conc_low) < 1.0e-15) {
        if (old_cal) {
            *new_cal = *old_cal;
        } else {
            linear_calibration_init(new_cal);
        }
        return;
    }

    new_cal->slope     = (sig_high - sig_low) / (conc_high - conc_low);
    new_cal->intercept = sig_low - new_cal->slope * conc_low;
    new_cal->r_squared = 1.0;  /* perfectly fits two points */
    new_cal->standard_error_slope     = 0.0;
    new_cal->standard_error_intercept = 0.0;
    new_cal->prediction_interval_95   = 0.0;
}

/**
 * Estimate remaining sensor lifetime from sensitivity decay.
 *
 * S(t) = S₀ · exp(-t/τ)
 *
 * From two sensitivity measurements at t=0 and t = days_in_use:
 *   τ = -days_in_use / ln(S_current / S_initial)
 *
 * Days remaining until S < S_min:
 *   t_remaining = -τ · ln(S_min / S_current)
 */
double estimate_sensor_remaining_life(double initial_sensitivity,
                                      double current_sensitivity,
                                      int days_in_use,
                                      double min_sensitivity)
{
    if (initial_sensitivity <= 0.0 || current_sensitivity <= 0.0 ||
        days_in_use <= 0 || min_sensitivity <= 0.0) return 0.0;

    if (current_sensitivity >= initial_sensitivity) return 365.0;  /* no decay */
    if (current_sensitivity <= min_sensitivity) return 0.0;         /* already dead */

    double ratio = current_sensitivity / initial_sensitivity;
    if (ratio <= 0.0 || ratio >= 1.0) return 0.0;

    double tau = -(double)days_in_use / log(ratio);
    if (tau <= 0.0) return 0.0;

    double remaining = -tau * log(min_sensitivity / current_sensitivity);
    return (remaining > 0.0) ? remaining : 0.0;
}

/* ============================================================================
 * L6: Calibration quality metrics
 * ============================================================================ */

/**
 * Mandel fitting test: compares linear model error to pure replicate error.
 *
 * SS_LOF = Σ n_i · (ȳ_i - ŷ_i)²  (lack-of-fit sum of squares)
 * SS_PE  = Σ Σ (y_ij - ȳ_i)²    (pure error sum of squares)
 *
 * F = (SS_LOF / (k - 2)) / (SS_PE / (N - k))
 *
 * If F > F_crit, the linear model is inadequate.
 */
bool mandel_fitting_test(const CalibrationStandards *standards,
                         const LinearCalibration *lc,
                         double f_critical)
{
    if (!standards || !lc || standards->replicates <= 1 ||
        standards->standard_count < 3) return true;  /* cannot test, assume OK */

    /* For simplicity, we use the signal_stds to estimate pure error.
     * If signal_stds are available, compute:
     * SS_PE = Σ (n_i - 1) · s_i²
     * SS_LOF = Σ n_i · (signal_i - ŷ_i)²
     */
    double ss_pe = 0.0;
    for (int i = 0; i < standards->standard_count; i++) {
        if (standards->signal_stds[i] > 0.0) {
            ss_pe += (standards->replicates - 1.0) *
                     standards->signal_stds[i] * standards->signal_stds[i];
        }
    }

    double ss_lof = 0.0;
    for (int i = 0; i < standards->standard_count; i++) {
        double pred = lc->slope * standards->concentrations[i] + lc->intercept;
        double res  = standards->signals[i] - pred;
        ss_lof += standards->replicates * res * res;
    }

    /* This is a simplified Mandel test using aggregate statistics.
     * Full implementation requires individual replicate data. */
    if (ss_pe <= 0.0) return true;  /* cannot compute pure error */

    double df_lof = standards->standard_count - 2.0;  /* k - 2 */
    double df_pe  = standards->standard_count * (standards->replicates - 1.0);

    double f_stat = (ss_lof / df_lof) / (ss_pe / df_pe);
    return (f_stat <= f_critical);
}

/**
 * Durbin-Watson test for residual autocorrelation.
 *
 * DW = Σ(e_i - e_{i-1})² / Σ e_i²
 *
 * Interpretation:
 *   DW ≈ 2: no autocorrelation (ideal)
 *   DW < 2: positive autocorrelation
 *   DW > 2: negative autocorrelation
 *   DW < d_L: significant positive autocorrelation at 5% level
 */
bool durbin_watson_test(const double *residuals,
                        int n_points,
                        double *dw_statistic)
{
    if (!residuals || !dw_statistic || n_points < 3) {
        if (dw_statistic) *dw_statistic = 2.0;
        return true;
    }

    double sum_diff_sq = 0.0;
    double sum_res_sq  = 0.0;

    for (int i = 0; i < n_points; i++) {
        sum_res_sq += residuals[i] * residuals[i];
        if (i > 0) {
            double diff = residuals[i] - residuals[i - 1];
            sum_diff_sq += diff * diff;
        }
    }

    *dw_statistic = (sum_res_sq > 0.0) ? (sum_diff_sq / sum_res_sq) : 2.0;

    /* For typical calibration sizes (k = 5-10), d_L ≈ 0.6 to 1.0 */
    double d_lower;
    if (n_points == 5)  d_lower = 0.61;
    else if (n_points == 6)  d_lower = 0.70;
    else if (n_points == 7)  d_lower = 0.78;
    else if (n_points == 8)  d_lower = 0.84;
    else if (n_points == 9)  d_lower = 0.90;
    else if (n_points == 10) d_lower = 0.95;
    else d_lower = 1.0;

    /* DW > d_upper (≈ 1.7 for small k): no evidence of autocorrelation */
    return (*dw_statistic > d_lower);
}

/* ============================================================================
 * L8: Calibration transfer — PDS
 * ============================================================================ */

/**
 * Piecewise Direct Standardization (PDS).
 *
 * For each channel j on the slave, a window of channels
 * [j-k, j+k] on the master is used to predict the slave
 * response via local PLS (here simplified to local OLS):
 *
 * s_j = M_{j, j-k:j+k} · b_j
 *
 * b_j = (M_window^T · M_window)⁻¹ · M_window^T · s_j
 *
 * The transfer matrix F transforms master spectra to slave-like spectra:
 *   X_slave ≈ X_master · F
 *
 * L8 Advanced Topic: Chemometric calibration transfer enables
 * deploying a calibration model built on a laboratory instrument
 * to field-deployable sensors without full recalibration.
 */
double pds_calibration_transfer(const double *master_responses,
                                const double *slave_responses,
                                int n_standards,
                                int n_channels,
                                int window_half_width,
                                double *transfer_matrix)
{
    if (!master_responses || !slave_responses || !transfer_matrix ||
        n_standards < 3 || n_channels < 2 || window_half_width < 0) return -1.0;

    /* Initialize transfer matrix as identity */
    for (int i = 0; i < n_channels; i++) {
        for (int j = 0; j < n_channels; j++) {
            transfer_matrix[i * n_channels + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    /* For each slave channel j, perform local regression */
    double total_error = 0.0;
    int valid_channels = 0;

    for (int j = 0; j < n_channels; j++) {
        int win_start = j - window_half_width;
        int win_end   = j + window_half_width;
        if (win_start < 0) win_start = 0;
        if (win_end >= n_channels) win_end = n_channels - 1;
        int win_size = win_end - win_start + 1;

        /* Build local regression: Y = slave responses at channel j,
         * X = master responses in window around j */
        double sum_x[5] = {0};  /* first 5 window channels */
        double sum_xx[5][5] = {{0}};
        double sum_xy[5] = {0};
        double sum_y = 0.0, sum_yy = 0.0;

        for (int s = 0; s < n_standards; s++) {
            double y = slave_responses[s * n_channels + j];
            sum_y  += y;
            sum_yy += y * y;

            for (int w = 0; w < win_size; w++) {
                double x = master_responses[s * n_channels + win_start + w];
                sum_x[w] += x;
                sum_xy[w] += x * y;
                for (int v = 0; v < win_size; v++) {
                    sum_xx[w][v] += master_responses[s * n_channels + win_start + w] *
                                    master_responses[s * n_channels + win_start + v];
                }
            }
        }

        /* Solve local regression: b = (X^T X)⁻¹ X^T y */
        /* For simplicity, store the regression coefficients in the
         * transfer_matrix column j, rows [win_start, win_end].
         * As an approximation, use simple ratio of means. */
        double y_mean = sum_y / n_standards;

        for (int i = 0; i < n_channels; i++) {
            transfer_matrix[i * n_channels + j] = 0.0;  /* clear column */
        }

        /* Simplified: use univariate correlation for each window element */
        for (int w = 0; w < win_size; w++) {
            double x_mean = sum_x[w] / n_standards;
            double cov = (sum_xy[w] / n_standards) - x_mean * y_mean;
            double var_x = (sum_xx[w][w] / n_standards) - x_mean * x_mean;

            if (var_x > 1.0e-15 && n_standards > 2) {
                double beta = cov / var_x;
                transfer_matrix[(win_start + w) * n_channels + j] = beta / win_size;
            }
        }

        /* Compute prediction error for this channel */
        double ch_error = 0.0;
        for (int s = 0; s < n_standards; s++) {
            double pred = 0.0;
            for (int w = 0; w < win_size; w++) {
                pred += transfer_matrix[(win_start + w) * n_channels + j] *
                        master_responses[s * n_channels + win_start + w];
            }
            double err = slave_responses[s * n_channels + j] - pred;
            ch_error += err * err;
        }
        total_error += sqrt(ch_error / n_standards);
        valid_channels++;
    }

    return (valid_channels > 0) ? (total_error / valid_channels) : -1.0;
}
