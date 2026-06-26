#include "measurement_calibration.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>

/* ─── Internal: Gaussian Elimination with Partial Pivoting ───────────────
 * Solves A x = b where A is n x n, b is n x 1.
 * A and b are modified in-place. Solution stored in x.
 * Returns 0 on success, -1 if singular.
 */

static int gauss_elimination(double *A, double *b, double *x, size_t n) {
    for (size_t col = 0; col < n; col++) {
        /* Partial pivoting: find max in column */
        size_t max_row = col;
        double max_val = fabs(A[col * n + col]);
        for (size_t row = col + 1; row < n; row++) {
            double val = fabs(A[row * n + col]);
            if (val > max_val) { max_val = val; max_row = row; }
        }
        if (max_val < 1e-15) return -1;  /* singular */

        /* Swap rows if needed */
        if (max_row != col) {
            for (size_t j = 0; j < n; j++) {
                double tmp = A[col * n + j];
                A[col * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }
            double tmp = b[col]; b[col] = b[max_row]; b[max_row] = tmp;
        }

        /* Eliminate below */
        double pivot = A[col * n + col];
        for (size_t row = col + 1; row < n; row++) {
            double factor = A[row * n + col] / pivot;
            A[row * n + col] = 0.0;
            for (size_t j = col + 1; j < n; j++) {
                A[row * n + j] -= factor * A[col * n + j];
            }
            b[row] -= factor * b[col];
        }
    }

    /* Back substitution */
    for (int i = (int)n - 1; i >= 0; i--) {
        double sum = b[i];
        for (size_t j = i + 1; j < n; j++) {
            sum -= A[i * n + j] * x[j];
        }
        x[i] = sum / A[i * n + i];
    }
    return 0;
}

/* ─── L5: Linear Least Squares ───────────────────────────────────────────
 * y = a0 + a1 * x
 * Normal equations:
 *   [ n      SUM x  ] [a0] = [SUM y  ]
 *   [ SUM x  SUM x^2] [a1]   [SUM x*y]
 */

int cal_linear_fit(const CalibrationPoint *points, size_t n,
                    CalibrationResult *result) {
    if (!points || !result || n < 2) return -1;
    memset(result, 0, sizeof(*result));
    result->model_type = CAL_LINEAR;
    result->order = 1;
    result->num_points = n;

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum_x  += points[i].x;
        sum_y  += points[i].y;
        sum_xy += points[i].x * points[i].y;
        sum_x2 += points[i].x * points[i].x;
    }

    double det = n * sum_x2 - sum_x * sum_x;
    if (fabs(det) < 1e-15) return -1;

    double a0 = (sum_x2 * sum_y - sum_x * sum_xy) / det;
    double a1 = (n * sum_xy - sum_x * sum_y) / det;

    result->coefficients[0] = a0;
    result->coefficients[1] = a1;
    result->num_coeffs = 2;

    /* Compute fit statistics */
    double ss_res = 0.0, ss_tot = 0.0;
    double y_mean = sum_y / (double)n;
    double max_res = 0.0, sum_res = 0.0;
    for (size_t i = 0; i < n; i++) {
        double y_pred = a0 + a1 * points[i].x;
        double res = points[i].y - y_pred;
        ss_res += res * res;
        ss_tot += (points[i].y - y_mean) * (points[i].y - y_mean);
        sum_res += res;
        if (fabs(res) > max_res) max_res = fabs(res);
    }
    result->r_squared = (ss_tot > 1e-15) ? (1.0 - ss_res / ss_tot) : 1.0;
    result->rmse = sqrt(ss_res / (double)n);
    result->max_residual = max_res;
    result->mean_residual = sum_res / (double)n;

    return 0;
}

/* ─── L5: Polynomial Least Squares ────────────────────────────────────────
 * y = a0 + a1*x + a2*x^2 + ... + ak*x^k
 * Uses Vandermonde normal equations: (V^T V) a = V^T y
 * Solved via Gaussian elimination.
 */

int cal_polynomial_fit(const CalibrationPoint *points, size_t n,
                        size_t order, CalibrationResult *result) {
    if (!points || !result || n < order + 1 || order < 1 || order > 4)
        return -1;
    memset(result, 0, sizeof(*result));
    result->model_type = CAL_POLYNOMIAL;
    result->order = order;
    result->num_points = n;

    size_t m = order + 1;  /* number of coefficients */

    /* Build normal equations: A_ij = SUM x^{i+j}, b_i = SUM x^i * y */
    double *A = (double*)calloc(m * m, sizeof(double));
    double *b = (double*)calloc(m, sizeof(double));
    if (!A || !b) { free(A); free(b); return -1; }

    for (size_t i = 0; i < n; i++) {
        double xi = points[i].x;
        double yi = points[i].y;
        double xi_pow = 1.0;
        for (size_t p = 0; p < m; p++) {
            b[p] += xi_pow * yi;
            /* xj_pow must start at x^{2p} since A[p][q] accumulates Σ x^{p+q} */
            double xj_pow = xi_pow * xi_pow;
            for (size_t q = p; q < m; q++) {
                A[p * m + q] += xj_pow;
                xj_pow *= xi;
            }
            xi_pow *= xi;
        }
    }
    /* Symmetrize A */
    for (size_t i = 0; i < m; i++)
        for (size_t j = i + 1; j < m; j++)
            A[j * m + i] = A[i * m + j];

    double *x = (double*)calloc(m, sizeof(double));
    int ret = gauss_elimination(A, b, x, m);
    if (ret == 0) {
        for (size_t i = 0; i < m; i++) {
            result->coefficients[i] = x[i];
        }
        result->num_coeffs = m;

        /* Compute statistics */
        double ss_res = 0.0, ss_tot = 0.0;
        double sum_y = 0.0;
        double max_res = 0.0, sum_res = 0.0;
        for (size_t i = 0; i < n; i++) sum_y += points[i].y;
        double y_mean = sum_y / (double)n;

        for (size_t i = 0; i < n; i++) {
            double y_pred = 0.0;
            double xi_pow = 1.0;
            for (size_t j = 0; j < m; j++) {
                y_pred += result->coefficients[j] * xi_pow;
                xi_pow *= points[i].x;
            }
            double res = points[i].y - y_pred;
            ss_res += res * res;
            ss_tot += (points[i].y - y_mean) * (points[i].y - y_mean);
            sum_res += res;
            if (fabs(res) > max_res) max_res = fabs(res);
        }
        result->r_squared = (ss_tot > 1e-15) ? (1.0 - ss_res / ss_tot) : 1.0;
        result->rmse = sqrt(ss_res / (double)n);
        result->max_residual = max_res;
        result->mean_residual = sum_res / (double)n;
    }

    free(A); free(b); free(x);
    return ret;
}

/* ─── L5: Weighted Least Squares ─────────────────────────────────────────
 * Minimizes SUM w_i * (y_i - f(x_i))^2
 * W should be diagonal: W_ii = w_i
 * Normal equations: (V^T W V) a = V^T W y
 */

int cal_weighted_poly_fit(const CalibrationPoint *points,
                          const double *weights, size_t n,
                          size_t order, CalibrationResult *result) {
    if (!points || !result || n < order + 1 || order < 1 || order > 4)
        return -1;

    /* If no weights provided, fall back to OLS */
    if (!weights) return cal_polynomial_fit(points, n, order, result);

    memset(result, 0, sizeof(*result));
    result->model_type = CAL_POLYNOMIAL;
    result->order = order;
    result->num_points = n;

    size_t m = order + 1;
    double *A = (double*)calloc(m * m, sizeof(double));
    double *b = (double*)calloc(m, sizeof(double));
    if (!A || !b) { free(A); free(b); return -1; }

    for (size_t i = 0; i < n; i++) {
        double xi = points[i].x;
        double yi = points[i].y;
        double wi = (weights[i] > 0.0) ? weights[i] : 1.0;
        double xi_pow = 1.0;
        for (size_t p = 0; p < m; p++) {
            b[p] += wi * xi_pow * yi;
            double xj_pow = xi_pow * xi_pow;
            for (size_t q = p; q < m; q++) {
                A[p * m + q] += wi * xj_pow;
                xj_pow *= xi;
            }
            xi_pow *= xi;
        }
    }
    for (size_t i = 0; i < m; i++)
        for (size_t j = i + 1; j < m; j++)
            A[j * m + i] = A[i * m + j];

    double *x = (double*)calloc(m, sizeof(double));
    int ret = gauss_elimination(A, b, x, m);
    if (ret == 0) {
        for (size_t i = 0; i < m; i++)
            result->coefficients[i] = x[i];
        result->num_coeffs = m;

        double ss_res = 0.0, sum_w = 0.0, sum_wy = 0.0;
        double max_res = 0.0;
        for (size_t i = 0; i < n; i++) {
            double wi = (weights[i] > 0.0) ? weights[i] : 1.0;
            sum_w += wi;
            sum_wy += wi * points[i].y;

            double y_pred = 0.0;
            double xi_pow = 1.0;
            for (size_t j = 0; j < m; j++) {
                y_pred += result->coefficients[j] * xi_pow;
                xi_pow *= points[i].x;
            }
            double res = points[i].y - y_pred;
            ss_res += wi * res * res;
            if (fabs(res) > max_res) max_res = fabs(res);
        }

        double y_wmean = (sum_w > 0.0) ? sum_wy / sum_w : 0.0;
        double ss_tot = 0.0;
        for (size_t i = 0; i < n; i++) {
            double wi = (weights[i] > 0.0) ? weights[i] : 1.0;
            ss_tot += wi * (points[i].y - y_wmean) * (points[i].y - y_wmean);
        }
        result->r_squared = (ss_tot > 1e-15) ? (1.0 - ss_res / ss_tot) : 1.0;
        result->rmse = sqrt(ss_res / (double)n);
        result->max_residual = max_res;
    }

    free(A); free(b); free(x);
    return ret;
}

/* ─── L5: Exponential Fit: y = a * exp(b * x) ────────────────────────────
 * Linearized: ln(y) = ln(a) + b * x
 * This minimizes relative errors in y (log-space residuals).
 */
int cal_exponential_fit(const CalibrationPoint *points, size_t n,
                         CalibrationResult *result) {
    if (!points || !result || n < 2) return -1;

    /* Check all y > 0 for log */
    for (size_t i = 0; i < n; i++) {
        if (points[i].y <= 0.0) return -1;
    }

    /* Transform to log space */
    CalibrationPoint *log_points = (CalibrationPoint*)malloc(
        n * sizeof(CalibrationPoint));
    if (!log_points) return -1;
    for (size_t i = 0; i < n; i++) {
        log_points[i].x = points[i].x;
        log_points[i].y = log(points[i].y);
    }

    int ret = cal_linear_fit(log_points, n, result);
    if (ret == 0) {
        result->model_type = CAL_EXPONENTIAL;
        /* a = exp(a0), b = a1 */
        result->coefficients[0] = exp(result->coefficients[0]);
        /* result->coefficients[1] stays as b */
        /* Recompute fit stats in original space */
        double ss_res = 0.0, sum_y = 0.0, max_res = 0.0;
        for (size_t i = 0; i < n; i++) sum_y += points[i].y;
        double y_mean = sum_y / (double)n;
        double ss_tot = 0.0;
        for (size_t i = 0; i < n; i++) {
            double y_pred = result->coefficients[0]
                            * exp(result->coefficients[1] * points[i].x);
            double res = points[i].y - y_pred;
            ss_res += res * res;
            ss_tot += (points[i].y - y_mean) * (points[i].y - y_mean);
            if (fabs(res) > max_res) max_res = fabs(res);
        }
        result->r_squared = (ss_tot > 1e-15) ? (1.0 - ss_res / ss_tot) : 1.0;
        result->rmse = sqrt(ss_res / (double)n);
        result->max_residual = max_res;
    }
    free(log_points);
    return ret;
}

/* ─── L5: Power-Law Fit: y = a * x^b ────────────────────────────────────
 * Linearized: ln(y) = ln(a) + b * ln(x)
 * Requires x > 0, y > 0.
 */
int cal_power_law_fit(const CalibrationPoint *points, size_t n,
                       CalibrationResult *result) {
    if (!points || !result || n < 2) return -1;
    for (size_t i = 0; i < n; i++) {
        if (points[i].x <= 0.0 || points[i].y <= 0.0) return -1;
    }

    CalibrationPoint *log_points = (CalibrationPoint*)malloc(
        n * sizeof(CalibrationPoint));
    if (!log_points) return -1;
    for (size_t i = 0; i < n; i++) {
        log_points[i].x = log(points[i].x);
        log_points[i].y = log(points[i].y);
    }

    int ret = cal_linear_fit(log_points, n, result);
    if (ret == 0) {
        result->model_type = CAL_POWER_LAW;
        result->coefficients[0] = exp(result->coefficients[0]);
        /* result->coefficients[1] = b */
        double ss_res = 0.0, sum_y = 0.0, max_res = 0.0;
        for (size_t i = 0; i < n; i++) sum_y += points[i].y;
        double y_mean = sum_y / (double)n;
        double ss_tot = 0.0;
        for (size_t i = 0; i < n; i++) {
            double y_pred = result->coefficients[0]
                            * pow(points[i].x, result->coefficients[1]);
            double res = points[i].y - y_pred;
            ss_res += res * res;
            ss_tot += (points[i].y - y_mean) * (points[i].y - y_mean);
            if (fabs(res) > max_res) max_res = fabs(res);
        }
        result->r_squared = (ss_tot > 1e-15) ? (1.0 - ss_res / ss_tot) : 1.0;
        result->rmse = sqrt(ss_res / (double)n);
        result->max_residual = max_res;
    }
    free(log_points);
    return ret;
}

/* ─── L5: Logarithmic Fit: y = a + b * ln(x) ────────────────────────────
 * Requires x > 0.
 */
int cal_logarithmic_fit(const CalibrationPoint *points, size_t n,
                         CalibrationResult *result) {
    if (!points || !result || n < 2) return -1;
    for (size_t i = 0; i < n; i++) {
        if (points[i].x <= 0.0) return -1;
    }

    CalibrationPoint *ln_points = (CalibrationPoint*)malloc(
        n * sizeof(CalibrationPoint));
    if (!ln_points) return -1;
    for (size_t i = 0; i < n; i++) {
        ln_points[i].x = log(points[i].x);
        ln_points[i].y = points[i].y;
    }

    int ret = cal_linear_fit(ln_points, n, result);
    if (ret == 0) {
        result->model_type = CAL_LOGARITHMIC;
    }
    free(ln_points);
    return ret;
}

/* ─── L5: Model Evaluation and Inverse ──────────────────────────────────── */

double cal_evaluate(const CalibrationResult *result, double x) {
    if (!result || result->num_coeffs == 0) return 0.0;
    const double *c = result->coefficients;

    switch (result->model_type) {
        case CAL_LINEAR:
            return c[0] + c[1] * x;
        case CAL_POLYNOMIAL: {
            double y = 0.0, xp = 1.0;
            for (size_t i = 0; i < result->num_coeffs; i++) {
                y += c[i] * xp;
                xp *= x;
            }
            return y;
        }
        case CAL_EXPONENTIAL:
            return c[0] * exp(c[1] * x);
        case CAL_LOGARITHMIC:
            if (x <= 0.0) return -INFINITY;
            return c[0] + c[1] * log(x);
        case CAL_POWER_LAW:
            if (x < 0.0) return 0.0;
            return c[0] * pow(x, c[1]);
        default:
            return 0.0;
    }
}

/* Inverse evaluation: given y, find x such that f(x) = y.
 * For linear: x = (y - a0) / a1
 * For others: uses Newton-Raphson.
 */
int cal_inverse_evaluate(const CalibrationResult *result, double y,
                          double x0, double *x_est) {
    if (!result || !x_est || result->num_coeffs < 2) return -1;
    const double *c = result->coefficients;

    if (result->model_type == CAL_LINEAR) {
        if (fabs(c[1]) < DBL_EPSILON) return -1;
        *x_est = (y - c[0]) / c[1];
        return 0;
    }

    /* Newton-Raphson: x_{k+1} = x_k - (f(x_k) - y) / f'(x_k) */
    double xk = x0;
    for (int iter = 0; iter < 50; iter++) {
        double f_val, f_prime;
        switch (result->model_type) {
            case CAL_EXPONENTIAL:
                f_val = c[0] * exp(c[1] * xk) - y;
                f_prime = c[0] * c[1] * exp(c[1] * xk);
                break;
            case CAL_POWER_LAW:
                if (xk <= 0.0) xk = 1e-6;
                f_val = c[0] * pow(xk, c[1]) - y;
                f_prime = c[0] * c[1] * pow(xk, c[1] - 1.0);
                break;
            case CAL_LOGARITHMIC:
                if (xk <= 0.0) xk = 1e-6;
                f_val = c[0] + c[1] * log(xk) - y;
                f_prime = c[1] / xk;
                break;
            case CAL_POLYNOMIAL: {
                f_val = -y; f_prime = 0.0;
                double xp = 1.0, xp_prime = 0.0;
                for (size_t i = 0; i < result->num_coeffs; i++) {
                    f_val += c[i] * xp;
                    f_prime += c[i] * (double)i * xp_prime;
                    xp_prime = xp;
                    xp *= xk;
                }
                break;
            }
            default:
                return -1;
        }
        double dx = (fabs(f_prime) > 1e-15) ? f_val / f_prime : 0.0;
        xk -= dx;
        if (fabs(dx) < 1e-10) break;
    }
    *x_est = xk;
    return 0;
}

/* ─── L5: Fit Statistics ────────────────────────────────────────────────── */

double cal_compute_r_squared(const CalibrationPoint *points, size_t n,
                              CalibrationResult *result) {
    if (!points || !result || n == 0) return 0.0;
    double ss_res = 0.0, sum_y = 0.0;
    for (size_t i = 0; i < n; i++) sum_y += points[i].y;
    double y_mean = sum_y / (double)n;
    double ss_tot = 0.0;
    for (size_t i = 0; i < n; i++) {
        double y_pred = cal_evaluate(result, points[i].x);
        double res = points[i].y - y_pred;
        ss_res += res * res;
        ss_tot += (points[i].y - y_mean) * (points[i].y - y_mean);
    }
    return (ss_tot > 1e-15) ? (1.0 - ss_res / ss_tot) : 1.0;
}

double cal_compute_rmse(const CalibrationPoint *points, size_t n,
                         CalibrationResult *result) {
    if (!points || !result || n == 0) return 0.0;
    double sse = 0.0;
    for (size_t i = 0; i < n; i++) {
        double y_pred = cal_evaluate(result, points[i].x);
        double res = points[i].y - y_pred;
        sse += res * res;
    }
    return sqrt(sse / (double)n);
}

double cal_lack_of_fit_f_test(const CalibrationPoint *points, size_t n,
                               CalibrationResult *result) {
    if (!points || !result || n < 3 || result->num_coeffs >= n)
        return 0.0;
    double sse_full = 0.0;
    for (size_t i = 0; i < n; i++) {
        double y_pred = cal_evaluate(result, points[i].x);
        double res = points[i].y - y_pred;
        sse_full += res * res;
    }
    /* Pure error: SSE of full model (fit each distinct x) approximated
     * by total variance times appropriate df. For simple case, use
     * lack-of-fit F = (SSE_model / (n-p)) / (SSE_pure / (n_unique-1)).
     * Here we return the F-statistic; comparison with F_critical is
     * left to the caller. */
    size_t p = result->num_coeffs;
    double ms_lof = sse_full / (double)(n - p);
    /* Approximate pure error as variance of y at repeated x */
    double ss_pe = 0.0;
    size_t n_pe = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (fabs(points[i].x - points[j].x) < 1e-10) {
                double diff = points[i].y - points[j].y;
                ss_pe += diff * diff * 0.5;
                n_pe++;
            }
        }
    }
    if (n_pe < 1 || ss_pe < 1e-15) return 0.0;
    double ms_pe = ss_pe / (double)n_pe;
    return ms_lof / ms_pe;
}

/* ─── L5: Prediction Uncertainty (ISO 11095) ─────────────────────────────
 *
 * For linear model y = a0 + a1*x:
 * u_pred^2 = s^2 * (1 + 1/n + (x - xbar)^2 / Sxx)
 * where s^2 = SSE / (n - 2), xbar = mean(x), Sxx = SUM (xi - xbar)^2
 */

double cal_prediction_uncertainty(const CalibrationResult *result,
                                   double x, const CalibrationPoint *points,
                                   size_t n) {
    if (!result || !points || n < 3) return INFINITY;
    if (result->model_type != CAL_LINEAR) return INFINITY;

    /* Compute xbar and Sxx */
    double sum_x = 0.0;
    for (size_t i = 0; i < n; i++) sum_x += points[i].x;
    double xbar = sum_x / (double)n;

    double sxx = 0.0;
    for (size_t i = 0; i < n; i++) {
        double dx = points[i].x - xbar;
        sxx += dx * dx;
    }

    /* Compute s^2 = SSE / (n - 2) */
    double sse = 0.0;
    for (size_t i = 0; i < n; i++) {
        double y_pred = result->coefficients[0]
                        + result->coefficients[1] * points[i].x;
        double res = points[i].y - y_pred;
        sse += res * res;
    }
    double s2 = (n > 2) ? sse / (double)(n - 2) : 0.0;

    double dx = x - xbar;
    return sqrt(s2 * (1.0 + 1.0 / (double)n + (dx * dx) / sxx));
}
