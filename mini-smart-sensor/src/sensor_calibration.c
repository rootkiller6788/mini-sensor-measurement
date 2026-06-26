/**
 * @file    sensor_calibration.c
 * @brief   Sensor calibration algorithms: regression, thermistor, thermocouple
 *
 * Knowledge Coverage:
 *   L5 (Algorithms): Linear regression, polynomial regression (normal equations),
 *                    Steinhart-Hart 3-point calibration, Beta parameter thermistor,
 *                    ITS-90 thermocouple inversion, two-point calibration,
 *                    multi-point LUT interpolation, polynomial evaluation
 *   L4 (Fundamental Laws): Gauss-Markov theorem (OLS is BLUE),
 *                          Propagation of uncertainty (GUM),
 *                          Standard error of the estimate
 *   L6 (Canonical Problems): Thermistor calibration, thermocouple calibration,
 *                            pressure sensor linearization
 *
 * Reference:
 *   Press, W.H. et al. (2007), "Numerical Recipes", 3rd ed., Cambridge Univ. Press
 *   Steinhart & Hart (1968), Deep Sea Research 15(4):497-503
 *   NIST Monograph 175 (1993) — ITS-90 Thermocouple Database
 *   JCGM 100:2008 (GUM) — Evaluation of measurement data
 *   Gauss, C.F. (1809), Theoria Motus Corporum Coelestium (least squares)
 *
 * Course Mapping:
 *   MIT 6.003 — Signals and Systems (least squares)
 *   Berkeley EE123 — Digital Signal Processing (numerical methods)
 *   Stanford EE102A — Signal Processing (parameter estimation)
 *   ETH 227-0427 — Signal Processing (statistical signal processing)
 */

#include "smart_sensor.h"
#include "sensor_calibration.h"
#include <string.h>
#include <math.h>
#include <stdio.h>

/* =========================================================================
 * L5: Linear Regression (Ordinary Least Squares)
 *
 * Gauss-Markov Theorem (1821): Under the assumptions of
 *   (1) linearity in parameters
 *   (2) E[epsilon_i] = 0 (zero-mean errors)
 *   (3) Var[epsilon_i] = sigma^2 (homoscedasticity)
 *   (4) Cov[epsilon_i, epsilon_j] = 0 for i != j (uncorrelated errors)
 * the OLS estimator is the Best Linear Unbiased Estimator (BLUE).
 *
 * Closed-form solution:
 *   Let Sx = sum(x_i), Sy = sum(y_i), Sxx = sum(x_i^2), Sxy = sum(x_i*y_i)
 *   slope a = (n*Sxy - Sx*Sy) / (n*Sxx - Sx^2)
 *   intercept b = (Sy - a*Sx) / n
 *
 * R-squared (coefficient of determination):
 *   R^2 = 1 - SS_res / SS_tot
 *   where SS_res = sum((y_i - (a*x_i+b))^2), SS_tot = sum((y_i - mean(y))^2)
 *
 * Complexity: O(n) — single pass for all running sums
 * ========================================================================= */

int ss_cal_linear_regression(const double *x, const double *y,
                             size_t n, ss_linear_model_t *model_out) {
    size_t i;
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    double denom, a, b, mean_y, ss_res, ss_tot;

    if (!x || !y || !model_out || n < 2) return -1;

    /* Compute running sums (single pass for numerical efficiency) */
    for (i = 0; i < n; i++) {
        sx  += x[i];
        sy  += y[i];
        sxx += x[i] * x[i];
        sxy += x[i] * y[i];
    }

    denom = (double)n * sxx - sx * sx;
    if (fabs(denom) < 1e-30) return -1; /* Near-degenerate: all x same */

    /* Slope and intercept */
    a = ((double)n * sxy - sx * sy) / denom;
    b = (sy - a * sx) / (double)n;

    /* Compute R-squared and RMS residual */
    mean_y = sy / (double)n;
    ss_res = 0.0;
    ss_tot = 0.0;
    for (i = 0; i < n; i++) {
        double y_pred = a * x[i] + b;
        double resid = y[i] - y_pred;
        double dev   = y[i] - mean_y;
        ss_res += resid * resid;
        ss_tot += dev * dev;
    }

    model_out->sensitivity = a;
    model_out->offset = b;
    model_out->r_squared = (ss_tot > 0.0) ? 1.0 - ss_res / ss_tot : 1.0;
    model_out->residuals_rms = sqrt(ss_res / (double)n);

    return 0;
}

/**
 * @brief Apply linear model forward: y = sensitivity * x + offset
 *
 * @param model  Calibrated linear model
 * @param x      Sensor reading (raw value)
 * @return       Estimated true value in engineering units
 */
double ss_cal_linear_apply(const ss_linear_model_t *model, double x) {
    if (!model) return x;
    return model->sensitivity * x + model->offset;
}

/**
 * @brief Apply linear model inverse: x = (y - offset) / sensitivity
 *
 * Used when you have a target output and need to compute the
 * corresponding input (e.g., calibration verification).
 *
 * @param model  Calibrated linear model
 * @param y      Target value in engineering units
 * @return       Estimated raw sensor reading
 */
double ss_cal_linear_inverse(const ss_linear_model_t *model, double y) {
    if (!model || fabs(model->sensitivity) < 1e-30) return y;
    return (y - model->offset) / model->sensitivity;
}

/* =========================================================================
 * L5: Polynomial Regression via Normal Equations
 *
 * Model: y = a_0 + a_1*x + a_2*x^2 + ... + a_d*x^d
 *
 * Normal equations: (X^T * X) * a = X^T * y
 *
 * X is the n x (d+1) Vandermonde matrix:
 *   X[i][j] = x_i^j
 *
 * System solved via Gaussian elimination with partial pivoting.
 * The normal matrix (X^T * X) is symmetric positive-definite
 * (assuming n > d and x_i are distinct), ensuring a unique solution.
 *
 * Complexity: O(n * d^2) to form normal matrix + O(d^3) to solve
 *
 * Caution: For high-degree polynomials (d > 5) with limited data,
 * the normal matrix can be ill-conditioned. Consider regularization
 * (ridge regression) or orthogonal polynomial bases for such cases.
 * ========================================================================= */

/* Solve A*x = b for d x d system via Gaussian elimination with pivoting */
static int gauss_solve(double *A, double *b, int d) {
    int i, j, k, max_row;
    double max_val, factor, tmp;

    /* Forward elimination with partial pivoting */
    for (i = 0; i < d; i++) {
        /* Find pivot row */
        max_row = i;
        max_val = fabs(A[i * d + i]);
        for (k = i + 1; k < d; k++) {
            if (fabs(A[k * d + i]) > max_val) {
                max_val = fabs(A[k * d + i]);
                max_row = k;
            }
        }

        /* Check for singularity */
        if (max_val < 1e-30) return -1;

        /* Swap rows if needed */
        if (max_row != i) {
            for (j = 0; j < d; j++) {
                tmp = A[i * d + j];
                A[i * d + j] = A[max_row * d + j];
                A[max_row * d + j] = tmp;
            }
            tmp = b[i];
            b[i] = b[max_row];
            b[max_row] = tmp;
        }

        /* Eliminate below */
        for (k = i + 1; k < d; k++) {
            factor = A[k * d + i] / A[i * d + i];
            for (j = i; j < d; j++) {
                A[k * d + j] -= factor * A[i * d + j];
            }
            b[k] -= factor * b[i];
        }
    }

    /* Back substitution */
    for (i = d - 1; i >= 0; i--) {
        for (j = i + 1; j < d; j++) {
            b[i] -= A[i * d + j] * b[j];
        }
        b[i] /= A[i * d + i];
    }

    return 0;
}

int ss_cal_polynomial_regression(const double *x, const double *y,
                                 size_t n, uint8_t order,
                                 ss_polynomial_model_t *model_out) {
    int d;
    size_t i;
    double *ATA, *ATy;
    double x_min_val, x_max_val;
    int result;

    if (!x || !y || !model_out || n < 2 || order < 1 || order > 7) {
        return -1;
    }
    /* Need at least order+1 points for a unique solution */
    if (n < (size_t)(order + 1)) return -1;

    d = (int)order + 1;

    /* Allocate normal matrix and RHS on stack for small systems
     * (max d=8 for order=7, so 64 doubles for ATA, 8 for ATy — safe) */
    {
        double ATA_local[64], ATy_local[8];
        int j, k;

        ATA = ATA_local;
        ATy = ATy_local;
        memset(ATA, 0, d * d * sizeof(double));
        memset(ATy, 0, d * sizeof(double));

        /* Form ATA (X^T * X) and ATy (X^T * y) */
        for (i = 0; i < n; i++) {
            double x_pow[8];  /* x^0, x^1, ..., x^d */
            x_pow[0] = 1.0;
            for (j = 1; j < d; j++) {
                x_pow[j] = x_pow[j - 1] * x[i];
            }
            for (j = 0; j < d; j++) {
                for (k = 0; k < d; k++) {
                    ATA[j * d + k] += x_pow[j] * x_pow[k];
                }
                ATy[j] += x_pow[j] * y[i];
            }
        }

        /* Solve ATA * a = ATy via Gaussian elimination */
        result = gauss_solve(ATA, ATy, d);
        if (result != 0) return -1;

        /* Populate output model (coefficients: highest degree first) */
        model_out->order = order;
        for (j = 0; j < d; j++) {
            model_out->coeffs[order - j] = ATy[j]; /* [a_n, ..., a_1, a_0] */
        }

        /* Find x range */
        x_min_val = x[0];
        x_max_val = x[0];
        for (i = 1; i < n; i++) {
            if (x[i] < x_min_val) x_min_val = x[i];
            if (x[i] > x_max_val) x_max_val = x[i];
        }
        model_out->x_min = x_min_val;
        model_out->x_max = x_max_val;

        /* Compute R-squared and max residual */
        {
            double ss_res = 0.0, ss_tot = 0.0, mean_y = 0.0;
            model_out->max_residual = 0.0;

            for (i = 0; i < n; i++) mean_y += y[i];
            mean_y /= (double)n;

            for (i = 0; i < n; i++) {
                double y_pred = ss_cal_polynomial_apply(model_out, x[i]);
                double resid = fabs(y[i] - y_pred);
                ss_res += resid * resid;
                ss_tot += (y[i] - mean_y) * (y[i] - mean_y);
                if (resid > model_out->max_residual) {
                    model_out->max_residual = resid;
                }
            }

            model_out->r_squared = (ss_tot > 0.0)
                                   ? 1.0 - ss_res / ss_tot : 1.0;
        }
    }

    return 0;
}

/**
 * @brief Evaluate polynomial via Horner's method
 *
 * Horner's method: p(x) = a_0 + x*(a_1 + x*(a_2 + ... + x*a_n))
 *
 * Advantages over direct evaluation:
 *   - Fewer multiplications: O(d) vs O(d^2)
 *   - Better numerical stability (reduced roundoff)
 *   - Basis of synthetic division
 *
 * Reference: Horner, W.G. (1819), Phil. Trans. R. Soc. Lond. 109:308-335
 *
 * @param model  Polynomial model (coeffs stored highest-degree first)
 * @param x      Point at which to evaluate
 * @return       p(x)
 */
double ss_cal_polynomial_apply(const ss_polynomial_model_t *model, double x) {
    int i;
    double result;
    if (!model) return x;

    /* coeffs[0] = a_n (highest degree), evaluate from highest to lowest */
    result = model->coeffs[0];
    for (i = 1; i <= (int)model->order; i++) {
        result = result * x + model->coeffs[i];
    }
    return result;
}

/**
 * @brief Inverse polynomial via Newton-Raphson
 *
 * Finds x such that p(x) = y_target.
 * Newton iteration: x_{k+1} = x_k - (p(x_k) - y_target) / p'(x_k)
 *
 * The derivative p'(x) is evaluated using the Horner scheme for
 * both p(x) and p'(x) simultaneously (synthetic division).
 *
 * Convergence: quadratic for simple roots (typical for monotonic
 * calibration curves).
 *
 * @param model     Polynomial model
 * @param y_target  Target output value
 * @param x_guess   Initial guess for the root
 * @param tol       Convergence tolerance (absolute error in y)
 * @param max_iter  Maximum iterations (suggest 20-50)
 * @return          Converged x value, NaN if not converged
 */
double ss_cal_polynomial_inverse(const ss_polynomial_model_t *model,
                                 double y_target, double x_guess,
                                 double tol, int max_iter) {
    int iter;
    double x = x_guess;
    double fx, fpx, dx;

    if (!model || max_iter <= 0) return nan("");

    for (iter = 0; iter < max_iter; iter++) {
        /* Evaluate p(x) using Horner */
        {
            int i;
            fx = model->coeffs[0];
            fpx = 0.0;
            for (i = 1; i <= (int)model->order; i++) {
                fpx = fpx * x + fx;  /* Derivative via synthetic division */
                fx = fx * x + model->coeffs[i];
            }
        }

        fx -= y_target;

        if (fabs(fx) < tol) return x;
        if (fabs(fpx) < 1e-30) return nan("");  /* Derivative vanished */

        dx = fx / fpx;
        x -= dx;

        if (fabs(dx) < tol * 1e-3) return x;
    }

    return nan("");  /* Failed to converge */
}

/* =========================================================================
 * L5: Steinhart-Hart Thermistor Calibration
 *
 * The Steinhart-Hart equation provides a highly accurate (±0.01 degC)
 * model for NTC thermistor resistance vs. temperature:
 *
 *   1/T = A + B*ln(R) + C*[ln(R)]^3
 *
 * where T is absolute temperature in Kelvin and R is resistance in Ohms.
 *
 * To calibrate, we need 3 (R, T) data points to solve for A, B, C:
 *   1/T_i = A + B*ln(R_i) + C*[ln(R_i)]^3    for i = 1, 2, 3
 *
 * This is a linear system in A, B, C solved via Cramer's rule.
 * ========================================================================= */

int ss_cal_steinhart_hart_3pt(double r1, double t1,
                              double r2, double t2,
                              double r3, double t3,
                              ss_steinhart_hart_t *model) {
    double L1, L2, L3;
    double Y1, Y2, Y3;
    double det, detA, detB, detC;
    double A, B, C;

    if (!model || r1 <= 0.0 || r2 <= 0.0 || r3 <= 0.0
        || t1 <= 0.0 || t2 <= 0.0 || t3 <= 0.0) return -1;

    /* Compute ln(R) values */
    L1 = log(r1);
    L2 = log(r2);
    L3 = log(r3);

    /* Compute 1/T values */
    Y1 = 1.0 / t1;
    Y2 = 1.0 / t2;
    Y3 = 1.0 / t3;

    /* Linear system:
     *   A + B*L1 + C*L1^3 = Y1
     *   A + B*L2 + C*L2^3 = Y2
     *   A + B*L3 + C*L3^3 = Y3
     *
     * Solve via Cramer's rule for the 3x3 Vandermonde-like system.
     */
    {
        double L1_3 = L1 * L1 * L1;
        double L2_3 = L2 * L2 * L2;
        double L3_3 = L3 * L3 * L3;

        /* Determinant of coefficient matrix */
        det = L1 * (L2_3 - L3_3) - L2 * (L1_3 - L3_3) + L3 * (L1_3 - L2_3);

        if (fabs(det) < 1e-30) return -1;  /* Degenerate points */

        detA = Y1 * (L2 * L3_3 - L3 * L2_3)
             - L1 * (Y2 * L3_3 - Y3 * L2_3)
             + L1_3 * (Y2 * L3 - Y3 * L2);
        detB = 1.0 * (Y2 * L3_3 - Y3 * L2_3)
             - Y1 * (1.0 * L3_3 - 1.0 * L2_3)
             + L1_3 * (1.0 * Y3 - 1.0 * Y2);
        detC = 1.0 * (L2 * Y3 - L3 * Y2)
             - L1 * (1.0 * Y3 - 1.0 * Y2)
             + Y1 * (1.0 * L3 - 1.0 * L2);

        A = detA / det;
        B = detB / det;
        C = detC / det;
    }

    memset(model, 0, sizeof(*model));
    model->coeff_A = A;
    model->coeff_B = B;
    model->coeff_C = C;
    model->r_at_25c = 0.0;  /* Not computed from 3-point calibration */
    model->beta_value = 0.0;

    return 0;
}

/**
 * @brief Simplified Beta-parameter model from 2 calibration points
 *
 * Beta model (simplified Steinhart-Hart):
 *   R(T) = R0 * exp(Beta * (1/T - 1/T0))
 *
 * Inverse:
 *   T = 1 / (1/T0 + ln(R/R0)/Beta)
 *
 * This is less accurate than the full Steinhart-Hart equation but
 * requires only 2 calibration points and is widely used in
 * microcontroller-based thermistor applications.
 *
 * Beta value from two points (R0,T0) and (R1,T1):
 *   Beta = ln(R1/R0) / (1/T1 - 1/T0)
 *
 * @param r0, t0  Reference resistance (Ohm) and temperature (K)
 * @param r1, t1  Second calibration point
 * @param model    Output model with beta_value set
 * @return         0 on success, -1 on invalid input
 */
int ss_cal_thermistor_beta(double r0, double t0,
                           double r1, double t1,
                           ss_steinhart_hart_t *model) {
    double beta;
    if (!model || r0 <= 0.0 || r1 <= 0.0 || t0 <= 0.0 || t1 <= 0.0
        || t0 == t1) {
        return -1;
    }

    beta = log(r1 / r0) / (1.0 / t1 - 1.0 / t0);

    memset(model, 0, sizeof(*model));
    model->r_at_25c = r0;
    model->beta_value = beta;
    model->coeff_A = 1.0 / t0 - log(r0) / beta;

    return 0;
}

/**
 * @brief Convert thermistor resistance to temperature (Steinhart-Hart)
 *
 * Uses the full Steinhart-Hart equation.
 * Falls back to Beta model if C coefficient is zero.
 *
 * @param model       Calibrated Steinhart-Hart model
 * @param resistance  Thermistor resistance (Ohm)
 * @return            Temperature in Kelvin
 */
double ss_cal_steinhart_hart_temp(const ss_steinhart_hart_t *model,
                                  double resistance) {
    double L, inv_T;
    if (!model || resistance <= 0.0) return 0.0;

    L = log(resistance);

    if (fabs(model->coeff_C) > 1e-30) {
        /* Full Steinhart-Hart equation */
        inv_T = model->coeff_A
                + model->coeff_B * L
                + model->coeff_C * L * L * L;
    } else if (fabs(model->beta_value) > 1e-30) {
        /* Beta model fallback */
        inv_T = 1.0 / 298.15 + log(resistance / model->r_at_25c)
                / model->beta_value;
    } else {
        /* Simple linear approximation */
        inv_T = model->coeff_A + model->coeff_B * L;
    }

    if (fabs(inv_T) < 1e-30) return 0.0;
    return 1.0 / inv_T;
}

/**
 * @brief Convert temperature to thermistor resistance
 *
 * Inverse of Steinhart-Hart using Newton-Raphson iteration.
 *
 * @param model   Steinhart-Hart model
 * @param temp_k  Temperature in Kelvin
 * @return        Resistance in Ohms
 */
double ss_cal_steinhart_hart_resistance(const ss_steinhart_hart_t *model,
                                        double temp_k) {
    double L, inv_T;
    int iter;

    if (!model || temp_k <= 0.0) return 0.0;
    inv_T = 1.0 / temp_k;

    /* Initial guess: use log-linear approximation */
    L = (inv_T - model->coeff_A) / model->coeff_B;

    /* Newton-Raphson refinement on ln(R) */
    for (iter = 0; iter < 20; iter++) {
        double L2 = L * L;
        double L3 = L2 * L;
        double f  = model->coeff_A + model->coeff_B * L
                    + model->coeff_C * L3 - inv_T;
        double fp = model->coeff_B + 3.0 * model->coeff_C * L2;

        if (fabs(f) < 1e-12) break;
        if (fabs(fp) < 1e-30) break;

        L -= f / fp;
        if (L < -50.0) L = -50.0;  /* Clamp extreme values */
        if (L > 20.0) L = 20.0;
    }

    return exp(L);
}

/* =========================================================================
 * L5: ITS-90 Thermocouple Calibration
 *
 * NIST ITS-90 provides polynomial coefficients for converting
 * thermocouple EMF (in mV) to temperature (in degC).
 *
 * The inverse polynomial form is:
 *   T = sum_{i=0}^{n} c_i * E^i
 *
 * where E is the thermoelectric EMF in mV.
 *
 * Cold-junction compensation: the measured voltage is relative to
 * the terminal temperature. We must add the EMF that would be
 * produced by a thermocouple at the cold-junction temperature
 * relative to 0 degC.
 *
 *   E_corrected = E_measured + E(cold_junction_temp)
 *
 * Then T = f(E_corrected) where f is the ITS-90 inverse polynomial.
 *
 * Reference: NIST Monograph 175 (1993), "Temperature-Electromotive
 *            Force Reference Functions and Tables for the
 *            Letter-Designated Thermocouple Types Based on the
 *            ITS-90"
 * ========================================================================= */

/**
 * @brief NIST ITS-90 inverse polynomial coefficients for Type K
 *
 * Temperature range: 0 to 1372 degC
 * These coefficients produce T in degC from EMF in mV.
 */
static const double tc_K_coeffs[] = {
     0.000000000000E+00,
     2.517346200000E+01,
    -1.166287800000E+00,
    -1.083363800000E+00,
    -8.977354000000E-01,
    -3.734237700000E-01,
    -8.663264300000E-02,
    -1.045059800000E-02,
    -5.192057700000E-04,
     0.000000000000E+00
};
static const uint8_t tc_K_ncoeffs = 10;

/* Type K forward coefficients (0 to 1372 degC) for CJC */
static const double tc_K_fwd_coeffs[] = {
     -0.176004136860E-01,
      0.389212049750E-01,
      0.185587700320E-04,
     -0.994575928740E-07,
      0.318409457190E-09,
     -0.560728448890E-12,
      0.560750590590E-15,
     -0.320207200030E-18,
      0.971511471520E-22,
     -0.121047212750E-25
};

/* Type J inverse coefficients (0 to 760 degC) */
static const double tc_J_coeffs[] = {
     0.000000000000E+00,
     1.952826800000E+01,
    -1.228618500000E+00,
    -1.075217800000E+00,
    -5.908693300000E-01,
    -1.725671300000E-01,
    -2.813151300000E-02,
    -2.396337000000E-03,
    -8.382332100000E-05,
     0.000000000000E+00
};
static const uint8_t tc_J_ncoeffs = 10;

static const double tc_J_fwd_coeffs[] = {
     -0.048868252E-01,
      0.198734129E-01,
     -0.208090012E-04,
     -0.115715128E-06,
      0.258387250E-09,
     -0.308921740E-12,
      0.225782520E-15,
     -0.100608150E-18,
      0.252032330E-22,
     -0.271282880E-26
};

int ss_cal_thermocouple_init(ss_thermocouple_model_t *model,
                             uint8_t tc_type, double cold_junction) {
    const double *coeffs;
    uint8_t n;

    if (!model || tc_type > 7) return -1;

    memset(model, 0, sizeof(*model));

    /* Select coefficients based on thermocouple type */
    switch (tc_type) {
    case 0: /* Type K */
        coeffs = tc_K_coeffs;
        n = tc_K_ncoeffs;
        model->temp_min = 0.0;
        model->temp_max = 1372.0;
        break;
    case 1: /* Type J */
        coeffs = tc_J_coeffs;
        n = tc_J_ncoeffs;
        model->temp_min = 0.0;
        model->temp_max = 760.0;
        break;
    default:
        /* For other types, use Type K as a reasonable default */
        coeffs = tc_K_coeffs;
        n = tc_K_ncoeffs;
        model->temp_min = 0.0;
        model->temp_max = 1372.0;
        break;
    }

    model->thermocouple_type = tc_type;
    model->cold_junction_temp = cold_junction;
    model->n_coeffs = n;
    memcpy(model->coeffs, coeffs, n * sizeof(double));

    return 0;
}

/**
 * @brief Compute cold-junction compensation EMF using forward polynomial
 *
 * Evaluates the thermocouple EMF at the cold-junction temperature
 * using the forward (temperature-to-EMF) polynomial.
 *
 * @param model      Thermocouple model (must be initialized)
 * @param cj_temp_c  Cold junction temperature in degC
 * @return           EMF correction in mV
 */
double ss_cal_thermocouple_cj_comp(const ss_thermocouple_model_t *model,
                                   double cj_temp_c) {
    const double *coeffs;
    double emf = 0.0;
    double t_pow = 1.0;
    int i;

    if (!model) return 0.0;
    if (cj_temp_c <= 0.0) return 0.0;

    /* Select forward coefficients */
    if (model->thermocouple_type == 0) {
        coeffs = tc_K_fwd_coeffs;
    } else {
        coeffs = tc_J_fwd_coeffs;
    }

    /* Evaluate forward polynomial: EMF = sum(c_i * T^i) */
    for (i = 0; i < 10; i++) {
        emf += coeffs[i] * t_pow;
        t_pow *= cj_temp_c;
    }

    return emf * 1000.0;  /* Convert V to mV */
}

/**
 * @brief Convert measured EMF to temperature with CJC
 *
 * 1. Compute cold-junction compensation voltage
 * 2. Add to measured voltage
 * 3. Evaluate ITS-90 inverse polynomial
 *
 * @param model   Initialized thermocouple model
 * @param emf_mv  Measured EMF at the terminals (mV)
 * @return        Temperature in degC
 */
double ss_cal_thermocouple_to_temp(const ss_thermocouple_model_t *model,
                                   double emf_mv) {
    double emf_total, temp, e_pow;
    int i;

    if (!model) return 0.0;

    /* Cold-junction compensation */
    emf_total = emf_mv + ss_cal_thermocouple_cj_comp(model,
                                                      model->cold_junction_temp);

    if (emf_total <= 0.0) return model->cold_junction_temp;

    /* Evaluate ITS-90 inverse polynomial */
    temp = 0.0;
    e_pow = 1.0;
    for (i = 0; i < (int)model->n_coeffs; i++) {
        temp += model->coeffs[i] * e_pow;
        e_pow *= emf_total;
    }

    return temp;
}

/* =========================================================================
 * L5: Two-Point and Multi-Point Calibration
 * ========================================================================= */

/**
 * @brief Two-point calibration (zero + span)
 *
 * Given two known reference points, computes a linear correction
 * that adjusts both offset (zero) and gain (span).
 *
 * Method:
 *   1. Fit line through (y1_raw, x1_true) and (y2_raw, x2_true)
 *   2. Output model: x_true = sensitivity * y_raw + offset
 *
 * This is the most common field calibration method. For highest
 * accuracy, the two calibration points should bracket the expected
 * measurement range.
 *
 * @param x1_true    True value at point 1
 * @param y1_raw     Raw sensor reading at point 1
 * @param x2_true    True value at point 2
 * @param y2_raw     Raw sensor reading at point 2
 * @param model_out   Output linear model
 * @return           0 on success, -1 if points are degenerate
 */
int ss_cal_two_point(double x1_true, double y1_raw,
                     double x2_true, double y2_raw,
                     ss_linear_model_t *model_out) {
    double dx, dy;
    if (!model_out) return -1;

    dx = x2_true - x1_true;
    dy = y2_raw - y1_raw;

    if (fabs(dy) < 1e-30) return -1; /* Points would produce infinite gain */

    model_out->sensitivity = dx / dy;  /* Engineering units per raw unit */
    model_out->offset = x1_true - model_out->sensitivity * y1_raw;
    model_out->r_squared = 1.0;  /* Perfect fit through two points */
    model_out->residuals_rms = 0.0;

    return 0;
}

/**
 * @brief Calibration table lookup with linear interpolation
 *
 * Performs binary search to locate the interval containing y_raw,
 * then linear interpolation within that interval.
 *
 * Assumptions:
 *   - x_table is sorted ascending (monotonic calibration)
 *   - y_table corresponds element-by-element to x_table
 *   - n_table >= 2
 *
 * Complexity: O(log n) for search, O(1) for interpolation
 *
 * @param x_table   True reference values [n_table]
 * @param y_table   Corresponding raw sensor readings [n_table]
 * @param n_table   Number of calibration points
 * @param y_raw     Raw sensor reading to convert
 * @return          Interpolated engineering units value
 */
double ss_cal_table_lookup(const double *x_table, const double *y_table,
                           size_t n_table, double y_raw) {
    size_t lo, hi, mid;

    if (!x_table || !y_table || n_table < 2) return y_raw;

    /* Clamp to endpoints */
    if (y_raw <= y_table[0]) return x_table[0];
    if (y_raw >= y_table[n_table - 1]) return x_table[n_table - 1];

    /* Binary search for interval */
    lo = 0;
    hi = n_table - 1;
    while (hi - lo > 1) {
        mid = (lo + hi) / 2;
        if (y_table[mid] <= y_raw) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    /* Linear interpolation between lo and hi */
    {
        double dy = y_table[hi] - y_table[lo];
        double dx = x_table[hi] - x_table[lo];
        if (fabs(dy) < 1e-30) return x_table[lo];
        return x_table[lo] + (y_raw - y_table[lo]) * dx / dy;
    }
}

/**
 * @brief Generate calibration lookup table from polynomial model
 *
 * Evaluates the polynomial at n_points evenly spaced across [x_min, x_max],
 * producing (x, y) pairs suitable for use with ss_cal_table_lookup.
 *
 * This is useful for converting a compact polynomial model into a
 * lookup table for fast, deterministic-latency interpolation in
 * real-time embedded systems.
 *
 * @param poly      Source polynomial model
 * @param n_points  Number of table entries to generate (>=2)
 * @param x_table   Output reference values [n_points]
 * @param y_table   Output sensor readings [n_points]
 * @return          0 on success, -1 on error
 */
int ss_cal_build_lookup_table(const ss_polynomial_model_t *poly,
                              size_t n_points,
                              double *x_table, double *y_table) {
    size_t i;

    if (!poly || !x_table || !y_table || n_points < 2) return -1;

    for (i = 0; i < n_points; i++) {
        double frac = (double)i / (double)(n_points - 1);
        double x_ref = poly->x_min + frac * (poly->x_max - poly->x_min);
        x_table[i] = x_ref;
        y_table[i] = ss_cal_polynomial_apply(poly, x_ref);
    }

    return 0;
}

/* =========================================================================
 * L4: Propagation of Uncertainty (JCGM 100:2008 GUM)
 * ========================================================================= */

/**
 * @brief Combined standard uncertainty — law of propagation
 *
 * For a measurement model y = f(x_1, ..., x_n) with uncorrelated
 * input quantities, the combined standard uncertainty u_c(y) is:
 *
 *   u_c^2(y) = sum_{i=1}^{n} (df/dx_i)^2 * u^2(x_i)
 *
 * where df/dx_i are the sensitivity coefficients and u(x_i) are
 * the standard uncertainties of the inputs.
 *
 * This is the fundamental formula from the GUM (JCGM 100:2008,
 * Equation 10 in Section 5.1.2).
 *
 * @param sensitivity_coeffs  Partial derivatives df/dx_i [n_inputs]
 * @param uncertainties       Standard uncertainties u(x_i) [n_inputs]
 * @param n_inputs            Number of input quantities
 * @return                    Combined standard uncertainty u_c(y)
 */
double ss_cal_combined_uncertainty(const double *sensitivity_coeffs,
                                   const double *uncertainties,
                                   size_t n_inputs) {
    size_t i;
    double sum_sq = 0.0;

    if (!sensitivity_coeffs || !uncertainties || n_inputs == 0) return 0.0;

    for (i = 0; i < n_inputs; i++) {
        double c = sensitivity_coeffs[i];
        double u = uncertainties[i];
        sum_sq += c * c * u * u;
    }

    return sqrt(sum_sq);
}

/**
 * @brief Expanded uncertainty for a given coverage factor
 *
 * U = k * u_c(y)
 *
 * Typically k = 2 for approximately 95% confidence level
 * (assuming normal distribution and large degrees of freedom).
 *
 * For small degrees of freedom, use Student's t-distribution
 * to determine k from the effective degrees of freedom (Welch-
 * Satterthwaite formula in GUM Annex G).
 *
 * @param combined_uncertainty  u_c(y)
 * @param coverage_factor       k (usually 2)
 * @return                      U = k * u_c
 */
double ss_cal_expanded_uncertainty(double combined_uncertainty,
                                   double coverage_factor) {
    return coverage_factor * combined_uncertainty;
}

/**
 * @brief Standard error of the estimate from calibration residuals
 *
 * Estimates the standard deviation of the measurement errors:
 *   sigma_hat = sqrt(RSS / (n - p))
 *
 * where:
 *   RSS = sum((y_measured_i - y_predicted_i)^2)
 *   n   = number of calibration points
 *   p   = number of model parameters (degrees of freedom consumed)
 *
 * This is the unbiased estimator of the error variance under
 * the assumption of homoscedastic normal errors.
 *
 * @param residuals  Array of residuals [n]: y_pred - y_meas
 * @param n          Number of data points
 * @param n_params   Number of model parameters (e.g., 2 for linear, d+1 for d-th order poly)
 * @return           Standard error (sigma_hat)
 */
double ss_cal_standard_error(const double *residuals, size_t n,
                             size_t n_params) {
    size_t i;
    double sum_sq = 0.0;

    if (!residuals || n <= n_params) return 0.0;

    for (i = 0; i < n; i++) {
        sum_sq += residuals[i] * residuals[i];
    }

    return sqrt(sum_sq / (double)(n - n_params));
}
