#ifndef MEASUREMENT_CALIBRATION_H
#define MEASUREMENT_CALIBRATION_H

#include <stddef.h>

/**
 * @file measurement_calibration.h
 * @brief L2, L5: Sensor calibration, regression analysis, curve fitting
 *
 * L2 Core Concepts:
 *   - Calibration: relationship between sensor output and measurand
 *   - Traceability: chain of comparisons to SI units
 *   - Calibration curve: y = f(x) mapping measured to true
 *
 * L5 Algorithms:
 *   - Ordinary Least Squares (OLS) linear regression
 *   - Polynomial regression (up to 4th order)
 *   - Weighted Least Squares (WLS)
 *   - Exponential and logarithmic curve fitting
 *   - Goodness-of-fit: R-squared, RMSE, max residual
 *
 * Reference: ISO 11095:1996 (linear calibration), NIST/SEMATECH e-Handbook
 */

typedef enum {
    CAL_LINEAR       = 0,
    CAL_POLYNOMIAL   = 1,
    CAL_EXPONENTIAL  = 2,
    CAL_LOGARITHMIC  = 3,
    CAL_POWER_LAW    = 4,
    CAL_RATIONAL     = 5
} CalibrationModelType;

#define MAX_CAL_COEFFS 8

typedef struct {
    double x;
    double y;
} CalibrationPoint;

typedef struct {
    CalibrationModelType model_type;
    size_t               order;
    double               coefficients[MAX_CAL_COEFFS];
    size_t               num_coeffs;
    double               r_squared;
    double               rmse;
    double               max_residual;
    double               mean_residual;
    size_t               num_points;
} CalibrationResult;

/* --- L5: Regression Algorithms --- */

/** Linear least squares: y = a0 + a1 * x */
int cal_linear_fit(const CalibrationPoint *points, size_t n,
                    CalibrationResult *result);

/** Polynomial least squares: y = a0 + a1*x + ... + ak*x^k */
int cal_polynomial_fit(const CalibrationPoint *points, size_t n,
                        size_t order, CalibrationResult *result);

/** Weighted least squares with diagonal weight matrix */
int cal_weighted_poly_fit(const CalibrationPoint *points,
                          const double *weights, size_t n,
                          size_t order, CalibrationResult *result);

/** Exponential fit: y = a * exp(b * x), via log-space linearization */
int cal_exponential_fit(const CalibrationPoint *points, size_t n,
                         CalibrationResult *result);

/** Power-law fit: y = a * x^b, via log-log linearization */
int cal_power_law_fit(const CalibrationPoint *points, size_t n,
                       CalibrationResult *result);

/** Logarithmic fit: y = a + b * ln(x), requires x > 0 */
int cal_logarithmic_fit(const CalibrationPoint *points, size_t n,
                         CalibrationResult *result);

/** Evaluate calibration model at x: return f(x) */
double cal_evaluate(const CalibrationResult *result, double x);

/** Inverse evaluation: given y, find x such that f(x) = y */
int cal_inverse_evaluate(const CalibrationResult *result, double y,
                          double x0, double *x_est);

/** Compute R-squared: 1 - SS_res / SS_tot */
double cal_compute_r_squared(const CalibrationPoint *points, size_t n,
                              CalibrationResult *result);

/** Compute RMSE: sqrt(SSE / n) */
double cal_compute_rmse(const CalibrationPoint *points, size_t n,
                         CalibrationResult *result);

/** Lack-of-fit F-test */
double cal_lack_of_fit_f_test(const CalibrationPoint *points, size_t n,
                               CalibrationResult *result);

/** Prediction uncertainty per ISO 11095 */
double cal_prediction_uncertainty(const CalibrationResult *result,
                                   double x, const CalibrationPoint *points,
                                   size_t n);

#endif /* MEASUREMENT_CALIBRATION_H */
