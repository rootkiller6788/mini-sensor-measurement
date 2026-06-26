#ifndef SENSOR_CALIBRATION_H
#define SENSOR_CALIBRATION_H
#include "smart_sensor.h"
#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * Sensor Calibration Module
 *
 * Implements calibration algorithms for smart sensors:
 *   L5: Linear regression, polynomial regression, Steinhart-Hart,
 *       thermocouple ITS-90 inversion, two-point calibration,
 *       multi-point calibration, temperature compensation
 *   L4: Least-squares estimation (Gauss-Markov theorem),
 *       Propagation of uncertainty
 *   L6: Thermistor calibration, thermocouple calibration,
 *       pressure sensor linearization
 *
 * Reference:
 *   JCGM 100:2008 (GUM) — Evaluation of measurement data
 *   Press, W.H. et al. (2007), "Numerical Recipes", 3rd ed., Cambridge
 *   NIST Monograph 175 — ITS-90 Thermocouple Database
 *   Steinhart & Hart (1968), Deep Sea Research 15(4):497-503
 *
 * Course Mapping:
 *   MIT 6.003 — Signals and Systems (least-squares estimation)
 *   Berkeley EE123 — Digital Signal Processing (numerical methods)
 *   ETH 227-0427 — Signal Processing (parameter estimation)
 * ========================================================================= */

/* =========================================================================
 * L5: Linear Regression Calibration
 * ========================================================================= */

/**
 * @brief Perform linear regression: y = a*x + b
 *
 * Implements ordinary least-squares (OLS) linear regression using
 * the closed-form solution:
 *   a = (n*sum(xy) - sum(x)*sum(y)) / (n*sum(x^2) - sum(x)^2)
 *   b = (sum(y) - a*sum(x)) / n
 *
 * This is the Gauss-Markov best linear unbiased estimator (BLUE)
 * under the assumption of homoscedastic uncorrelated errors.
 *
 * Complexity: O(n) time, O(1) space (single pass with running sums)
 *
 * @param x          Input (known reference) values [n]
 * @param y          Output (sensor reading) values [n]
 * @param n          Number of calibration points
 * @param model_out  Output linear model with computed parameters
 * @return 0 on success, -1 on invalid input (n<2 or null pointers)
 */
int ss_cal_linear_regression(const double *x, const double *y,
                             size_t n, ss_linear_model_t *model_out);

/**
 * @brief Apply linear model: y = sensitivity * x + offset
 * @param model  Linear calibration model
 * @param x      Input value (sensor reading)
 * @return       Calibrated output in engineering units
 */
double ss_cal_linear_apply(const ss_linear_model_t *model, double x);

/**
 * @brief Inverse linear model: x = (y - offset) / sensitivity
 * @param model  Linear calibration model
 * @param y      Calibrated value (engineering units)
 * @return       Sensor reading value
 */
double ss_cal_linear_inverse(const ss_linear_model_t *model, double y);

/* =========================================================================
 * L5: Polynomial Regression Calibration
 * ========================================================================= */

/**
 * @brief Perform polynomial regression: y = a_n*x^n + ... + a_1*x + a_0
 *
 * Uses the normal equations approach: (X^T * X) * a = X^T * y
 * where X is the Vandermonde matrix of order+1 columns.
 * Inversion via Gaussian elimination with partial pivoting.
 *
 * The normal equations approach minimizes the residual sum of squares:
 *   RSS = sum_i (y_i - sum_j a_j * x_i^j)^2
 *
 * Complexity: O(n * d^2 + d^3) where d = order+1
 *
 * @param x          Input (known reference) values [n]
 * @param y          Output (sensor reading) values [n]
 * @param n          Number of calibration points
 * @param order      Polynomial order (1-7)
 * @param model_out  Output polynomial model
 * @return 0 on success, -1 on invalid input or singular matrix
 */
int ss_cal_polynomial_regression(const double *x, const double *y,
                                 size_t n, uint8_t order,
                                 ss_polynomial_model_t *model_out);

/**
 * @brief Apply polynomial model using Horner's method
 *
 * Horner's method: p(x) = a_0 + x*(a_1 + x*(a_2 + ... + x*(a_n)...))
 * Reduces multiplication count from O(d^2) to O(d).
 *
 * @param model  Polynomial calibration model
 * @param x      Input value (sensor reading)
 * @return       Calibrated output in engineering units
 */
double ss_cal_polynomial_apply(const ss_polynomial_model_t *model, double x);

/**
 * @brief Inverse polynomial model using Newton-Raphson iteration
 *
 * Finds x such that poly_model(x) = y_target.
 * Iteration: x_{k+1} = x_k - f(x_k)/f'(x_k)
 *
 * @param model    Polynomial calibration model
 * @param y_target Target output value
 * @param x_guess  Initial guess for x
 * @param tol      Convergence tolerance
 * @param max_iter Maximum iterations
 * @return         Solution x, NaN if not converged
 */
double ss_cal_polynomial_inverse(const ss_polynomial_model_t *model,
                                 double y_target, double x_guess,
                                 double tol, int max_iter);

/* =========================================================================
 * L5: Steinhart-Hart Thermistor Calibration
 * ========================================================================= */

/**
 * @brief Compute Steinhart-Hart coefficients from 3 calibration points
 *
 * Given three (R, T) pairs, solves the linear system:
 *   1/T_i = A + B*ln(R_i) + C*[ln(R_i)]^3
 *
 * Three equations in three unknowns (A, B, C) solved directly
 * via Cramer's rule for the linear system in A, B, C.
 *
 * @param r1, t1  Resistance (Ohm) and Temperature (K) at point 1
 * @param r2, t2  Resistance (Ohm) and Temperature (K) at point 2
 * @param r3, t3  Resistance (Ohm) and Temperature (K) at point 3
 * @param model   Output Steinhart-Hart model
 * @return 0 on success, -1 on degenerate points
 */
int ss_cal_steinhart_hart_3pt(double r1, double t1,
                              double r2, double t2,
                              double r3, double t3,
                              ss_steinhart_hart_t *model);

/**
 * @brief Compute simplified Beta-parameter model from 2 points
 *
 * Beta model: R(T) = R0 * exp(Beta * (1/T - 1/T0))
 * Inverse:    T = 1 / (1/T0 + (1/Beta)*ln(R/R0))
 *
 * @param r0, t0  Reference resistance at reference temperature
 * @param r1, t1  Second calibration point
 * @param model   Output model with beta_value set
 * @return 0 on success
 */
int ss_cal_thermistor_beta(double r0, double t0,
                           double r1, double t1,
                           ss_steinhart_hart_t *model);

/**
 * @brief Convert thermistor resistance to temperature
 *
 * Full Steinhart-Hart equation:
 *   T_K = 1 / (A + B*ln(R) + C*[ln(R)]^3)
 *
 * @param model       Steinhart-Hart model
 * @param resistance  Thermistor resistance in Ohms
 * @return            Temperature in Kelvin
 */
double ss_cal_steinhart_hart_temp(const ss_steinhart_hart_t *model,
                                  double resistance);

/**
 * @brief Convert temperature back to thermistor resistance
 *
 * Solves: 1/T = A + B*ln(R) + C*[ln(R)]^3 for R
 * Uses Newton-Raphson iteration on ln(R).
 *
 * @param model  Steinhart-Hart model
 * @param temp_k Temperature in Kelvin
 * @return       Resistance in Ohms
 */
double ss_cal_steinhart_hart_resistance(const ss_steinhart_hart_t *model,
                                        double temp_k);

/* =========================================================================
 * L5: Thermocouple ITS-90 Calibration
 * ========================================================================= */

/**
 * @brief Initialize ITS-90 coefficients for a thermocouple type
 *
 * Loads NIST Monograph 175 polynomial coefficients for the specified
 * thermocouple type (K, J, T, E, N, R, S, B).
 *
 * @param model          Output model structure
 * @param tc_type        0=K, 1=J, 2=T, 3=E, 4=N, 5=R, 6=S, 7=B
 * @param cold_junction  Reference junction temperature in degC
 * @return 0 on success, -1 on invalid type
 */
int ss_cal_thermocouple_init(ss_thermocouple_model_t *model,
                             uint8_t tc_type, double cold_junction);

/**
 * @brief Convert thermocouple EMF to temperature
 *
 * Uses the NIST ITS-90 inverse polynomial:
 *   T = sum_{i=0}^{n} c_i * E^i
 * where E is the thermoelectric voltage in mV.
 *
 * Cold-junction compensation:
 *   E_actual = E_measured + E_cj
 * where E_cj = forward model(cold_junction_temp)
 *
 * @param model    Thermocouple model with cold-junction temp
 * @param emf_mv   Measured EMF in millivolts
 * @return         Temperature in degC
 */
double ss_cal_thermocouple_to_temp(const ss_thermocouple_model_t *model,
                                   double emf_mv);

/**
 * @brief Compute cold-junction compensation voltage
 *
 * Uses the forward polynomial to compute EMF at the cold junction
 * temperature, which is then added to the measured EMF.
 *
 * @param model        Thermocouple model
 * @param cj_temp_c    Cold junction temperature in degC
 * @return             EMF correction in mV
 */
double ss_cal_thermocouple_cj_comp(const ss_thermocouple_model_t *model,
                                   double cj_temp_c);

/* =========================================================================
 * L5: Multi-Point Calibration and Interpolation
 * ========================================================================= */

/**
 * @brief Two-point calibration (zero and span adjustment)
 *
 * Given two known reference points, computes linear correction
 * to adjust zero offset and span (gain).
 *
 * y_corrected = (y_raw - offset) * (span_ref / span_measured)
 *
 * @param x1, y1_raw  First reference: true value, raw reading
 * @param x2, y2_raw  Second reference: true value, raw reading
 * @param model_out    Output linear model
 * @return 0 on success, -1 if x1 == x2
 */
int ss_cal_two_point(double x1_true, double y1_raw,
                     double x2_true, double y2_raw,
                     ss_linear_model_t *model_out);

/**
 * @brief Multi-point piecewise-linear calibration table lookup
 *
 * Uses binary search to find the correct segment in a calibration
 * lookup table, then performs linear interpolation between the
 * two nearest calibration points.
 *
 * Complexity: O(log n_table) for binary search
 *
 * @param x_table     Known reference values [n_table] (sorted ascending)
 * @param y_table     Corresponding sensor readings [n_table]
 * @param n_table     Number of table entries
 * @param y_raw       Raw sensor reading
 * @return            Interpolated engineering units value
 */
double ss_cal_table_lookup(const double *x_table, const double *y_table,
                           size_t n_table, double y_raw);

/**
 * @brief Build a calibration lookup table from polynomial model
 *
 * Generates (x, y) pairs by evaluating the polynomial at evenly
 * spaced points across the valid range.
 *
 * @param poly      Polynomial model
 * @param n_points  Number of table points to generate
 * @param x_table   Output: reference values [n_points]
 * @param y_table   Output: sensor readings [n_points]
 * @return 0 on success
 */
int ss_cal_build_lookup_table(const ss_polynomial_model_t *poly,
                              size_t n_points,
                              double *x_table, double *y_table);

/* =========================================================================
 * L4: Propagation of Uncertainty (JCGM 100:2008)
 * ========================================================================= */

/**
 * @brief Compute combined standard uncertainty via propagation
 *
 * For uncorrelated input quantities:
 *   u_c^2(y) = sum_i (df/dx_i)^2 * u^2(x_i)
 *
 * This implements the "law of propagation of uncertainty" from
 * the GUM (JCGM 100:2008, Section 5).
 *
 * @param sensitivity_coeffs  Partial derivatives df/dx_i [n_inputs]
 * @param uncertainties       Standard uncertainties u(x_i) [n_inputs]
 * @param n_inputs            Number of input quantities
 * @return                    Combined standard uncertainty u_c(y)
 */
double ss_cal_combined_uncertainty(const double *sensitivity_coeffs,
                                   const double *uncertainties,
                                   size_t n_inputs);

/**
 * @brief Compute expanded uncertainty U = k * u_c
 *
 * Coverage factor k=2 corresponds to approximately 95% confidence
 * interval for a normal distribution.
 *
 * @param combined_uncertainty  u_c (combined standard uncertainty)
 * @param coverage_factor       k (typically 2 for 95% CI)
 * @return                      U (expanded uncertainty)
 */
double ss_cal_expanded_uncertainty(double combined_uncertainty,
                                   double coverage_factor);

/**
 * @brief Evaluate calibration uncertainty from residuals
 *
 * Computes the standard error of the estimate:
 *   sigma_est = sqrt(RSS / (n - p))
 * where RSS = sum(residual_i^2), n = number of points,
 * p = number of model parameters.
 *
 * @param residuals  Array of (y_predicted - y_measured) [n]
 * @param n          Number of calibration points
 * @param n_params   Number of model parameters
 * @return           Standard error of the estimate
 */
double ss_cal_standard_error(const double *residuals, size_t n,
                             size_t n_params);

#ifdef __cplusplus
}
#endif
#endif /* SENSOR_CALIBRATION_H */
