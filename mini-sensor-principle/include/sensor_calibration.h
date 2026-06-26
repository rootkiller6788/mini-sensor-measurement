/**
 * @file sensor_calibration.h
 * @brief Sensor calibration algorithms — least-squares, polynomial fit, temperature compensation
 *
 * Level: L5 — Algorithms/Methods
 *
 * Sensor calibration is the process of determining the relationship
 * between sensor output and the true value of the measurand,
 * usually by comparison with a reference standard.
 *
 * Key Algorithms:
 *
 * 1. Linear Least-Squares (one-point and two-point calibration):
 *    y = m·x + b, where (m, b) found by minimizing Σ(y_i - m·x_i - b)²
 *
 * 2. Polynomial Least-Squares (multi-point):
 *    y = a₀ + a₁·x + a₂·x² + ... + aₙ·xⁿ
 *    Solved via normal equations or QR decomposition.
 *
 * 3. Weighted Least-Squares:
 *    Minimizes Σ w_i·(y_i - f(x_i))² when uncertainty varies with x.
 *
 * 4. Temperature Compensation:
 *    Corrects for sensor output variation with temperature.
 *    Separates measurand response from thermal response.
 *
 * 5. Cross-Sensitivity Compensation:
 *    Corrects for sensor response to interfering quantities
 *    (e.g., humidity sensor sensitive to temperature).
 *
 * Calibration Traceability Chain (per ISO/IEC 17025):
 *   SI Definition → Primary Standard → Secondary Standard →
 *   Working Standard → Sensor Under Test
 *
 * References:
 *   - ISO/IEC 17025:2017 "General requirements for the competence of
 *     testing and calibration laboratories"
 *   - JCGM 100:2008 (GUM) "Evaluation of measurement data"
 *   - Bevington & Robinson, "Data Reduction and Error Analysis",
 *     3rd Ed, McGraw-Hill, 2002
 */

#ifndef SENSOR_CALIBRATION_H
#define SENSOR_CALIBRATION_H

#include <stddef.h>

/* ──── L5: Linear Least-Squares Calibration ──── */

/**
 * @brief One-point (offset) calibration.
 *
 * Compares one reference point y_ref at known input x_ref to raw output y_raw.
 * Assumes linear sensor with known sensitivity.
 *
 * Corrected output: y_cal = y_raw + (y_ref - y_raw_at_x_ref)
 *
 * This is the simplest calibration — corrects offset only.
 */
int calibration_offset(const double *x_raw_before, const double *x_ref,
                        size_t n_cal_points, double *offset);

/**
 * @brief Two-point (gain + offset / slope + intercept) calibration.
 *
 * Uses two reference points to determine both sensitivity and offset.
 *
 * m = (y_ref₂ - y_ref₁) / (y_raw₂ - y_raw₁)
 * b = y_ref₁ - m·y_raw₁
 *
 * Then: y_cal = m·y_raw + b
 *
 * Accuracy is limited by linearity of the sensor between the two points.
 */
int calibration_two_point(double x_raw_1, double x_ref_1,
                           double x_raw_2, double x_ref_2,
                           double *slope, double *intercept);

/**
 * @brief Multi-point linear least-squares regression.
 *
 * Given n calibration pairs (x_raw_i, y_ref_i), find (m, b) minimizing:
 *   Σ (y_ref_i - m·x_raw_i - b)²
 *
 * Solution (normal equations):
 *   m = [n·Σ(x_i·y_i) - Σx_i·Σy_i] / [n·Σ(x_i²) - (Σx_i)²]
 *   b = [Σy_i - m·Σx_i] / n
 *
 * Also computes R² (coefficient of determination) and RMS residual.
 */
int calibration_linear_ls(const double *x_raw, const double *y_ref,
                           size_t n, double *slope, double *intercept,
                           double *r_squared, double *rms_residual);

/* ──── L5: Polynomial Least-Squares Calibration ──── */

/**
 * @brief Polynomial least-squares regression (Vandermonde normal equations).
 *
 * Fits a polynomial of given order to calibration data:
 *   y = a₀ + a₁·x + a₂·x² + ... + a_order·x^order
 *
 * Solves the normal equations: (X^T·X)·a = X^T·y
 * where X is the Vandermonde matrix X_{i,j} = x_i^j.
 *
 * Uses Cholesky decomposition since X^T·X is symmetric positive definite.
 *
 * @param x input values [n]
 * @param y reference values [n]
 * @param n number of calibration points
 * @param order polynomial degree (1-5 recommended, n > order+1 required)
 * @param coeff output coefficients [order+1]
 * @param rms_residual output RMS residual
 * @return 0 on success, -1 if singular/malconditioned
 */
int calibration_polynomial_ls(const double *x, const double *y,
                               size_t n, int order, double *coeff,
                               double *rms_residual);

/**
 * @brief Compute calibration polynomial fit quality metrics.
 *
 * R² = 1 - SS_res / SS_tot
 * Adjusted R² = 1 - (1-R²)·(n-1)/(n-p-1)  where p = order
 *
 * RMS residual = √(SS_res / n)
 * Max residual = max_i |y_i - f(x_i)|
 */
int calibration_fit_metrics(const double *x, const double *y,
                             size_t n, const double *coeff, int order,
                             double *r_squared, double *adj_r_squared,
                             double *rms_residual, double *max_residual);

/* ──── L5: Weighted Least-Squares ──── */

/**
 * @brief Weighted linear least-squares.
 *
 * Minimizes Σ w_i·(y_i - m·x_i - b)²
 *
 * Solution:
 *   m = [W·Σw_i·x_i·y_i - (Σw_i·x_i)(Σw_i·y_i)] / [W·Σw_i·x_i² - (Σw_i·x_i)²]
 *   b = [Σw_i·y_i - m·Σw_i·x_i] / W
 *   where W = Σw_i
 *
 * Weights are typically w_i = 1/σ_i² where σ_i is the uncertainty of
 * reference point i. This ensures points with smaller uncertainty
 * have greater influence on the fit.
 */
int calibration_weighted_linear_ls(const double *x, const double *y,
                                    const double *weights, size_t n,
                                    double *slope, double *intercept,
                                    double *rms_residual);

/* ──── L5: Temperature Compensation ──── */

/**
 * @brief Compute temperature compensation coefficients.
 *
 * Many sensors exhibit temperature-dependent sensitivity and offset:
 *   y(T) = [S₀ + α·(T-T₀)] · x  +  [b₀ + β·(T-T₀)]
 *
 * where:
 *   S₀ = sensitivity at reference temperature T₀
 *   α  = temperature coefficient of sensitivity [ppm/°C or %/°C]
 *   b₀ = offset at T₀
 *   β  = temperature coefficient of offset [unit/°C]
 *
 * This function determines α, β from measurements at multiple temperatures.
 *
 * @param temps    temperature values [n_temps]
 * @param slopes   measured sensitivity at each temperature [n_temps]
 * @param offsets  measured offset at each temperature [n_temps]
 * @param n_temps  number of temperature calibration points
 * @param T0       reference temperature
 * @param S0_out   output: sensitivity at T0
 * @param alpha_out output: temp coefficient of sensitivity
 * @param b0_out   output: offset at T0
 * @param beta_out output: temp coefficient of offset
 * @return 0 on success
 */
int calibration_temp_compensation(const double *temps,
                                   const double *slopes,
                                   const double *offsets,
                                   size_t n_temps, double T0,
                                   double *S0_out, double *alpha_out,
                                   double *b0_out, double *beta_out);

/**
 * @brief Apply temperature compensation to a sensor reading.
 *
 * y_compensated = (y_raw - [b₀ + β·(T-T₀)]) / [S₀ + α·(T-T₀)]
 *               = (y_raw - b_eff) / S_eff
 */
double calibration_apply_temp_comp(double y_raw, double T_actual,
                                    double T0, double S0, double alpha,
                                    double b0, double beta);

/* ──── L5: Cross-Sensitivity Compensation ──── */

/**
 * @brief Two-parameter cross-sensitivity compensation.
 *
 * Sensor responds to measurand x and interference quantity z:
 *   y = S_x·x + S_z·z + b
 *
 * Given measurements with known x and z, compute compensated x:
 *   x = (y - b - S_z·z) / S_x
 *
 * This requires:
 *   1. Known cross-sensitivity S_z (from datasheet or calibration)
 *   2. Independent measurement of z (from another sensor)
 *
 * @param y_raw       raw sensor reading
 * @param z_measured  independently measured interference quantity
 * @param S_x         sensitivity to measurand
 * @param S_z         cross-sensitivity to interference
 * @param b           offset
 * @return compensated measurand value
 */
double calibration_cross_sensitivity_comp(double y_raw, double z_measured,
                                           double S_x, double S_z, double b);

/**
 * @brief Compute cross-sensitivity matrix for multi-sensor systems.
 *
 * For N sensors each responding to N measurands:
 *   y = S · x + b
 *   where S is the N×N sensitivity matrix.
 *
 * Compensation: x = S⁻¹ · (y - b)
 *
 * This is the foundation of multi-sensor data fusion at the calibration level.
 */
int calibration_cross_sensitivity_matrix(const double *y, const double *S,
                                          const double *b, int n,
                                          double *x_compensated);

/* ──── L5: Auto-Zero and Offset Correction ──── */

/**
 * @brief Simulate auto-zero (chopper stabilization) cycle.
 *
 * Phase 1 (zero): measure output with zero input → store offset
 * Phase 2 (measure): measure output with signal → subtract offset
 *
 * Auto-zeroing eliminates DC offset and low-frequency drift (1/f noise
 * below the chopping frequency).
 *
 * Effective after N cycles: offset reduced to ~1/√N of original noise.
 */
int calibration_auto_zero(const double *signal_with_offset,
                           const double *zero_readings,
                           size_t n_cycles, double *corrected_signal);

/**
 * @brief Three-point piecewise linearization.
 *
 * For sensors with known nonlinearity shape, using 3 calibration points
 * to create 2 linear segments can reduce nonlinearity error significantly
 * compared to a single straight line.
 *
 * Segment 1: (x0, y0) to (x1, y1), valid for x ∈ [x0, x1]
 * Segment 2: (x1, y1) to (x2, y2), valid for x ∈ [x1, x2]
 */
int calibration_piecewise_3point(double x0, double y0,
                                  double x1, double y1,
                                  double x2, double y2,
                                  double *m1, double *b1,
                                  double *m2, double *b2);

#endif /* SENSOR_CALIBRATION_H */
