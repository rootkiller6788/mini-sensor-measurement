/**
 * biosensor_calibration.h — Calibration methods and data reduction
 *
 * L5 Algorithms: Linear regression, 4PL/5PL nonlinear fitting,
 *   weighted least squares, standard addition, Grubbs outlier test.
 * L6 Canonical Problems: ELISA calibration, ion-selective electrode
 *   standard addition, sensor aging compensation.
 * L8 Advanced Topics: Multivariate calibration (PLS, PCR),
 *   transfer calibration across sensor batches.
 *
 * References:
 *   - Finney, "Statistical Method in Biological Assay" (1978)
 *   - Rodbard & Hutt, "Radioimmunoassay Data Processing" (1974)
 *   - ISO 11095 "Linear calibration using reference materials"
 *   - Grubbs, Technometrics 11 (1969) 1-21
 *
 * L9 alignment: ETH 227-0427, Caltech BE/APh 161
 */

#ifndef BIOSENSOR_CALIBRATION_H
#define BIOSENSOR_CALIBRATION_H

#include "biosensor_types.h"

/* ============================================================================
 * L5: Linear regression calibration
 * ============================================================================ */

/**
 * Ordinary least squares linear calibration (L5).
 *
 * y = a + b · x
 *
 * Minimizes Σ (y_i - (a + b·x_i))².
 * Returns slope (b, sensitivity), intercept (a, background),
 * R², and standard errors.
 *
 * L5 Algorithm: OLS regression with full uncertainty propagation.
 *
 * @param standards     Calibration standards data (x = conc, y = signal)
 * @param result        [out] Fitted LinearCalibration structure
 * @return              R² goodness of fit
 */
double linear_calibration_fit(const CalibrationStandards *standards,
                              LinearCalibration *result);

/**
 * Weighted least squares calibration (L5).
 *
 * Minimizes Σ w_i · (y_i - (a + b·x_i))²
 * where w_i = 1/σ_i² (inverse variance weighting).
 *
 * Preferred when heteroscedasticity is present (common in biosensors
 * where signal variance increases with concentration).
 *
 * @param standards     Calibration standards with signal_stds
 * @param result        [out] Weighted LinearCalibration structure
 * @return              Weighted R²
 */
double weighted_linear_calibration_fit(const CalibrationStandards *standards,
                                       LinearCalibration *result);

/**
 * Predict concentration from calibration with prediction interval (L5).
 *
 * x_pred = (y - intercept) / slope
 *
 * With 95% prediction interval accounting for calibration uncertainty:
 *   PI = t_0.95,ν · s_yx/slope · √(1 + 1/N + (x_pred - x̄)²/Σ(x_i - x̄)²)
 *
 * @param lc           Fitted calibration
 * @param signal       Measured sample signal
 * @param standards    Original calibration standards (for leverage calc)
 * @param concentration [out] Predicted concentration
 * @param lower_bound  [out] Lower 95% prediction interval
 * @param upper_bound  [out] Upper 95% prediction interval
 */
void linear_calibration_predict(const LinearCalibration *lc,
                                double signal,
                                const CalibrationStandards *standards,
                                double *concentration,
                                double *lower_bound,
                                double *upper_bound);

/* ============================================================================
 * L5: Nonlinear calibration — 4PL and 5PL Logistic Models
 * ============================================================================ */

/**
 * Fit 4-Parameter Logistic (4PL) model by iterative least squares (L5).
 *
 * Model: y = a + (d - a) / (1 + (x/c)^b)
 *
 * Uses Levenberg-Marquardt algorithm for robust convergence.
 * The 4PL is the NIH/WHO gold standard for immunoassay calibration.
 *
 * L5 Algorithm: Levenberg-Marquardt nonlinear regression.
 *
 * @param standards     Calibration standards
 * @param fit           [in/out] Initial parameter guesses, returns final fit
 * @param max_iter      Maximum LM iterations
 * @param lambda_init   Initial LM damping factor (suggest 0.01)
 * @return              Number of iterations used (0 if diverged)
 */
int fourpl_fit(const CalibrationStandards *standards,
               FourPLLogisticFit *fit,
               int max_iter,
               double lambda_init);

/**
 * Fit 5-Parameter Logistic (5PL) model (L5).
 *
 * Model: y = a + (d - a) / (1 + (x/c)^b)^g
 *
 * The asymmetry parameter g accounts for non-symmetric calibration
 * curves common in heterogeneous immunoassays.
 *
 * @param standards     Calibration standards
 * @param fit           [in/out] Initial guesses, returns final fit
 * @param max_iter      Maximum LM iterations
 * @param lambda_init   Initial LM damping factor
 * @return              Number of iterations used
 */
int fivepl_fit(const CalibrationStandards *standards,
               FivePLLogisticFit *fit,
               int max_iter,
               double lambda_init);

/**
 * Log-logistic linearization for initial 4PL parameter estimation (L5).
 *
 * The 4PL can be linearized as:
 *   log((d - y)/(y - a)) = b · log(x) - b · log(c)
 *
 * This provides excellent initial guesses for the full nonlinear fit.
 *
 * @param standards     Calibration standards
 * @param fit           [in] a and d must be set (asymptote estimates)
 *                      [out] b and c are estimated
 * @return              R² of the linearized fit
 */
double fourpl_linearized_initial_guess(const CalibrationStandards *standards,
                                       FourPLLogisticFit *fit);

/* ============================================================================
 * L5: Standard addition method
 * ============================================================================ */

/**
 * Standard addition calibration for complex sample matrices (L5).
 *
 * When the sample matrix causes signal interference (matrix effect),
 * standard addition eliminates the bias by building the calibration
 * curve within the sample itself.
 *
 * The unknown concentration is |-intercept/slope| from the standard
 * addition regression line.
 *
 * L6 Canonical Problem: Heavy metal detection in wastewater,
 * blood lead analysis with ISEs.
 *
 * @param sample_volume_ml     Volume of sample in each aliquot
 * @param spike_concentrations Added concentrations in each aliquot
 * @param spike_volumes_ml     Volume of standard spike added
 * @param total_volumes_ml     Final total volume of each aliquot
 * @param measured_signals     Measured response for each aliquot
 * @param n_aliquots           Number of aliquots (≥ 5 recommended)
 * @param original_conc        [out] Estimated original sample concentration
 * @param r_squared            [out] R² of standard addition line
 * @return                     true if fit is reliable (R² > 0.95)
 */
bool standard_addition_method(const double *sample_volume_ml,
                              const double *spike_concentrations,
                              const double *spike_volumes_ml,
                              const double *total_volumes_ml,
                              const double *measured_signals,
                              int n_aliquots,
                              double *original_conc,
                              double *r_squared);

/* ============================================================================
 * L5: Outlier detection and quality control
 * ============================================================================ */

/**
 * Grubbs' test for outlier detection (L5).
 *
 * G = |x_suspect - x̄| / s
 *
 * Compares against critical G values from Student's t-distribution:
 *   G_crit = (N-1)/√N · √(t²/(N-2+t²))
 *
 * L5 Algorithm: Statistical outlier detection for biosensor QC.
 *
 * @param values            Array of measurements
 * @param n_values          Number of replicates
 * @param alpha             Significance level (typically 0.05)
 * @param outlier_index     [out] Index of the outlier, or -1 if none found
 * @return                  true if an outlier was detected
 */
bool grubbs_outlier_test(const double *values,
                         int n_values,
                         double alpha,
                         int *outlier_index);

/**
 * Dixon's Q-test for small sample outlier detection (L5).
 *
 * Q = |x_suspect - x_nearest| / (x_max - x_min)
 *
 * Recommended for N = 2-7 replicates. Simpler than Grubbs.
 *
 * @param values            Sorted array of measurements
 * @param n_values          Number of replicates (2 ≤ N ≤ 7)
 * @param alpha             Significance level (typically 0.05)
 * @param is_high_outlier   true to test maximum, false to test minimum
 * @return                  true if the extreme value is an outlier
 */
bool dixon_q_test(const double *values,
                  int n_values,
                  double alpha,
                  bool is_high_outlier);

/**
 * Compute coefficient of variation (CV%) for replicate measurements (L5).
 *
 * CV = (s / x̄) × 100%
 *
 * Standard QC metric in clinical biosensor labs. FDA requires
 * CV < 15% for bioanalytical method validation (except at LLOQ
 * where CV < 20% is acceptable).
 *
 * @param values            Replicate measurements
 * @param n_values          Number of replicates
 * @return                  CV as percentage (0 to >100)
 */
double coefficient_of_variation(const double *values, int n_values);

/* ============================================================================
 * L6: Sensor drift and aging compensation
 * ============================================================================ */

/**
 * Two-point calibration adjustment for sensor aging (L6).
 *
 * When a biosensor ages, sensitivity decreases. Two-point
 * recalibration with low and high QC standards adjusts the
 * calibration curve:
 *
 *   slope_new = (sig_high - sig_low) / (conc_high - conc_low)
 *   intercept_new = sig_low - slope_new · conc_low
 *
 * L6 Canonical Problem: Daily QC recalibration of blood gas analyzers.
 *
 * @param conc_low          Concentration of low QC standard
 * @param sig_low           Measured signal at low QC
 * @param conc_high         Concentration of high QC standard
 * @param sig_high          Measured signal at high QC
 * @param old_cal           Existing (aged) calibration
 * @param new_cal           [out] Adjusted calibration
 */
void two_point_calibration_adjust(double conc_low, double sig_low,
                                  double conc_high, double sig_high,
                                  const LinearCalibration *old_cal,
                                  LinearCalibration *new_cal);

/**
 * Estimate sensor lifetime remaining based on sensitivity decay (L6).
 *
 * Exponential decay model for sensor sensitivity:
 *   S(t) = S₀ · exp(-t / τ)
 *
 * The remaining lifetime factor is S_current / S_min where S_min
 * is the minimum acceptable sensitivity.
 *
 * @param initial_sensitivity    S₀: factory sensitivity
 * @param current_sensitivity    S_current: measured at latest QC
 * @param days_in_use            Days since first use
 * @param min_sensitivity        S_min: minimum acceptable sensitivity
 * @return                       Estimated remaining days before S < S_min
 */
double estimate_sensor_remaining_life(double initial_sensitivity,
                                      double current_sensitivity,
                                      int days_in_use,
                                      double min_sensitivity);

/* ============================================================================
 * L6: Calibration quality metrics
 * ============================================================================ */

/**
 * Compute the Mandel's fitting test statistic (L5: Lack-of-fit).
 *
 * Compares the linear model error against the pure error from
 * replicates. If F > F_crit, the linear model is inadequate and
 * a nonlinear model should be considered.
 *
 * @param standards     Calibration standards with replicates
 * @param lc            Fitted linear calibration
 * @param f_critical    F_crit at chosen α (e.g., =3.5 for typical)
 * @return              true if linear model is adequate
 */
bool mandel_fitting_test(const CalibrationStandards *standards,
                         const LinearCalibration *lc,
                         double f_critical);

/**
 * Compute Durbin-Watson statistic for autocorrelation in residuals (L5).
 *
 * DW = Σ(e_i - e_{i-1})² / Σ e_i²
 *
 * DW ≈ 2: no autocorrelation (independent residuals).
 * DW < 1: positive autocorrelation (suspicious calibration).
 *
 * @param residuals     Residuals from calibration fit
 * @param n_points      Number of calibration standards
 * @param dw_statistic  [out] Durbin-Watson statistic
 * @return              true if residuals appear independent (DW ≈ 2)
 */
bool durbin_watson_test(const double *residuals,
                        int n_points,
                        double *dw_statistic);

/* ============================================================================
 * L8: Advanced — Transfer calibration between sensor batches
 * ============================================================================ */

/**
 * Piecewise direct standardization (PDS) for calibration transfer (L8).
 *
 * Transfers a calibration model from a master sensor (batch A)
 * to a slave sensor (batch B) using a small set of transfer standards.
 *
 * For each wavelength/channel j on slave, a window of master channels
 * [j-k, j+k] is used in a local PLS model to predict the slave response.
 *
 * L8 Advanced Topic: Chemometric calibration transfer.
 *
 * @param master_responses  Master sensor responses to transfer standards
 * @param slave_responses   Slave sensor responses to same standards
 * @param n_standards       Number of transfer standards
 * @param n_channels        Number of sensor channels/wavelengths
 * @param window_half_width Half-width of PDS window (k)
 * @param transfer_matrix   [out] PDS transfer matrix (n_channels × n_channels)
 * @return                  RMS prediction error of the transfer
 */
double pds_calibration_transfer(const double *master_responses,
                                const double *slave_responses,
                                int n_standards,
                                int n_channels,
                                int window_half_width,
                                double *transfer_matrix);

#endif /* BIOSENSOR_CALIBRATION_H */
