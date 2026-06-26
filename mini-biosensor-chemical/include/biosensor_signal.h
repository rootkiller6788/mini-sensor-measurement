/**
 * biosensor_signal.h — Signal processing and data analysis for biosensors
 *
 * L5 Algorithms: LOD/LOQ estimation, noise characterization, baseline
 *   correction, drift compensation, peak detection, Savitzky-Golay filtering,
 *   multivariate data analysis.
 * L6 Canonical Problems: Sensor baseline drift removal, signal-to-noise
 *   optimization, multi-electrode array denoising.
 *
 * References:
 *   - IUPAC Orange Book, "Compendium of Analytical Nomenclature" (1997)
 *   - Savitzky & Golay, Analytical Chemistry 36 (1964) 1627-1639
 *   - Massart et al., "Handbook of Chemometrics and Qualimetrics" (1998)
 *
 * L9 alignment: ETH 227-0427 Signal Processing, Michigan EECS 351 DSP
 */

#ifndef BIOSENSOR_SIGNAL_H
#define BIOSENSOR_SIGNAL_H

#include "biosensor_types.h"

/* ============================================================================
 * L5: Limit of Detection / Limit of Quantification estimation
 * ============================================================================ */

/**
 * Compute Limit of Detection (LOD) by the 3σ method (L5, IUPAC).
 *
 * LOD = 3.3 · σ_blank / S
 *
 * Where σ_blank is the standard deviation of blank (zero-analyte)
 * measurements, and S is the analytical sensitivity (calibration slope).
 *
 * The factor 3.3 (= √2 · t_0.05,∞ ≈ 2 · 1.645) accounts for the
 * 95% confidence interval assuming normally distributed blank errors.
 *
 * IUPAC recommendations: LOD = x̄_blank + 3·σ_blank as signal,
 * then converted to concentration via the calibration function.
 *
 * @param blank_std          σ_blank: standard deviation of N≥20 blank replicates
 * @param sensitivity        S: calibration curve slope [signal/concentration]
 * @return                   LOD [concentration units]
 */
double lod_three_sigma(double blank_std, double sensitivity);

/**
 * Compute Limit of Quantification (LOQ) by the 10σ method (L5).
 *
 * LOQ = 10 · σ_blank / S
 *
 * The signal at LOQ has RSD ≈ 10%, which is the threshold for
 * quantitative (vs qualitative) measurement. IUPAC recommends
 * LOQ = 10σ for standard analytical methods.
 *
 * @param blank_std          σ_blank
 * @param sensitivity        S [signal/concentration]
 * @return                   LOQ [concentration units]
 */
double loq_ten_sigma(double blank_std, double sensitivity);

/**
 * Compute LOD from calibration curve using ISO 11843 method (L5).
 *
 * LOD = (t_0.95,ν · σ_residual / S) · √(1 + 1/N + x̄²/Σ(xᵢ-x̄)²)
 *
 * This is the more rigorous ISO-recommended method that accounts
 * for both intercept uncertainty and calibration leverage.
 *
 * @param concentrations     Calibration standard concentrations
 * @param residuals          Residuals = measured - predicted signal
 * @param sensitivity        S: calibration slope
 * @param n_standards        Number of calibration standards
 * @param n_total            Total measurements (n_standards × replicates)
 * @return                   ISO 11843 LOD [concentration units]
 */
double lod_iso_11843(const double *concentrations,
                     const double *residuals,
                     double sensitivity,
                     int n_standards,
                     int n_total);

/* ============================================================================
 * L5: Baseline and drift correction
 * ============================================================================ */

/**
 * Asymmetric Least Squares baseline correction (Eilers, 2003) — L5.
 *
 * Minimizes: Σ w_i · (y_i - z_i)² + λ · Σ (Δ²z_i)²
 *
 * Where w_i = p for y_i > z_i (positive residuals, signal) and
 * w_i = 1-p for y_i ≤ z_i (negative residuals, baseline).
 *
 * Widely used in chromatography, Raman spectroscopy, and biosensor
 * drift removal. The asymmetry parameter p (0.001 to 0.1) determines
 * the "stiffness" toward positive peaks.
 *
 * L5 Algorithm: Penalized least squares baseline estimation.
 *
 * @param signal            Raw sensor signal array
 * @param n_points          Number of data points
 * @param lambda_smooth     Second-derivative smoothing penalty (10² to 10⁹)
 * @param p_asymmetry       Asymmetry parameter (0.001 = stiff, 0.1 = flexible)
 * @param max_iterations    Maximum ALS iterations (typically 5-20)
 * @param baseline          [out] Estimated baseline (must be pre-allocated)
 * @return                  Number of ALS iterations performed
 */
int als_baseline_correction(const double *signal,
                            int n_points,
                            double lambda_smooth,
                            double p_asymmetry,
                            int max_iterations,
                            double *baseline);

/**
 * Simple moving-average baseline estimation for sensor drift (L5).
 *
 * Computes baseline as the running minimum over a sliding window,
 * followed by interpolation. Suitable for slow sensor drift.
 *
 * @param signal            Raw signal array
 * @param n_points          Number of data points
 * @param window_size       Sliding window width for minimum tracking
 * @param baseline          [out] Estimated baseline
 */
void moving_minimum_baseline(const double *signal,
                             int n_points,
                             int window_size,
                             double *baseline);

/**
 * Exponential drift compensation (L5).
 *
 * Many biosensors exhibit exponential drift toward equilibrium:
 *   s_corrected = s_raw - A·(1 - exp(-t/τ))
 *
 * This function models and removes the exponential drift component.
 *
 * @param signal            Raw signal array
 * @param timestamps        Time array [sec]
 * @param n_points          Number of data points
 * @param drift_amplitude   A: asymptotic drift magnitude
 * @param drift_time_constant τ: drift time constant [sec]
 * @param corrected         [out] Drift-corrected signal (pre-allocated)
 */
void exponential_drift_correction(const double *signal,
                                  const double *timestamps,
                                  int n_points,
                                  double drift_amplitude,
                                  double drift_time_constant,
                                  double *corrected);

/* ============================================================================
 * L5: Smoothing and filtering
 * ============================================================================ */

/**
 * Savitzky-Golay smoothing filter (L5).
 *
 * Fits a polynomial of degree 'order' within a sliding window of
 * size 'window' (odd) using least squares, then evaluates the
 * polynomial at the center point. Preserves higher moments of the
 * signal better than moving-average filters.
 *
 * Reference: Savitzky & Golay, Anal. Chem. 36 (1964) 1627-1639.
 *
 * L5 Algorithm: Polynomial convolution smoothing.
 *
 * @param signal            Input signal array
 * @param n_points          Number of data points
 * @param window            Window size (odd, e.g., 5, 7, 9, 11)
 * @param order             Polynomial order (typically 2-4)
 * @param smoothed          [out] Smoothed signal (pre-allocated)
 * @return                  true if parameters are valid, false otherwise
 */
bool savitzky_golay_smooth(const double *signal,
                           int n_points,
                           int window,
                           int order,
                           double *smoothed);

/**
 * Median filter for spike removal in biosensor signals (L5).
 *
 * Replaces each point with the median of a sliding window.
 * Excellent for removing impulse noise without blurring edges.
 *
 * @param signal            Input signal array
 * @param n_points          Number of data points
 * @param window            Window size (odd, typically 3-7)
 * @param filtered          [out] Filtered signal (pre-allocated)
 */
void median_filter_spike_removal(const double *signal,
                                 int n_points,
                                 int window,
                                 double *filtered);

/**
 * Compute signal-to-noise ratio from replicate blank measurements (L5).
 *
 * SNR = μ_signal / σ_blank
 *
 * @param signal_mean        Mean signal at a given concentration
 * @param blank_std          Standard deviation of blank replicates
 * @return                   SNR (dimensionless)
 */
double compute_snr(double signal_mean, double blank_std);

/* ============================================================================
 * L6: Electrochemical sensor signal processing
 * ============================================================================ */

/**
 * Detect and quantify peaks in cyclic voltammogram (L6).
 *
 * Scans for local maxima (oxidation peaks) and local minima
 * (reduction peaks) in a voltammogram I vs V curve.
 *
 * L6 Canonical Problem: Cyclic voltammetry (CV) is the primary
 * characterization method for electrochemical biosensors.
 *
 * @param current           Measured current array [A]
 * @param voltage           Applied potential array [V]
 * @param n_points          Number of data points
 * @param peak_voltages     [out] Array of peak potentials (min size: n_points/2)
 * @param peak_currents     [out] Array of peak currents
 * @param max_peaks         Maximum number of peaks to detect
 * @return                  Number of peaks found
 */
int detect_voltammetry_peaks(const double *current,
                             const double *voltage,
                             int n_points,
                             double *peak_voltages,
                             double *peak_currents,
                             int max_peaks);

/**
 * Integrate current over time for coulometric analysis (L6).
 *
 * Q = ∫ I(t) dt   (trapezoidal integration)
 *
 * The total charge passed relates to the number of moles
 * reacted via Faraday's law: Q = n · F · N.
 *
 * @param current           Current array [A]
 * @param timestamps        Time array [s]
 * @param n_points          Number of data points
 * @return                  Total charge [C]
 */
double coulometric_charge(const double *current,
                          const double *timestamps,
                          int n_points);

/* ============================================================================
 * L7: Sensor array / electronic nose preprocessing
 * ============================================================================ */

/**
 * Normalize sensor array response pattern (L7).
 *
 * Fractional baseline manipulation:
 *   R_norm = (R - R_0) / R_0 = ΔR / R_0
 *
 * This is the standard preprocessing for electronic nose data,
 * compensating for different baseline resistances across sensors.
 *
 * L7 Application: Pattern recognition for gas/odor classification.
 *
 * @param raw_responses      Raw sensor resistances [Ω] (n_sensors long)
 * @param baseline_resistances Baseline R_0 [Ω] (n_sensors long)
 * @param normalized         [out] Normalized response (ΔR/R_0) (pre-allocated)
 * @param n_sensors          Number of sensors in array
 */
void sensor_array_normalize(const double *raw_responses,
                            const double *baseline_resistances,
                            double *normalized,
                            int n_sensors);

/**
 * Compute Euclidean distance between two sensor array patterns (L7).
 *
 * d = √(Σ (p_i - q_i)²)
 *
 * Used in k-NN classification of odors/volatiles in electronic nose systems.
 *
 * @param pattern_a          First normalized pattern
 * @param pattern_b          Second normalized pattern
 * @param n_sensors          Dimensionality of pattern
 * @return                   Euclidean distance
 */
double sensor_array_euclidean_distance(const double *pattern_a,
                                       const double *pattern_b,
                                       int n_sensors);

#endif /* BIOSENSOR_SIGNAL_H */
