/**
 * @file quantization.h
 * @brief Quantization theory: uniform and non-uniform quantizers,
 *        Lloyd-Max optimal quantization, dithering, noise shaping.
 *
 * Knowledge Coverage:
 *   L1 - Definitions: LSB, quantization step, quantization error,
 *        granular noise, overload distortion, quantization noise power
 *   L2 - Core Concepts: uniform vs non-uniform quantization,
 *        companding, mid-tread vs mid-rise quantizer
 *   L3 - Mathematical Structures: Quantization error as additive white
 *        noise model, probability density function of error (uniform),
 *        correlation structure of quantization error sequence
 *   L4 - Fundamental Laws: Quantization noise power = q²/12
 *        (derived from U[-q/2, q/2]), Bennett's conditions for
 *        white-noise approximation
 *   L5 - Algorithms: Lloyd-Max iterative algorithm for optimal
 *        non-uniform quantization, subtractive and non-subtractive
 *        dithering, noise shaping transfer function design
 *
 * References:
 *   - Bennett, "Spectra of Quantized Signals" BSTJ (1948)
 *   - Lloyd, "Least Squares Quantization in PCM" (1957/1982)
 *   - Max, "Quantizing for Minimum Distortion" IRE Trans (1960)
 *   - Gray & Neuhoff, "Quantization" IEEE Trans Info Theory (1998)
 *   - Oppenheim & Schafer (2010), Ch. 4
 *
 * Course Alignment:
 *   MIT 6.450 Digital Communications (Ch. 2, quantization)
 *   Stanford EE368 Digital Image Processing (Lloyd-Max)
 *   Berkeley EE123 DSP (quantization effects)
 */

#ifndef QUANTIZATION_H
#define QUANTIZATION_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1 — Quantizer Types
 * ========================================================================= */

/** Quantizer characteristic types. */
typedef enum {
    QUANTIZER_MID_TREAD,        /**< Zero at center of a quantization step (used in ADCs) */
    QUANTIZER_MID_RISE,          /**< Zero at boundary between steps */
    QUANTIZER_UNIFORM,           /**< Equal step sizes */
    QUANTIZER_LLOYD_MAX,         /**< Minimum MSE non-uniform quantizer */
    QUANTIZER_MU_LAW,            /**< μ-law companded quantizer (telephony) */
    QUANTIZER_A_LAW,             /**< A-law companded quantizer (European telephony) */
    QUANTIZER_TYPE_COUNT
} quantizer_type_t;

/* =========================================================================
 * L1 — Quantizer Specification
 * ========================================================================= */

/**
 * @brief Uniform scalar quantizer definition.
 *
 * Maps continuous input x to discrete output Q(x) = q * round((x - offset) / q)
 * clamped to [q * min_code, q * max_code].
 */
typedef struct {
    double          step_size;      /**< Quantization step q (V_LSB) */
    double          offset;         /**< DC offset in quantization levels */
    int64_t         min_code;       /**< Minimum output code */
    int64_t         max_code;       /**< Maximum output code */
    quantizer_type_t type;          /**< Quantizer characteristic type */
    double          overload_threshold; /**< Input level where overload begins */
} uniform_quantizer_t;

/**
 * @brief Non-uniform quantizer defined by decision and reconstruction levels.
 *
 * Decision levels d[0..N] partition input range: [d[i], d[i+1]) → code i.
 * Reconstruction levels r[0..N-1] give the output value for code i.
 *
 * For optimal (Lloyd-Max) quantizer:
 *   r[i] = E[x | x in [d[i], d[i+1])]  (centroid condition)
 *   d[i] = (r[i-1] + r[i]) / 2         (nearest-neighbor condition)
 */
typedef struct {
    uint32_t        num_levels;     /**< N = number of output levels */
    double         *decision_bounds;/**< Array of N+1 decision boundaries */
    double         *reconstruction; /**< Array of N reconstruction values */
    double          mse;            /**< Mean squared quantization error */
    quantizer_type_t type;          /**< Quantizer type */
    double          source_variance; /**< Variance of source distribution */
} nonuniform_quantizer_t;

/* =========================================================================
 * L1 — Quantization Error Statistics
 * ========================================================================= */

/**
 * @brief Quantization error statistics.
 *
 * Under Bennett's conditions (fine quantization, smooth PDF of input):
 * - Error e is uniformly distributed on [-q/2, q/2]
 * - Error is uncorrelated with input (and often white)
 * - Variance of error = q²/12
 */
typedef struct {
    double          mean_error;     /**< E[e] (should ≈ 0 for good quantizer) */
    double          var_error;      /**< Var[e] (≈ q²/12 for uniform) */
    double          rms_error;      /**< sqrt(var_error) */
    double          max_abs_error;  /**< ||e||_∞ */
    double          sqnr_db;        /**< Signal-to-quantization-noise ratio [dB] */
    uint64_t        overload_count; /**< Number of overload (clipping) events */
    uint64_t        total_samples;  /**< Total quantized samples */
    double          overload_rate;  /**< Fraction of samples that overloaded */
} quantization_error_stats_t;

/* =========================================================================
 * L1 — Uniform Quantization Functions
 * ========================================================================= */

/**
 * @brief Quantize a scalar input to nearest integer level.
 *
 * For mid-tread (standard ADC): code = round(input / step_size)
 * For mid-rise: code = floor(input / step_size + 0.5)
 *
 * @param x          Input value.
 * @param step_size  Quantization step (V_LSB).
 * @param type       Mid-tread or mid-rise.
 * @return Quantized code (integer).
 *
 * Complexity: O(1)
 */
int64_t quantize_uniform(double x, double step_size, quantizer_type_t type);

/**
 * @brief Reconstruct analog value from quantizer code.
 *
 * @param code       Quantized integer code.
 * @param step_size  Quantization step.
 * @param type       Quantizer type.
 * @return Reconstructed analog value.
 *
 * Complexity: O(1)
 */
double dequantize_uniform(int64_t code, double step_size, quantizer_type_t type);

/* =========================================================================
 * L4 — Quantization Noise Power
 * ========================================================================= */

/**
 * @brief Compute quantization noise power for uniform quantizer.
 *
 * For step size q, noise power = q² / 12.
 * This is derived from the uniform distribution of quantization error
 * over [-q/2, q/2]:
 *   ∫_{-q/2}^{q/2} e² * (1/q) de = [e³/(3q)]_{-q/2}^{q/2} = q²/12
 *
 * @param step_size Quantization step.
 * @return Noise power (variance) [units²].
 *
 * Complexity: O(1)
 * Reference: Bennett (1948), Widrow & Kollar "Quantization Noise" (2008)
 */
double quantization_noise_power(double step_size);

/**
 * @brief Compute SQNR (Signal-to-Quantization-Noise Ratio).
 *
 * For a full-scale sinusoid with amplitude A, sampled uniformly:
 *   SQNR = 10*log10( A²/2 / (q²/12) )
 *        = 6.02*b + 1.76 + 10*log10(A²/(V_fs/2)²)  [dB]
 *
 * @param signal_rms  RMS value of the input signal.
 * @param step_size   Quantization step.
 * @return SQNR in dB.
 *
 * Complexity: O(1)
 */
double quantization_sqnr(double signal_rms, double step_size);

/* =========================================================================
 * L5 — Lloyd-Max Optimal Quantizer Design
 * ========================================================================= */

/**
 * @brief Design a Lloyd-Max optimal quantizer for a given PDF.
 *
 * Minimizes MSE: E[(x - Q(x))²] for a given source distribution.
 * Iteratively alternates between centroid and nearest-neighbor conditions.
 *
 * The algorithm:
 *   1. Initialize decision boundaries
 *   2. Compute reconstruction levels as centroids: r_k = E[x | d_k ≤ x < d_{k+1}]
 *   3. Update decision boundaries: d_k = (r_{k-1} + r_k) / 2
 *   4. Repeat until convergence
 *
 * @param pdf_values    Array of PDF values over a grid.
 * @param grid_points   Corresponding grid points.
 * @param grid_len      Number of grid points.
 * @param num_levels    Desired output levels N.
 * @param[out] quant    Resulting quantizer (caller must free internal arrays).
 * @param max_iter      Maximum iterations.
 * @param tolerance     Convergence tolerance on MSE change.
 * @return 0 on success, -1 on error.
 *
 * Complexity: O(num_levels * grid_len * iterations)
 * Reference: Lloyd (1957), Max (1960)
 */
int lloyd_max_quantizer_design(const double *pdf_values,
                                const double *grid_points,
                                uint32_t grid_len,
                                uint32_t num_levels,
                                nonuniform_quantizer_t *quant,
                                uint32_t max_iter,
                                double tolerance);

/**
 * @brief Free memory allocated by lloyd_max_quantizer_design.
 */
void nonuniform_quantizer_free(nonuniform_quantizer_t *quant);

/* =========================================================================
 * L5 — Dithering
 * ========================================================================= */

/**
 * @brief Add subtractive (Schuchman) dither to a quantizer input.
 *
 * Subtractive dither: adds noise before quantization, subtracts same
 * noise from output to decorrelate error from signal without adding noise.
 *
 * @param input     Input sample.
 * @param dither    Dither noise sample (uniform [-q/2, q/2] recommended).
 * @param step_size Quantization step.
 * @param[out] dithered_input   Dithered value to quantize.
 * @param[out] subtractive_correction Correction to subtract after dequantization.
 *
 * Complexity: O(1)
 * Reference: Schuchman (1964), Gray & Stockham (1993)
 */
void dither_subtractive(double input, double dither, double step_size,
                         double *dithered_input, double *subtractive_correction);

/**
 * @brief Add non-subtractive dither (TPDF — triangular PDF dither).
 *
 * Two independent uniform random values summed to create triangular PDF.
 * This renders quantization error moments conditionally independent of input.
 *
 * @param input     Input sample.
 * @param dither1   First uniform random [-q, q] or [-q/2, q/2].
 * @param dither2   Second uniform random (same range).
 * @param step_size Quantization step.
 * @return Dithered input (quantize this).
 *
 * Complexity: O(1)
 */
double dither_nonsubtractive_tpdf(double input, double dither1,
                                   double dither2, double step_size);

/* =========================================================================
 * L5 — Noise Shaping (Sigma-Delta Foundation)
 * ========================================================================= */

/**
 * @brief First-order noise shaping.
 *
 * The quantizer error e[n] is filtered by NTF(z) = 1 - z^{-1}
 * (a high-pass filter), pushing noise to higher frequencies.
 *
 * @param input_sample    Current input sample.
 * @param step_size       Quantization step.
 * @param[in,out] error_prev Previous quantization error (state variable).
 * @return Quantized output code.
 *
 * Complexity: O(1)
 * Reference: Norsworthy, Schreier, Temes "Delta-Sigma Data Converters" (1997)
 */
int64_t noise_shape_first_order(double input_sample, double step_size,
                                 double *error_prev);

/**
 * @brief Second-order noise shaping with NTF(z) = (1 - z^{-1})².
 *
 * @param input_sample    Current input.
 * @param step_size       Quantization step.
 * @param[in,out] e1     First state (previous error).
 * @param[in,out] e2     Second state (error from two samples ago).
 * @return Quantized output.
 *
 * Complexity: O(1)
 */
int64_t noise_shape_second_order(double input_sample, double step_size,
                                  double *e1, double *e2);

/* =========================================================================
 * L5 — Quantization Error Statistics Collection
 * ========================================================================= */

/**
 * @brief Compute quantization error statistics for a batch of samples.
 *
 * Feed in original and quantized values to get aggregate statistics.
 *
 * @param original  Array of original analog values.
 * @param quantized Array of dequantized (reconstructed) values.
 * @param n         Number of samples.
 * @param step_size Quantization step.
 * @param[out] stats Output statistics (caller allocates).
 *
 * Complexity: O(n)
 */
void quantization_error_stats(const double *original,
                               const double *quantized,
                               uint64_t n,
                               double step_size,
                               quantization_error_stats_t *stats);

/* =========================================================================
 * L5 — μ-Law and A-Law Companding (Telephony)
 * ========================================================================= */

/**
 * @brief μ-law compression (North American telephony, μ=255).
 *
 *    y = sign(x) * ln(1 + μ*|x|) / ln(1 + μ)   for |x| ≤ 1
 *
 * @param x  Input in [-1, 1].
 * @param mu Compression parameter (standard: 255).
 * @return Compressed value in [-1, 1].
 *
 * Complexity: O(1)
 * Reference: ITU-T G.711
 */
double mu_law_compress(double x, double mu);

/**
 * @brief μ-law expansion (decode).
 */
double mu_law_expand(double y, double mu);

/**
 * @brief A-law compression (European telephony, A=87.6).
 *
 *    y = sign(x) * { A*|x| / (1+ln(A))           if |x| ≤ 1/A
 *                   { (1 + ln(A*|x|)) / (1+ln(A)) if 1/A < |x| ≤ 1
 *
 * @param x Input in [-1, 1].
 * @param A Compression parameter (standard: 87.6).
 * @return Compressed value in [-1, 1].
 *
 * Complexity: O(1)
 * Reference: ITU-T G.711
 */
double a_law_compress(double x, double A);

/**
 * @brief A-law expansion (decode).
 */
double a_law_expand(double y, double A);

/* =========================================================================
 * L6 — Quantizer Performance Evaluation
 * ========================================================================= */

/**
 * @brief Evaluate quantizer SNR vs input amplitude sweep.
 *
 * Sweeps input sine amplitude from -60 dBFS to 0 dBFS and computes
 * SQNR at each point. Classic characterization plot.
 *
 * @param quant        Quantizer specification.
 * @param amplitudes_db Array of input amplitudes in dBFS.
 * @param num_amps     Number of amplitude points.
 * @param f_points_per_amp Points per amplitude for averaging.
 * @param[out] sqnr_db Output SQNR array [num_amps].
 *
 * Complexity: O(num_amps * points_per_amp)
 */
void quantizer_snr_sweep(const uniform_quantizer_t *quant,
                          const double *amplitudes_db,
                          uint32_t num_amps,
                          uint32_t f_points_per_amp,
                          double *sqnr_db);

/* =========================================================================
 * L2 — Overload / Granular Distortion
 * ========================================================================= */

/**
 * @brief Compute overload distortion power.
 *
 * Overload occurs when |x| > V_fs/2. For Gaussian source with
 * variance σ², the overload probability is 2*Q(V_fs/(2σ)).
 *
 * @param v_fs        Full-scale voltage range.
 * @param source_rms  RMS of the input signal.
 * @return Overload distortion power (variance contribution).
 */
double overload_distortion_power(double v_fs, double source_rms);

/**
 * @brief Compute optimal step size that minimizes total distortion
 *        (granular + overload) for a given source distribution.
 *
 * @param source_rms   RMS of the input signal.
 * @param n_bits       Number of quantization bits.
 * @return Optimal step size.
 */
double optimal_step_size(double source_rms, uint32_t n_bits);

#ifdef __cplusplus
}
#endif

#endif /* QUANTIZATION_H */
