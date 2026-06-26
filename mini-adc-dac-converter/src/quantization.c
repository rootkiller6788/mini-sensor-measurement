/**
 * @file quantization.c
 * @brief Quantization theory implementation: uniform/non-uniform,
 *        Lloyd-Max, dithering, noise shaping, mu-law/A-law.
 *
 * Each function represents an independent knowledge point from
 * quantization theory (Bennett 1948, Lloyd 1957, Max 1960).
 */

#include "quantization.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =========================================================================
 * L1 — Uniform Quantization
 * ========================================================================= */

int64_t quantize_uniform(double x, double step_size, quantizer_type_t type)
{
    if (step_size <= 0.0) return 0;
    if (type == QUANTIZER_MID_TREAD) {
        /* Mid-tread: zero is at the center of a quantization step.
         * code = round(x / q)
         * This is the standard ADC quantizer characteristic. */
        return (int64_t)llround(x / step_size);
    } else {
        /* Mid-rise: zero is at the transition between two steps.
         * code = floor(x / q + 0.5) — equivalent to round but shifted. */
        return (int64_t)floor(x / step_size + 0.5);
    }
}

double dequantize_uniform(int64_t code, double step_size, quantizer_type_t type)
{
    (void)type; /* Type affects only quantization, not reconstruction */
    if (step_size <= 0.0) return NAN;
    /* For mid-tread: x_reconstructed = code * q */
    return (double)code * step_size;
}

/* =========================================================================
 * L4 — Quantization Noise Power (Fundamental Result)
 * ========================================================================= */

double quantization_noise_power(double step_size)
{
    if (step_size <= 0.0) return NAN;
    /* Derivation: error e ~ Uniform[-q/2, q/2]
     *   E[e²] = ∫_{-q/2}^{q/2} e² * (1/q) de
     *         = [e³/(3q)]_{-q/2}^{q/2}
     *         = (q³/24 + q³/24) / q
     *         = q²/12
     *
     * This result requires Bennett's conditions:
     * 1. Fine quantization (q small relative to signal range)
     * 2. Smooth input PDF
     * 3. No overload
     *
     * Under these conditions, the quantization error is well-approximated
     * as additive white noise uncorrelated with the input.
     */
    return (step_size * step_size) / 12.0;
}

double quantization_sqnr(double signal_rms, double step_size)
{
    if (signal_rms <= 0.0 || step_size <= 0.0) return NAN;
    double noise_power = quantization_noise_power(step_size);
    if (noise_power <= 0.0) return NAN;
    double signal_power = signal_rms * signal_rms;
    return 10.0 * log10(signal_power / noise_power);
}

/* =========================================================================
 * L5 — Lloyd-Max Optimal Quantizer Design
 * ========================================================================= */

/**
 * Compute centroid of a PDF interval [a, b].
 *
 * For discrete grid representation:
 *   centroid = Σ(x_i * pdf(x_i)) / Σ(pdf(x_i))  for x_i in [a, b]
 */
static double compute_centroid(const double *pdf_values,
                                const double *grid_points,
                                uint32_t start_idx, uint32_t end_idx)
{
    double numerator = 0.0;
    double denominator = 0.0;
    for (uint32_t i = start_idx; i <= end_idx && i < (uint32_t)-1; i++) {
        /* Use trapezoidal rule weighting */
        double weight = pdf_values[i];
        numerator += grid_points[i] * weight;
        denominator += weight;
    }
    if (denominator <= 0.0) {
        /* Empty interval: return midpoint */
        return (grid_points[start_idx] + grid_points[end_idx]) * 0.5;
    }
    return numerator / denominator;
}

int lloyd_max_quantizer_design(const double *pdf_values,
                                const double *grid_points,
                                uint32_t grid_len,
                                uint32_t num_levels,
                                nonuniform_quantizer_t *quant,
                                uint32_t max_iter,
                                double tolerance)
{
    if (!pdf_values || !grid_points || !quant || grid_len < 2 ||
        num_levels < 2 || num_levels > grid_len || max_iter == 0) {
        return -1;
    }

    /* Allocate decision boundaries (num_levels+1) and reconstruction (num_levels) */
    quant->num_levels = num_levels;
    quant->decision_bounds = (double *)malloc((num_levels + 1) * sizeof(double));
    quant->reconstruction = (double *)malloc(num_levels * sizeof(double));
    if (!quant->decision_bounds || !quant->reconstruction) {
        free(quant->decision_bounds);
        free(quant->reconstruction);
        quant->decision_bounds = NULL;
        quant->reconstruction = NULL;
        return -1;
    }

    /* Initialize decision boundaries uniformly over the grid range. */
    double x_min = grid_points[0];
    double x_max = grid_points[grid_len - 1];
    double range = x_max - x_min;

    /* Initial boundaries: uniform partitions */
    for (uint32_t k = 0; k <= num_levels; k++) {
        quant->decision_bounds[k] = x_min + range * (double)k / (double)num_levels;
    }
    /* Initial reconstruction levels: midpoints of decision intervals */
    for (uint32_t k = 0; k < num_levels; k++) {
        quant->reconstruction[k] = (quant->decision_bounds[k] +
                                     quant->decision_bounds[k + 1]) * 0.5;
    }

    double prev_mse = INFINITY;

    for (uint32_t iter = 0; iter < max_iter; iter++) {
        /* Step 1: Update decision boundaries as midpoints of adjacent
         * reconstruction levels (Nearest-Neighbor Condition). */
        for (uint32_t k = 1; k < num_levels; k++) {
            quant->decision_bounds[k] = (quant->reconstruction[k - 1] +
                                          quant->reconstruction[k]) * 0.5;
        }
        /* Outer boundaries stay fixed at min/max grid values */
        quant->decision_bounds[0] = x_min;
        quant->decision_bounds[num_levels] = x_max;

        /* Step 2: Update reconstruction levels as centroids
         * (Centroid Condition). */
        for (uint32_t k = 0; k < num_levels; k++) {
            /* Find grid indices that fall into interval [d[k], d[k+1]) */
            uint32_t start_idx = 0, end_idx = grid_len - 1;
            /* Search for start */
            while (start_idx < grid_len &&
                   grid_points[start_idx] < quant->decision_bounds[k]) {
                start_idx++;
            }
            /* Search for end */
            end_idx = start_idx;
            while (end_idx < grid_len &&
                   grid_points[end_idx] < quant->decision_bounds[k + 1]) {
                end_idx++;
            }
            if (end_idx > 0) end_idx--;

            if (start_idx <= end_idx && end_idx < grid_len) {
                quant->reconstruction[k] = compute_centroid(pdf_values,
                    grid_points, start_idx, end_idx);
            } else {
                /* Empty interval: keep midpoint */
                quant->reconstruction[k] = (quant->decision_bounds[k] +
                    quant->decision_bounds[k + 1]) * 0.5;
            }
        }

        /* Step 3: Compute MSE */
        double mse = 0.0;
        double total_weight = 0.0;
        for (uint32_t k = 0; k < num_levels; k++) {
            uint32_t start_idx = 0;
            while (start_idx < grid_len &&
                   grid_points[start_idx] < quant->decision_bounds[k])
                start_idx++;
            uint32_t end_idx = start_idx;
            while (end_idx < grid_len &&
                   grid_points[end_idx] < quant->decision_bounds[k + 1])
                end_idx++;
            if (end_idx > 0) end_idx--;

            for (uint32_t i = start_idx; i <= end_idx && i < grid_len; i++) {
                double diff = grid_points[i] - quant->reconstruction[k];
                mse += pdf_values[i] * diff * diff;
                total_weight += pdf_values[i];
            }
        }
        if (total_weight > 0.0) mse /= total_weight;

        /* Step 4: Check convergence */
        if (fabs(prev_mse - mse) < tolerance) {
            quant->mse = mse;
            break;
        }
        prev_mse = mse;
        quant->mse = mse;
    }

    quant->type = QUANTIZER_LLOYD_MAX;
    return 0;
}

void nonuniform_quantizer_free(nonuniform_quantizer_t *quant)
{
    if (!quant) return;
    free(quant->decision_bounds);
    free(quant->reconstruction);
    quant->decision_bounds = NULL;
    quant->reconstruction = NULL;
    quant->num_levels = 0;
}

/* =========================================================================
 * L5 — Dithering
 * ========================================================================= */

void dither_subtractive(double input, double dither, double step_size,
                         double *dithered_input, double *subtractive_correction)
{
    (void)step_size; /* Step size is inherent in quantizer, not used here */
    /* Add dither before quantization: x_d = input + dither */
    *dithered_input = input + dither;

    /* Store correction to subtract after dequantization:
     *   output = Q(input + dither) - dither
     * This makes the quantization error independent of the input
     * while not increasing total noise (dither is added then subtracted). */
    *subtractive_correction = dither;
}

double dither_nonsubtractive_tpdf(double input, double dither1,
                                   double dither2, double step_size)
{
    /* Triangular PDF dither: sum of two independent uniform random variables.
     *
     * TPDF dither of width 2q (each uniform over [-q/2, q/2]) renders
     * the first and second moments of the quantization error conditionally
     * independent of the input signal.
     *
     * This eliminates distortion (spurious tones) but adds more noise
     * than subtractive dither. */
    (void)step_size; /* Step size is inherent in the dither range */
    return input + dither1 + dither2;
}

/* =========================================================================
 * L5 — Noise Shaping (ΔΣ foundation)
 * ========================================================================= */

int64_t noise_shape_first_order(double input_sample, double step_size,
                                 double *error_prev)
{
    /* First-order noise shaping with NTF(z) = 1 - z^{-1}.
     *
     * Input to quantizer: v[n] = input[n] - e[n-1]
     * (The previous quantization error is subtracted from the input,
     * forming a negative feedback loop around the quantizer.)
     *
     * Output: y[n] = Q(v[n]) = v[n] + e[n]
     * So:     y[n] = input[n] + (e[n] - e[n-1])
     *
     * The input passes through with unity gain (STF = 1),
     * while the error is high-pass filtered (NTF = 1 - z^{-1}).
     *
     * This pushes quantization noise to higher frequencies,
     * reducing in-band noise when followed by a low-pass filter.
     */
    if (!error_prev) return 0;
    if (step_size <= 0.0) return 0;

    /* Shape: add negative of previous error to input */
    double shaped_input = input_sample - *error_prev;

    /* Quantize */
    int64_t code = (int64_t)llround(shaped_input / step_size);

    /* Compute new error */
    double reconstructed = (double)code * step_size;
    double new_error = reconstructed - shaped_input;

    /* Store for next sample */
    *error_prev = new_error;

    return code;
}

int64_t noise_shape_second_order(double input_sample, double step_size,
                                  double *e1, double *e2)
{
    /* Second-order noise shaping with NTF(z) = (1 - z^{-1})².
     *
     * Loop filter: H(z) = -2z^{-1} + z^{-2}
     *
     * v[n] = input[n] - 2*e[n-1] + e[n-2]
     * y[n] = v[n] + e[n]
     *      = input[n] + e[n] - 2*e[n-1] + e[n-2]
     *
     * NTF(z) = (1 - z^{-1})² = 1 - 2z^{-1} + z^{-2}
     *
     * This gives 15 dB/octave improvement in in-band noise
     * (vs 9 dB/octave for first-order).
     */
    if (!e1 || !e2) return 0;
    if (step_size <= 0.0) return 0;

    double shaped_input = input_sample - 2.0 * (*e1) + (*e2);

    int64_t code = (int64_t)llround(shaped_input / step_size);

    double reconstructed = (double)code * step_size;
    double new_error = reconstructed - shaped_input;

    /* Shift error history */
    *e2 = *e1;
    *e1 = new_error;

    return code;
}

/* =========================================================================
 * L5 — Quantization Error Statistics
 * ========================================================================= */

void quantization_error_stats(const double *original,
                               const double *quantized,
                               uint64_t n,
                               double step_size,
                               quantization_error_stats_t *stats)
{
    if (!original || !quantized || !stats || n == 0) return;
    memset(stats, 0, sizeof(*stats));
    stats->total_samples = n;

    double sum_error = 0.0;
    double sum_error_sq = 0.0;
    double max_abs = 0.0;
    double signal_power = 0.0;
    uint64_t overload_count = 0;

    for (uint64_t i = 0; i < n; i++) {
        double error = quantized[i] - original[i];
        sum_error += error;
        sum_error_sq += error * error;
        if (fabs(error) > max_abs) max_abs = fabs(error);
        signal_power += original[i] * original[i];

        /* Overload: error magnitude exceeds q/2 beyond what quantization
         * alone would produce (i.e., signal exceeded quantizer range). */
        if (fabs(error) > step_size * 1.5) {
            overload_count++;
        }
    }

    stats->mean_error = sum_error / (double)n;
    stats->var_error = sum_error_sq / (double)n - stats->mean_error * stats->mean_error;
    if (stats->var_error < 0.0) stats->var_error = 0.0;
    stats->rms_error = sqrt(stats->var_error);
    stats->max_abs_error = max_abs;
    stats->overload_count = overload_count;
    stats->overload_rate = (double)overload_count / (double)n;

    double signal_rms = sqrt(signal_power / (double)n);
    double noise_rms = stats->rms_error;
    if (noise_rms > 0.0 && signal_rms > 0.0) {
        stats->sqnr_db = 20.0 * log10(signal_rms / noise_rms);
    } else {
        stats->sqnr_db = INFINITY;
    }
}

/* =========================================================================
 * L5 — μ-Law Companding (ITU-T G.711)
 * ========================================================================= */

double mu_law_compress(double x, double mu)
{
    if (mu <= 0.0) return NAN;
    double sign_val = (x >= 0.0) ? 1.0 : -1.0;
    double abs_x = fabs(x);
    if (abs_x > 1.0) abs_x = 1.0;
    return sign_val * log1p(mu * abs_x) / log1p(mu);
}

double mu_law_expand(double y, double mu)
{
    if (mu <= 0.0) return NAN;
    double sign_val = (y >= 0.0) ? 1.0 : -1.0;
    double abs_y = fabs(y);
    if (abs_y > 1.0) abs_y = 1.0;
    return sign_val * (expm1(abs_y * log1p(mu))) / mu;
}

/* =========================================================================
 * L5 — A-Law Companding (ITU-T G.711)
 * ========================================================================= */

double a_law_compress(double x, double A)
{
    if (A <= 1.0) return NAN;
    double sign_val = (x >= 0.0) ? 1.0 : -1.0;
    double abs_x = fabs(x);
    if (abs_x > 1.0) abs_x = 1.0;

    double result;
    if (abs_x <= 1.0 / A) {
        /* Linear region: small signals amplified by A/(1+ln(A)) */
        result = A * abs_x / (1.0 + log(A));
    } else {
        /* Logarithmic region */
        result = (1.0 + log(A * abs_x)) / (1.0 + log(A));
    }
    return sign_val * result;
}

double a_law_expand(double y, double A)
{
    if (A <= 1.0) return NAN;
    double sign_val = (y >= 0.0) ? 1.0 : -1.0;
    double abs_y = fabs(y);
    if (abs_y > 1.0) abs_y = 1.0;

    double result;
    double threshold = 1.0 / (1.0 + log(A));
    if (abs_y <= threshold) {
        result = abs_y * (1.0 + log(A)) / A;
    } else {
        result = exp(abs_y * (1.0 + log(A)) - 1.0) / A;
    }
    if (result > 1.0) result = 1.0;
    return sign_val * result;
}

/* =========================================================================
 * L6 — Quantizer SNR Sweep
 * ========================================================================= */

void quantizer_snr_sweep(const uniform_quantizer_t *quant,
                          const double *amplitudes_db,
                          uint32_t num_amps,
                          uint32_t f_points_per_amp,
                          double *sqnr_db)
{
    if (!quant || !amplitudes_db || !sqnr_db || num_amps == 0 ||
        quant->step_size <= 0.0) {
        return;
    }

    for (uint32_t a = 0; a < num_amps; a++) {
        /* Convert dBFS to linear amplitude.
         * 0 dBFS = full-scale: amplitude = V_ref/2 (peak).
         * The quantizer range is [min_code*step_size, max_code*step_size]. */
        double fs_peak = (double)quant->max_code * quant->step_size * 0.5;
        double linear_amplitude = fs_peak * pow(10.0, amplitudes_db[a] / 20.0);
        if (linear_amplitude <= 0.0) {
            sqnr_db[a] = -INFINITY;
            continue;
        }

        double error_power = 0.0;
        double signal_power = 0.0;
        for (uint32_t k = 0; k < f_points_per_amp; k++) {
            double phase = 2.0 * M_PI * (double)k / (double)f_points_per_amp;
            double input = linear_amplitude * sin(phase);
            int64_t code = quantize_uniform(input, quant->step_size, quant->type);
            double reconstructed = dequantize_uniform(code, quant->step_size, quant->type);
            double error = reconstructed - input;
            error_power += error * error;
            signal_power += input * input;
        }
        error_power /= (double)f_points_per_amp;
        signal_power /= (double)f_points_per_amp;
        if (error_power > 0.0 && signal_power > 0.0) {
            sqnr_db[a] = 10.0 * log10(signal_power / error_power);
        } else {
            sqnr_db[a] = INFINITY;
        }
    }
}

/* =========================================================================
 * L2 — Overload and Granular Distortion
 * ========================================================================= */

/**
 * Q-function approximation using polynomial fit.
 *
 * Q(x) = (1/sqrt(2π)) * ∫_x^∞ exp(-t²/2) dt
 * ≈ 0.5 * erfc(x / sqrt(2))
 */
static double Q_function(double x)
{
    return 0.5 * erfc(x / sqrt(2.0));
}

double overload_distortion_power(double v_fs, double source_rms)
{
    if (v_fs <= 0.0 || source_rms <= 0.0) return NAN;
    /* For Gaussian source:
     * Overload probability: P_ol = 2 * Q(V_fs / (2 * σ))
     * Overload distortion: E[(x - V_fs/2)² | |x| > V_fs/2] * P_ol
     *
     * Approximate using tail integration. Simplified:
     * overload_power ≈ 2 * σ² * Q(α) * (1 + α²) where α = V_fs/(2σ)
     */
    double alpha = v_fs / (2.0 * source_rms);
    double q_alpha = Q_function(alpha);
    if (q_alpha <= 0.0) return 0.0;
    return 2.0 * source_rms * source_rms * q_alpha * (1.0 + alpha * alpha);
}

double optimal_step_size(double source_rms, uint32_t n_bits)
{
    if (source_rms <= 0.0 || n_bits == 0 || n_bits > 32) return NAN;
    /* For Gaussian source with variance σ², the optimal step size
     * balances granular and overload distortion.
     *
     * Rule of thumb: optimal loading factor γ = V_fs / (2σ) ≈ 4.
     * (The ADC clips at ~4σ, losing ~0.006% of samples to overload.)
     *
     * V_fs = 2^N * q
     * q_opt = 2σ * γ_opt / 2^N = σ * γ_opt / 2^{N-1}
     */
    double gamma_opt = 4.0; /* Loading factor for Gaussian signals */
    double v_fs = 2.0 * source_rms * gamma_opt;
    double n_levels = pow(2.0, (double)n_bits);
    return v_fs / n_levels;
}
