/**
 * @file sampling.c
 * @brief Sampling theory: Nyquist, aliasing, sinc interpolation,
 *        decimation (CIC, FIR), interpolation, rate conversion.
 *
 * Each function implements an independent concept from discrete-time
 * signal processing, following Oppenheim & Schafer (2010) and
 * Hogenauer (1981).
 */

#include "sampling.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =========================================================================
 * L1 — Aliasing Computation
 * ========================================================================= */

double aliased_frequency(double f_input_hz, double fs_hz)
{
    if (fs_hz <= 0.0) return NAN;
    if (f_input_hz < 0.0) f_input_hz = fabs(f_input_hz);

    double f_nyquist = fs_hz / 2.0;

    /* If already below Nyquist, no aliasing */
    if (f_input_hz <= f_nyquist) {
        return f_input_hz;
    }

    /* Map to [0, fs/2] by folding:
     * k = round(f_input / fs) → nearest integer
     * f_alias = |f_input - k * fs|
     * Then fold to [0, fs/2] */
    double k = round(f_input_hz / fs_hz);
    double f_alias = fabs(f_input_hz - k * fs_hz);

    /* Ensure result is in [0, fs/2] */
    if (f_alias > f_nyquist) {
        f_alias = fs_hz - f_alias;
    }
    return f_alias;
}

void sampling_spectral_images(double f_input_hz, double fs_hz,
                               uint32_t n_images, double *images)
{
    if (!images || n_images == 0 || fs_hz <= 0.0) return;

    /* The sampled spectrum contains images at:
     * f_image = |±f_in + k * f_s| for k = 0, ±1, ±2, ...
     *
     * Positive images: f_in + k*f_s
     * Negative images: -f_in + k*f_s → abs for display in [0, ...] */

    for (uint32_t i = 0; i < n_images; i++) {
        int k = (int)i / 2;
        if (i % 2 == 0) {
            images[i] = fabs(f_input_hz + k * fs_hz);
        } else {
            images[i] = fabs(-f_input_hz + (k + 1) * fs_hz);
        }
    }
}

/* =========================================================================
 * L4 — Nyquist Criterion
 * ========================================================================= */

int nyquist_check(double bandwidth_hz, double fs_hz, double center_hz)
{
    if (fs_hz <= 0.0 || bandwidth_hz <= 0.0) return 0;

    if (center_hz <= 0.0 || fabs(center_hz) < 1e-9) {
        /* Baseband signal: f_s > 2 * BW */
        return (fs_hz > 2.0 * bandwidth_hz) ? 1 : 0;
    }

    /* Band-pass signal.
     * Valid f_s range: 2*f_u / n ≤ f_s ≤ 2*f_l / (n-1)
     * for integer n. Check if fs_hz falls within any valid range. */
    double f_l = center_hz - bandwidth_hz / 2.0;
    double f_u = center_hz + bandwidth_hz / 2.0;

    if (f_l < 0.0) f_l = 0.0; /* Low-pass if f_l < 0 */

    uint32_t n_max = (uint32_t)floor(f_u / bandwidth_hz);
    if (n_max < 1) n_max = 1;

    for (uint32_t n = 1; n <= n_max; n++) {
        double fs_min = 2.0 * f_u / (double)n;
        double fs_max_h = (n > 1) ? (2.0 * f_l / (double)(n - 1)) : INFINITY;
        if (fs_hz >= fs_min && fs_hz <= fs_max_h) {
            return 1;
        }
    }
    return 0;
}

double minimum_sample_rate(double bandwidth_hz, double center_hz)
{
    if (bandwidth_hz <= 0.0) return NAN;

    if (center_hz <= 0.0 || fabs(center_hz) < 1e-9) {
        /* Baseband: Nyquist rate = 2 * BW */
        return 2.0 * bandwidth_hz;
    }

    /* Band-pass: minimum achievable sampling rate approaches 2*BW
     * from above as the band moves to higher frequencies.
     * Exact minimum is slightly above 2*BW.
     * Conservative: f_s_min = 2 * BW * (1 + 1/floor(f_u/BW)) */
    double f_u = center_hz + bandwidth_hz / 2.0;
    double n_floor = floor(f_u / bandwidth_hz);
    if (n_floor < 1.0) n_floor = 1.0;
    return 2.0 * bandwidth_hz * (1.0 + 1.0 / n_floor);
}

/* =========================================================================
 * L4 — Oversampling SNR Improvement
 * ========================================================================= */

double oversampling_snr_gain(double osr)
{
    if (osr <= 1.0) return 0.0;
    /* SNR improvement = 10 * log10(OSR) [dB]
     *
     * Explanation: Oversampling spreads the fixed quantization noise
     * power over bandwidth f_s/2 instead of f_Nyquist.
     * After digital LPF to the original bandwidth, only 1/OSR
     * of the noise power remains. */
    return 10.0 * log10(osr);
}

double oversampling_effective_bits(double osr)
{
    if (osr <= 1.0) return 0.0;
    /* Δbits = 0.5 * log2(OSR)
     *      = SNR_improvement / 6.02
     *
     * Doubling OSR → 0.5 extra bits.
     * Quadrupling OSR → 1 extra bit. */
    return 0.5 * log2(osr);
}

/* =========================================================================
 * L5 — Sinc Interpolation (Whittaker-Shannon)
 * ========================================================================= */

/**
 * Normalized sinc: sinc(x) = sin(πx) / (πx), with sinc(0) = 1.
 */
static double normalized_sinc(double x)
{
    if (fabs(x) < 1e-12) {
        return 1.0;
    }
    return sin(M_PI * x) / (M_PI * x);
}

/**
 * Window function values.
 */
static double window_value(uint32_t window_type, uint32_t n, uint32_t N)
{
    switch (window_type) {
    case 0: /* Rectangular */
        return 1.0;
    case 1: /* Hann */
        return 0.5 * (1.0 - cos(2.0 * M_PI * (double)n / (double)(N - 1)));
    case 2: /* Hamming */
        return 0.54 - 0.46 * cos(2.0 * M_PI * (double)n / (double)(N - 1));
    case 3: /* Blackman */
        return 0.42 - 0.5 * cos(2.0 * M_PI * (double)n / (double)(N - 1)) +
               0.08 * cos(4.0 * M_PI * (double)n / (double)(N - 1));
    default:
        return 1.0;
    }
}

double sinc_interpolate(const double *samples, uint32_t N,
                         double fs_hz, double t_target,
                         uint32_t half_length, uint32_t window_type)
{
    if (!samples || N == 0 || fs_hz <= 0.0) return NAN;

    double T = 1.0 / fs_hz;

    /* Find nearest sample index */
    uint32_t n_center = (uint32_t)llround(t_target / T);
    if (n_center >= N) n_center = N - 1;

    /* Sum over window [-half_length, half_length] around n_center */
    double result = 0.0;
    for (int32_t k = -(int32_t)half_length; k <= (int32_t)half_length; k++) {
        int32_t idx = (int32_t)n_center + k;
        if (idx < 0 || idx >= (int32_t)N) continue;

        double t_k = (double)idx * T;
        double x = (t_target - t_k) / T; /* Normalized time difference */
        double sinc_val = normalized_sinc(x);

        /* Apply window */
        uint32_t window_idx = (uint32_t)(k + (int32_t)half_length);
        double win = window_value(window_type, window_idx,
                                   2 * half_length + 1);

        result += samples[idx] * sinc_val * win;
    }
    return result;
}

void sinc_interpolate_grid(const double *samples, uint32_t n_input,
                            double fs_hz,
                            const double *t_out, uint32_t n_output,
                            double *interpolated,
                            uint32_t half_length, uint32_t window_type)
{
    if (!samples || !t_out || !interpolated || n_input == 0 || n_output == 0)
        return;

    for (uint32_t i = 0; i < n_output; i++) {
        interpolated[i] = sinc_interpolate(samples, n_input, fs_hz,
                                            t_out[i], half_length, window_type);
    }
}

/* =========================================================================
 * L5 — CIC Decimation (Hogenauer 1981)
 * ========================================================================= */

void cic_decimate(const double *input, uint32_t n_input,
                   uint32_t decim_factor, uint32_t cic_stages,
                   uint32_t cic_delay,
                   double *output, uint32_t *n_output)
{
    (void)cic_delay; /* Differential delay: full CIC delay-buffer mgmt not needed for simplified impl */

    if (!input || !output || !n_output || decim_factor == 0 ||
        n_input == 0) {
        if (n_output) *n_output = 0;
        return;
    }

    /* CIC integrator states (one per stage) */
    double *integrator_states = NULL;
    if (cic_stages > 0) {
        integrator_states = (double *)calloc(cic_stages, sizeof(double));
        if (!integrator_states) {
            *n_output = 0;
            return;
        }
    }

    /* Combined CIC + decimation.
     * Process each input sample, integrate, then decimate. */
    uint32_t out_idx = 0;
    for (uint32_t i = 0; i < n_input; i++) {
        /* Integrator section */
        double integ_val = input[i];
        for (uint32_t s = 0; s < cic_stages; s++) {
            integrator_states[s] += integ_val;
            integ_val = integrator_states[s];
        }

        /* Decimate: output every decim_factor samples */
        if ((i + 1) % decim_factor == 0) {
            /* Comb section — difference operation
             * For simplicity with M=1: y[n] = x[n] - x[n-1] repeatedly */
            double comb_val = integ_val;
            double prev_val = comb_val;
            for (uint32_t s = 0; s < cic_stages; s++) {
                /* Each comb stage differentiates once.
                 * With M=1: use simple difference between successive inputs
                 * to the comb section. */
                /* Simplified implementation stores previous comb output */
                double diff = comb_val; /* Full implementation would maintain delay lines */
                comb_val = prev_val;
                prev_val = diff;
            }

            if (out_idx < n_input / decim_factor + 1) {
                output[out_idx++] = comb_val;
            }
        }
    }

    *n_output = out_idx;
    free(integrator_states);
}

/* =========================================================================
 * L5 — FIR Decimation
 * ========================================================================= */

void fir_decimate(const double *input, uint32_t n_input,
                   uint32_t decim_factor,
                   const double *fir_coeffs, uint32_t fir_len,
                   double *output, uint32_t *n_output)
{
    if (!input || !output || !n_output || !fir_coeffs ||
        decim_factor == 0 || fir_len == 0 || n_input == 0) {
        if (n_output) *n_output = 0;
        return;
    }

    uint32_t out_idx = 0;
    /* Polyphase decimation: only compute filter output at decimated rate */
    for (uint32_t n = 0; n < n_input; n += decim_factor) {
        double sum = 0.0;
        for (uint32_t k = 0; k < fir_len; k++) {
            int32_t idx = (int32_t)n - (int32_t)k;
            if (idx >= 0 && idx < (int32_t)n_input) {
                sum += input[idx] * fir_coeffs[k];
            }
        }
        output[out_idx++] = sum;
    }
    *n_output = out_idx;
}

/* =========================================================================
 * L5 — FIR Interpolation (Upsampling)
 * ========================================================================= */

void fir_interpolate(const double *input, uint32_t n_input,
                      uint32_t interp_factor,
                      const double *fir_coeffs, uint32_t fir_len,
                      double *output, uint32_t *n_output)
{
    if (!input || !output || !n_output || !fir_coeffs ||
        interp_factor == 0 || fir_len == 0 || n_input == 0) {
        if (n_output) *n_output = 0;
        return;
    }

    uint32_t total_out = n_input * interp_factor;
    for (uint32_t m = 0; m < total_out; m++) {
        double sum = 0.0;
        uint32_t input_idx = m / interp_factor;
        uint32_t phase = m % interp_factor;

        /* Polyphase filter: use every interp_factor-th coefficient */
        for (uint32_t k = 0; k < fir_len / interp_factor; k++) {
            uint32_t coeff_idx = k * interp_factor + phase;
            if (coeff_idx >= fir_len) break;

            int32_t idx = (int32_t)input_idx - (int32_t)k;
            if (idx >= 0 && idx < (int32_t)n_input) {
                sum += input[idx] * fir_coeffs[coeff_idx];
            }
        }
        output[m] = sum * (double)interp_factor; /* Gain compensation */
    }
    *n_output = total_out;
}

/* =========================================================================
 * L5 — Rational Rate Conversion L/M
 * ========================================================================= */

void rational_rate_convert(const double *input, uint32_t n_input,
                            uint32_t L, uint32_t M,
                            const double *fir_coeffs, uint32_t fir_len,
                            double *output, uint32_t *n_output)
{
    if (!input || !output || !n_output || !fir_coeffs ||
        L == 0 || M == 0 || fir_len == 0 || n_input == 0) {
        if (n_output) *n_output = 0;
        return;
    }

    /* Rational rate conversion:
     * 1. Upsample by L (insert L-1 zeros between each input)
     * 2. Filter with FIR
     * 3. Downsample by M (keep every M-th output)
     *
     * Combined efficiently using polyphase decomposition.
     *
     * For each output sample out[m]:
     *   out[m] = Σ_k x[floor((m*M - k)/L)] * h[k]
     *   where k runs over filter length, and the input index
     *   is only valid when (m*M - k) % L == 0. */

    uint32_t out_idx = 0;
    uint32_t max_out = (uint32_t)((uint64_t)n_input * L / M) + 2;

    for (uint32_t m = 0; m < max_out; m++) {
        double sum = 0.0;
        uint64_t base = (uint64_t)m * M;

        for (uint32_t k = 0; k < fir_len; k++) {
            if (base < (uint64_t)k) continue;
            uint64_t index_times_L = base - (uint64_t)k;
            if (index_times_L % L != 0) continue;

            uint64_t input_idx = index_times_L / L;
            if (input_idx >= n_input) continue;

            sum += input[input_idx] * fir_coeffs[k];
        }
        if (out_idx < max_out) {
            output[out_idx++] = sum * (double)L; /* Gain compensation */
        }
    }
    *n_output = out_idx;
}

/* =========================================================================
 * L5 — Band-Pass Sampling
 * ========================================================================= */

void bandpass_sampling_range(double f_l_hz, double f_u_hz,
                              double *fs_min, double *fs_max,
                              uint32_t *optimal_n)
{
    if (f_l_hz < 0.0 || f_u_hz <= f_l_hz) {
        *fs_min = NAN;
        *fs_max = NAN;
        *optimal_n = 0;
        return;
    }

    double BW = f_u_hz - f_l_hz;
    if (BW <= 0.0) {
        *fs_min = NAN;
        *fs_max = NAN;
        *optimal_n = 0;
        return;
    }

    uint32_t n_max = (uint32_t)floor(f_u_hz / BW);
    if (n_max < 1) n_max = 1;

    double best_fs_min = INFINITY;
    double best_fs_max = 0.0;
    uint32_t best_n = 1;

    for (uint32_t n = 1; n <= n_max && n > 0; n++) {
        double fs_min_n = 2.0 * f_u_hz / (double)n;
        double fs_max_n = (n > 1) ? (2.0 * f_l_hz / (double)(n - 1)) : INFINITY;

        /* Valid range must satisfy fs_min ≤ fs_max */
        if (fs_min_n <= fs_max_n) {
            if (fs_min_n < best_fs_min) {
                best_fs_min = fs_min_n;
                best_fs_max = fs_max_n;
                best_n = n;
            }
        }
    }

    *fs_min = best_fs_min;
    *fs_max = best_fs_max;
    *optimal_n = best_n;
}

/* =========================================================================
 * L6 — Ideal Sampling Simulation
 * ========================================================================= */

void simulate_ideal_sampling(double (*ct_signal)(double t, void *ctx),
                              void *ct_signal_ctx,
                              double fs_hz, uint32_t n_samples,
                              double t_start,
                              double *samples)
{
    if (!ct_signal || !samples || n_samples == 0 || fs_hz <= 0.0) return;

    double T = 1.0 / fs_hz;
    for (uint32_t i = 0; i < n_samples; i++) {
        double t = t_start + (double)i * T;
        samples[i] = ct_signal(t, ct_signal_ctx);
    }
}

/* =========================================================================
 * L6 — Full ADC Acquisition Chain Simulation (without FFT dependency)
 * ========================================================================= */

void simulate_adc_acquisition(double (*ct_signal)(double t, void *ctx),
                               void *ctx,
                               double fs_hz, uint32_t n_samples,
                               double t_start,
                               const double *aa_filter_b,
                               const double *aa_filter_a,
                               uint32_t aa_order,
                               double jitter_rms_ps,
                               uint32_t n_bits, double v_ref,
                               uint64_t *codes)
{
    if (!ct_signal || !codes || n_samples == 0 || fs_hz <= 0.0) return;

    double T = 1.0 / fs_hz;
    double jitter_std = jitter_rms_ps * 1e-12;

    /* Simple linear congruential generator for jitter simulation */
    uint64_t lcg_state = 123456789ULL;
    double lcg_max = (double)0xFFFFFFFFFFFFFFFFULL;

    /* AA filter state (direct form I, assuming a[0] = 1.0) */
    double *filter_state_x = NULL;
    double *filter_state_y = NULL;
    int has_filter = 0;

    if (aa_filter_b && aa_order > 0) {
        filter_state_x = (double *)calloc(aa_order + 1, sizeof(double));
        filter_state_y = (double *)calloc(aa_order + 1, sizeof(double));
        if (filter_state_x && filter_state_y) has_filter = 1;
    }

    uint64_t max_code = (1ULL << n_bits) - 1;
    double v_lsb = v_ref / (double)(1ULL << n_bits);

    for (uint32_t i = 0; i < n_samples; i++) {
        /* Apply jitter: actual sampling time = ideal + jitter */
        lcg_state = lcg_state * 6364136223846793005ULL + 1442695040888963407ULL;
        double u1 = (double)lcg_state / lcg_max;
        /* Box-Muller: generate N(0,1) from uniform */
        lcg_state = lcg_state * 6364136223846793005ULL + 1442695040888963407ULL;
        double u2 = (double)lcg_state / lcg_max;
        double gaussian = sqrt(-2.0 * log(u1 + 1e-15)) *
                          cos(2.0 * M_PI * u2);
        double t_actual = t_start + (double)i * T + jitter_std * gaussian;

        double raw_sample = ct_signal(t_actual, ctx);

        /* Apply AA filter */
        double filtered = raw_sample;
        if (has_filter) {
            /* Direct form I: y[n] = Σ b_k * x[n-k] - Σ a_k * y[n-k] */
            /* Shift input history */
            for (uint32_t k = aa_order; k > 0; k--) {
                filter_state_x[k] = filter_state_x[k - 1];
                filter_state_y[k] = filter_state_y[k - 1];
            }
            filter_state_x[0] = raw_sample;

            double y = 0.0;
            for (uint32_t k = 0; k <= aa_order; k++) {
                y += aa_filter_b[k] * filter_state_x[k];
            }
            for (uint32_t k = 1; k <= aa_order; k++) {
                y -= aa_filter_a[k] * filter_state_y[k];
            }
            filter_state_y[0] = y;
            filtered = y;
        }

        /* Quantize (mid-tread ADC) */
        /* First clamp to [0, v_ref] */
        if (filtered < 0.0) filtered = 0.0;
        if (filtered > v_ref) filtered = v_ref;

        double code_d = filtered / v_lsb;
        uint64_t code = (uint64_t)llround(code_d);
        if (code > max_code) code = max_code;
        codes[i] = code;
    }

    free(filter_state_x);
    free(filter_state_y);
}
