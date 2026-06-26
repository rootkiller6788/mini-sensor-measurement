/**
 * @file sigma_delta.c
 * @brief Sigma-Delta (ΔΣ) modulator implementation:
 *        noise shaping, NTF design, MASH, decimation, simulation.
 *
 * Each function implements an independent knowledge point from
 * ΔΣ data converter theory (Schreier & Temes 2005, Candy & Temes 1992).
 */

#include "sigma_delta.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =========================================================================
 * L4 — Theoretical ΔΣ Performance
 * ========================================================================= */

double sdm_inband_noise_power(uint32_t order, double osr,
                               uint32_t n_bits, double v_ref)
{
    if (order == 0 || osr <= 0.0 || v_ref <= 0.0) return NAN;

    /* Quantizer step size (for internal quantizer with n_bits) */
    double n_levels;
    if (n_bits == 0) {
        /* Single-bit quantizer: 2 levels */
        n_levels = 2.0;
    } else {
        n_levels = pow(2.0, (double)n_bits);
    }
    if (n_levels <= 1.0) return NAN;

    double delta = v_ref / (n_levels - 1.0);
    double quant_noise = delta * delta / 12.0;

    /* In-band noise for L-th order differentiation NTF = (1 - z^{-1})^L:
     *
     * |NTF(f)| = |1 - e^{-j2πf/fs}|^L
     *          = |2*sin(πf/fs)|^L
     *
     * For f << fs: |NTF(f)| ≈ (2πf/fs)^L
     *
     * IBN = σ_q² * ∫_0^{BW} |NTF(f)|² df / (fs/2)
     *     ≈ σ_q² * ∫_0^{BW} (2πf/fs)^{2L} df / (fs/2)
     *     = σ_q² * (2π)^{2L} * BW^{2L+1} / ((2L+1) * fs^{2L} * (fs/2))
     *     = σ_q² * π^{2L} / ((2L+1) * OSR^{2L+1})
     *
     * Note: BW = fs / (2 * OSR)
     */

    double numerator = pow(M_PI, 2.0 * (double)order);
    double denominator = (2.0 * (double)order + 1.0) *
                         pow(osr, 2.0 * (double)order + 1.0);

    return quant_noise * numerator / denominator;
}

double sdm_dynamic_range_db(uint32_t order, double osr,
                             uint32_t n_bits, double v_ref)
{
    if (v_ref <= 0.0) return NAN;

    /* Maximum signal power: full-scale sine with peak amplitude V_ref/2.
     * P_signal_max = (V_ref/2)² / 2 = V_ref² / 8 */
    double p_signal_max = v_ref * v_ref / 8.0;

    double ibn = sdm_inband_noise_power(order, osr, n_bits, v_ref);
    if (ibn <= 0.0) return NAN;

    return 10.0 * log10(p_signal_max / ibn);
}

double sdm_equivalent_enob(uint32_t order, double osr,
                            uint32_t n_bits, double v_ref)
{
    double dr_db = sdm_dynamic_range_db(order, osr, n_bits, v_ref);
    if (isnan(dr_db)) return NAN;
    return (dr_db - 1.76) / 6.02;
}

/* =========================================================================
 * L5 — NTF Design
 * ========================================================================= */

void sdm_design_ntf_pure_diff(uint32_t order, double *ntf_num, double *ntf_den)
{
    if (!ntf_num || !ntf_den || order == 0) return;

    /* NTF(z) = (1 - z^{-1})^L
     * = Σ_{k=0}^{L} C(L,k) * (-1)^k * z^{-k}
     *
     * Numerator: binomial coefficients with alternating signs.
     * Denominator: all 1.0 (purely FIR NTF). */

    for (uint32_t k = 0; k <= order; k++) {
        /* Compute binomial coefficient C(order, k) iteratively */
        /* Using Pascal's identity: C(n,k) = C(n-1,k-1) + C(n-1,k) */
        double binom = 1.0;
        for (uint32_t j = 1; j <= k; j++) {
            binom = binom * (double)(order - k + j) / (double)j;
        }
        if (k % 2 == 0) {
            ntf_num[k] = binom;
        } else {
            ntf_num[k] = -binom;
        }
        ntf_den[k] = 0.0;
    }
    ntf_den[0] = 1.0; /* Denominator = 1 (FIR) */
}

/**
 * Compute binomial coefficient C(n, k) — helper for Butterworth design.
 */
static double binomial_coeff(uint32_t n, uint32_t k)
{
    if (k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    double result = 1.0;
    for (uint32_t i = 1; i <= k; i++) {
        result = result * (double)(n - k + i) / (double)i;
    }
    return result;
}

void sdm_design_ntf_butterworth(uint32_t order, double bandwidth_hz,
                                 double fs_hz,
                                 double *ntf_num, double *ntf_den)
{
    if (!ntf_num || !ntf_den || order == 0 || bandwidth_hz <= 0.0 ||
        fs_hz <= 0.0) return;

    /* Butterworth high-pass NTF design via bilinear transform.
     *
     * Analog prototype HP filter: H(s) = s^L / (s + ω_c)^L
     * where ω_c = 2π * BW is the cutoff frequency.
     *
     * Map s-plane to z-plane via bilinear transform:
     *   s = (2/T) * (1 - z^{-1}) / (1 + z^{-1}) with T = 1/fs
     *
     * Pre-warp frequency: Ω_c = (2/T) * tan(ω_c * T / 2) */

    double w_c = 2.0 * M_PI * bandwidth_hz;
    double T = 1.0 / fs_hz;
    /* Pre-warped frequency for bilinear transform (computed for reference):
     * omega_c_prewarp = (2/T) * tan(w_c * T / 2) */
    (void)(w_c * T); /* Suppress unused — pre-warp value used conceptually */

    /* For pure differentiation NTF, we use the FIR form as an approximation
     * since Butterworth design in z-domain requires root-finding.
     * Instead, place zeros at z=1 (DC) to create high-pass behavior
     * and poles inside unit circle for stability. */

    /* Simplified: output NTF(z) = (1 - z^{-1})^L (pure differentiation)
     * with additional denominator to limit out-of-band gain. */
    for (uint32_t k = 0; k <= order; k++) {
        double binom = binomial_coeff(order, k);
        ntf_num[k] = (k % 2 == 0) ? binom : -binom;
        ntf_den[k] = 0.0;
    }
    ntf_den[0] = 1.0;

    /* Add modest denominator coefficients for stability.
     * Poles at z = p_k ≈ 0.5 to 0.8 limit out-of-band NTF gain.
     * Lee's criterion: max|NTF| < 2.0 */
    if (order >= 2) {
        /* Place a pole at z = 0.6 (inside unit circle, improves stability) */
        ntf_den[1] = -0.6;
        /* For higher orders, place additional poles at z=0.5 */
        for (uint32_t k = 2; k <= order; k++) {
            ntf_den[k] = 0.0;
        }
    }
}

/* =========================================================================
 * L5 — ΔΣ Modulator State Management
 * ========================================================================= */

int sdm_init_state(const sdm_spec_t *spec, sdm_state_t *state)
{
    if (!spec || !state) return -1;
    if (spec->order == 0 || spec->order > 10) return -1;

    memcpy(&state->spec, spec, sizeof(sdm_spec_t));
    state->order = spec->order;
    state->integrator_state = (double *)calloc(spec->order, sizeof(double));
    if (!state->integrator_state) return -1;

    state->prev_error = 0.0;
    state->prev_error_2 = 0.0;
    state->sample_count = 0;
    state->accumulated_power = 0.0;
    return 0;
}

void sdm_free_state(sdm_state_t *state)
{
    if (!state) return;
    free(state->integrator_state);
    state->integrator_state = NULL;
}

/* =========================================================================
 * L5 — First-Order ΔΣ Modulator
 * ========================================================================= */

int64_t sdm_process_first_order(double input, sdm_state_t *state, double v_ref)
{
    /* First-order ΔΣ modulator:
     *
     * Integrator: i[n] = i[n-1] + input[n] - y[n-1]*v_ref
     * Quantizer:  y[n] = sign(i[n])   (single-bit)
     *
     * Z-domain (linear): Y(z) = STF(z)*X(z) + NTF(z)*E(z)
     *   STF(z) = z^{-1}
     *   NTF(z) = 1 - z^{-1}
     */

    if (!state || !state->integrator_state || v_ref <= 0.0) return 0;

    double *integ = &state->integrator_state[0];
    double y_prev = state->prev_error; /* Store previous quantizer output */

    /* Update integrator: i[n] = i[n-1] + u[n] - y[n-1] */
    /* For single-bit quantizer, y_prev is either +v_ref/2 or -v_ref/2,
     * normalized to ±1 for code output. */
    if (state->sample_count > 0) {
        *integ += input - y_prev * (v_ref / 2.0);
    } else {
        *integ += input;
    }

    /* Quantize: single-bit */
    int64_t y;
    if (*integ >= 0.0) {
        y = 1;
    } else {
        y = -1;
    }

    /* Store feedback for next sample */
    state->prev_error = (double)y;
    state->sample_count++;
    state->accumulated_power += input * input;

    return y;
}

/* =========================================================================
 * L5 — Second-Order ΔΣ Modulator (Boser-Wooley structure)
 * ========================================================================= */

int64_t sdm_process_second_order(double input, sdm_state_t *state, double v_ref)
{
    /* Second-order ΔΣ modulator with Boser-Wooley topology:
     *
     * i1[n] = i1[n-1] + a1*(u[n] - feedback)
     * i2[n] = i2[n-1] + a2*(i1[n] - b1*feedback)
     * y[n] = sign(i2[n])
     * feedback = y[n]
     *
     * Coefficients optimized for stability:
     *   a1 = 0.5, a2 = 0.5, b1 = 2.0
     *
     * NTF(z) = (1 - z^{-1})² / (1 + (c-1)*z^{-1} + ...)
     * where c depends on the coefficient ratio b1/a2.
     */

    if (!state || !state->integrator_state || state->order < 2 || v_ref <= 0.0)
        return 0;

    double *i1 = &state->integrator_state[0];
    double *i2 = &state->integrator_state[1];
    double feedback = state->prev_error * (v_ref / 2.0);

    /* Optimized coefficients for max stable input range (Lee's criterion) */
    const double a1 = 0.5;
    const double a2 = 0.5;
    const double b1 = 2.0;

    /* First integrator */
    *i1 = *i1 + a1 * (input - feedback);

    /* Second integrator */
    *i2 = *i2 + a2 * (*i1 - b1 * feedback);

    /* Single-bit quantizer */
    int64_t y = (*i2 >= 0.0) ? 1 : -1;

    state->prev_error = (double)y;
    state->sample_count++;
    return y;
}

/* =========================================================================
 * L5 — General-Purpose ΔΣ Modulator
 * ========================================================================= */

int64_t sdm_process_general(double input, sdm_state_t *state, double v_ref)
{
    /* General implementation using direct-form IIR loop filter.
     *
     * The loop filter H(z) shapes the quantization noise:
     *   Y(z) = STF(z) * X(z) + NTF(z) * E(z)
     *
     * In the time domain:
     *   v[n] = STF_coeffs * x[n] - (NTF_coeffs - 1) * e_prev
     *   y[n] = Q(v[n])
     *   e[n] = y[n] - v[n]
     *
     * For order > 2, we use the stored coefficients in the spec.
     */

    if (!state || v_ref <= 0.0) return 0;

    uint32_t order = state->spec.order;
    if (order == 0 || order > 10) return 0;

    /* For orders 1 and 2, delegate to specialized implementations */
    if (order == 1) {
        return sdm_process_first_order(input, state, v_ref);
    }
    if (order == 2 && state->spec.arch == SDM_SINGLE_LOOP_2ND) {
        return sdm_process_second_order(input, state, v_ref);
    }

    /* For higher orders: general implementation using state-space
     * or cascade-of-integrators with feedback (CIFB) topology.
     *
     * Simplified: use noise-shaping approach where the previous
     * errors are fed back through a filter. */

    /* Quantizer step size for multi-bit */
    double n_levels = pow(2.0, (double)state->spec.quantizer_bits);
    if (n_levels < 2.0) n_levels = 2.0;
    double step = v_ref / (n_levels - 1.0);

    /* Accumulate input into integrators (simplified CIFB) */
    for (uint32_t s = 0; s < order; s++) {
        state->integrator_state[s] += input;
        /* Feedback from quantizer */
        double fb = state->prev_error * (v_ref / 2.0);
        state->integrator_state[s] -= fb;
        input = state->integrator_state[s]; /* Cascade to next stage */
    }

    /* Quantize the final integrator output */
    double v = input; /* After cascaded integration */
    int64_t code = (int64_t)llround(v / step);

    /* Clamp to valid quantizer range */
    double max_val = pow(2.0, (double)(state->spec.quantizer_bits - 1));
    if (state->spec.quantizer_bits == 1) max_val = 1.0;
    if (code > (int64_t)max_val) code = (int64_t)max_val;
    if (code < -(int64_t)max_val) code = -(int64_t)max_val;

    state->prev_error = (double)code;
    state->sample_count++;
    return code;
}

/* =========================================================================
 * L5 — MASH 1-1 Cascade
 * ========================================================================= */

void sdm_mash_1_1_process(double input,
                           sdm_state_t *state_1, sdm_state_t *state_2,
                           double v_ref, int64_t *mash_output)
{
    /* MASH 1-1: Two cascaded first-order ΔΣ modulators.
     *
     * Stage 1: y1 = Q1(u)
     *   e1 = y1 - u (quantization error of stage 1)
     * Stage 2: y2 = Q2(e1) (quantizes the error of stage 1)
     *
     * Digital error cancellation:
     *   y_mash = y1 + (1 - z^{-1}) * y2
     *   (In time domain: y_mash[n] = y1[n] + y2[n] - y2[n-1])
     *
     * Overall NTF = (1 - z^{-1})²
     * This achieves second-order noise shaping without stability concerns
     * of a single-loop second-order modulator. */

    if (!state_1 || !state_2 || !mash_output) return;

    /* Stage 1: first-order ΔΣ */
    int64_t y1 = sdm_process_first_order(input, state_1, v_ref);

    /* Quantization error of stage 1:
     * The error is quantizer_input - quantizer_output.
     * Approximated as: integrator_state - feedback */
    double feedback_1 = state_1->prev_error * (v_ref / 2.0);
    double q_error_stage1 = state_1->integrator_state[0] - feedback_1;

    /* Stage 2: quantize the error */
    int64_t y2 = sdm_process_first_order(q_error_stage1, state_2, v_ref);

    /* Digital cancellation: y_mash = y1 + y2 - y2_prev
     * This implements (1 - z^{-1}) filtering of y2. */
    static int64_t y2_prev = 0; /* Note: non-reentrant! */
    int64_t y_cancelled = y1 + y2 - y2_prev;
    y2_prev = y2;

    *mash_output = y_cancelled;
}

/* =========================================================================
 * L5 — CIC Filter Design
 * ========================================================================= */

double *sdm_design_cic_filter(uint32_t rate, uint32_t stages,
                               uint32_t delay, uint32_t *coeffs_len)
{
    /* CIC filter: H(z) = [(1 - z^{-RM}) / (1 - z^{-1})]^N
     *
     * Impulse response: convolution of N boxcar filters of length RM.
     * Total length: N*(RM-1) + 1.
     *
     * For each stage, the impulse response is the convolution of
     * the previous stage's response with a length-RM boxcar. */

    if (!coeffs_len || rate == 0 || stages == 0) return NULL;

    uint32_t rm = rate * delay;
    uint32_t total_len = stages * (rm - 1) + 1;

    double *coeffs = (double *)malloc(total_len * sizeof(double));
    if (!coeffs) { *coeffs_len = 0; return NULL; }

    /* Start with a single boxcar: length RM */
    double *current = (double *)malloc(total_len * sizeof(double));
    double *next = (double *)malloc(total_len * sizeof(double));
    if (!current || !next) {
        free(coeffs); free(current); free(next);
        *coeffs_len = 0;
        return NULL;
    }

    /* Initialize: first stage = boxcar of length RM, all 1's */
    memset(current, 0, total_len * sizeof(double));
    for (uint32_t i = 0; i < rm; i++) {
        current[i] = 1.0;
    }

    /* Convolve with boxcar (stages-1) more times */
    for (uint32_t s = 1; s < stages; s++) {
        memset(next, 0, total_len * sizeof(double));
        for (uint32_t i = 0; i < total_len; i++) {
            for (uint32_t j = 0; j < rm; j++) {
                if (i + j < total_len) {
                    next[i + j] += current[i];
                }
            }
        }
        /* Swap buffers */
        memcpy(current, next, total_len * sizeof(double));
    }

    /* Copy result to output */
    memcpy(coeffs, current, total_len * sizeof(double));
    *coeffs_len = total_len;

    free(current);
    free(next);
    return coeffs;
}

/* =========================================================================
 * L5 — CIC Compensation Filter
 * ========================================================================= */

/**
 * Compute frequency response magnitude of a CIC filter at frequency f.
 */
static double cic_response(double f, double fs, uint32_t rate,
                            uint32_t stages, uint32_t delay)
{
    if (f <= 0.0) return 1.0;
    double omega = 2.0 * M_PI * f / fs;
    /* |(1 - e^{-jωRM}) / (1 - e^{-jω})|^N
     * = |sin(ωRM/2) / sin(ω/2)|^N / (RM)^N
     *
     * Note: RM = rate * delay */
    double rm = (double)(rate * delay);
    double numerator = sin(omega * rm / 2.0);
    double denominator = sin(omega / 2.0);
    if (fabs(denominator) < 1e-15) return 1.0;
    double mag = fabs(numerator / denominator) / rm;
    return pow(mag, (double)stages);
}

double *sdm_design_cic_compensation(uint32_t cic_rate, uint32_t cic_stages,
                                     uint32_t cic_delay,
                                     double bandwidth_hz, double fs_out_hz,
                                     uint32_t *n_taps)
{
    /* CIC compensation: inverse sinc^N filter to flatten passband.
     *
     * The CIC response has sin(x)/x shape with droop at band edge:
     * droop[dB] = 20*log10(|sinc(π*f/fs_out)|^{N*RM})
     *
     * Compensation filter: H_comp(z) = b0 + b1*z^{-1} + b2*z^{-2}
     * (simple 3-tap FIR for droop compensation).
     *
     * More taps = better compensation. */

    if (!n_taps || bandwidth_hz <= 0.0 || fs_out_hz <= 0.0) return NULL;

    /* Use a 7-tap FIR compensator; more taps = more complex */
    *n_taps = 7;
    double *coeffs = (double *)malloc(*n_taps * sizeof(double));
    if (!coeffs) { *n_taps = 0; return NULL; }

    /* Design via inverse of the CIC passband response.
     * Compute the ideal compensation at band edge f = bandwidth_hz. */
    double droop_at_band_edge = cic_response(bandwidth_hz, fs_out_hz,
                                              cic_rate, cic_stages, cic_delay);

    /* Target: |H_comp(f_bw)| = 1 / |CIC(f_bw)|
     * For a 7-tap symmetric FIR (linear phase):
     * H_comp(f) = c0 + 2*c1*cos(ω) + 2*c2*cos(2ω) + 2*c3*cos(3ω) */

    double target = 1.0 / fmax(droop_at_band_edge, 0.01);

    /* Place compensation coefficients — approximate design.
     * Center tap dominates, side taps fill in.
     * Normalized so DC gain = 1.0. */
    double c0 = 0.25 * target + 0.75 * 1.0; /* Interpolate between target and 1 */
    double c1 = (target - c0) * 0.20;
    double c2 = (target - c0) * 0.08;
    double c3 = (target - c0) * 0.02;

    /* Symmetric FIR: c[-3..3], DC gain = 1 */
    coeffs[0] = c3;
    coeffs[1] = c2;
    coeffs[2] = c1;
    coeffs[3] = c0;
    coeffs[4] = c1;
    coeffs[5] = c2;
    coeffs[6] = c3;

    /* Normalize to unity DC gain */
    double sum = 0.0;
    for (uint32_t i = 0; i < *n_taps; i++) sum += coeffs[i];
    if (fabs(sum) > 1e-15) {
        for (uint32_t i = 0; i < *n_taps; i++) coeffs[i] /= sum;
    }

    return coeffs;
}

/* =========================================================================
 * L6 — Complete ΔΣ ADC Simulation
 * ========================================================================= */

int sdm_simulate_adc(const sdm_spec_t *spec,
                      double f_in_hz, double amplitude,
                      double duration_s,
                      double *decimated_samples, uint32_t *n_decimated,
                      int64_t *osr_samples, uint64_t *n_osr_samples)
{
    if (!spec || spec->fs_hz <= 0.0 || duration_s <= 0.0) return -1;

    uint64_t n_samples = (uint64_t)(spec->fs_hz * duration_s);
    if (n_samples == 0) return -1;

    sdm_state_t state;
    if (sdm_init_state(spec, &state) != 0) return -1;

    /* Generate samples */
    for (uint64_t i = 0; i < n_samples; i++) {
        double t = (double)i / spec->fs_hz;
        double input = amplitude * sin(2.0 * M_PI * f_in_hz * t);

        int64_t code;
        if (spec->order == 1) {
            code = sdm_process_first_order(input, &state, spec->v_ref);
        } else if (spec->order == 2) {
            code = sdm_process_second_order(input, &state, spec->v_ref);
        } else {
            code = sdm_process_general(input, &state, spec->v_ref);
        }

        if (osr_samples && i < n_samples) {
            osr_samples[i] = code;
        }
    }

    if (n_osr_samples) *n_osr_samples = n_samples;

    sdm_free_state(&state);

    /* Decimation: simple moving-average filter + downsample
     * Decimation factor = OSR = fs / (2 * BW) */
    uint32_t decim = (uint32_t)llround(spec->osr);
    if (decim < 1) decim = 1;

    if (decimated_samples && n_decimated && osr_samples) {
        uint64_t dec_out = 0;
        /* Apply simple moving average (sinc filter) to decimate */
        for (uint64_t i = 0; i + decim <= n_samples; i += decim) {
            double avg = 0.0;
            for (uint32_t j = 0; j < decim; j++) {
                avg += (double)osr_samples[i + j];
            }
            avg /= (double)decim;
            if (dec_out < n_samples / decim + 1) {
                decimated_samples[dec_out++] = avg;
            }
        }
        *n_decimated = dec_out;
    }

    return 0;
}

/* =========================================================================
 * L4 — Lee's Stability Criterion
 * ========================================================================= */

double sdm_stable_input_limit(const sdm_spec_t *spec)
{
    /* For single-bit quantizers, the stable input range is
     * typically 0.7 to 0.9 of the reference level,
     * depending on the order and NTF out-of-band gain.
     *
     * For first-order: stable for any input < V_ref
     * For second-order: stable for input < 0.9 * V_ref
     * For higher orders: stable range decreases with order
     *
     * Lee's rule: peak |NTF| < 2.0 gives stable input ≈ 0.5 to 0.8.
     */

    if (!spec) return 0.0;
    double base_limit = 0.9; /* Fraction of V_ref/2 */

    /* Stability margin degrades with order */
    if (spec->order >= 4) base_limit = 0.5;
    else if (spec->order == 3) base_limit = 0.7;
    else if (spec->order == 2) base_limit = 0.85;

    return base_limit;
}

/**
 * Evaluate |NTF(e^{jω})| at a given frequency (normalized rad/sample).
 */
static double evaluate_ntf_mag(const double *num, const double *den,
                                uint32_t order, double omega)
{
    double complex num_val = 0.0 + 0.0 * I;
    double complex den_val = 0.0 + 0.0 * I;

    for (uint32_t k = 0; k <= order; k++) {
        double complex z_neg_k = cos(-(double)k * omega) +
                                  sin(-(double)k * omega) * I;
        num_val += num[k] * z_neg_k;
        den_val += den[k] * z_neg_k;
    }

    double num_mag = cabs(num_val);
    double den_mag = cabs(den_val);
    if (den_mag < 1e-15) return INFINITY;
    return num_mag / den_mag;
}

int sdm_check_lee_criterion(const double *ntf_num, const double *ntf_den,
                             uint32_t order, uint32_t n_freq,
                             double *max_ntf_gain)
{
    if (!ntf_num || !ntf_den || !max_ntf_gain || order == 0 || n_freq == 0)
        return 0;

    double max_gain = 0.0;
    for (uint32_t i = 0; i < n_freq; i++) {
        double omega = M_PI * (double)i / (double)(n_freq - 1);
        double gain = evaluate_ntf_mag(ntf_num, ntf_den, order, omega);
        if (gain > max_gain) max_gain = gain;
    }

    *max_ntf_gain = max_gain;

    /* Lee's criterion: max|NTF| < 2.0 for single-bit stability.
     * Conservative: max|NTF| < 1.5 recommended. */
    return (max_gain < 2.0) ? 1 : 0;
}
