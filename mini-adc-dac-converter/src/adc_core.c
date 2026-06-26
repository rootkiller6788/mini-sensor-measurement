/**
 * @file adc_core.c
 * @brief Core ADC/DAC implementation: ideal behavior, calibration, linearity.
 *
 * Knowledge Coverage:
 *   L2 - Nyquist-rate ADC modeling
 *   L4 - Ideal SNR formula derivation and verification
 *   L5 - DNL/INL computation algorithms
 *   L6 - Transfer characteristic modeling
 *
 * Each function implements a distinct, independent knowledge point
 * derived from IEEE Std 1241-2010 and graduate converter textbooks.
 */

#include "adc_dac_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

/* Boltzmann constant [J/K] */
#define K_BOLTZMANN 1.380649e-23

/* =========================================================================
 * L4 — Ideal ADC SNR: Derivation and Implementation
 * ========================================================================= */

/**
 * Ideal SNR for N-bit ADC with full-scale sinusoidal input.
 *
 * Derivation:
 *   Quantization step: q = V_ref / 2^N
 *   Quantization noise power (uniform [-q/2, q/2]): σ_q² = q² / 12
 *   Full-scale sine amplitude: A = V_ref/2 = 2^{N-1} * q
 *   Signal power (RMS): P_s = A² / 2 = 2^{2N-3} * q²
 *   SNR = 10*log10(P_s / σ_q²) = 10*log10(3 * 2^{2N-1})
 *       = 10*log10(3) + (2N-1)*10*log10(2)
 *       = 4.7712 + (2N-1)*3.0103
 *       = 6.0206*N + 1.7609
 *
 * Verified against: IEEE Std 1241-2010 Eq. (36).
 * This is invariant: for any V_ref, the ratio depends only on N.
 */
double adc_ideal_snr_db(uint32_t n_bits)
{
    if (n_bits == 0 || n_bits > 32) {
        return NAN;
    }
    /* Derivation inline for educational value:
     * SNR = 10*log10(1.5) + 20*N*log10(2)
     *     = 1.7609 + 6.0206*N
     */
    return 6.0206 * (double)n_bits + 1.7609;
}

/* =========================================================================
 * L1 — SNR from power ratio
 * ========================================================================= */

double adc_snr_from_powers(double p_signal, double p_noise)
{
    if (p_signal <= 0.0 || p_noise <= 0.0) {
        return NAN;
    }
    return 10.0 * log10(p_signal / p_noise);
}

/* =========================================================================
 * L2 — ENOB from SINAD
 * ========================================================================= */

double adc_enob_from_sinad(double sinad_db)
{
    if (sinad_db <= 1.76) {
        return NAN;
    }
    return (sinad_db - 1.76) / 6.02;
}

double adc_enob_from_sndr(double sndr_db)
{
    /* SNDR and SINAD are mathematically equivalent for the ENOB formula. */
    return adc_enob_from_sinad(sndr_db);
}

/* =========================================================================
 * L4 — Jitter-Limited SNR
 * ========================================================================= */

double adc_jitter_snr_limit(double freq_hz, double jitter_rms_ps)
{
    if (freq_hz <= 0.0 || jitter_rms_ps <= 0.0) {
        return NAN;
    }
    /* For sinusoid Vin(t) = A*sin(2πft), the slope is dV/dt = 2πfA*cos(2πft).
     * Max error magnitude: 2πfA * t_jitter
     * RMS error: 2πfA * t_jitter / sqrt(2)
     * SNR = A²/2 / ((2πfA * t_jitter)²/2) = 1 / (2πf * t_jitter)²
     * SNR_dB = -20*log10(2π * f * t_jitter)
     */
    double jitter_rms_s = jitter_rms_ps * 1e-12;
    double arg = 2.0 * M_PI * freq_hz * jitter_rms_s;
    if (arg <= 0.0) return INFINITY;
    return -20.0 * log10(arg);
}

/* =========================================================================
 * L2 — Total Input-Referred Noise
 * ========================================================================= */

double adc_total_input_noise(double v_lsb, double cap_pF, double temp_celsius)
{
    if (v_lsb <= 0.0 || cap_pF <= 0.0) {
        return NAN;
    }
    /* Quantization noise RMS: V_LSB / sqrt(12) */
    double v_quant_rms = v_lsb / sqrt(12.0);

    /* Thermal (kT/C) noise RMS.
     * Temperature in Kelvin = Celsius + 273.15. */
    double temp_kelvin = temp_celsius + 273.15;
    if (temp_kelvin <= 0.0) temp_kelvin = 300.0; /* clamp to room temp if nonsense */
    double cap_F = cap_pF * 1e-12;
    double v_thermal_rms = sqrt(K_BOLTZMANN * temp_kelvin / cap_F);

    /* Total RMS: root-sum-of-squares (independent noise sources) */
    double v_total = sqrt(v_quant_rms * v_quant_rms +
                          v_thermal_rms * v_thermal_rms);
    return v_total;
}

/* =========================================================================
 * L5 — Two-Point Calibration
 * ========================================================================= */

uint64_t adc_two_point_calibrate(uint64_t raw_code,
                                  uint64_t code_at_min_input,
                                  uint64_t code_at_max_input,
                                  uint32_t n_bits)
{
    if (n_bits == 0 || n_bits > 32) return raw_code;
    uint64_t max_code = (1ULL << n_bits) - 1;

    if (code_at_max_input <= code_at_min_input) return raw_code;
    if (code_at_min_input >= max_code || code_at_max_input > max_code)
        return raw_code;

    /* Compute correction:
     *   offset = code_at_min_input
     *   gain_correction = (max_code) / (code_at_max_input - code_at_min_input)
     *   corrected = (raw - offset) * gain_correction
     */
    double offset_code = (double)code_at_min_input;
    double gain_factor = (double)max_code /
                         (double)(code_at_max_input - code_at_min_input);
    double corrected_d = ((double)raw_code - offset_code) * gain_factor;
    if (corrected_d < 0.0) corrected_d = 0.0;
    if (corrected_d > (double)max_code) corrected_d = (double)max_code;
    return (uint64_t)llround(corrected_d);
}

void adc_calibration_params(uint64_t code_at_min,
                            uint64_t code_at_max,
                            uint32_t n_bits,
                            double *offset_code,
                            double *gain_factor)
{
    if (n_bits == 0 || n_bits > 32 ||
        code_at_max <= code_at_min) {
        *offset_code = 0.0;
        *gain_factor = 1.0;
        return;
    }
    uint64_t max_code = (1ULL << n_bits) - 1;
    *offset_code = (double)code_at_min;
    *gain_factor = (double)max_code /
                   (double)(code_at_max - code_at_min);
}

/* =========================================================================
 * L5 — ADC Code ↔ Voltage Conversion
 * ========================================================================= */

double adc_code_to_voltage(uint64_t code, double v_ref, uint32_t n_bits)
{
    if (n_bits == 0 || n_bits > 32 || v_ref <= 0.0) {
        return NAN;
    }
    uint64_t max_code = (1ULL << n_bits) - 1;
    if (code > max_code) code = max_code;
    return ((double)code * v_ref) / (double)(1ULL << n_bits);
}

uint64_t adc_voltage_to_code(double v_in, double v_ref, uint32_t n_bits)
{
    if (n_bits == 0 || n_bits > 32 || v_ref <= 0.0) {
        return 0;
    }
    uint64_t max_code = (1ULL << n_bits) - 1;
    double step = v_ref / (double)(1ULL << n_bits);

    /* Mid-tread quantizer: code = floor(V_in / step + 0.5) */
    double code_d = floor(v_in / step + 0.5);
    if (code_d < 0.0) code_d = 0.0;
    if (code_d > (double)max_code) code_d = (double)max_code;
    return (uint64_t)code_d;
}

/* =========================================================================
 * L5 — DNL Computation (IEEE 1241 §4.3.3)
 * ========================================================================= */

/**
 * Compute DNL from measured code transition voltages.
 *
 * Algorithm (ramp histogram method):
 *   For each transition k = 1..2^N-1:
 *     1. Compute actual code width: W[k] = T[k+1] - T[k]
 *     2. Compute average LSB: V_LSB_avg = (T[2^N] - T[0]) / 2^N
 *     3. DNL[k] = W[k] / V_LSB_avg - 1
 *
 * DNL is unitless, expressed in LSB.
 * DNL = 0 → ideal. DNL > 0 → wide code. DNL < 0 → narrow code.
 * DNL < -1 → missing code.
 */
void adc_compute_dnl(const double *transitions, uint32_t n_bits,
                     double *dnl_out)
{
    if (!transitions || !dnl_out || n_bits == 0 || n_bits > 32) {
        return;
    }
    uint64_t n_codes = 1ULL << n_bits;
    uint64_t n_transitions = n_codes + 1; /* T[0]..T[2^N] */

    /* End-point correction: remove offset and gain errors.
     * Compute average LSB from end points. */
    double v_lsb_avg = (transitions[n_transitions - 1] - transitions[0]) /
                       (double)n_codes;
    if (v_lsb_avg <= 0.0) {
        /* Degenerate case: flat transfer curve. Set all DNL = 0. */
        for (uint64_t i = 0; i < n_codes - 1; i++) {
            dnl_out[i] = 0.0;
        }
        return;
    }

    for (uint64_t k = 0; k < n_codes - 1; k++) {
        double code_width = transitions[k + 1] - transitions[k];
        dnl_out[k] = (code_width / v_lsb_avg) - 1.0;
    }

    /* DNL[0] is for code 1 (code 0 has no DNL definition) */
}

/* =========================================================================
 * L5 — INL from DNL
 * ========================================================================= */

/**
 * Compute INL by cumulative summation of DNL.
 *
 * The INL at a given code is the cumulative sum of DNL values
 * up to that code, representing the deviation of the actual
 * transfer function from a best-fit straight line.
 *
 * INL[0] = 0  (by definition: end-point corrected)
 * INL[k] = sum_{i=0}^{k-1} DNL[i]
 *
 * Relationship to integral:
 *   INL(code) = (V_actual(code) - V_ideal(code)) / V_LSB
 *
 * Interpreting INL:
 *   |INL| < 0.5 LSB → monotonic for DAC
 *   |INL| < 1.0 LSB → no missing codes for ADC
 */
void adc_compute_inl_from_dnl(const double *dnl, uint32_t n_bits,
                               double *inl_out)
{
    if (!dnl || !inl_out || n_bits == 0 || n_bits > 32) {
        return;
    }
    uint64_t n_codes = 1ULL << n_bits;

    inl_out[0] = 0.0; /* End-point INL is zero by definition. */
    double running_sum = 0.0;
    for (uint64_t k = 1; k < n_codes; k++) {
        /* dnl[k-1] corresponds to the code width between transition k-1 and k */
        running_sum += dnl[k - 1];
        inl_out[k] = running_sum;
    }
}

/* =========================================================================
 * L3 — Harmonic Distortion Model
 * ========================================================================= */

/**
 * Evaluate the polynomial transfer characteristic model.
 *
 *   V_out(V_in) = gain * V_in + offset + Σ a_i * Chebyshev_i(V_in/V_fs)
 *
 * Using Chebyshev polynomials of the first kind captures harmonic
 * distortion in a numerically stable manner.
 *
 * T_0(x) = 1, T_1(x) = x, T_{n+1}(x) = 2x*T_n(x) - T_{n-1}(x)
 */
static double chebyshev_T(int n, double x)
{
    if (n == 0) return 1.0;
    if (n == 1) return x;
    double T_n_minus_2 = 1.0;
    double T_n_minus_1 = x;
    double T_n = 0.0;
    for (int i = 2; i <= n; i++) {
        T_n = 2.0 * x * T_n_minus_1 - T_n_minus_2;
        T_n_minus_2 = T_n_minus_1;
        T_n_minus_1 = T_n;
    }
    return T_n;
}

/**
 * Compute THD from harmonic amplitudes in the transfer function.
 *
 * For a sinusoidal input V_in(t) = A*sin(2πf*t), the output
 * contains harmonics at integer multiples of f.
 *
 * Harmonic amplitudes are computed from the polynomial coefficients
 * using the Chebyshev expansion of sin^n(x).
 *
 * This function evaluates THD analytically (no FFT needed).
 */
double transfer_char_thd(const transfer_char_t *tc, double amplitude)
{
    if (!tc || amplitude <= 0.0 || tc->harmonic_order == 0) {
        return 0.0;
    }
    /* Fundamental amplitude through linear + nonlinear path:
     *
     * For input x(t) = A·sin(θ) where θ = 2πf·t:
     *   y(t) = gain·x + offset + Σ a_k·[x(t)]^k
     *
     * Using Chebyshev polynomials of the first kind:
     *   sin^k(θ) can be expressed as a linear combination of
     *   sin(θ), sin(2θ), ..., sin(kθ) for odd k, and
     *   DC + cos(2θ) + cos(4θ) + ... for even k.
     *
     * Specifically: sin^n(θ) = Σ c_m·T_m(sin(θ))
     * where c_m are Fourier coefficients of sin^n.
     *
     * For simplicity, the dominant term is:
     *   sin^k(θ) ≈ (1/2^{k-1}) · sin(kθ) + lower-order terms. */

    /* Fundamental: linear term + odd-order nonlinear contributions */
    double fundamental_amplitude = fabs(tc->gain) * amplitude;
    if (tc->harmonic_coeffs) {
        /* Odd-order terms also contribute to fundamental via sin(θ) component */
        for (uint32_t k = 3; k <= tc->harmonic_order; k += 2) {
            double a_k = tc->harmonic_coeffs[k - 1];
            if (fabs(a_k) < 1e-15) continue;
            /* Chebyshev coefficient: sin^k(θ) has fundamental component
             * proportional to 1/2^{k-1} * C(k, (k-1)/2) */
            double binom = 1.0;
            uint32_t mid = (k - 1) / 2;
            for (uint32_t j = 1; j <= mid; j++) {
                binom = binom * (double)(k - mid + j) / (double)j;
            }
            fundamental_amplitude += fabs(a_k) * binom *
                pow(amplitude, (double)k) / pow(2.0, (double)(k - 1));
        }
    }

    if (fundamental_amplitude <= 0.0) return 0.0;

    double harmonic_power_sum = 0.0;
    for (uint32_t k = 2; k <= tc->harmonic_order; k++) {
        if (!tc->harmonic_coeffs) break;
        double a_k = tc->harmonic_coeffs[k - 1];
        if (fabs(a_k) < 1e-15) continue;

        /* Evaluate harmonic amplitude using Chebyshev decomposition.
         * T_n(sin(θ)) = (-1)^{(n-1)/2}·sin(nθ)  for odd n.
         *
         * The k-th order term generates harmonics up to the k-th.
         * Dominant contribution at harmonic k: */
        double amplitude_k = fabs(a_k) * pow(amplitude, (double)k) /
                             pow(2.0, (double)(k - 1));

        /* Use Chebyshev evaluation to verify harmonic structure:
         * For k=2: sin²(θ) = ½ - ½·T₂(sin(θ)) → 2nd harmonic
         * For k=3: sin³(θ) = ¾·sin(θ) - ¼·sin(3θ) → 3rd harmonic
         * Chebyshev: T₃(x) = 4x³ - 3x  → captures 3rd harmonic */
        double x_norm = amplitude / tc->harmonic_coeffs[0]; /* normalize */
        if (fabs(x_norm) <= 1.0) {
            double T_k = chebyshev_T((int)k, x_norm);
            /* T_k isolates the k-th harmonic contribution */
            amplitude_k *= fabs(T_k) / (double)k;
        }

        harmonic_power_sum += amplitude_k * amplitude_k;
    }

    return sqrt(harmonic_power_sum) / fundamental_amplitude * 100.0;
}

/* =========================================================================
 * L1 — String Names for Enumerations
 * ========================================================================= */

const char *adc_arch_name(adc_architecture_t arch)
{
    static const char *names[] = {
        "Flash",
        "SAR",
        "Pipeline",
        "Sigma-Delta",
        "Dual-Slope",
        "Folding-Interpolating",
        "Time-Interleaved",
        "Subranging"
    };
    if ((uint32_t)arch >= ADC_ARCH_COUNT) return "Unknown";
    return names[arch];
}

const char *dac_arch_name(dac_architecture_t arch)
{
    static const char *names[] = {
        "Binary-Weighted",
        "R-2R Ladder",
        "String",
        "Current-Steering",
        "Sigma-Delta",
        "Charge-Scaling",
        "Hybrid"
    };
    if ((uint32_t)arch >= DAC_ARCH_COUNT) return "Unknown";
    return names[arch];
}

const char *adc_coding_name(adc_coding_t coding)
{
    static const char *names[] = {
        "Offset Binary",
        "Two's Complement",
        "Gray",
        "Thermometer"
    };
    if ((uint32_t)coding >= ADC_CODING_COUNT) return "Unknown";
    return names[coding];
}
