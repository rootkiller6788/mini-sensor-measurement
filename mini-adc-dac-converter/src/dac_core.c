/**
 * @file dac_core.c
 * @brief DAC core implementations: R-2R ladder, binary-weighted,
 *        string DAC, glitch energy, settling behavior.
 *
 * Knowledge Coverage:
 *   L1 - DAC specifications and architectures
 *   L2 - R-2R ladder operation principle, string DAC monotonicity
 *   L3 - Transfer function models for each DAC topology
 *   L5 - Glitch energy computation, settling time modeling
 *   L6 - R-2R ladder linearity analysis
 *
 * References:
 *   - Razavi, "Principles of Data Conversion System Design" (1995), Ch. 3-5
 *   - Kester, "Data Conversion Handbook" (2005), Ch. 3
 *   - IEEE Std 1658-2011: Terminology and Test Methods for DAC
 */

#include "adc_dac_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =========================================================================
 * L2 — R-2R Ladder DAC
 * ========================================================================= */

/**
 * R-2R ladder DAC voltage computation.
 *
 * Principle:
 * An R-2R ladder is a resistive network with only two resistor values:
 * R (series) and 2R (shunt). For N bits, the ladder has N stages.
 *
 * Each bit b_k (k = 0 for LSB, k = N-1 for MSB) controls a switch
 * that connects a 2R resistor to either V_ref (b_k = 1) or GND (b_k = 0).
 *
 * By Thévenin's theorem, the output voltage is:
 *   V_out = V_ref * (b_{N-1}/2 + b_{N-2}/4 + ... + b_0/2^N)
 *         = V_ref * code / 2^N
 *
 * R-2R advantage: only two resistor values needed regardless of N.
 * R-2R disadvantage: resistor matching limits linearity (~10-12 bits).
 *
 * @param code         Digital input code.
 * @param n_bits       Resolution.
 * @param v_ref        Reference voltage.
 * @return Output voltage.
 */
double r2r_dac_voltage(uint64_t code, uint32_t n_bits, double v_ref)
{
    if (n_bits == 0 || n_bits > 32 || v_ref <= 0.0) {
        return NAN;
    }
    uint64_t max_code = (1ULL << n_bits) - 1;
    if (code > max_code) code = max_code;
    return v_ref * (double)code / (double)(1ULL << n_bits);
}

/**
 * R-2R ladder output impedance.
 *
 * The output impedance looking back into the R-2R ladder
 * at the output node equals R (for any N, assuming ideal switches).
 *
 * Proof: Each stage has R to ground (via one 2R), plus 2R in
 * parallel with the Thevenin equivalent of the rest, recursively
 * yielding R_th = R.
 */
double r2r_output_impedance(double r_ohms)
{
    return r_ohms;
}

/**
 * Compute the contribution of each bit in an R-2R ladder.
 *
 * The MSB (bit N-1) contributes V_ref/2.
 * Bit k contributes V_ref / 2^{N-k}.
 *
 * @param bit_index  Bit index (0 = LSB, N-1 = MSB).
 * @param n_bits     Total resolution.
 * @param v_ref      Reference voltage.
 * @return Voltage contribution of this bit.
 */
double r2r_bit_weight(uint32_t bit_index, uint32_t n_bits, double v_ref)
{
    if (bit_index >= n_bits || n_bits == 0 || n_bits > 32) return NAN;
    if (v_ref <= 0.0) return NAN;
    uint32_t bit_position_from_msb = n_bits - 1 - bit_index;
    return v_ref / pow(2.0, (double)(bit_position_from_msb + 1));
}

/* =========================================================================
 * L2 — Binary-Weighted Resistor DAC
 * ========================================================================= */

/**
 * Binary-weighted resistor DAC voltage.
 *
 * Uses resistors: R_k = R / 2^k for bit k.
 * All resistors connected to a common summing junction (virtual ground).
 *
 * I_out = V_ref * Σ(b_k / R_k) = (V_ref/R) * Σ(b_k * 2^k)
 * V_out = -I_out * R_fb = -V_ref * (R_fb/R) * code
 *
 * For R_fb = R: V_out = -V_ref * code / 2^N (but inverted).
 *
 * Issues: Wide resistor spread (R_max/R_min = 2^{N-1}:1).
 * This limits practical resolution to ~8 bits.
 *
 * @param code    Digital code.
 * @param n_bits  Resolution.
 * @param v_ref   Reference voltage.
 * @param r_fb    Feedback resistor value.
 * @param r_base  Base resistor (R for MSB).
 * @return Output voltage (inverting).
 */
double binary_weighted_dac_voltage(uint64_t code, uint32_t n_bits,
                                    double v_ref, double r_fb, double r_base)
{
    if (n_bits == 0 || n_bits > 32 || v_ref <= 0.0 ||
        r_base <= 0.0 || r_fb <= 0.0) {
        return NAN;
    }
    uint64_t max_code = (1ULL << n_bits) - 1;
    if (code > max_code) code = max_code;
    /* V_out = -V_ref * (R_fb / R) * Σ(b_k * 2^k) / 2^N */
    return -v_ref * (r_fb / r_base) * (double)code / (double)(1ULL << n_bits);
}

/**
 * Compute the maximum resistor ratio in a binary-weighted DAC.
 *
 * ratio = R_MSB / R_LSB = 1 : (1 / 2^{N-1}) = 2^{N-1}
 * This is the key practical limitation.
 */
double binary_weighted_resistor_ratio(uint32_t n_bits)
{
    if (n_bits < 2 || n_bits > 32) return NAN;
    return pow(2.0, (double)(n_bits - 1));
}

/* =========================================================================
 * L2 — String DAC (Kelvin Divider)
 * ========================================================================= */

/**
 * String DAC (resistor string / Kelvin divider) output voltage.
 *
 * A chain of 2^N equal resistors between V_ref and GND.
 * A tree of switches (or thermometer decoder) selects one of
 * the 2^N tap points.
 *
 * V_out(k) = k * V_ref / 2^N   for k = 0..2^N-1
 *
 * Key advantage: INHERENTLY MONOTONIC.
 *   The output voltage always increases with increasing code.
 *   DNL ≥ -1 LSB is guaranteed.
 *
 * Disadvantage: 2^N resistors and 2^N switches → area explodes.
 *
 * @param code    Digital code.
 * @param n_bits  Resolution.
 * @param v_ref   Reference voltage.
 * @return Output voltage.
 */
double string_dac_voltage(uint64_t code, uint32_t n_bits, double v_ref)
{
    if (n_bits == 0 || n_bits > 32 || v_ref <= 0.0) return NAN;
    uint64_t max_code = (1ULL << n_bits) - 1;
    if (code > max_code) code = max_code;
    /* Divide V_ref across 2^N equal resistors, tap at code-th junction. */
    return v_ref * (double)code / (double)(1ULL << n_bits);
}

/**
 * Compute resistor value for string DAC.
 *
 * Each of the 2^N resistors: R_string = V_ref / (2^N * I_bias)
 * Total resistance: R_total = 2^N * R_string = V_ref / I_bias
 *
 * Settling time: τ ≈ R_string * C_load (worst case at midpoint)
 * R_string_at_worst = (R_total/2) || (R_total/2) = R_total/4
 */
double string_dac_resistor_value(double v_ref, double i_bias, uint32_t n_bits)
{
    if (n_bits == 0 || n_bits > 32 || v_ref <= 0.0 || i_bias <= 0.0)
        return NAN;
    uint64_t n_resistors = 1ULL << n_bits;
    return v_ref / ((double)n_resistors * i_bias);
}

/* =========================================================================
 * L5 — DAC Glitch Energy
 * ========================================================================= */

/**
 * Compute glitch energy for a DAC code transition.
 *
 * Glitch energy is the time integral of the voltage error during
 * a code transition, typically caused by timing skew between switches.
 *
 * E_glitch = ∫ |V_out(t) - V_ideal(t)| dt   [V·s]
 *
 * Worst-case glitch occurs at the MSB transition:
 * code 011...1 → 100...0, where all bits change simultaneously.
 * Any switch timing mismatch produces a large transient.
 *
 * Simplified model: glitch modeled as a triangular pulse
 * with amplitude A_glitch and duration T_glitch.
 *
 * @param amp_glitch_v Peak glitch amplitude [V].
 * @param dur_glitch_s Glitch duration [s].
 * @return Glitch energy in V·s (pico: multiply by 1e12 for pV·s).
 */
double dac_glitch_energy_triangular(double amp_glitch_v, double dur_glitch_s)
{
    if (amp_glitch_v < 0.0) amp_glitch_v = fabs(amp_glitch_v);
    /* Area of triangle: 0.5 * base * height */
    return 0.5 * dur_glitch_s * amp_glitch_v;
}

/**
 * Estimate glitch energy from MSB transition timing skew.
 *
 * When code transitions from 01...1 to 10...0:
 * - MSB turns ON (delay t_on)
 * - All other bits turn OFF (delay t_off)
 * If t_on ≠ t_off, a temporary state occurs:
 *   - t_on < t_off: briefly at ~V_ref (both MSB and all others on)
 *   - t_on > t_off: briefly at ~0 (all bits off)
 *
 * E_glitch ≈ V_ref * |t_on - t_off| / (2 * R_out)
 *
 * @param v_ref      Reference voltage.
 * @param t_skew_s   Timing skew between MSB and other bits.
 * @param r_out_ohm  Output resistance.
 * @return Estimated glitch energy in V·s.
 */
double dac_glitch_from_skew(double v_ref, double t_skew_s, double r_out_ohm)
{
    if (v_ref <= 0.0 || t_skew_s < 0.0 || r_out_ohm <= 0.0) return NAN;
    return v_ref * fabs(t_skew_s) / (2.0 * r_out_ohm);
}

/* =========================================================================
 * L5 — DAC Settling Behavior
 * ========================================================================= */

/**
 * Model DAC settling as a first-order RC circuit.
 *
 * V_out(t) = V_final * (1 - exp(-t / τ))
 *
 * Settling time to within 0.5 LSB:
 *   error = V_final * exp(-t / τ) < 0.5 * V_LSB
 *   t_settle > τ * ln(V_final / (0.5 * V_LSB))
 *           ≈ τ * (N * ln(2) + ln(2))   for V_final ≈ V_ref
 *           = τ * (N+1) * ln(2)
 *
 * @param t_current    Current time since transition start [s].
 * @param tau_s        Time constant τ = R*C [s].
 * @param v_final      Final settled voltage.
 * @return Output voltage at time t_current.
 */
double dac_settling_rc(double t_current, double tau_s, double v_final)
{
    if (t_current < 0.0 || tau_s <= 0.0) return NAN;
    return v_final * (1.0 - exp(-t_current / tau_s));
}

/**
 * Compute time to settle within a specified error.
 *
 * t = -τ * ln(V_error / V_step)
 *
 * @param tau_s        Time constant.
 * @param v_step       Magnitude of voltage step.
 * @param v_error_max  Maximum allowed error.
 * @return Settling time [s].
 */
double dac_settling_time(double tau_s, double v_step, double v_error_max)
{
    if (tau_s <= 0.0 || v_step <= 0.0 || v_error_max <= 0.0) return NAN;
    return -tau_s * log(v_error_max / v_step);
}

/**
 * Compute worst-case settling time for N-bit DAC.
 *
 * Full-scale step: V_step = V_ref.
 * Error budget = 0.5 * V_LSB = 0.5 * V_ref / 2^N.
 * t_settle = τ * ln(V_ref / (0.5*V_ref/2^N)) = τ * ln(2^{N+1}) = τ*(N+1)*ln(2)
 *
 * @param tau_s   DAC time constant.
 * @param n_bits  Resolution.
 * @return Worst-case settling time.
 */
double dac_worst_case_settling_time(double tau_s, uint32_t n_bits)
{
    if (tau_s <= 0.0 || n_bits == 0 || n_bits > 32) return NAN;
    return tau_s * (double)(n_bits + 1) * log(2.0);
}

/* =========================================================================
 * L5 — DAC Transfer Function (Segmented / Hybrid)
 * ========================================================================= */

/**
 * Compute V_out for a segmented (coarse + fine) DAC.
 *
 * Hybrid architecture: coarse DAC (M bits, thermometer) + fine DAC (K bits, R-2R).
 * Total resolution N = M + K.
 *
 * V_out = V_ref * (code_coarse * 2^K + code_fine) / 2^N
 *
 * Segmentation reduces glitch energy (coarse bits change less frequently)
 * and improves DNL (thermometer coarse guarantees monotonicity).
 *
 * @param code          Full N-bit code.
 * @param m_bits        Coarse DAC bits.
 * @param k_bits        Fine DAC bits.
 * @param v_ref         Reference voltage.
 * @return Output voltage.
 */
double segmented_dac_voltage(uint64_t code, uint32_t m_bits,
                              uint32_t k_bits, double v_ref)
{
    uint32_t n_bits = m_bits + k_bits;
    if (n_bits == 0 || n_bits > 32 || v_ref <= 0.0) return NAN;
    uint64_t max_code = (1ULL << n_bits) - 1;
    if (code > max_code) code = max_code;
    return v_ref * (double)code / (double)(1ULL << n_bits);
}

/**
 * Split a code into coarse and fine components for segmented DAC.
 *
 * @param code     Full code.
 * @param m_bits   Coarse bits.
 * @param k_bits   Fine bits.
 * @param[out] coarse_code  Upper M bits.
 * @param[out] fine_code    Lower K bits.
 */
void segmented_split_code(uint64_t code, uint32_t m_bits, uint32_t k_bits,
                           uint64_t *coarse_code, uint64_t *fine_code)
{
    *fine_code = code & ((1ULL << k_bits) - 1);
    *coarse_code = (code >> k_bits) & ((1ULL << m_bits) - 1);
}

/* =========================================================================
 * L5 — Current-Steering DAC Model
 * ========================================================================= */

/**
 * Compute current-steering DAC differential output voltage.
 *
 * Uses a differential pair of current sources. For code k:
 * I_plus = I_unit * k, I_minus = I_unit * (2^N - 1 - k)
 * V_out_diff = (I_plus - I_minus) * R_load
 *            = I_unit * (2k - 2^N + 1) * R_load
 *
 * @param code      Digital code.
 * @param n_bits    Resolution.
 * @param i_unit    Unit current source [A].
 * @param r_load    Load resistor [Ω].
 * @return Differential output voltage.
 */
double current_steering_dac_voltage(uint64_t code, uint32_t n_bits,
                                     double i_unit, double r_load)
{
    if (n_bits == 0 || n_bits > 32 || i_unit <= 0.0 || r_load <= 0.0)
        return NAN;
    uint64_t max_code = (1ULL << n_bits) - 1;
    if (code > max_code) code = max_code;
    double i_diff = i_unit * (2.0 * (double)code - (double)max_code);
    return i_diff * r_load;
}

/* =========================================================================
 * L5 — DAC SFDR from Amplitude Mismatch
 * ========================================================================= */

/**
 * Estimate SFDR from current source amplitude mismatch in current-steering DAC.
 *
 * For N unit current sources with standard deviation σ_I/I:
 * The mismatch creates harmonic distortion.
 *
 * SFDR[dB] ≈ 20*log10(1 / (σ_I/I)) + 20*log10(N/2) - 3
 *
 * Simplified empirical model:
 * SFDR ≈ 20*log10(1 / mismatch_std) + 10*log10(2^N) - 6  [dB]
 *
 * @param n_bits         Resolution.
 * @param mismatch_std   Current source standard deviation [fraction].
 * @return Estimated SFDR [dB].
 */
double current_steering_sfdr_estimate(uint32_t n_bits, double mismatch_std)
{
    if (n_bits == 0 || n_bits > 32 || mismatch_std <= 0.0) return NAN;
    /* SFDR limited by mismatch: ≈ -20*log10(σ_I/I) + 6*(N-1) [dB] */
    return -20.0 * log10(mismatch_std) + 6.0 * (double)(n_bits - 1);
}

/* =========================================================================
 * L5 — DAC Integral and Differential Nonlinearity
 * ========================================================================= */

/**
 * Compute DNL for a DAC from a set of measured output voltages.
 *
 * DNL[k] = (V_out[k+1] - V_out[k]) / V_LSB_ideal - 1   for k = 0..2^N-2
 *
 * @param v_out_measured  Array of measured output voltages (length ≥ 2^N).
 * @param n_bits          Resolution.
 * @param v_ref           Reference voltage.
 * @param[out] dnl        Output DNL array (length = 2^N - 1).
 */
void dac_compute_dnl(const double *v_out_measured, uint32_t n_bits,
                      double v_ref, double *dnl)
{
    if (!v_out_measured || !dnl || n_bits == 0 || n_bits > 32 ||
        v_ref <= 0.0) return;
    uint64_t n_codes = 1ULL << n_bits;
    double v_lsb_ideal = v_ref / (double)n_codes;
    if (v_lsb_ideal <= 0.0) return;

    for (uint64_t k = 0; k < n_codes - 1; k++) {
        double step = v_out_measured[k + 1] - v_out_measured[k];
        dnl[k] = (step / v_lsb_ideal) - 1.0;
    }
}

/**
 * Compute INL for a DAC from DNL (cumulative sum).
 */
void dac_compute_inl_from_dnl(const double *dnl, uint32_t n_bits,
                               double *inl)
{
    if (!dnl || !inl || n_bits == 0 || n_bits > 32) return;
    uint64_t n_codes = 1ULL << n_bits;
    inl[0] = 0.0;
    double running_sum = 0.0;
    for (uint64_t k = 1; k < n_codes; k++) {
        running_sum += dnl[k - 1];
        inl[k] = running_sum;
    }
}

/**
 * Determine if a DAC transfer function is monotonic from DNL.
 *
 * A DAC is monotonic if for every code increase, the output voltage
 * increases (or stays equal). This requires:
 *   step[k] = V_out[k+1] - V_out[k] ≥ 0  for all k.
 * Equivalently: DNL[k] > -1 LSB.
 *
 * @param dnl      DNL array.
 * @param n_codes  Number of output codes.
 * @return 1 if monotonic, 0 otherwise.
 */
int dac_is_monotonic(const double *dnl, uint64_t n_codes)
{
    if (!dnl || n_codes < 2) return 0;
    for (uint64_t k = 0; k < n_codes - 1; k++) {
        if (dnl[k] <= -1.0) return 0; /* Missing code / non-monotonic step */
    }
    return 1;
}
