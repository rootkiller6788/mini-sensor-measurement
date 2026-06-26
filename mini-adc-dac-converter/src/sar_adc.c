/**
 * @file sar_adc.c
 * @brief SAR ADC: CDAC model, binary search, redundant SAR,
 *        kT/C noise, comparator metastability, capacitor mismatch.
 *
 * Each function implements an independent knowledge point from
 * SAR ADC design (McCreary & Gray 1975, Kuttner 2002, Harpe 2018).
 */

#include "sar_adc.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define K_BOLTZMANN 1.380649e-23

/* =========================================================================
 * L2 — CDAC Array Management
 * ========================================================================= */

int cdac_init(uint32_t num_bits, double c_unit, cdac_array_t *cdac)
{
    if (!cdac || num_bits == 0 || num_bits > 24 || c_unit <= 0.0)
        return -1;

    cdac->num_bits = num_bits;
    cdac->c_unit_fF = c_unit;

    /* Allocate capacitors: one for each bit + terminal capacitor */
    cdac->capacitors_fF = (double *)calloc(num_bits + 1, sizeof(double));
    if (!cdac->capacitors_fF) return -1;

    /* Binary-weighted array: C_k = 2^k * C_unit
     * Terminal capacitor: C_term = C_unit (makes C_total = 2^{N+1} * C_unit) */
    cdac->c_total_fF = 0.0;
    for (uint32_t k = 0; k < num_bits; k++) {
        double cap_k = c_unit * pow(2.0, (double)k);
        cdac->capacitors_fF[k] = cap_k;
        cdac->c_total_fF += cap_k;
    }
    /* Terminal capacitor = C_unit */
    cdac->capacitors_fF[num_bits] = c_unit;
    cdac->c_total_fF += c_unit;

    cdac->switch_state = 0;
    cdac->top_plate_voltage = 0.0;
    return 0;
}

void cdac_free(cdac_array_t *cdac)
{
    if (!cdac) return;
    free(cdac->capacitors_fF);
    cdac->capacitors_fF = NULL;
    cdac->c_total_fF = 0.0;
}

double cdac_output_voltage(const cdac_array_t *cdac, uint64_t code)
{
    if (!cdac || cdac->c_total_fF <= 0.0) return NAN;

    /* V_out = V_ref * Σ(C_i * b_i) / C_total
     *
     * For binary-weighted array: code bits b_{N-1}...b_0 correspond
     * to capacitors C_{N-1}...C_0.
     *
     * LSB is bit 0, MSB is bit N-1. */

    double weighted_sum = 0.0;
    for (uint32_t k = 0; k < cdac->num_bits; k++) {
        if (code & (1ULL << k)) {
            weighted_sum += cdac->capacitors_fF[k];
        }
    }

    return weighted_sum / cdac->c_total_fF; /* Fraction of V_ref */
}

void cdac_sample(cdac_array_t *cdac, double v_in, double v_cm)
{
    if (!cdac) return;

    /* During sampling:
     * - Top plate is connected to V_cm (or V_in for top-plate sampling)
     * - Bottom plates are connected to V_in (input tracking)
     *
     * After sampling, the top plate is floating and holds the charge:
     *   Q_total = Σ C_i * (V_in - V_cm)
     *
     * When bottom plates switch to V_ref or GND during conversion:
     *   V_top = V_cm - V_in + V_ref * Σ(b_i * C_i) / C_total
     *
     * For the simplified model, we just compute the sampled voltage
     * and store it as the top plate's reference. */

    cdac->top_plate_voltage = v_in - v_cm; /* Relative to conversion result */
    cdac->switch_state = 0; /* Reset switches */
}

void cdac_set_code(cdac_array_t *cdac, uint64_t code)
{
    if (!cdac) return;
    cdac->switch_state = code;
}

/* =========================================================================
 * L5 — SAR Binary Search Algorithm
 * ========================================================================= */

void sar_bit_trial(sar_adc_state_t *state, cdac_array_t *cdac,
                    double v_in, double comparator_offset)
{
    if (!state || !cdac) return;

    uint32_t bit = state->current_bit;
    if (bit >= cdac->num_bits) return;

    /* Set the trial bit to '1' */
    uint64_t trial_code = state->dac_code | (1ULL << bit);

    /* Compute resulting DAC voltage */
    double dac_frac = cdac_output_voltage(cdac, trial_code);
    double v_ref = 1.0; /* Normalized V_ref = 1.0 for internal computation */
    double v_dac = dac_frac * v_ref;

    /* Compare: is V_in ≥ V_dac + offset? */
    double v_diff = v_in - v_dac - comparator_offset;
    state->comparator_input = v_diff;

    if (v_diff >= 0.0) {
        /* Keep the bit: V_in ≥ V_dac */
        state->dac_code = trial_code;
        state->comparator_result = 1;
    } else {
        /* Clear the bit: V_in < V_dac */
        state->comparator_result = -1;
        /* dac_code unchanged (bit stays 0) */
    }

    /* Move to next bit (less significant) */
    if (bit > 0) {
        state->current_bit = bit - 1;
    } else {
        state->current_bit = 0; /* This was the last bit */
    }
    state->cycles_remaining--;

    /* Update DAC voltage for the new code */
    double dac_frac_final = cdac_output_voltage(cdac, state->dac_code);
    state->dac_voltage = dac_frac_final;
}

void sar_full_conversion(sar_adc_state_t *state, cdac_array_t *cdac,
                          double v_in, const sar_adc_spec_t *spec,
                          uint64_t *code)
{
    if (!state || !cdac || !spec || !code) return;

    /* Initialize SAR state */
    uint32_t N = spec->resolution_bits;
    state->current_bit = N - 1; /* Start with MSB */
    state->dac_code = 0;
    state->dac_voltage = 0.0;
    state->cycles_remaining = N;
    state->sample_voltage = v_in;

    /* Run N bit trials */
    for (uint32_t i = 0; i < N; i++) {
        sar_bit_trial(state, cdac, v_in, spec->comparator_offset_mV * 1e-3);

        /* Model decision delay */
        state->cycles_remaining--;
        if (state->current_bit == 0 && i < N - 1) break; /* Completion */
    }

    *code = state->dac_code;
}

/* =========================================================================
 * L5 — Sub-Radix-2 (Redundant) SAR
 * ========================================================================= */

void sar_subradix_weights(uint32_t resolution_bits, double radix,
                           double *weights, uint32_t *total_steps)
{
    if (!weights || !total_steps || resolution_bits == 0) return;

    /* For radix r < 2, weights are w_k = r^k.
     * The total range must cover 2^N LSB:
     *   Σ_{k=0}^{P-1} r^k = (r^P - 1) / (r - 1) ≥ 2^N
     * So: P ≥ log_r(2^N * (r-1) + 1)
     *
     * Example: N=10, r=1.85 → P=15 steps to cover 1024 codes. */

    double target_range = pow(2.0, (double)resolution_bits);
    uint32_t P = (uint32_t)ceil(log(target_range * (radix - 1.0) + 1.0) /
                                 log(radix));
    if (P < resolution_bits) P = resolution_bits;

    /* Compute weights */
    double weight = 1.0;
    for (uint32_t k = 0; k < P; k++) {
        weights[k] = weight;
        weight *= radix;
    }

    /* Normalize so total sum = 2^N */
    double total = 0.0;
    for (uint32_t k = 0; k < P; k++) total += weights[k];
    double scale = target_range / total;
    for (uint32_t k = 0; k < P; k++) weights[k] *= scale;

    *total_steps = P;
}

void sar_redundant_conversion(sar_adc_state_t *state, cdac_array_t *cdac,
                               double v_in,
                               const double *weights, uint32_t n_steps,
                               uint64_t *code)
{
    if (!state || !cdac || !weights || !code || n_steps == 0) return;

    /* Redundant SAR conversion.
     *
     * Instead of strict binary search, we traverse the weight array
     * from largest to smallest. Each step:
     *   - Add weight[k] to DAC
     *   - If DAC > V_in: subtract weight[k] (undo the step)
     *
     * Because total weight sum > 2^N (redundancy), a "wrong" decision
     * at step k can be corrected by later steps.
     *
     * This is the basis for self-correcting SAR ADCs. */

    double dac_value = 0.0;
    uint64_t raw_code = 0;

    for (uint32_t k = 0; k < n_steps; k++) {
        dac_value += weights[k];
        raw_code |= (1ULL << k); /* Record that this weight was tried */

        /* Compare */
        if (dac_value > v_in) {
            /* Overshoot: remove this weight */
            dac_value -= weights[k];
            raw_code &= ~(1ULL << k);
        }
    }

    /* Convert sub-radix-2 code to standard binary by scaling */
    double total_range = 0.0;
    for (uint32_t k = 0; k < n_steps; k++) total_range += weights[k];
    double target_range = pow(2.0, (double)cdac->num_bits);

    /* Compute standard code: code_std = dac_value * 2^N / total_range */
    double code_std = dac_value * target_range / total_range;
    if (code_std < 0.0) code_std = 0.0;
    uint64_t max_code = (1ULL << cdac->num_bits) - 1;
    if (code_std > (double)max_code) code_std = (double)max_code;

    *code = (uint64_t)llround(code_std);
}

/* =========================================================================
 * L4 — kT/C Noise Limit
 * ========================================================================= */

double ktc_noise_rms(double temp_kelvin, double capacitance_F)
{
    if (temp_kelvin <= 0.0 || capacitance_F <= 0.0) return NAN;

    /* kT/C thermal noise on a capacitor.
     *
     * v_n_rms = sqrt(k * T / C)
     *
     * This is the fundamental noise limit of any sampled-data system.
     * It arises from the thermal noise of the switch resistance
     * being sampled onto the capacitor.
     *
     * For a 12-bit ADC with V_ref = 3.3 V and noise budget of 0.5 LSB:
     *   V_LSB = 3.3 / 4096 = 806 μV
     *   v_n_budget = 403 μV
     *   C_min = 1.38e-23 * 300 / (403e-6)² ≈ 25 fF
     *
     * In practice, margins require C >> C_min. */

    return sqrt(K_BOLTZMANN * temp_kelvin / capacitance_F);
}

double ktc_min_capacitance(double temp_kelvin, double noise_target_vrms)
{
    if (temp_kelvin <= 0.0 || noise_target_vrms <= 0.0) return NAN;

    /* C_min = k * T / (v_noise_rms)² */
    return K_BOLTZMANN * temp_kelvin / (noise_target_vrms * noise_target_vrms);
}

/* =========================================================================
 * L5 — Comparator Model
 * ========================================================================= */

int comparator_decision(double v_diff, double offset,
                         double noise_rms, double noise_val)
{
    /* Comparator with offset and Gaussian noise.
     *
     * Effective input: v_eff = v_diff + V_os + v_noise
     * where v_noise ~ N(0, σ²) and V_os is systematic offset.
     *
     * Decision: sign(v_eff) */

    double v_eff = v_diff + offset + noise_rms * noise_val;
    return (v_eff > 0.0) ? 1 : -1;
}

double comparator_metastability_prob(double decision_time_ns,
                                      double tau_ns, double v_diff)
{
    /* Comparator metastability probability.
     *
     * A regenerative latch comparator has exponential growth:
     * V_out(t) = V_in * exp(t/τ)
     *
     * If V_out doesn't reach logic threshold V_th within decision
     * time t_d, the output is metastable (unresolved).
     *
     * P(meta) = P(V_in * exp(t_d/τ) < V_th)
     *         = P(V_in < V_th * exp(-t_d/τ))
     *
     * For V_in uniformly distributed in [-V_range, V_range]:
     * P(meta) = V_th * exp(-t_d/τ) / V_range
     *
     * With V_th ≈ 1 V, V_range ≈ 1 V (normalized):
     * P(meta) ≈ exp(-t_d / τ)   for v_diff ≈ 0.
     *
     * For finite v_diff:
     * P(meta) = exp(-t_d/τ) * exp(-|v_diff|*exp(t_d/τ)/V_range) */

    if (decision_time_ns <= 0.0 || tau_ns <= 0.0) return 1.0;

    double p_base = exp(-decision_time_ns / tau_ns);

    /* For small v_diff, probability dominated by base rate.
     * For larger v_diff, probability drops exponentially. */
    double v_range = 1.0; /* Normalized comparator input range */
    double correction = exp(-fabs(v_diff) * exp(decision_time_ns / tau_ns) / v_range);

    return p_base * correction;
}

/* =========================================================================
 * L5 — Capacitor Mismatch Model (Pelgrom's Law)
 * ========================================================================= */

double capacitor_with_mismatch(double c_nominal, double area_um2,
                                double pelgrom_coeff, double random_val)
{
    if (c_nominal <= 0.0 || area_um2 <= 0.0) return NAN;

    /* Pelgrom's law for capacitor matching:
     *
     * σ(ΔC/C) = A_C / sqrt(W * L) = A_C / sqrt(Area)
     *
     * where A_C is the matching coefficient [%·μm].
     *
     * Typical values:
     *   A_C ≈ 0.5 to 2.0%·μm for MIM capacitors
     *   A_C ≈ 0.1 to 0.5%·μm for MOM capacitors
     *
     * C_actual = C_nominal + ΔC
     * ΔC = C_nominal * σ_C * random_val
     * σ_C = pelgrom_coeff / sqrt(Area) */

    double sigma = (pelgrom_coeff / 100.0) / sqrt(area_um2);
    double delta = c_nominal * sigma * random_val;
    double c_actual = c_nominal + delta;
    if (c_actual < 0.0) c_actual = 0.0; /* Physical constraint */
    return c_actual;
}

double sar_dnl_from_mismatch(uint32_t n_bits, double mismatch_pct)
{
    if (n_bits == 0 || n_bits > 24 || mismatch_pct <= 0.0) return NAN;

    /* DNL from binary-weighted capacitor mismatch.
     *
     * For binary-weighted CDAC, the DNL at the MSB transition is
     * approximately:
     *
     *   DNL_max ≈ sqrt( Σ(σ_i² * 2^i) ) / C_unit * 2^{N-1} / 2^N
     *           = sqrt( Σ(2^i) ) * σ_C * 2^{N-1} / 2^N
     *           ≈ sqrt(2^N - 1) * σ_C * 2^{N-1} / 2^N
     *           ≈ 2^{(N-1)/2} * σ_C * 0.5  (for large N)
     *
     * where σ_C = mismatch_pct / 100 * C_unit / sqrt(Area).
     * Simplified for unit capacitor with sqrt(Area) = 1: σ_C = mismatch_pct/100. */

    double sigma_c = mismatch_pct / 100.0;
    double dnl = sqrt(pow(2.0, (double)n_bits) - 1.0) * sigma_c *
                 pow(2.0, (double)(n_bits - 1)) / pow(2.0, (double)n_bits);

    return dnl;
}

/* =========================================================================
 * L6 — SAR ADC Full Simulation with Non-Idealities
 * ========================================================================= */

/* Simple pseudo-random number generator (xorshift64) */
static uint64_t xorshift64(uint64_t *state)
{
    uint64_t x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;
    return x;
}

/* Generate Gaussian random number using Box-Muller */
static double rand_gaussian(uint64_t *state)
{
    double u1 = (double)xorshift64(state) / (double)UINT64_MAX;
    double u2 = (double)xorshift64(state) / (double)UINT64_MAX;
    if (u1 < 1e-15) u1 = 1e-15;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

void sar_simulate_conversions(const sar_adc_spec_t *spec,
                               const double *v_in_array,
                               uint32_t n_conversions,
                               uint64_t *codes,
                               uint64_t *noise_seed)
{
    if (!spec || !v_in_array || !codes || n_conversions == 0) return;

    uint64_t seed = noise_seed ? *noise_seed : 12345ULL;
    uint32_t N = spec->resolution_bits;
    uint64_t max_code = (1ULL << N) - 1;

    /* Create CDAC for simulation */
    cdac_array_t cdac;
    if (cdac_init(N, spec->c_unit_fF, &cdac) != 0) return;

    /* Compute kT/C noise */
    double ktc_rms = ktc_noise_rms(300.0, cdac.c_total_fF * 1e-15);

    for (uint32_t conv = 0; conv < n_conversions; conv++) {
        double v_in = v_in_array[conv];

        /* Add kT/C noise */
        double ktc_noise_val = ktc_rms * rand_gaussian(&seed);
        double v_sampled = v_in + ktc_noise_val;

        /* Add comparator noise */
        double comp_noise = spec->comparator_noise_uVrms * 1e-6 *
                            rand_gaussian(&seed);

        /* Run binary search */
        sar_adc_state_t state;
        memset(&state, 0, sizeof(state));
        state.current_bit = N - 1;
        state.cycles_remaining = N;

        for (uint32_t bit = 0; bit < N; bit++) {
            uint32_t bit_idx = N - 1 - bit;

            /* Set trial bit */
            uint64_t trial_code = state.dac_code | (1ULL << bit_idx);
            double dac_frac = cdac_output_voltage(&cdac, trial_code);
            double dac_v = dac_frac; /* Normalized (V_ref = 1) */

            /* Comparison with offset and noise */
            double v_diff = v_sampled - dac_v - spec->comparator_offset_mV * 1e-3;
            double v_eff = v_diff + comp_noise * rand_gaussian(&seed);

            if (v_eff >= 0.0) {
                state.dac_code = trial_code;
            }
            /* else: bit stays 0 */

            state.current_bit = (bit_idx > 0) ? bit_idx - 1 : 0;
        }

        codes[conv] = state.dac_code;
        if (codes[conv] > max_code) codes[conv] = max_code;
    }

    cdac_free(&cdac);
    if (noise_seed) *noise_seed = seed;
}

/* =========================================================================
 * L5 — SAR Calibration
 * ========================================================================= */

void sar_calibrate_weights(const double *measured_weights,
                            uint32_t n_bits, double *corrections)
{
    if (!measured_weights || !corrections || n_bits == 0) return;

    /* Compute ideal binary weights: w_ideal[k] = 2^k.
     * Calibration factor = w_ideal[k] / w_measured[k]. */

    for (uint32_t k = 0; k < n_bits; k++) {
        double w_ideal = pow(2.0, (double)k);
        if (fabs(measured_weights[k]) > 1e-15) {
            corrections[k] = w_ideal / measured_weights[k];
        } else {
            corrections[k] = 1.0;
        }
    }
}

uint64_t sar_apply_calibration(uint64_t raw_code,
                                const double *corrections,
                                uint32_t n_bits)
{
    if (!corrections || n_bits == 0 || n_bits > 32) return raw_code;

    /* Apply per-bit calibration:
     * code_cal = Σ(b_k * w_ideal_k * correction_k)
     *          = Σ(b_k * 2^k * correction_k)
     *
     * Then scale to [0, 2^N - 1]. */

    double calibrated_value = 0.0;
    double total_range = 0.0;

    for (uint32_t k = 0; k < n_bits; k++) {
        if (raw_code & (1ULL << k)) {
            double w_calibrated = pow(2.0, (double)k) * corrections[k];
            calibrated_value += w_calibrated;
        }
        total_range += pow(2.0, (double)k) * corrections[k];
    }

    /* Scale to standard code range */
    double max_code_d = (double)((1ULL << n_bits) - 1);
    double code_d = calibrated_value * max_code_d / total_range;
    if (code_d < 0.0) code_d = 0.0;
    if (code_d > max_code_d) code_d = max_code_d;

    return (uint64_t)llround(code_d);
}
