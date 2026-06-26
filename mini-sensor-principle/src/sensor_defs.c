/**
 * @file sensor_defs.c
 * @brief Core sensor data structure operations, error budget computation
 *
 * Level: L1 — Definitions implementation
 */

#include "sensor_defs.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ──── L1: Sensor Specification Validation ──── */

/**
 * @brief Validate sensor specifications for physical consistency.
 *
 * Checks:
 *   1. Sensitivity > 0 (or < 0, but not 0)
 *   2. Full-scale range positive
 *   3. Resolution ≤ full scale
 *   4. Offset ≤ full scale (within reason)
 *   5. Bandwidth positive
 *   6. Temperature range physically meaningful
 *
 * Returns bitmask of violations:
 *   bit 0: zero sensitivity
 *   bit 1: invalid full scale
 *   bit 2: resolution exceeds FS
 *   bit 3: excessive offset
 *   bit 4: invalid bandwidth
 *   bit 5: invalid temperature range
 *   0 = all valid
 */
int sensor_specs_validate(const sensor_specs_t *s)
{
    int violations = 0;

    if (fabs(s->sensitivity) < 1e-30) {
        violations |= 0x01;
    }
    if (s->full_scale_input <= 0.0 || s->full_scale_output <= 0.0) {
        violations |= 0x02;
    }
    if (s->resolution > s->full_scale_input) {
        violations |= 0x04;
    }
    if (fabs(s->offset_error) > 2.0 * s->full_scale_output) {
        violations |= 0x08;
    }
    if (s->bandwidth <= 0.0) {
        violations |= 0x10;
    }
    if (s->temp_range_min >= s->temp_range_max ||
        s->temp_range_min < -273.15 || s->temp_range_max > 3000.0) {
        violations |= 0x20;
    }

    return violations;
}

/**
 * @brief Compute total static error (root-sum-square of all error sources).
 *
 * Per GUM: Total error cannot be simply summed because error sources
 * are independent. The RSS combination gives the most probable total error.
 *
 * ε_total = √(ε_nonlin² + ε_hysteresis² + ε_repeat² + ε_offset² + ε_gain²)
 *
 * This is expressed as % full scale.
 */
double sensor_total_error_percent(const sensor_specs_t *s)
{
    double e_nonlin  = s->nonlinearity;
    double e_hyst    = s->hysteresis;
    double e_repeat  = s->repeatability;
    double e_offset  = fabs(s->offset_error) / s->full_scale_output * 100.0;
    double e_gain    = s->gain_error;

    return sqrt(e_nonlin*e_nonlin + e_hyst*e_hyst + e_repeat*e_repeat
                + e_offset*e_offset + e_gain*e_gain);
}

/**
 * @brief Compute the sensor's dynamic error for a given input frequency.
 *
 * For 1st-order systems:
 *   Magnitude error = |1 - 1/√(1 + (f/f_c)²)| · 100%
 *
 * This represents the attenuation error when measuring a sinusoid
 * without compensating for the sensor's frequency response.
 *
 * @param f_input input signal frequency [Hz]
 * @param f_bandwidth sensor -3 dB bandwidth [Hz]
 * @return magnitude error [%]
 */
double sensor_dynamic_error_percent(double f_input, double f_bandwidth)
{
    if (f_bandwidth <= 0.0) return 100.0;
    double ratio = f_input / f_bandwidth;
    double magnitude = 1.0 / sqrt(1.0 + ratio * ratio);
    return fabs(1.0 - magnitude) * 100.0;
}

/**
 * @brief Simulate sensor response including all error sources.
 *
 * Applies to an ideal input value x_true:
 *   x_raw  = x_true + offset + noise
 *   x_gain = nonlinearity_error(x_raw) + hysteresis_error(x_raw)
 *   x_out  = x_raw · (1 + gain_error/100) + x_gain
 *
 * The noise is assumed Gaussian with given standard deviation.
 * Nonlinearity and hysteresis are modeled as worst-case ± values
 * applied with sinusoidal dependence on input (realistic smooth shape).
 *
 * @param x_true   true measurand value
 * @param s        sensor specifications
 * @param s_noise  noise standard deviation (0 = no noise)
 * @param ascending 1 = approaching from lower values, 0 = from higher
 * @return simulated sensor output including errors
 */
double sensor_simulate_reading(double x_true, const sensor_specs_t *s,
                                double s_noise, int ascending)
{
    /* Offset */
    double x = x_true + s->offset_error;

    /* Add Gaussian noise (Box-Muller approximation using uniform) */
    /* Simple approximation: use ±3σ rectangular distribution */
    if (s_noise > 0.0) {
        /* Generate pseudo-random from input value (deterministic simulation) */
        double pseudo_rand = fmod(fabs(x_true * 7919.0 + 104729.0), 1.0);
        x += s_noise * (pseudo_rand * 6.0 - 3.0);
    }

    /* Nonlinearity: model as cubic shape NL(x) = ε_nl · [(2x/FS)³ - (2x/FS)] */
    double x_norm = 2.0 * x / s->full_scale_input - 1.0;  /* [-1, 1] */
    double nl_shape = x_norm * x_norm * x_norm - x_norm;
    double nl_offset = (s->nonlinearity / 100.0) * s->full_scale_output
                       * nl_shape;

    /* Hysteresis: depends on direction */
    double hyst_offset = 0.0;
    if (ascending) {
        hyst_offset = -(s->hysteresis / 200.0) * s->full_scale_output
                      * (1.0 - x_norm * x_norm);
    } else {
        hyst_offset = (s->hysteresis / 200.0) * s->full_scale_output
                      * (1.0 - x_norm * x_norm);
    }

    /* Gain error */
    double x_with_gain = x * (1.0 + s->gain_error / 100.0);

    /* Combine: offset + gain + nonlinearity + hysteresis */
    double y = x_with_gain + nl_offset + hyst_offset;

    /* Clamp to plausible range */
    if (y > s->full_scale_output * 1.2) y = s->full_scale_output * 1.2;
    if (y < -s->full_scale_output * 0.2) y = -s->full_scale_output * 0.2;

    return y;
}

/* ──── L1: Uncertainty Budget Computation per GUM ──── */

/**
 * @brief Compute combined standard uncertainty from budget components.
 *
 * u_c = √(Σ u_i²)
 *
 * Each u_i represents one standard uncertainty component.
 * Type A: statistical (repeatability)
 * Type B: systematic (reference, resolution, thermal, loading, noise)
 *
 * The expanded uncertainty U = k · u_c uses:
 *   k = 2 for ~95% confidence (normal distribution, large ν_eff)
 *   k = t_(ν_eff, 0.95) otherwise per Student's t-distribution
 *
 * Effective degrees of freedom (Welch-Satterthwaite):
 *   ν_eff = u_c⁴ / Σ(u_i⁴/ν_i)
 */
int sensor_uncertainty_budget_compute(sensor_uncertainty_budget_t *b)
{
    if (!b) return -1;

    /* Sum all Type B components in quadrature */
    double sum_B2 = b->type_b_reference * b->type_b_reference
                  + b->type_b_resolution * b->type_b_resolution
                  + b->type_b_thermal * b->type_b_thermal
                  + b->type_b_loading * b->type_b_loading
                  + b->type_b_noise * b->type_b_noise
                  + b->type_b_calibration * b->type_b_calibration;

    /* Type A (statistical) */
    double sum_A2 = b->type_a_repeatability * b->type_a_repeatability;

    /* Combined standard uncertainty */
    b->combined_standard = sqrt(sum_A2 + sum_B2);

    /* Expanded uncertainty with k=2 */
    b->expanded_uncertainty = 2.0 * b->combined_standard;

    /* Effective degrees of freedom (simplified Welch-Satterthwaite):
     * Type A: ν_A = n - 1 (for n repeated readings, assume n=10 → ν=9)
     * Type B: ν_B = ∞ (infinite, as they come from datasheet limits)
     * With dominant Type B components, ν_eff ≈ ∞ → use k=2 */
    size_t n_repeat = b->degrees_of_freedom + 1;
    if (n_repeat <= 1) n_repeat = 2;
    double nu_A = (double)(n_repeat - 1);

    /* If Type A dominates, adjust; otherwise default to large ν_eff */
    if (sum_A2 > sum_B2 && nu_A > 1) {
        b->degrees_of_freedom = (size_t)nu_A;
    } else {
        b->degrees_of_freedom = 1000000;  /* effectively ∞, k=2 valid */
    }

    return 0;
}

/* ──── L1: Sensor Identity Utilities ──── */

/**
 * @brief Check if two sensor identities describe the same sensor model.
 *
 * Compares manufacturer, model, and domain (not serial number).
 */
int sensor_identity_same_model(const sensor_identity_t *a,
                                const sensor_identity_t *b)
{
    if (!a || !b) return 0;
    return (strncmp(a->manufacturer, b->manufacturer, 31) == 0
            && strncmp(a->model, b->model, 31) == 0
            && a->domain == b->domain);
}

/**
 * @brief Estimate power consumption of a passive sensor with given excitation.
 *
 * For passive sensors (RTD, strain gauge, thermistor):
 *   P = V_ex² / R_sensor  or  P = I_ex² · R_sensor
 *
 * For active sensors: P = V_supply · I_supply (from datasheet)
 */
double sensor_estimate_power_passive(double excitation, double resistance,
                                      int is_current_excitation)
{
    if (resistance <= 0.0) return 0.0;
    if (is_current_excitation) {
        return excitation * excitation * resistance;
    } else {
        return excitation * excitation / resistance;
    }
}

/**
 * @brief Compute thermal noise-limited resolution.
 *
 * The fundamental limit on resolution is set by Johnson noise:
 *   σ_V = √(4·k_B·T·R·Δf)
 *   Resolution = σ_V / sensitivity
 *
 * At room temperature (300 K), for R=1 kΩ, Δf=1 Hz:
 *   σ_V ≈ 4.07 nV → with sensitivity 1 V/unit: resolution = 4.07e-9 units
 */
double sensor_thermal_noise_limit(double temperature_K, double resistance,
                                   double bandwidth, double sensitivity)
{
    /* k_B = 1.380649e-23 J/K */
    const double kB = 1.380649e-23;
    double vn_rms = sqrt(4.0 * kB * temperature_K * resistance * bandwidth);
    if (fabs(sensitivity) < 1e-30) return INFINITY;
    return vn_rms / fabs(sensitivity);
}

/* ──── L1: Reading Validity Checks ──── */

/**
 * @brief Check if sensor reading is within expected range.
 *
 * Returns:
 *   0 = valid
 *   1 = out of range (high)
 *  -1 = out of range (low)
 *   2 = saturated
 *  -2 = invalid flag set
 */
int sensor_reading_check_valid(const sensor_reading_t *r,
                                double x_min, double x_max)
{
    if (!r) return -2;
    if (!r->valid) return -2;
    if (r->saturation) return 2;
    if (r->value > x_max) return 1;
    if (r->value < x_min) return -1;
    return 0;
}

/**
 * @brief Compute measurement rate from sample indices.
 *
 * Measurement rate = (idx_1 - idx_0) / (t_1 - t_0)
 *
 * Useful for detecting dropped samples or timing irregularities.
 */
double sensor_measurement_rate(const sensor_reading_t *r0,
                                const sensor_reading_t *r1)
{
    if (!r0 || !r1) return 0.0;
    uint64_t dt_us = r1->timestamp_us - r0->timestamp_us;
    if (dt_us == 0) return 0.0;
    double dt_s = (double)dt_us * 1e-6;
    double d_idx = (double)(r1->sample_index - r0->sample_index);
    return d_idx / dt_s;
}

/* ──── L1: Noise Type Classification ──── */

/**
 * @brief Classify the dominant noise type based on PSD frequency dependence.
 *
 * If the noise PSD S(f) ∝ f^α, then:
 *   α = 0  → white (thermal, shot)
 *   α = -1 → flicker / 1/f
 *   α = -2 → random walk (Brownian)
 *   α = +1 → blue noise
 *   α = +2 → violet noise
 *
 * This function estimates α from two PSD measurements at different frequencies.
 * Assumes S(f) ≈ K · f^α over the frequency range.
 *
 * @param S1 PSD at frequency f1 [V²/Hz]
 * @param f1 first frequency [Hz]
 * @param S2 PSD at frequency f2 [Hz]
 * @param f2 second frequency [Hz]
 * @return estimated exponent α
 */
double sensor_noise_exponent_estimate(double S1, double f1,
                                       double S2, double f2)
{
    if (S1 <= 0.0 || S2 <= 0.0 || f1 <= 0.0 || f2 <= 0.0) return 0.0;
    if (f1 == f2) return 0.0;
    return log(S2 / S1) / log(f2 / f1);
}

/**
 * @brief Compute noise bandwidth from single-pole filter.
 *
 * For a 1st-order low-pass filter with time constant τ:
 *   Noise equivalent bandwidth: B_n = 1/(4τ) = π/2 · f_c
 *
 * This is the bandwidth of an ideal brick-wall filter that passes
 * the same total noise power as the real filter.
 */
double sensor_noise_equivalent_bandwidth_1st_order(double tau)
{
    if (tau <= 0.0) return INFINITY;
    return 1.0 / (4.0 * tau);
}

/**
 * @brief Noise factor of a sensor system relative to thermal noise.
 *
 * F = total_output_noise_power / (G · k·T₀·B)
 *
 * Where G is the system gain, T₀ = 290 K.
 * F = 1 (0 dB) for an ideal noiseless system.
 * F > 1 means the system adds noise beyond the thermal floor.
 *
 * @param total_noise_Vrms total RMS output noise [V]
 * @param gain system voltage gain [V/V]
 * @param bandwidth noise bandwidth [Hz]
 * @param source_resistance source resistance [Ω] (use 0 if not applicable)
 * @return noise factor [linear, ≥ 1]
 */
double sensor_noise_factor(double total_noise_Vrms, double gain,
                            double bandwidth, double source_resistance)
{
    const double kB = 1.380649e-23;
    const double T0 = 290.0;  /* standard noise temperature */

    /* Input-referred noise */
    double input_noise_Vrms = total_noise_Vrms / gain;

    /* Thermal noise of source resistance */
    double thermal_Vrms = sqrt(4.0 * kB * T0 * source_resistance * bandwidth);

    /* Noise factor: F = (total input noise)^2 / (thermal noise)^2 */
    double total_input_noise2 = input_noise_Vrms * input_noise_Vrms;
    double thermal_noise2 = thermal_Vrms * thermal_Vrms;

    if (thermal_noise2 < 1e-30) {
        /* No source resistance specified, cannot compute noise factor.
         * Return total noise relative to 1 nV/√Hz floor */
        return total_noise_Vrms / (gain * 1e-9 * sqrt(bandwidth));
    }

    return total_input_noise2 / thermal_noise2;
}
