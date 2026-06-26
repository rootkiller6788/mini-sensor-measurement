/**
 * sigcond_excitation.c - Sensor Excitation Source Implementation
 *
 * Implements constant voltage/current source design, Howland current pump,
 * AC excitation for capacitive/inductive sensors, ratiometric measurement,
 * and bridge excitation optimization.
 *
 * L2: excitation stability, load regulation
 * L3: Howland current pump analysis, voltage divider models
 * L5: bridge excitation optimization, noise budget
 */

#include "sigcond_excitation.h"
#include <stdlib.h>

/* Load regulation: V_load = V_exc * R_load / (R_load + R_source)
 * Error fraction = (V_exc - V_load) / V_exc = R_source / (R_load + R_source) */
double excitation_load_error(double v_exc, double r_source_ohm, double r_load_ohm)
{
    if (r_load_ohm <= 0.0) return 1.0;
    return r_source_ohm / (r_source_ohm + r_load_ohm);
}

double excitation_max_source_impedance(double r_load_min_ohm, double accuracy_pct)
{
    if (r_load_min_ohm <= 0.0 || accuracy_pct <= 0.0) return 0.0;
    return r_load_min_ohm * accuracy_pct / 100.0;
}

double ratiometric_divider_output(double v_exc, double r_sensor, double r_ref)
{
    if (r_sensor + r_ref <= 0.0) return 0.0;
    return v_exc * r_sensor / (r_sensor + r_ref);
}

double ratiometric_sensor_resistance(double v_exc, double v_out, double r_ref)
{
    if (fabs(v_exc - v_out) < 1e-15) return INFINITY;
    if (v_out >= v_exc) return INFINITY;
    return r_ref * v_out / (v_exc - v_out);
}

/* Improved Howland current pump:
 * I_load = Vin / R_set
 * Compliance: V_load_max = V_supply_pos - V_sat - I_load * R_set
 *            V_load_min = V_supply_neg + V_sat + I_load * R_set
 */
double howland_current_source(double i_target_a, double *r_set_ohm,
                               double *vin_required_v, double v_supply_v,
                               double v_sat_v, double r_load_max_ohm)
{
    if (!r_set_ohm || !vin_required_v || i_target_a <= 0.0) return -1.0;

    /* Choose R_set such that Vin is in a reasonable range (~1-5V) */
    *r_set_ohm = 2.5 / i_target_a;
    if (*r_set_ohm < 100.0) *r_set_ohm = 100.0;
    if (*r_set_ohm > 1.0e6) *r_set_ohm = 1.0e6;

    *vin_required_v = i_target_a * (*r_set_ohm);

    /* Check compliance: V_load_max = V_supply - V_sat - I*R_set */
    double v_compliance = v_supply_v - v_sat_v - i_target_a * (*r_set_ohm);
    if (v_compliance < i_target_a * r_load_max_ohm) return -1.0;

    return v_compliance;
}

double howland_max_load(double v_supply_v, double v_sat_v, double i_target_a)
{
    if (i_target_a <= 0.0) return 0.0;
    return (v_supply_v - v_sat_v) / i_target_a;
}

double howland_output_impedance(double r_set, double cmrr_linear,
                                 double aol_open_loop_gain)
{
    if (cmrr_linear <= 0.0 || aol_open_loop_gain <= 0.0) return INFINITY;
    /* Z_out ~ R_set * CMRR_linear * AOL (first-order model) */
    return r_set * cmrr_linear * aol_open_loop_gain;
}

double ac_excitation_frequency_capacitive(double capacitance_f,
                                           double r_min_ohm, double r_max_ohm)
{
    if (capacitance_f <= 0.0 || r_min_ohm <= 0.0 || r_max_ohm <= r_min_ohm)
        return 0.0;

    /* Choose f such that Xc = sqrt(R_min * R_max) (geometric mean)
     * for best sensitivity */
    double r_opt = sqrt(r_min_ohm * r_max_ohm);
    return 1.0 / (2.0 * M_PI * capacitance_f * r_opt);
}

double capacitive_reactance(double frequency_hz, double capacitance_f)
{
    if (frequency_hz <= 0.0 || capacitance_f <= 0.0) return INFINITY;
    return 1.0 / (2.0 * M_PI * frequency_hz * capacitance_f);
}

double inductive_reactance(double frequency_hz, double inductance_h)
{
    if (inductance_h < 0.0) return 0.0;
    return 2.0 * M_PI * frequency_hz * inductance_h;
}

double ac_amplitude_rms_from_pp(double vpp)
{
    return vpp / (2.0 * sqrt(2.0));
}

/* Bridge excitation optimization:
 * Sensitivity: V_signal = V_exc * GF * strain / (4*active_factor)
 * Self-heating per gauge: P = V_exc^2 / (4 * R_gauge)
 * Optimal V_exc balances SNR (higher V_exc) vs self-heating (lower V_exc)
 *
 * For practical strain gauges:
 *   V_exc = sqrt(4 * R_gauge * P_max_per_gauge)
 *   where P_max depends on substrate thermal conductivity
 */
double bridge_optimal_excitation(double r_gauge_ohm,
                                  double max_power_mw_per_gauge,
                                  double *resulting_power_mw)
{
    if (r_gauge_ohm <= 0.0) return 0.0;

    double max_power = max_power_mw_per_gauge * 1e-3; /* mW to W */
    double v_exc = sqrt(4.0 * r_gauge_ohm * max_power);

    if (resulting_power_mw)
        *resulting_power_mw = v_exc * v_exc / (4.0 * r_gauge_ohm) * 1000.0;

    return v_exc;
}

double rtd_self_heating_error(double excitation_current_a,
                               double r0_ohm, double theta_c_per_w)
{
    double power = excitation_current_a * excitation_current_a * r0_ohm;
    return power * theta_c_per_w;
}

double excitation_noise_budget(double adc_fs_v, unsigned adc_bits,
                                double gain_vv, double target_lsb_fraction)
{
    double lsb_v = adc_fs_v / pow(2.0, (double)adc_bits);
    double max_noise_rti = lsb_v * target_lsb_fraction / gain_vv;
    return max_noise_rti;
}

double bridge_excitation_for_sensitivity(double target_v_per_unit,
                                          double gauge_factor,
                                          double strain_max)
{
    if (fabs(gauge_factor * strain_max) < 1e-15) return INFINITY;
    return target_v_per_unit * 4.0 / (gauge_factor * strain_max);
}
