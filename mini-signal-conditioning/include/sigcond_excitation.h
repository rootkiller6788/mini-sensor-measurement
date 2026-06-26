/**
 * sigcond_excitation.h - Sensor Excitation Sources
 *
 * Constant voltage, constant current, AC excitation, ratiometric
 * measurement, Howland current pump, bridge excitation optimization.
 *
 * L2-L5: excitation stability, load regulation, Howland pump design
 * Reference: Pallas-Areny & Webster (2001), Kester (2004)
 */
#ifndef SIGCOND_EXCITATION_H
#define SIGCOND_EXCITATION_H
#include "sigcond_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

/** Load regulation error: V_sensor = V_exc * R_sensor/(R_sensor+R_source) */
double excitation_load_error(double v_exc, double r_source_ohm, double r_load_ohm);

/** Max source impedance for target accuracy: R_source <= R_load * accuracy%/100 */
double excitation_max_source_impedance(double r_load_min_ohm, double accuracy_pct);

/** Ratiometric divider: V_out = V_exc * R_sensor / (R_sensor + R_ref) */
double ratiometric_divider_output(double v_exc, double r_sensor, double r_ref);

/** Sensor resistance from ratiometric measurement: R_sensor = R_ref * V_out / (V_exc - V_out) */
double ratiometric_sensor_resistance(double v_exc, double v_out, double r_ref);

/** Improved Howland current pump design: I_load = Vin / R_set */
double howland_current_source(double i_target_a, double *r_set_ohm, double *vin_required_v, double v_supply_v, double v_sat_v, double r_load_max_ohm);

/** Howland max load resistance within compliance */
double howland_max_load(double v_supply_v, double v_sat_v, double i_target_a);

/** Howland output impedance limited by op-amp CMRR and open-loop gain */
double howland_output_impedance(double r_set, double cmrr_linear, double aol_open_loop_gain);

/** Optimal AC excitation freq for capacitive sensor: Xc between R_min and R_max */
double ac_excitation_frequency_capacitive(double capacitance_f, double r_min_ohm, double r_max_ohm);

/** Capacitive reactance: Xc = 1/(2*pi*f*C) */
double capacitive_reactance(double frequency_hz, double capacitance_f);

/** Inductive reactance: Xl = 2*pi*f*L */
double inductive_reactance(double frequency_hz, double inductance_h);

/** AC amplitude: Vrms = Vpp / (2*sqrt(2)) */
double ac_amplitude_rms_from_pp(double vpp);

/** Optimal bridge excitation: trades SNR vs self-heating. Returns optimal V_exc. */
double bridge_optimal_excitation(double r_gauge_ohm, double max_power_mw_per_gauge, double *resulting_power_mw);

/** RTD self-heating temperature error: delta_T = P * theta */
double rtd_self_heating_error(double excitation_current_a, double r0_ohm, double theta_c_per_w);

/** Excitation noise must be << 1 LSB: V_noise_exc < V_LSB / gain */
double excitation_noise_budget(double adc_fs_v, unsigned adc_bits, double gain_vv, double target_lsb_fraction);

/** Bridge excitation for target sensitivity: V_exc = (target_V/unit) * 4 / GF */
double bridge_excitation_for_sensitivity(double target_v_per_unit, double gauge_factor, double strain_max);

#ifdef __cplusplus
}
#endif
#endif /* SIGCOND_EXCITATION_H */
