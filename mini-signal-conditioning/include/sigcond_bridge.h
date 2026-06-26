/**
 * sigcond_bridge.h - Wheatstone Bridge Completion and Balance
 *
 * Bridge transfer functions, completion networks, balance,
 * lead wire compensation, sensitivity optimization.
 *
 * L2-L4: Wheatstone bridge, lead compensation, bridge nonlinearity
 * Reference: Window & Holister (1982), Hoffmann (HBM, 1974)
 */
#ifndef SIGCOND_BRIDGE_H
#define SIGCOND_BRIDGE_H
#include "sigcond_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

/** Exact bridge output: V_out = V_exc * (R3/(R2+R3) - R4/(R1+R4)) */
double bridge_output_exact(double v_exc, double r1, double r2, double r3, double r4);

/** Linearized bridge output: V_out = V_exc * dR/(4R) for quarter bridge */
double bridge_output_linear(double v_exc, double r_nominal, double delta_r, unsigned active_arms);

/** Bridge nonlinearity error: NL = (V_linear - V_exact)/V_exact * 100 */
double bridge_nonlinearity_error(double r_nominal, double delta_r, unsigned active_arms);

/** Correct bridge nonlinearity for full bridge with known GF */
double bridge_nonlinearity_correct(double v_measured, double v_exc, double gauge_factor);

/** Design bridge completion network for quarter/half bridge */
void bridge_completion_design(bridge_topology_t topology, double r_sensor_nominal, bridge_completion_t *completion);

/** Error due to completion resistor tolerance mismatch */
double bridge_completion_error(double r_sensor_nominal, double r_completion_actual);

/** Check bridge balance: |R1*R3 - R2*R4| / R_nom^2 < tolerance */
bool bridge_is_balanced(double r1, double r2, double r3, double r4, double tolerance_ratio);

/** Balance potentiometer value for given max imbalance */
double bridge_balance_pot_value(double r_nominal, double max_imbalance_pct);

/** Software zero-offset subtraction: V_corrected = V_measured - V_zero */
double bridge_software_zero(double v_measured, double v_zero_offset);

/** Lead wire resistance: R = rho*L/A. Copper rho=1.68e-8 Ohm*m at 20C. */
double bridge_lead_resistance(double length_m, double awg, double temperature_c);

/** 3-wire RTD compensation (assumes equal lead resistances) */
double bridge_3wire_compensate(double v_exc, double v_sense_high, double v_sense_low, double r_ref, double r_lead_estimated);

/** 4-wire Kelvin measurement: R_sensor = V_sensed / I_excitation */
double bridge_4wire_kelvin(double i_excitation_a, double v_sensed_v);

/** Lead wire temperature error (uncompensated, copper TCR ~3900ppm/C) */
double bridge_lead_temp_error(double r_lead_20c, double temp_c, double r_sensor_nominal);

/** Bridge sensitivity for strain gauge: S = V_exc * GF / (4*active_arms_factor) */
double bridge_sensitivity_strain(double v_exc, double gauge_factor, unsigned active_arms);

/** Load cell sensitivity: V_out = V_exc * rated_mv_v * load/fs */
double bridge_sensitivity_loadcell(double v_exc, double rated_output_mv_per_v, double full_scale_load);

/** Minimum detectable strain from noise: eps_min = noise_RTI * 4 / (V_exc * GF * gain) */
double bridge_min_detectable_strain(double v_exc, double gauge_factor, double noise_rti_uv, double gain);

/** Bridge output for given strain */
double bridge_output_strain(double v_exc, double gauge_factor, double strain, unsigned active_arms);

/** Strain from bridge output (linear approximation) */
double bridge_strain_from_output(double v_out, double v_exc, double gauge_factor, unsigned active_arms);

#ifdef __cplusplus
}
#endif
#endif /* SIGCOND_BRIDGE_H */
