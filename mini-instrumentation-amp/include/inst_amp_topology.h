/**
 * inst_amp_topology.h - IA Topology Function Declarations
 * See inst_amp_defs.h for type definitions.
 */
#ifndef INST_AMP_TOPOLOGY_H
#define INST_AMP_TOPOLOGY_H
#include "inst_amp_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ==================================================================
 * L2: Topology Configuration Structures
 * ================================================================== */

typedef struct {
    double r1_ohm;
    double r2_ohm;
    double r3_ohm;
    double rg_ohm;
    double resistor_tol_pct;
} ia_three_opamp_config_t;

typedef struct {
    double r1_ohm;
    double r2_ohm;
    double r3_ohm;
    double r4_ohm;
    double rg_ohm;
} ia_two_opamp_config_t;

typedef struct {
    double rg_ohm;
    double rl_ohm;
    double gm_ms;
} ia_current_mode_config_t;

typedef struct {
    double cin_pf;
    double cfb_pf;
    double fsw_khz;
    double charge_injection_fc;
} ia_flying_cap_config_t;

typedef struct {
    double gm_input_ms;
    double rf_ohm;
} ia_indirect_current_config_t;

/* ==================================================================
 * L2: Function Declarations
 * ================================================================== */

/* 3-op-amp */
double ia_three_opamp_gain(const ia_three_opamp_config_t *cfg);
double ia_three_opamp_cmrr_resistor_limit(const ia_three_opamp_config_t *cfg);
double ia_three_opamp_rg_for_gain(const ia_three_opamp_config_t *cfg, double target_gain);
double ia_three_opamp_output(const ia_three_opamp_config_t *cfg, double v_in_plus, double v_in_minus, double v_ref, double cmrr_db);
bool ia_output_in_range(double v_out, double v_ee, double v_cc, double headroom_mv);

/* 2-op-amp */
double ia_two_opamp_gain(const ia_two_opamp_config_t *cfg);
double ia_two_opamp_cmrr_zero_freq(const ia_two_opamp_config_t *cfg, double gbp_mhz);
void ia_two_opamp_input_impedance(const ia_two_opamp_config_t *cfg, double *zin_plus, double *zin_minus);

/* Current-mode */
double ia_current_mode_gain(const ia_current_mode_config_t *cfg);
double ia_current_mode_bandwidth(const ia_current_mode_config_t *cfg, double c_comp_pf);
double ia_current_mode_output_max(const ia_current_mode_config_t *cfg, double v_cc, double v_sat);

/* Flying-capacitor */
double ia_flying_cap_gain(const ia_flying_cap_config_t *cfg);
double ia_flying_cap_cmrr_ideal(void);
double ia_flying_cap_charge_injection_error(const ia_flying_cap_config_t *cfg);

/* Indirect current feedback */
double ia_indirect_current_gain(const ia_indirect_current_config_t *cfg);
double ia_indirect_current_bandwidth(const ia_indirect_current_config_t *cfg, double c_parasitic_pf);

/* Topology comparison */
double ia_topology_cmrr_at_freq(ia_topology_t topo, double cmrr_dc_db, double f_hz, double gbp_mhz, double gain);
double ia_topology_bandwidth(ia_topology_t topo, double gain, double gbp_mhz, double c_comp_pf);

/* Auto-range */
double ia_autorange_gain(double v_in_rms, double v_fs_adc, double g_min, double g_max, int num_ranges);
bool ia_autorange_should_increase_gain(double v_current, double v_max_for_gain, double hysteresis);
bool ia_autorange_should_decrease_gain(double v_current, double v_max_for_gain, double hysteresis);

/* L2: Topology recommendation */
ia_topology_t ia_recommend_topology(sensor_type_t sensor, double supply_voltage);
double ia_total_input_noise(double en_amp_nv_rhz, double in_amp_pa_rhz, double rs_ohm, double bw_hz, double temp_c);

/* L1: Spec and state management */
void ia_spec_init(ia_spec_t *spec);
void ia_state_init(ia_state_t *state);
ia_error_code_t ia_spec_validate(const ia_spec_t *spec);

/* L2: Core utility functions */
double ia_gain_error_percent(double g_nominal, double g_actual);
double ia_effective_bits(double v_signal_rms, double v_noise_rms);
double ia_cmrr_at_frequency(double cmrr_dc_db, double f_hz, double f_pole_hz);
double ia_psrr_output_error(double v_supply_ripple, double psrr_db);
double ia_gain_for_sensor_span(double v_sensor_fs, double v_adc_fs);
double ia_e96_nearest(double r_ohm);
double ia_e24_nearest(double r_ohm);

/* L2: Calibration helpers */
void ia_calibrate_two_point(double v_cal1, double v_out1, double v_cal2, double v_out2, double *actual_gain, double *input_offset);
double ia_apply_calibration(double v_out_raw, double actual_gain, double nominal_gain, double input_offset);

/* L3: Bridge sensor */
double bridge_output_voltage(const bridge_sensor_t *sensor, bool exact);
double bridge_nonlinearity_pct(const bridge_sensor_t *sensor);

#ifdef __cplusplus
}
#endif
#endif