/**
 * inst_amp_analysis.h - Noise and Error Budget Analysis
 * L3/L4: Johnson noise, shot noise, 1/f noise, error budget, CMRR analysis
 */
#ifndef INST_AMP_ANALYSIS_H
#define INST_AMP_ANALYSIS_H
#include "inst_amp_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Johnson-Nyquist thermal noise */
double noise_thermal_nv_per_rhz(double r_ohm, double temp_c);
double noise_thermal_rms_uv(double r_ohm, double bw_hz, double temp_c);

/* Shot noise */
double noise_shot_fa_per_rhz(double i_dc_a);

/* 1/f noise */
double noise_one_over_f_nv_rms(double k_f_nv2, double f_low_hz, double f_high_hz);
double noise_corner_frequency(double k_f_nv2, double en_thermal_nv2_per_hz);

/* Total noise */
double noise_total_rti_nv_per_rhz(double en_amp_nv_rhz, double in_amp_pa_rhz, double rs_ohm, double rg_ohm, double r1_ohm, double gain, double temp_c);

/* Noise figure */
double noise_figure_db(double en_amp_total_nv_rhz, double rs_ohm, double temp_c);
double noise_figure_cascade_db(const double *nf_db, const double *gain_linear, int n_stages);

/* Error budget */
double error_budget_rti_uv(const ia_spec_t *spec, double rs_ohm, double v_cm_v, double v_ripple_v);
double gain_error_budget_pct(double rg_tol_pct, double G, double A_OL);

/* SNR and dynamic range */
double signal_to_noise_ratio_db(double v_signal_rms, double v_noise_rms);
double minimum_detectable_signal(double v_noise_rms);
double dynamic_range_db(double v_max, double v_noise_floor);

/* CMRR analysis */
double cmrr_from_impedance_imbalance(double zin_cm_ohm, double delta_z_ohm);
double cmrr_combined_db(double cmrr_amp_db, double cmrr_res_db, double cmrr_impedance_db);

/* Sensitivity analysis */
double gain_sensitivity_to_rg(double G);
double gain_tempco_ppm_per_c(double G, double r1_ppm, double rg_ppm);
double offset_drift_uv(double vos_25c_uv, double tc_vos_nv_per_c, double temp_c);
double bias_current_offset_uv(double ib_pa, double rs_ohm);
double total_unadjusted_error_mv(const ia_spec_t *spec, double G, double v_rti_error_uv);
double optimal_source_resistance(double en_nv_rhz, double in_pa_rhz);
double rg_noise_contribution_rti_nv_rhz(double rg_ohm, double r1_ohm, double G, double temp_c);

#ifdef __cplusplus
}
#endif
#endif