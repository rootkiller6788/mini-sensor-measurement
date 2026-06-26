/**
 * sigcond_selfcal.h - Self-Calibration and Auto-Zero
 *
 * L5: Two-point calibration, auto-zero, gain error compensation,
 * temperature drift compensation.
 */
#ifndef SIGCOND_SELFCAL_H
#define SIGCOND_SELFCAL_H
#include "sigcond_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

int selfcal_two_point(const calibration_point_t *p1, const calibration_point_t *p2, double *gain, double *offset);
double selfcal_apply_correction(double measured_output, double gain, double offset);
double selfcal_predict_output(double input_value, double gain, double offset);
double selfcal_autozero_measure(double measured_zero_output);
double selfcal_autozero_correct(double measured_output, double zero_offset);
double selfcal_offset_drift_estimate(double last_zero_offset, double drift_rate_v_per_sec, double time_since_zero_sec);
double selfcal_gain_calibrate(double v_input_known, double v_output_measured, double gain_nominal);
double selfcal_gain_correction_factor(double gain_actual, double gain_nominal);
double selfcal_gain_drift(double gain_at_t0, double gain_drift_ppm_per_c, double delta_temp_c);
void selfcal_multichannel_calibrate(unsigned num_channels, const sigcond_channel_config_t channels[], const double cal_p1_inputs[], const double cal_p1_outputs[], const double cal_p2_inputs[], const double cal_p2_outputs[], double gains[], double offsets[]);
double selfcal_temp_compensate_gain(double gain_at_t0, double alpha_gain_ppm_per_c, double delta_temp_c);
double selfcal_temp_compensate_offset(double offset_at_t0, double alpha_offset_uv_per_c, double beta_offset_uv_per_c2, double delta_temp_c);
bool selfcal_is_calibration_valid(double cal_temp_c, double current_temp_c, double cal_time_hours, double current_time_hours, double max_temp_drift_c, double max_time_hours);
double selfcal_total_uncertainty(double reference_uncertainty, double measurement_noise, double drift_uncertainty, double linearity_error);

#ifdef __cplusplus
}
#endif
#endif
