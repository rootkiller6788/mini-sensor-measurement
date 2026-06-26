/**
 * inst_amp_calibration.h - Calibration and Error Compensation
 * L5: Two-point, multi-point, polynomial calibration, auto-zero, chopper
 */
#ifndef INST_AMP_CALIBRATION_H
#define INST_AMP_CALIBRATION_H
#include "inst_amp_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Two-point calibration */
void ia_cal_two_point(double v_cal1, double v_out1, double v_cal2, double v_out2, double *actual_gain, double *output_offset);
double ia_cal_apply_correction(double v_out_raw, double actual_gain, double output_offset);

/* Multi-point calibration */
void ia_cal_linear_regression(const double *v_cal, const double *v_out, int n_points, double *gain, double *offset);
void ia_cal_polynomial_fit(const double *v_out, const double *v_in, int n_points, int degree, double *coeffs);
double ia_cal_eval_polynomial(const double *coeffs, int degree, double x);

/* Temperature compensation */
double ia_temp_compensate_offset(double vos_at_t0_uv, double tc1_nv_per_c, double tc2_nv_per_c2, double t_c, double t0_c);

/* Auto-zero and chopper */
double ia_autozero_measure_offset(double v_out_zero_input, double gain);
double ia_chopper_ripple_amplitude(double vos_input, double gain, double lpf_attenuation_db);
double ia_chopper_effective_corner_freq(double f_corner_original, double f_chop, double f_lpf);

/* Digital trim */
double ia_digital_trim_factor(double g_actual, double g_target);

/* Calibration utilities */
bool ia_cal_is_valid(const ia_calibration_t *cal);
double ia_ratiometric_correction(double reading_raw, double v_exc_nominal, double v_exc_actual);
double ia_cal_average_readings(const double *readings, int n);
double ia_cal_std_deviation(const double *readings, int n, double mean);

/* Self-test */
typedef struct { bool output_in_range; bool gain_plausible; bool not_oscillating; double measured_gain; double measured_offset; } ia_self_test_result_t;
void ia_self_test(double v_test_input, double v_out_measured, double expected_gain, ia_self_test_result_t *result);

#ifdef __cplusplus
}
#endif
#endif