#ifndef SENSOR_APPLICATIONS_H
#define SENSOR_APPLICATIONS_H
#include "smart_sensor.h"
#ifdef __cplusplus
extern "C" {
#endif

/* L7: Application APIs */

/* Industrial IoT Vibration */
double ss_app_vibration_severity(double accel_rms_g, double dominant_freq_hz);
int ss_app_vibration_zone(double velocity_rms_mm_s);
double ss_app_bpfo_frequency(int n_balls, double rpm, double ball_diameter, double pitch_diameter, double contact_angle);
double ss_app_bpfi_frequency(int n_balls, double rpm, double ball_diameter, double pitch_diameter, double contact_angle);

/* Automotive TPMS */
double ss_app_tpms_compensate_pressure(double pressure_kpa, double temp_meas_c, double temp_ref_c);
double ss_app_tpms_deviation_percent(double measured_pressure_kpa, double placard_pressure_kpa);
int ss_app_tpms_detect_slow_leak(const double *pressure_history, const double *time_hours, size_t n_samples, double *leak_rate_out);

/* Medical Pulse Oximeter */
double ss_app_spo2_estimate(double ac_red, double dc_red, double ac_ir, double dc_ir, double cal_a, double cal_b);
double ss_app_heart_rate_from_ppg(const double *peak_intervals_ms, size_t n_intervals);
double ss_app_perfusion_index(double ac, double dc);

/* Environmental Air Quality */
double ss_app_aqi_pm25(double pm25_ug_m3);
int ss_app_co2_comfort_level(double co2_ppm);
double ss_app_dew_point(double temp_c, double rel_humidity);
double ss_app_heat_index(double temp_c, double rel_humidity);

#ifdef __cplusplus
}
#endif
#endif
