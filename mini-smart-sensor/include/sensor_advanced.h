#ifndef SENSOR_ADVANCED_H
#define SENSOR_ADVANCED_H
#include "smart_sensor.h"
#ifdef __cplusplus
extern "C" {
#endif

/* L8: Advanced Topic APIs */

/* Energy Harvesting */
double ss_advanced_pv_power_uw(double cell_efficiency, double cell_area_cm2, double irradiance_uw_cm2);
double ss_advanced_teg_power(double seebeck_coeff, double delta_temp_k, double internal_r);
double ss_advanced_avg_power_uw(double p_active_mw, double p_sleep_uw, double duty_cycle, double wakeup_energy_uj, double wakeup_rate_hz);
double ss_advanced_max_duty_cycle(double p_harvest_uw, double p_active_mw, double p_sleep_uw);

/* Allan Variance */
int ss_advanced_allan_variance(const double *data, size_t n_points, double sample_rate_hz, const double *tau_values, size_t n_tau, double *allan_dev_out);
int ss_advanced_allan_fit(const double *tau_values, const double *allan_dev, size_t n_tau, ss_allan_variance_t *avar_out);

/* Edge ML Feature Extraction */
double ss_advanced_signal_rms(const double *data, size_t n);
double ss_advanced_crest_factor(const double *data, size_t n);
double ss_advanced_skewness(const double *data, size_t n);
double ss_advanced_spectral_centroid(const double *magnitudes, const double *freqs, size_t n_bins);
double ss_advanced_zero_crossing_rate(const double *data, size_t n);

#ifdef __cplusplus
}
#endif
#endif
