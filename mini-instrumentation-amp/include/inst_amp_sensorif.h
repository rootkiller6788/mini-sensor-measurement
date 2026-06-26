/**
 * inst_amp_sensorif.h - Sensor Interface Implementations
 * L6/L7: Strain gauge, thermocouple, RTD, load cell, ECG/EEG/EMG, current shunt
 */
#ifndef INST_AMP_SENSORIF_H
#define INST_AMP_SENSORIF_H
#include "inst_amp_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Strain gauge */
typedef struct { double gauge_factor; double max_strain_ue; double v_excitation; double v_adc_fullscale; double bridge_r_ohm; } strain_gauge_config_t;
double strain_gauge_bridge_output_uv(const strain_gauge_config_t *cfg, double strain_ue);
double strain_gauge_required_gain(const strain_gauge_config_t *cfg);
double strain_gauge_adc_to_ue(double v_adc_v, double gain, const strain_gauge_config_t *cfg);

/* Thermocouple */
typedef struct { thermocouple_type_t type; double temp_ref_c; double temp_meas_c; } thermocouple_config_t;
double thermocouple_seebeck_uv_per_c(thermocouple_type_t type);
double thermocouple_voltage_uv(const thermocouple_config_t *cfg);
double thermocouple_temp_from_voltage(double v_tc_uv, thermocouple_type_t type, double temp_ref_c);
double thermocouple_cjc_compensated(double v_out_amplified, double gain, thermocouple_type_t type, double temp_ref_c);

/* RTD */
typedef struct { double r0_ohm; double i_excitation_ma; double alpha; } rtd_config_t;
double rtd_resistance_at_temp(const rtd_config_t *cfg, double temp_c);
double rtd_temp_from_resistance(const rtd_config_t *cfg, double r_ohm);
double rtd_voltage_at_temp(const rtd_config_t *cfg, double temp_c);

/* Load cell */
typedef struct { double rated_output_mv_per_v; double rated_capacity_kg; double v_excitation; double v_zero_offset; } load_cell_config_t;
double load_cell_fullscale_voltage_mv(const load_cell_config_t *cfg);
double load_cell_required_gain(const load_cell_config_t *cfg, double v_adc_fullscale);
double load_cell_vout_to_kg(double v_out, const load_cell_config_t *cfg, double actual_gain);

/* ECG */
typedef struct { double gain; double cmrr_db; double highpass_hz; double lowpass_hz; double electrode_offset_mv; } ecg_frontend_config_t;
double ecg_required_gain(double ecg_amplitude_mv, double adc_fullscale_v);
double ecg_cm_noise_rti_uv(double v_cm_interference_v, double cmrr_db);
double ecg_highpass_attenuation(double f_hz, double f_cutoff_hz);

/* EEG */
typedef struct { double gain; double cmrr_db; double bw_low_hz; double bw_high_hz; } eeg_frontend_config_t;
double eeg_noise_floor_rti_uv(double en_rti_nv_rhz, double bw_hz);

/* EMG */
typedef struct { double gain; double bw_low_hz; double bw_high_hz; double cmrr_db; } emg_frontend_config_t;
double emg_signal_bandwidth_rms(const emg_frontend_config_t *cfg);

/* Current shunt */
typedef struct { double r_shunt_ohm; double i_max_a; double v_cm_v; bool high_side; } current_shunt_config_t;
double current_shunt_vout_at_current(const current_shunt_config_t *cfg, double i_load_a);
double current_shunt_required_gain(const current_shunt_config_t *cfg, double v_adc_fullscale);
double current_shunt_power_loss(const current_shunt_config_t *cfg, double i_a);

/* Hall sensor */
typedef struct { double sensitivity_mv_per_mt; double quiescent_output_v; double bias_current_ma; } hall_sensor_config_t;
double hall_field_from_voltage(double v_out, double v_quiescent, double sensitivity_mv_per_mt);

/* Ratiometric */
double ratiometric_adc_code(double v_bridge, double v_exc, int adc_bits);
double ratiometric_measurand_from_code(int adc_code, int adc_bits, double fullscale_value);

/* Sample buffer */
typedef struct { double *buffer; int capacity; int count; } sensor_sample_buffer_t;
void sensor_buffer_init(sensor_sample_buffer_t *buf, double *storage, int cap);
void sensor_buffer_add(sensor_sample_buffer_t *buf, double sample);
double sensor_buffer_mean(const sensor_sample_buffer_t *buf);
void sensor_buffer_clear(sensor_sample_buffer_t *buf);

/* Mains notch */
double mains_notch_filter(double current_sample, double prev_average, int period_samples, int *sample_index);

#ifdef __cplusplus
}
#endif
#endif