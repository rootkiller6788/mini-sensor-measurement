/**
 * inst_amp_sensorif.c - Sensor Interface and Application Implementations
 *
 * L6 Canonical Problems: Complete sensor readout implementations for
 * strain gauge bridges, thermocouples, RTDs, load cells, biopotential
 * electrodes (ECG/EEG/EMG), and current shunts.
 *
 * L7 Applications:
 *   - Medical instrumentation (ECG front-end per IEC 60601)
 *   - Industrial weighing (load cell signal conditioning)
 *   - Temperature measurement (thermocouple cold-junction compensation)
 *   - Automotive current sensing (shunt amplifier)
 *
 * Courses: Stanford EE359, Michigan EECS 455, Georgia Tech ECE 6601
 * Reference:
 *   Pall++s-Areny & Webster, "Sensors and Signal Conditioning" (2nd ed., 2001)
 *   Webster, "Medical Instrumentation: Application and Design" (4th ed., 2010)
 *   Kitchin & Counts (2006), Ch. 5, 7
 *   NIST ITS-90 Thermocouple Reference Tables
 */

#include "inst_amp_defs.h"
#include "inst_amp_topology.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ==================================================================
 * L6: Strain Gauge Bridge Readout
 *
 * Complete signal chain for a Wheatstone bridge strain gauge:
 *   Bridge -> IA -> LP filter -> ADC
 *
 * Strain calculation:
 *   epsilon = deltaL/L = (GF)^(-1) * deltaR/R
 *   where GF = gauge factor (~2 for metal foil, ~100 for semiconductor)
 *
 * Bridge output (quarter bridge, linear):
 *   V_bridge = V_exc * GF * epsilon / 4
 *
 * Required IA gain to map full-scale strain to ADC:
 *   G = V_adc_fs / (V_exc * GF * epsilon_max / 4)
 *
 * Example: GF=2, V_exc=5V, epsilon_max=1000 uStrain, V_adc_fs=5V
 *   V_bridge_fs = 5 * 2 * 0.001 / 4 = 2.5 mV
 *   G = 5 / 0.0025 = 2000
 *
 * Reference: Vishay Micro-Measurements, "Strain Gage Selection" (TN-505)
 * Complexity: O(1)
 * ================================================================== */

typedef struct {
    double gauge_factor;        /* GF, typically 2.0 for Constantan */
    double max_strain_ue;       /* Maximum strain in microstrain (ue = 1e-6) */
    double v_excitation;        /* Bridge excitation (V) */
    double v_adc_fullscale;     /* ADC full-scale input (V) */
    double bridge_r_ohm;        /* Nominal bridge resistance (Ohm) */
} strain_gauge_config_t;

double strain_gauge_bridge_output_uv(const strain_gauge_config_t *cfg,
                                      double strain_ue) {
    if (!cfg) return 0.0;
    /* V_bridge = V_exc * GF * epsilon / 4 */
    return cfg->v_excitation * cfg->gauge_factor * strain_ue * 1e-6 / 4.0 * 1e6;
}

double strain_gauge_required_gain(const strain_gauge_config_t *cfg) {
    if (!cfg) return 0.0;
    double v_fs = strain_gauge_bridge_output_uv(cfg, cfg->max_strain_ue) * 1e-6;
    if (v_fs <= 0.0) return 0.0;
    return cfg->v_adc_fullscale / v_fs;
}

/* ==================================================================
 * L6: Strain to Engineering Units Conversion
 *
 * Converts ADC reading back to microstrain:
 *   epsilon_ue = V_adc * 4 / (G * V_exc * GF) * 1e6
 *
 * Includes bridge nonlinearity correction for large strains:
 *   epsilon_true = epsilon_linear / (1 - epsilon_linear/2)
 *   (Taylor correction for quarter-bridge)
 *
 * Complexity: O(1)
 * ================================================================== */
double strain_gauge_adc_to_ue(double v_adc_v, double gain,
                               const strain_gauge_config_t *cfg) {
    if (!cfg || gain <= 0.0) return 0.0;

    double v_bridge = v_adc_v / gain;
    double epsilon_linear = v_bridge * 4.0 / (cfg->v_excitation * cfg->gauge_factor);

    /* Nonlinearity correction for quarter bridge (Taylor) */
    double epsilon_true = epsilon_linear / (1.0 - epsilon_linear / 2.0);

    return epsilon_true * 1e6;  /* Convert to microstrain */
}

/* ==================================================================
 * L6: Thermocouple Signal Conditioning
 *
 * Thermocouples produce a small differential voltage proportional
 * to the temperature difference between the hot junction (T_meas)
 * and the cold junction (T_ref).
 *
 * Seebeck effect: V_tc = S(T) * (T_meas - T_ref)
 * where S(T) is the temperature-dependent Seebeck coefficient.
 *
 * Key implementation challenges:
 *   1. Cold junction compensation (CJC): measure T_ref with RTD/thermistor
 *   2. Microvolt-level signals: need high-gain, low-noise IA
 *   3. Low source impedance (~10 Ohm): current noise less critical
 *   4. Common-mode filtering: industrial noise pickup
 *
 * Reference: NIST ITS-90 Thermocouple Database
 *   IEC 60584-1: Thermocouples - EMF specifications and tolerances
 * Complexity: O(1)
 * ================================================================== */

typedef struct {
    thermocouple_type_t type;
    double temp_ref_c;          /* Cold junction (reference) temperature (C) */
    double temp_meas_c;         /* Measured (hot junction) temperature (C) */
} thermocouple_config_t;

/**
 * Get nominal Seebeck coefficient at 25C for each thermocouple type.
 * Values from NIST ITS-90 Monograph 175.
 *
 * Complexity: O(1)
 */
double thermocouple_seebeck_uv_per_c(thermocouple_type_t type) {
    switch (type) {
        case TC_TYPE_K: return 40.6;   /* Chromel-Alumel */
        case TC_TYPE_J: return 51.9;   /* Iron-Constantan */
        case TC_TYPE_T: return 40.6;   /* Copper-Constantan */
        case TC_TYPE_E: return 60.9;   /* Chromel-Constantan */
        case TC_TYPE_N: return 26.0;   /* Nicrosil-Nisil */
        case TC_TYPE_R: return 10.4;   /* Pt13%Rh-Pt */
        case TC_TYPE_S: return 9.4;    /* Pt10%Rh-Pt */
        case TC_TYPE_B: return 5.3;    /* Pt30%Rh-Pt6%Rh */
        default: return 40.0;
    }
}

/**
 * Compute thermocouple output voltage using linear approximation.
 *
 * V_tc = S * (T_meas - T_ref)
 *
 * Note: Real thermocouple response is nonlinear. For precision work,
 * NIST polynomial coefficients (degree 8-10) should be used.
 * This linear model has ~0.5-2% error over a 100C range.
 *
 * Complexity: O(1)
 */
double thermocouple_voltage_uv(const thermocouple_config_t *cfg) {
    if (!cfg) return 0.0;

    double S = thermocouple_seebeck_uv_per_c(cfg->type);
    return S * (cfg->temp_meas_c - cfg->temp_ref_c);
}

/**
 * Compute temperature from thermocouple voltage (linear inverse).
 *
 * T_meas = T_ref + V_tc / S
 *
 * Complexity: O(1)
 */
double thermocouple_temp_from_voltage(double v_tc_uv,
                                       thermocouple_type_t type,
                                       double temp_ref_c) {
    double S = thermocouple_seebeck_uv_per_c(type);
    if (fabs(S) < 1e-9) return temp_ref_c;
    return temp_ref_c + v_tc_uv / S;
}

/**
 * Cold junction compensation using measured reference temperature.
 *
 * The thermocouple measures T_meas - T_ref. The CJC sensor (RTD,
 * thermistor, or IC sensor) provides T_ref. The IA amplifies the
 * thermocouple voltage V_tc. The digital system adds T_ref to get
 * T_meas.
 *
 * T_meas = T_ref + V_tc_amplified / (G * S)
 *
 * Complexity: O(1)
 */
double thermocouple_cjc_compensated(double v_out_amplified, double gain,
                                     thermocouple_type_t type, double temp_ref_c) {
    if (gain <= 0.0) return temp_ref_c;

    double v_tc_uv = v_out_amplified / gain * 1e6;  /* Convert V to uV */
    return thermocouple_temp_from_voltage(v_tc_uv, type, temp_ref_c);
}

/* ==================================================================
 * L6: RTD (Pt100) Signal Conditioning
 *
 * RTD resistance varies with temperature:
 *   R(T) = R_0 * [1 + A*T + B*T^2 + C*(T-100)*T^3]
 *   (Callendar-van Dusen equation, IEC 60751)
 *
 * For Pt100: R_0 = 100 Ohm at 0C
 *   A = 3.9083e-3 /C
 *   B = -5.775e-7 /C^2
 *   C = -4.183e-12 /C^4 (for T < 0C only)
 *
 * Linear approximation (acceptable for -50 to 150C):
 *   R(T) = R_0 * [1 + alpha*(T - T_0)]
 *   alpha = 0.00385 /C (DIN 43760)
 *
 * RTD readout uses a constant current source to generate
 * a voltage proportional to resistance:
 *   V_rtd = I_exc * R(T)
 *
 * The IA amplifies the difference between R(T) and a fixed
 * reference resistor.
 *
 * Reference: IEC 60751, Industrial platinum resistance thermometers
 * Complexity: O(1)
 * ================================================================== */

typedef struct {
    double r0_ohm;              /* Resistance at 0C (100 for Pt100) */
    double i_excitation_ma;     /* Excitation current (mA) */
    double alpha;               /* Temp coefficient (0.00385 for Pt) */
} rtd_config_t;

double rtd_resistance_at_temp(const rtd_config_t *cfg, double temp_c) {
    if (!cfg) return 0.0;

    /* Linear approximation R(T) = R0 * (1 + alpha * T) */
    return cfg->r0_ohm * (1.0 + cfg->alpha * temp_c);
}

double rtd_temp_from_resistance(const rtd_config_t *cfg, double r_ohm) {
    if (!cfg || cfg->alpha <= 0.0 || cfg->r0_ohm <= 0.0) return 0.0;
    return (r_ohm - cfg->r0_ohm) / (cfg->r0_ohm * cfg->alpha);
}

double rtd_voltage_at_temp(const rtd_config_t *cfg, double temp_c) {
    if (!cfg) return 0.0;
    double R = rtd_resistance_at_temp(cfg, temp_c);
    return cfg->i_excitation_ma * 1e-3 * R;  /* V = I * R */
}

/* ==================================================================
 * L6: Load Cell Signal Conditioning
 *
 * Load cells use a fully active Wheatstone bridge (4 strain gauges)
 * to measure force/weight. Key specifications:
 *   - Rated output: 2 mV/V (or 3 mV/V) at full scale
 *   - Excitation: 5-10 VDC
 *   - Bridge resistance: 350 Ohm or 1000 Ohm
 *
 * Full-scale output at 10V excitation, 2 mV/V:
 *   V_fs = 10V * 2 mV/V = 20 mV
 *
 * Required IA gain for 5V ADC range:
 *   G = 5V / 20mV = 250
 *
 * Load cell calibration uses a linear model:
 *   F = (V_out - V_zero) / sensitivity
 *
 * Reference: Vishay Revere, "Load Cell Technology" (TN-507)
 * Complexity: O(1)
 * ================================================================== */

typedef struct {
    double rated_output_mv_per_v; /* mV/V at full scale */
    double rated_capacity_kg;     /* Full-scale capacity (kg) */
    double v_excitation;          /* Excitation voltage (V) */
    double v_zero_offset;         /* Zero-load output voltage */
} load_cell_config_t;

double load_cell_fullscale_voltage_mv(const load_cell_config_t *cfg) {
    if (!cfg) return 0.0;
    return cfg->rated_output_mv_per_v * cfg->v_excitation;
}

double load_cell_required_gain(const load_cell_config_t *cfg,
                                double v_adc_fullscale) {
    double v_fs = load_cell_fullscale_voltage_mv(cfg) * 1e-3;
    if (v_fs <= 0.0) return 0.0;
    return v_adc_fullscale / v_fs;
}

double load_cell_vout_to_kg(double v_out, const load_cell_config_t *cfg,
                             double actual_gain) {
    if (!cfg || actual_gain <= 0.0) return 0.0;

    double v_bridge = v_out / actual_gain;
    double v_fs_mv = load_cell_fullscale_voltage_mv(cfg);
    if (v_fs_mv <= 0.0) return 0.0;

    double fraction = (v_bridge * 1000.0 - cfg->v_zero_offset) / v_fs_mv;
    return fraction * cfg->rated_capacity_kg;
}

/* ==================================================================
 * L6/L7: ECG/EKG Biopotential Front-End
 *
 * The electrocardiogram (ECG) measures the electrical activity of
 * the heart via surface electrodes. Typical signal characteristics:
 *   - Amplitude: 0.5 - 5 mV (R-wave peak)
 *   - Frequency range: 0.05 - 150 Hz (diagnostic)
 *   - Electrode offset: up to +/-300 mV DC
 *   - Common-mode noise: 50/60 Hz mains, up to several volts
 *
 * Key IA requirements for ECG (IEC 60601-2-25):
 *   - CMRR > 90 dB at 50/60 Hz
 *   - Input impedance > 10 MOhm
 *   - Gain: 100-1000
 *   - Noise RTI < 30 uVpp (0.05-150 Hz)
 *   - Patient leakage current < 10 uA (IEC 60601)
 *
 * Driven Right Leg (DRL) circuit:
 *   The common-mode voltage is sensed, inverted, and fed back
 *   to the patient's right leg via a resistor. This actively
 *   cancels common-mode interference, improving effective CMRR
 *   by 20-40 dB.
 *
 * Reference:
 *   Webster, "Medical Instrumentation: Application and Design" (4th ed., 2010)
 *   IEC 60601-2-25: Medical electrical equipment - ECG
 * Complexity: O(1)
 * ================================================================== */

typedef struct {
    double gain;                /* IA gain (100-1000 typical) */
    double cmrr_db;             /* IA CMRR at mains frequency */
    double highpass_hz;         /* AC coupling cutoff (0.05 Hz typical) */
    double lowpass_hz;          /* Anti-aliasing filter (150 Hz typical) */
    double electrode_offset_mv; /* Maximum electrode DC offset */
} ecg_frontend_config_t;

double ecg_required_gain(double ecg_amplitude_mv, double adc_fullscale_v) {
    if (ecg_amplitude_mv <= 0.0) return 0.0;
    return adc_fullscale_v / (ecg_amplitude_mv * 1e-3);
}

double ecg_cm_noise_rti_uv(double v_cm_interference_v, double cmrr_db) {
    double cmrr_linear = pow(10.0, cmrr_db / 20.0);
    if (cmrr_linear <= 0.0) return v_cm_interference_v * 1e6;
    return (v_cm_interference_v / cmrr_linear) * 1e6;
}

/**
 * DC blocking (highpass) filter for electrode offset removal.
 *
 * A first-order highpass filter with cutoff f_c:
 *   |H(f)| = 1 / sqrt(1 + (f_c/f)^2)
 *
 * At f = f_c: -3 dB attenuation
 * At f = 0.05 Hz with f_c = 0.05 Hz: ECG attenuated by ~3 dB at lowest freq
 * Electrode offset (DC) blocked completely.
 *
 * Complexity: O(1)
 */
double ecg_highpass_attenuation(double f_hz, double f_cutoff_hz) {
    if (f_cutoff_hz <= 0.0) return 1.0;
    return 1.0 / sqrt(1.0 + (f_cutoff_hz / f_hz) * (f_cutoff_hz / f_hz));
}

/* ==================================================================
 * L6/L7: EEG Biopotential Front-End
 *
 * Electroencephalogram (EEG) signals are 10-50x smaller than ECG
 * and require lower noise and higher gain:
 *   - Amplitude: 10-200 uV
 *   - Frequency: 0.5-50 Hz (clinical), DC-200 Hz (research)
 *   - Electrode impedance: 1-10 kOhm (wet), 10-100 kOhm (dry)
 *
 * IA requirements for EEG:
 *   - CMRR > 110 dB
 *   - Input impedance > 100 MOhm
 *   - Noise RTI < 1 uVpp (0.5-50 Hz)
 *   - Gain: 1000-10000
 *
 * Reference: Niedermeyer & Lopes da Silva,
 *   "Electroencephalography: Basic Principles" (5th ed., 2005)
 * Complexity: O(1)
 * ================================================================== */

typedef struct {
    double gain;                /* IA gain (1000-10000 typical) */
    double cmrr_db;             /* CMRR at 50/60 Hz */
    double bw_low_hz;           /* Low-frequency cutoff (0.5 Hz) */
    double bw_high_hz;          /* High-frequency cutoff (50 Hz for clinical) */
} eeg_frontend_config_t;

double eeg_noise_floor_rti_uv(double en_rti_nv_rhz, double bw_hz) {
    /* Integrated RTI noise over bandwidth */
    return en_rti_nv_rhz * sqrt(bw_hz) * 1e-3;
}

/* ==================================================================
 * L6/L7: EMG Biopotential Front-End
 *
 * Electromyogram (EMG) signals from skeletal muscle:
 *   - Amplitude: 0.1-5 mV
 *   - Frequency: 20-500 Hz (surface EMG)
 *   - Source impedance: 1-10 kOhm
 *
 * IA requirements for EMG:
 *   - Gain: 100-1000
 *   - CMRR > 90 dB
 *   - Bandwidth: 20-500 Hz
 *   - Input impedance > 100 MOhm
 *
 * Reference: Merletti & Parker,
 *   "Electromyography: Physiology, Engineering, and Noninvasive Applications" (2004)
 * Complexity: O(1)
 * ================================================================== */

typedef struct {
    double gain;                /* IA gain (100-1000 typical) */
    double bw_low_hz;           /* Highpass cutoff (20 Hz) */
    double bw_high_hz;          /* Lowpass cutoff (500 Hz) */
    double cmrr_db;             /* CMRR */
} emg_frontend_config_t;

double emg_signal_bandwidth_rms(const emg_frontend_config_t *cfg) {
    if (!cfg) return 0.0;
    return cfg->bw_high_hz - cfg->bw_low_hz;
}

/* ==================================================================
 * L6: Current Shunt Measurement
 *
 * Current sensing using a shunt resistor:
 *   V_shunt = I_load * R_shunt
 *
 * Shunt selection trade-off:
 *   - Larger R_shunt: higher signal, less gain needed, but more
 *     power dissipation and voltage drop
 *   - Smaller R_shunt: lower power loss, but need higher IA gain
 *
 * For a 10A load with R_shunt = 1 mOhm:
 *   V_shunt = 10 * 0.001 = 10 mV (at full load)
 *   Power loss = I^2 * R = 100 * 0.001 = 0.1 W
 *
 * Common-mode voltage can be very high (e.g., 48V bus).
 * The IA must reject this large CM voltage. High-side sensing
 * requires high CMRR IAs or dedicated current-sense amplifiers.
 *
 * Reference: Zetex AN39, "Current Measurement Applications Handbook"
 * Complexity: O(1)
 * ================================================================== */

typedef struct {
    double r_shunt_ohm;         /* Shunt resistance (Ohm) */
    double i_max_a;             /* Maximum load current (A) */
    double v_cm_v;              /* Common-mode voltage (V) */
    bool high_side;             /* High-side (true) or low-side sensing */
} current_shunt_config_t;

double current_shunt_vout_at_current(const current_shunt_config_t *cfg,
                                      double i_load_a) {
    if (!cfg) return 0.0;
    return i_load_a * cfg->r_shunt_ohm;
}

double current_shunt_required_gain(const current_shunt_config_t *cfg,
                                    double v_adc_fullscale) {
    if (!cfg) return 0.0;
    double v_fs = current_shunt_vout_at_current(cfg, cfg->i_max_a);
    if (v_fs <= 0.0) return 0.0;
    return v_adc_fullscale / v_fs;
}

double current_shunt_power_loss(const current_shunt_config_t *cfg, double i_a) {
    if (!cfg) return 0.0;
    return i_a * i_a * cfg->r_shunt_ohm;
}

/* ==================================================================
 * L6: Hall Effect Sensor Interface
 *
 * Hall sensors produce a voltage proportional to magnetic field:
 *   V_hall = S * B * I_bias
 *   where S = sensitivity (mV/mT at 1 mA bias)
 *
 * Typical Hall sensor (Allegro A1324):
 *   - Sensitivity: 5 mV/Gauss at 5V supply
 *   - Quiescent output: 2.5V (half supply)
 *   - Bandwidth: DC-20 kHz
 *
 * The IA must handle a large common-mode DC offset (the quiescent
 * output) while amplifying small field-dependent variations.
 *
 * Complexity: O(1)
 * ================================================================== */

typedef struct {
    double sensitivity_mv_per_mt;  /* Sensitivity (mV/mT) */
    double quiescent_output_v;     /* Zero-field output voltage */
    double bias_current_ma;        /* Sensor bias current */
} hall_sensor_config_t;

double hall_field_from_voltage(double v_out, double v_quiescent,
                                double sensitivity_mv_per_mt) {
    if (sensitivity_mv_per_mt <= 0.0) return 0.0;
    return (v_out - v_quiescent) / (sensitivity_mv_per_mt * 1e-3);
}

/* ==================================================================
 * L6: Sensor Excitation and Ratiometric Configuration
 *
 * Many bridge sensors are used ratiometrically: both the bridge
 * excitation and the ADC reference are derived from the same
 * voltage source. This cancels excitation drift:
 *
 *   ADC_reading = (V_bridge / V_exc) * (2^N - 1)
 *               = k * measurand * (2^N - 1)
 *
 * The measurand is directly proportional to the ADC code,
 * independent of V_exc variations.
 *
 * This is the gold standard for precision sensor measurements
 * and eliminates the need for a precision voltage reference
 * for the bridge excitation.
 *
 * Reference: Analog Devices, "Ratiometric Measurement Techniques"
 * Complexity: O(1)
 * ================================================================== */
double ratiometric_adc_code(double v_bridge, double v_exc, int adc_bits) {
    if (v_exc <= 0.0) return 0.0;
    return (v_bridge / v_exc) * ((1 << adc_bits) - 1);
}

double ratiometric_measurand_from_code(int adc_code, int adc_bits,
                                        double fullscale_value) {
    int max_code = (1 << adc_bits) - 1;
    if (max_code <= 0) return 0.0;
    return ((double)adc_code / max_code) * fullscale_value;
}

/* ==================================================================
 * L6/L7: Complete Sensor Readout Pipeline
 *
 * Generic sensor readout sequence:
 *   1. Configure IA gain based on sensor type and range
 *   2. Apply excitation (if bridge sensor)
 *   3. Wait for settling
 *   4. Acquire N samples from ADC
 *   5. Apply digital filtering (moving average, notch 50/60Hz)
 *   6. Apply calibration correction
 *   7. Convert to engineering units
 *
 * Complexity: O(N) for N-sample acquisition
 * ================================================================== */

typedef struct {
    double *buffer;             /* Sample buffer */
    int capacity;               /* Buffer capacity */
    int count;                  /* Current number of samples */
} sensor_sample_buffer_t;

void sensor_buffer_init(sensor_sample_buffer_t *buf, double *storage, int cap) {
    if (!buf || !storage || cap <= 0) return;
    buf->buffer = storage;
    buf->capacity = cap;
    buf->count = 0;
}

void sensor_buffer_add(sensor_sample_buffer_t *buf, double sample) {
    if (!buf || buf->count >= buf->capacity) return;
    buf->buffer[buf->count++] = sample;
}

double sensor_buffer_mean(const sensor_sample_buffer_t *buf) {
    if (!buf || buf->count <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < buf->count; i++) sum += buf->buffer[i];
    return sum / buf->count;
}

void sensor_buffer_clear(sensor_sample_buffer_t *buf) {
    if (!buf) return;
    buf->count = 0;
}

/* ==================================================================
 * L6: Notch Filter for 50/60 Hz Power Line Rejection
 *
 * A simple digital notch filter to remove mains interference:
 *   y[n] = (x[n] + x[n-2]) / 2  (simplified)
 *
 * For a proper notch at f_notch with sampling rate f_s:
 *   y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2]
 *        - a1*y[n-1] - a2*y[n-2]
 *
 * Coefficients for 50 Hz notch at f_s = 500 Hz:
 *   w0 = 2*pi*50/500 = 0.6283
 *   Q = 30 (quality factor, narrow notch)
 *   b0 = 1, b1 = -2*cos(w0), b2 = 1
 *   a1 = -2*cos(w0), a2 = 1 - sin(w0)/(2*Q) ... (simplified)
 *
 * For many applications, a simple moving average of 10-20 samples
 * at the mains period (20 ms for 50 Hz, 16.67 ms for 60 Hz)
 * provides excellent rejection of mains harmonics.
 *
 * Complexity: O(1) per sample
 * ================================================================== */

/**
 * Synchronous moving average filter for mains rejection.
 * If f_sample is an integer multiple of f_mains, averaging over
 * one full mains period cancels the fundamental and all harmonics.
 *
 * For 50 Hz mains, sample at 500 Hz: average 10 samples -> 20 ms period
 * For 60 Hz mains, sample at 600 Hz: average 10 samples -> 16.67 ms period
 */
double mains_notch_filter(double current_sample, double prev_average,
                           int period_samples, int *sample_index) {
    if (!sample_index || period_samples <= 0) return current_sample;

    /* Simple exponential moving average with window = 1 period */
    double alpha = 2.0 / (period_samples + 1.0);
    (void)sample_index;
    return prev_average + alpha * (current_sample - prev_average);
}