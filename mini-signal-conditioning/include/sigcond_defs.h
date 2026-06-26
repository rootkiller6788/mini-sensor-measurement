/**
 * sigcond_defs.h - Signal Conditioning Core Type Definitions
 *
 * Defines all core data types, enumerations, and structures for
 * the signal conditioning pipeline: amplification, filtering,
 * linearization, excitation, isolation, and CJC.
 *
 * L1 Definitions: signal conditioning stages, filter types,
 * excitation modes, isolation topologies, linearization methods.
 *
 * Courses: MIT 6.003, Stanford EE247, Berkeley EE105, Michigan EECS 411,
 *   Georgia Tech ECE 4270, TU Munich, ETH 227-0427, Tsinghua
 *
 * Reference:
 *   Pallas-Areny & Webster, "Sensors and Signal Conditioning" (2001)
 *   Kester, "Analog-Digital Conversion" (Analog Devices, 2004)
 *   Franco, "Design with Operational Amplifiers and Analog ICs" (2015)
 *   Horowitz & Hill, "The Art of Electronics" (3rd ed., 2015), Ch. 15
 */
#ifndef SIGCOND_DEFS_H
#define SIGCOND_DEFS_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==================================================================
 * L1: Pipeline Stages (10 stages)
 *
 * Sensor >> Protection >> Excitation >> Amplification >> Filtering >>
 * Linearization >> Isolation >> Offset Null >> CJC >> ADC Buffer
 * ================================================================== */
typedef enum {
    SIGCOND_STAGE_PROTECTION    = 0,
    SIGCOND_STAGE_EXCITATION    = 1,
    SIGCOND_STAGE_AMPLIFICATION = 2,
    SIGCOND_STAGE_FILTERING     = 3,
    SIGCOND_STAGE_LINEARIZATION = 4,
    SIGCOND_STAGE_ISOLATION     = 5,
    SIGCOND_STAGE_OFFSET_NULL   = 6,
    SIGCOND_STAGE_CJC           = 7,
    SIGCOND_STAGE_BRIDGE_BAL    = 8,
    SIGCOND_STAGE_ADC_BUFFER    = 9,
    SIGCOND_STAGE_COUNT           = 10
} sigcond_stage_t;

/* ==================================================================
 * L1: Filter Approximation Families (10 types)
 *
 * Butterworth (1930): |H(jw)|^2 = 1/(1 + (w/wc)^(2n))
 *   Maximally flat passband, monotonic roll-off.
 * Chebyshev I: |H(jw)|^2 = 1/(1 + e^2 * Tn^2(w/wc))
 *   Passband ripple, faster roll-off than Butterworth.
 * Chebyshev II: monotonic passband, stopband ripple.
 * Elliptic (Cauer): ripple in both bands, steepest transition.
 * Bessel (Thomson): maximally flat group delay.
 * Gaussian: time-domain step response, no overshoot.
 * Moving Average: boxcar FIR, sin(pi*f/fs)/sin(pi*f/fs/N) response.
 * Notch/Bandstop: narrow rejection band (e.g., 50/60 Hz).
 *
 * Reference: Zverev (1967), Van Valkenburg (1982), Sedra & Smith (2020)
 * ================================================================== */
typedef enum {
    FILTER_NONE            = 0,
    FILTER_BUTTERWORTH     = 1,
    FILTER_CHEBYSHEV_I     = 2,
    FILTER_CHEBYSHEV_II    = 3,
    FILTER_ELLIPTIC        = 4,
    FILTER_BESSEL          = 5,
    FILTER_GAUSSIAN        = 6,
    FILTER_SYNC            = 7,
    FILTER_MOVING_AVERAGE  = 8,
    FILTER_NOTCH_BANDSTOP  = 9,
    FILTER_COUNT             = 10
} filter_family_t;

/* ==================================================================
 * L1: Active-RC Filter Topologies (8 types)
 *
 * Sallen-Key: positive feedback, single op-amp, low Z_out.
 * MFB (Rauch): inverting, low component sensitivity.
 * State-Variable (KHN/Tow-Thomas): three op-amps, tunable Q.
 * Twin-T: passive notch + op-amp for Q enhancement.
 * GIC: inductor simulation without magnetics.
 * FDNR: frequency-dependent negative resistance.
 * Switched-Cap: capacitor ratios replace resistors.
 * Gm-C: transconductance-C, GHz operation.
 *
 * Reference: Huelsman & Allen (1980)
 * ================================================================== */
typedef enum {
    TOPO_SALLEN_KEY        = 0,
    TOPO_MULTIPLE_FEEDBACK = 1,
    TOPO_STATE_VARIABLE    = 2,
    TOPO_TWIN_T            = 3,
    TOPO_GIC_BASED         = 4,
    TOPO_FDNR              = 5,
    TOPO_SWITCHED_CAP      = 6,
    TOPO_GM_C              = 7,
    TOPO_COUNT               = 8
} filter_topology_t;

/** Filter response type: Lowpass, Highpass, Bandpass, Bandstop, Allpass */
typedef enum {
    FILTER_RESP_LOWPASS     = 0,
    FILTER_RESP_HIGHPASS    = 1,
    FILTER_RESP_BANDPASS    = 2,
    FILTER_RESP_BANDSTOP    = 3,
    FILTER_RESP_ALLPASS     = 4,
    FILTER_RESP_COUNT         = 5
} filter_response_t;

/* ==================================================================
 * L1: Linearization Methods (10 methods)
 *
 * Steinhart-Hart: 1/T = A + B*ln(R) + C*(ln(R))^3 (NTC thermistor)
 * Callendar-Van Dusen: R(T) = R0*(1 + A*T + B*T^2 + ...) (RTD)
 * NIST ITS-90: published 8th-14th degree polynomials for TCs
 * Rational (Pade): y = P(x)/Q(x), better extrapolation
 *
 * Reference: Pallas-Areny & Webster (2001), Ch. 3
 * ================================================================== */
typedef enum {
    LINEARIZE_NONE             = 0,
    LINEARIZE_LOOKUP_TABLE     = 1,
    LINEARIZE_POLYNOMIAL       = 2,
    LINEARIZE_PIECEWISE_LINEAR = 3,
    LINEARIZE_PIECEWISE_POLY   = 4,
    LINEARIZE_STEINHART_HART   = 5,
    LINEARIZE_CALLENDAR_VAN_DUSEN = 6,
    LINEARIZE_NIST_ITS90       = 7,
    LINEARIZE_RATIONAL         = 8,
    LINEARIZE_PHYSICAL_MODEL   = 9,
    LINEARIZE_COUNT              = 10
} linearize_method_t;

/* ==================================================================
 * L1: Excitation Types (7 types)
 *
 * Ratiometric measurement: Vout = k * Vref * measurand
 *   Cancels Vref drift from the measurement equation.
 * AC excitation needed for: capacitive sensors, LVDT, inductive.
 * Pulsed excitation: low-power, duty-cycled sensors.
 *
 * Reference: Pallas-Areny & Webster (2001), Ch. 2
 * ================================================================== */
typedef enum {
    EXCITATION_NONE           = 0,
    EXCITATION_CONST_VOLTAGE  = 1,
    EXCITATION_CONST_CURRENT  = 2,
    EXCITATION_AC_SINE        = 3,
    EXCITATION_AC_PULSE       = 4,
    EXCITATION_RATIOMETRIC    = 5,
    EXCITATION_CURRENT_PULSE  = 6,
    EXCITATION_COUNT            = 7
} excitation_type_t;

/* ==================================================================
 * L1: Isolation Barrier Types (7 types)
 *
 * Galvanic isolation blocks DC/low-freq current.
 * Required for: patient safety (IEC 60601-1, 4kV),
 *   ground loop elimination, high CMRR.
 * iCoupler (ADuM): chip-scale micro-transformers, >100kV/us CMTI.
 *
 * Reference: Kester (1999), Ch. 10
 * ================================================================== */
typedef enum {
    ISOLATION_NONE            = 0,
    ISOLATION_OPTICAL         = 1,
    ISOLATION_TRANSFORMER     = 2,
    ISOLATION_CAPACITIVE      = 3,
    ISOLATION_MAGNETIC        = 4,
    ISOLATION_HALL            = 5,
    ISOLATION_RF              = 6,
    ISOLATION_COUNT             = 7
} isolation_type_t;

/* ==================================================================
 * L1: Bridge Topologies (7 types)
 *
 * Quarter bridge: Vout = Vexc * dR/(4R), single active arm
 * Half bridge: Vout = Vexc * dR/(2R), adjacent arms
 * Full bridge: Vout = Vexc * dR/R, max sensitivity
 * 3-wire: single lead compensation for RTD
 * 4-wire Kelvin: eliminates lead resistance entirely
 * ================================================================== */
typedef enum {
    BRIDGE_QUARTER_SINGLE     = 0,
    BRIDGE_QUARTER_3WIRE      = 1,
    BRIDGE_QUARTER_4WIRE      = 2,
    BRIDGE_HALF_ADJACENT      = 3,
    BRIDGE_HALF_OPPOSITE      = 4,
    BRIDGE_FULL               = 5,
    BRIDGE_NONE                 = 6
} bridge_topology_t;

/* ==================================================================
 * L1: Error Classification (13 error sources)
 *
 * RSS combination per ISO GUM:
 *   u_total = sqrt(sum(u_i^2))
 * Expanded (k=2): U = 2 * u_total
 * ================================================================== */
typedef enum {
    COND_ERR_OFFSET             = 0,
    COND_ERR_GAIN               = 1,
    COND_ERR_NONLINEARITY       = 2,
    COND_ERR_THERMAL_DRIFT      = 3,
    COND_ERR_NOISE              = 4,
    COND_ERR_CMRR               = 5,
    COND_ERR_PSRR               = 6,
    COND_ERR_EXCITATION_DRIFT   = 7,
    COND_ERR_FILTER_ATTEN       = 8,
    COND_ERR_FILTER_PHASE       = 9,
    COND_ERR_ISOLATION_LEAK     = 10,
    COND_ERR_QUANTIZATION       = 11,
    COND_ERR_ALIASING           = 12,
    COND_ERR_COUNT                = 13
} sigcond_error_type_t;

/* ==================================================================
 * L1: Sensor Input Types (17 types)
 * ================================================================== */
typedef enum {
    SENSOR_TEMPERATURE_TC      = 0,
    SENSOR_TEMPERATURE_RTD     = 1,
    SENSOR_TEMPERATURE_NTC     = 2,
    SENSOR_TEMPERATURE_IC      = 3,
    SENSOR_STRAIN              = 4,
    SENSOR_PRESSURE_WHEATSTONE = 5,
    SENSOR_PRESSURE_PIEZO      = 6,
    SENSOR_FORCE_LOADCELL      = 7,
    SENSOR_TORQUE              = 8,
    SENSOR_POSITION_LVDT       = 9,
    SENSOR_POSITION_POT        = 10,
    SENSOR_PHOTO_CONDUCTIVE    = 11,
    SENSOR_PHOTO_VOLTAIC       = 12,
    SENSOR_HALL                = 13,
    SENSOR_CURRENT_SHUNT       = 14,
    SENSOR_ACCEL_MEMS          = 15,
    SENSOR_HUMIDITY_CAP        = 16
} sensor_input_type_t;

/* ==================================================================
 * L1: Analog Filter Specification
 *
 * Q-factor: Q = 1/(2*zeta) = f0 / BW_3dB
 * Butterworth 2nd-order: Q = 0.7071 (max flat magnitude)
 * Bessel 2nd-order:      Q = 0.5774 (max flat delay)
 * Chebyshev 2nd-order:   Q depends on ripple
 * ================================================================== */
typedef struct {
    filter_family_t     family;
    filter_response_t   response;
    filter_topology_t   topology;
    unsigned            order;
    double              fc_hz;
    double              fc2_hz;
    double              passband_ripple_db;
    double              stopband_attn_db;
    double              stopband_freq_hz;
    double              q_factor;
    double              group_delay_us;
    double              dc_gain_vv;
} analog_filter_spec_t;

/* ==================================================================
 * L1: Digital Filter Specification
 *
 * IIR difference equation:
 *   y[n] = sum(b_k * x[n-k]) - sum(a_k * y[n-k])
 * FIR: all a_k = 0, linear phase achievable with symmetric coefficients.
 *
 * Kaiser window FIR length estimate:
 *   N = (A - 7.95) / (14.36 * delta_f/fs)
 *   where A = stopband_attenuation_dB, delta_f = transition width
 * ================================================================== */
typedef struct {
    filter_family_t     family;
    filter_response_t   response;
    unsigned            num_taps;
    double              fs_hz;
    double              fc_hz;
    double              fc2_hz;
    double              passband_ripple_db;
    double              stopband_attn_db;
    unsigned            filter_length;
    double              coefficients[32];
    double              denominator[32];
    unsigned            num_numerator;
    unsigned            num_denominator;
} digital_filter_spec_t;

/* ==================================================================
 * L1: Channel Configuration (master config struct)
 * ================================================================== */
typedef struct {
    unsigned            channel_id;
    sensor_input_type_t sensor_type;
    const char*         channel_name;
    excitation_type_t   excitation;
    double              excitation_v;
    double              excitation_current_a;
    bool                ratiometric;
    bool                ac_excitation;
    double              ac_excitation_freq;
    double              gain_vv;
    double              v_ref;
    analog_filter_spec_t pre_filter;
    analog_filter_spec_t post_filter;
    bool                notch_50hz_enable;
    bool                notch_60hz_enable;
    linearize_method_t  linearize_method;
    bool                linearize_enable;
    isolation_type_t    isolation;
    double              isolation_withstand_v;
    bool                cjc_enable;
    bool                cjc_sensor_internal;
    double              adc_fs;
    unsigned            adc_resolution_bits;
    double              sample_rate_hz;
    bridge_topology_t   bridge_topo;
    double              bridge_r_nominal;
    double              bridge_v_excitation;
    double              ovp_threshold_v;
    double              input_z_min;
} sigcond_channel_config_t;

/* ==================================================================
 * L1: Support Structures
 * ================================================================== */
typedef struct {
    bool                channels_configured;
    bool                calibration_valid;
    bool                self_test_passed;
    bool                overtemperature;
    bool                power_good;
    unsigned            active_channels;
    double              board_temp_c;
    double              supply_voltage_v;
} sigcond_system_status_t;

/** Biquadratic filter section: H(s) = k*(s^2 + (wz/Qz)s + wz^2)/(s^2 + (wp/Qp)s + wp^2) */
typedef struct {
    double section_gain, pole_freq_hz, pole_q;
    double zero_freq_hz, zero_q;
    double coeff_b0, coeff_b1, coeff_b2;
    double coeff_a1, coeff_a2;
} filter_biquad_t;

/** Single calibration data point */
typedef struct {
    double measurand_value, sensor_output_raw;
    double conditioned_output, temperature_c, humidity_pct;
} calibration_point_t;

/** Polynomial: y = c[0] + c[1]*x + ... + c[N]*x^N */
typedef struct {
    unsigned degree;
    double coefficients[15];
    double valid_min, valid_max, residual_error_max;
} polynomial_coeffs_t;

/** Steinhart-Hart: 1/T = A + B*ln(R) + C*(ln(R))^3 */
typedef struct {
    double a_coeff, b_coeff, c_coeff, d_coeff;
    double r_ref, t_ref_c, beta_k;
    bool use_4param;
} steinhart_hart_params_t;

/** RTD Callendar-Van Dusen (IEC 60751:2008) */
typedef struct {
    double r0, a_coeff, b_coeff, c_coeff;
    double alpha, delta, beta;
    bool above_zero;
    double nominal_current_a;
} rtd_cvd_params_t;

/** ITS-90 thermocouple range model */
typedef struct {
    unsigned type_code;
    double seebeck_25c, temp_min_c, temp_max_c;
    double voltage_min_uv, voltage_max_uv;
    polynomial_coeffs_t forward_poly, inverse_poly;
    double error_max_c, error_typical_c;
} thermocouple_nist_model_t;

/** Isolation specification */
typedef struct {
    isolation_type_t barrier_type;
    double withstand_vrms, working_vrms, surge_kv;
    double cmrr_db, cmrr_at_50hz_db;
    double barrier_cap_pf, leakage_ua;
    double cmti_kv_per_us, mtbf_hours;
    unsigned creepage_mm, clearance_mm;
} isolation_spec_t;

/** Error budget (ISO GUM aggregation) */
typedef struct {
    double offset_error_uv, gain_error_ppm, nonlinearity_ppm;
    double noise_rti_uvrms, noise_rti_uvpp;
    double cmrr_error_uv, psrr_error_uv;
    double thermal_drift_uv, excitation_drift_ppm;
    double filter_atten_error_db, isolation_leakage_uv;
    double quantization_error_uv, aliasing_error_uv;
    double total_uncertainty_uv, expanded_uncertainty_uv;
} error_budget_t;

/** Wheatstone bridge completion network */
typedef struct {
    double r1, r2, r3;
    double r_balance_min, r_balance_max;
    double r_lead_per_wire;
    unsigned num_wires;
    double tcr_ppm_per_c;
} bridge_completion_t;

/** Cold junction compensation (thermocouple) */
typedef struct {
    bool enabled;
    bool internal_sensor;
    double cj_temp_c;
    double cj_temp_accuracy_c;
    double cj_voltage_uv;
    polynomial_coeffs_t cjc_sensor_poly;
} cjc_config_t;

#ifdef __cplusplus
}
#endif
#endif /* SIGCOND_DEFS_H */
