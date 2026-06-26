/**
 * inst_amp_core.c - Three-Op-Amp Instrumentation Amplifier Core
 *
 * L2 Core Concepts: Complete implementation of the classic 3-op-amp
 * instrumentation amplifier. Includes gain computation, R_G selection,
 * CMRR analysis, offset voltage correction, and output swing validation.
 *
 * L4 Fundamental Laws:
 *   - Gain: G = (1 + 2*R1/R_G) * (R3/R2)
 *   - CMRR_resistors = (1 + G_diff) / (4 * epsilon)
 *   - V_out = G * V_diff + V_ref
 *   - Output swing: V_out in [V_ee+headroom, V_cc-headroom]
 *
 * Courses: Stanford EE247, Berkeley EE105, MIT 6.630
 * Reference:
 *   Kitchin & Counts (2006) "A Designer's Guide to Instrumentation Amplifiers"
 *   Sedra & Smith (2020) "Microelectronic Circuits" Ch. 2.4
 *   Pall++s-Areny & Webster (1991) "CMRR in Differential Amplifiers", IEEE TIM
 */

#include "inst_amp_defs.h"
#include "inst_amp_topology.h"
#include <stdlib.h>
#include <string.h>

/* ==================================================================
 * L2/L4: Three-Op-Amp Gain Computation
 *
 * Theorem (L4 - Ideal 3-op-amp gain):
 *   The output voltage of a 3-op-amp IA is:
 *     V_out = (V_in+ - V_in-) * (1 + 2*R1/R_G) * (R3/R2) + V_ref
 *
 *   When R2 = R3 (unity-gain diff amp):
 *     V_out = (V_in+ - V_in-) * (1 + 2*R1/R_G) + V_ref
 *
 * Proof outline:
 *   Stage 1 (A1, A2): V_a = V_in-*(1+R1/R_G) - V_in+*(R1/R_G)
 *                      V_b = V_in+*(1+R1/R_G) - V_in-*(R1/R_G)
 *                      V_a - V_b = (V_in+ - V_in-)*(1 + 2*R1/R_G)
 *   Stage 2 (A3):      V_out = (V_b - V_a)*(R3/R2) + V_ref
 *   Substitute:        V_out = (V_in+ - V_in-)*(1 + 2*R1/R_G)*(R3/R2) + V_ref
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_three_opamp_gain(const ia_three_opamp_config_t *cfg) {
    if (!cfg) return 0.0;
    if (cfg->rg_ohm <= 0.0) return 0.0;
    if (cfg->r2_ohm <= 0.0 || cfg->r3_ohm <= 0.0) return 0.0;

    double g_input_stage = 1.0 + (2.0 * cfg->r1_ohm / cfg->rg_ohm);
    double g_diff_amp = cfg->r3_ohm / cfg->r2_ohm;
    return g_input_stage * g_diff_amp;
}

/* ==================================================================
 * L4: CMRR Limited by Resistor Mismatch
 *
 * Theorem (Pall++s-Areny & Webster, 1991):
 *   For a difference amplifier with resistors R2 (input), R3 (feedback),
 *   the CMRR due to resistor mismatch is:
 *     CMRR_res = (1 + R3/R2) / (4 * epsilon)
 *   where epsilon = |deltaR|/R is the relative resistor tolerance.
 *
 *   CMRR_total = [1/CMRR_opamp + 1/CMRR_res]^(-1)
 *
 *   Typical values: 0.1% resistors -> CMRR_res = 1.01/(4*0.001) = 252.5 = 48 dB
 *                   0.01% resistors -> CMRR_res = 68 dB
 *
 *   For 3-op-amp IA with R2=R3: CMRR_total ~ CMRR_A3 / (1 + CMRR_A3/CMRR_res)
 *   The input buffers (A1, A2) see common-mode signal equally, so their
 *   CMRR contributions are largely rejected by the diff amp.
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_three_opamp_cmrr_resistor_limit(const ia_three_opamp_config_t *cfg) {
    if (!cfg) return 0.0;

    double g_diff = cfg->r3_ohm / cfg->r2_ohm;
    double epsilon = cfg->resistor_tol_pct / 100.0;
    if (epsilon <= 0.0) return INFINITY;

    double cmrr_linear = (1.0 + g_diff) / (4.0 * epsilon);
    return 20.0 * log10(cmrr_linear);  /* Return in dB */
}

/* ==================================================================
 * L5: R_G Selection for Target Gain
 *
 * Algorithm (L5 - Resistor selection for 3-op-amp IA):
 *   For R2 = R3 and desired gain G_target:
 *     R_G = 2 * R1 / (G_target - 1)
 *
 *   This is derived from G = 1 + 2*R1/R_G (for R2=R3).
 *   R_G must be positive. For G_target arbitrarily close to 1,
 *   R_G approaches infinity (open circuit).
 *
 *   Practical constraints:
 *     - R_G_min: limited by op-amp output drive (typically > 100 Ohm)
 *     - R_G_max: limited by input bias current * R_G dropping across R_G
 *     - Gain range: typically 1 to 10000
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_three_opamp_rg_for_gain(const ia_three_opamp_config_t *cfg, double target_gain) {
    if (!cfg) return -1.0;
    if (target_gain <= 1.0) return -1.0;

    return (2.0 * cfg->r1_ohm) / (target_gain - 1.0);
}

/* ==================================================================
 * L2: Output Voltage Computation for 3-Op-Amp IA
 *
 * Computes the output voltage given differential and common-mode
 * input voltages, including non-ideal effects.
 *
 * Total output (L4):
 *   V_out = G * V_diff + V_ref + G * V_cm / CMRR
 *   where CMRR is in linear units (10^(CMRR_dB/20)).
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_three_opamp_output(const ia_three_opamp_config_t *cfg,
                              double v_in_plus, double v_in_minus,
                              double v_ref, double cmrr_db) {
    if (!cfg) return 0.0;

    double G = ia_three_opamp_gain(cfg);
    double v_diff = v_in_plus - v_in_minus;
    double v_cm = (v_in_plus + v_in_minus) / 2.0;

    /* Common-mode error: V_cm / CMRR_linear referred to input */
    double cmrr_linear = pow(10.0, cmrr_db / 20.0);
    double v_cm_error = v_cm / cmrr_linear;

    return G * (v_diff + v_cm_error) + v_ref;
}

/* ==================================================================
 * L2: Check if output is within valid swing range
 *
 * Most rail-to-rail IAs need some headroom (typically 100-200mV
 * from each rail). This function checks if the computed output
 * can be delivered without saturation/clipping.
 *
 * Complexity: O(1)
 * ================================================================== */
bool ia_output_in_range(double v_out, double v_ee, double v_cc,
                         double headroom_mv) {
    double hr = headroom_mv / 1000.0;
    return (v_out >= (v_ee + hr)) && (v_out <= (v_cc - hr));
}

/* ==================================================================
 * L5: Input Offset Voltage Compensation
 *
 * Two-point calibration method:
 *   1. Apply known input V_cal1, measure V_out1
 *   2. Apply known input V_cal2, measure V_out2
 *   3. Solve: V_out = G * (V_in + V_os)
 *      V_os = (V_out1*V_cal2 - V_out2*V_cal1) / (G*(V_out1 - V_out2)) ... no
 *
 *   Correct method:
 *     slope = (V_out2 - V_out1) / (V_cal2 - V_cal1)  = actual gain
 *     offset_at_input = V_cal1 - V_out1/slope
 *
 *   Corrected output: V_out_corr = (V_out_raw - offset_at_output) * (G_nominal/G_actual)
 *   where offset_at_output = offset_at_input * G_actual
 *
 * Complexity: O(1)
 * ================================================================== */
void ia_calibrate_two_point(double v_cal1, double v_out1,
                             double v_cal2, double v_out2,
                             double *actual_gain, double *input_offset) {
    if (!actual_gain || !input_offset) return;

    double dv_in = v_cal2 - v_cal1;
    double dv_out = v_out2 - v_out1;

    if (fabs(dv_in) < 1e-12) {
        *actual_gain = 0.0;
        *input_offset = 0.0;
        return;
    }

    *actual_gain = dv_out / dv_in;
    *input_offset = v_cal1 - v_out1 / (*actual_gain);
}

/* ==================================================================
 * L5: Apply calibration to raw output reading
 *
 * V_out_corrected = (V_out_raw - V_os_output) * (G_nominal / G_actual)
 * where V_os_output = V_os_input * G_actual
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_apply_calibration(double v_out_raw, double actual_gain,
                             double nominal_gain, double input_offset) {
    double offset_at_output = input_offset * actual_gain;
    double v_corrected = v_out_raw - offset_at_output;
    if (actual_gain > 0.0) {
        v_corrected *= (nominal_gain / actual_gain);
    }
    return v_corrected;
}

/* ==================================================================
 * L3: Total Input-Referred Noise Computation
 *
 * Three noise sources combine:
 *   1. Amplifier voltage noise: en_amp [nV/rtHz]
 *   2. Amplifier current noise flowing through source: in_amp * R_s [nV/rtHz]
 *   3. Johnson-Nyquist (thermal) noise of source resistance:
 *      en_thermal = sqrt(4 * k * T * R_s) [nV/rtHz]
 *      where k = 1.380649e-23 J/K, T in Kelvin
 *
 *   en_total = sqrt(en_amp^2 + (in_amp*R_s)^2 + 4*k*T*R_s)
 *
 *   For total RMS noise over bandwidth BW:
 *   Vn_rms = en_total * sqrt(BW)
 *
 * Reference: Motchenbacher & Connelly, "Low-Noise Electronic System Design" (1993)
 * Complexity: O(1)
 * ================================================================== */
double ia_total_input_noise(double en_amp_nv_rhz, double in_amp_pa_rhz,
                             double rs_ohm, double bw_hz, double temp_c) {
    double k_B = 1.380649e-23;
    double T_K = temp_c + 273.15;

    /* Convert current noise from pA/rtHz to nV/rtHz through R_s */
    double en_current = (in_amp_pa_rhz * 1e-12) * rs_ohm * 1e9;

    /* Thermal (Johnson) noise of source resistance */
    double en_thermal_nv_rhz = sqrt(4.0 * k_B * T_K * rs_ohm) * 1e9;

    /* Total noise density */
    double en_total_nv_rhz = sqrt(en_amp_nv_rhz * en_amp_nv_rhz
                                + en_current * en_current
                                + en_thermal_nv_rhz * en_thermal_nv_rhz);

    /* RMS noise over bandwidth */
    return en_total_nv_rhz * sqrt(bw_hz) * 1e-9;  /* Return in Vrms */
}

/* ==================================================================
 * L2: Gain Error Calculation
 *
 * Given nominal and actual gains, compute the percentage error.
 *   gain_error_pct = (G_actual - G_nominal) / G_nominal * 100
 *
 * This includes contributions from:
 *   - R_G tolerance (+/- 0.1% to 1%)
 *   - Internal resistor matching
 *   - Op-amp finite open-loop gain
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_gain_error_percent(double g_nominal, double g_actual) {
    if (g_nominal <= 0.0) return 0.0;
    return (g_actual - g_nominal) / g_nominal * 100.0;
}

/* ==================================================================
 * L2: Effective Number of Bits (ENOB) at IA Output
 *
 * ENOB = (SNR_dB - 1.76) / 6.02
 * where SNR = 20*log10(V_signal_rms / V_noise_rms)
 *
 * This determines the resolution limited by analog front-end noise.
 * For a 16-bit ADC with 10V FS: 1 LSB = 153 uV.
 * IA noise should be < 1 LSB for the LSB to be meaningful.
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_effective_bits(double v_signal_rms, double v_noise_rms) {
    if (v_noise_rms <= 0.0) return 24.0;  /* Theoretical limit */
    double snr_db = 20.0 * log10(v_signal_rms / v_noise_rms);
    return (snr_db - 1.76) / 6.02;
}

/* ==================================================================
 * L2: CMRR Frequency Roll-off Model
 *
 * CMRR degrades with frequency. A simple single-pole model:
 *   CMRR(f) = CMRR_DC / sqrt(1 + (f/f_cmrr_pole)^2)
 *   CMRR_dB(f) = CMRR_DC_dB - 20*log10(sqrt(1 + (f/f_cmrr_pole)^2))
 *
 * For the 2-op-amp topology, the CMRR pole is:
 *   f_pole = GBW / G  (approximately)
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_cmrr_at_frequency(double cmrr_dc_db, double f_hz,
                             double f_pole_hz) {
    if (f_pole_hz <= 0.0) return cmrr_dc_db;
    double rolloff = 10.0 * log10(1.0 + (f_hz / f_pole_hz) * (f_hz / f_pole_hz));
    return cmrr_dc_db - rolloff;
}

/* ==================================================================
 * L2: Power Supply Rejection Ratio (PSRR) Error
 *
 * Supply ripple appears at output attenuated by PSRR:
 *   V_out_ripple = V_supply_ripple / 10^(PSRR_dB/20)
 *
 * Reference to input: divide by gain G.
 *   V_rti_ripple = V_out_ripple / G
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_psrr_output_error(double v_supply_ripple, double psrr_db) {
    double psrr_linear = pow(10.0, psrr_db / 20.0);
    return v_supply_ripple / psrr_linear;
}

/* ==================================================================
 * L5: Select R_G from Standard E96 Resistor Values
 *
 * Converts computed R_G to the nearest standard E96 (1%) value.
 * E96 series: 1.00, 1.02, 1.05, 1.07, ..., 9.76 (96 values/decade).
 *
 * Algorithm: round to nearest standard value using decade mantissa.
 *
 * Complexity: O(log(E96)) = O(1)  (binary search in 96-element table)
 * ================================================================== */
static const double e96_mantissas[96] = {
    1.00, 1.02, 1.05, 1.07, 1.10, 1.13, 1.15, 1.18, 1.21, 1.24,
    1.27, 1.30, 1.33, 1.37, 1.40, 1.43, 1.47, 1.50, 1.54, 1.58,
    1.62, 1.65, 1.69, 1.74, 1.78, 1.82, 1.87, 1.91, 1.96, 2.00,
    2.05, 2.10, 2.15, 2.21, 2.26, 2.32, 2.37, 2.43, 2.49, 2.55,
    2.61, 2.67, 2.74, 2.80, 2.87, 2.94, 3.01, 3.09, 3.16, 3.24,
    3.32, 3.40, 3.48, 3.57, 3.65, 3.74, 3.83, 3.92, 4.02, 4.12,
    4.22, 4.32, 4.42, 4.53, 4.64, 4.75, 4.87, 4.99, 5.11, 5.23,
    5.36, 5.49, 5.62, 5.76, 5.90, 6.04, 6.19, 6.34, 6.49, 6.65,
    6.81, 6.98, 7.15, 7.32, 7.50, 7.68, 7.87, 8.06, 8.25, 8.45,
    8.66, 8.87, 9.09, 9.31, 9.53, 9.76
};

double ia_e96_nearest(double r_ohm) {
    if (r_ohm <= 0.0) return 0.0;

    /* Find decade and mantissa */
    double decade = 1.0;
    while (r_ohm >= 10.0) { r_ohm /= 10.0; decade *= 10.0; }
    while (r_ohm < 1.0)  { r_ohm *= 10.0; decade /= 10.0; }

    /* Binary search for nearest E96 mantissa */
    int lo = 0, hi = 95;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (e96_mantissas[mid] < r_ohm) lo = mid;
        else hi = mid;
    }

    double best = fabs(e96_mantissas[lo] - r_ohm) < fabs(e96_mantissas[hi] - r_ohm)
                  ? e96_mantissas[lo] : e96_mantissas[hi];
    return best * decade;
}

/* ==================================================================
 * L2: Select R_G from Standard E24 Resistor Values
 *
 * E24 series (5% tolerance): 24 values per decade.
 *
 * Complexity: O(1)
 * ================================================================== */
static const double e24_mantissas[24] = {
    1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0,
    2.2, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.3,
    4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1
};

double ia_e24_nearest(double r_ohm) {
    if (r_ohm <= 0.0) return 0.0;

    double decade = 1.0;
    while (r_ohm >= 10.0) { r_ohm /= 10.0; decade *= 10.0; }
    while (r_ohm < 1.0)  { r_ohm *= 10.0; decade /= 10.0; }

    double best_diff = INFINITY;
    double best_val = 1.0;

    for (int i = 0; i < 24; i++) {
        double diff = fabs(e24_mantissas[i] - r_ohm);
        if (diff < best_diff) {
            best_diff = diff;
            best_val = e24_mantissas[i];
        }
    }
    return best_val * decade;
}

/* ==================================================================
 * L2: Load Cell Sensitivity Converted to IA Input Range
 *
 * Load cell specification: 2 mV/V at full scale.
 * With V_exc = 5V, full-scale output = 10 mV.
 * For a 10 kg load cell, sensitivity = 1 mV/kg.
 *
 * This function computes the required IA gain to map the sensor
 * full-scale to the ADC input range.
 *
 *   G_required = V_adc_fs / (V_exc * sensitivity_mv_per_v)
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_gain_for_sensor_span(double v_sensor_fs, double v_adc_fs) {
    if (v_sensor_fs <= 0.0) return 0.0;
    return v_adc_fs / v_sensor_fs;
}

/* ==================================================================
 * L3: Bridge Sensor Differential Output Voltage
 *
 * Exact Wheatstone bridge equation:
 *   V_diff = V_exc * [R3/(R2+R3) - R4/(R1+R4)]
 *
 * For a quarter bridge with one active arm R1 = R + dR:
 *   V_diff = V_exc * [R/(R+R) - R/(R+dR+R)]
 *          = V_exc * dR / (4R + 2*dR)
 *
 * Linear approximation (dR << R):
 *   V_diff ~ V_exc * dR / (4R)
 *
 * Nonlinearity error for quarter bridge:
 *   err = (exact - approx) / exact = dR/(2R + dR) ~ dR/(2R)
 *   For dR/R = 0.01 (1% strain): err = 0.5%
 *
 * Reference: The Strain Gage Primer, Vishay Micro-Measurements
 * Complexity: O(1)
 * ================================================================== */
double bridge_output_voltage(const bridge_sensor_t *sensor, bool exact) {
    if (!sensor) return 0.0;

    double R = sensor->r_nominal;
    double dR = sensor->delta_r;
    double V_exc = sensor->v_excitation;

    if (exact) {
        /* Full exact formula for one-arm quarter bridge */
        if (R <= 0.0) return 0.0;
        return V_exc * dR / (4.0 * R + 2.0 * dR);
    } else {
        /* Linear approximation */
        if (R <= 0.0) return 0.0;
        return V_exc * dR / (4.0 * R);
    }
}

/* ==================================================================
 * L3: Bridge Nonlinearity Error Estimation
 *
 * For a quarter bridge:
 *   nonlinearity_pct = dR / (2R + dR) * 100
 *
 * For a half bridge (adjacent arms in compression/tension):
 *   V_diff = V_exc * dR / (2R)  (inherently linear!)
 *   nonlinearity = 0 (cancellation in bridge topology)
 *
 * Complexity: O(1)
 * ================================================================== */
double bridge_nonlinearity_pct(const bridge_sensor_t *sensor) {
    if (!sensor || sensor->r_nominal <= 0.0) return 0.0;

    if (sensor->active_arms >= 2) {
        return 0.0;  /* Half and full bridges are inherently linear */
    }

    double R = sensor->r_nominal;
    double dR = sensor->delta_r;
    return dR / (2.0 * R + dR) * 100.0;
}

/* ==================================================================
 * L1: Initialize IA Specification with Default Values
 *
 * Sets reasonable defaults for a general-purpose IA.
 * Complexity: O(1)
 * ================================================================== */
void ia_spec_init(ia_spec_t *spec) {
    if (!spec) return;
    memset(spec, 0, sizeof(ia_spec_t));
    spec->gain = 1.0;
    spec->gain_min = 1.0;
    spec->gain_max = 1000.0;
    spec->gain_error_pct = 0.1;
    spec->gain_drift_ppm = 10.0;
    spec->gain_nonlinearity = 10.0;
    spec->vos_uv = 50.0;
    spec->vos_drift_nv_per_c = 500.0;
    spec->ib_pa = 500.0;
    spec->ios_pa = 100.0;
    spec->bw_khz = 100.0;
    spec->sr_v_per_us = 0.5;
    spec->settling_time_us = 10.0;
    spec->gbp_mhz = 1.0;
    spec->cmrr_db = 90.0;
    spec->cmrr_at_50hz_db = 85.0;
    spec->psrr_db = 80.0;
    spec->en_rti_nv_per_rhz = 10.0;
    spec->en_0p1_10hz_uvpp = 0.5;
    spec->en_corner_hz = 10.0;
    spec->zin_diff_mohm = 10.0;
    spec->zin_cm_mohm = 10.0;
    spec->input_range_v = 10.0;
    spec->output_swing_v = 10.0;
    spec->supply_min_v = 2.7;
    spec->supply_max_v = 36.0;
    spec->iq_ua = 750.0;
}

/* ==================================================================
 * L1: Initialize State
 * Complexity: O(1)
 * ================================================================== */
void ia_state_init(ia_state_t *state) {
    if (!state) return;
    memset(state, 0, sizeof(ia_state_t));
    state->powered = false;
    state->calibrated = false;
}

/* ==================================================================
 * L1: Validate IA Specification
 *
 * Checks that all parameters are within physically reasonable bounds.
 * Returns IA_OK if valid, specific error code otherwise.
 * Complexity: O(1)
 * ================================================================== */
ia_error_code_t ia_spec_validate(const ia_spec_t *spec) {
    if (!spec) return IA_ERR_NULL_POINTER;

    if (spec->gain_min < 0.1 || spec->gain_max > 100000.0)
        return IA_ERR_INVALID_GAIN;
    if (spec->supply_min_v < 1.0 || spec->supply_max_v > 100.0)
        return IA_ERR_PARAMETER;
    if (spec->supply_min_v >= spec->supply_max_v)
        return IA_ERR_PARAMETER;
    if (spec->gain_min > spec->gain_max)
        return IA_ERR_INVALID_GAIN;

    return IA_OK;
}

/* ==================================================================
 * L2: Topology Recommendation Engine
 *
 * Selects best IA topology based on sensor type and supply voltage.
 *
 * Decision rules (from Kitchin & Counts, Ch. 5):
 *   - Strain gauge / load cell (low Z, bridge): 3-op-amp
 *   - Thermocouple (uV, high Z): flying-cap/chopper for low drift
 *   - Biopotential (AC, high Z): 3-op-amp with AC coupling
 *   - Low voltage < 3V: current-mode
 *   - High bandwidth needed: indirect current feedback
 *   - Hall effect (mV, low Z): 3-op-amp or 2-op-amp
 *   - Current shunt (uOhm-mOhm): 3-op-amp or current-mode
 *
 * Complexity: O(1)
 * ================================================================== */
ia_topology_t ia_recommend_topology(sensor_type_t sensor, double supply_voltage) {
    if (supply_voltage < 3.0) {
        return IA_TOPO_CURRENT_MODE;
    }

    switch (sensor) {
        case SENSOR_STRAIN_GAUGE:
        case SENSOR_LOAD_CELL:
        case SENSOR_PRESSURE_BRIDGE:
            return IA_TOPO_THREE_OPAMP;

        case SENSOR_THERMOCOUPLE:
            return IA_TOPO_FLYING_CAP;  /* nV drift needed */

        case SENSOR_RTD:
        case SENSOR_THERMISTOR:
            return IA_TOPO_THREE_OPAMP;

        case SENSOR_ECG_ELECTRODE:
        case SENSOR_EEG_ELECTRODE:
        case SENSOR_EMG_ELECTRODE:
            return IA_TOPO_THREE_OPAMP;

        case SENSOR_CURRENT_SHUNT:
            return IA_TOPO_CURRENT_MODE;

        case SENSOR_HALL_EFFECT:
        case SENSOR_ACCELEROMETER:
            return IA_TOPO_TWO_OPAMP;

        case SENSOR_PHOTODIODE:
            return IA_TOPO_INDIRECT_CURRENT;  /* BW needed */

        case SENSOR_LVDT:
            return IA_TOPO_THREE_OPAMP;  /* AC excitation */

        default:
            return IA_TOPO_THREE_OPAMP;
    }
}