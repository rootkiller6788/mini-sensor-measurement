/**
 * sigcond_selfcal.c - Self-Calibration and Auto-Zero Implementation
 *
 * Two-point calibration, auto-zero offset correction, gain error
 * calibration, system calibration matrix, and temperature drift
 * compensation for signal conditioning chains.
 *
 * L2: calibration concepts, offset/gain error models
 * L3: linear calibration models, matrix formulation
 * L5: two-point calibration, multi-point least-squares
 *
 * A signal conditioning channel's transfer function:
 *   V_out = G * (V_in + V_os) + V_os_out
 * where G is gain, V_os is input offset, V_os_out is output offset.
 *
 * Two-point calibration solves for G and V_os:
 *   V_out1 = G * (V_in1 + V_os), V_out2 = G * (V_in2 + V_os)
 *   G = (V_out2 - V_out1) / (V_in2 - V_in1)
 *   V_os = V_out1/G - V_in1
 *
 * Reference:
 *   Pallas-Areny & Webster (2001), Ch. 5
 *   Kester (2004), Ch. 7
 */

#include "sigcond_selfcal.h"
#include <stdlib.h>
#include <string.h>

/* ==================================================================
 * L5: Two-Point Calibration
 *
 * Uses two known reference inputs to determine gain and offset.
 * Most common calibration method for linear systems.
 *
 * Gain:     G = (Y2 - Y1) / (X2 - X1)
 * Offset:   Offset = Y1 - G * X1
 *
 * Result:   X_compensated = Y_measured * (1/G) - Offset/G
 *          or equivalently: X = (Y - Offset) / G
 *
 * Accuracy depends on:
 *   - Precision of reference inputs (X1, X2)
 *   - Repeatability of measurements
 *   - Linearity between calibration points
 *   - Temperature stability during calibration
 * ================================================================== */

/** Perform two-point calibration.
 *  Input:  two (known_input, measured_output) pairs
 *  Output: gain and offset for the linear model
 *  Returns 0 on success, -1 on error (identical inputs, null pointers)
 */
int selfcal_two_point(const calibration_point_t *p1,
                       const calibration_point_t *p2,
                       double *gain, double *offset)
{
    if (!p1 || !p2 || !gain || !offset) return -1;

    double dx = p2->measurand_value - p1->measurand_value;
    double dy = p2->conditioned_output - p1->conditioned_output;

    if (fabs(dx) < 1e-15) return -1;  /* Inputs too close */

    *gain = dy / dx;
    *offset = p1->conditioned_output - *gain * p1->measurand_value;

    return 0;
}

/** Apply calibration correction.
 *  corrected_value = (measured_output - offset) / gain
 */
double selfcal_apply_correction(double measured_output,
                                 double gain, double offset)
{
    if (fabs(gain) < 1e-15) return measured_output;
    return (measured_output - offset) / gain;
}

/** Reverse calibration: predict output for a given input */
double selfcal_predict_output(double input_value,
                               double gain, double offset)
{
    return gain * input_value + offset;
}

/* ==================================================================
 * L5: Auto-Zero (Offset Nulling)
 *
 * Auto-zero periodically shorts the input to measure and subtract
 * the offset voltage. Implementations:
 *   - Chopper: modulates input to AC, demodulates after amplification
 *   - Ping-pong: alternates between two amplifiers
 *   - Software: periodic short-input measurement + subtraction
 *
 * Software auto-zero algorithm:
 *   1. Short input (or connect to known reference)
 *   2. Measure output: V_zero = G * V_os_total
 *   3. For subsequent measurements: V_corrected = V_measured - V_zero
 *
 * Residual offset after auto-zero:
 *   V_os_residual = V_os_drift_rate * time_since_zero
 *   plus V_os_noise (1/f component at averaging frequency)
 * ================================================================== */

/** Perform auto-zero measurement.
 *  When input is grounded/shorted, measure the output voltage
 *  which represents the total offset referred to output.
 */
double selfcal_autozero_measure(double measured_zero_output)
{
    return measured_zero_output;  /* This is the offset voltage */
}

/** Apply auto-zero correction to subsequent measurements */
double selfcal_autozero_correct(double measured_output,
                                 double zero_offset)
{
    return measured_output - zero_offset;
}

/** Estimate residual offset drift since last auto-zero */
double selfcal_offset_drift_estimate(double last_zero_offset,
                                      double drift_rate_v_per_sec,
                                      double time_since_zero_sec)
{
    return last_zero_offset + drift_rate_v_per_sec * time_since_zero_sec;
}

/* ==================================================================
 * L5: Gain Error Calibration
 *
 * Gain error sources:
 *   - Resistor tolerance in gain-setting network
 *   - Finite open-loop gain of op-amps
 *   - Temperature coefficient of gain-setting resistors
 *   - Reference voltage drift (ADC)
 *
 * For in-amp with external R_G:
 *   G = 1 + (49.4 kOhm / R_G_kohm)  ... AD620 formula
 *   dG/G = -dR_G/R_G                ... gain inversely proportional
 *
 * For PGA with digital gain:
 *   G_error = G_nominal * (1 + gain_error_ppm/1e6)
 *
 * System-level calibration:
 *   Apply known precision input, measure output.
 *   G_actual = V_out / V_in
 *   Correction factor = G_nominal / G_actual
 * ================================================================== */

/** Calibrate gain using a known precision input */
double selfcal_gain_calibrate(double v_input_known,
                               double v_output_measured,
                               double gain_nominal)
{
    if (fabs(v_input_known) < 1e-15) return gain_nominal;
    return v_output_measured / v_input_known;
}

/** Compute gain correction factor */
double selfcal_gain_correction_factor(double gain_actual,
                                       double gain_nominal)
{
    if (fabs(gain_actual) < 1e-15) return 1.0;
    return gain_nominal / gain_actual;
}

/** Estimate gain temperature drift */
double selfcal_gain_drift(double gain_at_t0, double gain_drift_ppm_per_c,
                           double delta_temp_c)
{
    return gain_at_t0 * (1.0 + gain_drift_ppm_per_c * 1e-6 * delta_temp_c);
}

/* ==================================================================
 * L5: System Calibration Matrix
 *
 * For multi-channel systems, calibration accounts for:
 *   - Individual channel gain/offset
 *   - Channel-to-channel crosstalk
 *   - Reference voltage errors
 *   - ADC transfer function
 *
 * Calibration matrix M (NxN for N channels):
 *   Y_corrected = M * (Y_measured - b)
 *   where b is offset vector, M is inverse gain+crosstalk matrix.
 *
 * For independent channels (no crosstalk):
 *   M is diagonal: M_ii = 1/G_i
 *   b_i = offset_i
 * ================================================================== */

/** Multi-channel two-point calibration.
 *  Independent channel assumption: no crosstalk.
 *  channels: array of channel configs
 *  num_channels: number of channels
 *  cal_p1, cal_p2: calibration points at two reference levels
 *  gains: output array of calibrated gains [num_channels]
 *  offsets: output array of calibrated offsets [num_channels]
 */
void selfcal_multichannel_calibrate(unsigned num_channels,
                                     const sigcond_channel_config_t channels[],
                                     const double cal_p1_inputs[],
                                     const double cal_p1_outputs[],
                                     const double cal_p2_inputs[],
                                     const double cal_p2_outputs[],
                                     double gains[],
                                     double offsets[])
{
    if (!channels || !gains || !offsets || num_channels == 0) return;

    for (unsigned i = 0; i < num_channels; i++) {
        double dx = cal_p2_inputs[i] - cal_p1_inputs[i];
        double dy = cal_p2_outputs[i] - cal_p1_outputs[i];

        if (fabs(dx) < 1e-15) {
            gains[i] = 1.0;
            offsets[i] = 0.0;
        } else {
            gains[i] = dy / dx;
            offsets[i] = cal_p1_outputs[i] - gains[i] * cal_p1_inputs[i];
        }
    }
}

/** Temperature compensation using pre-characterized drift coefficients.
 *
 *  G(T) = G_0 * (1 + alpha_G * delta_T)
 *  V_os(T) = V_os,0 + alpha_Vos * delta_T + beta_Vos * delta_T^2
 *
 *  where alpha and beta are drift coefficients from characterization.
 */
double selfcal_temp_compensate_gain(double gain_at_t0,
                                     double alpha_gain_ppm_per_c,
                                     double delta_temp_c)
{
    return gain_at_t0 * (1.0 + alpha_gain_ppm_per_c * 1e-6 * delta_temp_c);
}

double selfcal_temp_compensate_offset(double offset_at_t0,
                                       double alpha_offset_uv_per_c,
                                       double beta_offset_uv_per_c2,
                                       double delta_temp_c)
{
    return offset_at_t0
         + alpha_offset_uv_per_c * delta_temp_c
         + beta_offset_uv_per_c2 * delta_temp_c * delta_temp_c;
}

/* ==================================================================
 * L5: Calibration Validity Checks
 * ================================================================== */

/** Check if calibration is still valid based on time and temperature */
bool selfcal_is_calibration_valid(double cal_temp_c,
                                   double current_temp_c,
                                   double cal_time_hours,
                                   double current_time_hours,
                                   double max_temp_drift_c,
                                   double max_time_hours)
{
    double temp_drift = fabs(current_temp_c - cal_temp_c);
    double time_elapsed = fabs(current_time_hours - cal_time_hours);

    return (temp_drift <= max_temp_drift_c) && (time_elapsed <= max_time_hours);
}

/** Compute total system uncertainty after calibration (RSS) */
double selfcal_total_uncertainty(double reference_uncertainty,
                                  double measurement_noise,
                                  double drift_uncertainty,
                                  double linearity_error)
{
    return sqrt(reference_uncertainty * reference_uncertainty
              + measurement_noise * measurement_noise
              + drift_uncertainty * drift_uncertainty
              + linearity_error * linearity_error);
}
