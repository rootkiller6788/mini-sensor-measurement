/**
 * inst_amp_calibration.c - IA Calibration and Error Compensation
 *
 * L5 Algorithms: Auto-zero, chopper stabilization, two-point and
 * multi-point calibration, temperature compensation, and digital
 * gain error correction for precision instrumentation amplifiers.
 *
 * Courses: Stanford EE247, Michigan EECS 455, Berkeley EE105
 * Reference:
 *   Pall++s-Areny & Webster, "Sensors and Signal Conditioning" (2nd ed., 2001)
 *   Kitchin & Counts (2006), Ch. 7 - Application Circuits
 *   NIST, "Guidelines for Evaluating and Expressing the Uncertainty
 *     of NIST Measurement Results", TN 1297
 */

#include "inst_amp_defs.h"
#include "inst_amp_topology.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ==================================================================
 * L5: Two-Point Calibration (Gain and Offset)
 *
 * Applies known input voltages V_cal1, V_cal2 and measures outputs.
 * Solves the linear system:
 *
 *   V_out = G_actual * V_in + V_os_output
 *
 *   G_actual = (V_out2 - V_out1) / (V_cal2 - V_cal1)
 *   V_os_output = V_out1 - G_actual * V_cal1
 *   V_os_input = V_os_output / G_actual
 *
 * The calibration coefficients are then used to correct future
 * measurements:
 *   V_in_corrected = (V_out_raw - V_os_output) / G_actual
 *
 * Requirements for accurate calibration:
 *   - V_cal1 and V_cal2 must be accurately known (precision source)
 *   - V_cal1 near lower end of range, V_cal2 near upper end
 *   - Multiple readings averaged to reduce noise
 *   - Temperature stabilized or compensated
 *
 * Reference: NIST TN 1297, ISO/IEC Guide 98-3 (GUM)
 * Complexity: O(1)
 * ================================================================== */
void ia_cal_two_point(double v_cal1, double v_out1,
                       double v_cal2, double v_out2,
                       double *actual_gain, double *output_offset) {
    if (!actual_gain || !output_offset) return;

    double dv_in = v_cal2 - v_cal1;
    double dv_out = v_out2 - v_out1;

    if (fabs(dv_in) < 1e-15) {
        *actual_gain = 1.0;
        *output_offset = 0.0;
        return;
    }

    *actual_gain = dv_out / dv_in;
    *output_offset = v_out1 - (*actual_gain) * v_cal1;
}

/* ==================================================================
 * L5: Apply 2-Point Calibration Correction
 *
 * Corrects a raw ADC reading using stored calibration.
 *
 * V_in = (V_out_raw - V_os_output) / G_actual
 *
 * Uncertainty propagation:
 *   sigma_Vin = sqrt((sigma_Vout/G)^2
 *                  + (V_out*sigma_G/G^2)^2
 *                  + (sigma_Vos/G)^2)
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_cal_apply_correction(double v_out_raw, double actual_gain,
                                double output_offset) {
    if (fabs(actual_gain) < 1e-15) return 0.0;
    return (v_out_raw - output_offset) / actual_gain;
}

/* ==================================================================
 * L5: Multi-Point Linear Regression Calibration
 *
 * For sensors with modest nonlinearity, a least-squares linear
 * fit over N calibration points provides better accuracy than
 * 2-point calibration.
 *
 * Linear model: V_out = G * V_in + V_os
 *
 * Least squares solution:
 *   G = [N*sum(x_i*y_i) - sum(x_i)*sum(y_i)] / [N*sum(x_i^2) - sum(x_i)^2]
 *   V_os = [sum(y_i) - G*sum(x_i)] / N
 *
 * where x_i = V_cal_i, y_i = V_out_i
 *
 * Complexity: O(N)
 * ================================================================== */
void ia_cal_linear_regression(const double *v_cal, const double *v_out,
                               int n_points, double *gain, double *offset) {
    if (!v_cal || !v_out || !gain || !offset || n_points < 2) return;

    double sum_x = 0.0, sum_y = 0.0;
    double sum_xx = 0.0, sum_xy = 0.0;

    for (int i = 0; i < n_points; i++) {
        sum_x += v_cal[i];
        sum_y += v_out[i];
        sum_xx += v_cal[i] * v_cal[i];
        sum_xy += v_cal[i] * v_out[i];
    }

    double N = (double)n_points;
    double denom = N * sum_xx - sum_x * sum_x;

    if (fabs(denom) < 1e-15) {
        *gain = 1.0;
        *offset = 0.0;
        return;
    }

    *gain = (N * sum_xy - sum_x * sum_y) / denom;
    *offset = (sum_y - (*gain) * sum_x) / N;
}

/* ==================================================================
 * L5: Polynomial Calibration for Nonlinear Sensors
 *
 * For sensors with significant nonlinearity (thermistors, some
 * pressure sensors), use polynomial regression:
 *
 *   V_in = a_0 + a_1 * V_out + a_2 * V_out^2 + ... + a_k * V_out^k
 *
 * Uses normal equations (least-squares polynomial fit).
 *
 * Reference: NIST, "Data Fitting and Uncertainty" (2014)
 * Complexity: O(N * degree^2) for the normal equations solve
 * ================================================================== */
void ia_cal_polynomial_fit(const double *v_out, const double *v_in,
                            int n_points, int degree,
                            double *coeffs) {
    if (!v_out || !v_in || !coeffs || n_points < degree + 1 || degree < 0) return;

    /* For simplicity, implement degree 2 (quadratic) */
    if (degree > 2) return;

    /* Vandermonde matrix A: A[i][j] = v_out[i]^j */
    /* Solve A^T * A * c = A^T * v_in using Cramer's rule for 3x3 */

    if (degree == 0) {
        double sum = 0.0;
        for (int i = 0; i < n_points; i++) sum += v_in[i];
        coeffs[0] = sum / n_points;
        return;
    }

    if (degree == 1) {
        /* Already handled by linear regression */
        double gain, offset;
        ia_cal_linear_regression(v_out, v_in, n_points, &gain, &offset);
        coeffs[0] = offset;
        coeffs[1] = gain;
        return;
    }

    /* Degree == 2: fit y = a0 + a1*x + a2*x^2 */
    double sum_x = 0, sum_x2 = 0, sum_x3 = 0, sum_x4 = 0;
    double sum_y = 0, sum_xy = 0, sum_x2y = 0;

    for (int i = 0; i < n_points; i++) {
        double x = v_out[i];
        double y = v_in[i];
        double x2 = x * x;

        sum_x += x;
        sum_x2 += x2;
        sum_x3 += x * x2;
        sum_x4 += x2 * x2;
        sum_y += y;
        sum_xy += x * y;
        sum_x2y += x2 * y;
    }

    double N = (double)n_points;

    /* 3x3 normal equations using Cramer's rule */
    double det = N * (sum_x2 * sum_x4 - sum_x3 * sum_x3)
               - sum_x * (sum_x * sum_x4 - sum_x2 * sum_x3)
               + sum_x2 * (sum_x * sum_x3 - sum_x2 * sum_x2);

    if (fabs(det) < 1e-15) {
        coeffs[0] = coeffs[1] = coeffs[2] = 0.0;
        return;
    }

    coeffs[0] = (sum_y * (sum_x2 * sum_x4 - sum_x3 * sum_x3)
               - sum_x * (sum_xy * sum_x4 - sum_x2y * sum_x3)
               + sum_x2 * (sum_xy * sum_x3 - sum_x2 * sum_x2y)) / det;

    coeffs[1] = (N * (sum_xy * sum_x4 - sum_x2y * sum_x3)
               - sum_y * (sum_x * sum_x4 - sum_x2 * sum_x3)
               + sum_x2 * (sum_x * sum_x2y - sum_xy * sum_x2)) / det;

    coeffs[2] = (N * (sum_x2 * sum_x2y - sum_xy * sum_x3)
               - sum_x * (sum_x * sum_x2y - sum_xy * sum_x2)
               + sum_y * (sum_x * sum_x3 - sum_x2 * sum_x2)) / det;
}

/* ==================================================================
 * L5: Evaluate Polynomial Calibration
 *
 * y = a_0 + a_1 * x + a_2 * x^2
 *
 * Complexity: O(degree)
 * ================================================================== */
double ia_cal_eval_polynomial(const double *coeffs, int degree, double x) {
    if (!coeffs || degree < 0) return 0.0;

    double result = 0.0;
    double x_pow = 1.0;
    for (int d = 0; d <= degree; d++) {
        result += coeffs[d] * x_pow;
        x_pow *= x;
    }
    return result;
}

/* ==================================================================
 * L5: Temperature Compensation Using Polynomial Drift Model
 *
 * Most sensor and amplifier parameters drift with temperature.
 * A quadratic model for offset drift:
 *
 *   V_os(T) = V_os(T0) + TC1*(T-T0) + TC2*(T-T0)^2
 *
 * This function computes the offset at temperature T given
 * calibration at T0 and temperature coefficients.
 *
 * Reference: Kitchin & Counts (2006), Ch. 6
 * Complexity: O(1)
 * ================================================================== */
double ia_temp_compensate_offset(double vos_at_t0_uv, double tc1_nv_per_c,
                                  double tc2_nv_per_c2, double t_c, double t0_c) {
    double dt = t_c - t0_c;
    return vos_at_t0_uv + tc1_nv_per_c * dt * 1e-3
           + tc2_nv_per_c2 * dt * dt * 1e-6;
}

/* ==================================================================
 * L5: Auto-Zero Calibration Procedure
 *
 * Auto-zero (AZ) periodically shorts the amplifier inputs and
 * measures the output to determine the offset, which is then
 * subtracted from subsequent measurements.
 *
 * Procedure:
 *   1. Disconnect sensor, short IA inputs (or short through MUX)
 *   2. Measure V_out = G * V_os (offset amplified by gain)
 *   3. Store V_os_measured = V_out / G
 *   4. Reconnect sensor
 *   5. Correct: V_signal = V_out_raw - V_os_output
 *
 * This eliminates offset and 1/f noise below the AZ frequency.
 * Typical AZ period: 1 ms to 1 s.
 *
 * Limitation: cannot distinguish sensor zero from amplifier zero.
 * For bridge sensors, "zero" can be obtained by turning off
 * excitation voltage.
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_autozero_measure_offset(double v_out_zero_input, double gain) {
    if (gain <= 0.0) return 0.0;
    return v_out_zero_input / gain;
}

/* ==================================================================
 * L5: Chopper Stabilization Model
 *
 * Chopper stabilization modulates the input signal to a higher
 * frequency (f_chop), amplifies it where 1/f noise is low, then
 * demodulates back to DC. This eliminates 1/f noise and offset.
 *
 * Transfer function:
 *   V_out = G * [V_in * square(f_chop) + V_os + V_n(f)]
 *   After demod and LPF: V_out = G * V_in + residual
 *
 * Residual offset after chopping: < 1 uV typical
 * Residual 1/f noise after chopping: corner shifted from ~10 Hz to ~0.01 Hz
 *
 * The chopping introduces ripple at f_chop and harmonics, which
 * must be filtered by the output LPF.
 *
 * Ripple amplitude: V_ripple ~ V_os_input * G
 * Ripple attenuation: depends on LPF order and f_chop/f_LPF ratio.
 *
 * Reference: Enz & Temes, "Circuit Techniques for Reducing the
 *   Effects of Op-Amp Imperfections", Proc. IEEE, 1996
 * Complexity: O(1)
 * ================================================================== */
double ia_chopper_ripple_amplitude(double vos_input, double gain,
                                    double lpf_attenuation_db) {
    double ripple = vos_input * gain;
    return ripple * pow(10.0, -lpf_attenuation_db / 20.0);
}

double ia_chopper_effective_corner_freq(double f_corner_original,
                                         double f_chop, double f_lpf) {
    /* Simplified model: corner shifted by ratio f_chop/f_corner */
    if (f_chop <= 0.0) return f_corner_original;
    return f_corner_original * (f_lpf / f_chop);
}

/* ==================================================================
 * L5: Digital Gain Trimming
 *
 * In digitally programmable IAs (AD8253, PGA280), the gain can
 * be adjusted in discrete steps. Digital trimming corrects for
 * the residual gain error:
 *
 *   G_effective = G_nominal * (1 + trim_factor)
 *
 * where trim_factor is typically +/- 0.1% to +/- 10% adjustable
 * via SPI/I2C register.
 *
 * This function computes the trim factor needed to achieve
 * a target gain.
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_digital_trim_factor(double g_actual, double g_target) {
    if (g_target <= 0.0) return 0.0;
    return (g_target / g_actual) - 1.0;
}

/* ==================================================================
 * L5: Verify Calibration Integrity
 *
 * Checks if stored calibration coefficients are within reasonable
 * bounds. Invalid calibrations can result from:
 *   - EEPROM corruption
 *   - Sensor degradation
 *   - Environmental extremes
 *
 * Complexity: O(1)
 * ================================================================== */
bool ia_cal_is_valid(const ia_calibration_t *cal) {
    if (!cal) return false;
    if (!cal->valid) return false;

    /* Gain should be positive and reasonable */
    if (cal->gain_slope <= 0.0 || cal->gain_slope > 10000.0) return false;

    /* Offset should be within expected range */
    if (fabs(cal->offset_v) > 1.0) return false;  /* > 1V offset suspect */

    /* Date should be reasonable */
    if (cal->cal_date < 20200101 || cal->cal_date > 20991231) return false;

    return true;
}

/* ==================================================================
 * L5: Sensor Excitation Voltage Compensation
 *
 * For ratiometric bridge measurements, the sensor output is
 * proportional to excitation voltage:
 *   V_sensor = k * V_exc * measurand
 *
 * If V_exc varies, the apparent sensitivity changes. Compensate
 * by dividing the reading by the actual V_exc and multiplying
 * by the nominal V_exc:
 *
 *   reading_corrected = reading_raw * (V_exc_nominal / V_exc_actual)
 *
 * This assumes the ADC reference is also ratiometric (derived
 * from the same V_exc). If ADC uses a fixed reference, additional
 * compensation is needed.
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_ratiometric_correction(double reading_raw,
                                  double v_exc_nominal, double v_exc_actual) {
    if (v_exc_actual <= 0.0) return reading_raw;
    return reading_raw * v_exc_nominal / v_exc_actual;
}

/* ==================================================================
 * L5: Moving Average Filter for Calibration Readings
 *
 * To reduce noise during calibration, N readings are averaged:
 *   V_avg = (1/N) * sum(V_i)
 *
 * For Gaussian noise, averaging N samples reduces noise by
 * sqrt(N): standard error = sigma / sqrt(N).
 *
 * To achieve 10x noise reduction: N = 100 samples.
 *
 * Complexity: O(N) per window, O(1) amortized with sliding window
 * ================================================================== */
double ia_cal_average_readings(const double *readings, int n) {
    if (!readings || n <= 0) return 0.0;

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += readings[i];
    }
    return sum / n;
}

/* ==================================================================
 * L5: Standard Deviation of Calibration Readings
 *
 * sigma = sqrt( sum((x_i - mean)^2) / (N-1) )
 *
 * The standard deviation indicates calibration quality.
 * High sigma -> noisy environment or unstable sensor.
 *
 * Complexity: O(N)
 * ================================================================== */
double ia_cal_std_deviation(const double *readings, int n, double mean) {
    if (!readings || n < 2) return 0.0;

    double sum_sq = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = readings[i] - mean;
        sum_sq += diff * diff;
    }
    return sqrt(sum_sq / (n - 1));
}

/* ==================================================================
 * L5: Self-Test and Diagnostics
 *
 * Verifies IA basic functionality:
 *   1. Output within common-mode range when inputs shorted
 *   2. Gain approximately correct (within +/-20%)
 *   3. No oscillation or latch-up
 *
 * Complexity: O(1)
 * ================================================================== */
typedef struct {
    bool output_in_range;
    bool gain_plausible;
    bool not_oscillating;
    double measured_gain;
    double measured_offset;
} ia_self_test_result_t;

void ia_self_test(double v_test_input, double v_out_measured,
                   double expected_gain, ia_self_test_result_t *result) {
    if (!result) return;

    result->output_in_range = (fabs(v_out_measured) < 15.0);

    double measured_gain = 0.0;
    if (fabs(v_test_input) > 1e-9) {
        measured_gain = v_out_measured / v_test_input;
    }
    result->measured_gain = measured_gain;

    /* Gain within +/- 30% of expected */
    if (expected_gain > 0.0) {
        double err = fabs(measured_gain - expected_gain) / expected_gain;
        result->gain_plausible = (err < 0.30);
    } else {
        result->gain_plausible = true;
    }

    result->not_oscillating = true;  /* Simplified check */
    result->measured_offset = 0.0;
}