#ifndef SENSOR_CONDITIONING_H
#define SENSOR_CONDITIONING_H
#include "smart_sensor.h"
#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * Sensor Signal Conditioning Module
 *
 * Implements analog and digital signal conditioning algorithms:
 *   L3: Wheatstone bridge analysis, transfer function evaluation,
 *       noise analysis, effective number of bits (ENOB)
 *   L5: Digital filtering (moving average, median filter, IIR),
 *       linearization, temperature compensation, peak detection
 *   L6: Thermistor linearization, strain gauge bridge conditioning
 *
 * Reference:
 *   Sedra, A.S. & Smith, K.C. (2020), "Microelectronic Circuits", 8th ed., Oxford
 *   Kester, W. (2005), "Data Conversion Handbook", Analog Devices/Newnes
 *   Pallas-Areny, R. & Webster, J.G. (2001), "Sensors and Signal Conditioning",
 *     2nd ed., Wiley
 *
 * Course Mapping:
 *   Berkeley EE105 — Analog Circuits (Wheatstone bridge, amplifiers)
 *   MIT 6.003 — Signal Processing (filtering)
 *   Stanford EE102A — Signal Processing (digital filter design)
 * ========================================================================= */

/* =========================================================================
 * L3: Wheatstone Bridge Analysis
 *
 * The Wheatstone bridge is the most common signal conditioning
 * circuit for resistive sensors (strain gauges, RTDs, thermistors,
 * pressure sensors).
 *
 * Bridge output voltage (quarter-bridge):
 *   V_out = V_ex * (Delta_R / (4R_0 + 2*Delta_R))
 *          ≈ V_ex * Delta_R / (4R_0)   (for small Delta_R << R_0)
 *
 * For a full bridge with equal initial resistances:
 *   V_out = V_ex * Delta_R / R_0
 * ========================================================================= */

/**
 * @brief Compute Wheatstone bridge output voltage
 *
 * Quarter-bridge configuration (one active arm, three fixed resistors):
 *   V_out = V_ex * (R1/(R1+R2) - R3/(R3+R4))
 *
 * @param v_ex   Bridge excitation voltage (V)
 * @param r1     Resistance of active gauge (Ohm) = R0 + Delta_R
 * @param r2     Fixed resistor R2 (Ohm)
 * @param r3     Fixed resistor R3 (Ohm)
 * @param r4     Fixed resistor R4 (Ohm)
 * @return       Bridge output voltage (V)
 */
double ss_cond_wheatstone_voltage(double v_ex, double r1, double r2,
                                  double r3, double r4);

/**
 * @brief Wheatstone bridge linearized strain approximation
 *
 * For small strains (Delta_R << R0), uses the linear approximation:
 *   strain = (4 / GF) * (V_out / V_ex)   (quarter-bridge)
 *   strain = (2 / GF) * (V_out / V_ex)   (half-bridge)
 *   strain = (1 / GF) * (V_out / V_ex)   (full-bridge)
 *
 * @param v_out      Measured bridge output voltage (V)
 * @param v_ex       Bridge excitation voltage (V)
 * @param gf         Gauge factor (typically 2.0 for metal foil gauges)
 * @param bridge_cfg 0=quarter, 1=half, 2=full bridge
 * @return           Estimated strain (unitless, microstrain*1e-6)
 */
double ss_cond_strain_from_bridge(double v_out, double v_ex,
                                  double gf, int bridge_cfg);

/**
 * @brief Compute required amplifier gain for bridge signal
 *
 * Determines amplifier gain needed to map sensor full-scale output
 * to ADC full-scale input range:
 *   gain = adc_vref / (bridge_vout_at_full_scale)
 *
 * @param bridge_fs_output  Bridge output at full-scale sensor input (V)
 * @param adc_vref          ADC reference voltage (V)
 * @return                  Required amplifier gain (V/V)
 */
double ss_cond_amplifier_gain_needed(double bridge_fs_output,
                                     double adc_vref);

/* =========================================================================
 * L3: ADC Effective Number of Bits (ENOB) and Noise Analysis
 *
 * ENOB quantifies the real-world resolution of an ADC accounting
 * for noise and distortion:
 *   ENOB = (SINAD_dB - 1.76) / 6.02
 *
 * SINAD = Signal-to-Noise-and-Distortion ratio in dB.
 *
 * Reference: IEEE Std 1241-2010, "Standard for Terminology and Test
 *            Methods for Analog-to-Digital Converters"
 * ========================================================================= */

/**
 * @brief Compute ADC effective number of bits (ENOB)
 *
 * @param sinad_db  Signal-to-Noise-and-Distortion ratio in dB
 * @return          ENOB (fractional bits)
 */
double ss_cond_adc_enob(double sinad_db);

/**
 * @brief Compute ADC quantization noise RMS
 *
 * For an ideal N-bit ADC with full-scale range V_FS:
 *   sigma_q = V_FS / (2^N * sqrt(12))
 *
 * @param resolution_bits  ADC resolution (bits)
 * @param vref             ADC reference voltage (V)
 * @return                 RMS quantization noise (V)
 */
double ss_cond_adc_quantization_noise(uint8_t resolution_bits, double vref);

/**
 * @brief Compute ADC signal-to-quantization-noise ratio (SQNR)
 *
 * For a full-scale sine wave input:
 *   SQNR_dB = 6.02 * N + 1.76
 *
 * @param resolution_bits  ADC resolution (bits)
 * @return                 SQNR in dB
 */
double ss_cond_adc_sqnr(uint8_t resolution_bits);

/**
 * @brief Compute total sensor chain noise (RTI)
 *
 * Root-sum-square of all referred-to-input noise sources:
 *   total_noise = sqrt(sensor_noise^2 + amp_noise^2 + adc_noise^2)
 *
 * @param sensor_noise_rti   Sensor noise referred to input
 * @param amplifier_noise_rti Amplifier input-referred noise
 * @param adc_noise_rti      ADC noise referred to input
 * @return                   Total RMS noise RTI
 */
double ss_cond_total_noise_rti(double sensor_noise_rti,
                               double amplifier_noise_rti,
                               double adc_noise_rti);

/* =========================================================================
 * L5: Digital Filtering for Sensor Signals
 * ========================================================================= */

/**
 * @brief Simple moving average filter
 *
 * Computes the arithmetic mean of the last N samples.
 * Effective for reducing uncorrelated zero-mean noise.
 * Noise reduction factor: sqrt(N) (for white noise).
 *
 * Complexity: O(N) per call (direct implementation)
 *
 * @param buffer      Sliding window of last N samples (updated, oldest shifted out)
 * @param window_size Number of samples in the averaging window
 * @param new_sample  New sample to incorporate
 * @return            Filtered (averaged) output
 */
double ss_cond_moving_average(double *buffer, size_t window_size,
                              double new_sample);

/**
 * @brief Running moving average (O(1) per sample)
 *
 * Maintains running sum to avoid re-summing the window each time.
 * When the window is not yet full, computes partial average.
 *
 * @param window       Circular buffer of last N samples [window_size]
 * @param window_size  Max window size
 * @param count        Sample count up to window_size (updated in place)
 * @param idx          Current insert index in circular buffer (updated in place)
 * @param running_sum  Running sum accumulator (updated in place)
 * @param new_sample   New sample
 * @return             Running average
 */
double ss_cond_running_moving_avg(double *window, size_t window_size,
                                  size_t *count, size_t *idx,
                                  double *running_sum, double new_sample);

/**
 * @brief Median filter for impulse noise removal
 *
 * Sorts a small window (3, 5, or 7 samples) and returns the middle value.
 * Excellent for removing salt-and-pepper noise and single-sample spikes
 * while preserving edges.
 *
 * Complexity: O(N log N) for sorting, but N is small (3-7).
 * The window_size must be odd (3, 5, or 7).
 *
 * @param buffer      Sliding window of samples [window_size] (updated in place)
 * @param window_size  Must be odd: 3, 5, or 7
 * @param new_sample  New sample to incorporate
 * @return            Median value of the current window
 */
double ss_cond_median_filter(double *buffer, size_t window_size,
                             double new_sample);

/**
 * @brief First-order IIR lowpass filter
 *
 * Transfer function:
 *   H(z) = alpha / (1 - (1-alpha)*z^{-1})
 *
 * Time-domain recursion:
 *   y[n] = alpha * x[n] + (1-alpha) * y[n-1]
 *
 * -3 dB cutoff frequency for a given alpha and sample rate fs:
 *   f_c ≈ alpha * fs / (2*pi)   (for alpha << 1)
 *   alpha = 2*pi*f_c / fs
 *
 * @param prev_output  Previous filter output (updated in place)
 * @param input        Current input sample
 * @param alpha        Filter coefficient (0-1)
 *                     0 < alpha <= 1: alpha=1 = no filtering (output=input)
 *                     alpha small = more filtering, slower response
 * @return             Filtered output
 */
double ss_cond_iir_lowpass(double *prev_output, double input, double alpha);

/**
 * @brief First-order IIR highpass filter
 *
 * Transfer function:
 *   H(z) = (1-alpha) * (1 - z^{-1}) / (1 - (1-alpha)*z^{-1})
 *
 * Time-domain recursion:
 *   y[n] = (1-alpha) * (y[n-1] + x[n] - x[n-1])
 *
 * -3 dB cutoff frequency:
 *   alpha = 2*pi*f_c / fs
 *
 * @param prev_output  Previous filter output (updated in place)
 * @param prev_input   Previous input sample (updated in place)
 * @param input        Current input sample
 * @param alpha        Filter coefficient (0-1)
 * @return             Filtered output
 */
double ss_cond_iir_highpass(double *prev_output, double *prev_input,
                            double input, double alpha);

/**
 * @brief Peak detector with exponential decay
 *
 * Tracks the maximum value with controlled decay. Useful for
 * envelope detection in vibration/impact sensors.
 *
 *   if (x[n] > peak): peak = x[n]
 *   else:             peak = peak * decay_factor
 *
 * @param peak          Current peak value (updated in place)
 * @param input         Current input sample
 * @param decay_factor  Per-sample decay multiplier (0-1, e.g., 0.999)
 * @return              Updated peak value
 */
double ss_cond_peak_detector(double *peak, double input, double decay_factor);

/**
 * @brief Rate-of-change (derivative) estimator
 *
 * Simple backward-difference derivative approximation:
 *   dx/dt ≈ (x[n] - x[n-1]) / dt
 *
 * @param prev_value  Previous sample value (updated in place)
 * @param curr_value  Current sample value
 * @param dt          Time between samples (seconds)
 * @return            Estimated rate of change (units/sec)
 */
double ss_cond_rate_of_change(double *prev_value, double curr_value, double dt);

/**
 * @brief Threshold detector with hysteresis (Schmitt trigger)
 *
 * Prevents oscillation when signal crosses threshold repeatedly due
 * to noise. Uses separate rising and falling thresholds.
 *
 * State machine:
 *   Below lower threshold: state = 0 (inactive)
 *   Above upper threshold: state = 1 (active)
 *   Between thresholds:    state maintains previous value (hysteresis)
 *
 * @param state         Current output state (0 or 1, updated in place)
 * @param input         Current signal value
 * @param low_thresh    Low threshold (to transition 1->0)
 * @param high_thresh   High threshold (to transition 0->1)
 * @return              Current state (0 or 1)
 */
int ss_cond_threshold_hysteresis(int *state, double input,
                                 double low_thresh, double high_thresh);

/* =========================================================================
 * L5: Sensor Linearization and Compensation
 * ========================================================================= */

/**
 * @brief Piecewise linearization using a lookup table
 *
 * Given a monotonic calibration table, uses linear interpolation
 * to convert a raw sensor reading to engineering units.
 *
 * @param x_table   Reference (known input) values [n] — sorted ascending
 * @param y_table   Raw sensor output values [n]
 * @param n         Number of table entries (>=2)
 * @param y_input   Raw sensor reading to linearize
 * @return          Linearized engineering units value
 */
double ss_cond_piecewise_linearize(const double *x_table,
                                   const double *y_table,
                                   size_t n, double y_input);

/**
 * @brief Temperature compensation: zero and span correction
 *
 * Adjusts sensor reading for ambient temperature effects:
 *   reading_corrected = (reading_raw - zero_shift) / span_factor
 *
 * where:
 *   zero_shift = TC_zero * (T - T_ref)
 *   span_factor = 1 + TC_span * (T - T_ref)
 *
 * @param raw_reading      Uncompensated sensor reading
 * @param current_temp     Current ambient temperature (degC)
 * @param ref_temp         Reference calibration temperature (degC)
 * @param tc_zero          Zero drift coefficient (units/degC)
 * @param tc_span          Span drift coefficient (1/degC)
 * @return                 Temperature-compensated reading
 */
double ss_cond_temp_compensate(double raw_reading, double current_temp,
                               double ref_temp, double tc_zero,
                               double tc_span);

/**
 * @brief Non-linearity correction using polynomial
 *
 * Applies an inverse polynomial correction to linearize sensor output.
 * If the sensor forward model is y = f(x) where f is nonlinear,
 * this applies the inverse: x_estimate = g(y) where g ≈ f^{-1}.
 *
 * @param y_raw     Raw (nonlinear) sensor output
 * @param coeffs    Inverse polynomial coefficients [order+1], highest degree first
 * @param order     Polynomial order
 * @return          Linearized estimate
 */
double ss_cond_poly_nonlinearity_correct(double y_raw,
                                         const double *coeffs,
                                         uint8_t order);

/* =========================================================================
 * L3: Sensor Transfer Function Analysis
 * ========================================================================= */

/**
 * @brief Evaluate a generic sensor transfer function
 *
 * For sensors characterized by a mathematical transfer function,
 * this evaluates the output given an input.
 *
 * Supports: linear, polynomial, exponential, logarithmic, and
 * power-law sensor models.
 *
 * Model types: 0=linear, 1=exp, 2=log, 3=power, 4=sigmoid
 *
 * @param input       Physical input value
 * @param params      Model parameters [4]: [a, b, c, d]
 * @param model_type  0=ax+b, 1=a*exp(bx)+c, 2=a*log(bx)+c, 3=a*x^b+c, 4=a/(1+exp(-b(x-c)))+d
 * @return            Sensor output
 */
double ss_cond_transfer_function(double input, const double params[4],
                                 int model_type);

/**
 * @brief Sensor sensitivity analysis at operating point
 *
 * Computes the local sensitivity (derivative) of the sensor
 * transfer function at a given operating point using finite
 * differences.
 *
 *   sensitivity(x0) ≈ [f(x0+h) - f(x0-h)] / (2h)
 *
 * @param input       Operating point
 * @param params      Transfer function parameters [4]
 * @param model_type  Transfer function type (0-4)
 * @param h           Step size for finite difference
 * @return            Local sensitivity df/dx at operating point
 */
double ss_cond_sensitivity_analysis(double input, const double params[4],
                                    int model_type, double h);

/**
 * @brief Compute sensor dynamic range in dB
 *
 *   DR_dB = 20 * log10(full_scale / resolution)
 *
 * @param full_scale  Full-scale input range
 * @param resolution  Smallest detectable input
 * @return            Dynamic range in dB
 */
double ss_cond_dynamic_range_db(double full_scale, double resolution);

/* =========================================================================
 * L4: Signal-to-Noise Ratio Analysis
 * ========================================================================= */

/**
 * @brief Compute SNR from signal and noise powers
 *
 *   SNR_dB = 10 * log10(P_signal / P_noise)
 *
 * @param signal_rms  RMS signal amplitude
 * @param noise_rms   RMS noise amplitude
 * @return            SNR in dB
 */
double ss_cond_snr_db(double signal_rms, double noise_rms);

/**
 * @brief Compute minimum detectable signal (MDS)
 *
 * The smallest signal that can be reliably detected,
 * typically defined at SNR = 0 dB (signal power = noise power):
 *   MDS = noise_rms
 *
 * Or at a specified SNR threshold:
 *   MDS = noise_rms * 10^(SNR_min_dB / 20)
 *
 * @param noise_rms     RMS noise amplitude
 * @param snr_min_db    Minimum required SNR for detection (dB)
 * @return              Minimum detectable signal
 */
double ss_cond_min_detectable_signal(double noise_rms, double snr_min_db);

#ifdef __cplusplus
}
#endif
#endif /* SENSOR_CONDITIONING_H */
