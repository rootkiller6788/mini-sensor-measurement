/**
 * @file    sensor_conditioning.c
 * @brief   Sensor signal conditioning: analog and digital processing
 *
 * Knowledge Coverage:
 *   L3 (Math Structures): Wheatstone bridge analysis, transfer function
 *                         models, noise analysis, ADC performance metrics
 *   L5 (Algorithms): Digital filtering (moving average, median, IIR),
 *                    peak detection, threshold hysteresis, temperature
 *                    compensation, piecewise linearization
 *   L4 (Fundamental Laws): Bridge balance equation, ENOB formula,
 *                          quantization noise, Nyquist sampling
 *
 * Reference:
 *   Sedra & Smith (2020), "Microelectronic Circuits", 8th ed., Oxford
 *   Kester, W. (2005), "Data Conversion Handbook", Analog Devices
 *   Pallas-Areny & Webster (2001), "Sensors and Signal Conditioning", 2nd ed., Wiley
 *   IEEE Std 1241-2010 — Terminology and Test Methods for ADCs
 *   Oppenheim & Schafer (2010), "Discrete-Time Signal Processing", 3rd ed., Pearson
 *
 * Course Mapping:
 *   Berkeley EE105 — Analog Circuits (bridge circuits, amplifiers)
 *   MIT 6.003 — Signals and Systems (filtering, noise analysis)
 *   Stanford EE102A — Signal Processing (digital filter design)
 *   Illinois ECE 310 — DSP (digital filtering algorithms)
 */

#include "smart_sensor.h"
#include "sensor_conditioning.h"
#include <string.h>
#include <math.h>
#include <stdio.h>

/* =========================================================================
 * L3: Wheatstone Bridge Analysis
 * ========================================================================= */

/**
 * @brief Compute Wheatstone bridge differential output voltage
 *
 * The Wheatstone bridge is a four-arm bridge circuit used to measure
 * small resistance changes with high precision. The output voltage
 * is the difference between two voltage dividers:
 *
 *   V_out = V_ex * (R1/(R1+R2) - R3/(R3+R4))
 *
 * When balanced (R1*R4 = R2*R3), V_out = 0.
 * A change in any resistor unbalances the bridge, producing an output
 * proportional to the change.
 *
 * Applications: strain gauges, RTDs, pressure sensors, load cells.
 *
 * Reference: Wheatstone, C. (1843), Phil. Trans. R. Soc. Lond. 133:303-327
 *
 * @param v_ex  Bridge excitation voltage (V)
 * @param r1    Top-left resistor (Ohm)
 * @param r2    Bottom-left resistor (Ohm)
 * @param r3    Top-right resistor (Ohm)
 * @param r4    Bottom-right resistor (Ohm)
 * @return      Differential output voltage (V)
 */
double ss_cond_wheatstone_voltage(double v_ex, double r1, double r2,
                                  double r3, double r4) {
    double v_left, v_right;
    if (r1 + r2 == 0.0 || r3 + r4 == 0.0) return 0.0;
    v_left  = v_ex * r2 / (r1 + r2);
    v_right = v_ex * r4 / (r3 + r4);
    return v_left - v_right;
}

/**
 * @brief Estimate strain from bridge output using linearized formula
 *
 * For a quarter-bridge with one active gauge:
 *   V_out ≈ V_ex * GF * strain / 4
 *   strain ≈ 4 * V_out / (V_ex * GF)
 *
 * For a half-bridge (two active gauges in adjacent arms):
 *   V_out ≈ V_ex * GF * strain / 2
 *   strain ≈ 2 * V_out / (V_ex * GF)
 *
 * For a full-bridge (four active gauges):
 *   V_out ≈ V_ex * GF * strain
 *   strain ≈ V_out / (V_ex * GF)
 *
 * Gauge factor (GF) = (Delta_R/R) / strain, typically ~2.0 for metals.
 *
 * @param v_out      Measured bridge output voltage (V)
 * @param v_ex       Excitation voltage (V)
 * @param gf         Gauge factor
 * @param bridge_cfg 0=quarter, 1=half, 2=full bridge
 * @return           Estimated strain (epsilon, dimensionless, i.e., m/m)
 */
double ss_cond_strain_from_bridge(double v_out, double v_ex,
                                  double gf, int bridge_cfg) {
    double factor;
    if (v_ex == 0.0 || gf == 0.0) return 0.0;

    switch (bridge_cfg) {
    case 0: factor = 4.0; break;  /* Quarter bridge */
    case 1: factor = 2.0; break;  /* Half bridge */
    case 2: factor = 1.0; break;  /* Full bridge */
    default: factor = 4.0; break;
    }

    return factor * v_out / (v_ex * gf);
}

/**
 * @brief Compute required amplifier gain for ADC full-scale matching
 *
 * Determines the gain needed to map the bridge's full-scale output
 * to the ADC's full-scale input range:
 *
 *   gain_required = ADC_Vref / Bridge_Vout_at_FS
 *
 * This ensures optimal use of ADC dynamic range.
 *
 * @param bridge_fs_output  Bridge output at sensor full-scale (V)
 * @param adc_vref          ADC reference voltage (V)
 * @return                  Required amplifier gain (V/V)
 */
double ss_cond_amplifier_gain_needed(double bridge_fs_output,
                                     double adc_vref) {
    if (bridge_fs_output == 0.0 || adc_vref == 0.0) return 1.0;
    return adc_vref / bridge_fs_output;
}

/* =========================================================================
 * L3: ADC Performance Metrics
 * ========================================================================= */

/**
 * @brief Compute ADC Effective Number of Bits (ENOB)
 *
 * ENOB quantifies the real-world dynamic performance of an ADC,
 * accounting for noise and distortion beyond ideal quantization.
 *
 *   ENOB = (SINAD_dB - 1.76) / 6.02
 *
 * Where 1.76 dB comes from the theoretical SNR of an ideal N-bit ADC
 * with a full-scale sinusoidal input:
 *   SNR_ideal = 6.02*N + 1.76 dB
 *
 * Reference: IEEE Std 1241-2010, Section 4.5.1
 *
 * A 16-bit ADC with SINAD = 84 dB:
 *   ENOB = (84 - 1.76) / 6.02 = 13.66 bits (effective bits)
 *
 * @param sinad_db  Signal-to-Noise-and-Distortion ratio (dB)
 * @return          ENOB (fractional number of bits)
 */
double ss_cond_adc_enob(double sinad_db) {
    return (sinad_db - 1.76) / 6.02;
}

/**
 * @brief Compute ADC RMS quantization noise
 *
 * For an ideal ADC with uniform quantization error distribution
 * over [-LSB/2, +LSB/2], the RMS quantization noise is:
 *
 *   sigma_q = LSB / sqrt(12) = V_FS / (2^N * sqrt(12))
 *
 * This assumes no missing codes, no DNL errors, and uniform
 * distribution of quantization errors (Bennett's approximation).
 *
 * Reference: Bennett, W.R. (1948), Bell System Tech. J. 27:446-472
 *
 * @param resolution_bits  ADC resolution (bits)
 * @param vref             Reference voltage (V), equals full-scale range
 * @return                 RMS quantization noise (V)
 */
double ss_cond_adc_quantization_noise(uint8_t resolution_bits, double vref) {
    if (resolution_bits == 0 || vref <= 0.0) return 0.0;
    return vref / (pow(2.0, (double)resolution_bits) * sqrt(12.0));
}

/**
 * @brief Compute ideal ADC SQNR (Signal-to-Quantization-Noise Ratio)
 *
 * For a full-scale sinusoidal input:
 *   SQNR_dB = 6.02 * N + 1.76
 *
 * This is the theoretical upper bound for SNR of an ideal N-bit ADC.
 * Real ADCs achieve less due to thermal noise, jitter, DNL, and distortion.
 *
 * Derivation:
 *   Signal RMS: V_FS / (2*sqrt(2))
 *   Noise RMS:  V_FS / (2^N * sqrt(12))
 *   SQNR = (Signal RMS / Noise RMS)^2 = (3/2) * 2^(2N)
 *   SQNR_dB = 10*log10(3/2) + 20*N*log10(2) = 1.76 + 6.02*N
 *
 * @param resolution_bits  ADC resolution (bits)
 * @return                 SQNR in dB
 */
double ss_cond_adc_sqnr(uint8_t resolution_bits) {
    return 6.02 * (double)resolution_bits + 1.76;
}

/**
 * @brief Compute total sensor chain noise (referred to input)
 *
 * Combines multiple independent noise sources using root-sum-of-squares
 * (assuming uncorrelated Gaussian noise sources):
 *
 *   total_noise_RTI = sqrt(sensor_noise^2 + amp_noise^2 + adc_noise^2)
 *
 * All noise values must be referred to the same point — typically
 * the sensor input — for meaningful combination.
 *
 * @param sensor_noise_rti     Sensor element noise (units RTI)
 * @param amplifier_noise_rti  Amplifier input-referred noise
 * @param adc_noise_rti        ADC noise referred to sensor input
 * @return                     Total RMS noise RTI
 */
double ss_cond_total_noise_rti(double sensor_noise_rti,
                               double amplifier_noise_rti,
                               double adc_noise_rti) {
    return sqrt(sensor_noise_rti * sensor_noise_rti
              + amplifier_noise_rti * amplifier_noise_rti
              + adc_noise_rti * adc_noise_rti);
}

/* =========================================================================
 * L5: Digital Filtering — Moving Average
 * ========================================================================= */

/**
 * @brief Simple moving average (SMA) filter
 *
 * Computes the mean of the last N samples by shifting a window buffer
 * and summing the contents. The SMA filter has:
 *   - Linear phase response (symmetric — no phase distortion)
 *   - Frequency response: H(f) = sin(pi*N*f/fs) / (N*sin(pi*f/fs))
 *   - First null at f = fs/N
 *   - -3 dB cutoff at approximately f_c ≈ 0.443 * fs / N
 *   - Noise reduction: sigma_out = sigma_in / sqrt(N)
 *
 * Complexity: O(N) per call (direct method — sums entire window)
 *
 * @param buffer       Sliding window [window_size], shifted left, new_sample at end
 * @param window_size  Number of samples in the window (>0)
 * @param new_sample   Newest sample to incorporate
 * @return             Filtered (averaged) output
 */
double ss_cond_moving_average(double *buffer, size_t window_size,
                              double new_sample) {
    size_t i;
    double sum = 0.0;
    if (!buffer || window_size == 0) return new_sample;

    /* Shift buffer left and append new sample */
    for (i = 0; i < window_size - 1; i++) {
        buffer[i] = buffer[i + 1];
    }
    buffer[window_size - 1] = new_sample;

    /* Sum window */
    for (i = 0; i < window_size; i++) {
        sum += buffer[i];
    }

    return sum / (double)window_size;
}

/**
 * @brief Running moving average (O(1) per call)
 *
 * Maintains a running sum to achieve constant-time averaging.
 * Until the window fills, uses partial averaging.
 *
 * This is the preferred implementation for real-time systems where
 * per-sample computation must be bounded.
 *
 * @param window       Circular buffer [window_size]
 * @param window_size  Maximum window size
 * @param count        Running sample counter (updated in place)
 * @param idx          Circular buffer insert index (updated in place)
 * @param running_sum  Running sum accumulator (updated in place)
 * @param new_sample   New sample to include
 * @return             Running average (partial average if count < window_size)
 */
double ss_cond_running_moving_avg(double *window, size_t window_size,
                                  size_t *count, size_t *idx,
                                  double *running_sum, double new_sample) {
    if (!window || window_size == 0 || !count || !idx || !running_sum) {
        return new_sample;
    }

    if (*count < window_size) {
        /* Window not yet full: simple append */
        window[*count] = new_sample;
        *running_sum += new_sample;
        (*count)++;
        return *running_sum / (double)(*count);
    } else {
        /* Window full: remove oldest, add newest */
        double oldest = window[*idx];
        window[*idx] = new_sample;
        *running_sum = *running_sum - oldest + new_sample;
        *idx = (*idx + 1) % window_size;
        return *running_sum / (double)window_size;
    }
}

/* =========================================================================
 * L5: Digital Filtering — Median Filter
 * ========================================================================= */

/**
 * @brief Median filter for impulse noise removal
 *
 * Sorts the current window and returns the middle element.
 * Highly effective against salt-and-pepper noise (isolated outliers)
 * while preserving step edges — unlike linear filters which blur edges.
 *
 * Properties:
 *   - Nonlinear: output cannot exceed input range
 *   - Edge-preserving: sharp transitions survive filtering
 *   - Impulse rejection: single-sample outliers are completely removed
 *   - Window must be odd: 3, 5, or 7 samples
 *
 * Sorting: Insertion sort for small N (efficient for N <= 10)
 *
 * Complexity: O(N^2) for insertion sort, but N is very small (3-7)
 *
 * @param buffer       Sliding window [window_size] (shifted, updated in place)
 * @param window_size  Odd number: 3, 5, or 7
 * @param new_sample   Newest sample
 * @return             Median value of window
 */
double ss_cond_median_filter(double *buffer, size_t window_size,
                             double new_sample) {
    double sorted[7];
    size_t i, j;

    if (!buffer || window_size < 3 || window_size > 7
        || window_size % 2 == 0) {
        return new_sample;
    }

    /* Shift buffer left */
    for (i = 0; i < window_size - 1; i++) {
        buffer[i] = buffer[i + 1];
    }
    buffer[window_size - 1] = new_sample;

    /* Copy to sort buffer */
    for (i = 0; i < window_size; i++) {
        sorted[i] = buffer[i];
    }

    /* Insertion sort (efficient for small N) */
    for (i = 1; i < window_size; i++) {
        double key = sorted[i];
        j = i;
        while (j > 0 && sorted[j - 1] > key) {
            sorted[j] = sorted[j - 1];
            j--;
        }
        sorted[j] = key;
    }

    /* Return median: middle element of sorted array */
    return sorted[window_size / 2];
}

/* =========================================================================
 * L5: Digital Filtering — First-Order IIR Filters
 * ========================================================================= */

/**
 * @brief First-order IIR lowpass filter
 *
 * Difference equation:
 *   y[n] = alpha * x[n] + (1 - alpha) * y[n-1]
 *
 * Transfer function:
 *   H(z) = alpha / (1 - (1-alpha) * z^{-1})
 *
 * Frequency response magnitude at f:
 *   |H(f)| = alpha / sqrt(1 + (1-alpha)^2 - 2*(1-alpha)*cos(2*pi*f/fs))
 *
 * Design equation for desired cutoff frequency f_c:
 *   omega_c = 2*pi*f_c / fs
 *   alpha = -cos(omega_c) + sqrt(cos^2(omega_c) - 4*cos(omega_c) + 3) — 1
 *          ≈ 2*pi*f_c/fs              (for f_c << fs, small-alpha approximation)
 *
 * Time constant: tau = -1 / (fs * ln(1-alpha)) ≈ 1/(alpha*fs) for small alpha
 *
 * This is the simplest possible digital lowpass — one multiply, one add,
 * one subtraction per sample. Ubiquitous in embedded sensor systems.
 *
 * Complexity: O(1)
 *
 * @param prev_output  Previous filter output y[n-1] (updated in place)
 * @param input        Current input x[n]
 * @param alpha        Filter coefficient (0 < alpha <= 1)
 * @return             Filtered output y[n]
 */
double ss_cond_iir_lowpass(double *prev_output, double input, double alpha) {
    double output;
    if (!prev_output) return input;
    if (alpha <= 0.0) return *prev_output;  /* No update */
    if (alpha >= 1.0) { *prev_output = input; return input; } /* No filtering */

    output = alpha * input + (1.0 - alpha) * (*prev_output);
    *prev_output = output;
    return output;
}

/**
 * @brief First-order IIR highpass filter
 *
 * Difference equation:
 *   y[n] = (1-alpha) * (y[n-1] + x[n] - x[n-1])
 *
 * Transfer function:
 *   H(z) = (1-alpha) * (1 - z^{-1}) / (1 - (1-alpha) * z^{-1})
 *
 * DC gain: 0 (blocks DC completely)
 * High-frequency gain: 1 (passes through unchanged)
 *
 * Used for removing DC offset/bias from sensor signals.
 * For example, removing gravity from an accelerometer to isolate
 * dynamic acceleration.
 *
 * Complexity: O(1)
 *
 * @param prev_output  Previous filter output y[n-1] (updated in place)
 * @param prev_input   Previous input x[n-1] (updated in place)
 * @param input        Current input x[n]
 * @param alpha        Filter coefficient (0 < alpha <= 1)
 * @return             Filtered output y[n]
 */
double ss_cond_iir_highpass(double *prev_output, double *prev_input,
                            double input, double alpha) {
    double output;
    if (!prev_output || !prev_input) return 0.0;
    if (alpha <= 0.0) return 0.0;
    if (alpha >= 1.0) { output = input - *prev_input; goto update; }

    output = (1.0 - alpha) * (*prev_output + input - *prev_input);

update:
    *prev_input = input;
    *prev_output = output;
    return output;
}

/* =========================================================================
 * L5: Peak Detection with Exponential Decay
 * ========================================================================= */

/**
 * @brief Peak detector with controlled decay rate
 *
 * Tracks the maximum signal value with an exponential decay that
 * allows the peak to "forget" old maxima over time. Useful for:
 *   - Vibration envelope detection
 *   - Impact/shock event detection
 *   - Signal envelope tracking in demodulation
 *
 * Algorithm:
 *   if input > peak: peak = input          (instant attack)
 *   else:            peak *= decay_factor  (exponential decay)
 *
 * Attack time: 0 (instantaneous)
 * Decay time constant: tau = -1 / (fs * ln(decay_factor))
 *   For decay_factor=0.999 and fs=1000Hz: tau ≈ 1 second
 *
 * @param peak          Current peak value (updated in place)
 * @param input         Current input sample
 * @param decay_factor  Per-sample decay (0-1, 0.999 typical)
 * @return              Updated peak value
 */
double ss_cond_peak_detector(double *peak, double input, double decay_factor) {
    if (!peak) return input;
    if (input > *peak) {
        *peak = input;
    } else {
        *peak *= decay_factor;
    }
    return *peak;
}

/* =========================================================================
 * L5: Rate of Change Estimation
 * ========================================================================= */

/**
 * @brief Estimate first derivative using backward difference
 *
 * Simple derivative approximation:
 *   dx/dt[n] ≈ (x[n] - x[n-1]) / dt
 *
 * This is a first-order finite difference with truncation error O(dt).
 * For noisy signals, consider combining with a lowpass filter first.
 *
 * Amplifies high-frequency noise with gain = 1/dt.
 * For a sensor with 10 Hz sampling (dt=0.1s), noise is amplified by 10x.
 *
 * @param prev_value  Previous sample x[n-1] (updated in place to x[n])
 * @param curr_value  Current sample x[n]
 * @param dt          Sampling period (seconds)
 * @return            Estimated dx/dt (units/sec)
 */
double ss_cond_rate_of_change(double *prev_value, double curr_value,
                              double dt) {
    double rate;
    if (!prev_value || dt <= 0.0) return 0.0;

    rate = (curr_value - *prev_value) / dt;
    *prev_value = curr_value;
    return rate;
}

/* =========================================================================
 * L5: Threshold Detection with Hysteresis
 * ========================================================================= */

/**
 * @brief Schmitt-trigger style threshold detector
 *
 * Implements hysteresis to prevent rapid toggling near the threshold
 * boundary due to noise. Uses separate rising (on) and falling (off)
 * thresholds with a deadband between them.
 *
 * State machine:
 *   state=0 (inactive): if input > high_thresh -> state=1 (turn on)
 *   state=1 (active):   if input < low_thresh  -> state=0 (turn off)
 *   otherwise:          state unchanged (hysteresis band)
 *
 * Hysteresis width = high_thresh - low_thresh. Must be positive
 * and should be at least 2-3x the peak-to-peak noise amplitude.
 *
 * @param state        Current state (0 or 1, updated in place)
 * @param input        Current signal value
 * @param low_thresh   Falling threshold (state 1->0)
 * @param high_thresh  Rising threshold (state 0->1)
 * @return             Current state (0 or 1)
 */
int ss_cond_threshold_hysteresis(int *state, double input,
                                 double low_thresh, double high_thresh) {
    if (!state) return 0;

    if (*state == 0 && input >= high_thresh) {
        *state = 1;
    } else if (*state == 1 && input <= low_thresh) {
        *state = 0;
    }
    /* else: stay in current state (hysteresis) */

    return *state;
}

/* =========================================================================
 * L5: Sensor Linearization and Temperature Compensation
 * ========================================================================= */

/**
 * @brief Piecewise linear interpolation from calibration table
 *
 * Given a calibrated lookup table (x_i = known reference, y_i = raw reading),
 * performs binary search + linear interpolation to convert raw sensor
 * reading to engineering units.
 *
 * Interpolation formula for y_raw between y_table[i] and y_table[i+1]:
 *   x = x_table[i] + (y_raw - y_table[i])
 *       * (x_table[i+1] - x_table[i]) / (y_table[i+1] - y_table[i])
 *
 * Assumes y_table is monotonically increasing with x_table.
 *
 * Complexity: O(log n) for binary search, O(1) for interpolation
 *
 * @param x_table  Known reference values [n] (sorted ascending)
 * @param y_table  Raw sensor readings at each reference [n]
 * @param n        Number of table entries
 * @param y_input  Raw sensor reading to linearize
 * @return         Interpolated engineering units value
 */
double ss_cond_piecewise_linearize(const double *x_table,
                                   const double *y_table,
                                   size_t n, double y_input) {
    size_t lo, hi, mid;

    if (!x_table || !y_table || n < 2) return y_input;

    /* Handle out-of-range: clamp to endpoints */
    if (y_input <= y_table[0]) return x_table[0];
    if (y_input >= y_table[n - 1]) return x_table[n - 1];

    /* Binary search for the correct segment */
    lo = 0;
    hi = n - 1;
    while (hi - lo > 1) {
        mid = (lo + hi) / 2;
        if (y_table[mid] <= y_input) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    /* Linear interpolation between lo and hi */
    {
        double dy = y_table[hi] - y_table[lo];
        double dx = x_table[hi] - x_table[lo];
        if (dy == 0.0) return x_table[lo];
        return x_table[lo] + (y_input - y_table[lo]) * dx / dy;
    }
}

/**
 * @brief Temperature compensation: zero and span correction
 *
 * Corrects sensor reading for ambient temperature effects using
 * linear temperature coefficient model.
 *
 * The sensor output at temperature T is:
 *   reading(T) = (true_value * nominal_span) * (1 + TC_span * (T-T_ref))
 *                + zero_offset * (1 + TC_zero * (T-T_ref))
 *
 * Solving for true_value:
 *   zero_shift    = TC_zero * (T - T_ref)
 *   span_factor   = 1 + TC_span * (T - T_ref)
 *   true_value    = (reading - zero_shift) / span_factor
 *
 * Zero drift (offset drift): TC_zero in units/degC
 * Span drift (gain drift):   TC_span in fraction/degC (often ppm/degC)
 *
 * @param raw_reading    Uncompensated sensor output
 * @param current_temp   Current ambient temperature (degC)
 * @param ref_temp       Reference temperature from calibration (degC)
 * @param tc_zero        Zero drift coefficient (units/degC)
 * @param tc_span        Span drift coefficient (1/degC, fraction per degC)
 * @return               Temperature-compensated reading
 */
double ss_cond_temp_compensate(double raw_reading, double current_temp,
                               double ref_temp, double tc_zero,
                               double tc_span) {
    double dt = current_temp - ref_temp;
    double zero_shift = tc_zero * dt;
    double span_factor = 1.0 + tc_span * dt;

    if (span_factor == 0.0) return raw_reading; /* Degenerate */
    return (raw_reading - zero_shift) / span_factor;
}

/**
 * @brief Polynomial non-linearity correction
 *
 * Applies an inverse polynomial to linearize sensor output.
 * Given coefficients c_0...c_n (increasing degree), computes:
 *
 *   x_corrected = c_0 + c_1*y + c_2*y^2 + ... + c_n*y^n
 *
 * using Horner's method for numerical efficiency and stability.
 *
 * @param y_raw   Raw nonlinear sensor output
 * @param coeffs   Polynomial coefficients [order+1], lowest degree first
 * @param order    Polynomial order
 * @return         Linearized estimate
 */
double ss_cond_poly_nonlinearity_correct(double y_raw,
                                         const double *coeffs,
                                         uint8_t order) {
    int i;
    double result;

    if (!coeffs || order == 0) return y_raw;

    /* Horner's method: evaluate polynomial at y_raw */
    result = coeffs[order];
    for (i = (int)order - 1; i >= 0; i--) {
        result = result * y_raw + coeffs[i];
    }
    return result;
}

/* =========================================================================
 * L3: Generic Sensor Transfer Function Evaluation
 * ========================================================================= */

/**
 * @brief Evaluate a generic sensor transfer function
 *
 * Supports five common sensor model types with 4 parameters each:
 *   0: Linear      f(x) = a*x + b                         (params: [a, b, -, -])
 *   1: Exponential f(x) = a * exp(b*x) + c                (params: [a, b, c, -])
 *   2: Logarithmic f(x) = a * log(b*x) + c                (params: [a, b, c, -])
 *   3: Power law   f(x) = a * x^b + c                     (params: [a, b, c, -])
 *   4: Sigmoid     f(x) = a / (1 + exp(-b*(x-c))) + d     (params: [a, b, c, d])
 *
 * Used for sensor modeling when the physical principle dictates
 * a known functional form (e.g., exponential for thermistor,
 * power-law for semiconductor gas sensors).
 *
 * @param input       Physical input value x
 * @param params      Model parameters [4]: [a, b, c, d]
 * @param model_type  0=linear, 1=exp, 2=log, 3=power, 4=sigmoid
 * @return            Sensor output f(x)
 */
double ss_cond_transfer_function(double input, const double params[4],
                                 int model_type) {
    if (!params) return input;

    switch (model_type) {
    case 0: /* Linear: a*x + b */
        return params[0] * input + params[1];

    case 1: /* Exponential: a*exp(b*x) + c */
        return params[0] * exp(params[1] * input) + params[2];

    case 2: /* Logarithmic: a*log(b*x) + c */
        if (params[1] * input <= 0.0) return params[2];
        return params[0] * log(params[1] * input) + params[2];

    case 3: /* Power law: a*x^b + c */
        if (input < 0.0 && params[1] != floor(params[1])) return params[2];
        return params[0] * pow(input, params[1]) + params[2];

    case 4: /* Sigmoid (logistic): a/(1+exp(-b*(x-c))) + d */
        return params[0] / (1.0 + exp(-params[1] * (input - params[2])))
               + params[3];

    default:
        return input;
    }
}

/**
 * @brief Compute local sensitivity via central finite difference
 *
 * Sensitivity (local derivative) at operating point x0:
 *   S(x0) = df/dx|_{x=x0} ≈ [f(x0+h) - f(x0-h)] / (2h)
 *
 * Uses central difference for O(h^2) accuracy (vs. O(h) for forward diff).
 * The optimal step size h balances truncation error against roundoff:
 *   h_opt ≈ cbrt(eps) * |x0| ≈ 6e-6 * |x0|  (for double precision)
 *
 * @param input       Operating point x0
 * @param params      Transfer function parameters [4]
 * @param model_type  Model type (0-4)
 * @param h           Step size for finite difference
 * @return            Local sensitivity df/dx at x0
 */
double ss_cond_sensitivity_analysis(double input, const double params[4],
                                    int model_type, double h) {
    double f_plus, f_minus;
    if (h == 0.0) h = fabs(input) * 1e-6 + 1e-9;

    f_plus  = ss_cond_transfer_function(input + h, params, model_type);
    f_minus = ss_cond_transfer_function(input - h, params, model_type);

    return (f_plus - f_minus) / (2.0 * h);
}

/**
 * @brief Compute sensor dynamic range in dB
 *
 * Dynamic range is the ratio of largest to smallest detectable signal:
 *   DR_dB = 20 * log10(full_scale / resolution)
 *
 * For a sensor with 100 degC range and 0.01 degC resolution:
 *   DR_dB = 20 * log10(100/0.01) = 20 * log10(10000) = 80 dB
 *
 * @param full_scale   Maximum measurable input
 * @param resolution   Minimum detectable input change
 * @return             Dynamic range in dB
 */
double ss_cond_dynamic_range_db(double full_scale, double resolution) {
    if (resolution <= 0.0 || full_scale <= 0.0) return 0.0;
    return 20.0 * log10(full_scale / resolution);
}

/* =========================================================================
 * L4: Signal-to-Noise Ratio Fundamentals
 * ========================================================================= */

/**
 * @brief Compute Signal-to-Noise Ratio in dB
 *
 * SNR_dB = 10 * log10(P_signal / P_noise)
 *        = 20 * log10(V_signal_rms / V_noise_rms)
 *
 * Note: Factor of 10 for power ratio, 20 for amplitude ratio
 * because power is proportional to amplitude squared.
 *
 * @param signal_rms  RMS signal amplitude
 * @param noise_rms   RMS noise amplitude
 * @return            SNR in dB
 */
double ss_cond_snr_db(double signal_rms, double noise_rms) {
    if (noise_rms <= 0.0) return INFINITY;
    return 20.0 * log10(signal_rms / noise_rms);
}

/**
 * @brief Compute Minimum Detectable Signal (MDS)
 *
 * The smallest signal amplitude that can be detected above the noise
 * floor at a specified SNR threshold:
 *
 *   MDS = noise_rms * 10^(SNR_min_dB / 20)
 *
 * For SNR_min = 0 dB: MDS = noise_rms (signal = noise floor)
 * For SNR_min = 10 dB: MDS = noise_rms * 3.16
 * For SNR_min = 20 dB: MDS = noise_rms * 10
 *
 * In radar/sonar, SNR_min for detection is typically 12-15 dB
 * (corresponding to Pd=0.9, Pfa=1e-6 for a single pulse).
 *
 * @param noise_rms   RMS noise level
 * @param snr_min_db  Minimum SNR required for reliable detection
 * @return            Minimum detectable signal amplitude
 */
double ss_cond_min_detectable_signal(double noise_rms, double snr_min_db) {
    return noise_rms * pow(10.0, snr_min_db / 20.0);
}
