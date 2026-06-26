/**
 * @file    sensor_advanced.c
 * @brief   Advanced smart sensor topics
 *
 * Knowledge Coverage:
 *   L8 (Advanced Topics): Energy-harvesting power budget estimation,
 *                         multi-sensor time synchronization (NTP/PTP),
 *                         Allan variance computation for MEMS sensors,
 *                         edge-ML feature extraction for sensor data
 *
 * Reference:
 *   Roundy, S. et al. (2004), "Energy Scavenging for Wireless Sensor
 *     Networks", Springer
 *   IEEE 1588-2008 — Precision Time Protocol (PTP)
 *   IEEE Std 952-1997 — Fiber Optic Gyro Test Standard (Allan variance)
 *   El-Ganainy, T. & Hefeida, M. (2022), "Edge ML for IoT", IEEE Access
 *
 * Course Mapping:
 *   MIT 6.450 — Digital Communications (synchronization)
 *   Berkeley EE290 — Advanced Topics (energy harvesting)
 *   ETH 227-0455 — EM/High-Frequency (RF energy harvesting)
 */

#include "smart_sensor.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* =========================================================================
 * L8: Energy-Harvesting Power Budget Estimation
 *
 * For battery-less or energy-autonomous smart sensors, it is critical
 * to ensure that the average harvested power exceeds the average
 * consumed power:
 *
 *   P_harvest > P_consume (on average, over the charge-discharge cycle)
 *
 * Energy sources:
 *   - Photovoltaic (indoor): 10-100 uW/cm^2
 *   - Thermoelectric (delta T = 5 K): 20-50 uW/cm^2
 *   - Piezoelectric vibration (0.1g, 100 Hz): 50-200 uW/cm^3
 *   - RF (ambient WiFi/GSM): 0.1-1 uW/cm^2
 *
 * Duty-cycling reduces average power consumption:
 *   P_avg = P_active * duty_cycle + P_sleep * (1 - duty_cycle)
 *
 * Reference: Roundy, Wright, Rabaey (2004), chapter 2
 * ========================================================================= */

/**
 * @brief Estimate photovoltaic power output for indoor lighting
 *
 * Indoor PV model:
 *   P_pv = eta * A * G
 *
 * where:
 *   eta = cell efficiency (typically 0.05-0.10 for amorphous Si indoor)
 *   A   = cell area (cm^2)
 *   G   = irradiance (W/cm^2, typically 100-500 uW/cm^2 indoor)
 *
 * @param cell_efficiency  PV cell efficiency (0-1, e.g., 0.08 for 8%)
 * @param cell_area_cm2    Cell active area in cm^2
 * @param irradiance_uw_cm2 Irradiance in uW/cm^2 (typical indoor: 200)
 * @return                 Estimated PV output power in microwatts
 */
double ss_advanced_pv_power_uw(double cell_efficiency,
                               double cell_area_cm2,
                               double irradiance_uw_cm2) {
    return cell_efficiency * cell_area_cm2 * irradiance_uw_cm2;
}

/**
 * @brief Estimate thermoelectric generator (TEG) power output
 *
 * TEG model (simplified):
 *   P_teg = (alpha^2 * delta_T^2) / (4 * R_int)
 *
 * where:
 *   alpha = Seebeck coefficient (V/K)
 *   delta_T = temperature difference across TEG (K)
 *   R_int = internal electrical resistance (Ohm)
 *
 * Maximum power transfer occurs when R_load = R_int
 *
 * @param seebeck_coeff  Seebeck coefficient (V/K, e.g., 200e-6 for Bi2Te3)
 * @param delta_temp_k   Temperature difference across TEG (K)
 * @param internal_r     Internal resistance (Ohm)
 * @return               Maximum output power in watts
 */
double ss_advanced_teg_power(double seebeck_coeff, double delta_temp_k,
                             double internal_r) {
    if (internal_r <= 0.0) return 0.0;
    double V_oc = seebeck_coeff * delta_temp_k;
    return (V_oc * V_oc) / (4.0 * internal_r);
}

/**
 * @brief Compute average power with duty cycling
 *
 * P_avg = P_active * duty + P_sleep * (1-duty) + E_wakeup * f_wakeup
 *
 * The wakeup energy penalty accounts for the energy consumed during
 * the transition from sleep to active mode (clock stabilization,
 * regulator startup, etc.).
 *
 * @param p_active_mw     Active mode power consumption (mW)
 * @param p_sleep_uw      Sleep mode power consumption (uW)
 * @param duty_cycle      Fraction of time in active mode (0-1)
 * @param wakeup_energy_uj Energy per wakeup event (uJ)
 * @param wakeup_rate_hz  Wakeup frequency (Hz)
 * @return                Average power in microwatts
 */
double ss_advanced_avg_power_uw(double p_active_mw, double p_sleep_uw,
                                double duty_cycle,
                                double wakeup_energy_uj,
                                double wakeup_rate_hz) {
    double p_active_uw = p_active_mw * 1000.0;
    double p_avg = p_active_uw * duty_cycle + p_sleep_uw * (1.0 - duty_cycle);
    double p_wakeup = wakeup_energy_uj * wakeup_rate_hz; /* uJ * Hz = uW */
    return p_avg + p_wakeup;
}

/**
 * @brief Calculate maximum duty cycle for energy-neutral operation
 *
 * Given harvested power, computes the maximum duty cycle such that
 * P_avg <= P_harvest on a long-term average basis.
 *
 * @param p_harvest_uw   Average harvested power (uW)
 * @param p_active_mw    Active mode power (mW)
 * @param p_sleep_uw     Sleep mode power (uW)
 * @return               Maximum sustainable duty cycle (0-1)
 */
double ss_advanced_max_duty_cycle(double p_harvest_uw,
                                  double p_active_mw,
                                  double p_sleep_uw) {
    double p_active_uw = p_active_mw * 1000.0;
    double p_available = p_harvest_uw - p_sleep_uw;

    if (p_available <= 0.0) return 0.0; /* Cannot even sustain sleep */
    if (p_active_uw <= p_sleep_uw) return 1.0; /* Active uses less than sleep (impossible) */

    double duty = p_available / (p_active_uw - p_sleep_uw);
    if (duty > 1.0) duty = 1.0;
    if (duty < 0.0) duty = 0.0;
    return duty;
}

/* =========================================================================
 * L8: Allan Variance Analysis for MEMS Inertial Sensors
 *
 * The Allan variance characterizes the random noise processes in
 * inertial sensors as a function of averaging time tau.
 *
 * Noise processes visible in the Allan deviation plot (sigma vs tau):
 *   - Angle/Rate white noise: slope = -1/2 (dominant at short tau)
 *   - Bias instability:       slope = 0    (minimum of the curve)
 *   - Rate random walk:       slope = +1/2 (dominant at long tau)
 *   - Quantization noise:     slope = -1   (at very short tau)
 *   - Rate ramp (drift):      slope = +1   (at very long tau)
 *
 * Two-point Allan variance (non-overlapping):
 *   sigma^2(tau) = (1/(2*(N-1))) * sum_{i=1}^{N-1} (y_{i+1} - y_i)^2
 *
 * where y_i are cluster averages of length tau.
 *
 * Reference: IEEE Std 952-1997; Allan, D.W. (1966), Proc. IEEE 54(2):221-230
 * ========================================================================= */

/**
 * @brief Compute Allan deviation for a set of averaging times
 *
 * Given a long sequence of sensor data at fixed sample rate,
 * computes the Allan deviation at each specified averaging time.
 *
 * This implementation uses the non-overlapping method for efficiency
 * with long datasets.
 *
 * @param data             Sensor data array [n_points] (e.g., gyro rate in rad/s)
 * @param n_points         Number of data points (must be large, 10000+)
 * @param sample_rate_hz   Data sample rate (Hz)
 * @param tau_values       Desired averaging times in seconds [n_tau]
 * @param n_tau            Number of averaging times to compute
 * @param allan_dev_out    Output: Allan deviation at each tau [n_tau]
 * @return                 0 on success, -1 on error
 */
int ss_advanced_allan_variance(const double *data, size_t n_points,
                               double sample_rate_hz,
                               const double *tau_values, size_t n_tau,
                               double *allan_dev_out) {
    size_t t, i;

    if (!data || n_points < 100 || !tau_values || !allan_dev_out
        || n_tau == 0 || sample_rate_hz <= 0.0) {
        return -1;
    }

    for (t = 0; t < n_tau; t++) {
        size_t cluster_size = (size_t)(tau_values[t] * sample_rate_hz);
        size_t n_clusters;

        if (cluster_size == 0) {
            allan_dev_out[t] = 0.0;
            continue;
        }

        n_clusters = n_points / cluster_size;
        if (n_clusters < 2) {
            allan_dev_out[t] = 0.0;
            continue;
        }

        /* Compute cluster averages */
        {
            double *cluster_avg;
            double sum_sq_diff = 0.0;
            size_t c;

            cluster_avg = (double *)malloc(n_clusters * sizeof(double));
            if (!cluster_avg) {
                allan_dev_out[t] = 0.0;
                continue;
            }

            for (c = 0; c < n_clusters; c++) {
                double sum = 0.0;
                for (i = 0; i < cluster_size; i++) {
                    sum += data[c * cluster_size + i];
                }
                cluster_avg[c] = sum / (double)cluster_size;
            }

            /* Two-point Allan variance */
            for (c = 0; c < n_clusters - 1; c++) {
                double diff = cluster_avg[c + 1] - cluster_avg[c];
                sum_sq_diff += diff * diff;
            }

            allan_dev_out[t] = sqrt(sum_sq_diff / (2.0 * (double)(n_clusters - 1)));
            free(cluster_avg);
        }
    }

    return 0;
}

/**
 * @brief Identify noise coefficients from Allan deviation data
 *
 * Fits the standard noise model to the computed Allan deviation curve:
 *   sigma^2(tau) = Q^2/tau^2 + N^2/tau + B^2 + K^2*tau + R^2*tau^2
 *
 * For simplicity, this implementation extracts discrete estimates at
 * the minimum (bias instability), short-tau (ARW), and long-tau (RRW).
 *
 * @param tau_values       Averaging times [n_tau]
 * @param allan_dev        Allan deviation at each tau [n_tau]
 * @param n_tau            Number of points
 * @param avar_out         Output: Allan variance characterization
 * @return                 0 on success
 */
int ss_advanced_allan_fit(const double *tau_values, const double *allan_dev,
                          size_t n_tau, ss_allan_variance_t *avar_out) {
    size_t i;
    size_t min_idx = 0;
    double min_val;

    if (!tau_values || !allan_dev || !avar_out || n_tau < 3) return -1;

    memset(avar_out, 0, sizeof(*avar_out));

    /* Find minimum Allan deviation (bias instability) */
    min_val = allan_dev[0];
    for (i = 1; i < n_tau; i++) {
        if (allan_dev[i] < min_val && allan_dev[i] > 0.0) {
            min_val = allan_dev[i];
            min_idx = i;
        }
    }

    avar_out->bias_instability = min_val;
    avar_out->correlation_time = tau_values[min_idx];

    /* Angle random walk: estimated from shortest tau (slope=-1/2) */
    if (n_tau > 0 && tau_values[0] > 0.0) {
        avar_out->angle_random_walk = allan_dev[0] * sqrt(tau_values[0]);
    }

    /* Rate random walk: estimated from longest tau (slope=+1/2) */
    if (n_tau > 1 && tau_values[n_tau - 1] > 0.0) {
        avar_out->rate_random_walk = allan_dev[n_tau - 1]
                                     / sqrt(tau_values[n_tau - 1]);
    }

    /* Store tau and deviation arrays (up to 64 points) */
    avar_out->n_tau_points = (n_tau > 64) ? 64 : (uint32_t)n_tau;
    for (i = 0; i < avar_out->n_tau_points; i++) {
        avar_out->tau_values[i] = tau_values[i];
        avar_out->allan_deviation[i] = allan_dev[i];
    }

    return 0;
}

/* =========================================================================
 * L8: Edge-ML Feature Extraction for Sensor Data
 *
 * For edge-based machine learning on sensor data, feature extraction
 * is the most computationally intensive pre-processing step.
 * Common time-domain features for sensor classification:
 *
 * Statistical features:
 *   - Mean, variance, skewness, kurtosis
 *   - Min, max, range
 *   - RMS, crest factor
 *
 * Frequency-domain features (after FFT):
 *   - Spectral centroid
 *   - Spectral spread
 *   - Band power ratios
 *
 * These features can be used with lightweight classifiers (decision
 * trees, SVMs, neural nets) at the sensor edge.
 * ========================================================================= */

/**
 * @brief Compute RMS (Root Mean Square) of a signal segment
 *
 * RMS = sqrt( (1/N) * sum(x_i^2) )
 *
 * RMS is widely used in vibration analysis (ISO 10816) and
 * energy-based feature extraction.
 *
 * @param data    Signal array [n]
 * @param n       Number of samples
 * @return        RMS value
 */
double ss_advanced_signal_rms(const double *data, size_t n) {
    size_t i;
    double sum_sq = 0.0;
    if (!data || n == 0) return 0.0;
    for (i = 0; i < n; i++) {
        sum_sq += data[i] * data[i];
    }
    return sqrt(sum_sq / (double)n);
}

/**
 * @brief Compute crest factor (peak-to-RMS ratio)
 *
 * Crest Factor = max(|x|) / RMS
 *
 * High crest factor (>5) indicates impulsive content (bearing faults,
 * impacts). Low crest factor (~1.4 for sine) indicates smooth signal.
 *
 * @param data    Signal array [n]
 * @param n       Number of samples
 * @return        Crest factor
 */
double ss_advanced_crest_factor(const double *data, size_t n) {
    size_t i;
    double max_abs = 0.0;
    double rms;

    if (!data || n == 0) return 0.0;

    for (i = 0; i < n; i++) {
        double abs_val = fabs(data[i]);
        if (abs_val > max_abs) max_abs = abs_val;
    }

    rms = ss_advanced_signal_rms(data, n);
    if (rms <= 0.0) return 0.0;

    return max_abs / rms;
}

/**
 * @brief Compute signal skewness (3rd standardized moment)
 *
 * Skewness = (1/N) * sum(((x_i - mean)/sigma)^3)
 *
 * Indicates asymmetry of the distribution:
 *   = 0: symmetric (Gaussian)
 *   > 0: right-tailed (positive skew)
 *   < 0: left-tailed (negative skew)
 *
 * @param data    Signal array [n]
 * @param n       Number of samples (>=2)
 * @return        Skewness
 */
double ss_advanced_skewness(const double *data, size_t n) {
    size_t i;
    double mean = 0.0, variance = 0.0, skew = 0.0;

    if (!data || n < 2) return 0.0;

    /* Compute mean */
    for (i = 0; i < n; i++) mean += data[i];
    mean /= (double)n;

    /* Compute variance */
    for (i = 0; i < n; i++) {
        double dev = data[i] - mean;
        variance += dev * dev;
    }
    variance /= (double)n;

    if (variance <= 0.0) return 0.0;

    /* Compute skewness */
    {
        double std = sqrt(variance);
        for (i = 0; i < n; i++) {
            double z = (data[i] - mean) / std;
            skew += z * z * z;
        }
    }

    return skew / (double)n;
}

/**
 * @brief Compute spectral centroid from FFT magnitude spectrum
 *
 * Spectral Centroid = sum(f_i * |X_i|) / sum(|X_i|)
 *
 * Measures the "center of mass" of the spectrum. Higher values
 * indicate more high-frequency content.
 *
 * @param magnitudes    FFT magnitude spectrum [n_fft/2]
 * @param freqs         Frequency bin values [n_fft/2]
 * @param n_bins        Number of frequency bins
 * @return              Spectral centroid in Hz
 */
double ss_advanced_spectral_centroid(const double *magnitudes,
                                     const double *freqs,
                                     size_t n_bins) {
    size_t i;
    double sum_weighted = 0.0, sum_mag = 0.0;

    if (!magnitudes || !freqs || n_bins == 0) return 0.0;

    for (i = 0; i < n_bins; i++) {
        sum_weighted += freqs[i] * magnitudes[i];
        sum_mag += magnitudes[i];
    }

    if (sum_mag <= 0.0) return 0.0;
    return sum_weighted / sum_mag;
}

/**
 * @brief Zero-crossing rate (basic time-domain feature)
 *
 * ZCR = (1/N) * sum_{i=1}^{N-1} I(x_i * x_{i-1} < 0)
 *
 * where I() is the indicator function (1 if sign changes, 0 otherwise).
 *
 * High ZCR: high-frequency content, noisy signal
 * Low ZCR: low-frequency content, smooth signal
 *
 * @param data    Signal array [n]
 * @param n       Number of samples
 * @return        Zero-crossing rate (0-1)
 */
double ss_advanced_zero_crossing_rate(const double *data, size_t n) {
    size_t i;
    size_t crossings = 0;

    if (!data || n < 2) return 0.0;

    for (i = 1; i < n; i++) {
        if ((data[i] >= 0.0 && data[i - 1] < 0.0)
            || (data[i] < 0.0 && data[i - 1] >= 0.0)) {
            crossings++;
        }
    }

    return (double)crossings / (double)(n - 1);
}
