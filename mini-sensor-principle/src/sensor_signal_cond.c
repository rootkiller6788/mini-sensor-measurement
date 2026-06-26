/**
 * @file sensor_signal_cond.c
 * @brief Sensor signal conditioning — amplification, filtering, linearization, excitation
 *
 * Levels: L2 — Core Concepts · L5 — Algorithms/Methods
 */

#include "sensor_signal_cond.h"
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ═══════════════════════════════════════════════════════════════════════
 * L2: Operational Amplifier Gain Calculations
 * ═══════════════════════════════════════════════════════════════════════ */

double amp_non_inverting_gain(double Rf, double Rg)
{
    if (Rg <= 0.0) return 1.0;  /* unity-gain buffer (voltage follower) */
    return 1.0 + Rf / Rg;
}

double amp_inverting_gain(double Rf, double Rin)
{
    if (Rin <= 0.0) return 0.0;
    return -Rf / Rin;
}

double amp_differential_gain(double Rf, double Rin)
{
    if (Rin <= 0.0) return 0.0;
    return Rf / Rin;
}

/**
 * @brief Estimate CMRR degradation due to resistor mismatch.
 *
 * For a differential amplifier with ideal op-amp but mismatched resistors:
 *   CMRR ≈ (1 + G) / (4·δ)
 *   where δ = ΔR/R is the resistor tolerance, G = Rf/Rin.
 *
 * This shows that high-gain diff amps are very sensitive to resistor matching.
 * At G=100 with 1% resistors: CMRR ≈ 101/(4·0.01) = 2525 ≈ 68 dB
 * → 68 dB CMRR is marginal for bridge sensors → need 0.1% resistors
 *   or integrated instrumentation amplifier.
 */
double amp_differential_cmrr_estimate(double gain, double resistor_tolerance)
{
    if (resistor_tolerance <= 0.0) return INFINITY;
    return (1.0 + gain) / (4.0 * resistor_tolerance);
}

/**
 * @brief Three-op-amp instrumentation amplifier gain.
 *
 * Input stage: two non-inverting amplifiers sharing Rg.
 *   G1 = 1 + 2·Rf/Rg
 *
 * Output stage: unity-gain difference amplifier (R2 = R1 typically).
 *   G2 = 1 (or R2/R1 if not unity)
 *
 * Total: G = (1 + 2·Rf/Rg) · (R2/R1)
 *
 * The key advantage is that Rg sets the entire gain with a single resistor.
 * Input impedance: > 10⁹ Ω at each input.
 * CMRR: typically > 100 dB for integrated solutions (laser-trimmed resistors).
 *
 * Reference: Kitchin & Counts, "A Designer's Guide to Instrumentation
 *            Amplifiers", 3rd Ed, Analog Devices, 2006
 */
double amp_instrumentation_gain(double Rf, double Rg, double R2, double R1)
{
    if (Rg <= 0.0) return 0.0;
    if (R1 <= 0.0) return 0.0;
    double G1 = 1.0 + 2.0 * Rf / Rg;
    double G2 = R2 / R1;
    return G1 * G2;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L5: Anti-Alias Filter Design
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Compute cutoff frequency for anti-alias filter.
 *
 * Rule of thumb: f_c ≈ f_s / (2·OSR), where OSR ≥ 10.
 *
 * For f_s = 1 kHz, OSR=10 → f_c = 50 Hz.
 * At Nyquist (500 Hz): attenuation ≈ 20·log₁₀(1/√(1+10²)) ≈ -20 dB.
 *
 * For more aggressive alias protection (OSR=100):
 *   f_c = 5 Hz, Nyquist attenuation ≈ -40 dB.
 */
double anti_alias_fc_from_fs(double fs_Hz, double oversampling_factor)
{
    if (fs_Hz <= 0.0 || oversampling_factor <= 1.0) return fs_Hz / 2.0;
    return fs_Hz / (2.0 * oversampling_factor);
}

/**
 * @brief Compute anti-alias filter attenuation at Nyquist frequency.
 *
 * For an Nth-order Butterworth low-pass:
 *   |H(f)|² = 1 / (1 + (f/f_c)^{2N})
 *
 * At f = f_s/2: H_dB = -10·log₁₀(1 + (f_s/(2f_c))^{2N})
 *
 * This determines whether the alias rejection is sufficient
 * for the target ADC effective number of bits.
 */
double anti_alias_attenuation_at_nyquist(double fc_Hz, double fs_Hz, int order)
{
    if (fc_Hz <= 0.0 || fs_Hz <= 0.0 || order < 1) return 0.0;
    double ratio = fs_Hz / (2.0 * fc_Hz);
    return -10.0 * log10(1.0 + pow(ratio, 2.0 * order));
}

/**
 * @brief Choose standard E24 resistor value and compute capacitor.
 *
 * Picks R from standard E24 values (10% tolerance series):
 *   1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7,
 *   3.0, 3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1
 *
 * Then: C = 1/(2π·R·f_c)
 *
 * For typical sensor bandwidths (DC to 1 kHz),
 * R = 1-100 kΩ, C = 1.5 nF to 0.1 µF.
 */
int anti_alias_rc_values(double fc_target, double *R_out, double *C_out)
{
    if (!R_out || !C_out || fc_target <= 0.0) return -1;

    /* E24 standard values (×10^n) */
    const double e24[] = {1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7,
                          3.0, 3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8,
                          7.5, 8.2, 9.1};

    /* Try standard R values from 1k to 100k */
    double best_R = 10000.0, best_C = 0.0;
    double best_error = INFINITY;

    for (int decade = 3; decade <= 5; decade++) {
        double scale = pow(10.0, decade);
        for (int j = 0; j < 24; j++) {
            double R = e24[j] * scale;
            double C = 1.0 / (2.0 * M_PI * R * fc_target);
            /* Check if C is in a practical range (100 pF to 10 µF) */
            if (C >= 100e-12 && C <= 10e-6) {
                /* Compute actual fc with this R and C */
                double fc_actual = 1.0 / (2.0 * M_PI * R * C);
                double error = fabs(fc_actual - fc_target) / fc_target;
                if (error < best_error) {
                    best_error = error;
                    best_R = R;
                    best_C = C;
                }
            }
        }
    }

    *R_out = best_R;
    *C_out = best_C;
    return (best_error < 0.1) ? 0 : 1;  /* 0 = good match, 1 = approximate */
}

/**
 * @brief Compute required filter order for given alias rejection.
 *
 * From Butterworth LP attenuation formula:
 *   N ≥ log₁₀(10^{A_dB/10} - 1) / [2·log₁₀(f_s/(2f_c))]
 *
 * Example: f_s=1kHz, f_c=100Hz, want A=40dB
 *          N ≥ log₁₀(10^4-1) / [2·log₁₀(5)] = log₁₀(9999)/(2·0.699)
 *            = 4/(1.398) = 2.86 → need 3rd order
 */
int anti_alias_filter_order_required(double fs_hz, double fc_hz,
                                      double attenuation_db)
{
    if (fs_hz <= 0.0 || fc_hz <= 0.0 || fc_hz >= fs_hz / 2.0
        || attenuation_db <= 0.0) {
        return 1;
    }

    double ratio = fs_hz / (2.0 * fc_hz);
    if (ratio <= 1.0) return 1;

    double num = log10(pow(10.0, attenuation_db / 10.0) - 1.0);
    double den = 2.0 * log10(ratio);

    if (den <= 0.0) return 1;
    int N = (int)ceil(num / den);
    if (N < 1) N = 1;
    return N;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L2: 4-20 mA Current Loop
 * ═══════════════════════════════════════════════════════════════════════ */

double current_loop_sensor_to_current(double sensor_value,
                                       double x_min, double x_max)
{
    if (x_max <= x_min) return 4.0;
    double fraction = (sensor_value - x_min) / (x_max - x_min);
    if (fraction < 0.0) fraction = 0.0;
    if (fraction > 1.0) fraction = 1.0;
    return 4.0 + 16.0 * fraction;
}

double current_loop_current_to_sensor(double I_loop,
                                       double x_min, double x_max)
{
    if (x_max <= x_min) return x_min;
    if (I_loop < 4.0) return x_min;
    if (I_loop > 20.0) return x_max;
    return x_min + (x_max - x_min) * (I_loop - 4.0) / 16.0;
}

double current_loop_sense_resistor(double V_adc_fullscale)
{
    if (V_adc_fullscale <= 0.0) return 0.0;
    return V_adc_fullscale / 0.020;  /* 250 Ω for 5V ADC */
}

/**
 * @brief Maximum loop distance (cable length) based on supply voltage.
 *
 * The loop transmitter drops some minimum voltage (typically 7-12V).
 * The remaining voltage is available for the sensing resistor and cable drops.
 *
 * R_max_total = (V_supply - V_min_transmitter) / 0.020
 * R_cable = R_max_total - R_sense
 * L = R_cable / (2·R_per_meter)
 *
 * This is important for industrial installations where the sensor
 * may be hundreds of meters from the control room.
 */
double current_loop_max_distance(double V_supply, double V_min_transmitter,
                                  double R_sense, double R_per_meter)
{
    if (V_supply <= V_min_transmitter || R_per_meter <= 0.0) return 0.0;

    double R_max_total = (V_supply - V_min_transmitter) / 0.020;
    double R_cable = R_max_total - R_sense;
    if (R_cable <= 0.0) return 0.0;
    return R_cable / (2.0 * R_per_meter);
}

/* ═══════════════════════════════════════════════════════════════════════
 * L5: Sensor Linearization Functions
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Square-root extraction for flow sensors.
 *
 * For orifice plate, Venturi, and Pitot tube flow measurements:
 *   ΔP ∝ Q²  →  Q = K·√(ΔP)
 *
 * The square-root relationship means that at low flows, small pressure
 * measurement errors cause large flow errors. At 10% of full flow,
 * ΔP is at 1% of full scale.
 *
 * K is typically calibrated empirically following ISO 5167.
 */
double linearize_sqrt_law(double differential_pressure, double K)
{
    if (differential_pressure < 0.0) return 0.0;
    return K * sqrt(differential_pressure);
}

/**
 * @brief Logarithmic compression for wide dynamic range.
 *
 * Using a diode-connected transistor in the feedback loop:
 *   V_out = -V_T · ln(I_in / I_s)
 *
 * where V_T = kT/q ≈ 26 mV at 300 K, I_s ≈ 10⁻¹⁴ to 10⁻¹² A.
 *
 * This compresses many decades of input current into a manageable
 * voltage range. Commonly used in photodiode amplifiers for
 * optical power measurement spanning nW to mW.
 *
 * The log amp needs temperature compensation since V_T ∝ T.
 * This is typically done by a matched reference transistor.
 */
double linearize_log_compression(double I_signal, double I_reference, double V_T)
{
    if (I_signal <= 0.0 || I_reference <= 0.0 || V_T <= 0.0) return 0.0;
    return -V_T * log(I_signal / I_reference);
}

/**
 * @brief Active bridge linearization using op-amp feedback.
 *
 * By placing the bridge's active arm as the feedback element of an
 * inverting amplifier:
 *   V_out = -V_ref · (ΔR / R₀)
 *
 * This completely eliminates the bridge nonlinearity inherent in the
 * passive voltage-divider approach. The op-amp forces a constant current
 * through the active arm, making V_out linearly proportional to ΔR.
 *
 * Reference: Williams, "Bridge Circuits: Marrying Gain and Balance",
 *            Linear Technology Application Note 43, 1990
 */
double linearize_bridge_active(double delta_R, double R0, double V_ref)
{
    if (R0 <= 0.0) return 0.0;
    return -V_ref * delta_R / R0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L2: Sensor Excitation
 * ═══════════════════════════════════════════════════════════════════════ */

double excitation_stability_required(double total_accuracy_percent)
{
    return total_accuracy_percent / 3.0;
}

/**
 * @brief Maximum excitation current to avoid excessive self-heating.
 *
 * I_max = √(ΔT_max / (R_sensor · θ_ja))
 *
 * where θ_ja is thermal resistance junction-to-ambient [K/W].
 *
 * For RTD Pt100 with θ_ja ≈ 50 K/W in still air,
 * limiting self-heating to 0.01 K:
 *   I_max = √(0.01 / (100·50)) = √(2e-6) = 1.4 mA
 *
 * For 10 kΩ thermistor with θ_ja ≈ 200 K/W,
 * limiting self-heating to 0.1 K:
 *   I_max = √(0.1/(10000·200)) = √(5e-8) = 0.22 mA
 */
double excitation_current_max_for_self_heating(double R_sensor,
                                                double thermal_resistance,
                                                double max_temp_rise)
{
    if (R_sensor <= 0.0 || thermal_resistance <= 0.0 || max_temp_rise <= 0.0)
        return 0.0;

    return sqrt(max_temp_rise / (R_sensor * thermal_resistance));
}

/**
 * @brief Compute reactance for AC excitation.
 *
 * Capacitive: X_C = 1/(2π·f·C)  →  Z decreases with frequency
 * Inductive:  X_L = 2π·f·L      →  Z increases with frequency
 *
 * AC excitation is used for:
 *   - Capacitive sensors (displacement, humidity, touch)
 *   - Inductive sensors (LVDT, proximity, eddy current)
 *   - Avoiding DC errors (thermal EMF, offset drift)
 *
 * For capacitive MEMS accelerometers: C varies ~0.1-1 pF
 * at f=100 kHz: X_C ≈ 1/(2π·10⁵·10⁻¹²) ≈ 1.6 MΩ
 */
double ac_excitation_reactance(double f_hz, double C_or_L, int is_capacitive)
{
    if (f_hz <= 0.0 || C_or_L <= 0.0) return INFINITY;
    if (is_capacitive) {
        return 1.0 / (2.0 * M_PI * f_hz * C_or_L);
    } else {
        return 2.0 * M_PI * f_hz * C_or_L;
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * L5: Digital Filtering for Sensor Signal Processing
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Simple moving average (boxcar) filter.
 *
 * y[k] = (1/N)·Σ_{i=0}^{N-1} x[k-i]
 *
 * Transfer function (z-domain):
 *   H(z) = (1/N)·(1 - z^{-N})/(1 - z^{-1})
 *
 * Frequency response nulls at f = k·fs/N, k = 1, 2, ...
 *
 * Reduces white noise RMS by factor √N.
 * Impulse response: rectangular (hence "boxcar").
 *
 * Efficient implementation using running sum:
 *   sum[k] = sum[k-1] + x[k] - x[k-N]
 *   y[k] = sum[k]/N
 */
int signal_moving_average(const double *input, double *output,
                           size_t n_samples, size_t window_size)
{
    if (!input || !output || n_samples == 0 || window_size == 0
        || window_size > n_samples) {
        return -1;
    }

    /* Running sum for efficiency */
    double sum = 0.0;
    for (size_t k = 0; k < n_samples; k++) {
        sum += input[k];
        if (k >= window_size) {
            sum -= input[k - window_size];
        }
        size_t valid_count = (k + 1 < window_size) ? (k + 1) : window_size;
        output[k] = sum / (double)valid_count;
    }
    return 0;
}

/**
 * @brief First-order exponential (IIR) low-pass filter.
 *
 * y[k] = α·x[k] + (1-α)·y[k-1]
 *
 * This is the discrete-time equivalent of an RC filter.
 *
 * Time constant: τ = T_s / α  (for small α)
 * Cutoff: f_c = α / (2π·T_s)
 *
 * α = 0: no filtering (y[k] = y[k-1], DC only)
 * α = 1: no filtering (y[k] = x[k], passthrough)
 * α = 0.1: moderate filtering
 *
 * Equivalent to: y[k] = y[k-1] + α·(x[k] - y[k-1])
 * The error-correction form shows the filter tracks changes.
 */
double signal_exponential_filter(double x_new, double y_prev, double alpha)
{
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;
    return alpha * x_new + (1.0 - alpha) * y_prev;
}

/**
 * @brief Median filter for impulse/outlier removal.
 *
 * For each sample, computes median of the window [k-w/2, k+w/2].
 *
 * Properties:
 *   - Excellent at removing impulse noise (salt-and-pepper)
 *   - Preserves edges (unlike moving average which blurs them)
 *   - Destroys small signal features smaller than w/2
 *   - Not linear → frequency-domain analysis not applicable
 *   - Computationally more expensive than moving average
 *
 * For streaming applications, a running median can be computed
 * using two heaps (min-heap + max-heap) in O(log w) per sample.
 */
int signal_median_filter(const double *input, double *output,
                          size_t n, size_t window)
{
    if (!input || !output || n == 0 || window < 3 || window % 2 == 0)
        return -1;

    size_t half = window / 2;
    double *buf = (double *)calloc(window, sizeof(double));
    if (!buf) return -1;

    for (size_t k = 0; k < n; k++) {
        /* Collect window samples */
        size_t count = 0;
        for (size_t j = (k >= half ? k - half : 0);
             j <= (k + half < n ? k + half : n - 1) && count < window;
             j++) {
            buf[count++] = input[j];
        }

        /* Simple selection sort to find median (OK for small windows ≤ 9) */
        /* Find median using quickselect-like approach: nth-element */
        for (size_t i = 0; i <= count / 2; i++) {
            size_t min_idx = i;
            for (size_t j = i + 1; j < count; j++) {
                if (buf[j] < buf[min_idx]) min_idx = j;
            }
            double tmp = buf[i];
            buf[i] = buf[min_idx];
            buf[min_idx] = tmp;
        }
        output[k] = buf[count / 2];
    }

    free(buf);
    return 0;
}
