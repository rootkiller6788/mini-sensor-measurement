/**
 * @file sensor_signal_cond.h
 * @brief Sensor signal conditioning — amplification, filtering, linearization, excitation
 *
 * Level: L2 — Core Concepts · L5 — Algorithms/Methods
 *
 * Signal conditioning transforms raw sensor output into a form suitable
 * for analog-to-digital conversion (ADC) or direct use.
 *
 * Core Functions:
 *   1. Amplification: raise mV-level signals to ADC input range
 *   2. Filtering: remove noise, prevent aliasing
 *   3. Linearization: compensate for sensor nonlinearity in analog or digital domain
 *   4. Excitation: provide stable voltage/current to passive sensors
 *   5. Isolation: galvanically separate sensor from measurement system
 *   6. Cold junction compensation: for thermocouples
 *   7. Current loop: 4-20 mA industrial standard
 *
 * The 4-20 mA Current Loop Standard (ISA-50.1):
 *   - 4 mA = zero/live-zero (distinguishes from broken wire = 0 mA)
 *   - 20 mA = full scale
 *   - 2-wire: power + signal on same pair (loop-powered)
 *   - 3-wire: separate power ground
 *   - 4-wire: separate power and signal pairs
 *
 *   Supply voltage: typically 24 VDC
 *   Maximum loop resistance: R_max = (V_supply - V_min_sensor) / 0.020
 *   e.g., V_supply=24V, V_min=8V → R_max = 800 Ω
 *
 * References:
 *   - Pallas-Areny & Webster, "Sensors and Signal Conditioning", 2nd Ed, 2000
 *   - Kester, "Sensor Signal Conditioning", Analog Devices, 1999
 *   - ISA-50.1 "Compatibility of Analog Signals for Electronic Industrial
 *     Process Instruments"
 */

#ifndef SENSOR_SIGNAL_COND_H
#define SENSOR_SIGNAL_COND_H

#include <stddef.h>
#include <stdint.h>

/* ──── L2: Voltage Amplification ──── */

/**
 * @brief Compute non-inverting amplifier gain.
 *
 *   G = 1 + Rf/Rg
 *
 * Input impedance: very high (~GΩ for FET-input op-amp)
 * Output impedance: very low (~mΩ at DC)
 *
 * Preferred for high-impedance sensors (piezoelectric, pH, ISFET).
 */
double amp_non_inverting_gain(double Rf, double Rg);

/**
 * @brief Compute inverting amplifier gain.
 *
 *   G = -Rf/Rin
 *
 * Input impedance = Rin (typically low, loads sensor).
 * Preferred when virtual ground summing of multiple sensor inputs is needed.
 */
double amp_inverting_gain(double Rf, double Rin);

/**
 * @brief Compute differential amplifier gain.
 *
 *   V_out = (Rf/Rin) · (V⁺ - V⁻)   [for R1=Rin, R2=Rf on both sides]
 *
 * CMRR depends on resistor matching:
 *   CMRR ≈ (1 + G) / (4 · ΔR/R)   where ΔR/R is resistor tolerance
 *
 * For G=1 with 0.1% resistors: CMRR ≈ 500 (~54 dB)
 * For G=100 with 0.1% resistors: CMRR ≈ 25000 (~88 dB)
 */
double amp_differential_gain(double Rf, double Rin);
double amp_differential_cmrr_estimate(double gain, double resistor_tolerance);

/**
 * @brief Instrumentation amplifier gain (3-op-amp topology).
 *
 * Stage 1 (input buffers): G1 = 1 + 2·Rf/Rg
 * Stage 2 (difference amp): G2 = R2/R1
 * Total: G = (1 + 2·Rf/Rg) · (R2/R1)
 *
 * Key advantage: very high input impedance, gain set by single resistor Rg.
 * CMRR > 100 dB achievable with integrated resistor matching.
 */
double amp_instrumentation_gain(double Rf, double Rg, double R2, double R1);

/* ──── L5: Anti-Alias Filter Design ──── */

/**
 * @brief Compute first-order RC low-pass filter cutoff frequency.
 *
 *   f_c = 1/(2π·R·C)
 *
 * For anti-alias filtering before ADC with sampling rate f_s:
 *   f_c ≤ f_s / (2 · k)  where k ≥ 5-10 for reasonable alias rejection
 *
 * At f = f_s/2 (Nyquist), a 1st-order filter attenuates by:
 *   |H| = 1/√(1 + (f_s/(2f_c))²) = 1/√(1 + k²) ≈ 1/k
 *
 * For 10-bit ADC (60 dB SNR), k=31 gives ~30 dB alias rejection → may need
 * 2nd-order or higher for demanding applications.
 */
double anti_alias_fc_from_fs(double fs_Hz, double oversampling_factor);
double anti_alias_attenuation_at_nyquist(double fc_Hz, double fs_Hz, int order);

/**
 * @brief Compute RC filter component values for target cutoff.
 *
 * Chooses standard E12/E24 resistor values and computes required capacitor.
 */
int anti_alias_rc_values(double fc_target, double *R_out, double *C_out);

/**
 * @brief Compute required filter order for given alias rejection.
 *
 * For Butterworth low-pass with cutoff f_c and stopband at f_s/2:
 *   Attenuation(dB) = 10·log₁₀(1 + (f_s/(2f_c))^{2N})
 *   Solve for N: N ≥ log₁₀(10^{A/10} - 1) / (2·log₁₀(f_s/(2f_c)))
 *
 * @param fs_hz       sampling frequency [Hz]
 * @param fc_hz       filter cutoff frequency [Hz]
 * @param attenuation_db required stopband attenuation at Nyquist
 * @return required filter order (rounded up to integer)
 */
int anti_alias_filter_order_required(double fs_hz, double fc_hz,
                                      double attenuation_db);

/* ──── L2: 4-20 mA Current Loop ──── */

/**
 * @brief Convert sensor value to 4-20 mA loop current.
 *
 *   I_loop = 4 mA + (x - x_min) / (x_max - x_min) · 16 mA
 *
 * Live zero (4 mA) provides:
 *   - Broken wire detection (I < 3.6 mA or I > 21 mA per NAMUR NE43)
 *   - Power for 2-wire loop-powered transmitters
 *
 * NAMUR NE43 fault thresholds:
 *   I < 3.6 mA  → under-range / broken wire
 *   3.8–20.5 mA  → valid measurement range
 *   I > 21.0 mA → over-range / short circuit
 */
double current_loop_sensor_to_current(double sensor_value,
                                       double x_min, double x_max);

/**
 * @brief Convert 4-20 mA loop current back to sensor value.
 */
double current_loop_current_to_sensor(double I_loop,
                                       double x_min, double x_max);

/**
 * @brief Compute sensing resistor value for current-to-voltage conversion.
 *
 *   R_sense = V_adc_max / 0.020   [for 20 mA → full ADC range]
 *
 * e.g., V_adc_max=3.3V → R_sense = 165 Ω
 *       V_adc_max=5.0V → R_sense = 250 Ω (standard)
 *
 * Power dissipated: P = I²·R = 0.02²·R_sense ≤ 100 mW typically.
 */
double current_loop_sense_resistor(double V_adc_fullscale);

/**
 * @brief Compute maximum loop cable length given wire gauge.
 *
 * R_max = (V_supply - V_min_transmitter - V_sense) / 0.020
 * L_max = R_max / (2 · R_per_meter)   [round-trip resistance]
 *
 * e.g., AWG 24 (84 Ω/km), V_supply=24V, V_min=8V, V_sense=5V (250Ω)
 *       R_max = (24-8-5)/0.02 = 550 Ω → L_max = 550/(2·0.084) = 3274 m
 */
double current_loop_max_distance(double V_supply, double V_min_transmitter,
                                  double R_sense, double R_per_meter);

/* ──── L5: Sensor Linearization ──── */

/**
 * @brief Square root linearization (e.g., orifice plate flow sensor).
 *
 * For differential-pressure flow sensors:
 *   Q = K · √(ΔP)
 *
 * The square-root extraction can be done in analog (log-antilog circuit)
 * or digitally (lookup table, polynomial expansion).
 */
double linearize_sqrt_law(double differential_pressure, double K);

/**
 * @brief Logarithmic linearization (e.g., thermistor with op-amp feedback).
 *
 * Using a logarithmic amplifier:
 *   V_out = -V_T · ln(I_in / I_s)
 *   where V_T = kT/q ≈ 26 mV at 300 K
 *
 * Useful for photodiodes (log compression of wide dynamic range signals).
 */
double linearize_log_compression(double I_signal, double I_reference, double V_T);

/**
 * @brief Bridge linearization via op-amp feedback.
 *
 * Placing the active bridge arm in an op-amp feedback loop makes
 * V_out linearly proportional to ΔR (instead of the nonlinear bridge equation).
 *
 * With ideal op-amp: V_out = -V_ref · (ΔR / R₀)
 *
 * This technique is used in precision RTD and strain gauge signal conditioners
 * to avoid digital linearization.
 */
double linearize_bridge_active(double delta_R, double R0, double V_ref);

/* ──── L2: Sensor Excitation ──── */

/**
 * @brief Compute excitation voltage stability requirement.
 *
 * For a sensor whose output is proportional to excitation:
 *   Δy/y = ΔV_ex/V_ex   (sensitivity is a ratio metric)
 *
 * If measurement target accuracy is ε_total, then:
 *   ε_excitation ≤ ε_total / 3  (assign 1/3 of error budget to excitation)
 *
 * For 0.1% total accuracy → V_ex must be stable to ±0.033%.
 */
double excitation_stability_required(double total_accuracy_percent);

/**
 * @brief Constant current source for RTD or thermistor excitation.
 *
 * V_out = I_ex · R_sensor
 *
 * Self-heating must be considered:
 *   P = I_ex² · R_sensor  →  ΔT = P · θ_ja  [thermal resistance]
 *
 * For Pt100 with 1 mA excitation:
 *   P = 0.001² · 100 = 100 µW
 *   ΔT ≈ 100µW · 50 K/W = 5 mK  → acceptable
 *
 * For 10 kΩ thermistor with 1 mA:
 *   P = 0.001² · 10000 = 10 mW
 *   ΔT ≈ 0.01 · 200 K/W = 2 K  → unacceptable
 *   → Use 100 µA instead: P = 0.1 µW, ΔT = 0.02 K
 */
double excitation_current_max_for_self_heating(double R_sensor,
                                                double thermal_resistance,
                                                double max_temp_rise);

/**
 * @brief AC excitation for capacitive/inductive sensors.
 *
 * Capacitive reactance: X_C = 1/(2π·f·C)
 * Inductive reactance: X_L = 2π·f·L
 *
 * Carrier frequency f_carrier should be >> signal bandwidth
 * but << parasitic effects become dominant.
 *
 * Synchronous demodulation recovers amplitude and phase.
 */
double ac_excitation_reactance(double f_hz, double C_or_L, int is_capacitive);

/* ──── L5: Digital Filtering for Sensor Noise ──── */

/**
 * @brief Simple moving average filter.
 *
 * y[k] = (1/N) · Σ_{i=0}^{N-1} x[k-i]
 *
 * Frequency response: H(f) = sin(πfN/fs) / [N·sin(πf/fs)] · e^{-jπf(N-1)/fs}
 * First null at f = fs/N.
 *
 * Reduces white noise by √N.
 * Effective for removing periodic interference at multiples of fs/N.
 */
int signal_moving_average(const double *input, double *output,
                           size_t n_samples, size_t window_size);

/**
 * @brief Exponential moving average (first-order IIR low-pass).
 *
 * y[k] = α·x[k] + (1-α)·y[k-1]
 *
 * Time constant: τ = dt / (1-α) ≈ dt/α for small α
 * Cutoff frequency: f_c = α / (2π·dt) for small α
 *
 * This is the digital equivalent of an RC filter.
 * Memory efficient: only stores previous output.
 */
double signal_exponential_filter(double x_new, double y_prev, double alpha);

/**
 * @brief Median filter for impulse noise (salt-and-pepper) removal.
 *
 * Replaces each sample with the median of its N neighbors.
 * Preserves edges better than moving average.
 * Effective for removing occasional outlier spikes from EMI or arcing.
 *
 * @param input  input signal
 * @param output output signal
 * @param n      number of samples
 * @param window median window size (odd, ≥ 3)
 */
int signal_median_filter(const double *input, double *output,
                          size_t n, size_t window);

#endif /* SENSOR_SIGNAL_COND_H */
