/**
 * @file sensor_noise.h
 * @brief Sensor noise models — Johnson, shot, flicker, Allan variance, NEP, D*
 *
 * Level: L3 — Mathematical Structures · L4 — Fundamental Laws
 *
 * Fundamental Noise Theorems:
 *
 * 1. Johnson-Nyquist (Thermal) Noise (Nyquist, 1928):
 *    V_n_rms = √(4·k_B·T·R·Δf)
 *    P_n = k_B·T·Δf   [available noise power, independent of R]
 *
 * 2. Shot Noise (Schottky, 1918):
 *    I_n_rms = √(2·q·I_DC·Δf)
 *    Requires DC current flow across a potential barrier.
 *
 * 3. 1/f (Flicker) Noise (Hooge, 1969):
 *    S_V(f) = K / f^α   with α ≈ 0.8–1.3
 *    Arises from conductance fluctuations. Corner frequency = flicker/white intercept.
 *
 * 4. Allan Variance (Allan, 1966):
 *    σ²(τ) = 1/2 ⟨(ȳ_{k+1} - ȳ_k)²⟩
 *    Reveals noise processes by their τ-dependence in log-log plot.
 *
 * 5. Noise Equivalent Power (NEP):
 *    NEP = V_n_rms / R_V   [W/√Hz]
 *    Minimum detectable power for SNR=1 in 1 Hz bandwidth.
 *
 * 6. Specific Detectivity D*:
 *    D* = √(A·Δf) / NEP   [cm·√Hz/W]  (Jones, 1959)
 *    Normalized detectivity allowing comparison across detector sizes.
 *
 * References:
 *   - Nyquist, Phys. Rev. 32:110-113, 1928
 *   - Schottky, Ann. Phys. 362:541-567, 1918
 *   - Allan, Proc. IEEE 54(2):221-230, 1966
 *   - Jones, Proc. IRE 47(9):1495-1502, 1959
 *   - Hooge, IEEE Trans. Electron Devices, 41(11):1926-1935, 1994
 */

#ifndef SENSOR_NOISE_H
#define SENSOR_NOISE_H

#include "sensor_defs.h"
#include <stddef.h>

/* ──── L4: Johnson-Nyquist Noise ──── */

/**
 * @brief Compute Johnson-Nyquist thermal noise voltage.
 *
 * Theorem (Nyquist, 1928):
 *   V_n_rms = √(4·k_B·T·R·Δf)
 *
 * k_B = 1.380649×10⁻²³ J/K (Boltzmann constant, exact per SI 2019)
 *
 * At T = 300 K (27°C):
 *   √(4k_BT) ≈ 4.07×10⁻¹¹ √(W/Hz)² = 4.07 nV/√Hz into 1 Ω
 *   Into 1 kΩ: 129 nV/√Hz
 *   Into 1 MΩ: 4.07 µV/√Hz
 *
 * @param T temperature [K]
 * @param R resistance [Ω]
 * @param delta_f bandwidth [Hz]
 * @return V_n_rms [V]
 */
double johnson_noise_voltage(double T, double R, double delta_f);

/**
 * @brief Available noise power from a resistor.
 *
 * P_n_available = k_B·T·Δf   [W]
 *
 * This is the maximum noise power that can be delivered to a matched load.
 * Independent of resistance value — a fundamental result.
 */
double johnson_noise_power_available(double T, double delta_f);

/**
 * @brief Noise figure contribution of a resistor.
 *
 * When a resistor is used in series with a signal source, it adds noise.
 * The noise temperature contribution is:
 *   T_noise = T_physical (for a passive resistor at thermal equilibrium)
 *
 * Excess noise ratio: ENR = (T_noise - T₀) / T₀ where T₀ = 290 K
 */
double johnson_noise_temperature(double T_physical);

/* ──── L4: Shot Noise ──── */

/**
 * @brief Compute shot noise RMS current.
 *
 * Theorem (Schottky, 1918):
 *   I_n_rms = √(2·q·I_DC·Δf)
 *
 * q = 1.602176634×10⁻¹⁹ C (elementary charge, exact per SI 2019)
 *
 * Shot noise is white (constant PSD) up to f ≈ 1/(2π·τ_transit).
 * Present in: p-n junctions, Schottky diodes, photodiodes, vacuum tubes.
 * Absent in: metallic conductors (electrons are not independent).
 *
 * @param I_dc DC current [A]
 * @param delta_f bandwidth [Hz]
 * @return I_n_rms [A]
 */
double shot_noise_current(double I_dc, double delta_f);

/**
 * @brief Detect shot noise in a measurement by checking if current
 *        exceeds pure thermal noise prediction.
 *
 * If I_measured_noise² >> 4k_BTΔf/R, shot noise is significant.
 */
int shot_noise_is_significant(double I_dc, double R, double T,
                               double delta_f, double threshold_ratio);

/* ──── L3: 1/f (Flicker) Noise ──── */

/**
 * @brief 1/f noise PSD and integrated noise.
 *
 * Power spectral density: S_V(f) = K_f / f^α
 * where K_f is the flicker coefficient and α ≈ 1.
 *
 * Integrated noise from f_low to f_high:
 *   V_n_rms² = K_f · ln(f_high / f_low)    [for α = 1]
 *   V_n_rms² = K_f · (f_high^{1-α} - f_low^{1-α}) / (1-α)  [for α ≠ 1]
 *
 * Corner frequency f_c: where 1/f PSD equals white noise PSD.
 *   S_white = S_flicker(f_c) → f_c = K_f / S_white  [for α=1]
 */
typedef struct {
    double K_f;     /* flicker noise coefficient [V²/Hz or A²/Hz at 1 Hz] */
    double alpha;   /* frequency exponent, typically 0.8–1.3 */
    double f_low;   /* lower integration limit [Hz] */
    double f_high;  /* upper integration limit [Hz] */
} flicker_noise_spec_t;

double flicker_noise_psd(const flicker_noise_spec_t *fn, double freq);
double flicker_noise_integrated_rms(const flicker_noise_spec_t *fn);
double flicker_noise_corner_freq(const flicker_noise_spec_t *fn,
                                  double white_noise_psd);
double flicker_noise_alpha_from_slope(double log_slope);

/* ──── L3: Total Noise in Sensor Systems ──── */

/**
 * @brief Combine uncorrelated noise sources (root-sum-square).
 *
 * For N independent noise sources with RMS values v_i:
 *   V_total_rms = √(Σ v_i²)
 *
 * This follows from: E[(Σ X_i)²] = Σ E[X_i²] for independent zero-mean X_i.
 */
double noise_rss_combine(const double *noise_rms, size_t n);

/**
 * @brief Compute effective number of bits (ENOB) limited by noise.
 *
 * ENOB = log₂(FSR / (σ_noise · √12))
 *
 * where FSR = full-scale range, σ_noise = RMS noise.
 * Derivation: quantization noise σ_q = FS/(2^N · √12).
 * Setting σ_q = σ_noise yields the effective bits.
 *
 * @param full_scale_range peak-to-peak signal range
 * @param noise_rms RMS noise voltage
 * @return ENOB [bits]
 */
double enob_from_noise(double full_scale_range, double noise_rms);

/**
 * @brief Compute signal-to-noise ratio from signal and noise RMS.
 *
 * SNR_dB = 20 · log₁₀(V_signal_rms / V_noise_rms)
 */
double snr_db_from_rms(double signal_rms, double noise_rms);

/* ──── L3: Noise Equivalent Power & Detectivity ──── */

/**
 * @brief Compute Noise Equivalent Power (NEP).
 *
 * NEP = V_noise_rms / R_V   [W] or [W/√Hz]
 *
 * where R_V [V/W] is the responsivity (voltage output per watt of input).
 * NEP is the minimum detectable optical/radiant power for SNR = 1.
 *
 * Normalized NEP per √Hz:
 *   NEP_norm = NEP / √Δf   [W/√Hz]
 */
double nep_compute(double noise_voltage_rms, double responsivity,
                    double bandwidth);

/**
 * @brief Compute Specific Detectivity D* (D-star).
 *
 * D* = √(A_detector · Δf) / NEP   [cm·√Hz/W]
 *      = R_V · √(A) / V_n_rms  [normalized form]
 *
 * D* allows comparison of detectors with different areas.
 * Higher D* = better detector.
 *
 * Typical D* values (500 K blackbody, 1 kHz):
 *   - Thermal (bolometer):  ~10⁸
 *   - PbS (photoconductive): ~10¹¹
 *   - InSb (photovoltaic):  ~10¹¹
 *   - HgCdTe (photovoltaic): ~10¹⁰–10¹¹
 *   - Ideal BLIP: ~10¹²
 */
double detectivity_d_star(double NEP, double area_cm2, double bandwidth);

/**
 * @brief Background-Limited Infrared Photodetector (BLIP) D* limit.
 *
 * D*_BLIP = λ / (2·h·c·√(Q_B))   where Q_B is the background photon flux.
 *
 * This is the fundamental limit set by photon noise from the background.
 */
double detectivity_blip_limit(double wavelength_um, double background_flux,
                               double quantum_efficiency);

/* ──── L3: Allan Variance Computation ──── */

/**
 * @brief Compute overlapping Allan variance from a time series.
 *
 * Input: equally-spaced time series y[0..N-1] with sampling period dt.
 *
 * Overlapping estimator (preferred, more efficient):
 *   σ²(τ) = (1 / [2·m²·(N-2m+1)]) · Σ_{j=0}^{N-2m} [Σ_{i=j}^{j+m-1} (y_{i+m} - y_i)]²
 *   where τ = m·dt is the cluster time.
 *
 * The function computes σ²(τ) for logarithmically-spaced τ values
 * from dt to N·dt/3 (for statistical significance).
 *
 * Noise slope identification (log-log Allan deviation plot):
 *   white noise:       slope = -1/2  (σ ∝ τ⁻⁰·⁵)
 *   flicker noise:     slope =  0    (σ ∝ τ⁰, horizontal)
 *   random walk:       slope = +1/2  (σ ∝ τ⁰·⁵)
 *   quantization:      slope = -1    (σ ∝ τ⁻¹)
 *   drift/ramp:        slope = +1    (σ ∝ τ¹)
 *
 * References:
 *   - IEEE Std 952-1997, "IEEE Standard Specification Format Guide
 *     and Test Procedure for Single-Axis Interferometric Fiber Optic Gyros"
 *   - Allan, D.W., "Statistics of Atomic Frequency Standards",
 *     Proc. IEEE, 54(2):221-230, 1966
 *
 * @param data      input time series [input_unit]
 * @param N         number of data points
 * @param dt        sampling period [s]
 * @param av        output: populated allan_variance_t structure
 * @param max_clusters maximum number of τ values to compute
 * @return 0 on success, -1 on error
 */
int allan_variance_compute_overlapping(const double *data, size_t N,
                                        double dt, allan_variance_t *av,
                                        size_t max_clusters);

/**
 * @brief Identify noise coefficients from Allan variance data.
 *
 * Fits power-law model: σ²(τ) = Σ h_i · τ^{μ_i}
 * where μ_i ∈ {-2, -1, 0, 1, 2} correspond to noise types.
 *
 * Uses log-log linear regression on segments of the Allan deviation curve.
 */
int allan_noise_coefficients_identify(allan_variance_t *av);

void allan_variance_free(allan_variance_t *av);

/* ──── L3: Signal-to-Noise and Detection Metrics ──── */

/**
 * @brief Minimum detectable signal (MDS) given noise floor.
 *
 * MDS = k_detection · σ_noise
 * where k_detection depends on required detection probability
 * (k=3 gives ~99.7% confidence for Gaussian noise, k=6 for low false-alarm rate).
 */
double minimum_detectable_signal(double noise_rms, double k_detection);

/**
 * @brief Compute sensor resolution from noise floor.
 *
 * Resolution = 3 · σ_noise / sensitivity
 * (3-sigma rule: the smallest input change distinguishable from noise with 99.7% confidence)
 */
double resolution_from_noise(double noise_rms, double sensitivity);

/**
 * @brief Rose criterion for imaging sensors (SNR threshold).
 *
 * For human observers to detect a feature with 50% probability:
 * SNR > 5 (Rose criterion, also known as "Rose model")
 *
 * For 100% probability: SNR > 7-10
 *
 * Reference: Rose, "Vision: Human and Electronic", Plenum Press, 1973
 */
int rose_criterion_check(double snr, double probability_threshold);

/* ──── L3: Noise Figure and Noise Temperature ──── */

/**
 * @brief Noise factor and noise figure for cascaded stages (Friis formula).
 *
 * Friis formula for cascaded stages:
 *   F_total = F₁ + (F₂-1)/G₁ + (F₃-1)/(G₁·G₂) + ...
 *
 * where F_i = linear noise factor, G_i = linear gain of stage i.
 *
 * Noise figure NF = 10·log₁₀(F)  [dB]
 *
 * This is critical for sensor preamplifier design: the first stage's
 * noise figure dominates the cascade if G₁ is sufficiently large.
 *
 * References:
 *   - Friis, Proc. IRE, 32(7):419-422, 1944
 */
double noise_factor_cascade(const double *noise_factors,
                             const double *gains, size_t n_stages);
double noise_figure_db_from_factor(double noise_factor);
double noise_factor_from_figure_db(double noise_figure_db);
double noise_temperature_from_factor(double noise_factor, double T_ref);
double noise_factor_from_temperature(double T_noise, double T_ref);

#endif /* SENSOR_NOISE_H */
