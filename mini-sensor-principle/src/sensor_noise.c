/**
 * @file sensor_noise.c
 * @brief Sensor noise models — Johnson, shot, flicker, Allan variance, NEP, D*
 *
 * Levels: L3 — Mathematical Structures · L4 — Fundamental Laws
 */

#include "sensor_noise.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Physical constants (SI 2019 exact values) */
#define K_B  (1.380649e-23)  /* Boltzmann constant [J/K] */
#define Q_E  (1.602176634e-19) /* elementary charge [C] */

/* ═══════════════════════════════════════════════════════════════════════
 * L4: Johnson-Nyquist (Thermal) Noise
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Johnson-Nyquist noise voltage.
 *
 * Theorem: V_n_rms = √(4·k_B·T·R·Δf)
 *
 * This is a consequence of the fluctuation-dissipation theorem
 * (Callen & Welton, 1951). The available noise power is independent
 * of resistance: P_n_avail = k_B·T·Δf.
 *
 * At T = 290 K (standard IEEE noise temperature):
 *   √(4k_BT) ≈ 126.5 pV/√(Ω·Hz)
 *
 * Example: R=50 Ω, Δf=1 MHz → V_n = 126.5e-12·√(50·1e6) = 0.89 µVrms
 */
double johnson_noise_voltage(double T, double R, double delta_f)
{
    if (T < 0.0 || R < 0.0 || delta_f < 0.0) return 0.0;
    return sqrt(4.0 * K_B * T * R * delta_f);
}

double johnson_noise_power_available(double T, double delta_f)
{
    if (T < 0.0 || delta_f < 0.0) return 0.0;
    return K_B * T * delta_f;
}

double johnson_noise_temperature(double T_physical)
{
    return T_physical;  /* for a passive resistor at equilibrium */
}

/* ═══════════════════════════════════════════════════════════════════════
 * L4: Shot Noise
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Shot noise RMS current.
 *
 * Theorem (Schottky, 1918): I_n_rms = √(2·q·I_DC·Δf)
 *
 * Shot noise arises from the discrete nature of charge carriers crossing
 * a potential barrier independently. Each carrier arrival is a Poisson process.
 *
 * For Poisson processes with rate λ, the variance equals the mean.
 * The noise PSD is S_I(f) = 2·q·I_DC (double-sided), white up to
 * f ≈ 1/(2πτ_transit) where τ_transit is the transit time.
 *
 * Example: I_DC=1 µA, Δf=1 kHz → I_n=√(2·1.6e-19·1e-6·1e3)=17.9 pArms
 */
double shot_noise_current(double I_dc, double delta_f)
{
    if (I_dc < 0.0 || delta_f < 0.0) return 0.0;
    return sqrt(2.0 * Q_E * I_dc * delta_f);
}

int shot_noise_is_significant(double I_dc, double R, double T,
                               double delta_f, double threshold_ratio)
{
    double shot_power = 2.0 * Q_E * I_dc * delta_f * R;   /* I² * R */
    double thermal_power = 4.0 * K_B * T * delta_f * R;    /* (V_n²/R) * R */
    if (thermal_power < 1e-30) return (shot_power > 0.0);
    return (shot_power / thermal_power) > threshold_ratio;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L3: 1/f (Flicker) Noise
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief 1/f noise power spectral density at a given frequency.
 *
 * S(f) = K_f / f^α
 *
 * Origin: conductivity fluctuations due to traps in semiconductors,
 * contact resistance variations, and other slow processes.
 * Hooge's empirical formula: S_R(f)/R² = α_H/(N·f) where N = carrier count.
 */
double flicker_noise_psd(const flicker_noise_spec_t *fn, double freq)
{
    if (!fn || freq <= 0.0) return 0.0;
    return fn->K_f / pow(freq, fn->alpha);
}

/**
 * @brief Integrated RMS noise from 1/f noise over [f_low, f_high].
 *
 * V² = ∫_{f_low}^{f_high} K_f / f^α df
 *
 * For α = 1: V² = K_f·ln(f_high/f_low)
 * For α ≠ 1: V² = K_f·(f_high^{1-α} - f_low^{1-α})/(1-α)
 *
 * This integral diverges if:
 *   - α ≥ 1 and f_low → 0 (infrared catastrophe)
 *   - α ≤ 1 and f_high → ∞
 *
 * In practice, f_low is set by observation time: f_low ≈ 1/T_obs.
 */
double flicker_noise_integrated_rms(const flicker_noise_spec_t *fn)
{
    if (!fn) return 0.0;
    double alpha = fn->alpha;

    if (fabs(alpha - 1.0) < 1e-6) {
        if (fn->f_high <= fn->f_low || fn->f_low <= 0.0) return 0.0;
        return sqrt(fn->K_f * log(fn->f_high / fn->f_low));
    } else {
        double term1 = pow(fn->f_high, 1.0 - alpha);
        double term2 = pow(fn->f_low, 1.0 - alpha);
        double integral = fn->K_f * (term1 - term2) / (1.0 - alpha);
        if (integral < 0.0) integral = 0.0;
        return sqrt(integral);
    }
}

/**
 * @brief Compute the corner frequency where flicker noise PSD equals white noise PSD.
 *
 * K_f / f_c^α = S_white → f_c = (K_f / S_white)^{1/α}
 *
 * Below f_c, 1/f noise dominates; above f_c, white noise dominates.
 * Low f_c is desirable — indicates low 1/f noise.
 */
double flicker_noise_corner_freq(const flicker_noise_spec_t *fn,
                                  double white_noise_psd)
{
    if (!fn || white_noise_psd <= 0.0 || fn->K_f <= 0.0) return 0.0;
    return pow(fn->K_f / white_noise_psd, 1.0 / fn->alpha);
}

/**
 * @brief Estimate α (frequency exponent) from log-log PSD slope.
 *
 * On a log-log plot: log S(f) = log K_f - α·log f
 * The slope is -α, so α = -slope.
 */
double flicker_noise_alpha_from_slope(double log_slope)
{
    return -log_slope;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L3: Total Noise in Sensor Systems
 * ═══════════════════════════════════════════════════════════════════════ */

double noise_rss_combine(const double *noise_rms, size_t n)
{
    if (!noise_rms || n == 0) return 0.0;
    double sum_sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum_sq += noise_rms[i] * noise_rms[i];
    }
    return sqrt(sum_sq);
}

double enob_from_noise(double full_scale_range, double noise_rms)
{
    if (full_scale_range <= 0.0 || noise_rms <= 0.0) return 0.0;
    /* ENOB = log₂(FSR / (noise_rms · √12))
     * Derivation: ideal quantization σ_q = FS/(2^N · √12)
     * Setting σ_q = σ_noise → 2^N = FS/(σ_noise·√12) */
    double bits = full_scale_range / (noise_rms * sqrt(12.0));
    if (bits <= 1.0) return 0.0;
    return log2(bits);
}

double snr_db_from_rms(double signal_rms, double noise_rms)
{
    if (noise_rms <= 0.0) return INFINITY;
    if (signal_rms <= 0.0) return -INFINITY;
    return 20.0 * log10(signal_rms / noise_rms);
}

/* ═══════════════════════════════════════════════════════════════════════
 * L3: Noise Equivalent Power (NEP) & Detectivity (D*)
 * ═══════════════════════════════════════════════════════════════════════ */

double nep_compute(double noise_voltage_rms, double responsivity,
                    double bandwidth)
{
    if (responsivity <= 0.0 || bandwidth <= 0.0) return INFINITY;
    return noise_voltage_rms / (responsivity * sqrt(bandwidth));
}

double detectivity_d_star(double NEP, double area_cm2, double bandwidth)
{
    if (NEP <= 0.0 || bandwidth <= 0.0) return 0.0;
    return sqrt(area_cm2 * bandwidth) / NEP;
}

/**
 * @brief BLIP detectivity limit.
 *
 * D*_BLIP = η·λ / (h·c·√(2·η·Φ_B))   → simplified for ideal
 *         ≈ λ / (2·h·c·√(Φ_B·η))
 *
 * where:
 *   λ = wavelength [m]
 *   h = 6.62607015e-34 J·s (Planck constant)
 *   c = 2.99792458e8 m/s (speed of light)
 *   Φ_B = background photon flux [photons/(s·cm²)]
 *   η = quantum efficiency
 *
 * This is the fundamental detectivity limit for photon detectors
 * when performance is limited by photon noise from the background.
 */
double detectivity_blip_limit(double wavelength_um, double background_flux,
                               double quantum_efficiency)
{
    if (wavelength_um <= 0.0 || background_flux <= 0.0
        || quantum_efficiency <= 0.0 || quantum_efficiency > 1.0) {
        return 0.0;
    }

    const double h = 6.62607015e-34;
    const double c = 2.99792458e8;
    double lambda_m = wavelength_um * 1e-6;

    /* Photon noise variance from background (Poisson): σ² = 2·η·Φ_B·A·t
     * Per unit area, unit bandwidth: σ² = 2·η·Φ_B
     * But the standard expression for ideal photovoltaic detector is:
     * D*_BLIP = η / (h·c/λ · √(2·η·Φ_B)) → simplified below */
    double photon_energy = h * c / lambda_m;
    double noise_per_root_area_bw = photon_energy * sqrt(2.0 * background_flux
                                                          / quantum_efficiency);
    if (noise_per_root_area_bw < 1e-30) return INFINITY;
    return 1.0 / noise_per_root_area_bw;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L3: Allan Variance
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Compute overlapping Allan variance.
 *
 * The overlapping estimator uses all possible clusters of length m,
 * providing better confidence than the non-overlapping estimator.
 *
 * Algorithm:
 *   1. For each cluster length m (τ = m·dt), compute the average
 *      of each cluster: ȳ_k = (1/m)·Σ_{i=k}^{k+m-1} y_i
 *   2. σ²(τ) = (1/[2·(N-2m+1)]) · Σ_{k=0}^{N-2m-1} (ȳ_{k+m} - ȳ_k)²
 *
 * The cluster times τ are logarithmically spaced from dt to N·dt/3
 * to ensure statistical significance (at least 3 independent clusters).
 */
int allan_variance_compute_overlapping(const double *data, size_t N,
                                        double dt, allan_variance_t *av,
                                        size_t max_clusters)
{
    if (!data || N < 4 || dt <= 0.0 || !av || max_clusters < 2) return -1;

    /* Generate logarithmically spaced cluster sizes */
    size_t m_max = N / 3;
    if (m_max < 2) m_max = 2;

    /* Determine number of clusters */
    size_t n_tau = max_clusters;
    if (n_tau > m_max) n_tau = m_max;
    if (n_tau > 100) n_tau = 100;

    av->num_clusters = n_tau;
    av->tau_values = (double *)calloc(n_tau, sizeof(double));
    av->allan_var = (double *)calloc(n_tau, sizeof(double));
    av->allan_dev = (double *)calloc(n_tau, sizeof(double));

    if (!av->tau_values || !av->allan_var || !av->allan_dev) {
        allan_variance_free(av);
        return -1;
    }

    /* Generate cluster sizes: powers of 2 */
    size_t m = 1;
    for (size_t i = 0; i < n_tau; i++) {
        if (m > m_max) m = m_max;
        double tau = m * dt;
        av->tau_values[i] = tau;

        /* Compute overlapping Allan variance for this cluster size */
        size_t num_terms = N - 2 * m + 1;
        if (num_terms < 1) num_terms = 1;

        double sum_sq_diff = 0.0;
        for (size_t k = 0; k < num_terms; k++) {
            /* Compute mean of cluster k and cluster k+m */
            double mean_k = 0.0, mean_km = 0.0;
            for (size_t j = 0; j < m; j++) {
                mean_k += data[k + j];
                mean_km += data[k + m + j];
            }
            mean_k /= (double)m;
            mean_km /= (double)m;
            double diff = mean_km - mean_k;
            sum_sq_diff += diff * diff;
        }

        av->allan_var[i] = sum_sq_diff / (2.0 * (double)num_terms);
        if (av->allan_var[i] > 0.0) {
            av->allan_dev[i] = sqrt(av->allan_var[i]);
        } else {
            av->allan_dev[i] = 0.0;
        }

        m *= 2;
    }

    /* Find minimum Allan deviation (bias instability floor) */
    av->min_allan_dev = INFINITY;
    av->tau_at_min = 0.0;
    for (size_t i = 0; i < n_tau; i++) {
        if (av->allan_dev[i] > 0.0 && av->allan_dev[i] < av->min_allan_dev) {
            av->min_allan_dev = av->allan_dev[i];
            av->tau_at_min = av->tau_values[i];
        }
    }
    if (av->min_allan_dev >= INFINITY) av->min_allan_dev = 0.0;

    /* Identify noise coefficients */
    allan_noise_coefficients_identify(av);

    return 0;
}

/**
 * @brief Identify noise coefficients by analyzing log-log Allan deviation slope.
 *
 * On a log-log plot of σ(τ) vs τ:
 *   slope = -1/2  →  white noise (angle random walk) [dominant at small τ]
 *   slope = 0     →  flicker noise (bias instability) [dominant at medium τ]
 *   slope = +1/2  →  random walk (rate random walk)    [dominant at large τ]
 *   slope = +1    →  drift / ramp
 *   slope = -1    →  quantization noise
 *
 * The coefficients are extracted by fitting to each segment.
 */
int allan_noise_coefficients_identify(allan_variance_t *av)
{
    if (!av || av->num_clusters < 3) return -1;

    /* Reset coefficients */
    av->angle_random_walk = 0.0;
    av->bias_instability = av->min_allan_dev;
    av->rate_random_walk = 0.0;
    av->quantization_noise = 0.0;
    av->drift_ramp = 0.0;

    size_t n = av->num_clusters;

    /* Estimate ARW from smallest τ: σ(τ) ≈ ARW / √τ → ARW ≈ σ(τ₁)·√τ₁ */
    if (av->tau_values[0] > 0.0) {
        av->angle_random_walk = av->allan_dev[0] * sqrt(av->tau_values[0]);
    }

    /* Estimate RRW from largest τ: σ(τ) ≈ RRW·√τ → RRW ≈ σ(τ_n)/√τ_n */
    if (av->tau_values[n-1] > 0.0) {
        av->rate_random_walk = av->allan_dev[n-1] / sqrt(av->tau_values[n-1]);
    }

    /* Estimate quantization noise from the τ⁻¹ slope region (smallest τ).
     * σ_Q ≈ Q / τ → Q ≈ σ(τ₁)·τ₁ (if slope ≈ -1) */
    if (n >= 2 && av->tau_values[0] > 0.0) {
        double slope = log(av->allan_dev[1] / av->allan_dev[0])
                       / log(av->tau_values[1] / av->tau_values[0]);
        if (slope < -0.75) {
            av->quantization_noise = av->allan_dev[0] * av->tau_values[0];
        }
    }

    /* Estimate drift ramp from large τ slope ≈ +1 */
    if (n >= 2) {
        double slope = log(av->allan_dev[n-1] / av->allan_dev[n-2])
                       / log(av->tau_values[n-1] / av->tau_values[n-2]);
        if (slope > 0.75) {
            av->drift_ramp = av->allan_dev[n-1] / av->tau_values[n-1];
        }
    }

    return 0;
}

void allan_variance_free(allan_variance_t *av)
{
    if (av) {
        free(av->tau_values);
        free(av->allan_var);
        free(av->allan_dev);
        memset(av, 0, sizeof(allan_variance_t));
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * L3: Detection Metrics
 * ═══════════════════════════════════════════════════════════════════════ */

double minimum_detectable_signal(double noise_rms, double k_detection)
{
    if (noise_rms < 0.0 || k_detection < 0.0) return 0.0;
    return k_detection * noise_rms;
}

double resolution_from_noise(double noise_rms, double sensitivity)
{
    if (fabs(sensitivity) < 1e-30) return INFINITY;
    /* 3-sigma rule for 99.7% confidence */
    return 3.0 * noise_rms / fabs(sensitivity);
}

/**
 * @brief Rose criterion for human/electronic detection.
 *
 * The Rose criterion states that a feature can be detected with ~50%
 * probability when SNR ≥ k, where k depends on the observer/model:
 *   k = 3-5 for human observers under ideal conditions
 *   k ≈ 7-10 for automated detection with low false-alarm rate
 *
 * This function checks if SNR meets the threshold for given probability.
 *
 * Reference: Rose, A., "The sensitivity performance of the human eye
 * on an absolute scale", J. Opt. Soc. Am., 38(2):196-208, 1948
 */
int rose_criterion_check(double snr, double probability_threshold)
{
    if (probability_threshold >= 1.0) return (snr >= 7.0);
    if (probability_threshold >= 0.5) return (snr >= 5.0);
    return (snr >= 3.0);
}

/* ═══════════════════════════════════════════════════════════════════════
 * L3: Noise Figure and Noise Temperature
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Cascaded noise factor (Friis formula, 1944).
 *
 * F_total = F₁ + (F₂-1)/G₁ + (F₃-1)/(G₁·G₂) + ...
 *
 * This shows that the first stage's noise factor dominates:
 * if G₁ is large, subsequent contributions are divided down.
 *
 * Example: LNA with F₁=1.5 (1.76 dB), G₁=100 (20 dB)
 *          Mixer with F₂=10 (10 dB)
 *          F_total = 1.5 + (10-1)/100 = 1.59 (2.01 dB)
 *          → The mixer adds only 0.25 dB to the noise figure!
 */
double noise_factor_cascade(const double *noise_factors,
                             const double *gains, size_t n_stages)
{
    if (!noise_factors || !gains || n_stages == 0) return 1.0;

    double F_total = noise_factors[0];
    double G_product = gains[0];

    for (size_t i = 1; i < n_stages; i++) {
        if (G_product < 1e-30) break;
        F_total += (noise_factors[i] - 1.0) / G_product;
        G_product *= gains[i];
    }

    return F_total;
}

double noise_figure_db_from_factor(double noise_factor)
{
    if (noise_factor <= 0.0) return 0.0;
    return 10.0 * log10(noise_factor);
}

double noise_factor_from_figure_db(double noise_figure_db)
{
    return pow(10.0, noise_figure_db / 10.0);
}

/**
 * @brief Noise temperature from noise factor.
 *
 * T_e = T₀ · (F - 1)
 * where T₀ = 290 K (standard IEEE reference temperature).
 *
 * This represents the equivalent input noise temperature:
 * an ideal (noiseless) system with a resistor at T_e at its input
 * would produce the same output noise as the actual system.
 */
double noise_temperature_from_factor(double noise_factor, double T_ref)
{
    return T_ref * (noise_factor - 1.0);
}

double noise_factor_from_temperature(double T_noise, double T_ref)
{
    if (T_ref <= 0.0) return 1.0;
    return 1.0 + T_noise / T_ref;
}
