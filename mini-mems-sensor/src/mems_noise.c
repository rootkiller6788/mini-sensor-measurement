/**
 * @file mems_noise.c
 * @brief L3-L4 MEMS Noise Analysis — Thermo-mechanical noise limit,
 *        Allan variance computation, power spectral density estimation,
 *        noise source identification and simulation
 *
 * Reference: IEEE Std 1139-2008 (Allan Variance)
 *            IEEE Std 952-2020 (Inertial Sensor Noise)
 *            Allan, "Statistics of Atomic Frequency Standards" (1966)
 *            Gabrielson, "Mechanical-Thermal Noise in MEMS" (1993)
 *            Kasdin, "Discrete Simulation of Colored Noise" (1995)
 */

#include "mems_noise.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* Boltzmann constant [J/K] */
#define K_B 1.380649e-23
/* Electron charge [C] */
#define Q_E 1.602176634e-19

/* ─── L3: Noise Process Name ─────────────────────────────────────────────── */

const char *mems_noise_process_name(MemsNoiseProcess np) {
    static const char *names[] = {
        [MEMS_NOISE_WHITE]          = "White Noise (Angle/Velocity RW)",
        [MEMS_NOISE_FLICKER_FM]     = "Flicker FM (Bias Instability)",
        [MEMS_NOISE_FLICKER_PM]     = "Flicker Phase Noise",
        [MEMS_NOISE_RANDOM_WALK_FM] = "Random Walk FM (Rate RW)",
        [MEMS_NOISE_QUANTIZATION]   = "Quantization Noise",
        [MEMS_NOISE_RATE_RAMP]      = "Rate Ramp / Drift",
        [MEMS_NOISE_SINUSOIDAL]     = "Sinusoidal Noise",
        [MEMS_NOISE_WHITE_PM]       = "White Phase Noise"
    };
    if (np <= MEMS_NOISE_WHITE_PM) return names[np];
    return "Unknown";
}

/* ─── L4: Thermo-Mechanical Noise ────────────────────────────────────────────
 *
 * The fundamental noise floor of MEMS inertial sensors arises from
 * Brownian motion of the proof mass due to molecular collisions with
 * the surrounding gas.
 *
 * The fluctuation-dissipation theorem (Callen & Welton, 1951) relates
 * the damping coefficient to the noise force PSD:
 *   S_F = 4·k_B·T·c   [N²/Hz]
 *
 * Dividing by mass gives acceleration noise:
 *   S_a = 4·k_B·T·c / m² = (4·k_B·T·ω_n) / (m·Q)  since c = m·ω_n/Q
 *
 * NEA = sqrt(S_a) = sqrt(4·k_B·T·ω_n / (m·Q))  [m/s²/√Hz]
 *
 * For gyroscopes, the equivalent rate noise (NER):
 *   NER = NEA / (ω_d·A_d) = NEA / v_drive  [rad/s/√Hz]
 *
 * Key insight: larger mass, higher Q, and lower temperature reduce NEA.
 * A 1μg proof mass at room temperature with Q=10⁴ has NEA ≈ 10 μg/√Hz.
 */

double mems_thermomech_nea(double mass, double natural_freq_Hz,
                            double q_factor, double temperature_k) {
    if (mass <= DBL_EPSILON || natural_freq_Hz <= 0.0 ||
        q_factor <= 0.0 || temperature_k <= 0.0) return 0.0;

    double omega_n = 2.0 * M_PI * natural_freq_Hz;
    /* S_a = 4·k_B·T·ω_n / (m·Q) */
    double s_a = 4.0 * K_B * temperature_k * omega_n / (mass * q_factor);
    if (s_a < 0.0) return 0.0;
    return sqrt(s_a);
}

double mems_thermomech_ner(double nea, double v_drive, double sf_norm) {
    if (fabs(v_drive) < DBL_EPSILON || sf_norm <= 0.0) return INFINITY;
    return nea / (v_drive * sf_norm);
}

/* ─── L4: Johnson-Nyquist Noise ────────────────────────────────────────────
 *
 * Johnson (thermal) noise in resistive elements:
 *   v_n_rms = sqrt(4·k_B·T·R·Δf)
 *   PSD: S_v = 4·k_B·T·R  [V²/Hz]
 *
 * This noise source is present in MEMS readout electronics (transimpedance
 * amplifiers, resistive bridges, TIA feedback resistors).
 */

double mems_johnson_noise_density(double resistance, double temperature_k) {
    if (resistance < 0.0 || temperature_k < 0.0) return 0.0;
    return sqrt(4.0 * K_B * temperature_k * resistance);
}

double mems_johnson_noise_rms(double resistance, double temperature_k,
                               double bandwidth_hz) {
    if (resistance < 0.0 || bandwidth_hz < 0.0) return 0.0;
    return sqrt(4.0 * K_B * temperature_k * resistance * bandwidth_hz);
}

/* ─── L4: Shot Noise ──────────────────────────────────────────────────────
 *
 * Shot noise arises from discrete charge carriers crossing a barrier.
 * It is Poisson-distributed and has white spectrum:
 *   i_n_density = sqrt(2·q_e·I_DC)  [A/√Hz]
 *
 * Relevant for MEMS tunneling sensors and photodetectors.
 */

double mems_shot_noise_density(double dc_current) {
    if (dc_current < 0.0) return 0.0;
    return sqrt(2.0 * Q_E * dc_current);
}

/* ─── L3: Simple PRNG — xorshift64* ───────────────────────────────────────── */

double mems_rand64(uint64_t *seed) {
    /* xorshift64* — Marsaglia (2003), Vigna (2014)
     * Period: 2^64 - 1, fast, good statistical properties */
    *seed ^= *seed >> 12;
    *seed ^= *seed << 25;
    *seed ^= *seed >> 27;
    /* Multiply by 0x2545F4914F6CDD1D and take upper 53 bits for double */
    uint64_t x = *seed * 0x2545F4914F6CDD1DULL;
    /* Return [0, 1) using top 53 bits (double precision mantissa) */
    return (double)(x >> 11) * 0x1.0p-53;
}

/* ─── L3: White Gaussian Noise Generation (Box-Muller) ────────────────────── */

void mems_generate_white_noise(double *output, size_t n,
                                double std_dev, uint64_t *seed) {
    if (!output || n == 0 || std_dev < 0.0) return;

    for (size_t i = 0; i < n; i += 2) {
        /* Box-Muller: two independent U(0,1) → two N(0,1) */
        double u1 = mems_rand64(seed);
        double u2 = mems_rand64(seed);
        /* Avoid log(0) */
        if (u1 < 1e-15) u1 = 1e-15;
        double r = sqrt(-2.0 * log(u1));
        double theta = 2.0 * M_PI * u2;
        double z1 = r * cos(theta);
        output[i] = z1 * std_dev;
        if (i + 1 < n) {
            double z2 = r * sin(theta);
            output[i + 1] = z2 * std_dev;
        }
    }
}

/* ─── L3: 1/f Noise Generation via Fractional Integration ───────────────────
 *
 * 1/f noise (flicker noise) is ubiquitous in MEMS — bias instability
 * appears as 1/f frequency modulation.
 *
 * Kasdin (1995) method: sum filtered white noise with 1/sqrt(f) response.
 * Approximated by direct filtering of white noise:
 *   y[n] = y[n-1] + α·(w[n] - y[n-1]/tau)
 *
 * More accurate: use the fractional differencing filter:
 *   y[n] = w[n] + Σ_{k=1}^{n} h[k]·y[n-k]
 *   where h[k] = -Γ(k-α/2) / (Γ(-α/2)·Γ(k+1)) ≈ (α/2)·k^{-(1+α/2)}
 *   for large k.
 *
 * Here we use a simplified but physically motivated approximation.
 */

void mems_generate_flicker_noise(double *output, size_t n,
                                  double std_dev, double alpha,
                                  uint64_t *seed) {
    if (!output || n == 0) return;
    if (alpha < 0.1) alpha = 0.1;  /* clamp to near-white */
    if (alpha > 2.0) alpha = 2.0;  /* clamp to near-random-walk */

    /* Generate driving white noise */
    double *white = (double *)malloc(n * sizeof(double));
    if (!white) return;
    mems_generate_white_noise(white, n, std_dev, seed);

    /* Apply leaky integrator per pole to approximate 1/f^alpha */
    /* For alpha = 1: 1/f, use a pole at low frequency */
    size_t num_poles = (size_t)(alpha * 5.0 + 0.5);  /* more poles = closer to ideal */
    if (num_poles < 1) num_poles = 1;
    if (num_poles > 10) num_poles = 10;

    for (size_t p = 0; p < num_poles; p++) {
        double tau = (double)(1 << (p + 2));  /* geometrically spaced time constants */
        double a = 1.0 / (tau + 1.0);  /* filter pole */
        double state = 0.0;
        for (size_t i = 0; i < n; i++) {
            state = (1.0 - a) * white[i] + a * state;
            white[i] = state;
        }
    }

    /* Scale to preserve approximate variance */
    double scale = 1.0 / sqrt((double)num_poles);
    for (size_t i = 0; i < n; i++) {
        output[i] = white[i] * scale;
    }

    free(white);
}

/* ─── L3: Allan Variance Computation ────────────────────────────────────────
 *
 * Allan variance σ²(τ) measures time-domain stability.
 * For a dataset x[0..N-1] with sample interval τ₀:
 *
 * Overlapping Allan variance (fully overlapped estimator):
 *   σ²(τ) = (1/(2·(N-2m))) · Σ_{i=0}^{N-2m-1} (x[i+2m] - 2·x[i+m] + x[i])²
 *
 * where m = τ/τ₀ is the cluster size factor.
 *
 * Each noise process has a distinct slope on the Allan deviation plot:
 *   White PM:        σ(τ) ∝ τ⁻¹   (slope -1)
 *   Flicker PM:      σ(τ) ∝ τ⁻¹   (slope -1)  — often grouped
 *   White FM (ARW):  σ(τ) ∝ τ⁻½   (slope -1/2)
 *   Flicker FM (BI): σ(τ) ∝ τ⁰    (slope 0, minimum)
 *   Random Walk FM:  σ(τ) ∝ τ⁺½   (slope +1/2)
 *   Rate Ramp:       σ(τ) ∝ τ⁺¹   (slope +1)
 *
 * The bias instability is the minimum of the Allan deviation curve.
 */

double mems_allan_variance(const double *data, size_t n,
                            double sample_rate, double tau) {
    if (!data || n < 4 || sample_rate <= 0.0 || tau <= 0.0) return -1.0;

    /* Cluster size: m = tau / tau_0 */
    double tau_0 = 1.0 / sample_rate;
    size_t m = (size_t)(tau / tau_0 + 0.5);
    if (m < 1) m = 1;
    if (m > n / 3) return -1.0;  /* need at least 3 independent clusters */

    /* Overlapping Allan variance */
    double sum_sq = 0.0;
    size_t count = 0;
    for (size_t i = 0; i + 2 * m < n; i++) {
        double diff = data[i + 2 * m] - 2.0 * data[i + m] + data[i];
        sum_sq += diff * diff;
        count++;
    }

    if (count < 2) return -1.0;
    double avar = sum_sq / (2.0 * (double)count * tau * tau);
    return avar;
}

MemsAllanResult mems_allan_analysis(const double *data, size_t n,
                                     double sample_rate, size_t num_points) {
    MemsAllanResult result;
    memset(&result, 0, sizeof(result));

    if (!data || n < 8 || sample_rate <= 0.0 || num_points < 3) {
        return result;
    }

    result.points = (MemsAllanPoint *)malloc(num_points * sizeof(MemsAllanPoint));
    if (!result.points) return result;

    result.sample_rate = sample_rate;
    result.total_time = (double)n / sample_rate;
    result.num_points = num_points;

    double tau_0 = 1.0 / sample_rate;
    double tau_min = tau_0;
    double tau_max = result.total_time / 3.0;  /* at least 3 clusters */

    for (size_t i = 0; i < num_points; i++) {
        /* Log-spaced τ values */
        double log_tau_min = log10(tau_min);
        double log_tau_max = log10(tau_max);
        double log_tau = log_tau_min +
                         (double)i / (double)(num_points - 1) *
                         (log_tau_max - log_tau_min);
        double tau = pow(10.0, log_tau);

        result.points[i].tau = tau;
        result.points[i].cluster_size = (size_t)(tau / tau_0 + 0.5);
        if (result.points[i].cluster_size < 1)
            result.points[i].cluster_size = 1;

        double avar = mems_allan_variance(data, n, sample_rate, tau);
        result.points[i].avar = avar;
        result.points[i].adev = (avar >= 0.0) ? sqrt(avar) : -1.0;

        /* Confidence intervals: ±1σ ≈ 1/sqrt(2·N_clusters)
         * For overlapping estimator, effective degrees of freedom ≈ N-2m+2 */
        size_t m = result.points[i].cluster_size;
        result.points[i].num_clusters = (n > 2 * m) ? (n - 2 * m + 2) : 0;
        if (result.points[i].num_clusters > 0 && avar >= 0.0) {
            double err = 1.0 / sqrt(2.0 * (double)result.points[i].num_clusters);
            result.points[i].conf_hi = (1.0 + err) * avar;
            result.points[i].conf_lo = (1.0 - err) * avar;
            if (result.points[i].conf_lo < 0.0)
                result.points[i].conf_lo = 0.0;
        }
    }

    /* Fit noise coefficients */
    mems_allan_fit_noise(&result);

    return result;
}

void mems_allan_fit_noise(MemsAllanResult *result) {
    if (!result || !result->points || result->num_points < 3) return;

    /* Identify noise coefficients from Allan deviation curve.
     * We search for characteristic slopes in log-log space:
     *
     * Slope -1/2 region → ARW: N = σ(τ=1)  (extrapolated)
     * Slope 0 region → Bias Instability: B = min(σ(τ))
     * Slope +1/2 region → Rate Random Walk: K = σ(τ=3)
     * Slope -1 region → Quantization Noise: Q = σ(τ=τ₀)·τ₀/√3
     * Slope +1 region → Rate Ramp: R = σ(τ)·τ for large τ
     */

    double min_adev = INFINITY;
    double sum_log_tau = 0.0, sum_log_adev = 0.0;
    double sum_log_tau2 = 0.0, sum_log_tau_log_adev = 0.0;
    size_t n_valid = 0;

    /* Find minimum (bias instability) */
    for (size_t i = 0; i < result->num_points; i++) {
        if (result->points[i].adev > 0.0) {
            if (result->points[i].adev < min_adev) {
                min_adev = result->points[i].adev;
            }
            double lt = log10(result->points[i].tau);
            double la = log10(result->points[i].adev);
            sum_log_tau += lt;
            sum_log_adev += la;
            sum_log_tau2 += lt * lt;
            sum_log_tau_log_adev += lt * la;
            n_valid++;
        }
    }

    result->bias_instability = (min_adev < INFINITY) ? min_adev : 0.0;

    /* ARW: Extrapolate from short τ — fit slope in first few valid points */
    if (n_valid >= 2) {
        /* Fit slope to estimate ARW: σ(τ) ≈ N / sqrt(τ)
         * → log σ = log N - 0.5·log τ
         * → N = σ(τ_short) · sqrt(τ_short)
         */
        for (size_t i = 0; i < result->num_points; i++) {
            if (result->points[i].adev > 0.0) {
                result->angle_rw = result->points[i].adev *
                                   sqrt(result->points[i].tau);
                break;
            }
        }

        /* Rate RW: from longer τ — σ(τ) ≈ K · sqrt(τ/3)
         * → K = σ(τ_long) · sqrt(3/τ_long) */
        double last_adev = 0.0, last_tau = 0.0;
        for (size_t i = 0; i < result->num_points; i++) {
            if (result->points[i].adev > 0.0) {
                last_adev = result->points[i].adev;
                last_tau = result->points[i].tau;
            }
        }
        if (last_tau > 0.0) {
            result->rate_rw = last_adev * sqrt(3.0 / last_tau);
        }
    }

    /* Quantization noise from shortest τ */
    if (result->points[0].adev > 0.0) {
        double tau_0 = 1.0 / result->sample_rate;
        result->quant_noise = result->points[0].adev * tau_0 / M_SQRT2;
    }
}

/* ─── L3: Welch's PSD Estimation ────────────────────────────────────────────
 *
 * Welch's method (1967) reduces variance of periodogram by averaging
 * windowed segments.
 *
 * Procedure:
 * 1. Divide N samples into K overlapping segments of length L
 * 2. Apply window (Hann) to each segment
 * 3. Compute FFT of each windowed segment
 * 4. Average periodograms: P̂(f) = (1/K)·Σ |X_k(f)|²
 *
 * The Hann window: w[n] = 0.5·(1 - cos(2π·n/(L-1)))
 * Its noise bandwidth: ENBW = 1.5 / L (bins normalized)
 */

/* ─── Internal FFT (radix-2 DIT Cooley-Tukey) ─────────────────────────── */

static void mems_fft_complex(double *real, double *imag, size_t n) {
    /* Bit-reverse */
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double tr = real[i], ti = imag[i];
            real[i] = real[j]; imag[i] = imag[j];
            real[j] = tr; imag[j] = ti;
        }
    }

    /* FFT butterfly */
    for (size_t len = 2; len <= n; len <<= 1) {
        double ang = -2.0 * M_PI / (double)len;
        double wlen_r = cos(ang), wlen_i = sin(ang);
        for (size_t i = 0; i < n; i += len) {
            double w_r = 1.0, w_i = 0.0;
            for (size_t j = 0; j < len / 2; j++) {
                double u_r = real[i + j], u_i = imag[i + j];
                size_t idx = i + j + len / 2;
                double t_r = w_r * real[idx] - w_i * imag[idx];
                double t_i = w_r * imag[idx] + w_i * real[idx];
                real[idx] = u_r - t_r;
                imag[idx] = u_i - t_i;
                real[i + j] = u_r + t_r;
                imag[i + j] = u_i + t_i;
                double nw_r = w_r * wlen_r - w_i * wlen_i;
                double nw_i = w_r * wlen_i + w_i * wlen_r;
                w_r = nw_r; w_i = nw_i;
            }
        }
    }
}

static size_t mems_next_pow2(size_t n) {
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

MemsPsdEstimate mems_welch_psd(const double *data, size_t n,
                                double sample_rate, size_t seg_len,
                                double overlap) {
    MemsPsdEstimate psd;
    memset(&psd, 0, sizeof(psd));

    if (!data || n < 16 || sample_rate <= 0.0) return psd;

    /* Auto segment length: ~n/8, rounded to power of 2 */
    if (seg_len == 0) {
        seg_len = mems_next_pow2(n / 8);
    } else {
        seg_len = mems_next_pow2(seg_len);
    }
    if (seg_len < 4 || seg_len > n) return psd;

    /* Default overlap: 50% */
    if (overlap <= 0.0) overlap = 0.5;
    if (overlap >= 1.0) overlap = 0.9;

    size_t step = (size_t)((double)seg_len * (1.0 - overlap));
    if (step < 1) step = 1;

    size_t num_segments = (n - seg_len) / step + 1;
    if (num_segments < 2) num_segments = 2;

    /* Only need positive frequencies: seg_len/2 + 1 bins */
    psd.num_bins = seg_len / 2 + 1;
    psd.freq = (double *)calloc(psd.num_bins, sizeof(double));
    psd.psd  = (double *)calloc(psd.num_bins, sizeof(double));
    if (!psd.freq || !psd.psd) {
        free(psd.freq); free(psd.psd);
        psd.freq = NULL; psd.psd = NULL;
        return psd;
    }

    double *segment_real = (double *)calloc(seg_len, sizeof(double));
    double *segment_imag = (double *)calloc(seg_len, sizeof(double));
    if (!segment_real || !segment_imag) {
        free(segment_real); free(segment_imag);
        return psd;
    }

    /* Hann window and window power */
    double *window = (double *)malloc(seg_len * sizeof(double));
    double window_power = 0.0;
    if (window) {
        for (size_t i = 0; i < seg_len; i++) {
            window[i] = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i /
                                         (double)(seg_len - 1)));
            window_power += window[i] * window[i];
        }
    }

    /* Process segments */
    for (size_t seg = 0; seg < num_segments; seg++) {
        size_t offset = seg * step;
        for (size_t i = 0; i < seg_len && (offset + i) < n; i++) {
            double val = data[offset + i];
            if (window) val *= window[i];
            segment_real[i] = val;
            segment_imag[i] = 0.0;
        }

        mems_fft_complex(segment_real, segment_imag, seg_len);

        /* Accumulate periodogram */
        for (size_t i = 0; i < psd.num_bins; i++) {
            psd.psd[i] += (segment_real[i] * segment_real[i] +
                           segment_imag[i] * segment_imag[i]);
        }
    }

    /* Normalize */
    double scale;
    if (window && window_power > 0.0) {
        scale = 1.0 / ((double)num_segments * window_power * sample_rate);
    } else {
        scale = 1.0 / ((double)num_segments * (double)seg_len * sample_rate);
    }

    /* Double the one-sided PSD (except DC and Nyquist) */
    for (size_t i = 0; i < psd.num_bins; i++) {
        psd.psd[i] *= scale;
        if (i > 0 && i < psd.num_bins - 1) psd.psd[i] *= 2.0;
        psd.freq[i] = (double)i * sample_rate / (double)seg_len;
    }

    psd.resolution = sample_rate / (double)seg_len;

    /* Compute total power */
    for (size_t i = 0; i < psd.num_bins; i++) {
        psd.total_power += psd.psd[i] * psd.resolution;
    }

    free(window);
    free(segment_real);
    free(segment_imag);
    return psd;
}

double mems_psd_rms(const MemsPsdEstimate *psd, double f_low, double f_hi) {
    if (!psd || !psd->freq || !psd->psd || f_low < 0.0 || f_hi <= f_low)
        return 0.0;

    double power = 0.0;
    for (size_t i = 0; i < psd->num_bins; i++) {
        if (psd->freq[i] >= f_low && psd->freq[i] <= f_hi) {
            double psd_val = psd->psd[i];
            if (psd_val < 0.0) psd_val = 0.0;
            power += psd_val;
        }
    }
    power *= psd->resolution;
    if (power < 0.0) power = 0.0;
    return sqrt(power);
}

void mems_psd_free(MemsPsdEstimate *psd) {
    if (psd) {
        free(psd->freq);
        free(psd->psd);
        psd->freq = NULL;
        psd->psd = NULL;
        psd->num_bins = 0;
    }
}

/* ─── L3: Noise Propagation Through 2nd-order System ────────────────────────
 *
 * For a white noise input with PSD S_0 passed through H(s) = ω_n²/(s²+2ζω_n·s+ω_n²),
 * the output variance is:
 *
 *   σ²_out = S_0 · ∫₀^∞ |H(f)|² df = S_0 · ω_n / (4·ζ)
 *
 * This is the noise equivalent bandwidth: NEB = ω_n / (4·ζ) [Hz, rad]
 * In Hz: NEB = f_n · π / (2·ζ)
 *
 * Reference: Brown & Hwang, "Introduction to Random Signals" (2012) §4.4
 */

double mems_noise_propagate_2nd_order(double input_psd_density,
                                       double natural_freq_Hz,
                                       double zeta) {
    if (natural_freq_Hz <= 0.0 || zeta <= 0.0) return 0.0;
    double omega_n = 2.0 * M_PI * natural_freq_Hz;
    /* Noise equivalent bandwidth (rad/s): ω_n / (4·ζ)
     * Output variance = input_PSD · NEB = input_PSD · ω_n / (4·ζ) */
    double output_variance = input_psd_density * omega_n / (4.0 * zeta);
    if (output_variance < 0.0) return 0.0;
    return sqrt(output_variance);
}
