/**
 * @file mems_noise.h
 * @brief L3-L4 MEMS Noise Analysis — Thermo-mechanical noise, Allan variance,
 *        power spectral density, noise source identification
 *
 * Covers:
 *   L3 - Power spectral density, Allan variance, noise processes (white, flicker,
 *        random walk, quantization, rate ramp, sinusoidal)
 *   L4 - Thermo-mechanical noise (Brownian noise limit), Johnson-Nyquist noise,
 *        shot noise, noise equivalent acceleration/rate
 *
 * Reference: IEEE Std 1139-2008 (Frequency Stability — Allan Variance)
 *            IEEE Std 952-2020 (Inertial Sensor Noise)
 *            Allan, "Statistics of Atomic Frequency Standards" (1966)
 *            Gabrielson, "Mechanical-Thermal Noise in MEMS" (1993)
 * Courses: MIT 6.630, Stanford EE247, Berkeley EE117,
 *          Illinois ECE 310 (DSP), ETH 227-0427 (Signal Processing)
 */

#ifndef MEMS_NOISE_H
#define MEMS_NOISE_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include "mems_basics.h"

/* ─── L3: Noise Process Types (IEEE 1139-2008) ────────────────────────────── */
typedef enum {
    MEMS_NOISE_WHITE         = 0,  /* Angle/velocity random walk */
    MEMS_NOISE_FLICKER_FM    = 1,  /* Bias instability (1/f FM) */
    MEMS_NOISE_FLICKER_PM    = 2,  /* 1/f phase noise */
    MEMS_NOISE_RANDOM_WALK_FM = 3, /* Rate/acceleration random walk */
    MEMS_NOISE_QUANTIZATION  = 4,  /* Quantization noise */
    MEMS_NOISE_RATE_RAMP     = 5,  /* Rate ramp / drift */
    MEMS_NOISE_SINUSOIDAL    = 6,  /* Periodic noise */
    MEMS_NOISE_WHITE_PM      = 7   /* White phase noise */
} MemsNoiseProcess;

const char *mems_noise_process_name(MemsNoiseProcess np);

/* ─── L3: Allan Variance Cluster ───────────────────────────────────────────── */
typedef struct {
    size_t   cluster_size;  /* number of samples per cluster */
    size_t   num_clusters;  /* number of independent clusters */
    double   tau;           /* cluster time [s] */
    double   avar;          /* Allan variance at this τ */
    double   adev;          /* Allan deviation (sqrt of avar) */
    double   conf_hi;       /* 68% upper confidence */
    double   conf_lo;       /* 68% lower confidence */
} MemsAllanPoint;

/* ─── L3: Allan Variance Analysis Results ─────────────────────────────────── */
typedef struct {
    MemsAllanPoint *points;
    size_t          num_points;
    double          sample_rate;    /* Hz */
    double          total_time;     /* s */
    /* Identified noise coefficients */
    double          angle_rw;       /* Angle/Velocity Random Walk */
    double          bias_instability; /* Bias Instability minimum */
    double          rate_rw;        /* Rate/Acceleration Random Walk */
    double          quant_noise;    /* Quantization noise */
    double          rate_ramp;      /* Rate ramp coefficient */
} MemsAllanResult;

/* ─── L3: Power Spectral Density Estimate ──────────────────────────────────── */
typedef struct {
    double *freq;        /* frequency bins [Hz] */
    double *psd;         /* power spectral density */
    size_t  num_bins;
    double  resolution;  /* frequency resolution [Hz] */
    double  total_power; /* integrated power */
} MemsPsdEstimate;

/* ─── L3: Noise Source Model ───────────────────────────────────────────────── */
typedef struct {
    MemsNoiseProcess type;
    double           coefficient;  /* noise coefficient in appropriate units */
    double           cutoff_freq;  /* Hz, for colored noise */
    double           corner_freq;  /* Hz, where 1/f meets white */
    char             description[64];
} MemsNoiseSource;

/* ─── L4: Thermo-mechanical (Brownian) Noise ──────────────────────────────────
 * Fundamental noise limit in MEMS due to molecular collisions.
 * Total noise equivalent acceleration (NEA):
 *   NEA = sqrt(4·k_B·T·ω_n / (m·Q))   [m/s²/√Hz]
 *
 * Reference: Gabrielson (1993), "Mechanical-Thermal Noise in MEMS"
 *            Levin, "Internal Thermal Noise in LIGO" (1998)
 */

/**
 * @brief Thermo-mechanical noise equivalent acceleration (NEA)
 *        NEA = sqrt(4·k_B·T·c) / m = sqrt(4·k_B·T·ω_n / (m·Q))
 *        where c = sqrt(k·m) / Q = m·ω_n / Q
 *
 * @param mass            Proof mass [kg]
 * @param natural_freq_Hz Natural frequency [Hz], ω_n / (2π)
 * @param q_factor        Quality factor
 * @param temperature_k   Temperature [K]
 * @return NEA in [m/s²/√Hz], returns 0 if parameters invalid
 *
 * Complexity: O(1)
 */
double mems_thermomech_nea(double mass, double natural_freq_Hz,
                            double q_factor, double temperature_k);

/**
 * @brief Thermo-mechanical noise equivalent rate (NER) for gyros
 *        Same physics as NEA but expressed as angular rate noise.
 *        NER = NEA / (v_drive · SF_normalized)
 *
 * @param nea         Acceleration noise [m/s²/√Hz]
 * @param v_drive     Drive velocity [m/s]
 * @param sf_norm     Normalized scale factor
 */
double mems_thermomech_ner(double nea, double v_drive, double sf_norm);

/* ─── L4: Johnson-Nyquist Noise (in MEMS readout circuits) ────────────────── */

/**
 * @brief Johnson-Nyquist voltage noise density
 *        v_n = sqrt(4·k_B·T·R)  [V/√Hz]
 *        Reference: Nyquist (1928), Johnson (1928)
 *
 * @param resistance   R [Ω]
 * @param temperature  T [K]
 * @param bandwidth    Δf [Hz], if 0 returns density
 */
double mems_johnson_noise_density(double resistance, double temperature_k);

/**
 * @brief Johnson noise RMS in a given bandwidth
 *        V_rms = sqrt(4·k_B·T·R·Δf)
 */
double mems_johnson_noise_rms(double resistance, double temperature_k,
                               double bandwidth_hz);

/* ─── L4: Shot Noise (in MEMS with small currents) ─────────────────────────── */

/**
 * @brief Shot noise current density
 *        i_n = sqrt(2·q_e·I_DC)  [A/√Hz]
 *        Reference: Schottky (1918)
 *
 * @param dc_current  I_DC [A]
 */
double mems_shot_noise_density(double dc_current);

/* ─── L3: Allan Variance Computation ───────────────────────────────────────── */

/**
 * @brief Compute Allan variance for a given cluster size τ
 *        σ²(τ) = (1/(2·τ²·(N-2n))) · Σ(x_{i+2n} - 2·x_{i+n} + x_i)²
 *        for overlapping clusters.
 *        Reference: IEEE 1139-2008, Allan (1966)
 *
 * @param data        Input time series
 * @param n           Number of data points
 * @param sample_rate Sample rate [Hz]
 * @param tau         Cluster time [s]
 * @return Allan variance value, or -1 on error
 *
 * Complexity: O(N) where N is n
 */
double mems_allan_variance(const double *data, size_t n,
                            double sample_rate, double tau);

/**
 * @brief Full Allan deviation analysis with log-spaced τ values
 *        τ from τ_min to τ_max with n_points log-spaced clusters.
 *
 * @param data        Input time series
 * @param n           Number of data points
 * @param sample_rate Sample rate [Hz]
 * @param num_points  Number of τ values to evaluate
 * @return MemsAllanResult (caller frees result.points)
 */
MemsAllanResult mems_allan_analysis(const double *data, size_t n,
                                     double sample_rate, size_t num_points);

/**
 * @brief Fit noise coefficients from Allan variance curve
 *        Using log-space linear regression:
 *        White noise: slope = -1/2
 *        Bias instability: slope = 0 (minimum)
 *        Random walk: slope = +1/2
 *        Rate ramp: slope = +1
 *        Quantization: slope = -1
 *        Reference: IEEE 952-2020 §6.3
 *
 * @param result Allan variance analysis result (modified in-place)
 */
void mems_allan_fit_noise(MemsAllanResult *result);

/* ─── L3: Power Spectral Density (Welch's Method) ──────────────────────────── */

/**
 * @brief Compute PSD using Welch's averaged periodogram
 *        Divides data into overlapping segments, applies Hann window,
 *        computes FFT, averages periodograms.
 *
 * @param data        Input time series
 * @param n           Number of points
 * @param sample_rate Sample rate [Hz]
 * @param seg_len     Segment length (0 = auto: n/8)
 * @param overlap     Overlap fraction [0, 1) (0 = 0.5 default)
 * @return MemsPsdEstimate (caller frees freq and psd arrays)
 *
 * Complexity: O(N·log(N)) with auto segment length
 */
MemsPsdEstimate mems_welch_psd(const double *data, size_t n,
                                double sample_rate, size_t seg_len,
                                double overlap);

/**
 * @brief Compute RMS noise from PSD in a frequency band
 *        RMS = sqrt(∫_{f1}^{f2} PSD(f) df)
 *
 * @param psd   PSD estimate
 * @param f_low Lower frequency [Hz]
 * @param f_hi  Upper frequency [Hz]
 */
double mems_psd_rms(const MemsPsdEstimate *psd, double f_low, double f_hi);

/**
 * @brief Free PSD estimate memory
 */
void mems_psd_free(MemsPsdEstimate *psd);

/* ─── L3: Noise Simulation ─────────────────────────────────────────────────── */

/**
 * @brief Generate white Gaussian noise samples using Box-Muller
 *        y[i] = σ · z[i] where z[i] ~ N(0,1)
 *
 * @param output     Pre-allocated output array
 * @param n          Number of samples
 * @param std_dev    Standard deviation σ
 * @param seed       PRNG seed value (modified in-place)
 */
void mems_generate_white_noise(double *output, size_t n,
                                double std_dev, uint64_t *seed);

/**
 * @brief Generate 1/f (flicker) noise via fractional integration
 *        Approximated by summing filtered Gaussian noise.
 *        Reference: Kasdin (1995), "Discrete Simulation of Colored Noise"
 *
 * @param output     Pre-allocated output array
 * @param n          Number of samples
 * @param std_dev    Standard deviation of driving white noise
 * @param alpha      1/f exponent (typically 0.5 to 1.5)
 * @param seed       PRNG seed value
 */
void mems_generate_flicker_noise(double *output, size_t n,
                                  double std_dev, double alpha,
                                  uint64_t *seed);

/**
 * @brief Simple PRNG: xorshift64*
 *        Returns uniform [0, 1) and updates seed.
 */
double mems_rand64(uint64_t *seed);

/* ─── L3: Noise Propagation Through Transfer Function ──────────────────────── */

/**
 * @brief Propagate input noise through a 2nd-order system
 *        Given input PSD S_in(f), output PSD = |H(f)|² · S_in(f)
 *        where H(f) = 1 / (1 - r² + j·2ζr) with r = f/ω_n
 *
 * @param input_psd_density   Input white noise density (constant PSD)
 * @param natural_freq_Hz     System natural frequency [Hz]
 * @param zeta                Damping ratio
 * @param output_rms          Output RMS noise [in same units as input]
 *
 * Complexity: O(1) via analytical integration
 * ∫₀^∞ |H(f)|² df = ω_n / (4·ζ) — the noise bandwidth
 */
double mems_noise_propagate_2nd_order(double input_psd_density,
                                       double natural_freq_Hz,
                                       double zeta);

#endif /* MEMS_NOISE_H */
