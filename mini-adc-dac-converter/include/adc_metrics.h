/**
 * @file adc_metrics.h
 * @brief ADC dynamic performance metrics: SNR, SINAD, SFDR, THD, ENOB, IMD.
 *
 * Knowledge Coverage:
 *   L1 - Definitions: SINAD, SFDR, THD, ENOB, SNR, two-tone IMD,
 *        noise floor, spurious-free dynamic range
 *   L2 - Core Concepts: spectral leakage, window functions for FFT-based
 *        metrics, coherent vs non-coherent sampling
 *   L4 - Fundamental Laws: Parseval's theorem for power calculation,
 *        relationship SINAD = SNR + THD (in linear, approx)
 *
 * All metrics follow IEEE Std 1241-2010 definitions.
 *
 * References:
 *   - IEEE Std 1241-2010: Terminology and Test Methods for ADC
 *   - Kester, "Data Conversion Handbook" (2005), Ch. 2,5
 *   - Oppenheim & Schafer, "Discrete-Time Signal Processing" (2010), Ch. 10
 *
 * Course Alignment:
 *   MIT 6.003 Signal Processing, Stanford EE315 Data Converters,
 *   Berkeley EE123 DSP, ETH 227-0427 Signal Processing
 */

#ifndef ADC_METRICS_H
#define ADC_METRICS_H

#include "adc_dac_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1 Definitions — Dynamic Metric Structs
 * ========================================================================= */

/**
 * @brief Result of a single-tone ADC dynamic test.
 *
 * Captures all IEEE 1241 standard dynamic parameters.
 * Typically obtained via FFT of a coherently sampled sine wave.
 */
typedef struct {
    double          sinad_db;       /**< Signal-to-Noise-and-Distortion [dB] */
    double          snr_db;         /**< Signal-to-Noise (excl. harmonics) [dB] */
    double          thd_db;         /**< Total Harmonic Distortion [dB] */
    double          thd_percent;    /**< THD as percentage */
    double          sfdr_db;        /**< Spurious-Free Dynamic Range [dB] */
    double          sfdr_dbFS;       /**< SFDR relative to full-scale */
    double          enob;           /**< Effective Number of Bits */
    double          noise_floor_dBFS; /**< RMS noise floor relative to FS */
    double          fundamental_power_dBFS; /**< Fundamental signal power */
    double          dc_offset_dBFS; /**< DC offset relative to FS */
    uint32_t        num_harmonics;  /**< Number of harmonics measured */
    double         *harmonic_powers;/**< Power in each harmonic [dBFS] */
} adc_dynamic_metrics_t;

/**
 * @brief Two-tone intermodulation distortion result.
 *
 * Tests linearity with two closely spaced tones at f1, f2.
 * Key products: IM2 at |f1±f2|, IM3 at |2f1±f2|, |f1±2f2|.
 */
typedef struct {
    double          imd2_db;        /**< 2nd-order IMD [dBc] */
    double          imd3_db;        /**< 3rd-order IMD [dBc] */
    double          imd2_dbFS;       /**< 2nd-order IMD [dBFS] */
    double          imd3_dbFS;       /**< 3rd-order IMD [dBFS] */
    double          ip2_dBm;        /**< 2nd-order intercept point */
    double          ip3_dBm;        /**< 3rd-order intercept point */
    double          f1_hz;          /**< First tone frequency */
    double          f2_hz;          /**< Second tone frequency */
    double          imd_freqs_hz[6]; /**< IM product frequencies */
    double          imd_levels_db[6];/**< IM product levels */
} adc_imd_metrics_t;

/* =========================================================================
 * L1 — SINAD (Signal-to-Noise and Distortion Ratio)
 * ========================================================================= */

/**
 * @brief Compute SINAD from time-domain samples.
 *
 * SINAD = 10 * log10( P_signal / P_noise_and_distortion )
 *
 * P_signal is computed from the fundamental frequency component
 * (fitted sine), and P_n+d = P_total - P_signal.
 *
 * Uses three-parameter least-squares sine fit (IEEE 1241 §4.1.4.1).
 *
 * @param samples       Array of time-domain ADC output samples.
 * @param num_samples   Number of samples.
 * @param fs_hz         Sampling frequency [Hz].
 * @param expected_f_hz Expected input signal frequency [Hz].
 * @return SINAD in dB, NAN on error.
 *
 * Complexity: O(N) where N = num_samples
 * Reference: IEEE Std 1241-2010 §4.1.4
 */
double adc_sinad_from_samples(const double *samples, uint32_t num_samples,
                               double fs_hz, double expected_f_hz);

/* =========================================================================
 * L5 — SNR (Excluding Harmonics)
 * ========================================================================= */

/**
 * @brief Compute SNR from time-domain samples (harmonics excluded).
 *
 * Removes the fundamental and its harmonics from total power,
 * leaving only quantization noise, thermal noise, and jitter.
 *
 * @param samples       Time-domain ADC output samples.
 * @param num_samples   Number of samples.
 * @param fs_hz         Sampling frequency.
 * @param expected_f_hz Input frequency.
 * @param n_harmonics   Number of harmonics to exclude.
 * @return SNR in dB.
 *
 * Complexity: O(N * H) where H = n_harmonics
 */
double adc_snr_excl_harmonics(const double *samples, uint32_t num_samples,
                               double fs_hz, double expected_f_hz,
                               uint32_t n_harmonics);

/* =========================================================================
 * L1/L5 — THD (Total Harmonic Distortion)
 * ========================================================================= */

/**
 * @brief Compute THD from harmonic amplitudes.
 *
 * THD = 20 * log10( sqrt(sum_{k=2..H} A_k²) / A_1 )
 * where A_1 is the fundamental amplitude, A_k are harmonic amplitudes.
 *
 * @param harmonics_db   Array of harmonic powers [dB] (index 0 = fundamental).
 * @param num_harmonics  Total number of harmonics (incl. fundamental).
 * @return THD in dB.
 *
 * Complexity: O(H)
 */
double adc_thd_from_harmonics(const double *harmonics_db,
                               uint32_t num_harmonics);

/**
 * @brief Compute THD in percent from dB value.
 */
double adc_thd_db_to_percent(double thd_db);

/**
 * @brief Compute THD in dB from percent value.
 */
double adc_thd_percent_to_db(double thd_percent);

/* =========================================================================
 * L5 — SFDR (Spurious-Free Dynamic Range)
 * ========================================================================= */

/**
 * @brief Compute SFDR from spectral data.
 *
 * SFDR = Power(fundamental) - Power(largest_spur) [dBc]
 * or
 * SFDR = 0 dBFS - Power(largest_spur) [dBFS]
 *
 * Only considers spurs not harmonically related to input.
 *
 * @param spectrum_dBFS  Array of spectrum values in dBFS.
 * @param spectrum_len   Length of single-sided spectrum.
 * @param fund_bin       FFT bin index of the fundamental.
 * @param harmonic_bins  Array of harmonic bin indices.
 * @param n_harmonics    Number of harmonics.
 * @param sfdr_dBc       Output SFDR in dBc (relative to carrier).
 * @param sfdr_dBFS       Output SFDR in dBFS (relative to full-scale).
 * @param spur_freq_bin  Output bin index of largest spur.
 *
 * Complexity: O(N)
 */
void adc_sfdr_from_spectrum(const double *spectrum_dBFS,
                             uint32_t spectrum_len,
                             uint32_t fund_bin,
                             const uint32_t *harmonic_bins,
                             uint32_t n_harmonics,
                             double *sfdr_dBc,
                             double *sfdr_dBFS,
                             uint32_t *spur_freq_bin);

/* =========================================================================
 * L5 — Spectral Averaging / Noise Floor
 * ========================================================================= */

/**
 * @brief Estimate noise floor from spectrum with signal bins excluded.
 *
 * Averages the power spectral density over "noise-only" bins
 * (excluding fundamental, harmonics, and DC).
 *
 * @param spectrum_dBFS  Single-sided spectrum in dBFS.
 * @param spectrum_len   Number of spectral bins.
 * @param exclude_bins   Array of bin indices to exclude.
 * @param n_exclude      Number of excluded bins.
 * @return RMS noise floor in dBFS.
 *
 * Complexity: O(N)
 */
double adc_noise_floor_dBFS(const double *spectrum_dBFS,
                             uint32_t spectrum_len,
                             const uint32_t *exclude_bins,
                             uint32_t n_exclude);

/* =========================================================================
 * L5 — Window Functions for Coherent/Non-Coherent Sampling
 * ========================================================================= */

/**
 * @brief Apply a window function to N samples in-place.
 *
 * Window types follow standard DSP naming (Oppenheim & Schafer Ch. 10).
 *
 * @param samples  Input/output sample array.
 * @param N        Number of samples.
 * @param window_type Window type: 0=none/rect, 1=Hann, 2=Hamming, 3=Blackman, 4=Blackman-Harris.
 *
 * Complexity: O(N)
 */
void adc_apply_window(double *samples, uint32_t N, uint32_t window_type);

/**
 * @brief Compute window processing gain (coherent gain) for normalization.
 *
 * @param window_type Window type (same enum as apply_window).
 * @return Processing gain factor.
 */
double adc_window_coherent_gain(uint32_t window_type);

/**
 * @brief Compute equivalent noise bandwidth (ENBW) of window in bins.
 */
double adc_window_enbw(uint32_t window_type);

/* =========================================================================
 * L5 — FFT-Based Tone Analysis
 * ========================================================================= */

/**
 * @brief Estimate fundamental frequency and amplitude from ADC samples.
 *
 * Uses interpolated FFT (IFFT) for improved accuracy when the
 * signal frequency is not perfectly coherent.
 *
 * @param samples       Time-domain samples.
 * @param num_samples   Number of samples.
 * @param fs_hz         Sampling frequency [Hz].
 * @param[out] est_freq Estimated fundamental frequency [Hz].
 * @param[out] est_ampl Estimated fundamental amplitude.
 * @param[out] est_phase Estimated phase [rad].
 *
 * Complexity: O(N log N) via FFT
 */
void adc_estimate_tone_params(const double *samples, uint32_t num_samples,
                               double fs_hz,
                               double *est_freq, double *est_ampl,
                               double *est_phase);

/* =========================================================================
 * L5 — Coherent Sampling Check
 * ========================================================================= */

/**
 * @brief Check if sampling parameters yield coherent sampling.
 *
 * Coherent condition: M * f_in = J * f_s, with M cycles, J samples,
 * both integers, GCD(J, M) = 1 (mutually prime).
 *
 * @param f_in_hz       Input signal frequency.
 * @param fs_hz         Sampling frequency.
 * @param num_samples   Number of samples acquired.
 * @param[out] cycles   Integer number of signal periods in record.
 * @return 1 if coherent, 0 if not.
 */
int adc_check_coherent_sampling(double f_in_hz, double fs_hz,
                                 uint32_t num_samples,
                                 uint32_t *cycles);

/* =========================================================================
 * L5 — Two-Tone IMD Computation
 * ========================================================================= */

/**
 * @brief Compute IM2 and IM3 from two-tone measurement.
 *
 * IM2 = Power_fundamental - Power_{f1±f2} [dBc]
 * IM3 = Power_fundamental - Power_{2f1±f2, 2f2±f1} [dBc]
 *
 * @param samples       Time-domain samples.
 * @param num_samples   Sample count.
 * @param fs_hz         Sampling rate.
 * @param f1_hz         First tone frequency.
 * @param f2_hz         Second tone frequency.
 * @param[out] metrics  IMD metrics output.
 *
 * Complexity: O(N log N) via FFT
 */
void adc_compute_two_tone_imd(const double *samples, uint32_t num_samples,
                               double fs_hz, double f1_hz, double f2_hz,
                               adc_imd_metrics_t *metrics);

/* =========================================================================
 * L6 — Full Dynamic Test Suite
 * ========================================================================= */

/**
 * @brief Run the complete IEEE 1241 dynamic test suite on ADC samples.
 *
 * Computes SNR, SINAD, THD, SFDR, ENOB, noise floor in one call.
 *
 * @param samples       Time-domain samples.
 * @param num_samples   Number of samples.
 * @param fs_hz         Sampling rate.
 * @param f_in_hz       Input sine frequency.
 * @param n_harmonics   Number of harmonics to consider.
 * @param[out] metrics  Full dynamic metrics.
 *
 * Complexity: O(N log N) dominated by FFT
 */
void adc_full_dynamic_test(const double *samples, uint32_t num_samples,
                            double fs_hz, double f_in_hz,
                            uint32_t n_harmonics,
                            adc_dynamic_metrics_t *metrics);

#ifdef __cplusplus
}
#endif

#endif /* ADC_METRICS_H */
