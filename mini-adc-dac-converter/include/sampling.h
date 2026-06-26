/**
 * @file sampling.h
 * @brief Nyquist-Shannon sampling theory, aliasing analysis,
 *        oversampling, decimation, and interpolation.
 *
 * Knowledge Coverage:
 *   L1 - Definitions: sampling rate (f_s), Nyquist frequency (f_s/2),
 *        aliasing, undersampling (band-pass sampling), oversampling ratio
 *   L2 - Core Concepts: impulse-train model of ideal sampling,
 *        periodic replication of spectrum, reconstruction via sinc interpolation
 *   L3 - Mathematical Structures: Fourier transform of sampled signal,
 *        Poisson summation formula for spectral copies, sampling of
 *        band-pass signals (intermediate frequency sampling)
 *   L4 - Fundamental Laws: Nyquist-Shannon sampling theorem,
 *        Whittaker-Shannon interpolation formula,
 *        Oversampling SNR improvement: 3 dB/octave (10*log10(OSR))
 *   L5 - Algorithms: decimation filter (CIC + FIR), interpolation by
 *        zero-insertion + filtering, polyphase decomposition for
 *        efficient rate conversion
 *
 * References:
 *   - Shannon, "Communication in the Presence of Noise" Proc. IRE (1949)
 *   - Nyquist, "Certain Topics in Telegraph Transmission Theory" (1928)
 *   - Oppenheim & Schafer (2010), Ch. 4, 10
 *   - Hogenauer, "An Economical Class of Digital Filters for
 *     Decimation and Interpolation" IEEE Trans ASSP (1981)
 *
 * Course Alignment:
 *   MIT 6.003 Signals and Systems (Ch. 7, sampling)
 *   Stanford EE264 Digital Signal Processing (Ch. 4, multirate)
 *   Berkeley EE123 DSP (Ch. 4, sampling and reconstruction)
 */

#ifndef SAMPLING_H
#define SAMPLING_H

#include <stdint.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1 — Sampling Parameters
 * ========================================================================= */

/** Anti-aliasing filter specification. */
typedef struct {
    double          passband_hz;        /**< Passband edge frequency */
    double          stopband_hz;        /**< Stopband edge frequency */
    double          passband_ripple_db; /**< Passband ripple [dB] */
    double          stopband_atten_db;  /**< Minimum stopband attenuation [dB] */
    uint32_t        filter_order;       /**< Filter order */
} antialias_filter_spec_t;

/**
 * @brief Sampling system specification.
 */
typedef struct {
    double          fs_hz;              /**< Sampling frequency */
    double          f_nyquist_hz;       /**< Nyquist frequency = fs/2 */
    double          signal_bandwidth_hz;/**< Bandwidth of input signal */
    double          oversampling_ratio; /**< OSR = fs / (2 * BW) */
    int             is_bandpass;        /**< 1 if band-pass sampling */
    double          bandpass_center_hz; /**< Center freq for band-pass sampling */
    antialias_filter_spec_t aa_filter; /**< Required anti-alias filter */
} sampling_system_t;

/* =========================================================================
 * L1 — Aliasing / Spectral Copy
 * ========================================================================= */

/**
 * @brief Compute the frequency at which a tone appears after aliasing.
 *
 * Given f_in > f_nyquist, compute f_alias = |f_in - k * f_s|
 * where k = round(f_in / f_s), with result folded to [0, f_s/2].
 *
 * @param f_input_hz  Input frequency.
 * @param fs_hz       Sampling frequency.
 * @return Aliased frequency in [0, fs/2].
 *
 * Complexity: O(1)
 */
double aliased_frequency(double f_input_hz, double fs_hz);

/**
 * @brief Compute all spectral image frequencies for a given tone.
 *
 * The sampled spectrum contains copies at f_image = |±f_in + k * f_s| for all k.
 * This computes the first N_images that fall within [0, n * fs/2].
 *
 * @param f_input_hz  Input tone frequency.
 * @param fs_hz       Sampling rate.
 * @param n_images    Number of image frequencies to compute.
 * @param[out] images Array of length n_images (filled with NAN for non-existent).
 */
void sampling_spectral_images(double f_input_hz, double fs_hz,
                               uint32_t n_images, double *images);

/* =========================================================================
 * L4 — Nyquist-Shannon Sampling Theorem
 * ========================================================================= */

/**
 * @brief Verify if a signal bandwidth satisfies the Nyquist criterion.
 *
 * Band-limited signal with bandwidth BW is perfectly reconstructable
 * from its samples if f_s > 2 * BW.
 *
 * For band-pass signals (BW centered at f_c), the condition is:
 *   2*f_u / floor(f_u / BW)  ≤  f_s  ≤  2*f_l / floor(f_l / BW)
 * where f_u = f_c + BW/2, f_l = f_c - BW/2.
 *
 * @param bandwidth_hz  Signal bandwidth (two-sided).
 * @param fs_hz         Sampling frequency.
 * @param center_hz     Center frequency (0 for baseband).
 * @return 1 if Nyquist satisfied, 0 if undersampling.
 *
 * Complexity: O(1)
 */
int nyquist_check(double bandwidth_hz, double fs_hz, double center_hz);

/**
 * @brief Compute the minimum sampling rate to avoid aliasing.
 *
 * For baseband: f_s_min = 2 * BW.
 * For bandpass: f_s_min = 2 * BW * (1 + floor(f_c / BW)).
 *
 * @param bandwidth_hz Signal bandwidth.
 * @param center_hz    Center frequency (0 for baseband).
 * @return Minimum sampling frequency.
 */
double minimum_sample_rate(double bandwidth_hz, double center_hz);

/* =========================================================================
 * L4 — Oversampling SNR Improvement
 * ========================================================================= */

/**
 * @brief Compute SNR improvement from oversampling.
 *
 * Oversampling by factor OSR spreads quantization noise over a wider
 * bandwidth. After digital low-pass filtering to the original bandwidth:
 *   SNR_improvement = 10 * log10(OSR)  [dB]
 *
 * For every factor-of-4 increase in OSR, SNR improves by ~6 dB,
 * equivalent to 1 additional bit of resolution.
 *
 * @param osr Oversampling ratio (f_s / (2 * BW_signal)).
 * @return SNR improvement in dB.
 *
 * Complexity: O(1)
 * Reference: Proakis & Manolakis "DSP" (2007) Ch. 9
 */
double oversampling_snr_gain(double osr);

/**
 * @brief Compute effective resolution increase from oversampling.
 *
 * Δbits = log2(OSR) / 2
 *
 * @param osr Oversampling ratio.
 * @return Additional effective bits.
 */
double oversampling_effective_bits(double osr);

/* =========================================================================
 * L5 — Sinc Interpolation (Ideal Reconstruction)
 * ========================================================================= */

/**
 * @brief Ideal band-limited interpolation using the sinc kernel.
 *
 * x(t) = sum_{n} x[n] * sinc((t - n*T) / T)
 * where T = 1/fs.
 *
 * In practice, a windowed sinc is used for finite-length interpolation.
 *
 * @param samples    Array of sample values.
 * @param N          Number of sample points.
 * @param fs_hz      Sampling rate.
 * @param t_target   Target time instant for interpolation.
 * @param half_length Half-length of the sinc window (number of samples each side).
 * @param window_type Window type for truncation (0=rect, 1=hann, 2=hamming).
 * @return Interpolated value at t_target.
 *
 * Complexity: O(half_length)
 * Reference: Whittaker (1915), Shannon (1949)
 */
double sinc_interpolate(const double *samples, uint32_t N,
                         double fs_hz, double t_target,
                         uint32_t half_length, uint32_t window_type);

/**
 * @brief Interpolate to an arbitrary grid using windowed sinc.
 *
 * @param samples      Input samples at rate fs_hz.
 * @param n_input      Number of input samples.
 * @param fs_hz        Sampling rate.
 * @param t_out        Output time grid.
 * @param n_output     Number of output points.
 * @param[out] interpolated Output array.
 * @param half_length  Sinc half-length.
 * @param window_type  Window type.
 *
 * Complexity: O(n_output * half_length)
 */
void sinc_interpolate_grid(const double *samples, uint32_t n_input,
                            double fs_hz,
                            const double *t_out, uint32_t n_output,
                            double *interpolated,
                            uint32_t half_length, uint32_t window_type);

/* =========================================================================
 * L5 — Decimation (Downsampling)
 * ========================================================================= */

/**
 * @brief Decimate a signal by factor M (keep every M-th sample).
 *
 * Requires prior anti-aliasing filtering! This function implements
 * a combined CIC (Cascaded Integrator-Comb) filter + decimation.
 *
 * CIC filter H(z) = [(1 - z^{-RM}) / (1 - z^{-1})]^N
 * where R is rate change, M is differential delay (typically 1 or 2),
 * N is number of stages.
 *
 * @param input        Input samples at rate fs_in.
 * @param n_input      Number of input samples.
 * @param decim_factor Decimation factor M.
 * @param cic_stages   Number of CIC stages N.
 * @param cic_delay    Differential delay (1 or 2).
 * @param[out] output  Decimated output (length ≈ n_input / M).
 * @param[out] n_output Number of output samples.
 *
 * Complexity: O(n_input * cic_stages)
 * Reference: Hogenauer (1981)
 */
void cic_decimate(const double *input, uint32_t n_input,
                   uint32_t decim_factor, uint32_t cic_stages,
                   uint32_t cic_delay,
                   double *output, uint32_t *n_output);

/**
 * @brief Simple FIR anti-alias filter + decimate by M.
 *
 * @param input         Input samples.
 * @param n_input       Number of input samples.
 * @param decim_factor  Decimation factor.
 * @param fir_coeffs    FIR anti-aliasing filter coefficients.
 * @param fir_len       Number of FIR taps.
 * @param[out] output   Decimated output.
 * @param[out] n_output Number of output samples.
 *
 * Complexity: O(n_input * fir_len / decim_factor)
 */
void fir_decimate(const double *input, uint32_t n_input,
                   uint32_t decim_factor,
                   const double *fir_coeffs, uint32_t fir_len,
                   double *output, uint32_t *n_output);

/* =========================================================================
 * L5 — Interpolation (Upsampling)
 * ========================================================================= */

/**
 * @brief Interpolate by factor L using zero-insertion + FIR filter.
 *
 * x_up[m] = x[m/L] if m % L == 0 else 0
 * Followed by low-pass FIR filter with gain L.
 *
 * @param input         Input samples.
 * @param n_input       Number of input samples.
 * @param interp_factor Interpolation factor L.
 * @param fir_coeffs    Interpolation FIR filter.
 * @param fir_len       Number of taps.
 * @param[out] output   Interpolated output (length = n_input * L).
 * @param[out] n_output Number of output samples.
 *
 * Complexity: O(n_input * L * fir_len)
 */
void fir_interpolate(const double *input, uint32_t n_input,
                      uint32_t interp_factor,
                      const double *fir_coeffs, uint32_t fir_len,
                      double *output, uint32_t *n_output);

/* =========================================================================
 * L5 — Rational Rate Conversion (L/M)
 * ========================================================================= */

/**
 * @brief Sample rate conversion by rational factor L/M.
 *
 * Combines interpolation by L with decimation by M using a single
 * common FIR filter (polyphase decomposition).
 *
 * @param input      Input samples at rate fs_in.
 * @param n_input    Number of input samples.
 * @param L          Interpolation factor (numerator).
 * @param M          Decimation factor (denominator).
 * @param fir_coeffs Combined interpolation/decimation filter.
 * @param fir_len    Filter length (must be multiple of L).
 * @param[out] output Converted samples (length ≈ n_input * L / M).
 * @param[out] n_output Number of output samples.
 *
 * Complexity: O(n_input * L * fir_len / (L*M)) = O(n_input * fir_len / M)
 */
void rational_rate_convert(const double *input, uint32_t n_input,
                            uint32_t L, uint32_t M,
                            const double *fir_coeffs, uint32_t fir_len,
                            double *output, uint32_t *n_output);

/* =========================================================================
 * L5 — Band-Pass Sampling (IF Sampling)
 * ========================================================================= */

/**
 * @brief Determine valid sampling frequency range for band-pass signal.
 *
 * For a band-pass signal occupying [f_l, f_u] with BW = f_u - f_l,
 * valid fs must satisfy (for integer n):
 *   2*f_u / n  ≤  f_s  ≤  2*f_l / (n-1)
 * where n_max = floor(f_u / BW).
 *
 * @param f_l_hz        Lower band edge.
 * @param f_u_hz        Upper band edge.
 * @param[out] fs_min   Minimum valid sampling rate.
 * @param[out] fs_max   Maximum valid sampling rate (at n=1).
 * @param[out] optimal_n Optimal value of n.
 *
 * Complexity: O(n_max) where n_max = floor(f_u/BW)
 */
void bandpass_sampling_range(double f_l_hz, double f_u_hz,
                              double *fs_min, double *fs_max,
                              uint32_t *optimal_n);

/* =========================================================================
 * L6 — Sampling System Simulation
 * ========================================================================= */

/**
 * @brief Simulate ideal sampling of a continuous-time signal.
 *
 * Apply the sampling system model: anti-alias filter → sample → quantize.
 *
 * @param ct_signal     Function pointer to continuous-time signal.
 * @param ct_signal_ctx User context for signal function.
 * @param fs_hz         Sampling rate.
 * @param n_samples     Number of samples.
 * @param t_start       Start time.
 * @param[out] samples  Output samples [n_samples].
 *
 * Complexity: O(n_samples) × cost of ct_signal()
 */
void simulate_ideal_sampling(double (*ct_signal)(double t, void *ctx),
                              void *ct_signal_ctx,
                              double fs_hz, uint32_t n_samples,
                              double t_start,
                              double *samples);

/* =========================================================================
 * L6 — Full Acquisition Chain Simulation
 * ========================================================================= */

/**
 * @brief Simulate complete ADC acquisition including AA filter,
 *        sampling jitter, and quantization.
 *
 * @param ct_signal     Continuous-time signal.
 * @param ctx           Signal context.
 * @param fs_hz         Sampling rate.
 * @param n_samples     Number of samples.
 * @param t_start       Start time.
 * @param aa_filter_b   FIR anti-alias filter numerator coefficients.
 * @param aa_filter_a   FIR anti-alias filter denominator (1.0 for FIR).
 * @param aa_order      Filter order.
 * @param jitter_rms_ps RMS aperture jitter.
 * @param n_bits        ADC resolution.
 * @param v_ref         ADC reference voltage.
 * @param[out] codes    Output ADC codes [n_samples].
 *
 * Complexity: O(n_samples * (aa_order + 1))
 */
void simulate_adc_acquisition(double (*ct_signal)(double t, void *ctx),
                               void *ctx,
                               double fs_hz, uint32_t n_samples,
                               double t_start,
                               const double *aa_filter_b,
                               const double *aa_filter_a,
                               uint32_t aa_order,
                               double jitter_rms_ps,
                               uint32_t n_bits, double v_ref,
                               uint64_t *codes);

#ifdef __cplusplus
}
#endif

#endif /* SAMPLING_H */
