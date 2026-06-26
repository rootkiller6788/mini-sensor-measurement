/**
 * @file adc_metrics.c
 * @brief ADC dynamic performance metrics: FFT-based SINAD, SNR, THD,
 *        SFDR, noise floor, window functions, two-tone IMD.
 *
 * Implements IEEE Std 1241-2010 dynamic test procedures.
 * Each function embodies an independent measurement concept.
 */

#include "adc_metrics.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =========================================================================
 * L5 — In-Place Radix-2 FFT (Cooley-Tukey, Decimation-in-Time)
 * ========================================================================= */

/**
 * Bit-reversal permutation for FFT reordering.
 */
static uint32_t bit_reverse(uint32_t x, uint32_t log2_n)
{
    uint32_t result = 0;
    for (uint32_t i = 0; i < log2_n; i++) {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    return result;
}

/**
 * Compute forward FFT of real-valued time-domain data.
 *
 * Uses complex-out radix-2 Cooley-Tukey FFT.
 * Output: X[0..N/2] = complex spectrum (single-sided magnitude).
 *
 * @param time_data  Real input samples, length N (power of 2).
 * @param N          Number of samples (must be power of 2).
 * @param[out] real_out Real part of spectrum [N/2+1].
 * @param[out] imag_out Imag part of spectrum [N/2+1].
 *
 * Complexity: O(N log N)
 */
static void fft_real_to_complex(const double *time_data, uint32_t N,
                                 double *real_out, double *imag_out)
{
    if (!time_data || !real_out || !imag_out || N < 2) return;

    /* Check N is power of 2 */
    uint32_t log2_n = 0;
    uint32_t temp = N;
    while (temp > 1) { temp >>= 1; log2_n++; }
    if ((1U << log2_n) != N) return;

    /* Copy to complex array with bit-reversed order */
    double complex *data = (double complex *)malloc(N * sizeof(double complex));
    if (!data) return;

    for (uint32_t i = 0; i < N; i++) {
        uint32_t br = bit_reverse(i, log2_n);
        data[br] = time_data[i] + 0.0 * I;
    }

    /* Cooley-Tukey FFT */
    for (uint32_t stage = 1; stage <= log2_n; stage++) {
        uint32_t step = 1U << stage;
        uint32_t half_step = step >> 1;
        double complex w_m = cexp(-2.0 * M_PI * I / (double)step);

        for (uint32_t k = 0; k < N; k += step) {
            double complex w = 1.0 + 0.0 * I;
            for (uint32_t j = 0; j < half_step; j++) {
                double complex t = w * data[k + j + half_step];
                double complex u = data[k + j];
                data[k + j] = u + t;
                data[k + j + half_step] = u - t;
                w *= w_m;
            }
        }
    }

    /* Extract single-sided spectrum: 0 to N/2 */
    for (uint32_t i = 0; i <= N / 2; i++) {
        real_out[i] = creal(data[i]);
        imag_out[i] = cimag(data[i]);
    }

    free(data);
}

/**
 * Compute the magnitude spectrum normalized to full-scale in dB.
 *
 * For ADC testing, full-scale (0 dBFS) corresponds to a sine wave
 * with peak amplitude = V_ref/2 (code = 2^{N-1}).
 *
 * @param real_fft   Real part of FFT.
 * @param imag_fft   Imag part of FFT.
 * @param N          FFT length.
 * @param[out] mag_dBFS Magnitude spectrum in dBFS [N/2+1].
 */
static void fft_magnitude_dBFS(const double *real_fft, const double *imag_fft,
                                uint32_t N, double *mag_dBFS)
{
    uint32_t half_N = N / 2 + 1;
    /* Normalization: divide by N (FFT scaling), then by N/2 for full-scale sine */
    double norm_factor = N; /* peak sine bin amplitude = A * N / 2 */

    for (uint32_t i = 0; i < half_N; i++) {
        double mag = sqrt(real_fft[i] * real_fft[i] +
                          imag_fft[i] * imag_fft[i]);
        /* DC and Nyquist bins use N normalization; others use N/2 for peak */
        if (i == 0 || i == half_N - 1) {
            mag /= norm_factor;
        } else {
            mag /= (norm_factor / 2.0);
        }
        /* Convert to dBFS: 0 dBFS = full-scale peak amplitude */
        if (mag > 1e-15) {
            mag_dBFS[i] = 20.0 * log10(mag);
        } else {
            mag_dBFS[i] = -300.0; /* Noise floor */
        }
    }
}

/* =========================================================================
 * L5 — Window Functions
 * ========================================================================= */

void adc_apply_window(double *samples, uint32_t N, uint32_t window_type)
{
    if (!samples || N < 2) return;

    for (uint32_t n = 0; n < N; n++) {
        double w;
        switch (window_type) {
        case 0: /* Rectangular */
            w = 1.0;
            break;
        case 1: /* Hann */
            w = 0.5 * (1.0 - cos(2.0 * M_PI * (double)n / (double)(N - 1)));
            break;
        case 2: /* Hamming */
            w = 0.54 - 0.46 * cos(2.0 * M_PI * (double)n / (double)(N - 1));
            break;
        case 3: /* Blackman */
            w = 0.42 - 0.5 * cos(2.0 * M_PI * (double)n / (double)(N - 1)) +
                0.08 * cos(4.0 * M_PI * (double)n / (double)(N - 1));
            break;
        case 4: /* Blackman-Harris (3-term) */
            w = 0.35875 - 0.48829 * cos(2.0 * M_PI * (double)n / (double)(N - 1)) +
                0.14128 * cos(4.0 * M_PI * (double)n / (double)(N - 1)) -
                0.01168 * cos(6.0 * M_PI * (double)n / (double)(N - 1));
            break;
        default:
            w = 1.0;
            break;
        }
        samples[n] *= w;
    }
}

double adc_window_coherent_gain(uint32_t window_type)
{
    /* Coherent gain = average of window values.
     * Used to correct the amplitude after windowing:
     * corrected_amplitude = measured_amplitude / coherent_gain. */
    switch (window_type) {
    case 0: return 1.0;                    /* Rectangular */
    case 1: return 0.5;                    /* Hann */
    case 2: return 0.54;                   /* Hamming */
    case 3: return 0.42;                   /* Blackman */
    case 4: return 0.35875;                /* Blackman-Harris */
    default: return 1.0;
    }
}

double adc_window_enbw(uint32_t window_type)
{
    /* Equivalent Noise BandWidth in FFT bins.
     * ENBW = N * Σ(w[n]²) / (Σ(w[n]))²
     *
     * This is the bandwidth of a rectangular filter that would
     * pass the same amount of noise as the window. */
    switch (window_type) {
    case 0: return 1.0;     /* Rectangular */
    case 1: return 1.5;     /* Hann */
    case 2: return 1.3628;  /* Hamming */
    case 3: return 1.7268;  /* Blackman */
    case 4: return 2.0044;  /* Blackman-Harris */
    default: return 1.0;
    }
}

/* =========================================================================
 * L5 — Fundamental Frequency Estimation
 * ========================================================================= */

void adc_estimate_tone_params(const double *samples, uint32_t num_samples,
                               double fs_hz,
                               double *est_freq, double *est_ampl,
                               double *est_phase)
{
    if (!samples || !est_freq || !est_ampl || !est_phase || num_samples < 4)
        return;

    /* Find the next power of 2 for FFT */
    uint32_t N_fft = 1;
    while (N_fft < num_samples) N_fft <<= 1;

    /* Allocate padded array */
    double *padded = (double *)calloc(N_fft, sizeof(double));
    if (!padded) return;
    memcpy(padded, samples, num_samples * sizeof(double));

    /* Apply Hann window to reduce spectral leakage */
    adc_apply_window(padded, N_fft, 1);

    /* FFT */
    uint32_t half_N = N_fft / 2 + 1;
    double *real_fft = (double *)malloc(half_N * sizeof(double));
    double *imag_fft = (double *)malloc(half_N * sizeof(double));
    if (!real_fft || !imag_fft) {
        free(padded); free(real_fft); free(imag_fft);
        return;
    }

    fft_real_to_complex(padded, N_fft, real_fft, imag_fft);

    /* Find peak bin (skip DC) */
    uint32_t peak_bin = 1;
    double peak_mag = 0.0;
    for (uint32_t i = 1; i < half_N; i++) {
        double mag = sqrt(real_fft[i] * real_fft[i] +
                          imag_fft[i] * imag_fft[i]);
        if (mag > peak_mag) {
            peak_mag = mag;
            peak_bin = i;
        }
    }

    /* Interpolated FFT for improved frequency estimation.
     * Use parabolic interpolation on magnitude:
     *   Δbin = 0.5 * (M_left - M_right) / (M_left - 2*M_peak + M_right) */

    double mag_left = 0.0, mag_right = 0.0;
    if (peak_bin > 0) {
        mag_left = sqrt(real_fft[peak_bin - 1] * real_fft[peak_bin - 1] +
                        imag_fft[peak_bin - 1] * imag_fft[peak_bin - 1]);
    }
    if (peak_bin + 1 < half_N) {
        mag_right = sqrt(real_fft[peak_bin + 1] * real_fft[peak_bin + 1] +
                         imag_fft[peak_bin + 1] * imag_fft[peak_bin + 1]);
    }

    double delta_bin = 0.0;
    double denom = mag_left - 2.0 * peak_mag + mag_right;
    if (fabs(denom) > 1e-15) {
        delta_bin = 0.5 * (mag_left - mag_right) / denom;
    }

    /* Estimated frequency */
    double bin_freq = fs_hz / (double)N_fft;
    *est_freq = ((double)peak_bin + delta_bin) * bin_freq;

    /* Amplitude: normalize by window coherent gain and FFT scaling */
    double coherent_gain = adc_window_coherent_gain(1); /* Hann */
    *est_ampl = (2.0 * peak_mag) / ((double)N_fft * coherent_gain);

    /* Phase at peak bin */
    *est_phase = atan2(imag_fft[peak_bin], real_fft[peak_bin]);

    free(padded);
    free(real_fft);
    free(imag_fft);
}

/* =========================================================================
 * L1/L5 — SINAD Computation
 * ========================================================================= */

double adc_sinad_from_samples(const double *samples, uint32_t num_samples,
                               double fs_hz, double expected_f_hz)
{
    if (!samples || num_samples < 8 || fs_hz <= 0.0) return NAN;

    /* FFT setup */
    uint32_t N_fft = 1;
    while (N_fft < num_samples) N_fft <<= 1;

    double *padded = (double *)calloc(N_fft, sizeof(double));
    if (!padded) return NAN;
    memcpy(padded, samples, num_samples * sizeof(double));

    /* Use minimum 4-term Blackman-Harris window for best dynamic range */
    adc_apply_window(padded, N_fft, 4);

    uint32_t half_N = N_fft / 2 + 1;
    double *real_fft = (double *)malloc(half_N * sizeof(double));
    double *imag_fft = (double *)malloc(half_N * sizeof(double));
    double *mag_dBFS = (double *)malloc(half_N * sizeof(double));

    if (!real_fft || !imag_fft || !mag_dBFS) {
        free(padded); free(real_fft); free(imag_fft); free(mag_dBFS);
        return NAN;
    }

    fft_real_to_complex(padded, N_fft, real_fft, imag_fft);
    fft_magnitude_dBFS(real_fft, imag_fft, N_fft, mag_dBFS);

    /* Find fundamental bin */
    double bin_res = fs_hz / (double)N_fft;
    uint32_t fund_bin = (uint32_t)llround(expected_f_hz / bin_res);
    if (fund_bin >= half_N) fund_bin = half_N - 1;

    /* Total power (dBFS to linear) */
    double total_power = 0.0;
    double signal_power = 0.0;
    for (uint32_t i = 1; i < half_N; i++) { /* Skip DC */
        double linear_power = pow(10.0, mag_dBFS[i] / 10.0);
        total_power += linear_power;

        /* Fundamental and its immediate neighbors carry signal power */
        uint32_t low = (fund_bin > 0) ? fund_bin - 1 : fund_bin;
        uint32_t high = fund_bin + 1;
        if (high >= half_N) high = half_N - 1;
        if (i >= low && i <= high) {
            signal_power += linear_power;
        }
    }

    /* Noise + Distortion power = total - signal */
    double nd_power = total_power - signal_power;
    if (nd_power <= 0.0) nd_power = 1e-20;

    double sinad_db = 10.0 * log10(signal_power / nd_power);

    free(padded);
    free(real_fft);
    free(imag_fft);
    free(mag_dBFS);

    return sinad_db;
}

/* =========================================================================
 * L5 — SNR Excluding Harmonics
 * ========================================================================= */

double adc_snr_excl_harmonics(const double *samples, uint32_t num_samples,
                               double fs_hz, double expected_f_hz,
                               uint32_t n_harmonics)
{
    if (!samples || num_samples < 8 || fs_hz <= 0.0) return NAN;

    uint32_t N_fft = 1;
    while (N_fft < num_samples) N_fft <<= 1;

    double *padded = (double *)calloc(N_fft, sizeof(double));
    if (!padded) return NAN;
    memcpy(padded, samples, num_samples * sizeof(double));

    adc_apply_window(padded, N_fft, 4);

    uint32_t half_N = N_fft / 2 + 1;
    double *real_fft = (double *)malloc(half_N * sizeof(double));
    double *imag_fft = (double *)malloc(half_N * sizeof(double));
    double *mag_dBFS = (double *)malloc(half_N * sizeof(double));
    int *exclude_mask = (int *)calloc(half_N, sizeof(int));

    if (!real_fft || !imag_fft || !mag_dBFS || !exclude_mask) {
        free(padded); free(real_fft); free(imag_fft);
        free(mag_dBFS); free(exclude_mask);
        return NAN;
    }

    fft_real_to_complex(padded, N_fft, real_fft, imag_fft);
    fft_magnitude_dBFS(real_fft, imag_fft, N_fft, mag_dBFS);

    double bin_res = fs_hz / (double)N_fft;
    uint32_t fund_bin = (uint32_t)llround(expected_f_hz / bin_res);

    /* Mark fundamental bins for exclusion */
    if (fund_bin > 0) exclude_mask[fund_bin - 1] = 1;
    exclude_mask[fund_bin] = 1;
    if (fund_bin + 1 < half_N) exclude_mask[fund_bin + 1] = 1;

    /* Mark harmonic bins for exclusion */
    for (uint32_t h = 2; h <= n_harmonics + 1; h++) {
        uint32_t harm_bin = (uint32_t)llround((double)h * expected_f_hz / bin_res);
        /* Allow for spectral leakage: ±2 bins around harmonic */
        for (int32_t offset = -2; offset <= 2; offset++) {
            int32_t bin = (int32_t)harm_bin + offset;
            if (bin > 0 && bin < (int32_t)half_N) {
                exclude_mask[bin] = 1;
            }
        }
    }
    /* Exclude DC */
    exclude_mask[0] = 1;

    /* Compute signal and noise power */
    double signal_power = 0.0;
    double noise_power = 0.0;

    for (uint32_t i = 1; i < half_N; i++) {
        double linear_power = pow(10.0, mag_dBFS[i] / 10.0);
        if (i >= fund_bin - 1 && i <= fund_bin + 1) {
            signal_power += linear_power;
        } else if (!exclude_mask[i]) {
            noise_power += linear_power;
        }
    }

    double snr_db = (noise_power > 0.0 && signal_power > 0.0) ?
                    10.0 * log10(signal_power / noise_power) : NAN;

    free(padded); free(real_fft); free(imag_fft);
    free(mag_dBFS); free(exclude_mask);

    return snr_db;
}

/* =========================================================================
 * L5 — THD Computation
 * ========================================================================= */

double adc_thd_from_harmonics(const double *harmonics_db,
                               uint32_t num_harmonics)
{
    if (!harmonics_db || num_harmonics < 2) return NAN;

    /* THD = 10 * log10( sum_{k=2..H} 10^{P_k/10} / 10^{P_1/10} )
     *     = 10 * log10( sum_{k=2..H} 10^{(P_k - P_1)/10} )
     *
     * where P_k is the power of the k-th harmonic in dB. */

    double fundamental_linear = pow(10.0, harmonics_db[0] / 10.0);
    double harmonic_sum_linear = 0.0;

    for (uint32_t k = 1; k < num_harmonics; k++) {
        harmonic_sum_linear += pow(10.0, harmonics_db[k] / 10.0);
    }

    if (fundamental_linear <= 0.0) return NAN;
    double thd_ratio = harmonic_sum_linear / fundamental_linear;

    return 10.0 * log10(thd_ratio);
}

double adc_thd_db_to_percent(double thd_db)
{
    return 100.0 * pow(10.0, thd_db / 20.0);
}

double adc_thd_percent_to_db(double thd_percent)
{
    if (thd_percent <= 0.0) return -INFINITY;
    return 20.0 * log10(thd_percent / 100.0);
}

/* =========================================================================
 * L5 — SFDR Computation
 * ========================================================================= */

void adc_sfdr_from_spectrum(const double *spectrum_dBFS,
                             uint32_t spectrum_len,
                             uint32_t fund_bin,
                             const uint32_t *harmonic_bins,
                             uint32_t n_harmonics,
                             double *sfdr_dBc,
                             double *sfdr_dBFS,
                             uint32_t *spur_freq_bin)
{
    if (!spectrum_dBFS || spectrum_len < 4 || fund_bin >= spectrum_len) {
        if (sfdr_dBc) *sfdr_dBc = NAN;
        if (sfdr_dBFS) *sfdr_dBFS = NAN;
        if (spur_freq_bin) *spur_freq_bin = 0;
        return;
    }

    /* Create exclusion mask: fundamental + harmonics + DC + their neighbors */
    int *exclude = (int *)calloc(spectrum_len, sizeof(int));
    if (!exclude) return;

    /* Exclude DC */
    exclude[0] = 1;
    /* Exclude fundamental and ±2 bins */
    for (int32_t off = -2; off <= 2; off++) {
        int32_t idx = (int32_t)fund_bin + off;
        if (idx > 0 && idx < (int32_t)spectrum_len) exclude[idx] = 1;
    }
    /* Exclude harmonics and ±1 bin */
    for (uint32_t h = 0; h < n_harmonics; h++) {
        uint32_t hb = harmonic_bins[h];
        for (int32_t off = -1; off <= 1; off++) {
            int32_t idx = (int32_t)hb + off;
            if (idx > 0 && idx < (int32_t)spectrum_len) exclude[idx] = 1;
        }
    }

    /* Find largest spur among non-excluded bins */
    double max_spur = -INFINITY;
    uint32_t max_spur_bin = 0;
    for (uint32_t i = 1; i < spectrum_len; i++) {
        if (!exclude[i] && spectrum_dBFS[i] > max_spur) {
            max_spur = spectrum_dBFS[i];
            max_spur_bin = i;
        }
    }

    double fundamental_power = spectrum_dBFS[fund_bin];

    if (sfdr_dBc) {
        *sfdr_dBc = (max_spur > -INFINITY) ?
                     fundamental_power - max_spur : NAN;
    }
    if (sfdr_dBFS) {
        *sfdr_dBFS = (max_spur > -INFINITY) ?
                      0.0 - max_spur : NAN;
    }
    if (spur_freq_bin) {
        *spur_freq_bin = max_spur_bin;
    }

    free(exclude);
}

/* =========================================================================
 * L5 — Noise Floor
 * ========================================================================= */

double adc_noise_floor_dBFS(const double *spectrum_dBFS,
                             uint32_t spectrum_len,
                             const uint32_t *exclude_bins,
                             uint32_t n_exclude)
{
    if (!spectrum_dBFS || spectrum_len == 0) return NAN;

    /* Create exclusion mask */
    int *exclude = (int *)calloc(spectrum_len, sizeof(int));
    if (!exclude) return NAN;

    exclude[0] = 1; /* Always exclude DC */
    for (uint32_t i = 0; i < n_exclude; i++) {
        if (exclude_bins[i] < spectrum_len) {
            exclude[exclude_bins[i]] = 1;
            /* Also exclude adjacent bins to be conservative */
            if (exclude_bins[i] > 0) exclude[exclude_bins[i] - 1] = 1;
            if (exclude_bins[i] + 1 < spectrum_len) exclude[exclude_bins[i] + 1] = 1;
        }
    }

    /* Average noise power over non-excluded bins */
    double noise_power_sum = 0.0;
    uint32_t noise_bin_count = 0;

    for (uint32_t i = 1; i < spectrum_len; i++) {
        if (!exclude[i]) {
            noise_power_sum += pow(10.0, spectrum_dBFS[i] / 10.0);
            noise_bin_count++;
        }
    }

    free(exclude);

    if (noise_bin_count == 0 || noise_power_sum <= 0.0) return NAN;

    double avg_noise_power = noise_power_sum / (double)noise_bin_count;
    return 10.0 * log10(avg_noise_power);
}

/* =========================================================================
 * L5 — Coherent Sampling Check
 * ========================================================================= */

int adc_check_coherent_sampling(double f_in_hz, double fs_hz,
                                 uint32_t num_samples,
                                 uint32_t *cycles)
{
    if (f_in_hz <= 0.0 || fs_hz <= 0.0 || num_samples == 0) return 0;

    /* Coherent sampling: M * f_in = J * f_s for integers M (cycles), J (samples).
     * Equivalently: M = J * f_in / f_s must be integer.
     *
     * With J = num_samples, we need M = num_samples * f_in / f_s ∈ ℤ.
     * And GCD(J, M) = 1 for unique samples. */

    double exact_cycles = (double)num_samples * f_in_hz / fs_hz;
    uint32_t M = (uint32_t)llround(exact_cycles);
    if (cycles) *cycles = M;

    /* Check if M is sufficiently close to integer */
    if (fabs(exact_cycles - (double)M) > 0.01) return 0;

    /* Check mutually prime: GCD(num_samples, M) == 1 */
    uint32_t a = num_samples;
    uint32_t b = M;
    while (b != 0) {
        uint32_t t = b;
        b = a % b;
        a = t;
    }
    return (a == 1) ? 1 : 0;
}

/* =========================================================================
 * L5 — Two-Tone IMD
 * ========================================================================= */

void adc_compute_two_tone_imd(const double *samples, uint32_t num_samples,
                               double fs_hz, double f1_hz, double f2_hz,
                               adc_imd_metrics_t *metrics)
{
    if (!samples || !metrics || num_samples < 16 || fs_hz <= 0.0) return;

    memset(metrics, 0, sizeof(*metrics));
    metrics->f1_hz = f1_hz;
    metrics->f2_hz = f2_hz;

    /* FFT */
    uint32_t N_fft = 1;
    while (N_fft < num_samples) N_fft <<= 1;

    double *padded = (double *)calloc(N_fft, sizeof(double));
    if (!padded) return;
    memcpy(padded, samples, num_samples * sizeof(double));
    adc_apply_window(padded, N_fft, 4);

    uint32_t half_N = N_fft / 2 + 1;
    double *real_fft = (double *)malloc(half_N * sizeof(double));
    double *imag_fft = (double *)malloc(half_N * sizeof(double));
    double *mag_dBFS = (double *)malloc(half_N * sizeof(double));

    if (!real_fft || !imag_fft || !mag_dBFS) {
        free(padded); free(real_fft); free(imag_fft); free(mag_dBFS);
        return;
    }

    fft_real_to_complex(padded, N_fft, real_fft, imag_fft);
    fft_magnitude_dBFS(real_fft, imag_fft, N_fft, mag_dBFS);

    double bin_res = fs_hz / (double)N_fft;
    uint32_t bin_f1 = (uint32_t)llround(f1_hz / bin_res);
    uint32_t bin_f2 = (uint32_t)llround(f2_hz / bin_res);

    if (bin_f1 >= half_N) bin_f1 = half_N - 1;
    if (bin_f2 >= half_N) bin_f2 = half_N - 1;

    /* Fundamental power: average of two tones */
    double fund_power_linear = (pow(10.0, mag_dBFS[bin_f1] / 10.0) +
                                 pow(10.0, mag_dBFS[bin_f2] / 10.0)) * 0.5;
    double fund_power_dB = 10.0 * log10(fund_power_linear);

    /* IM2 products: |f1 ± f2| */
    double im2_freqs[2] = { fabs(f1_hz - f2_hz), f1_hz + f2_hz };
    double im2_max_dBc = -INFINITY;

    for (int i = 0; i < 2; i++) {
        if (im2_freqs[i] < fs_hz / 2.0 && im2_freqs[i] > 0.0) {
            uint32_t bin = (uint32_t)llround(im2_freqs[i] / bin_res);
            if (bin < half_N) {
                double level = mag_dBFS[bin];
                double dBc = fund_power_dB - level;
                if (dBc > im2_max_dBc) im2_max_dBc = dBc;
            }
        }
    }
    metrics->imd2_db = (im2_max_dBc > -INFINITY) ? im2_max_dBc : NAN;

    /* IM3 products: |2f1 ± f2|, |f1 ± 2f2| */
    double im3_freqs[4] = {
        fabs(2.0 * f1_hz - f2_hz),
        2.0 * f1_hz + f2_hz,
        fabs(2.0 * f2_hz - f1_hz),
        2.0 * f2_hz + f1_hz
    };
    double im3_max_dBc = -INFINITY;

    for (int i = 0; i < 4; i++) {
        if (im3_freqs[i] < fs_hz / 2.0 && im3_freqs[i] > 0.0) {
            uint32_t bin = (uint32_t)llround(im3_freqs[i] / bin_res);
            if (bin < half_N) {
                double level = mag_dBFS[bin];
                double dBc = fund_power_dB - level;
                if (dBc > im3_max_dBc) im3_max_dBc = dBc;
            }
        }
    }
    metrics->imd3_db = (im3_max_dBc > -INFINITY) ? im3_max_dBc : NAN;

    /* Store all 6 IM product frequencies and levels */
    double all_freqs[6];
    all_freqs[0] = im2_freqs[0]; all_freqs[1] = im2_freqs[1];
    all_freqs[2] = im3_freqs[0]; all_freqs[3] = im3_freqs[1];
    all_freqs[4] = im3_freqs[2]; all_freqs[5] = im3_freqs[3];

    for (int i = 0; i < 6; i++) {
        metrics->imd_freqs_hz[i] = all_freqs[i];
        if (all_freqs[i] < fs_hz / 2.0 && all_freqs[i] > 0.0) {
            uint32_t bin = (uint32_t)llround(all_freqs[i] / bin_res);
            if (bin < half_N) {
                metrics->imd_levels_db[i] = mag_dBFS[bin];
            } else {
                metrics->imd_levels_db[i] = NAN;
            }
        } else {
            metrics->imd_levels_db[i] = NAN;
        }
    }

    free(padded); free(real_fft); free(imag_fft); free(mag_dBFS);
}

/* =========================================================================
 * L6 — Full Dynamic Test Suite (IEEE 1241)
 * ========================================================================= */

void adc_full_dynamic_test(const double *samples, uint32_t num_samples,
                            double fs_hz, double f_in_hz,
                            uint32_t n_harmonics,
                            adc_dynamic_metrics_t *metrics)
{
    if (!samples || !metrics || num_samples < 16 || fs_hz <= 0.0) return;

    memset(metrics, 0, sizeof(*metrics));

    /* FFT setup */
    uint32_t N_fft = 1;
    while (N_fft < num_samples) N_fft <<= 1;

    double *padded = (double *)calloc(N_fft, sizeof(double));
    if (!padded) return;
    memcpy(padded, samples, num_samples * sizeof(double));
    adc_apply_window(padded, N_fft, 4); /* Blackman-Harris for best isolation */

    uint32_t half_N = N_fft / 2 + 1;
    double *real_fft = (double *)malloc(half_N * sizeof(double));
    double *imag_fft = (double *)malloc(half_N * sizeof(double));
    double *mag_dBFS = (double *)malloc(half_N * sizeof(double));

    if (!real_fft || !imag_fft || !mag_dBFS) {
        free(padded); free(real_fft); free(imag_fft); free(mag_dBFS);
        return;
    }

    fft_real_to_complex(padded, N_fft, real_fft, imag_fft);
    fft_magnitude_dBFS(real_fft, imag_fft, N_fft, mag_dBFS);

    double bin_res = fs_hz / (double)N_fft;
    uint32_t fund_bin = (uint32_t)llround(f_in_hz / bin_res);
    if (fund_bin >= half_N) fund_bin = 0;

    /* --- Fundamental power --- */
    double fund_power_linear = 0.0;
    for (int32_t off = -1; off <= 1; off++) {
        int32_t idx = (int32_t)fund_bin + off;
        if (idx > 0 && idx < (int32_t)half_N) {
            fund_power_linear += pow(10.0, mag_dBFS[idx] / 10.0);
        }
    }
    metrics->fundamental_power_dBFS = (fund_power_linear > 0.0) ?
        10.0 * log10(fund_power_linear) : -300.0;

    /* --- DC offset --- */
    metrics->dc_offset_dBFS = mag_dBFS[0];

    /* --- Total power --- */
    double total_power_linear = 0.0;
    for (uint32_t i = 1; i < half_N; i++) {
        total_power_linear += pow(10.0, mag_dBFS[i] / 10.0);
    }

    /* --- Harmonic powers --- */
    metrics->num_harmonics = n_harmonics;
    metrics->harmonic_powers = (double *)malloc(n_harmonics * sizeof(double));
    uint32_t *harmonic_bins = (uint32_t *)malloc(n_harmonics * sizeof(uint32_t));

    double harmonic_power_sum = 0.0;
    int *exclude_mask = (int *)calloc(half_N, sizeof(int));

    /* Mark fundamental range */
    for (int32_t off = -2; off <= 2; off++) {
        int32_t idx = (int32_t)fund_bin + off;
        if (idx > 0 && idx < (int32_t)half_N) exclude_mask[idx] = 1;
    }
    exclude_mask[0] = 1; /* Exclude DC */

    if (metrics->harmonic_powers && harmonic_bins) {
        for (uint32_t h = 0; h < n_harmonics; h++) {
            uint32_t harm_order = h + 2; /* 2nd, 3rd, ... harmonic */
            double harm_freq = (double)harm_order * f_in_hz;
            double folded_freq = harm_freq;

            /* Handle aliasing (fold back to [0, fs/2]) */
            double k = round(harm_freq / fs_hz);
            folded_freq = fabs(harm_freq - k * fs_hz);
            if (folded_freq > fs_hz / 2.0) folded_freq = fs_hz - folded_freq;

            uint32_t harm_bin = (uint32_t)llround(folded_freq / bin_res);
            if (harm_bin >= half_N) harm_bin = half_N - 1;
            harmonic_bins[h] = harm_bin;

            /* Sum power in ±1 bin around harmonic */
            double harm_power = 0.0;
            for (int32_t off = -1; off <= 1; off++) {
                int32_t idx = (int32_t)harm_bin + off;
                if (idx > 0 && idx < (int32_t)half_N) {
                    double p = pow(10.0, mag_dBFS[idx] / 10.0);
                    harm_power += p;
                    exclude_mask[idx] = 1;
                }
            }
            metrics->harmonic_powers[h] = (harm_power > 0.0) ?
                10.0 * log10(harm_power) : -300.0;
            if (h > 0) harmonic_power_sum += harm_power;
        }
    }

    /* --- SNR (noise only, excl harmonics) --- */
    double noise_power_linear = 0.0;
    for (uint32_t i = 1; i < half_N; i++) {
        if (!exclude_mask[i]) {
            noise_power_linear += pow(10.0, mag_dBFS[i] / 10.0);
        }
    }
    metrics->snr_db = (noise_power_linear > 0.0 && fund_power_linear > 0.0) ?
        10.0 * log10(fund_power_linear / noise_power_linear) : NAN;

    /* --- SINAD = Signal / (Noise + Distortion) --- */
    double nd_power = noise_power_linear + harmonic_power_sum;
    metrics->sinad_db = (nd_power > 0.0 && fund_power_linear > 0.0) ?
        10.0 * log10(fund_power_linear / nd_power) : NAN;

    /* --- THD --- */
    double thd_linear = (fund_power_linear > 0.0) ?
        harmonic_power_sum / fund_power_linear : 0.0;
    metrics->thd_db = (thd_linear > 0.0) ?
        10.0 * log10(thd_linear) : -300.0;
    metrics->thd_percent = 100.0 * sqrt(thd_linear);

    /* --- SFDR --- */
    uint32_t spur_bin;
    adc_sfdr_from_spectrum(mag_dBFS, half_N, fund_bin,
                            harmonic_bins, n_harmonics,
                            &metrics->sfdr_db, &metrics->sfdr_dbFS, &spur_bin);

    /* --- Noise floor --- */
    /* Exclude fund + harmonics */
    uint32_t exc_count = n_harmonics + 3;
    uint32_t *exc = (uint32_t *)malloc(exc_count * sizeof(uint32_t));
    if (exc) {
        exc[0] = 0; /* DC */
        exc[1] = fund_bin;
        exc[2] = (fund_bin > 0) ? fund_bin - 1 : 0;
        for (uint32_t i = 0; i < n_harmonics; i++) {
            exc[3 + i] = harmonic_bins[i];
        }
        metrics->noise_floor_dBFS = adc_noise_floor_dBFS(mag_dBFS, half_N,
                                                          exc, exc_count);
        free(exc);
    }

    /* --- ENOB --- */
    metrics->enob = adc_enob_from_sinad(metrics->sinad_db);

    free(padded); free(real_fft); free(imag_fft); free(mag_dBFS);
    free(harmonic_bins); free(exclude_mask);
}
