/**
 * @file example_sampling.c
 * @brief End-to-end sampling and reconstruction example.
 *
 * Demonstrates:
 * - Aliasing visualization (sweep input frequency past Nyquist)
 * - Oversampling SNR improvement computation
 * - Band-pass sampling range calculation
 * - Sinc interpolation of a sampled signal
 */
#include "sampling.h"
#include "adc_dac_core.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("=== Sampling Theory Demonstration ===\n\n");

    /* ── Part 1: Aliasing Sweep ── */
    printf("--- Aliasing Sweep (fs = 100 Hz) ---\n");
    printf("f_in [Hz]  f_alias [Hz]  Aliased?\n");
    printf("-----------------------------------\n");

    double fs = 100.0;
    double f_inputs[] = {10, 40, 60, 90, 110, 150, 190, 210, 490};
    for (int i = 0; i < 9; i++) {
        double fi = f_inputs[i];
        double fa = aliased_frequency(fi, fs);
        printf("  %6.1f     %6.1f        %s\n",
               fi, fa, (fabs(fi - fa) > 0.1) ? "YES" : "no");
    }

    /* ── Part 2: Oversampling SNR Improvement ── */
    printf("\n--- Oversampling SNR Gain ---\n");
    printf("OSR    ΔSNR[dB]  ΔENOB[bits]\n");
    printf("--------------------------------\n");

    double osr_vals[] = {2, 4, 8, 16, 32, 64, 128, 256};
    for (int i = 0; i < 8; i++) {
        double gain = oversampling_snr_gain(osr_vals[i]);
        double bits = oversampling_effective_bits(osr_vals[i]);
        printf(" %-6.0f %-9.2f %-12.3f\n", osr_vals[i], gain, bits);
    }
    printf("→ Each 4× OSR increase gives ~1 extra bit of resolution.\n");

    /* ── Part 3: Band-Pass Sampling ── */
    printf("\n--- Band-Pass Sampling (70-90 MHz IF band) ---\n");

    double fs_min, fs_max;
    uint32_t n_opt;
    bandpass_sampling_range(70e6, 90e6, &fs_min, &fs_max, &n_opt);
    printf("Signal band: 70-90 MHz (BW = 20 MHz)\n");
    printf("Valid fs range: %.2f - %.2f MHz (optimal n = %u)\n",
           fs_min/1e6, fs_max/1e6, n_opt);

    /* Check a few candidate sampling rates */
    double candidates[] = {40e6, 45e6, 50e6, 56e6, 60e6};
    for (int i = 0; i < 5; i++) {
        int ok = nyquist_check(20e6, candidates[i], 80e6);
        printf("  fs = %.0f MHz: %s\n",
               candidates[i]/1e6, ok ? "VALID" : "INVALID");
    }

    /* ── Part 4: Sinc Interpolation ── */
    printf("\n--- Sinc Interpolation (4-point sine, 10× upsampling) ---\n");

    double samples[] = {0.0, 1.0, 0.0, -1.0};
    double t_out[21];
    double interp[21];

    for (int i = 0; i < 21; i++) {
        t_out[i] = (double)i * 0.15; /* 0 to 3.0 */
    }
    sinc_interpolate_grid(samples, 4, 4.0, t_out, 21, interp, 8, 1);

    printf("t[s]    x(t)\n");
    printf("-------------\n");
    for (int i = 0; i < 21; i++) {
        printf("%5.2f  %+7.4f\n", t_out[i], interp[i]);
    }

    return 0;
}
