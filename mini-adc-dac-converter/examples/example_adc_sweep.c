/**
 * @file example_adc_sweep.c
 * @brief End-to-end example: ADC SNR vs resolution characterization.
 *
 * Computes and displays ideal SNR, ENOB, and jitter limits
 * for ADC resolutions from 6 to 24 bits — demonstrating the
 * classic 6.02N+1.76 relationship and practical limitations.
 */
#include "adc_dac_core.h"
#include "adc_metrics.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("=== ADC Resolution Sweep ===\n");
    printf("%4s %10s %10s %10s %12s\n",
           "Bits", "SNR[dB]", "ENOB", "LSB[uV@3.3V]", "JitterLimit@10MHz");
    printf("----------------------------------------------------------\n");

    double f_in = 10e6;    /* 10 MHz input */
    double t_jitter = 1.0; /* 1 ps RMS jitter */

    for (uint32_t N = 6; N <= 24; N += 2) {
        double snr_ideal = adc_ideal_snr_db(N);
        double enob = adc_enob_from_sinad(snr_ideal);
        double v_lsb_uv = 3.3 / pow(2.0, (double)N) * 1e6;
        double snr_jitter = adc_jitter_snr_limit(f_in, t_jitter);

        printf("%4u %10.2f %10.2f %10.2f %12.2f\n",
               N, snr_ideal, enob, v_lsb_uv, snr_jitter);

        if (snr_jitter < snr_ideal) {
            printf("  → Jitter-limited! f_in=%g MHz creates noise floor\n",
                   f_in / 1e6);
        }
    }

    printf("\nKey insight: Beyond ~14 bits, 1 ps jitter at 10 MHz\n");
    printf("limits SNR to ~84 dB regardless of resolution.\n");

    return 0;
}
