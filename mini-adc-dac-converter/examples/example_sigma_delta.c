/**
 * @file example_sigma_delta.c
 * @brief End-to-end ΔΣ ADC design example: Audio-band 16-bit equivalent.
 *
 * Demonstrates:
 * - Theoretical performance prediction for 1st, 2nd, 3rd order ΔΣ
 * - Single-loop modulator simulation
 * - Comparison of ideal vs simulated performance
 */
#include "sigma_delta.h"
#include "adc_dac_core.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void)
{
    printf("=== Sigma-Delta ADC Design for Audio (16-bit target) ===\n\n");

    double f_bw = 20e3;        /* Audio bandwidth */
    double v_ref = 1.0;        /* Reference voltage */

    /* Design space exploration */
    printf("Order  OSR    f_s[MHz]  DR[dB]    ENOB    Stable?\n");
    printf("------------------------------------------------------\n");

    uint32_t orders[] = {1, 2, 3};
    double osrs[] = {32, 64, 128, 256};

    for (int oi = 0; oi < 3; oi++) {
        uint32_t L = orders[oi];
        for (int si = 0; si < 4; si++) {
            double osr = osrs[si];
            double fs = 2.0 * f_bw * osr;
            double dr = sdm_dynamic_range_db(L, osr, 1, v_ref);
            double enob = sdm_equivalent_enob(L, osr, 1, v_ref);

            sdm_spec_t spec = {
                .order = L,
                .quantizer_bits = 1,
                .osr = osr,
                .fs_hz = fs,
                .bandwidth_hz = f_bw,
                .v_ref = v_ref
            };
            double limit = sdm_stable_input_limit(&spec);

            printf("%-6u %-6.0f %-10.3f %-8.1f %-7.2f %s\n",
                   L, osr, fs/1e6, dr, enob,
                   (limit > 0.7) ? "YES" : "marginal");
        }
    }

    /* Simulate a 2nd-order ΔΣ with OSR=128 (fs=5.12 MHz) */
    printf("\n--- 2nd-Order ΔΣ Simulation (OSR=128) ---\n");

    sdm_spec_t spec = {
        .order = 2,
        .quantizer_bits = 1,
        .osr = 128,
        .fs_hz = 5.12e6,
        .bandwidth_hz = f_bw,
        .v_ref = 1.0
    };

    double amplitude = 0.5;
    double duration = 0.001; /* 1 ms */

    /* Allocate buffers for modulator output */
    uint64_t n_full = (uint64_t)(spec.fs_hz * duration);
    int64_t *bitstream = (int64_t *)malloc(n_full * sizeof(int64_t));
    double *decimated = (double *)malloc((uint64_t)(n_full / 128 + 1) * sizeof(double));

    if (bitstream && decimated) {
        uint32_t n_dec;
        uint64_t n_osr;
        sdm_simulate_adc(&spec, 1e3, amplitude, duration,
                          decimated, &n_dec, bitstream, &n_osr);

        /* Count ones in bitstream (density should track input) */
        uint64_t ones = 0;
        for (uint64_t i = 0; i < n_osr && i < 10000; i++) {
            if (bitstream[i] > 0) ones++;
        }
        double density = (double)ones / 10000.0;

        printf("Input amplitude: %.2f\n", amplitude);
        printf("Bitstream density (first 10k samples): %.4f (expect ~0.5)\n",
               density);
        printf("Bitstream average: ~%.3f (should track input)\n", 2.0*density - 1.0);
        printf("Decimated samples: %u\n", n_dec);
        if (n_dec > 0) {
            printf("First decimated sample: %.6f\n", decimated[0]);
        }
    }

    /* Compute theoretical DR and ENOB */
    double dr = sdm_dynamic_range_db(2, 128, 1, v_ref);
    double enob = sdm_equivalent_enob(2, 128, 1, v_ref);
    printf("\nPredicted DR: %.1f dB → ENOB: %.2f bits\n", dr, enob);

    free(bitstream);
    free(decimated);
    return 0;
}
