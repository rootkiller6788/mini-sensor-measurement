/**
 * @file example_sar_sim.c
 * @brief End-to-end SAR ADC simulation with non-idealities.
 *
 * Demonstrates a complete 10-bit SAR ADC conversion including:
 * - kT/C noise on CDAC
 * - Comparator offset and noise
 * - Capacitor mismatch effects
 * - DNL/INL computation from simulated data
 */
#include "sar_adc.h"
#include "adc_dac_core.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void)
{
    printf("=== SAR ADC 10-bit Simulation ===\n\n");

    /* Configure a 10-bit SAR ADC with practical parameters */
    sar_adc_spec_t spec = {
        .resolution_bits = 10,
        .sample_rate_hz = 1e6,
        .v_ref_volts = 3.3,
        .c_unit_fF = 50.0,
        .switching = SAR_SWITCHING_MONOTONIC,
        .comparator_offset_mV = 2.0,
        .comparator_noise_uVrms = 100.0
    };

    printf("Spec: %u-bit, fs=%.1f MS/s, Vref=%.1f V\n",
           spec.resolution_bits, spec.sample_rate_hz/1e6, spec.v_ref_volts);

    /* Compute kT/C noise limit */
    cdac_array_t cdac;
    cdac_init(spec.resolution_bits, spec.c_unit_fF, &cdac);
    double c_total = cdac.c_total_fF * 1e-15;
    double ktc_rms = ktc_noise_rms(300.0, c_total);
    double v_lsb = spec.v_ref_volts / pow(2.0, spec.resolution_bits);
    printf("C_total = %.1f pF, kT/C noise = %.1f uVrms, V_LSB = %.1f uV\n",
           cdac.c_total_fF * 1e-3, ktc_rms * 1e6, v_lsb * 1e6);
    printf("Noise/LSB ratio: %.2f (target < 0.5)\n\n", ktc_rms / v_lsb);

    /* Simulate a ramp input sweep */
    #define N_STEPS 64
    double v_in[N_STEPS];
    uint64_t codes[N_STEPS];

    for (int i = 0; i < N_STEPS; i++) {
        v_in[i] = spec.v_ref_volts * (double)i / (double)(N_STEPS - 1);
    }

    uint64_t seed = 12345;
    sar_simulate_conversions(&spec, v_in, N_STEPS, codes, &seed);

    printf("Transfer characteristic:\n");
    printf("  Vin[V]   Code   V_recon[V]   Error[mV]\n");
    printf("  ----------------------------------------\n");
    double max_error = 0.0;
    for (int i = 0; i < N_STEPS; i += 8) {
        double v_recon = adc_code_to_voltage(codes[i], spec.v_ref_volts,
                                              spec.resolution_bits);
        double error = (v_recon - v_in[i]) * 1e3;
        if (fabs(error) > max_error) max_error = fabs(error);
        printf("  %6.3f   %4llu   %7.4f     %+7.2f\n",
               v_in[i], (unsigned long long)codes[i], v_recon, error);
    }
    printf("  Max |error|: %.2f mV\n\n", max_error);

    /* Check monotonicity */
    int monotonic = 1;
    for (int i = 1; i < N_STEPS; i++) {
        if (codes[i] < codes[i-1]) { monotonic = 0; break; }
    }
    printf("Monotonic: %s (expected YES for SAR)\n", monotonic ? "YES" : "NO");
    printf("DNL from mismatch estimate: %.3f LSB\n",
           sar_dnl_from_mismatch(spec.resolution_bits, 0.1));

    cdac_free(&cdac);
    return 0;
}
