/**
 * ex1_thermocouple_chain.c - Complete Thermocouple Signal Chain
 *
 * Demonstrates a complete Type K thermocouple measurement chain
 * with CJC, amplification, anti-alias filtering, and linearization.
 *
 * L6: Thermocouple measurement problem
 * L7: Industrial temperature sensing application (Boeing, NASA)
 *
 * Signal chain:
 *   TC hot junction (T_hot unknown)
 *     -> TC-Cu cold junction at terminal block (T_cj measured)
 *        -> CJC sensor (NTC thermistor)
 *           -> Instrumentation amplifier (gain = 200)
 *              -> Anti-alias filter (Butterworth 4th-order LP, fc=10Hz)
 *                 -> 16-bit ADC (V_ref=5V, fs=60Hz)
 *                    -> Software: CJC compensation + ITS-90 linearization
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigcond_defs.h"
#include "sigcond_filter.h"
#include "sigcond_linearize.h"
#include "sigcond_cjc.h"

int main(void) {
    printf("=== Thermocouple Measurement Signal Chain ===\n\n");

    /* Configure system */
    thermocouple_nist_model_t tcK;
    nist_tc_init_typeK(&tcK);

    /* Cold junction thermistor (NTC, 10k at 25C) */
    steinhart_hart_params_t cjc_thermistor = {
        .a_coeff = 1.129241e-3, .b_coeff = 2.341077e-4,
        .c_coeff = 8.775468e-8, .use_4param = false,
        .r_ref = 10000.0, .t_ref_c = 25.0, .beta_k = 3977.0
    };

    /* Simulated measurement scenario: Boeing 747 engine exhaust monitoring
     * Target: T_hot around 800C with accuracy <= 2C */
    double t_hot_actual = 800.0;  /* degrees C */
    double t_cj_actual = 35.0;    /* ambient in engine nacelle */

    printf("Scenario: Jet engine exhaust temperature monitoring\n");
    printf("  Actual hot junction: %.1f C\n", t_hot_actual);
    printf("  Actual cold junction: %.1f C\n", t_cj_actual);

    /* Step 1: Simulate raw thermocouple voltage (T_hot to T_cj) */
    double v_hot = nist_tc_temp_to_voltage(&tcK, t_hot_actual);
    double v_cj = nist_tc_temp_to_voltage(&tcK, t_cj_actual);
    double v_measured = v_hot - v_cj;  /* Measured differential voltage */

    printf("\nStep 1 - Raw measurement:\n");
    printf("  V(T_hot ref 0C) = %.1f uV\n", v_hot);
    printf("  V(T_cj ref 0C) = %.1f uV\n", v_cj);
    printf("  V_measured = %.1f uV\n", v_measured);

    /* Step 2: Amplify */
    double gain = 200.0;
    double v_amplified = v_measured * 1e-6 * gain;  /* Convert uV to V, apply gain */
    double v_noise = 1e-6;  /* 1uV RTI noise * gain */
    double v_with_noise = v_amplified + v_noise * (rand() / (double)RAND_MAX - 0.5);

    printf("\nStep 2 - Amplification (G=%.0f):\n", gain);
    printf("  V_amplified = %.6f V\n", v_with_noise);

    /* Step 3: Anti-alias filter (Butterworth 4th order, fc=10Hz) */
    int filter_order = filter_butterworth_order(10.0, 30.0, 0.1, 80.0);
    bool aa_ok = filter_verify_antialias(10.0, 60.0, 16, 90.0);

    printf("\nStep 3 - Anti-alias filter:\n");
    printf("  Filter order: %d\n", filter_order);
    printf("  Nyquist compliant: %s\n", aa_ok ? "YES" : "NO");

    /* Step 4: ADC conversion */
    double v_ref = 5.0;
    unsigned adc_bits = 16;
    double lsb = v_ref / pow(2.0, adc_bits);
    unsigned adc_code = (unsigned)(v_with_noise / lsb);
    if (adc_code >= (1u << adc_bits)) adc_code = (1u << adc_bits) - 1;
    double v_digitized = adc_code * lsb;

    printf("\nStep 4 - ADC (16-bit, Vref=%.1fV):\n", v_ref);
    printf("  LSB = %.3f uV\n", lsb * 1e6);
    printf("  ADC code = %u\n", adc_code);
    printf("  V_digitized = %.6f V\n", v_digitized);

    /* Step 5: Software conversion back to uV at input */
    double v_recovered_uv = v_digitized / gain * 1e6;

    printf("\nStep 5 - Digital back to uV:\n");
    printf("  V_recovered = %.1f uV\n", v_recovered_uv);

    /* Step 6: CJC measurement (thermistor at cold junction) */
    /* Use Steinhart-Hart for accuracy (10k NTC at 35C => ~6530 Ohm) */
    double T_cj_measured = cjc_measure_thermistor(6530.0, &cjc_thermistor);

    printf("\nStep 6 - CJC sensor:\n");
    printf("  Thermistor R = %.1f Ohm\n", 6530.0);
    printf("  T_cj measured = %.2f C\n", T_cj_measured);

    /* Step 7: CJC compensation */
    double T_hot_calculated = cjc_compensate(v_recovered_uv, T_cj_measured, &tcK);

    double error = T_hot_calculated - t_hot_actual;

    printf("\nStep 7 - CJC compensation result:\n");
    printf("  T_hot calculated = %.2f C\n", T_hot_calculated);
    printf("  T_hot actual = %.2f C\n", t_hot_actual);
    printf("  Error = %.3f C\n", error);

    printf("\n=== Measurement Complete ===\n");

    if (fabs(error) < 5.0)
        printf("PASS: Error within 5C spec (NASA engine test tolerance)\n");
    else
        printf("WARN: Error exceeds tolerance\n");

    return 0;
}
