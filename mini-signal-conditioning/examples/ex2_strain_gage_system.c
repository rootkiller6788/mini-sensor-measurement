/**
 * ex2_strain_gage_system.c - Strain Gauge Measurement System
 *
 * Complete strain measurement chain: bridge excitation, amplification,
 * filtering, ADC, and strain calculation with nonlinearity correction.
 *
 * L6: Strain gauge measurement problem
 * L7: Structural health monitoring (Boeing 787 composite wing, F-35 airframe)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigcond_defs.h"
#include "sigcond_filter.h"
#include "sigcond_excitation.h"
#include "sigcond_bridge.h"

int main(void) {
    printf("=== Strain Gauge Measurement System ===\n\n");

    /* Application: Aerospace structural testing
     * Boeing 787 carbon fiber wing spar, strain gauge rosette
     * Gauge factor GF = 2.0, nominal R = 350 Ohm */

    double r_nominal = 350.0;
    double gauge_factor = 2.0;
    double max_power_mw = 2.0;  /* Composite: limit self-heating */

    printf("Scenario: Boeing 787 wing spar strain monitoring\n");
    printf("  Gauge: 350 Ohm, GF=2.0\n\n");

    /* Step 1: Optimize excitation voltage */
    double power_mw;
    double v_exc = bridge_optimal_excitation(r_nominal, max_power_mw, &power_mw);

    printf("Step 1 - Excitation optimization:\n");
    printf("  V_exc optimal = %.2f V\n", v_exc);
    printf("  Power per gauge = %.3f mW\n", power_mw);

    /* Step 2: Simulate strain measurement (1000 microstrain) */
    double strain_applied = 1000e-6;  /* 1000 ue */
    double v_bridge = bridge_output_strain(v_exc, gauge_factor, strain_applied, 1);

    printf("\nStep 2 - Bridge output:\n");
    printf("  Strain applied = %.0f ue\n", strain_applied * 1e6);
    printf("  V_bridge (linear) = %.6f V\n", v_bridge);

    /* Exact output with nonlinearity */
    double dR = r_nominal * gauge_factor * strain_applied;
    double v_exact = bridge_output_exact(v_exc, r_nominal, r_nominal,
                                          r_nominal + dR, r_nominal);
    double nl_error = bridge_nonlinearity_error(r_nominal, dR, 1);

    printf("  V_bridge (exact) = %.6f V\n", v_exact);
    printf("  Nonlinearity error = %.3f %%\n", nl_error);

    /* Step 3: Amplify (gain = 500) */
    double gain = 500.0;
    double v_amp = v_exact * gain;
    printf("\nStep 3 - Amplification (G=%.0f):\n", gain);
    printf("  V_amplified = %.4f V\n", v_amp);

    /* Step 4: Anti-alias filter */
    int order = filter_butterworth_order(100.0, 500.0, 0.1, 60.0);
    printf("\nStep 4 - Anti-alias filter:\n");
    printf("  Required order (fc=100Hz, fs=2kHz): %d\n", order);

    /* Step 5: ADC (24-bit sigma-delta) */
    double v_ref = 5.0;
    unsigned adc_bits = 24;
    double lsb = v_ref / pow(2.0, adc_bits);
    double noise_floor = 1e-6;  /* 1 uV RMS noise */
    double min_strain = bridge_min_detectable_strain(v_exc, gauge_factor,
                                                      noise_floor * 1e6, gain);

    printf("\nStep 5 - 24-bit ADC:\n");
    printf("  LSB = %.1f nV\n", lsb * 1e9);
    printf("  Min detectable strain = %.3f ue\n", min_strain * 1e6);

    /* Step 6: Calculate strain from digitized measurement */
    double v_digitized = round(v_amp / lsb) * lsb;
    double v_rti = v_digitized / gain;
    double strain_calc = bridge_strain_from_output(v_rti, v_exc,
                                                    gauge_factor, 1);
    double strain_corrected = bridge_nonlinearity_correct(v_rti, v_exc,
                                                           gauge_factor);

    printf("\nStep 6 - Strain calculation:\n");
    printf("  Strain (linear) = %.3f ue\n", strain_calc * 1e6);
    printf("  Strain (corrected) = %.3f ue\n", strain_corrected * 1e6);
    printf("  Error = %.3f ue\n", fabs(strain_corrected - strain_applied) * 1e6);

    printf("\n=== Airbus A350 / Boeing 787 structural test complete ===\n");
    return 0;
}
