/**
 * ex2_thermocouple.c - Thermocouple Measurement with Cold Junction Compensation
 *
 * L6 Canonical Problem: K-type thermocouple readout with CJC.
 * Demonstrates: thermocouple voltage -> IA amplification ->
 * cold junction compensation -> temperature in degrees C.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "inst_amp_defs.h"
#include "inst_amp_topology.h"
#include "inst_amp_calibration.h"
#include "inst_amp_sensorif.h"

int main(void) {
    printf("=== Thermocouple Measurement Example ===\n\n");

    /* 1. K-type thermocouple configuration */
    thermocouple_config_t tc = {
        .type = TC_TYPE_K,
        .temp_ref_c = 25.0,    /* Cold junction at room temp */
        .temp_meas_c = 350.0   /* Hot junction at 350C */
    };

    double seebeck = thermocouple_seebeck_uv_per_c(tc.type);
    printf("Thermocouple: Type K, Seebeck = %.1f uV/C at 25C\n", seebeck);

    /* 2. Expected thermocouple voltage */
    double v_tc_uv = thermocouple_voltage_uv(&tc);
    printf("Expected V_TC = %.1f uV (for dT = %.0f C)\n",
           v_tc_uv, tc.temp_meas_c - tc.temp_ref_c);

    /* 3. IA configuration for microvolt signals */
    double G_needed = 5.0 / (v_tc_uv * 1e-6);  /* Map to 5V ADC */
    if (G_needed > 1000.0) G_needed = 1000.0;   /* Cap at 1000 */
    printf("Required IA gain: %.0f V/V\n", G_needed);

    /* 4. Simulate amplified output */
    double v_amplified = v_tc_uv * 1e-6 * G_needed;
    printf("Amplified output: %.3f V\n", v_amplified);

    /* 5. Cold junction compensation: recover temperature */
    double temp_meas = thermocouple_cjc_compensated(v_amplified, G_needed,
                                                     tc.type, tc.temp_ref_c);
    printf("Measured temperature (CJC): %.1f C\n", temp_meas);
    printf("Error: %.1f C (%.1f%%)\n", temp_meas - tc.temp_meas_c,
           100.0 * fabs(temp_meas - tc.temp_meas_c) / tc.temp_meas_c);

    /* 6. Compare different thermocouple types */
    printf("\nThermocouple type comparison at dT=100C:\n");
    thermocouple_type_t types[] = {TC_TYPE_K, TC_TYPE_J, TC_TYPE_T, TC_TYPE_E};
    const char *names[] = {"K", "J", "T", "E"};
    for (int i = 0; i < 4; i++) {
        double S = thermocouple_seebeck_uv_per_c(types[i]);
        printf("  Type %s: S=%.1f uV/C, V(100C)=%.1f uV\n",
               names[i], S, S * 100.0);
    }

    printf("\n=== Thermocouple Measurement Complete ===\n");
    return 0;
}