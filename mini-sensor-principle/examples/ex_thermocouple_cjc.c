/**
 * @file ex_thermocouple_cjc.c
 * @brief Example: Type K thermocouple measurement with cold junction compensation
 *
 * This example demonstrates:
 *   1. Thermocouple EMF computation per NIST ITS-90
 *   2. Cold junction temperature measurement via thermistor
 *   3. Cold junction compensation to get true hot junction temperature
 *   4. Error analysis without CJC vs with CJC
 *
 * This is a canonical L6 problem: thermocouple-based temperature measurement
 * with analog cold junction compensation.
 */

#include <stdio.h>
#include <math.h>
#include "../include/sensor_transfer.h"

int main(void) {
    printf("=== Type K Thermocouple with Cold Junction Compensation ===\n\n");

    /* ──── Initialize Type K thermocouple model ──── */
    tc_model_t tc;
    if (tc_model_init(&tc, TC_TYPE_K) != 0) {
        printf("ERROR: Failed to initialize thermocouple model\n");
        return 1;
    }
    printf("Thermocouple Type K initialized\n");
    printf("  Range: %.0f °C to %.0f °C\n", tc.t_min_C, tc.t_max_C);
    printf("  Sensitivity at 25°C: %.1f µV/°C\n", tc.sensitivity_typical);

    /* ──── Initialize NTC thermistor for cold junction ──── */
    thermistor_model_t tm;
    thermistor_model_init_beta(&tm, 10000.0, 298.15, 3977.0);
    printf("\nCold junction sensor: NTC thermistor (10kΩ @ 25°C, β=3977K)\n");

    /* ──── Scenario: Measure hot junction at various temperatures ──── */
    printf("\n--- Measurement Simulation ---\n");
    printf("┌─────────┬──────────┬───────────┬────────────┬──────────┐\n");
    printf("│ T_hot   │ EMF(T_h) │ T_cj      │ V_measured │ T_calc   │\n");
    printf("│ (°C)    │ (µV)     │ (°C)      │ (µV)       │ (°C)     │\n");
    printf("├─────────┼──────────┼───────────┼────────────┼──────────┤\n");

    double test_temps[] = {0.0, 25.0, 100.0, 250.0, 500.0, 750.0, 1000.0};
    double T_cj = 22.5;  /* cold junction at room temperature */

    for (int i = 0; i < 7; i++) {
        double T_hot = test_temps[i];

        /* EMF that would be generated if cold junction were at 0°C */
        double E_0C = tc_emf_from_T(&tc, T_hot, 0.0);

        /* EMF that is actually measured (cold junction at T_cj) */
        double V_measured = tc_emf_from_T(&tc, T_hot, T_cj);

        /* Cold junction compensation: recover T_hot from V_measured */
        double T_calculated = tc_T_from_emf(&tc, V_measured, T_cj, 0.01);

        printf("│ %7.0f │ %8.1f │ %9.1f │ %10.1f │ %8.1f │\n",
               T_hot, E_0C, T_cj, V_measured, T_calculated);
    }
    printf("└─────────┴──────────┴───────────┴────────────┴──────────┘\n");

    /* ──── Demonstrate CJC importance ──── */
    printf("\n--- CJC Error Analysis ---\n");
    double T_true = 500.0;
    /* V_correct is the EMF measured with cold junction at T_cj */
    double V_correct = tc_emf_from_T(&tc, T_true, T_cj);

    double T_wrong = tc_T_from_emf(&tc, V_correct, 0.0, 0.01); /* skip CJC */
    double T_right = tc_T_from_emf(&tc, V_correct, T_cj, 0.01);

    printf("True temperature: %.1f °C\n", T_true);
    printf("If CJC ignored (T_ref assumed 0°C): %.1f °C (error: %.1f °C)\n",
           T_wrong, T_wrong - T_true);
    printf("With correct CJC at %.1f °C: %.1f °C (error: %.1f °C)\n",
           T_cj, T_right, T_right - T_true);
    printf("\nMoral: CJC is essential — ignoring it causes ~T_cj error!\n");

    /* ──── Law of Intermediate Metals ──── */
    printf("\n--- Law of Intermediate Metals Check ---\n");
    int law_ok = tc_law_intermediate_metals_check(V_correct, T_true, T_cj, T_cj);
    printf("Copper traces at same temperature (%.1f°C): %s\n",
           T_cj, law_ok == 0 ? "LAW SATISFIED" : "VIOLATION");
    law_ok = tc_law_intermediate_metals_check(V_correct, T_true, T_cj, T_cj + 10.0);
    printf("Copper traces at different temps (%.1f vs %.1f°C): %s\n",
           T_cj, T_cj + 10.0, law_ok == 0 ? "OK" : "SIGNIFICANT ERROR POSSIBLE");

    return 0;
}
