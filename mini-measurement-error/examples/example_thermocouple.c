#include "measurement_error.h"
#include "measurement_uncertainty.h"
#include "measurement_applications.h"
#include <stdio.h>

int main(void) {
    printf("=== Thermocouple Measurement Error Budget ===\n");
    printf("Real-world Type K thermocouple system analysis\n\n");

    ThermocoupleSystemSpec spec;
    thermocouple_spec_init(&spec);

    printf("System Configuration:\n");
    printf("  Seebeck coefficient: %.3f mV/C\n", spec.seebeck_coeff_mv_per_c);
    printf("  CJC accuracy: %.1f C\n", spec.cjc_sensor_accuracy_c);
    printf("  Amplifier offset: %.1f uV\n", spec.amplifier_offset_uv);
    printf("  ADC resolution: %u bits\n\n", (unsigned)spec.adc_resolution_bits);

    UncertaintyBudget budget;
    double temps[] = {0.0, 100.0, 500.0, 1000.0};

    for (int i = 0; i < 4; i++) {
        thermocouple_error_budget(&spec, temps[i], &budget);
        printf("Measured Temperature: %.0f C\n", temps[i]);
        printf("  Combined uncertainty u_c: %.2f C\n",
               budget.combined_standard_uncertainty);
        printf("  Expanded uncertainty U (k=2): %.2f C\n",
               budget.expanded_uncertainty);
        printf("  Relative uncertainty: %.2f %%\n\n",
               100.0 * budget.expanded_uncertainty / temps[i]);
    }

    return 0;
}
