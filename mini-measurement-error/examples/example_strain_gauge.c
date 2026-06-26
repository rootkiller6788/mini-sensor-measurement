#include "measurement_error.h"
#include "measurement_uncertainty.h"
#include "measurement_applications.h"
#include "measurement_compensation.h"
#include <stdio.h>

int main(void) {
    printf("=== Strain Gauge Bridge Uncertainty Analysis ===\n\n");

    StrainGaugeSystemSpec spec;
    straingauge_spec_init(&spec);

    printf("System Configuration:\n");
    printf("  Excitation: %.1f V (accuracy: %.2f%%)\n",
           spec.excitation_voltage, spec.excitation_accuracy_pct);
    printf("  Gauge Factor: %.2f (uncertainty: %.1f%%)\n",
           spec.gauge_factor, spec.gf_uncertainty_pct);
    printf("  ADC: %u bits\n\n", (unsigned)spec.adc_bits);

    /* Simulate a measurement at 1000 microstrain */
    double V_out = 1000e-6 * spec.excitation_voltage * spec.gauge_factor / 4.0;
    double measured_strain = straingauge_strain_from_voltage(
        V_out, spec.excitation_voltage, spec.gauge_factor);

    printf("Simulated measurement:\n");
    printf("  Bridge output: %.6f V\n", V_out);
    printf("  Apparent strain: %.1f microstrain\n\n", measured_strain * 1e6);

    UncertaintyBudget budget;
    straingauge_uncertainty_budget(&spec, measured_strain, &budget);

    printf("Uncertainty Budget:\n");
    for (size_t i = 0; i < budget.num_components; i++) {
        printf("  %-20s  u=%.3e  (Type %s)\n",
               budget.components[i].name,
               budget.components[i].standard_uncertainty,
               budget.components[i].eval_type == EVAL_TYPE_A ? "A" : "B");
    }

    printf("\n  Combined uncertainty: %.2f microstrain\n",
           budget.combined_standard_uncertainty * 1e6);
    printf("  Expanded uncertainty (k~2): %.2f microstrain\n",
           budget.expanded_uncertainty * 1e6);
    printf("  Result: %.1f +/- %.1f microstrain\n",
           measured_strain * 1e6, budget.expanded_uncertainty * 1e6);

    return 0;
}
