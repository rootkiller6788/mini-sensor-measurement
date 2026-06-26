/**
 * example_glucose_biosensor.c — End-to-end glucose biosensor simulation
 *
 * L6 Canonical Problem: Clark-type glucose biosensor operating across
 * the clinical range from hypoglycemic to diabetic thresholds.
 *
 * Demonstrates:
 * - Amperometric glucose sensing with Michaelis-Menten kinetics
 * - O₂ co-substrate limitation
 * - Calibration curve generation
 * - Clinical glucose classification from sensor current
 *
 * Reference: Clark & Lyons, Ann. NY Acad. Sci. 102 (1962) 29-45
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "biosensor_types.h"
#include "biosensor_electrochemical.h"
#include "biosensor_kinetics.h"
#include "biosensor_calibration.h"

int main(void)
{
    printf("=== Glucose Biosensor Simulation ===\n");
    printf("Model: GOx-based amperometric biosensor\n");
    printf("Clinical thresholds (fasting, plasma):\n");
    printf("  Hypoglycemic:  < 3.9 mM (70 mg/dL)\n");
    printf("  Normal:        3.9 - 5.6 mM\n");
    printf("  Prediabetic:   5.6 - 7.0 mM\n");
    printf("  Diabetic:      > 7.0 mM\n\n");

    /* Sensor configuration */
    double enzyme_activity  = 2.0;       /* U/cm² */
    double electrode_area   = 0.0707;    /* cm² (3 mm diameter) */
    double oxygen_conc      = 0.21;      /* mM (air-saturated buffer) */
    double km_apparent      = 8.0;       /* mM, immobilized GOx */

    /* Generate calibration curve */
    printf("--- Calibration Curve ---\n");
    printf("Glucose (mM) | Current (nA) | Classification\n");
    printf("-------------|-------------|----------------\n");

    double glucose_levels[] = {1.0, 2.5, 3.9, 5.0, 5.6, 7.0, 10.0, 15.0, 20.0};
    int n_levels = 9;
    double *currents = malloc((size_t)n_levels * sizeof(double));

    for (int i = 0; i < n_levels; i++) {
        currents[i] = glucose_biosensor_current(glucose_levels[i],
                                                 enzyme_activity,
                                                 electrode_area,
                                                 oxygen_conc);
        const char *class_str;
        if (glucose_levels[i] < 3.9)      class_str = "HYPO";
        else if (glucose_levels[i] <= 5.6) class_str = "NORMAL";
        else if (glucose_levels[i] <= 7.0) class_str = "PREDIABETIC";
        else                               class_str = "DIABETIC";

        printf("  %6.1f mM   |  %6.1f nA   | %s\n",
               glucose_levels[i], currents[i], class_str);
    }

    /* Fit linear calibration from low-glucose region */
    printf("\n--- Linear Calibration (0-7 mM range) ---\n");
    CalibrationStandards cs;
    calibration_standards_init(&cs, 5, 1);

    cs.concentrations[0] = 1.0;  cs.signals[0] = currents[0];
    cs.concentrations[1] = 2.5;  cs.signals[1] = currents[1];
    cs.concentrations[2] = 3.9;  cs.signals[2] = currents[2];
    cs.concentrations[3] = 5.0;  cs.signals[3] = currents[3];
    cs.concentrations[4] = 7.0;  cs.signals[4] = currents[5];

    LinearCalibration lc;
    double r2 = linear_calibration_fit(&cs, &lc);
    printf("  Slope:      %.2f nA/mM\n", lc.slope);
    printf("  Intercept:  %.2f nA\n", lc.intercept);
    printf("  R²:         %.4f\n", r2);
    printf("  Sensitivity: %.1f nA/mM\n", lc.slope);

    /* Estimate unknown sample */
    printf("\n--- Unknown Sample Analysis ---\n");
    double sample_current = 350.0;  /* nA, simulated measurement */
    double conc, lower, upper;
    linear_calibration_predict(&lc, sample_current, &cs, &conc, &lower, &upper);
    printf("  Measured current: %.1f nA\n", sample_current);
    printf("  Estimated glucose: %.1f mM [%.1f, %.1f]\n", conc, lower, upper);

    const char *diag;
    if (conc < 3.9)      diag = "CRITICAL LOW - Possible Hypoglycemia";
    else if (conc <= 5.6) diag = "NORMAL - Fasting glucose within range";
    else if (conc <= 7.0) diag = "ELEVATED - Prediabetic range, follow-up advised";
    else                  diag = "HIGH - Diabetic range, consult physician";
    printf("  Assessment: %s\n", diag);

    /* Demonstrate oxygen limitation */
    printf("\n--- O₂ Co-Substrate Limitation ---\n");
    printf("O₂ (mM) | Current at 10mM glucose (nA)\n");
    printf("--------|---------------------------\n");
    double o2_levels[] = {0.05, 0.10, 0.15, 0.21, 0.30};
    for (int i = 0; i < 5; i++) {
        double i_limited = glucose_biosensor_current(10.0, enzyme_activity,
                                                      electrode_area, o2_levels[i]);
        printf("  %.2f   |  %.1f\n", o2_levels[i], i_limited);
    }
    printf("\nNote: Low O₂ (venous blood ~0.05 mM) significantly reduces signal.\n");

    /* Michaelis-Menten fit to the full calibration */
    printf("\n--- Michaelis-Menten Kinetic Fit ---\n");
    MichaelisMentenParams mm;
    michaelis_menten_params_init(&mm, km_apparent, currents[8] * 2.0);
    double km_fit, vmax_fit;
    int iters = michaelis_menten_nls_fit(glucose_levels, currents, n_levels,
                                          8.0, 2000.0, &km_fit, &vmax_fit,
                                          100, 1.0e-8);
    printf("  Iterations:  %d\n", iters);
    printf("  Km (apparent): %.1f mM\n", km_fit);
    printf("  Imax:         %.1f nA\n", vmax_fit);
    printf("  Linear range: %.1f mM (≈10%% of Km)\n", km_fit * 0.1);

    /* Clean up */
    calibration_standards_free(&cs);
    free(currents);

    printf("\n=== Simulation complete ===\n");
    return 0;
}
