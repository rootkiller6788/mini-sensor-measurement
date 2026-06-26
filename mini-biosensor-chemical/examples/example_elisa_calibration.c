/**
 * example_elisa_calibration.c — ELISA immunoassay calibration workflow
 *
 * L6 Canonical Problem: 4PL and 5PL calibration for sandwich ELISA.
 * L7 Application: Clinical immunoassay data reduction per FDA guidelines.
 *
 * Demonstrates:
 * - 4PL nonlinear calibration fitting
 * - LOD/LOQ determination from calibration data
 * - Outlier detection in replicate wells
 * - Weighted vs unweighted regression comparison
 * - Inverse prediction (concentration from signal)
 *
 * Reference: Findlay et al., J. Pharm. Biomed. Anal. 21 (2000) 1249-1273
 */

#include <stdio.h>
#include <math.h>
#include "biosensor_types.h"
#include "biosensor_calibration.h"
#include "biosensor_optical.h"
#include "biosensor_signal.h"

int main(void)
{
    printf("=== ELISA Immunoassay Calibration Demo ===\n\n");

    /* Simulated human IL-6 ELISA data (pg/mL)
     * Signal in OD450, 8 standards, 2 replicates each */
    double conc[]  = {0.0, 3.1, 6.25, 12.5, 25.0, 50.0, 100.0, 200.0};
    double sig_mean[]  = {0.045, 0.098, 0.178, 0.342, 0.658, 1.182, 1.742, 2.051};
    double sig_std[]   = {0.005, 0.008, 0.012, 0.018, 0.031, 0.052, 0.071, 0.089};
    int n_standards = 8;

    printf("--- Calibration Standards (IL-6, pg/mL) ---\n");
    printf("Conc (pg/mL) | Mean OD450 | Std Dev | CV%%\n");
    printf("-------------|-----------|---------|------\n");
    for (int i = 0; i < n_standards; i++) {
        double cv = (sig_mean[i] > 0.0) ? (sig_std[i] / sig_mean[i] * 100.0) : 0.0;
        printf("  %6.1f     |  %6.3f   | %6.3f  | %4.1f%%\n",
               conc[i], sig_mean[i], sig_std[i], cv);
    }

    /* LOD estimation (3σ method) */
    double blank_std = sig_std[0];
    printf("\n--- Detection Limits ---\n");
    printf("  Blank SD:           %.4f OD\n", blank_std);

    /* Estimate sensitivity from lowest non-zero standard */
    double sensitivity = (sig_mean[1] - sig_mean[0]) / (conc[1] - conc[0]);
    printf("  Sensitivity:        %.4f OD/(pg/mL)\n", sensitivity);
    double lod = lod_three_sigma(blank_std, sensitivity);
    double loq = loq_ten_sigma(blank_std, sensitivity);
    printf("  LOD (3.3σ):         %.2f pg/mL\n", lod);
    printf("  LOQ (10σ):          %.2f pg/mL\n", loq);

    /* Fit 4PL model */
    printf("\n--- 4-Parameter Logistic Fit ---\n");
    CalibrationStandards cs;
    calibration_standards_init(&cs, n_standards, 2);
    for (int i = 0; i < n_standards; i++) {
        cs.concentrations[i] = conc[i];
        cs.signals[i]        = sig_mean[i];
        cs.signal_stds[i]    = sig_std[i];
    }

    FourPLLogisticFit fit;
    fourpl_fit_init(&fit);
    fit.a = sig_mean[0];      /* blank */
    fit.d = sig_mean[7] + 0.1; /* upper asymptote estimate */

    /* Get initial b, c from log-logistic linearization */
    double r2_lin = fourpl_linearized_initial_guess(&cs, &fit);
    printf("  Linearized init R²:  %.4f\n", r2_lin);
    printf("  Initial guess: a=%.3f, b=%.2f, c=%.1f, d=%.3f\n",
           fit.a, fit.b, fit.c, fit.d);

    int iters = fourpl_fit(&cs, &fit, 100, 0.01);
    printf("  LM iterations:       %d\n", iters);
    printf("  Final fit:\n");
    printf("    a (lower asymptote):  %.4f\n", fit.a);
    printf("    b (slope factor):     %.2f\n", fit.b);
    printf("    c (EC50, pg/mL):     %.2f\n", fit.c);
    printf("    d (upper asymptote):  %.4f\n", fit.d);
    printf("    R²:                   %.4f\n", fit.r_squared);
    printf("    RMSE:                 %.5f\n", fit.rmse);

    /* Back-fit: evaluate model at all standards */
    printf("\n--- Model Evaluation ---\n");
    printf("Conc | Measured | Predicted | Residual | %%Error\n");
    printf("-----|---------|----------|---------|-------\n");
    for (int i = 0; i < n_standards; i++) {
        double pred = elisa_fourpl_predict(&fit, conc[i]);
        double resid = sig_mean[i] - pred;
        double pct_err = (pred > 0.01) ? (resid / pred * 100.0) : 0.0;
        printf(" %5.1f | %6.3f  |  %6.3f  | %+6.4f | %+5.1f%%\n",
               conc[i], sig_mean[i], pred, resid, pct_err);
    }

    /* Unknown sample prediction */
    printf("\n--- Unknown Sample Analysis ---\n");
    double unknown_signal = 1.050;  /* OD450 */
    double unknown_conc  = elisa_fourpl_inverse(&fit, unknown_signal);
    printf("  Measured OD450:  %.3f\n", unknown_signal);
    printf("  Estimated IL-6:  %.1f pg/mL\n", unknown_conc);

    if (unknown_conc < lod) {
        printf("  Result: BELOW DETECTION LIMIT\n");
    } else if (unknown_conc < loq) {
        printf("  Result: DETECTED (below quantification limit)\n");
    } else {
        printf("  Result: QUANTIFIED\n");
        /* Clinical context: IL-6 normal < 5 pg/mL, sepsis > 100 pg/mL */
        if (unknown_conc < 5.0) {
            printf("  Clinical: Normal IL-6 level\n");
        } else if (unknown_conc < 50.0) {
            printf("  Clinical: Mildly elevated (inflammation)\n");
        } else if (unknown_conc < 100.0) {
            printf("  Clinical: Significantly elevated (possible sepsis)\n");
        } else {
            printf("  Clinical: Critically elevated (sepsis/septic shock)\n");
        }
    }

    /* Weighted mean of replicates (simulated) */
    printf("\n--- Replicate Analysis (Weighted Mean) ---\n");
    double rep_concs[] = {42.3, 38.7, 45.1, 40.2};  /* pg/mL, 4 wells */
    double rep_vars[]  = {2.5, 3.1, 1.8, 2.9};       /* variance */
    double wm = elisa_weighted_mean(rep_concs, rep_vars, 4);
    double cv_rep = coefficient_of_variation(rep_concs, 4);
    printf("  Replicates:  %.1f, %.1f, %.1f, %.1f pg/mL\n",
           rep_concs[0], rep_concs[1], rep_concs[2], rep_concs[3]);
    printf("  Weighted mean: %.2f pg/mL\n", wm);
    printf("  CV:            %.1f%%\n", cv_rep);

    /* Outlier test */
    int outlier_idx;
    bool has_outlier = grubbs_outlier_test(rep_concs, 4, 0.05, &outlier_idx);
    printf("  Outlier test:  %s", has_outlier ? "SUSPECT FOUND" : "PASSED");
    if (has_outlier) printf(" (index %d)", outlier_idx);
    printf("\n");

    /* Clean up */
    calibration_standards_free(&cs);

    printf("\n=== ELISA calibration complete ===\n");
    return 0;
}
