/**
 * example_sensor_array.c — Electronic nose / chemical sensor array
 *
 * L7 Application: Gas/odor classification using chemiresistor sensor arrays.
 * L8 Advanced Topic: Pattern recognition, calibration transfer across batches.
 *
 * Demonstrates:
 * - Chemiresistor array configuration and baseline measurement
 * - Response normalization
 * - Euclidean distance-based nearest-neighbor classification
 * - Calibration transfer between sensor batches
 * - Drift compensation and baseline tracking
 *
 * Reference: Persaud & Dodd, Nature 299 (1982) 352-355 (first e-nose)
 *            Gardner & Bartlett, "Electronic Noses" (1999)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "biosensor_types.h"
#include "biosensor_signal.h"
#include "biosensor_calibration.h"

/* Known odor profiles (normalized responses for 8-sensor array) */
static const double profile_ethanol[8]    = {0.15, 0.45, 0.08, 0.03, 0.20, 0.55, 0.12, 0.05};
static const double profile_acetone[8]    = {0.30, 0.15, 0.05, 0.02, 0.40, 0.25, 0.08, 0.03};
static const double profile_ammonia[8]    = {0.05, 0.02, 0.60, 0.35, 0.04, 0.01, 0.50, 0.42};
static const double profile_toluene[8]    = {0.08, 0.10, 0.03, 0.02, 0.55, 0.60, 0.10, 0.07};
static const char   *profile_names[4]     = {"Ethanol", "Acetone", "Ammonia", "Toluene"};

static const double *profiles[4] = {
    profile_ethanol, profile_acetone, profile_ammonia, profile_toluene
};

int main(void)
{
    printf("=== Electronic Nose — Chemical Sensor Array ===\n\n");

    /* Configure 8-sensor MOS chemiresistor array */
    SensorArray sa;
    sensor_array_init(&sa, 8);
    sa.sampling_rate_hz = 10.0;
    sa.array_name = "MOS-ENose-V2";

    /* Configure sensor elements */
    const char *materials[] = {
        "SnO2 (undoped)", "SnO2 (Pd-doped)", "WO3", "ZnO",
        "In2O3", "TiO2 (Pt-doped)", "MoO3", "CuO"
    };
    double baseline_r[] = {2.5e6, 1.2e6, 3.8e6, 5.1e6, 0.89e6, 4.2e6, 7.5e6, 1.8e6};
    double sensitivities[] = {0.15, 0.45, 0.08, 0.03, 0.20, 0.55, 0.12, 0.05};

    printf("--- Sensor Configuration ---\n");
    printf("Sensor | Material           | R0 (MΩ)   | Sensitivity (ΔR/R0 per ppm)\n");
    printf("-------|-------------------|-----------|-----------------------------\n");
    for (int i = 0; i < 8; i++) {
        sa.sensors[i].sensor_id = i;
        sa.sensors[i].sensor_material = materials[i];
        sa.sensors[i].baseline_resistance_ohm = baseline_r[i];
        sa.sensors[i].sensitivity_ppm = sensitivities[i];
        sa.sensors[i].operating_temp_celsius = 300.0 + i * 15.0;
        printf("  %2d   | %-17s | %6.1f   | %.2f\n",
               i, materials[i], baseline_r[i] / 1.0e6, sensitivities[i]);
    }
    sa.sensor_count = 8;

    /* Simulate response to an unknown vapor */
    printf("\n--- Unknown Vapor Exposure ---\n");

    /* Simulated raw responses (with noise) for a vapor closest to acetone */
    double raw_responses[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    double noise_level = 0.02;  /* 2% noise */
    srand(42);

    printf("Sensor | Raw R (MΩ) | Baseline (MΩ) | ΔR/R0  | Normalized\n");
    printf("-------|------------|---------------|--------|----------\n");
    for (int i = 0; i < 8; i++) {
        /* Generate response: shifts toward acetone profile + noise */
        double response = 1.0 - 0.3 * profile_acetone[i];
        double noise = (rand() / (double)RAND_MAX - 0.5) * 2.0 * noise_level;
        double r_raw = baseline_r[i] * (response + noise);
        raw_responses[i] = r_raw;

        double delta_r = (r_raw - baseline_r[i]) / baseline_r[i];
        printf("  %2d   |  %8.3f   |  %8.3f    | %+6.4f | %+7.4f\n",
               i, r_raw / 1.0e6, baseline_r[i] / 1.0e6, delta_r, delta_r);
    }

    /* Normalize to fractional baseline */
    double normalized[8];
    sensor_array_normalize(raw_responses, baseline_r, normalized, 8);

    printf("\n--- Normalized Response Pattern ---\n");
    for (int i = 0; i < 8; i++) {
        printf("  Sensor %d: %+6.4f", i, normalized[i]);
        if (i == 3) printf("\n");
    }
    if (8 % 4 != 0) printf("\n");

    /* Nearest-neighbor classification */
    printf("\n--- Nearest-Neighbor Classification ---\n");
    printf("Target Vapor    | Euclidean Distance | Match?\n");
    printf("----------------|-------------------|--------\n");

    int best_idx = -1;
    double best_dist = INFINITY;
    for (int p = 0; p < 4; p++) {
        double dist = sensor_array_euclidean_distance(normalized,
                                                       profiles[p], 8);
        const char *marker = "";
        printf("  %-13s |  %12.6f     | %s\n", profile_names[p], dist, marker);
        if (dist < best_dist) { best_dist = dist; best_idx = p; }
    }

    printf("\n  => Best match: %s (distance = %.6f)\n",
           profile_names[best_idx], best_dist);

    /* Demonstrate baseline drift and correction */
    printf("\n--- Baseline Drift Simulation ---\n");
    printf("Time (s) | Raw Signal | Drift Est. | Corrected\n");
    printf("---------|-----------|-----------|----------\n");

    double time[50], raw_sig[50], corrected[50];
    for (int t = 0; t < 10; t++) {
        time[t] = t * 5.0;  /* every 5 seconds */
        raw_sig[t] = 1.0 + 0.3 * (1.0 - exp(-time[t] / 30.0)) +
                     0.02 * sin(time[t] * 0.5);
    }
    exponential_drift_correction(raw_sig, time, 10, 0.3, 30.0, corrected);

    for (int t = 0; t < 10; t++) {
        double drift_est = 0.3 * (1.0 - exp(-time[t] / 30.0));
        printf("  %6.1f  |  %7.4f   | %7.4f   | %7.4f\n",
               time[t], raw_sig[t], drift_est, corrected[t]);
    }

    printf("\nNote: After drift correction, the true signal (~1.0) is recovered.\n");

    /* Demonstrate calibration transfer */
    printf("\n--- Calibration Transfer (Batch A → Batch B) ---\n");
    int n_chan = 4;
    double master_resp[] = {1.0, 0.5, 0.3, 0.2,  0.8, 0.9, 0.4, 0.3,  0.6, 1.2, 0.5, 0.4};
    double slave_resp[]  = {1.1, 0.55, 0.33, 0.22, 0.88, 0.99, 0.44, 0.33, 0.66, 1.32, 0.55, 0.44};
    double *transfer_mat = calloc(n_chan * n_chan, sizeof(double));

    double rmse = pds_calibration_transfer(master_resp, slave_resp,
                                            3, n_chan, 1, transfer_mat);
    printf("  Channels:          %d\n", n_chan);
    printf("  Transfer standards: 3\n");
    printf("  PDS RMSE:          %.4f\n", rmse);
    printf("  Transfer matrix (first 2 columns):\n");
    for (int r = 0; r < n_chan; r++) {
        printf("    [");
        for (int c = 0; c < 2; c++) {
            printf(" %+7.4f", transfer_mat[r * n_chan + c]);
        }
        printf(" ... ]\n");
    }

    /* Clean up */
    free(transfer_mat);
    sensor_array_free(&sa);

    printf("\n=== Electronic nose simulation complete ===\n");
    return 0;
}
