#include "measurement_filtering.h"
#include "measurement_statistics.h"
#include "measurement_error.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Sensor Fusion and Monte Carlo Uncertainty Demo ===\n\n");

    printf("--- Multi-Sensor Fusion ---\n");
    double readings[] = {10.05, 9.95, 10.12, 9.98};
    double vars[] = {0.01, 0.04, 0.09, 0.02};
    double fused_var;

    double fused = sensor_fuse_inverse_variance(readings, vars, 4, &fused_var);
    printf("Sensor readings: %.2f, %.2f, %.2f, %.2f\n",
           readings[0], readings[1], readings[2], readings[3]);
    printf("Individual uncertainties: %.3f, %.3f, %.3f, %.3f\n",
           sqrt(vars[0]), sqrt(vars[1]), sqrt(vars[2]), sqrt(vars[3]));
    printf("Fused estimate: %.4f +/- %.4f\n", fused, sqrt(fused_var));
    printf("Variance reduction: %.1fx\n\n", sensor_fuse_variance_reduction(4));

    printf("--- Grubbs Outlier Test ---\n");
    double data[] = {100.1, 100.3, 99.8, 100.0, 105.2, 100.2, 99.9};
    int outlier;
    double G = stats_grubbs_test(data, 7, &outlier);
    printf("Data: ");
    for (int i = 0; i < 7; i++) printf("%.1f ", data[i]);
    printf("\nGrubbs G = %.3f, suspected outlier at index %d (value=%.1f)\n",
           G, outlier, data[outlier]);
    printf("Critical G(0.05, 7) ~ 2.02\n\n");

    printf("--- Modified Z-Score Method ---\n");
    int flags[7];
    stats_modified_zscore_outliers(data, 7, flags);
    printf("Outlier flags: ");
    for (int i = 0; i < 7; i++) printf("%d ", flags[i]);
    printf("\n(|Modified Z| > 3.5 => outlier)\n\n");

    printf("--- Moving Average vs EMA Comparison ---\n");
    MovingAvgFilter mavg;
    mavg_init(&mavg, 5);
    EMAFilter ema;
    ema_init(&ema, 0.3);

    printf("Step  Raw     MA(5)   EMA(0.3)\n");
    double signal[] = {10.0, 10.5, 9.8, 10.2, 11.0, 10.1, 9.9, 10.3};
    for (int i = 0; i < 8; i++) {
        double ma = mavg_update(&mavg, signal[i]);
        double em = ema_update(&ema, signal[i]);
        printf("%4d  %6.2f  %6.2f  %6.2f\n", i+1, signal[i], ma, em);
    }
    mavg_free(&mavg);

    return 0;
}
