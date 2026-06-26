#include "measurement_error.h"
#include "measurement_uncertainty.h"
#include "measurement_applications.h"
#include "measurement_filtering.h"
#include "measurement_compensation.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== ADC Error Analysis and Sensor Fusion Demo ===\n\n");

    /* Part 1: ADC Error Analysis */
    printf("--- ADC System Analysis ---\n");
    ADCSpec adc;
    adc_spec_init(&adc);
    adc.bits = 12;
    adc.fsr_voltage = 3.3;
    adc.inl_lsb = 1.5;
    adc.noise_rms_lsb = 0.25;

    double enob = adc_enob(&adc);
    double sinad = adc_sinad_db(&adc);
    double noise_lsb = adc_total_noise_lsb(&adc);

    printf("ADC: %u-bit, %.1fV FSR\n", (unsigned)adc.bits, adc.fsr_voltage);
    printf("  Total RMS noise: %.3f LSB\n", noise_lsb);
    printf("  SINAD: %.1f dB\n", sinad);
    printf("  ENOB: %.1f bits\n\n", enob);

    /* Part 2: Sensor fusion with Kalman filter */
    printf("--- Sensor Fusion with Kalman Filter ---\n");

    KalmanFilter1D kf;
    kf1d_init(&kf, 25.0, 2.0, 0.01, 0.5, 1.0);

    double true_temp = 25.0;
    double noisy_readings[] = {25.2, 25.8, 24.5, 25.1, 25.9,
                                24.8, 25.3, 24.9, 25.5, 25.0};
    int n = sizeof(noisy_readings) / sizeof(noisy_readings[0]);

    printf("True temperature: %.1f C\n", true_temp);
    printf("Raw readings (noisy):");
    for (int i = 0; i < n; i++) printf(" %.1f", noisy_readings[i]);

    double raw_mean = 0.0;
    for (int i = 0; i < n; i++) raw_mean += noisy_readings[i];
    raw_mean /= n;
    printf("\nRaw mean: %.2f C\n", raw_mean);

    /* Apply Kalman filter */
    printf("Kalman filtered:     ");
    for (int i = 0; i < n; i++) {
        double filtered = kf1d_update(&kf, noisy_readings[i]);
        if (i >= n - 5) printf(" %.2f", filtered);
    }
    printf("\nFinal estimate: %.2f C\n\n", kf.x_hat);

    /* Part 3: Temperature compensation demo */
    printf("--- Temperature Compensation ---\n");
    TempCompensationModel tcm = {25.0, 0.002, 0.0005, 0.0, 0};

    double raw_readings_t[] = {10.0, 10.15, 10.08, 10.22, 10.12};
    double temps[] = {25.0, 40.0, 35.0, 50.0, 30.0};

    printf("Raw -> Compensated (at various temperatures):\n");
    for (int i = 0; i < 5; i++) {
        double corrected = tempcomp_apply(raw_readings_t[i], temps[i], &tcm);
        printf("  %.2f @ %.0fC -> %.4f\n", raw_readings_t[i], temps[i], corrected);
    }

    return 0;
}
