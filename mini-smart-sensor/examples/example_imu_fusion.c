/**
 * @file    example_imu_fusion.c
 * @brief   L6: IMU orientation estimation via sensor fusion
 *
 * Demonstrates sensor fusion for MEMS IMU orientation estimation:
 *   1. Kalman filter (gyro + accelerometer fusion)
 *   2. Complementary filter (computationally lighter alternative)
 *   3. Comparison of both approaches under simulated conditions
 *
 * This is a canonical problem in smart sensors: combining
 * complementary sensors (gyroscope = good short-term,
 * accelerometer = good long-term reference) to estimate
 * orientation angles.
 *
 * Reference:
 *   Valenti, R.G. et al. (2015), Sensors 15(8):19302-19330
 *   Mahony, R. et al. (2008), IEEE TAC 53(5):1203-1218
 *   MPU-6050 MEMS IMU Datasheet, InvenSense/TDK
 */

#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <math.h>
#include <string.h>
#include "smart_sensor.h"
#include "sensor_fusion.h"
#include <stdlib.h>

/* Simulate IMU data: constant tilt (15 deg roll, 10 deg pitch) with noise */
static void simulate_imu_sample(double t, double *gx, double *gy, double *gz,
                                double *ax, double *ay, double *az) {
    double roll_true  = 15.0 * M_PI / 180.0;  /* 15 deg */
    double pitch_true = 10.0 * M_PI / 180.0;  /* 10 deg */

    (void)t;

    /* Gyroscope: zero rate (static) + noise */
    *gx = 0.0 + 0.001 * ((double)rand() / RAND_MAX - 0.5);  /* rad/s noise */
    *gy = 0.0 + 0.001 * ((double)rand() / RAND_MAX - 0.5);
    *gz = 0.0 + 0.001 * ((double)rand() / RAND_MAX - 0.5);

    /* Accelerometer: gravity vector at true orientation + noise */
    /* Gravity components for roll=15 deg, pitch=10 deg */
    double g_roll  = sin(roll_true);
    double g_pitch = -sin(pitch_true) * cos(roll_true);
    double g_down  = cos(pitch_true) * cos(roll_true);

    *ax = g_pitch  + 0.005 * ((double)rand() / RAND_MAX - 0.5);  /* g noise */
    *ay = g_roll   + 0.005 * ((double)rand() / RAND_MAX - 0.5);
    *az = g_down   + 0.005 * ((double)rand() / RAND_MAX - 0.5);
}

int main(void) {
    printf("============================================================\n");
    printf("  L6: IMU Sensor Fusion — Orientation Estimation\n");
    printf("  Kalman Filter vs Complementary Filter\n");
    printf("============================================================\n\n");

    srand(12345);  /* Fixed seed for reproducibility */

    /* ============================================================
     * Part 1: Kalman Filter IMU
     * ============================================================ */

    printf("--- Part 1: Kalman Filter IMU Orientation ---\n\n");

    ss_kalman_3d_imu_t kalman;
    ss_fusion_kalman_imu_init(&kalman, 0.01,  /* dt = 10 ms */
                               0.001,         /* gyro noise (rad^2/s^2) */
                               0.01);         /* accel noise (g^2) */

    double roll_true  = 15.0 * M_PI / 180.0;
    double pitch_true = 10.0 * M_PI / 180.0;

    printf("True orientation: Roll = %.1f deg, Pitch = %.1f deg\n",
           roll_true * 180.0 / M_PI, pitch_true * 180.0 / M_PI);
    printf("Sampling: dt = 10 ms, 200 iterations = 2 seconds\n\n");

    printf("Kalman Filter convergence:\n");
    printf("  %-8s %-14s %-14s %-14s %-14s\n",
           "Iter", "Roll (deg)", "Roll Err", "Pitch (deg)", "Pitch Err");
    printf("  %-8s %-14s %-14s %-14s %-14s\n",
           "----", "----------", "----------", "-----------", "-----------");

    {
        int i;
        for (i = 0; i <= 200; i += 20) {
            /* Run several steps between prints */
            int j;
            double gx, gy, gz, ax, ay, az;
            for (j = 0; j < 20 && (i + j) < 200; j++) {
                simulate_imu_sample(0.01 * (i + j), &gx, &gy, &gz, &ax, &ay, &az);
                ss_fusion_kalman_imu_predict(&kalman, gx, gy, gz);
                ss_fusion_kalman_imu_update(&kalman, ax, ay, az);
            }

            double roll_est, pitch_est;
            ss_fusion_kalman_imu_get_angles(&kalman, &roll_est, &pitch_est);

            printf("  %-8d %-14.3f %-14.4f %-14.3f %-14.4f\n",
                   i + 20,
                   roll_est * 180.0 / M_PI,
                   fabs(roll_est - roll_true) * 180.0 / M_PI,
                   pitch_est * 180.0 / M_PI,
                   fabs(pitch_est - pitch_true) * 180.0 / M_PI);
        }
    }

    /* Final state */
    {
        double roll_final, pitch_final;
        ss_fusion_kalman_imu_get_angles(&kalman, &roll_final, &pitch_final);
        printf("\nFinal Kalman estimate:\n");
        printf("  Roll  = %.3f deg (error = %.4f deg)\n",
               roll_final * 180.0 / M_PI,
               fabs(roll_final - roll_true) * 180.0 / M_PI);
        printf("  Pitch = %.3f deg (error = %.4f deg)\n",
               pitch_final * 180.0 / M_PI,
               fabs(pitch_final - pitch_true) * 180.0 / M_PI);
    }

    /* ============================================================
     * Part 2: Complementary Filter IMU
     * ============================================================ */

    printf("\n--- Part 2: Complementary Filter IMU Orientation ---\n\n");

    ss_complementary_filter_t comp;
    ss_fusion_complementary_init(&comp, 0.98, 0.01);

    printf("Complementary filter (alpha=0.98) convergence:\n");
    printf("  %-8s %-14s %-14s %-14s %-14s\n",
           "Iter", "Roll (deg)", "Roll Err", "Pitch (deg)", "Pitch Err");
    printf("  %-8s %-14s %-14s %-14s %-14s\n",
           "----", "----------", "----------", "-----------", "-----------");

    {
        int i;
        for (i = 0; i <= 200; i += 20) {
            int j;
            double gx, gy, gz, ax, ay, az;
            for (j = 0; j < 20 && (i + j) < 200; j++) {
                simulate_imu_sample(0.01 * (i + j), &gx, &gy, &gz, &ax, &ay, &az);
                ss_fusion_complementary_update(&comp, gx, gy, gz, ax, ay, az);
            }

            printf("  %-8d %-14.3f %-14.4f %-14.3f %-14.4f\n",
                   i + 20,
                   comp.roll * 180.0 / M_PI,
                   fabs(comp.roll - roll_true) * 180.0 / M_PI,
                   comp.pitch * 180.0 / M_PI,
                   fabs(comp.pitch - pitch_true) * 180.0 / M_PI);
        }
    }

    printf("\nFinal Complementary estimate:\n");
    printf("  Roll  = %.3f deg (error = %.4f deg)\n",
           comp.roll * 180.0 / M_PI,
           fabs(comp.roll - roll_true) * 180.0 / M_PI);
    printf("  Pitch = %.3f deg (error = %.4f deg)\n",
           comp.pitch * 180.0 / M_PI,
           fabs(comp.pitch - pitch_true) * 180.0 / M_PI);

    /* ============================================================
     * Part 3: Multi-Sensor Statistical Fusion
     * ============================================================ */

    printf("\n--- Part 3: Multi-Sensor Statistical Fusion ---\n\n");

    /* Simulate 5 temperature sensors with different noise characteristics */
    double temps[] = {25.1, 24.8, 25.3, 25.0, 27.5};  /* last one is outlier */
    double vars[]  = {0.04, 0.09, 0.01, 0.16, 0.09};   /* variances */

    printf("Sensor readings:\n");
    {
        int i;
        for (i = 0; i < 5; i++) {
            printf("  Sensor %d: %.1f degC (sigma = %.2f degC)\n",
                   i + 1, temps[i], sqrt(vars[i]));
        }
    }

    /* Simple average (corrupted by outlier) */
    double simple_avg = 0.0;
    {
        int i;
        for (i = 0; i < 5; i++) simple_avg += temps[i];
        simple_avg /= 5.0;
    }

    /* Inverse-variance weighted average */
    double fused_var;
    double weighted_avg = ss_fusion_weighted_average(temps, vars, 5, &fused_var);

    /* Median */
    double median_val = ss_fusion_median(temps, 5);

    /* Trimmed mean (trim 1 from each end) */
    double trimmed_val = ss_fusion_trimmed_mean(temps, 5, 1);

    /* Outlier rejection */
    uint32_t mask;
    double robust_val;
    int n_accepted = ss_fusion_outlier_rejection(temps, vars, 5, 2.5, &mask, &robust_val);

    printf("\nFusion results (true value ≈ 25.0 degC):\n");
    printf("  Simple average:            %.3f degC  [corrupted by outlier]\n", simple_avg);
    printf("  Weighted average:          %.3f degC  (sigma = %.4f degC)\n",
           weighted_avg, sqrt(fused_var));
    printf("  Median:                    %.3f degC  [robust to outlier]\n", median_val);
    printf("  Trimmed mean (1 each):     %.3f degC\n", trimmed_val);
    printf("  Outlier rejection:         %.3f degC  (%d sensors accepted)\n",
           robust_val, n_accepted);

    printf("\n============================================================\n");
    printf("  IMU Sensor Fusion Example Complete\n");
    printf("============================================================\n");
    return 0;
}
