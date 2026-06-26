/**
 * @file example_imu_calib.c
 * @brief L6 Example: Full IMU calibration pipeline (accelerometer 6-position
 *        + gyroscope rate table + magnetometer ellipsoid fitting)
 *
 * This demonstrates a complete MEMS-IMU calibration workflow as used
 * in drone autopilots (ArduPilot, PX4), smartphone sensor calibration,
 * and automotive sensor characterization.
 *
 * Key real-world references:
 *   - DJI Phantom 4 IMU calibration
 *   - Apple iPhone motion calibration
 *   - Bosch BMI270 factory calibration
 *   - Tesla Autopilot IMU calibration
 */

#include <stdio.h>
#include <math.h>
#include "mems_calibration.h"
#include "mems_basics.h"

#define G 9.80665

int main(void) {
    printf("=== MEMS IMU Calibration Example ===\n\n");

    /* ─── Step 1: Accelerometer 6-position calibration ───────────────── */
    printf("1. Accelerometer 6-Position Calibration\n");
    printf("   (Simulating Bosch BMI270 at ±2g range)\n\n");

    MemsSixPosData accel_data;
    accel_data.num_positions = 6;

    /* Set expected gravity references */
    double g[6][3] = {
        { G, 0, 0},   /* +X up */
        {-G, 0, 0},   /* -X down */
        {0,  G, 0},   /* +Y up */
        {0, -G, 0},   /* -Y down */
        {0, 0,  G},   /* +Z up */
        {0, 0, -G}    /* -Z down */
    };

    /* Simulate measurements with imperfections:
     *   scale factor error: +2% on each axis
     *   bias: 30mg on X, -20mg on Y, 15mg on Z
     *   noise: ~5mg random */
    double sf_err = 1.02;  /* +2% scale error */
    double b_x = 0.03 * G, b_y = -0.02 * G, b_z = 0.015 * G;

    printf("   True parameters:\n");
    printf("     Scale factor error: +2%%\n");
    printf("     Bias X: +30mg, Bias Y: -20mg, Bias Z: +15mg\n\n");

    for (int i = 0; i < 6; i++) {
        double mx = g[i][0] * sf_err + b_x;
        double my = g[i][1] * sf_err + b_y;
        double mz = g[i][2] * sf_err + b_z;

        accel_data.expected[i].x = g[i][0];
        accel_data.expected[i].y = g[i][1];
        accel_data.expected[i].z = g[i][2];

        accel_data.measured[i].x = mx;
        accel_data.measured[i].y = my;
        accel_data.measured[i].z = mz;

        accel_data.temperature[i] = 25.0;
    }

    MemsAccelCalib accel_calib;
    double accel_rms = mems_calib_accel_6pos(&accel_data, &accel_calib);

    printf("   Calibration results:\n");
    printf("     Bias X: %.4f m/s^2 (%.2f mg)\n",
           accel_calib.bias.x, mems_ms2_to_g(accel_calib.bias.x) * 1000.0);
    printf("     Bias Y: %.4f m/s^2 (%.2f mg)\n",
           accel_calib.bias.y, mems_ms2_to_g(accel_calib.bias.y) * 1000.0);
    printf("     Bias Z: %.4f m/s^2 (%.2f mg)\n",
           accel_calib.bias.z, mems_ms2_to_g(accel_calib.bias.z) * 1000.0);
    printf("     Scale X: %.4f, Y: %.4f, Z: %.4f\n",
           accel_calib.scale_factor.x, accel_calib.scale_factor.y,
           accel_calib.scale_factor.z);
    printf("     RMS error: %.4f m/s^2\n\n", accel_rms);

    /* ─── Step 2: Gyroscope rate table calibration ──────────────────── */
    printf("2. Gyroscope Rate Table Calibration\n");
    printf("   (Simulating ST LSM6DSO on rate table ±200°/s)\n\n");

    const size_t n_rates = 7;
    double rates[7][3];
    double gyro_meas[7][3];
    double accels[7][3];

    /* Known rotation rates (rad/s) */
    double rate_vals[7] = {-3.49, -2.0, -1.0, 0.0, 1.0, 2.0, 3.49};
    double gyro_sf = 0.98;  /* -2% scale error */
    double gyro_bias = 0.05; /* 0.05 rad/s bias ≈ 2.86°/s */

    for (size_t i = 0; i < n_rates; i++) {
        rates[i][0] = rate_vals[i];
        rates[i][1] = 0.0;
        rates[i][2] = 0.0;
        gyro_meas[i][0] = gyro_sf * rate_vals[i] + gyro_bias;
        gyro_meas[i][1] = 0.0;
        gyro_meas[i][2] = 0.0;
        accels[i][0] = 0.0; accels[i][1] = 0.0; accels[i][2] = G;
    }

    MemsGyroCalib gyro_calib;
    double gyro_rms = mems_calib_gyro_rate(rates, gyro_meas, accels,
                                            n_rates, &gyro_calib);

    printf("   Calibration results:\n");
    printf("     Bias X: %.4f rad/s (%.2f deg/s)\n",
           gyro_calib.bias.x, mems_rad_to_deg(gyro_calib.bias.x));
    printf("     Scale X: %.4f\n", gyro_calib.scale_factor.x);
    printf("     RMS error: %.4f rad/s\n\n", gyro_rms);

    /* ─── Step 3: Magnetometer calibration ─────────────────────────── */
    printf("3. Magnetometer Ellipsoid Calibration\n");
    printf("   (Simulating ST LIS3MDL with hard-iron + soft-iron distortion)\n\n");

    const size_t n_mag = 20;
    double mag_meas[20][3];
    double mag_expected = 50.0;  /* Earth's field ≈ 50 μT */
    double hard_x = 12.0, hard_y = -8.0, hard_z = 5.0;  /* hard iron */

    for (size_t i = 0; i < n_mag; i++) {
        double theta = 2.0 * M_PI * (double)i / (double)n_mag;
        double phi = M_PI * (double)i / (double)n_mag;
        mag_meas[i][0] = hard_x + mag_expected * sin(phi) * cos(theta);
        mag_meas[i][1] = hard_y + mag_expected * sin(phi) * sin(theta);
        mag_meas[i][2] = hard_z + mag_expected * cos(phi);
    }

    MemsMagCalib mag_calib;
    double mag_rms = mems_calib_mag_ellipsoid(mag_meas, n_mag,
                                               mag_expected, &mag_calib);

    printf("   Calibration results:\n");
    printf("     Hard Iron X: %.2f, Y: %.2f, Z: %.2f μT\n",
           mag_calib.hard_iron.x, mag_calib.hard_iron.y, mag_calib.hard_iron.z);
    printf("     Scale Factors: %.4f, %.4f, %.4f\n",
           mag_calib.scale_factor.x, mag_calib.scale_factor.y,
           mag_calib.scale_factor.z);
    printf("     RMS error: %.4f μT\n\n", mag_rms);

    /* ─── Step 4: Verify calibration quality ──────────────────────── */
    printf("4. Calibration Verification\n\n");

    /* Apply calibrations to a known test point */
    MemsVector3 test_accel = {0.0, 0.0, G * sf_err + b_z};
    MemsVector3 corrected_accel = mems_apply_accel_calib(test_accel,
                                                          &accel_calib);
    printf("   Accel test: (0, 0, %.2f) → (%.3f, %.3f, %.3f) [expected: (0,0,%.2f)]\n",
           G * sf_err + b_z,
           corrected_accel.x, corrected_accel.y, corrected_accel.z, G);

    /* Test gyro at 1 rad/s */
    MemsVector3 test_gyro = {1.0 * gyro_sf + gyro_bias, 0.0, 0.0};
    MemsVector3 corrected_gyro = mems_apply_gyro_calib(test_gyro, &gyro_calib);
    printf("   Gyro test:  (%.3f, 0, 0) → (%.3f, %.3f, %.3f) [expected: (1,0,0)]\n\n",
           test_gyro.x, corrected_gyro.x, corrected_gyro.y, corrected_gyro.z);

    printf("=== IMU Calibration Complete ===");

    return 0;
}
