/**
 * @file test_calibration.c
 * @brief Tests for L5: 6-position accel calibration, gyro rate table,
 *        magnetometer ellipsoid fit, temperature compensation
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include "mems_calibration.h"

#define PASS(msg) printf("  PASS: %s\n", msg)
#define G 9.80665

static void test_accel_6pos_calib(void) {
    printf("test_accel_6pos_calib:\n");
    MemsSixPosData data;
    data.num_positions = 6;

    /* Expected gravity vectors (±g along each axis) */
    data.expected[0] = (MemsVector3){ G, 0.0, 0.0};
    data.expected[1] = (MemsVector3){-G, 0.0, 0.0};
    data.expected[2] = (MemsVector3){0.0,  G, 0.0};
    data.expected[3] = (MemsVector3){0.0, -G, 0.0};
    data.expected[4] = (MemsVector3){0.0, 0.0,  G};
    data.expected[5] = (MemsVector3){0.0, 0.0, -G};

    /* Simulated measurements: scale factor 1.02, bias 0.1g on x */
    double bias_x = 0.1 * G;
    double sf = 1.02;
    for (int i = 0; i < 6; i++) {
        data.measured[i].x = data.expected[i].x * sf + bias_x;
        data.measured[i].y = data.expected[i].y * sf;
        data.measured[i].z = data.expected[i].z * sf;
        data.temperature[i] = 25.0;
    }

    MemsAccelCalib calib;
    double rms = mems_calib_accel_6pos(&data, &calib);
    assert(rms >= 0.0);
    PASS("6-pos calibration converged");

    /* Bias should be approximately bias_x */
    assert(fabs(calib.bias.x - bias_x) < G * 0.02);
    PASS("Bias x identified correctly");

    /* Null pointer safety */
    assert(mems_calib_accel_6pos(NULL, &calib) < 0.0);
    PASS("NULL data handled");
}

static void test_apply_accel_calib(void) {
    printf("test_apply_accel_calib:\n");
    MemsAccelCalib calib = {0};
    calib.bias.x = 1.0;
    calib.bias.y = 2.0;
    calib.bias.z = 3.0;
    calib.scale_mis.m[0][0] = 2.0;
    calib.scale_mis.m[1][1] = 2.0;
    calib.scale_mis.m[2][2] = 2.0;

    MemsVector3 raw = {4.0, 5.0, 6.0};
    MemsVector3 corrected = mems_apply_accel_calib(raw, &calib);
    /* After bias subtraction: (3,3,3), scale ×2 → (6,6,6) */
    assert(fabs(corrected.x - 6.0) < 1e-10);
    assert(fabs(corrected.y - 6.0) < 1e-10);
    assert(fabs(corrected.z - 6.0) < 1e-10);
    PASS("Accel calibration applied correctly");
}

static void test_gyro_rate_calib(void) {
    printf("test_gyro_rate_calib:\n");
    /* Generate rate table data with known SF and bias */
    const size_t n = 5;
    double rates[n][3];
    double measured[n][3];
    double accel[n][3];

    for (size_t i = 0; i < n; i++) {
        rates[i][0] = (double)((int)i - 2) * 0.1;  /* -0.2, -0.1, 0, 0.1, 0.2 */
        rates[i][1] = 0.0;
        rates[i][2] = 0.0;
        /* Measurement: SF=1.1, bias=0.02 */
        measured[i][0] = 1.1 * rates[i][0] + 0.02;
        measured[i][1] = 0.0;
        measured[i][2] = 0.0;
        accel[i][0] = accel[i][1] = 0.0;
        accel[i][2] = G;
    }

    MemsGyroCalib calib;
    double rms = mems_calib_gyro_rate(rates, measured, accel, n, &calib);
    assert(rms >= 0.0);
    PASS("Gyro rate calibration converged");

    /* Verify SF and bias approximately correct */
    assert(fabs(calib.scale_factor.x - 1.1) < 0.1);
    PASS("X scale factor ≈ 1.1");
    assert(fabs(calib.bias.x - 0.02) < 0.05);
    PASS("X bias ≈ 0.02");
}

static void test_apply_gyro_calib(void) {
    printf("test_apply_gyro_calib:\n");
    MemsGyroCalib calib = {0};
    calib.bias.x = 0.1;
    calib.bias.y = 0.2;
    calib.bias.z = 0.3;
    calib.scale_mis.m[0][0] = 1.0;
    calib.scale_mis.m[1][1] = 1.0;
    calib.scale_mis.m[2][2] = 1.0;

    MemsVector3 raw = {1.1, 1.2, 1.3};
    MemsVector3 corrected = mems_apply_gyro_calib(raw, &calib);
    assert(fabs(corrected.x - 1.0) < 1e-10);
    PASS("Gyro cal applied");
}

static void test_mag_ellipsoid_calib(void) {
    printf("test_mag_ellipsoid_calib:\n");
    /* Generate points on a shifted sphere */
    const size_t n = 20;
    double measurements[n][3];
    double expected_mag = 50.0;  /* μT */
    double hx = 10.0, hy = -5.0, hz = 3.0;

    for (size_t i = 0; i < n; i++) {
        double theta = 2.0 * M_PI * (double)i / (double)n;
        double phi = M_PI * (double)i / (double)n;
        measurements[i][0] = hx + expected_mag * sin(phi) * cos(theta);
        measurements[i][1] = hy + expected_mag * sin(phi) * sin(theta);
        measurements[i][2] = hz + expected_mag * cos(phi);
    }

    MemsMagCalib calib;
    double rms = mems_calib_mag_ellipsoid(measurements, n, expected_mag, &calib);
    assert(rms >= 0.0);
    PASS("Mag ellipsoid fit converged");

    /* Hard-iron center should be approximately (hx, hy, hz) */
    assert(fabs(calib.hard_iron.x - hx) < expected_mag * 0.4);
    assert(fabs(calib.hard_iron.y - hy) < expected_mag * 0.4);
    assert(fabs(calib.hard_iron.z - hz) < expected_mag * 0.4);
    PASS("Hard-iron center identified");
}

static void test_apply_mag_calib(void) {
    printf("test_apply_mag_calib:\n");
    MemsMagCalib calib = {0};
    calib.hard_iron.x = 10.0;
    calib.hard_iron.y = -5.0;
    calib.hard_iron.z = 3.0;
    calib.soft_iron.m[0][0] = 1.0;
    calib.soft_iron.m[1][1] = 1.0;
    calib.soft_iron.m[2][2] = 1.0;

    MemsVector3 raw = {10.0, -5.0, 3.0};  /* at center → output (0,0,0) */
    MemsVector3 corrected = mems_apply_mag_calib(raw, &calib);
    assert(fabs(corrected.x) < 1e-10);
    assert(fabs(corrected.y) < 1e-10);
    assert(fabs(corrected.z) < 1e-10);
    PASS("Mag calibration: center → origin");
}

static void test_temp_compensation(void) {
    printf("test_temp_compensation:\n");
    /* Generate temperature-dependent errors:
     * error(T) = a1·ΔT + a2·ΔT² with a1=0.01, a2=-0.0005
     * T_ref = 25°C */
    const size_t n = 10;
    double temps[n], errors[n];

    for (size_t i = 0; i < n; i++) {
        temps[i] = -20.0 + (double)i * 10.0;
        double dt = temps[i] - 25.0;
        errors[i] = 0.01 * dt - 0.0005 * dt * dt;
    }

    MemsTempComp comp;
    double r2 = mems_fit_temp_comp(temps, errors, n, 2, &comp);
    assert(r2 > 0.9);
    PASS("Temp comp R² > 0.9");
    assert(fabs(comp.ref_temp - 25.0) < 10.0);
    PASS("Reference temperature identified");

    /* Apply compensation: should reduce error */
    double raw = 100.0;
    double dt = 30.0;  /* 30°C deviation from ref */
    double comp_val = mems_apply_temp_comp(raw, 25.0 + dt, &comp);
    /* Compensated value should not be NaN */
    assert(!isnan(comp_val));
    PASS("Temp compensation applied without NaN");
}

static void test_alignment_matrix(void) {
    printf("test_alignment_matrix:\n");
    /* Identity alignment test */
    const size_t n = 4;
    double measured[4][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 1.0, 0.0}
    };
    double reference[4][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 1.0, 0.0}
    };

    MemsMatrix33 S;
    double cond = mems_alignment_matrix(measured, reference, n, &S);
    assert(cond > 0.0);
    PASS("Alignment matrix computed");

    /* For identity, S should be near identity */
    assert(fabs(S.m[0][0] - 1.0) < 0.1);
    PASS("S[0][0] ≈ 1");
    assert(fabs(S.m[1][1] - 1.0) < 0.1);
    PASS("S[1][1] ≈ 1");
}

static void test_apply_matrix33(void) {
    printf("test_apply_matrix33:\n");
    MemsMatrix33 M = {{{2.0, 0.0, 0.0},
                        {0.0, 3.0, 0.0},
                        {0.0, 0.0, 4.0}}};
    MemsVector3 v = {1.0, 1.0, 1.0};
    MemsVector3 result = mems_apply_matrix33(v, &M);
    assert(fabs(result.x - 2.0) < 1e-10);
    assert(fabs(result.y - 3.0) < 1e-10);
    assert(fabs(result.z - 4.0) < 1e-10);
    PASS("Matrix-vector multiplication");
}

int main(void) {
    printf("=== test_calibration ===\n\n");
    test_accel_6pos_calib();
    test_apply_accel_calib();
    test_gyro_rate_calib();
    test_apply_gyro_calib();
    test_mag_ellipsoid_calib();
    test_apply_mag_calib();
    test_temp_compensation();
    test_alignment_matrix();
    test_apply_matrix33();
    printf("\n=== All calibration tests PASSED ===\n");
    return 0;
}
