/**
 * @file test_basics.c
 * @brief Tests for L1-L2: MEMS sensor types, vector math, quaternion, spec init
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include "mems_basics.h"

#define PASS(msg) printf("  PASS: %s\n", msg)

static void test_sensor_type_names(void) {
    printf("test_sensor_type_names:\n");
    assert(strcmp(mems_sensor_type_name(MEMS_TYPE_ACCELEROMETER), "Accelerometer") == 0);
    PASS("Accelerometer");
    assert(strcmp(mems_sensor_type_name(MEMS_TYPE_GYROSCOPE), "Gyroscope") == 0);
    PASS("Gyroscope");
    assert(strcmp(mems_sensor_type_name(MEMS_TYPE_MAGNETOMETER), "Magnetometer") == 0);
    PASS("Magnetometer");
    assert(strcmp(mems_sensor_type_name(MEMS_TYPE_MICROPHONE), "Microphone (MEMS)") == 0);
    PASS("Microphone");
    assert(strcmp(mems_sensor_type_name(MEMS_TYPE_RESONATOR), "Resonator/Oscillator") == 0);
    PASS("Resonator");
    assert(strcmp(mems_sensor_type_name(999), "Unknown") == 0);
    PASS("Invalid type -> Unknown");
}

static void test_transduce_type_names(void) {
    printf("test_transduce_type_names:\n");
    assert(strcmp(mems_transduce_type_name(MEMS_TRANSDUCE_CAPACITIVE), "Capacitive") == 0);
    PASS("Capacitive");
    assert(strcmp(mems_transduce_type_name(MEMS_TRANSDUCE_PIEZORESISTIVE), "Piezoresistive") == 0);
    PASS("Piezoresistive");
    assert(strcmp(mems_transduce_type_name(MEMS_TRANSDUCE_PIEZOELECTRIC), "Piezoelectric") == 0);
    PASS("Piezoelectric");
    assert(strcmp(mems_transduce_type_name(MEMS_TRANSDUCE_TUNNELING), "Electron Tunneling") == 0);
    PASS("Tunneling");
    assert(strcmp(mems_transduce_type_name(999), "Unknown") == 0);
    PASS("Invalid -> Unknown");
}

static void test_fab_type_names(void) {
    printf("test_fab_type_names:\n");
    assert(strcmp(mems_fab_type_name(MEMS_FAB_BULK_MICROMACHINED), "Bulk Micromachined") == 0);
    PASS("Bulk Micromachined");
    assert(strcmp(mems_fab_type_name(MEMS_FAB_SURFACE_MICROMACH), "Surface Micromachined") == 0);
    PASS("Surface Micromachined");
    assert(strcmp(mems_fab_type_name(MEMS_FAB_CMOS_MEMS), "CMOS-MEMS") == 0);
    PASS("CMOS-MEMS");
}

static void test_resonant_freq(void) {
    printf("test_resonant_freq:\n");
    /* Si proof mass: m = 1e-9 kg, k = 100 N/m → f_n = sqrt(100/1e-9)/(2π) ≈ 50329 Hz */
    double f1 = mems_resonant_freq(1e-9, 100.0);
    assert(f1 > 40000.0 && f1 < 60000.0);
    PASS("f_n ≈ 50 kHz for m=1ng, k=100 N/m");

    /* Handling of invalid inputs */
    assert(mems_resonant_freq(0.0, 100.0) == 0.0);
    PASS("Zero mass → 0");
    assert(mems_resonant_freq(1e-9, -1.0) == 0.0);
    PASS("Negative k → 0");
}

static void test_quality_factor(void) {
    printf("test_quality_factor:\n");
    /* m = 1e-9, k = 100, c = 1e-8
     * Q = sqrt(100*1e-9)/1e-8 = sqrt(1e-7)/1e-8 ≈ 3.162e-4/1e-8 = 31623 */
    double Q = mems_quality_factor(1e-9, 100.0, 1e-8);
    assert(Q > 30000.0 && Q < 33000.0);
    PASS("Q = sqrt(k·m)/c ≈ 31623");

    /* Invalid input */
    assert(mems_quality_factor(0.0, 100.0, 0.01) == 0.0);
    PASS("Zero mass → 0");
    assert(mems_quality_factor(1e-9, 100.0, 0.0) == 0.0);
    PASS("Zero damping → 0");
}

static void test_damping_ratio(void) {
    printf("test_damping_ratio:\n");
    /* m = 1e-9, k = 100, c = 1e-8 → ζ = 1e-8/(2*sqrt(100*1e-9)) = 1e-8/(2*3.16e-4) ≈ 1.58e-5 */
    double zeta = mems_damping_ratio(1e-9, 100.0, 1e-8);
    assert(zeta > 0.0 && zeta < 1.0);
    PASS("ζ = c/(2·sqrt(k·m)) < 1 (underdamped)");

    assert(mems_damping_ratio(1e-9, 100.0, 0.0) == 0.0);
    PASS("Zero damping → ζ=0 (undamped)");
}

static void test_capacitive_sensitivity(void) {
    printf("test_capacitive_sensitivity:\n");
    /* ε = 8.854e-12, A = 1e-8 m², gap = 2e-6 m, displacement = 1e-9 m
     * ΔC ≈ 8.854e-12 * 1e-8 * 1e-9 / (2e-6)² = 8.854e-29 / 4e-12 ≈ 2.2e-17 F */
    double dc = mems_capacitive_sensitivity(8.854e-12, 1e-8, 2e-6, 1e-9);
    assert(fabs(dc - 2.21e-17) < 1e-17);
    PASS("ΔC ≈ 2.2e-17 F");

    assert(mems_capacitive_sensitivity(8.854e-12, 1e-8, 0.0, 1e-9) == 0.0);
    PASS("Zero gap → 0");
}

static void test_vector3_ops(void) {
    printf("test_vector3_ops:\n");
    MemsVector3 a = {3.0, 4.0, 0.0};
    double mag = mems_vector3_mag(a);
    assert(fabs(mag - 5.0) < 1e-10);
    PASS("|v| = 5");

    MemsVector3 n = mems_vector3_normalize(a);
    assert(fabs(n.x - 0.6) < 1e-10 && fabs(n.y - 0.8) < 1e-10);
    PASS("Normalize → (0.6, 0.8, 0)");

    MemsVector3 b = {1.0, 0.0, 0.0};
    MemsVector3 c = mems_vector3_cross(a, b);
    /* (3,4,0)×(1,0,0) = (0,0,-4) */
    assert(fabs(c.x) < 1e-10 && fabs(c.y) < 1e-10 && fabs(c.z + 4.0) < 1e-10);
    PASS("Cross product");

    double dot = mems_vector3_dot(a, b);
    assert(fabs(dot - 3.0) < 1e-10);
    PASS("Dot product = 3");
}

static void test_unit_conversions(void) {
    printf("test_unit_conversions:\n");
    double rad = mems_deg_to_rad(180.0);
    assert(fabs(rad - M_PI) < 1e-10);
    PASS("180° → π rad");

    double deg = mems_rad_to_deg(M_PI);
    assert(fabs(deg - 180.0) < 1e-10);
    PASS("π rad → 180°");

    double g_val = mems_ms2_to_g(9.80665);
    assert(fabs(g_val - 1.0) < 1e-10);
    PASS("9.80665 m/s² → 1g");

    double ms2 = mems_g_to_ms2(1.0);
    assert(fabs(ms2 - 9.80665) < 1e-10);
    PASS("1g → 9.80665 m/s²");

    double radps = mems_dps_to_radps(57.2957795);
    assert(fabs(radps - 1.0) < 1e-6);
    PASS("57.3°/s → ≈1 rad/s");
}

static void test_spec_init(void) {
    printf("test_spec_init:\n");
    MemsSpec spec;
    mems_spec_typical_accel(&spec);
    assert(spec.type == MEMS_TYPE_ACCELEROMETER);
    PASS("Accel spec: correct type");
    assert(spec.fsr > 0.0);
    PASS("Accel spec: FSR > 0");
    assert(spec.noise_density > 0.0);
    PASS("Accel spec: noise density > 0");

    mems_spec_typical_gyro(&spec);
    assert(spec.type == MEMS_TYPE_GYROSCOPE);
    PASS("Gyro spec: correct type");
    assert(spec.g_sensitivity > 0.0);
    PASS("Gyro spec: g-sensitivity > 0");

    mems_spec_typical_mag(&spec);
    assert(spec.type == MEMS_TYPE_MAGNETOMETER);
    PASS("Mag spec: correct type");
    assert(spec.fsr == 1600.0);
    PASS("Mag spec: FSR = 1600 μT");
}

int main(void) {
    printf("=== test_basics ===\n\n");
    test_sensor_type_names();
    test_transduce_type_names();
    test_fab_type_names();
    test_resonant_freq();
    test_quality_factor();
    test_damping_ratio();
    test_capacitive_sensitivity();
    test_vector3_ops();
    test_unit_conversions();
    test_spec_init();
    printf("\n=== All basics tests PASSED ===\n");
    return 0;
}
