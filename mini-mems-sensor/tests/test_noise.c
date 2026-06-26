/**
 * @file test_noise.c
 * @brief Tests for L3-L4: Thermo-mechanical noise, Allan variance,
 *        PSD estimation, noise propagation, Johnson/Shot noise
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include "mems_noise.h"

#define PASS(msg) printf("  PASS: %s\n", msg)

static void test_noise_process_names(void) {
    printf("test_noise_process_names:\n");
    assert(strcmp(mems_noise_process_name(MEMS_NOISE_WHITE),
           "White Noise (Angle/Velocity RW)") == 0);
    PASS("White noise");
    assert(strcmp(mems_noise_process_name(MEMS_NOISE_FLICKER_FM),
           "Flicker FM (Bias Instability)") == 0);
    PASS("Flicker FM");
    assert(strcmp(mems_noise_process_name(MEMS_NOISE_RATE_RAMP),
           "Rate Ramp / Drift") == 0);
    PASS("Rate ramp");
}

static void test_thermomechanical_nea(void) {
    printf("test_thermomechanical_nea:\n");
    /* Typical consumer accel: m=2e-8 kg, fn=500 Hz, Q=200, T=300 K
     * NEA = sqrt(4*kB*T*2π*500 / (2e-8*200))
     *     = sqrt(4*1.38e-23*300*3142 / 4e-6)
     *     = sqrt(5.2e-17 / 4e-6) = sqrt(1.3e-11) ≈ 3.6e-6 g/√Hz? */
    double nea = mems_thermomech_nea(2e-8, 500.0, 200.0, 300.0);
    assert(nea > 1e-8 && nea < 1e-4);
    PASS("NEA in μg/√Hz range");

    /* Higher Q → lower NEA */
    double nea_hq = mems_thermomech_nea(2e-8, 500.0, 400.0, 300.0);
    assert(nea_hq < nea);
    PASS("Higher Q → lower NEA");

    /* Higher temperature → higher NEA */
    double nea_ht = mems_thermomech_nea(2e-8, 500.0, 200.0, 400.0);
    assert(nea_ht > nea);
    PASS("Higher T → higher NEA");

    /* Invalid inputs */
    assert(mems_thermomech_nea(0.0, 500.0, 200.0, 300.0) == 0.0);
    PASS("Zero mass → 0");
    assert(mems_thermomech_nea(2e-8, 0.0, 200.0, 300.0) == 0.0);
    PASS("Zero frequency → 0");
}

static void test_thermomechanical_ner(void) {
    printf("test_thermomechanical_ner:\n");
    double ner = mems_thermomech_ner(1e-3, 0.01, 1.0);
    /* NER = 1e-3 / 0.01 = 0.1 rad/s/√Hz */
    assert(fabs(ner - 0.1) < 1e-6);
    PASS("NER = NEA/v_drive");

    assert(isinf(mems_thermomech_ner(1e-3, 0.0, 1.0)));
    PASS("Zero v_drive → INF");
}

static void test_johnson_noise(void) {
    printf("test_johnson_noise:\n");
    /* R=10kΩ, T=300K → v_n_density = sqrt(4*1.38e-23*300*10000)
     * = sqrt(1.656e-16) ≈ 12.87 nV/√Hz */
    double vn = mems_johnson_noise_density(10000.0, 300.0);
    assert(fabs(vn - 12.87e-9) < 0.5e-9);
    PASS("v_n ≈ 12.87 nV/√Hz at 10kΩ");

    /* RMS in 1 kHz bandwidth */
    double vrms = mems_johnson_noise_rms(10000.0, 300.0, 1000.0);
    assert(fabs(vrms - 12.87e-9 * sqrt(1000.0)) < 1e-9);
    PASS("V_rms = v_n_density·√Δf");
}

static void test_shot_noise(void) {
    printf("test_shot_noise:\n");
    /* I_DC=1μA → i_n = sqrt(2*1.602e-19*1e-6) = sqrt(3.204e-25) ≈ 5.66e-13 A/√Hz */
    double in = mems_shot_noise_density(1e-6);
    assert(fabs(in - 5.66e-13) < 1e-14);
    PASS("i_n ≈ 0.566 pA/√Hz at 1μA");

    assert(mems_shot_noise_density(0.0) == 0.0);
    PASS("Zero current → 0");
}

static void test_allan_variance_basic(void) {
    printf("test_allan_variance_basic:\n");
    /* Constant data: Allan variance should be ~0 */
    double data[100];
    for (int i = 0; i < 100; i++) data[i] = 1.0;
    double avar = mems_allan_variance(data, 100, 1.0, 1.0);
    assert(fabs(avar) < 1e-12);
    PASS("Constant signal → AVAR ≈ 0");

    /* Not enough data */
    assert(mems_allan_variance(data, 3, 1.0, 1.0) < 0.0);
    PASS("n<4 → error");
}

static void test_allan_analysis(void) {
    printf("test_allan_analysis:\n");
    /* Generate white noise and compute Allan variance */
    double data[1024];
    uint64_t seed = 12345;
    for (int i = 0; i < 1024; i++) data[i] = 0.0;
    mems_generate_white_noise(data, 1024, 0.01, &seed);

    MemsAllanResult result = mems_allan_analysis(data, 1024, 100.0, 20);
    assert(result.num_points == 20);
    PASS("20-point Allan analysis");
    assert(result.points != NULL);
    PASS("Points allocated");

    /* White noise: ARW should be non-zero */
    assert(result.angle_rw > 0.0);
    PASS("ARW from white noise > 0");
    assert(result.bias_instability > 0.0);
    PASS("Bias instability identified");

    free(result.points);
}

static void test_rand64(void) {
    printf("test_rand64:\n");
    uint64_t seed = 42;
    double r = mems_rand64(&seed);
    assert(r >= 0.0 && r < 1.0);
    PASS("rand64 ∈ [0,1)");
}

static void test_white_noise_generation(void) {
    printf("test_white_noise_generation:\n");
    double output[1000];
    uint64_t seed = 99;

    for (int i = 0; i < 1000; i++) output[i] = 0.0;
    mems_generate_white_noise(output, 1000, 1.0, &seed);

    /* Check statistical properties: mean ≈ 0, std ≈ 1 */
    double sum = 0.0, sum_sq = 0.0;
    for (int i = 0; i < 1000; i++) {
        sum += output[i];
        sum_sq += output[i] * output[i];
    }
    double mean = sum / 1000.0;
    double std = sqrt(sum_sq / 1000.0 - mean * mean);
    assert(fabs(mean) < 0.1);
    PASS("Mean ≈ 0");
    assert(fabs(std - 1.0) < 0.2);
    PASS("Std ≈ 1");
}

static void test_flicker_noise_generation(void) {
    printf("test_flicker_noise_generation:\n");
    double output[500];
    uint64_t seed = 777;

    mems_generate_flicker_noise(output, 500, 1.0, 1.0, &seed);

    /* Ensure output is non-zero */
    int non_zero = 0;
    for (int i = 0; i < 500; i++)
        if (fabs(output[i]) > 1e-10) non_zero++;
    assert(non_zero > 0);
    PASS("Flicker noise has non-zero output");
}

static void test_noise_propagation(void) {
    printf("test_noise_propagation:\n");
    /* White noise PSD=1e-6 through system with fn=100Hz, ζ=0.1
     * Output RMS = sqrt(1e-6 * 2π*100 / (4*0.1)) = sqrt(1e-6 * 1570.8) ≈ 0.0396 */
    double rms = mems_noise_propagate_2nd_order(1e-6, 100.0, 0.1);
    assert(rms > 1e-3 && rms < 1e-1);
    PASS("Output RMS computed from NEB");

    /* Higher ζ → lower output noise */
    double rms_hz = mems_noise_propagate_2nd_order(1e-6, 100.0, 0.2);
    assert(rms_hz < rms);
    PASS("Higher ζ → lower output RMS");

    /* Invalid input */
    assert(mems_noise_propagate_2nd_order(1e-6, 0.0, 0.1) == 0.0);
    PASS("Zero ω_n → 0");
}

static void test_psd_welch(void) {
    printf("test_psd_welch:\n");
    /* Generate tone + noise */
    double data[512];
    uint64_t seed = 42;
    for (int i = 0; i < 512; i++) {
        data[i] = sin(2.0 * M_PI * 10.0 * (double)i / 100.0) +
                  0.1 * (mems_rand64(&seed) - 0.5);
    }

    MemsPsdEstimate psd = mems_welch_psd(data, 512, 100.0, 128, 0.5);
    assert(psd.num_bins > 0);
    PASS("PSD bins > 0");
    assert(psd.freq != NULL && psd.psd != NULL);
    PASS("Freq and PSD arrays allocated");
    assert(psd.total_power > 0.0);
    PASS("Total power > 0");

    /* RMS in band */
    double rms = mems_psd_rms(&psd, 5.0, 15.0);
    assert(rms > 0.0);
    PASS("RMS in [5,15] Hz band > 0");

    mems_psd_free(&psd);
    assert(psd.freq == NULL);
    PASS("PSD memory freed");
}

int main(void) {
    printf("=== test_noise ===\n\n");
    test_noise_process_names();
    test_thermomechanical_nea();
    test_thermomechanical_ner();
    test_johnson_noise();
    test_shot_noise();
    test_allan_variance_basic();
    test_allan_analysis();
    test_rand64();
    test_white_noise_generation();
    test_flicker_noise_generation();
    test_noise_propagation();
    test_psd_welch();
    printf("\n=== All noise tests PASSED ===\n");
    return 0;
}
