/**
 * @file test_sensor_core.c
 * @brief Tests for core sensor definitions, transfer functions, noise
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "../include/sensor_defs.h"
#include "../include/sensor_transfer.h"
#include "../include/sensor_noise.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  TEST %s ... ", #name); \
} while(0)

#define PASS() do { \
    tests_passed++; \
    printf("PASS\n"); \
} while(0)

#define FAIL(msg) do { \
    printf("FAIL: %s\n", msg); \
} while(0)

#define ASSERT_TRUE(cond, msg) do { \
    if (!(cond)) { FAIL(msg); return; } \
} while(0)

#define ASSERT_FLOAT_EQ(a, b, tol, msg) do { \
    if (fabs((a) - (b)) > (tol)) { \
        printf("FAIL: %s (got %.6g, expected %.6g)\n", msg, (a), (b)); \
        return; \
    } \
} while(0)

/* ──── L1: Sensor Specs Validation ──── */

static void test_specs_validation(void) {
    TEST(specs_validation);
    sensor_specs_t s = {0};
    s.sensitivity = 1.0;
    s.full_scale_input = 100.0;
    s.full_scale_output = 5.0;
    s.bandwidth = 1000.0;
    s.temp_range_min = -40.0;
    s.temp_range_max = 85.0;
    ASSERT_TRUE(sensor_specs_validate(&s) == 0, "valid specs should return 0");
    PASS();
}

static void test_specs_invalid_zero_sensitivity(void) {
    TEST(specs_zero_sensitivity);
    sensor_specs_t s = {0};
    s.sensitivity = 0.0;
    s.full_scale_input = 100.0;
    s.full_scale_output = 5.0;
    s.bandwidth = 1000.0;
    s.temp_range_min = -40.0;
    s.temp_range_max = 85.0;
    ASSERT_TRUE(sensor_specs_validate(&s) == 0x01, "zero sensitivity should set bit 0");
    PASS();
}

static void test_total_error_percent(void) {
    TEST(total_error_percent);
    sensor_specs_t s = {0};
    s.nonlinearity = 0.1;   /* 0.1% FS */
    s.hysteresis = 0.05;
    s.repeatability = 0.02;
    s.offset_error = 0.01;
    s.full_scale_output = 5.0;
    s.gain_error = 0.1;
    double err = sensor_total_error_percent(&s);
    /* Expected: sqrt(0.1² + 0.05² + 0.02² + (0.01/5*100)² + 0.1²)
     *         = sqrt(0.01 + 0.0025 + 0.0004 + 0.04 + 0.01)
     *         = sqrt(0.0629) ≈ 0.2508% */
    ASSERT_TRUE(err > 0.24 && err < 0.26, "total error should be ~0.25%");
    PASS();
}

/* ──── L3: Dynamic Error ──── */

static void test_dynamic_error_percent(void) {
    TEST(dynamic_error_percent);
    /* At f = f_c (bandwidth), magnitude = 1/√2 ≈ 0.707, error ≈ 29.3% */
    double err = sensor_dynamic_error_percent(100.0, 100.0);
    ASSERT_FLOAT_EQ(err, 100.0 * (1.0 - 1.0/sqrt(2.0)), 0.1, "at f=fc, error ~29.3%");
    /* At f = 0, error should be 0 */
    err = sensor_dynamic_error_percent(0.0, 100.0);
    ASSERT_FLOAT_EQ(err, 0.0, 0.001, "at DC, error=0");
    PASS();
}

/* ──── L3: Uncertainty Budget ──── */

static void test_uncertainty_budget(void) {
    TEST(uncertainty_budget);
    sensor_uncertainty_budget_t b = {0};
    b.type_a_repeatability = 0.01;   /* statistical */
    b.type_b_reference     = 0.02;   /* systematic */
    b.type_b_resolution    = 0.005;
    b.type_b_thermal       = 0.01;
    b.type_b_noise         = 0.002;
    b.degrees_of_freedom   = 9;       /* n=10 → ν=9 */
    ASSERT_TRUE(sensor_uncertainty_budget_compute(&b) == 0, "budget compute should succeed");
    ASSERT_TRUE(b.combined_standard > 0.02, "combined uncertainty exists");
    ASSERT_TRUE(b.expanded_uncertainty == 2.0 * b.combined_standard,
                "expanded should be k=2 * combined");
    PASS();
}

/* ──── L4: Johnson-Nyquist Noise ──── */

static void test_johnson_noise(void) {
    TEST(johnson_noise);
    /* At T=300K, R=1kΩ, Δf=1Hz: Vn ≈ 4.07 nV */
    double vn = johnson_noise_voltage(300.0, 1000.0, 1.0);
    double expected = sqrt(4.0 * 1.380649e-23 * 300.0 * 1000.0 * 1.0);
    ASSERT_FLOAT_EQ(vn, expected, 1e-15, "johnson noise voltage calculated correctly");

    /* Power available independent of resistance */
    double P1 = johnson_noise_power_available(300.0, 1000.0);
    double P2 = johnson_noise_power_available(300.0, 1000.0);  /* same param */
    ASSERT_FLOAT_EQ(P1, P2, 1e-30, "available noise power doesn't depend on R");
    PASS();
}

/* ──── L4: Shot Noise ──── */

static void test_shot_noise(void) {
    TEST(shot_noise);
    /* I_DC=1mA, Δf=1kHz: In ≈ 5.66e-10 A = 566 pA */
    double inoise = shot_noise_current(1e-3, 1000.0);
    double expected = sqrt(2.0 * 1.602176634e-19 * 1e-3 * 1000.0);
    ASSERT_FLOAT_EQ(inoise, expected, 1e-25, "shot noise calculated correctly");

    /* Shot/thermal ratio ≈ qI/(2kT):
     *   At I=1mA:  ratio ≈ 0.019 → NOT significant (thermal dominates)
     *   At I=100mA: ratio ≈ 1.93 → IS significant (shot dominates)
     * Threshold of 1.0 means shot power must exceed thermal power. */
    int sig_low = shot_noise_is_significant(1e-3, 50.0, 300.0, 1000.0, 1.0);
    ASSERT_TRUE(!sig_low, "shot noise NOT significant at 1mA (thermal dominates)");

    int sig_high = shot_noise_is_significant(0.1, 50.0, 300.0, 1000.0, 1.0);
    ASSERT_TRUE(sig_high, "shot noise IS significant at 100mA through 50Ω");
    PASS();
}

/* ──── L3: 1/f Noise ──── */

static void test_flicker_noise(void) {
    TEST(flicker_noise);
    flicker_noise_spec_t fn = {0};
    fn.K_f = 1e-12;     /* 1 pV²/Hz at 1 Hz */
    fn.alpha = 1.0;
    fn.f_low = 0.1;     /* 100 mHz */
    fn.f_high = 10.0;   /* 10 Hz */

    double vrms = flicker_noise_integrated_rms(&fn);
    /* For α=1: V² = K_f·ln(f_high/f_low) = 1e-12·ln(100) ≈ 4.6e-12 */
    double expected_var = 1e-12 * log(100.0);
    ASSERT_FLOAT_EQ(vrms, sqrt(expected_var), 1e-9, "1/f integrated noise correct");
    PASS();
}

/* ──── L3: Allan Variance ──── */

static void test_allan_variance(void) {
    TEST(allan_variance);
    /* Generate synthetic data: white noise + drift */
    size_t N = 1024;
    double *data = (double *)calloc(N, sizeof(double));
    assert(data);

    /* White noise with std=0.1, plus linear drift 0.001 per sample */
    for (size_t i = 0; i < N; i++) {
        /* Deterministic pseudo-random */
        double noise = sin((double)i * 7919.0) * 0.1;
        double drift = (double)i * 0.001;
        data[i] = noise + drift;
    }

    allan_variance_t av = {0};
    int ret = allan_variance_compute_overlapping(data, N, 0.01, &av, 10);
    ASSERT_TRUE(ret == 0, "allan variance computation should succeed");
    ASSERT_TRUE(av.num_clusters > 0, "should have some clusters");
    ASSERT_TRUE(av.min_allan_dev > 0.0, "should have non-zero allan deviation");

    allan_variance_free(&av);
    free(data);
    PASS();
}

/* ──── L2: Polynomial Transfer ──── */

static void test_polynomial_transfer(void) {
    TEST(polynomial_transfer);
    poly_transfer_t pt = {0};
    double coeff[] = {1.0, 2.0, 0.5};  /* y = 1 + 2x + 0.5x² */
    ASSERT_TRUE(poly_transfer_init(&pt, coeff, 2, -10.0, 10.0) == 0,
                "poly init should succeed");

    /* f(0) = 1.0 */
    ASSERT_FLOAT_EQ(poly_transfer_evaluate(&pt, 0.0), 1.0, 1e-10, "f(0)=1.0");
    /* f(2) = 1 + 4 + 2 = 7.0 */
    ASSERT_FLOAT_EQ(poly_transfer_evaluate(&pt, 2.0), 7.0, 1e-10, "f(2)=7.0");
    /* Sensitivity at x=0: f'(x) = 2 + x, f'(0) = 2 */
    ASSERT_FLOAT_EQ(poly_transfer_sensitivity(&pt, 0.0), 2.0, 1e-10, "S(0)=2.0");

    poly_transfer_free(&pt);
    PASS();
}

/* ──── L2: 1st Order Sensor ──── */

static void test_first_order_sensor(void) {
    TEST(first_order_sensor);
    sensor_1st_order_t s1;
    sensor_1st_order_init(&s1, 2.0, 0.5, 0.0, 0.0);

    /* Step response at t=tau: y = K·A·(1-1/e) = 2·1·0.6321 = 1.264 */
    double y_tau = sensor_1st_order_step_response(&s1, 1.0, 0.5);
    ASSERT_FLOAT_EQ(y_tau, 2.0 * (1.0 - 1.0/M_E), 0.01, "step response at t=tau");

    /* Bandwidth: f_c = 1/(2πτ) */
    double bw = sensor_1st_order_bandwidth(&s1);
    ASSERT_FLOAT_EQ(bw, 1.0/(2.0*M_PI*0.5), 0.001, "bandwidth = 1/(2πτ)");

    /* Rise time: t_r = τ·ln(9) */
    double tr = sensor_1st_order_rise_time(&s1);
    ASSERT_FLOAT_EQ(tr, 0.5 * log(9.0), 0.001, "rise time = τ·ln(9)");
    PASS();
}

/* ──── L2: 2nd Order Sensor ──── */

static void test_second_order_sensor(void) {
    TEST(second_order_sensor);
    sensor_2nd_order_t s2;
    sensor_2nd_order_init(&s2, 1.0, 10.0, 0.3, 0.0, 0.0);

    /* Overshoot: M_p = exp(-π·ζ/√(1-ζ²)) */
    double Mp = sensor_2nd_order_overshoot_percent(&s2);
    double Mp_expected = 100.0 * exp(-M_PI * 0.3 / sqrt(1.0 - 0.09));
    ASSERT_FLOAT_EQ(Mp, Mp_expected, 0.1, "overshoot percentage");

    /* Peak time: t_p = π/(ω_n·√(1-ζ²)) */
    double tp = sensor_2nd_order_peak_time(&s2);
    double tp_expected = M_PI / (10.0 * sqrt(1.0 - 0.09));
    ASSERT_FLOAT_EQ(tp, tp_expected, 0.001, "peak time");

    /* For ζ ≥ 1, no overshoot */
    sensor_2nd_order_t s2_crit;
    sensor_2nd_order_init(&s2_crit, 1.0, 10.0, 1.0, 0.0, 0.0);
    ASSERT_TRUE(sensor_2nd_order_overshoot_percent(&s2_crit) == 0.0,
                "critically damped has no overshoot");
    PASS();
}

/* ──── L2: Thermistor Model ──── */

static void test_thermistor_beta(void) {
    TEST(thermistor_beta);
    thermistor_model_t tm;
    /* Standard NTC: R25=10kΩ, β=3977K (25/85) */
    thermistor_model_init_beta(&tm, 10000.0, 298.15, 3977.0);

    /* At T=25°C (298.15K), R should be ~10kΩ */
    double R = thermistor_R_from_temp(&tm, 298.15);
    ASSERT_FLOAT_EQ(R, 10000.0, 10.0, "R(25°C) ≈ 10000Ω");

    /* Temperature from that resistance */
    double T = thermistor_temp_from_R(&tm, R);
    ASSERT_FLOAT_EQ(T, 298.15, 0.1, "T from R round-trip");

    /* At T=85°C (358.15K), R should be lower */
    double R_85 = thermistor_R_from_temp(&tm, 358.15);
    ASSERT_TRUE(R_85 < 10000.0, "NTC: R decreases with temperature");
    PASS();
}

/* ──── L2: RTD CVD Model ──── */

static void test_rtd_cvd(void) {
    TEST(rtd_cvd);
    rtd_cvd_model_t rtd;
    rtd_cvd_init_iec751(&rtd, 100.0);  /* Pt100 */

    /* R(0°C) = 100Ω */
    ASSERT_FLOAT_EQ(rtd_R_from_T(&rtd, 0.0), 100.0, 0.01, "R(0°C)=100Ω");

    /* R(100°C) ≈ 100·(1 + 0.0039083·100 - 5.775e-7·10000) ≈ 138.51 */
    double R100 = rtd_R_from_T(&rtd, 100.0);
    ASSERT_FLOAT_EQ(R100, 138.51, 0.1, "R(100°C) ≈ 138.51Ω");

    /* R(200°C) */
    double R200 = rtd_R_from_T(&rtd, 200.0);
    ASSERT_TRUE(R200 > 170.0, "R increases with temperature");
    ASSERT_TRUE(R200 < 180.0, "R at 200C is plausible");
    PASS();
}

/* ──── L3: Linearity Analysis ──── */

static void test_linearity_independent(void) {
    TEST(linearity_independent);
    /* Perfectly linear data */
    double x[] = {0, 1, 2, 3, 4, 5};
    double y[] = {0, 2, 4, 6, 8, 10};  /* y = 2x */
    double slope, intercept, max_err;

    sensor_linearity_independent(x, y, 6, &slope, &intercept, &max_err);
    ASSERT_FLOAT_EQ(slope, 2.0, 0.001, "slope = 2");
    ASSERT_FLOAT_EQ(intercept, 0.0, 0.001, "intercept = 0");
    ASSERT_FLOAT_EQ(max_err, 0.0, 0.001, "perfect linear: max error = 0%");
    PASS();
}

int main(void) {
    printf("=== Sensor Core Tests ===\n\n");

    test_specs_validation();
    test_specs_invalid_zero_sensitivity();
    test_total_error_percent();
    test_dynamic_error_percent();
    test_uncertainty_budget();
    test_johnson_noise();
    test_shot_noise();
    test_flicker_noise();
    test_allan_variance();
    test_polynomial_transfer();
    test_first_order_sensor();
    test_second_order_sensor();
    test_thermistor_beta();
    test_rtd_cvd();
    test_linearity_independent();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
