/**
 * @file test_sensor_bridge_cal.c
 * @brief Tests for Wheatstone bridge and calibration algorithms
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "../include/sensor_bridge.h"
#include "../include/sensor_calibration.h"
#include "../include/sensor_signal_cond.h"

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST %s ... ", #name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); return; } while(0)
#define ASSERT_TRUE(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)
#define ASSERT_FLOAT_EQ(a, b, tol, msg) do { \
    if (fabs((a) - (b)) > (tol)) { \
        printf("FAIL: %s (got %.6g, expected %.6g)\n", msg, (a), (b)); \
        return; \
    } \
} while(0)

/* ──── Bridge: Balanced Bridge Zero Output ──── */

static void test_bridge_balanced_zero(void) {
    TEST(bridge_balanced_zero);
    wheatstone_bridge_t b;
    bridge_init_balanced(&b, 350.0, 5.0);  /* 350Ω strain gauge, 5V excitation */

    double V_out = bridge_output_voltage_exact(&b);
    ASSERT_FLOAT_EQ(V_out, 0.0, 1e-10, "balanced bridge output ~0V");
    PASS();
}

/* ──── Bridge: Quarter Bridge Output ──── */

static void test_bridge_quarter(void) {
    TEST(bridge_quarter_output);
    wheatstone_bridge_t b;
    bridge_init_balanced(&b, 350.0, 5.0);
    bridge_set_quarter_active(&b, 350.0, 353.5, 1);  /* +1% ΔR (≡ ~5000 µε for GF=2) */

    double V_out = bridge_output_voltage_exact(&b);
    /* V_out ≈ V_ex·ΔR/(4R₀) ≈ 5·3.5/(4·350) = 12.5 mV */
    ASSERT_TRUE(V_out > 0.010 && V_out < 0.015, "quarter bridge output ~12.5mV");
    PASS();
}

/* ──── Bridge: Full Bridge Output ──── */

static void test_bridge_full(void) {
    TEST(bridge_full_output);
    wheatstone_bridge_t b;
    bridge_init_balanced(&b, 350.0, 5.0);
    bridge_set_full_pushpull(&b, 350.0, 7.0);  /* 2% ΔR */

    double V_out = bridge_output_voltage_exact(&b);
    /* V_out ≈ V_ex·ΔR/R₀ ≈ 5·7/350 = 100 mV */
    ASSERT_TRUE(V_out > 0.090 && V_out < 0.110, "full bridge output ~100mV");
    PASS();
}

/* ──── Bridge: Nonlinearity Error ──── */

static void test_bridge_nonlinearity_error(void) {
    TEST(bridge_nonlinearity_error);
    wheatstone_bridge_t b;
    bridge_init_balanced(&b, 350.0, 5.0);

    /* For ΔR/R₀ = 0.05 (5%), nonlinearity should be ~2.5% */
    double err = bridge_nonlinearity_error(&b, 350.0, 17.5);
    ASSERT_TRUE(err > 1.0 && err < 5.0, "nonlinearity error in expected range");
    PASS();
}

/* ──── Bridge: Strain Calculation ──── */

static void test_bridge_strain(void) {
    TEST(bridge_strain);
    /* For quarter bridge with GF=2.0, V_out/V_ex=0.0025:
     * ε = 4·0.0025/2.0 = 0.005 = 5000 µε */
    double strain = strain_from_bridge_output(0.0125, 5.0, 2.0, BRIDGE_QUARTER);
    ASSERT_FLOAT_EQ(strain, 0.005, 0.0001, "strain = 5000 µε");

    /* Nonlinearity correction: ε_true = ε_app/(1 - GF·ε_app/2) */
    double strain_app = 0.01;  /* 10000 µε apparent */
    double strain_corr = bridge_strain_nonlinearity_corrected(strain_app, 2.0);
    ASSERT_TRUE(strain_corr > strain_app, "corrected strain > apparent strain");
    PASS();
}

/* ──── Calibration: Two-Point ──── */

static void test_cal_two_point(void) {
    TEST(cal_two_point);
    double slope, intercept;
    /* y_ref = 2·y_raw + 3 */
    ASSERT_TRUE(calibration_two_point(0.0, 3.0, 5.0, 13.0, &slope, &intercept) == 0,
                "two-point calibration succeeds");
    ASSERT_FLOAT_EQ(slope, 2.0, 0.001, "slope = 2");
    ASSERT_FLOAT_EQ(intercept, 3.0, 0.001, "intercept = 3");
    PASS();
}

/* ──── Calibration: Linear Least-Squares ──── */

static void test_cal_linear_ls(void) {
    TEST(cal_linear_ls);
    double x[] = {0, 1, 2, 3, 4};
    double y[] = {1.1, 3.0, 5.2, 7.1, 8.9};  /* y ≈ 2x + 1 with noise */
    double slope, intercept, r2, rms;

    ASSERT_TRUE(calibration_linear_ls(x, y, 5, &slope, &intercept, &r2, &rms) == 0,
                "linear LS succeeds");
    ASSERT_TRUE(slope > 1.9 && slope < 2.1, "slope ≈ 2");
    ASSERT_TRUE(intercept > 0.9 && intercept < 1.3, "intercept ≈ 1");
    ASSERT_TRUE(r2 > 0.99, "R² close to 1");
    PASS();
}

/* ──── Calibration: Polynomial LS ──── */

static void test_cal_polynomial_ls(void) {
    TEST(cal_polynomial_ls);
    /* Data from y = 1 + 2x + 0.5x² with small noise */
    size_t n = 10;
    double *x = (double *)calloc(n, sizeof(double));
    double *y = (double *)calloc(n, sizeof(double));
    assert(x && y);

    for (size_t i = 0; i < n; i++) {
        x[i] = (double)i;
        y[i] = 1.0 + 2.0 * x[i] + 0.5 * x[i] * x[i] + ((double)i - 5.0) * 0.01;
    }

    double coeff[3];
    double rms;
    ASSERT_TRUE(calibration_polynomial_ls(x, y, n, 2, coeff, &rms) == 0,
                "polynomial LS succeeds");
    ASSERT_FLOAT_EQ(coeff[0], 1.0, 0.05, "a0 ≈ 1");
    ASSERT_FLOAT_EQ(coeff[1], 2.0, 0.05, "a1 ≈ 2");
    ASSERT_FLOAT_EQ(coeff[2], 0.5, 0.01, "a2 ≈ 0.5");

    free(x); free(y);
    PASS();
}

/* ──── Calibration: Temperature Compensation ──── */

static void test_cal_temp_comp(void) {
    TEST(cal_temp_comp);
    /* Sensitivity varies linearly with temperature:
     *   S(0°C)=1.0, S(50°C)=1.05, S(100°C)=1.10
     *   → α = 0.001 /°C
     * Offset: b(0°C)=0.0, b(50°C)=0.2, b(100°C)=0.4
     *   → β = 0.004 /°C */
    double temps[] = {0.0, 50.0, 100.0};
    double slopes[] = {1.0, 1.05, 1.10};
    double offsets[] = {0.0, 0.2, 0.4};
    double S0, alpha, b0, beta;

    ASSERT_TRUE(calibration_temp_compensation(temps, slopes, offsets, 3, 25.0,
                                               &S0, &alpha, &b0, &beta) == 0,
                "temp comp succeeds");
    /* Linear fit: b(T) = 0.004*T, at T0=25°C: b0 = 0.004*25 = 0.1 */
    ASSERT_FLOAT_EQ(S0, 1.025, 0.01, "S0 near 1.025 at 25°C");
    ASSERT_FLOAT_EQ(b0, 0.1, 0.01, "b0 near 0.1 at 25°C (0.004*25)");

    /* Apply compensation */
    double y_comp = calibration_apply_temp_comp(1.5, 75.0, 25.0, S0, alpha, b0, beta);
    ASSERT_TRUE(y_comp > 0.0, "compensated value is valid");
    PASS();
}

/* ──── Signal Conditioning: 4-20 mA ──── */

static void test_current_loop(void) {
    TEST(current_loop);
    /* Linear mapping: 0-100°C → 4-20 mA */
    double I = current_loop_sensor_to_current(50.0, 0.0, 100.0);
    ASSERT_FLOAT_EQ(I, 12.0, 0.01, "mid-scale = 12mA");

    double val = current_loop_current_to_sensor(I, 0.0, 100.0);
    ASSERT_FLOAT_EQ(val, 50.0, 0.01, "round-trip preserves value");

    /* Sense resistor for 5V ADC */
    double R = current_loop_sense_resistor(5.0);
    ASSERT_FLOAT_EQ(R, 250.0, 0.1, "R_sense = 250Ω for 5V ADC");
    PASS();
}

/* ──── Signal Conditioning: Anti-Alias Filter ──── */

static void test_anti_alias(void) {
    TEST(anti_alias);
    double fc = anti_alias_fc_from_fs(1000.0, 10.0);
    ASSERT_FLOAT_EQ(fc, 50.0, 0.1, "fc = fs/(2·OSR) = 50Hz");

    /* Attenuation at Nyquist: 1st order with OSR=10 → ~20 dB */
    double att = anti_alias_attenuation_at_nyquist(50.0, 1000.0, 1);
    ASSERT_TRUE(att < -15.0 && att > -25.0, "~20dB attenuation at Nyquist");

    /* Filter order needed for 60 dB attenuation with fc=100Hz, fs=1kHz */
    int order = anti_alias_filter_order_required(1000.0, 100.0, 60.0);
    ASSERT_TRUE(order >= 3, "need ≥3rd order for 60dB alias rejection");
    PASS();
}

/* ──── Signal Conditioning: Digital Filters ──── */

static void test_moving_average(void) {
    TEST(moving_average);
    double input[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    double output[8];
    ASSERT_TRUE(signal_moving_average(input, output, 8, 4) == 0, "moving avg succeeds");
    /* With window=4: y[7] = (5+6+7+8)/4 = 6.5 */
    ASSERT_FLOAT_EQ(output[7], 6.5, 0.01, "moving avg output correct");
    PASS();
}

static void test_exponential_filter(void) {
    TEST(exponential_filter);
    /* α=0.5: y[0]=0, x[0]=10 → y[1]=5 */
    double y = signal_exponential_filter(10.0, 0.0, 0.5);
    ASSERT_FLOAT_EQ(y, 5.0, 0.01, "α=0.5: y = 0.5·10+0.5·0 = 5");
    /* Next step: x=10, y_prev=5 → y=7.5 */
    y = signal_exponential_filter(10.0, y, 0.5);
    ASSERT_FLOAT_EQ(y, 7.5, 0.01, "second step converges toward 10");
    PASS();
}

/* ──── Calibration: Cross-Sensitivity ──── */

static void test_cross_sensitivity(void) {
    TEST(cross_sensitivity);
    /* Sensor reads 5.0V, but 1.0V is from temperature interference.
     * S_x=1.0 V/unit, S_z=0.1 V/°C, z=30°C, b=0.5V */
    double x = calibration_cross_sensitivity_comp(5.0, 30.0, 1.0, 0.1, 0.5);
    ASSERT_FLOAT_EQ(x, 1.5, 0.1, "cross-sensitivity compensation correct");
    PASS();
}

int main(void) {
    printf("=== Sensor Bridge & Calibration Tests ===\n\n");

    test_bridge_balanced_zero();
    test_bridge_quarter();
    test_bridge_full();
    test_bridge_nonlinearity_error();
    test_bridge_strain();
    test_cal_two_point();
    test_cal_linear_ls();
    test_cal_polynomial_ls();
    test_cal_temp_comp();
    test_current_loop();
    test_anti_alias();
    test_moving_average();
    test_exponential_filter();
    test_cross_sensitivity();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
