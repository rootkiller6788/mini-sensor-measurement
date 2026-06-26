/**
 * @file test_adc_core.c
 * @brief Tests for ADC core definitions and fundamental functions.
 */
#include "adc_dac_core.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  TEST %-40s", name); \
    fflush(stdout); \
} while(0)

#define PASS() do { tests_passed++; printf(" PASS\n"); } while(0)
#define FAIL(msg) do { printf(" FAIL: %s\n", msg); } while(0)

int main(void)
{
    printf("=== test_adc_core ===\n");

    /* L4: Ideal ADC SNR formula */
    TEST("ideal_snr_6bit");
    {
        double snr = adc_ideal_snr_db(6);
        double expected = 6.02 * 6.0 + 1.76;
        if (fabs(snr - expected) < 0.5) PASS();
        else FAIL("SNR mismatch");
    }

    TEST("ideal_snr_12bit");
    {
        double snr = adc_ideal_snr_db(12);
        double expected = 6.02 * 12.0 + 1.76;
        if (fabs(snr - expected) < 0.5) PASS();
        else FAIL("SNR mismatch");
    }

    TEST("ideal_snr_16bit");
    {
        double snr = adc_ideal_snr_db(16);
        double expected = 6.02 * 16.0 + 1.76;
        if (fabs(snr - expected) < 0.5) PASS();
        else FAIL("SNR mismatch");
    }

    TEST("ideal_snr_zero_bits_returns_nan");
    {
        double snr = adc_ideal_snr_db(0);
        if (isnan(snr)) PASS();
        else FAIL("Should return NaN for 0 bits");
    }

    /* SNR from powers */
    TEST("snr_from_powers");
    {
        double snr = adc_snr_from_powers(1.0, 0.01);
        if (fabs(snr - 20.0) < 0.1) PASS();
        else FAIL("SNR should be 20 dB for 100:1 power ratio");
    }

    TEST("snr_from_powers_invalid");
    {
        double snr = adc_snr_from_powers(0.0, 1.0);
        if (isnan(snr)) PASS();
        else FAIL("Should return NaN for zero signal");
    }

    /* ENOB */
    TEST("enob_from_sinad");
    {
        double enob = adc_enob_from_sinad(74.0);
        double expected = (74.0 - 1.76) / 6.02;
        if (fabs(enob - expected) < 0.01) PASS();
        else FAIL("ENOB mismatch");
    }

    /* Jitter SNR limit */
    TEST("jitter_snr_100MHz");
    {
        double snr = adc_jitter_snr_limit(100e6, 1.0);
        double expected = -20.0 * log10(2.0 * M_PI * 100e6 * 1e-12);
        if (fabs(snr - expected) < 0.1) PASS();
        else FAIL("Jitter SNR mismatch");
    }

    TEST("jitter_snr_invalid_freq");
    {
        double snr = adc_jitter_snr_limit(0.0, 1.0);
        if (isnan(snr)) PASS();
        else FAIL("Should return NaN");
    }

    /* Total input noise */
    TEST("total_input_noise");
    {
        double noise = adc_total_input_noise(0.001, 10.0, 25.0);
        if (noise > 0.0 && isfinite(noise)) PASS();
        else FAIL("Noise should be positive and finite");
    }

    /* Code to voltage */
    TEST("code_to_voltage_full_scale");
    {
        double v = adc_code_to_voltage(4095, 3.3, 12);
        double expected = 3.3 * 4095.0 / 4096.0;
        if (fabs(v - expected) < 1e-6) PASS();
        else FAIL("Full-scale voltage mismatch");
    }

    TEST("code_to_voltage_zero");
    {
        double v = adc_code_to_voltage(0, 3.3, 12);
        if (fabs(v - 0.0) < 1e-9) PASS();
        else FAIL("Zero code should give 0V");
    }

    /* Voltage to code */
    TEST("voltage_to_code_mid");
    {
        uint64_t code = adc_voltage_to_code(1.65, 3.3, 12);
        if (code == 2048) PASS();
        else FAIL("Mid-scale should be 2048");
    }

    TEST("voltage_to_code_clamped");
    {
        uint64_t code = adc_voltage_to_code(5.0, 3.3, 12);
        if (code == 4095) PASS();
        else FAIL("Should clamp to max code");
    }

    /* Two-point calibration */
    TEST("two_point_calibrate");
    {
        uint64_t cal = adc_two_point_calibrate(1000, 10, 4080, 12);
        if (cal > 0 && cal <= 4095) PASS();
        else FAIL("Calibration failed");
    }

    /* DNL computation */
    TEST("dnl_ideal_ramp");
    {
        double transitions[9]; /* 3-bit: 9 transitions */
        for (int i = 0; i < 9; i++) transitions[i] = 0.1 * (double)i;
        double dnl[7]; /* 2^3 - 1 = 7 */
        adc_compute_dnl(transitions, 3, dnl);
        int ok = 1;
        for (int i = 0; i < 7; i++) {
            if (fabs(dnl[i]) > 0.01) ok = 0;
        }
        if (ok) PASS();
        else FAIL("Ideal DNL should be near zero");
    }

    /* INL from DNL */
    TEST("inl_from_ideal_dnl");
    {
        double dnl[7] = {0, 0, 0, 0, 0, 0, 0};
        double inl[8];
        adc_compute_inl_from_dnl(dnl, 3, inl);
        int ok = 1;
        for (int i = 0; i < 8; i++) ok &= (fabs(inl[i]) < 0.001);
        if (ok) PASS();
        else FAIL("INL should be zero for ideal DNL");
    }

    /* String names */
    TEST("adc_arch_name");
    {
        const char *name = adc_arch_name(ADC_ARCH_SAR);
        if (name && name[0] == 'S') PASS();
        else FAIL("SAR name wrong");
    }

    TEST("dac_arch_name");
    {
        const char *name = dac_arch_name(DAC_ARCH_R2R_LADDER);
        if (name && strstr(name, "R-2R")) PASS();
        else FAIL("R-2R name wrong");
    }

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
