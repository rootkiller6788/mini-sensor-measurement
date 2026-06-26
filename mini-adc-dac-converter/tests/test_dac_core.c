/**
 * @file test_dac_core.c
 * @brief Tests for DAC core implementations.
 */
#include "adc_dac_core.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

static int passed = 0, total = 0;
#define T(desc) do { total++; printf("  %-40s", desc); fflush(stdout); } while(0)
#define OK() do { passed++; printf("PASS\n"); } while(0)
#define BAD(m) do { printf("FAIL: %s\n", m); } while(0)

/* Forward declarations from dac_core.c */
double r2r_dac_voltage(uint64_t code, uint32_t n_bits, double v_ref);
double r2r_output_impedance(double r_ohms);
double r2r_bit_weight(uint32_t bit_index, uint32_t n_bits, double v_ref);
double binary_weighted_dac_voltage(uint64_t code, uint32_t n_bits, double v_ref, double r_fb, double r_base);
double binary_weighted_resistor_ratio(uint32_t n_bits);
double string_dac_voltage(uint64_t code, uint32_t n_bits, double v_ref);
double string_dac_resistor_value(double v_ref, double i_bias, uint32_t n_bits);
double dac_glitch_energy_triangular(double amp_glitch_v, double dur_glitch_s);
double dac_glitch_from_skew(double v_ref, double t_skew_s, double r_out_ohm);
double dac_settling_rc(double t_current, double tau_s, double v_final);
double dac_settling_time(double tau_s, double v_step, double v_error_max);
double dac_worst_case_settling_time(double tau_s, uint32_t n_bits);
double segmented_dac_voltage(uint64_t code, uint32_t m_bits, uint32_t k_bits, double v_ref);
void segmented_split_code(uint64_t code, uint32_t m, uint32_t k, uint64_t *co, uint64_t *fi);
double current_steering_dac_voltage(uint64_t code, uint32_t n_bits, double i_unit, double r_load);
double current_steering_sfdr_estimate(uint32_t n_bits, double mismatch_std);
void dac_compute_dnl(const double *v_out, uint32_t n_bits, double v_ref, double *dnl);
void dac_compute_inl_from_dnl(const double *dnl, uint32_t n_bits, double *inl);
int dac_is_monotonic(const double *dnl, uint64_t n_codes);

int main(void)
{
    printf("=== test_dac_core ===\n");

    /* R-2R ladder */
    T("r2r_full_scale");
    {
        double v = r2r_dac_voltage(255, 8, 5.0);
        if (fabs(v - 5.0 * 255.0 / 256.0) < 1e-6) OK();
        else BAD("voltage mismatch");
    }

    T("r2r_zero_code");
    {
        double v = r2r_dac_voltage(0, 8, 5.0);
        if (fabs(v) < 1e-9) OK();
        else BAD("zero code not 0V");
    }

    T("r2r_output_impedance");
    {
        double z = r2r_output_impedance(1000.0);
        if (fabs(z - 1000.0) < 1e-6) OK();
        else BAD("impedance mismatch");
    }

    T("r2r_bit_weight_msb");
    {
        double w = r2r_bit_weight(7, 8, 5.0);
        if (fabs(w - 2.5) < 1e-6) OK();
        else BAD("MSB weight should be Vref/2");
    }

    T("r2r_bit_weight_lsb");
    {
        double w = r2r_bit_weight(0, 8, 5.0);
        double expected = 5.0 / 256.0;
        if (fabs(w - expected) < 1e-9) OK();
        else BAD("LSB weight mismatch");
    }

    /* Binary-weighted DAC */
    T("binary_weighted_voltage");
    {
        double v = binary_weighted_dac_voltage(128, 8, 5.0, 1000.0, 1000.0);
        double expected = -5.0 * 128.0 / 256.0;
        if (fabs(v - expected) < 1e-6) OK();
        else BAD("voltage wrong");
    }

    T("binary_weighted_resistor_ratio");
    {
        double ratio = binary_weighted_resistor_ratio(8);
        if (fabs(ratio - 128.0) < 1e-6) OK();
        else BAD("ratio should be 2^7 = 128");
    }

    /* String DAC */
    T("string_dac_voltage");
    {
        double v = string_dac_voltage(512, 10, 2.5);
        if (fabs(v - 1.25) < 1e-6) OK();
        else BAD("half-scale wrong");
    }

    T("string_dac_resistor_value");
    {
        double r = string_dac_resistor_value(5.0, 0.001, 8);
        if (r > 0.0) OK();
        else BAD("resistor value should be positive");
    }

    /* Glitch energy */
    T("glitch_energy_triangular");
    {
        double e = dac_glitch_energy_triangular(0.1, 1e-9);
        double expected = 0.5 * 1e-9 * 0.1;
        if (fabs(e - expected) < 1e-18) OK();
        else BAD("glitch energy wrong");
    }

    T("glitch_from_skew");
    {
        double e = dac_glitch_from_skew(3.3, 100e-12, 1000.0);
        if (e > 0.0) OK();
        else BAD("glitch energy should be positive");
    }

    /* Settling */
    T("settling_rc_mid");
    {
        double v = dac_settling_rc(0.693147, 1.0, 1.0);
        if (fabs(v - 0.5) < 0.01) OK();
        else BAD("settling at τ*ln(2) should be 0.5*V_final");
    }

    T("settling_time");
    {
        double t = dac_settling_time(1e-9, 1.0, 0.001);
        if (t > 0.0) OK();
        else BAD("settling time should be positive");
    }

    T("worst_case_settling");
    {
        double t = dac_worst_case_settling_time(1e-9, 12);
        double expected = 1e-9 * 13.0 * log(2.0);
        if (fabs(t - expected) < 1e-12) OK();
        else BAD("worst-case settling time wrong");
    }

    /* Segmented DAC */
    T("segmented_dac_voltage");
    {
        double v = segmented_dac_voltage(512, 4, 6, 5.0);
        if (fabs(v - 5.0 * 512.0 / 1024.0) < 1e-6) OK();
        else BAD("voltage mismatch");
    }

    T("segmented_split_code");
    {
        uint64_t coarse, fine;
        segmented_split_code(0b1100110, 4, 6, &coarse, &fine);
        /* 0b1100110: lower 6 bits = 0b100110 = 38, upper 4 bits = 0b1 = 1 */
        if (coarse == 1 && fine == 38) OK();
        else BAD("code splitting wrong");
    }

    /* Current-steering DAC */
    T("current_steering_dac");
    {
        double v = current_steering_dac_voltage(128, 8, 1e-6, 1000.0);
        if (isfinite(v)) OK();
        else BAD("voltage should be finite");
    }

    T("current_steering_sfdr");
    {
        double sfdr = current_steering_sfdr_estimate(10, 0.01);
        if (sfdr > 40.0) OK();
        else BAD("SFDR should be > 40 dB");
    }

    /* DNL/INL for DAC */
    T("dac_dnl_ideal");
    {
        double v_out[9];
        for (int i = 0; i < 9; i++) v_out[i] = 3.3 * i / 8.0;
        double dnl[7];
        dac_compute_dnl(v_out, 3, 3.3, dnl);
        int ok = 1;
        for (int i = 0; i < 7; i++) if (fabs(dnl[i]) > 0.01) ok = 0;
        if (ok) OK();
        else BAD("DNL not ideal");
    }

    T("dac_inl_from_dnl");
    {
        double dnl[7] = {0.1, -0.05, 0.2, -0.1, 0.05, 0.0, -0.05};
        double inl[8];
        dac_compute_inl_from_dnl(dnl, 3, inl);
        if (fabs(inl[0]) < 1e-9) OK();
        else BAD("INL[0] should be 0");
    }

    T("dac_monotonic_ideal");
    {
        double dnl[7] = {0, 0, 0, 0, 0, 0, 0};
        if (dac_is_monotonic(dnl, 8)) OK();
        else BAD("ideal DAC should be monotonic");
    }

    T("dac_non_monotonic");
    {
        double dnl[7] = {0, 0, -2.0, 0, 0, 0, 0};
        if (!dac_is_monotonic(dnl, 8)) OK();
        else BAD("DNL < -1 should be non-monotonic");
    }

    printf("\n=== %d/%d tests passed ===\n", passed, total);
    return (passed == total) ? 0 : 1;
}
