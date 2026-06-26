#include "sigma_delta.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

static int p = 0, t = 0;
#define T(d) do { t++; printf("  %-40s", d); fflush(stdout); } while(0)
#define OK() do { p++; printf("PASS\n"); } while(0)
#define BAD(m) do { printf("FAIL: %s\n", m); } while(0)

int main(void)
{
    printf("=== test_sigma_delta ===\n");

    T("sdm_inband_noise_power_1st");
    {
        double ibn = sdm_inband_noise_power(1, 64.0, 1, 1.0);
        if (ibn > 0.0 && isfinite(ibn)) OK();
        else BAD("IBN should be positive");
    }

    T("sdm_dynamic_range_1st_osr64");
    {
        double dr = sdm_dynamic_range_db(1, 64.0, 1, 1.0);
        if (dr > 30.0) OK(); else BAD("DR should be > 30 dB for OSR=64");
    }

    T("sdm_dynamic_range_2nd_osr64");
    {
        double dr = sdm_dynamic_range_db(2, 64.0, 1, 1.0);
        if (dr > 50.0) OK(); else BAD("DR should be > 50 dB for 2nd order");
    }

    T("sdm_equivalent_enob");
    {
        double enob = sdm_equivalent_enob(2, 64.0, 1, 1.0);
        if (enob > 8.0) OK(); else BAD("ENOB should be > 8 for 2nd/OSR64");
    }

    T("ntf_pure_diff_1st");
    {
        double num[2], den[2];
        sdm_design_ntf_pure_diff(1, num, den);
        if (fabs(num[0] - 1.0) < 1e-9 && fabs(num[1] + 1.0) < 1e-9) OK();
        else BAD("1st-order diff NTF should be [1, -1]");
    }

    T("ntf_pure_diff_2nd");
    {
        double num[3] = {0}, den[3] = {0};
        sdm_design_ntf_pure_diff(2, num, den);
        if (fabs(num[0] - 1.0) < 1e-9 && fabs(num[1] + 2.0) < 1e-9 && fabs(num[2] - 1.0) < 1e-9) OK();
        else BAD("2nd order NTF should be [1, -2, 1]");
    }

    T("ntf_butterworth");
    {
        double num[3] = {0}, den[3] = {0};
        sdm_design_ntf_butterworth(2, 20e3, 2.8224e6, num, den);
        if (fabs(num[0] - 1.0) < 1e-9) OK();
        else BAD("Butterworth NTF design failed");
    }

    T("sdm_first_order_process");
    {
        sdm_spec_t spec = {.order=1, .quantizer_bits=1, .osr=64, .fs_hz=2.8224e6, .bandwidth_hz=20e3, .v_ref=1.0};
        sdm_state_t state;
        if (sdm_init_state(&spec, &state) == 0) {
            int64_t y = sdm_process_first_order(0.1, &state, 1.0);
            if (y == 1 || y == -1) { sdm_free_state(&state); OK(); }
            else { sdm_free_state(&state); BAD("output should be ±1"); }
        } else BAD("init failed");
    }

    T("sdm_second_order_process");
    {
        sdm_spec_t spec = {.order=2, .quantizer_bits=1, .osr=64, .fs_hz=2.8224e6, .bandwidth_hz=20e3, .v_ref=1.0};
        sdm_state_t state;
        if (sdm_init_state(&spec, &state) == 0) {
            int64_t y = sdm_process_second_order(0.1, &state, 1.0);
            if (y == 1 || y == -1) { sdm_free_state(&state); OK(); }
            else { sdm_free_state(&state); BAD("output should be ±1"); }
        } else BAD("init failed");
    }

    T("sdm_general_process");
    {
        sdm_spec_t spec = {.order=1, .quantizer_bits=1, .osr=64, .fs_hz=1e6, .bandwidth_hz=10e3, .v_ref=1.0};
        sdm_state_t state;
        if (sdm_init_state(&spec, &state) == 0) {
            int64_t y = sdm_process_general(0.3, &state, 1.0);
            sdm_free_state(&state);
            if (y != 0) OK(); else BAD("general process should output nonzero");
        } else BAD("init failed");
    }

    T("sdm_mash_1_1_process");
    {
        sdm_spec_t spec1 = {.order=1, .quantizer_bits=1, .v_ref=1.0};
        sdm_spec_t spec2 = {.order=1, .quantizer_bits=1, .v_ref=1.0};
        sdm_state_t st1, st2;
        if (sdm_init_state(&spec1, &st1) == 0 && sdm_init_state(&spec2, &st2) == 0) {
            int64_t out;
            sdm_mash_1_1_process(0.2, &st1, &st2, 1.0, &out);
            sdm_free_state(&st1); sdm_free_state(&st2);
            /* MASH 1-1 output is multi-valued; just check it's within reasonable range */
            if (out >= -3 && out <= 3) OK(); else BAD("MASH should produce output");
        } else BAD("init failed");
    }

    T("cic_filter_design");
    {
        uint32_t len;
        double *coeffs = sdm_design_cic_filter(8, 2, 1, &len);
        if (coeffs && len > 0) { free(coeffs); OK(); }
        else BAD("CIC design failed");
    }

    T("cic_compensation");
    {
        uint32_t n_taps;
        double *coeffs = sdm_design_cic_compensation(64, 2, 1, 20e3, 44.1e3, &n_taps);
        if (coeffs && n_taps > 0) { free(coeffs); OK(); }
        else BAD("CIC compensation failed");
    }

    T("sdm_stable_input_limit");
    {
        sdm_spec_t spec = {.order=2, .quantizer_bits=1};
        double limit = sdm_stable_input_limit(&spec);
        if (limit > 0.0 && limit < 1.0) OK();
        else BAD("stable input limit should be in (0, 1)");
    }

    T("lee_criterion");
    {
        double num[] = {1.0, -1.0, 0.0};
        double den[] = {1.0, -0.5, 0.0};
        double max_gain;
        int stable = sdm_check_lee_criterion(num, den, 2, 100, &max_gain);
        if (stable && max_gain > 0.0) OK();
        else BAD("1st-order diff NTF should satisfy Lee's criterion");
    }

    T("sdm_simulate_adc");
    {
        sdm_spec_t spec = {.order=1, .quantizer_bits=1, .osr=16, .fs_hz=640e3, .bandwidth_hz=20e3, .v_ref=1.0};
        double decimated[1024];
        int64_t *osr = (int64_t *)malloc(640*sizeof(int64_t));
        uint32_t n_dec;
        uint64_t n_osr;
        int r = sdm_simulate_adc(&spec, 1e3, 0.5, 0.001, decimated, &n_dec, osr, &n_osr);
        if (r == 0 && n_osr > 0) OK(); else BAD("simulation failed");
        free(osr);
    }

    printf("\n=== %d/%d tests passed ===\n", p, t);
    return (p == t) ? 0 : 1;
}
