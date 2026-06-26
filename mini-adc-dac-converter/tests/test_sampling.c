/**
 * @file test_sampling.c
 * @brief Tests for sampling theory implementation.
 */
#include "sampling.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int p = 0, t = 0;
#define T(desc) do { t++; printf("  %-40s", desc); fflush(stdout); } while(0)
#define OK() do { p++; printf("PASS\n"); } while(0)
#define BAD(m) do { printf("FAIL: %s\n", m); } while(0)

static double test_sine(double time, void *ctx) {
    double freq = *(double *)ctx;
    return sin(2.0 * M_PI * freq * time);
}

int main(void)
{
    printf("=== test_sampling ===\n");

    /* Aliased frequency */
    T("aliased_above_nyquist");
    {
        double fa = aliased_frequency(70.0, 100.0);
        if (fabs(fa - 30.0) < 0.01) OK();
        else BAD("70 Hz at 100 Hz fs should alias to 30 Hz");
    }

    T("aliased_below_nyquist");
    {
        double fa = aliased_frequency(30.0, 100.0);
        if (fabs(fa - 30.0) < 0.01) OK();
        else BAD("30 Hz below Nyquist should not alias");
    }

    /* Nyquist check */
    T("nyquist_satisfied_baseband");
    {
        int ok = nyquist_check(10.0, 30.0, 0.0);
        if (ok) OK(); else BAD("fs=30 > 2*10=20 should satisfy Nyquist");
    }

    T("nyquist_violated");
    {
        int ok = nyquist_check(10.0, 15.0, 0.0);
        if (!ok) OK(); else BAD("fs=15 < 2*10=20 should violate Nyquist");
    }

    /* Minimum sample rate */
    T("minimum_sample_rate_baseband");
    {
        double fs = minimum_sample_rate(10.0, 0.0);
        if (fabs(fs - 20.0) < 1e-6) OK();
        else BAD("min fs should be 2*BW");
    }

    /* Oversampling SNR gain */
    T("oversampling_snr_gain");
    {
        double gain = oversampling_snr_gain(4.0);
        double expected = 10.0 * log10(4.0);
        if (fabs(gain - expected) < 0.01) OK();
        else BAD("OSR=4 should give ~6 dB gain");
    }

    T("oversampling_effective_bits");
    {
        double bits = oversampling_effective_bits(16.0);
        double expected = 0.5 * log2(16.0);
        if (fabs(bits - expected) < 0.01) OK();
        else BAD("OSR=16 should give 2 extra bits");
    }

    /* Sinc interpolation */
    T("sinc_interpolate_midpoint");
    {
        double samples[] = {0.0, 1.0, 0.0, -1.0, 0.0};
        double val = sinc_interpolate(samples, 5, 4.0, 0.25, 10, 1);
        if (isfinite(val)) OK(); else BAD("interpolation should be finite");
    }

    /* Bandpass sampling */
    T("bandpass_sampling_range");
    {
        double fs_min, fs_max;
        uint32_t n_opt;
        bandpass_sampling_range(70e6, 90e6, &fs_min, &fs_max, &n_opt);
        if (fs_min > 0.0 && fs_max > fs_min) OK();
        else BAD("bandpass sampling range invalid");
    }

    /* Ideal sampling simulation */
    T("simulate_ideal_sampling");
    {
        double freq = 10.0;
        double samples[8];
        simulate_ideal_sampling(test_sine, &freq, 80.0, 8, 0.0, samples);
        if (fabs(samples[0]) < 1e-9) OK();
        else BAD("sin(0) should be 0");
    }

    /* CIC decimation */
    T("cic_decimate");
    {
        double input[64], output[32];
        for (int i = 0; i < 64; i++) input[i] = 1.0;
        uint32_t n_out;
        cic_decimate(input, 64, 4, 2, 1, output, &n_out);
        if (n_out > 0) OK(); else BAD("CIC should produce output");
    }

    /* FIR decimation */
    T("fir_decimate");
    {
        double input[32], coeffs[] = {0.25, 0.5, 0.25}, output[32];
        for (int i = 0; i < 32; i++) input[i] = (double)i;
        uint32_t n_out;
        fir_decimate(input, 32, 4, coeffs, 3, output, &n_out);
        if (n_out > 0) OK(); else BAD("FIR decimate should produce output");
    }

    /* FIR interpolation */
    T("fir_interpolate");
    {
        double input[] = {1.0, 0.5}, output[64];
        double coeffs[] = {1.0, 0.0, 0.0, 0.0};
        uint32_t n_out;
        fir_interpolate(input, 2, 4, coeffs, 4, output, &n_out);
        if (n_out == 8) OK(); else BAD("interpolation should give n_input*L output");
    }

    /* Rational rate conversion */
    T("rational_rate_convert");
    {
        double input[] = {1.0, 2.0, 3.0, 4.0}, output[64];
        double coeffs[] = {1.0};
        uint32_t n_out;
        rational_rate_convert(input, 4, 3, 2, coeffs, 1, output, &n_out);
        if (n_out > 0) OK(); else BAD("rate conversion should produce output");
    }

    /* AA filter sampling */
    T("simulate_adc_acquisition_basic");
    {
        double freq = 1e3;
        double coeff_b[] = {1.0};
        double coeff_a[] = {1.0, 0.0};
        uint64_t codes[16];
        simulate_adc_acquisition(test_sine, &freq, 10e3, 16, 0.0,
                                  coeff_b, coeff_a, 0, 0.0, 8, 3.3, codes);
        if (codes[0] == 0 && codes[2] > 0) OK();
        else BAD("ADC simulation output invalid");
    }

    printf("\n=== %d/%d tests passed ===\n", p, t);
    return (p == t) ? 0 : 1;
}
