/**
 * @file test_quantization.c
 * @brief Tests for quantization theory implementation.
 */
#include "quantization.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

static int p = 0, t = 0;
#define T(desc) do { t++; printf("  %-40s", desc); fflush(stdout); } while(0)
#define OK() do { p++; printf("PASS\n"); } while(0)
#define BAD(m) do { printf("FAIL: %s\n", m); } while(0)

int main(void)
{
    printf("=== test_quantization ===\n");

    /* Uniform quantization */
    T("quantize_uniform_mid_tread");
    {
        int64_t c = quantize_uniform(0.75, 1.0, QUANTIZER_MID_TREAD);
        if (c == 1) OK(); else BAD("0.75 with step=1 should be code 1");
    }

    T("quantize_uniform_negative");
    {
        int64_t c = quantize_uniform(-0.6, 1.0, QUANTIZER_MID_TREAD);
        if (c == -1) OK(); else BAD("-0.6 with step=1 should be code -1");
    }

    T("dequantize_uniform");
    {
        double v = dequantize_uniform(3, 0.5, QUANTIZER_MID_TREAD);
        if (fabs(v - 1.5) < 1e-9) OK(); else BAD("dequantization wrong");
    }

    /* Quantization noise power */
    T("noise_power_q_sq_over_12");
    {
        double np = quantization_noise_power(0.1);
        double expected = 0.01 / 12.0;
        if (fabs(np - expected) < 1e-12) OK();
        else BAD("noise power should be q²/12");
    }

    T("sqnr_sinusoid");
    {
        /* Full-scale sine with step=q: amplitude=Vref/2 = (2^N)*q/2 */
        /* For N=1, amplitude = q, rms = q/√2 */
        double signal_rms = 0.1 / sqrt(2.0);
        double sqnr = quantization_sqnr(signal_rms, 0.1);
        if (sqnr > 0.0) OK(); else BAD("SQNR should be positive");
    }

    /* Lloyd-Max quantizer */
    T("lloyd_max_gaussian_design");
    {
        /* Simple 5-point Gaussian PDF */
        double grid[] = {-2.0, -1.0, 0.0, 1.0, 2.0};
        double pdf[]  = {0.054, 0.242, 0.399, 0.242, 0.054};
        nonuniform_quantizer_t q;
        int r = lloyd_max_quantizer_design(pdf, grid, 5, 2, &q, 50, 1e-6);
        if (r == 0 && q.num_levels == 2 && q.mse >= 0.0) {
            nonuniform_quantizer_free(&q);
            OK();
        } else BAD("Lloyd-Max design failed");
    }

    /* Dithering */
    T("subtractive_dither");
    {
        double di, corr;
        dither_subtractive(0.5, 0.02, 0.1, &di, &corr);
        if (fabs(di - 0.52) < 1e-9 && fabs(corr - 0.02) < 1e-9) OK();
        else BAD("subtractive dither wrong");
    }

    T("nonsubtractive_tpdf_dither");
    {
        double di = dither_nonsubtractive_tpdf(0.5, 0.02, -0.01, 0.1);
        if (fabs(di - 0.51) < 1e-9) OK();
        else BAD("TPDF dither wrong");
    }

    /* Noise shaping */
    T("noise_shaping_first_order");
    {
        double e_prev = 0.0;
        int64_t c = noise_shape_first_order(0.55, 0.1, &e_prev);
        /* Output should be quantized to nearest step */
        double expected = round(0.55 / 0.1) * 0.1;
        double actual = (double)c * 0.1;
        if (fabs(actual - expected) < 0.01) OK();
        else BAD("first-order noise shaping wrong");
    }

    T("noise_shaping_second_order");
    {
        double e1 = 0.0, e2 = 0.0;
        int64_t c = noise_shape_second_order(0.55, 0.1, &e1, &e2);
        if (c != 0) OK(); else BAD("second-order shaping should produce output");
    }

    /* Quantization error statistics */
    T("error_stats");
    {
        double orig[] = {0.1, 0.2, 0.3, 0.4, 0.5};
        double quant[] = {0.1, 0.2, 0.3, 0.4, 0.5};
        quantization_error_stats_t stats;
        quantization_error_stats(orig, quant, 5, 0.1, &stats);
        if (stats.mean_error < 1e-9 && stats.var_error < 1e-9) OK();
        else BAD("zero error stats wrong");
    }

    /* μ-law */
    T("mu_law_compress_zero");
    {
        double y = mu_law_compress(0.0, 255.0);
        if (fabs(y) < 1e-9) OK(); else BAD("mu-law(0) should be 0");
    }

    T("mu_law_roundtrip");
    {
        double x = 0.5;
        double y = mu_law_compress(x, 255.0);
        double x2 = mu_law_expand(y, 255.0);
        if (fabs(x - x2) < 0.01) OK(); else BAD("mu-law roundtrip");
    }

    /* A-law */
    T("a_law_compress_zero");
    {
        double y = a_law_compress(0.0, 87.6);
        if (fabs(y) < 1e-9) OK(); else BAD("A-law(0) should be 0");
    }

    T("a_law_roundtrip");
    {
        double x = 0.5;
        double y = a_law_compress(x, 87.6);
        double x2 = a_law_expand(y, 87.6);
        if (fabs(x - x2) < 0.05) OK(); else BAD("A-law roundtrip");
    }

    /* Overload / optimal step size */
    T("overload_distortion");
    {
        double ov = overload_distortion_power(1.0, 0.2);
        if (ov >= 0.0) OK(); else BAD("overload distortion should be ≥ 0");
    }

    T("optimal_step_size");
    {
        double q = optimal_step_size(1.0, 8);
        if (q > 0.0) OK(); else BAD("step size should be positive");
    }

    printf("\n=== %d/%d tests passed ===\n", p, t);
    return (p == t) ? 0 : 1;
}
