#include "adc_metrics.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int p = 0, t = 0;
#define T(d) do { t++; printf("  %-40s", d); fflush(stdout); } while(0)
#define OK() do { p++; printf("PASS\n"); } while(0)
#define BAD(m) do { printf("FAIL: %s\n", m); } while(0)

int main(void)
{
    printf("=== test_adc_metrics ===\n");

    T("sinad_from_samples");
    {
        /* Generate a clean 1 kHz sine sampled at 10 kHz, 8-bit */
        double samples[128];
        for (int i = 0; i < 128; i++) {
            double t = (double)i / 10000.0;
            samples[i] = sin(2.0 * M_PI * 1000.0 * t);
        }
        double sinad = adc_sinad_from_samples(samples, 128, 10000.0, 1000.0);
        if (sinad > 0.0) OK(); else BAD("SINAD should be > 0");
    }

    T("snr_excl_harmonics");
    {
        double samples[128];
        for (int i = 0; i < 128; i++) {
            double t = (double)i / 10000.0;
            samples[i] = sin(2.0 * M_PI * 1000.0 * t);
        }
        double snr = adc_snr_excl_harmonics(samples, 128, 10000.0, 1000.0, 5);
        if (snr > 0.0) OK(); else BAD("SNR should be > 0");
    }

    T("thd_from_harmonics");
    {
        double harmonics[] = {0.0, -40.0, -50.0, -60.0}; /* fund + 3 harmonics in dB */
        double thd = adc_thd_from_harmonics(harmonics, 4);
        if (thd < -30.0) OK(); else BAD("THD should be < -30 dB");
    }

    T("thd_db_to_percent_roundtrip");
    {
        double pct = adc_thd_db_to_percent(-40.0);
        double db = adc_thd_percent_to_db(pct);
        if (fabs(db + 40.0) < 0.01) OK(); else BAD("roundtrip failed");
    }

    T("sfdr_from_spectrum");
    {
        /* Simulated spectrum: fund at -1 dBFS, spur at -60 dBFS */
        double spectrum[64];
        for (int i = 0; i < 64; i++) spectrum[i] = -100.0;
        spectrum[10] = -1.0; /* Fundamental */
        spectrum[20] = -60.0; /* Spur */
        uint32_t harm_bins[] = {20, 30};
        double sfdr_dBc, sfdr_dBFS;
        uint32_t spur_bin;
        adc_sfdr_from_spectrum(spectrum, 64, 10, harm_bins, 0, &sfdr_dBc, &sfdr_dBFS, &spur_bin);
        /* Spur at bin 20 not excluded (not a harmonic), should be found */
        if (spur_bin == 20) OK(); else BAD("should find spur at bin 20");
    }

    T("noise_floor");
    {
        double spectrum[64];
        for (int i = 0; i < 64; i++) spectrum[i] = -80.0;
        spectrum[0] = 0.0;
        spectrum[10] =  0.0;
        spectrum[20] = -20.0;
        uint32_t exclude[] = {10, 20};
        double nf = adc_noise_floor_dBFS(spectrum, 64, exclude, 2);
        if (nf < -70.0) OK(); else BAD("noise floor should be ~ -80 dBFS");
    }

    T("window_apply_hann");
    {
        double samples[16];
        for (int i = 0; i < 16; i++) samples[i] = 1.0;
        adc_apply_window(samples, 16, 1);
        /* Hann: w[0]=0, w[8]=0.5*(1-cos(16π/15))≈0.989 */
        if (samples[0] < 1e-9 && fabs(samples[8] - 1.0) < 0.02) OK();
        else BAD("Hann window should taper to 0 at edges");
    }

    T("window_coherent_gain");
    {
        double cg = adc_window_coherent_gain(1);
        if (fabs(cg - 0.5) < 0.01) OK();
        else BAD("Hann coherent gain should be 0.5");
    }

    T("window_enbw");
    {
        double enbw = adc_window_enbw(1);
        if (fabs(enbw - 1.5) < 0.01) OK();
        else BAD("Hann ENBW should be 1.5");
    }

    T("estimate_tone_params");
    {
        double samples[128];
        for (int i = 0; i < 128; i++) {
            double t = (double)i / 10000.0;
            samples[i] = 0.8 * sin(2.0 * M_PI * 1000.0 * t + 0.5);
        }
        double freq, ampl, phase;
        adc_estimate_tone_params(samples, 128, 10000.0, &freq, &ampl, &phase);
        if (fabs(freq - 1000.0) < 100.0 && ampl > 0.5) OK();
        else BAD("tone estimation inaccurate");
    }

    T("coherent_sampling_check");
    {
        uint32_t cycles;
        int coherent = adc_check_coherent_sampling(1000.0, 8192.0, 8192, &cycles);
        /* 1000*8192/8192 = 1000 cycles, GCD(8192,1000)=8 not 1 */
        if (!coherent) OK(); else BAD("this should not be coherent (GCD != 1)");
    }

    T("coherent_sampling_check_good");
    {
        uint32_t cycles;
        /* f_in=997Hz, fs=8192Hz, N=8192 → M=997 cycles, GCD(8192,997)=1 */
        int coherent = adc_check_coherent_sampling(997.0, 8192.0, 8192, &cycles);
        if (coherent) OK(); else BAD("this should be coherent");
    }

    T("full_dynamic_test");
    {
        double samples[256];
        for (int i = 0; i < 256; i++) {
            double t = (double)i / 8192.0;
            samples[i] = 0.5 * sin(2.0 * M_PI * 997.0 * t);
        }
        adc_dynamic_metrics_t metrics;
        adc_full_dynamic_test(samples, 256, 8192.0, 997.0, 5, &metrics);
        free(metrics.harmonic_powers);
        if (metrics.sfdr_db > 0.0) OK();
        else BAD("full dynamic test should produce SFDR");
    }

    printf("\n=== %d/%d tests passed ===\n", p, t);
    return (p == t) ? 0 : 1;
}
