/**
 * demo_filter_response.c - Filter Response Visualization
 *
 * Demonstrates Butterworth, Chebyshev, Bessel filter responses
 * (magnitude, phase, group delay) for comparison.
 *
 * L6: Filter design trade-off visualization
 */

#include <stdio.h>
#include <math.h>
#include "sigcond_filter.h"

int main(void) {
    printf("=== Filter Response Comparison Demo ===\n\n");

    double fc = 1000.0;  /* 1 kHz cutoff */
    unsigned order = 4;

    /* Butterworth biquads */
    double bw_w0[10], bw_Q[10];
    unsigned bw_ns;
    butterworth_biquads(order, bw_w0, bw_Q, &bw_ns);
    printf("Butterworth %dth-order (maximally flat):\n", order);
    for (unsigned k = 0; k < bw_ns; k++)
        printf("  Section %u: w0=%.4f, Q=%.4f\n", k+1, bw_w0[k], bw_Q[k]);

    /* Chebyshev biquads */
    double ch_w0[10], ch_Q[10];
    unsigned ch_ns;
    chebyshev1_biquads(order, 1.0, ch_w0, ch_Q, &ch_ns);
    printf("\nChebyshev Type I %dth-order (1dB ripple):\n", order);
    for (unsigned k = 0; k < ch_ns; k++)
        printf("  Section %u: w0=%.4f, Q=%.4f\n", k+1, ch_w0[k], ch_Q[k]);

    /* Bessel biquads */
    double be_w0[10], be_Q[10];
    unsigned be_ns;
    bessel_biquads(order, be_w0, be_Q, &be_ns);
    printf("\nBessel %dth-order (max flat delay):\n", order);
    for (unsigned k = 0; k < be_ns; k++)
        printf("  Section %u: w0=%.4f, Q=%.4f\n", k+1, be_w0[k], be_Q[k]);

    /* Sallen-Key design for Butterworth LP */
    double r, c, rf_by_rg;
    double gain = sallen_key_lp_design(1000.0, 0.7071, &r, &c, &rf_by_rg);
    printf("\nSallen-Key LP (Butterworth, fc=1kHz):\n");
    printf("  R = %.1f kOhm\n", r / 1000.0);
    printf("  C = %.1f nF\n", c * 1e9);
    printf("  Rf/Rg = %.3f (Gain = %.3f)\n", rf_by_rg, gain);

    /* Sallen-Key sensitivity */
    double s_r1, s_r2, s_c1, s_c2;
    sallen_key_sensitivity(r, r, c, c, &s_r1, &s_r2, &s_c1, &s_c2);
    printf("  Sensitivity: S_R=%.1f, S_C=%.1f\n", s_r1, s_c1);

    printf("\n=== Key trade-offs ===\n");
    printf("  Butterworth: flat passband, moderate transition\n");
    printf("  Chebyshev I: ripple in passband, faster roll-off\n");
    printf("  Bessel: flat group delay, gradual magnitude roll-off\n");
    printf("  Elliptic: ripple everywhere, steepest transition\n");

    return 0;
}
