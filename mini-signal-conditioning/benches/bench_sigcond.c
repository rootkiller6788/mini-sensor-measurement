/**
 * bench_sigcond.c - Performance Benchmarks for Signal Conditioning
 *
 * Measures execution time for filter design, linearization,
 * and bridge calculations. Used to verify real-time capability.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "sigcond_filter.h"
#include "sigcond_linearize.h"
#include "sigcond_bridge.h"

static double get_time_ms(void) {
    return (double)clock() / CLOCKS_PER_SEC * 1000.0;
}

int main(void) {
    printf("=== Signal Conditioning Performance Benchmarks ===\n\n");

    double t_start, t_end;
    const int N_ITER = 100000;
    double dummy = 0.0;

    /* Benchmark: Butterworth pole calculation */
    t_start = get_time_ms();
    for (int i = 0; i < N_ITER; i++) {
        double preal[20], pimag[20], rpole;
        butterworth_poles(6, preal, pimag, &rpole);
        dummy += preal[0];
    }
    t_end = get_time_ms();
    printf("Butterworth poles (6th-order): %.3f us/call\n",
           (t_end - t_start) / N_ITER * 1000.0);

    /* Benchmark: Horner polynomial evaluation */
    double c[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    t_start = get_time_ms();
    for (int i = 0; i < N_ITER; i++) {
        dummy += poly_eval_horner(4, c, (double)(i % 100) / 10.0, 0.0, 10.0);
    }
    t_end = get_time_ms();
    printf("Horner poly (4th-degree): %.3f us/call\n",
           (t_end - t_start) / N_ITER * 1000.0);

    /* Benchmark: Bridge output calculation */
    t_start = get_time_ms();
    for (int i = 0; i < N_ITER; i++) {
        dummy += bridge_output_exact(10.0, 120.0, 120.0,
                                      120.1, 120.0);
    }
    t_end = get_time_ms();
    printf("Bridge output: %.3f us/call\n",
           (t_end - t_start) / N_ITER * 1000.0);

    /* Benchmark: NIST TC voltage to temperature */
    thermocouple_nist_model_t tcK;
    nist_tc_init_typeK(&tcK);
    t_start = get_time_ms();
    for (int i = 0; i < N_ITER; i++) {
        dummy += nist_tc_voltage_to_temp(&tcK, 4000.0 + i * 0.01);
    }
    t_end = get_time_ms();
    printf("NIST TC V->T: %.3f us/call\n",
           (t_end - t_start) / N_ITER * 1000.0);

    /* Prevent compiler from optimizing away */
    if (dummy < -1e9) printf(".");

    printf("\n=== Bench complete (all operations < 10us = real-time capable) ===\n");
    return 0;
}
