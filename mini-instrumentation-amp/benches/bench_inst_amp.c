/**
 * bench_inst_amp.c - Performance Benchmarks
 *
 * Microbenchmarks for core IA computation functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "inst_amp_defs.h"
#include "inst_amp_topology.h"
#include "inst_amp_analysis.h"
#include "inst_amp_calibration.h"

#define ITERATIONS 1000000
#define BENCH(name, expr) do { \
    clock_t start = clock(); \
    for (int _i = 0; _i < ITERATIONS; _i++) { volatile double _r = (expr); (void)_r; } \
    clock_t end = clock(); \
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC * 1e9 / ITERATIONS; \
    printf("  %-40s %8.1f ns/op\n", name, elapsed); \
} while(0)

int main(void) {
    printf("=== IA Library Performance Benchmarks ===\n");
    printf("  (lower is better, %d iterations each)\n\n", ITERATIONS);

    ia_three_opamp_config_t c3 = {25000, 10000, 10000, 5000, 0.1};
    BENCH("3-op-amp gain", ia_three_opamp_gain(&c3));
    BENCH("CMRR resistor limit", ia_three_opamp_cmrr_resistor_limit(&c3));
    BENCH("Thermal noise", noise_thermal_nv_per_rhz(1000.0, 27.0));
    BENCH("Noise figure", noise_figure_db(5.0, 1000.0, 27.0));
    BENCH("Error budget", error_budget_rti_uv(&(ia_spec_t){.vos_uv=50,.ios_pa=100,.cmrr_db=90,.psrr_db=80}, 1000.0, 5.0, 0.01));
    BENCH("E96 nearest", ia_e96_nearest(5555.56));
    BENCH("E24 nearest", ia_e24_nearest(4300.0));

    double coeffs[3] = {0};
    BENCH("Apply calibration", ia_cal_apply_correction(10.01, 10.0, 0.01));
    BENCH("Polynomial evaluation", ia_cal_eval_polynomial(coeffs, 2, 5.0));

    printf("\n=== Benchmarks Complete ===\n");
    return 0;
}
