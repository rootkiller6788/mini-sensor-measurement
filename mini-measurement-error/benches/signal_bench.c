#include "measurement_error.h"
#include "measurement_uncertainty.h"
#include "measurement_calibration.h"
#include "measurement_filtering.h"
#include <stdio.h>
#include <time.h>

#define ITER 100000
#define N_PTS 100

static double bench_time(const char *name, void (*fn)(void)) {
    clock_t start = clock();
    fn();
    clock_t end = clock();
    double ms = 1000.0 * (end - start) / CLOCKS_PER_SEC;
    printf("  %-30s %8.2f ms\n", name, ms);
    return ms;
}

/* Benchmark targets */
static void bm_absolute_error(void) {
    double sum = 0;
    for (int i = 0; i < ITER; i++)
        sum += meas_absolute_error((double)i, (double)(i-1));
    if (sum < 0) printf(""); /* prevent optimization */
}

static CalibrationPoint bench_pts[N_PTS];
static void bm_linear_fit(void) {
    CalibrationResult r;
    for (int i = 0; i < ITER/100; i++)
        cal_linear_fit(bench_pts, N_PTS, &r);
}

static void bm_moving_avg(void) {
    MovingAvgFilter f;
    mavg_init(&f, 32);
    double sum = 0;
    for (int i = 0; i < ITER; i++)
        sum += mavg_update(&f, (double)i);
    mavg_free(&f);
    if (sum < 0) printf("");
}

static void bm_gaussian_pdf(void) {
    double sum = 0;
    for (int i = 0; i < ITER; i++)
        sum += gaussian_pdf((double)i/ITER, 0.0, 1.0);
    if (sum < 0) printf("");
}

int main(void) {
    printf("=== Measurement Error Performance Benchmarks ===\n");
    printf("Iterations: %d\n\n", ITER);

    /* Setup calibration points */
    for (int i = 0; i < N_PTS; i++) {
        bench_pts[i].x = (double)i;
        bench_pts[i].y = 2.0 * (double)i + 1.0;
    }

    bench_time("absolute_error", bm_absolute_error);
    bench_time("linear_fit (100 pts)", bm_linear_fit);
    bench_time("moving_average (w=32)", bm_moving_avg);
    bench_time("gaussian_pdf", bm_gaussian_pdf);

    printf("\n=== Done ===\n");
    return 0;
}
