#include "measurement_error.h"
#include "measurement_uncertainty.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define EPS 1e-10
#define NEAR(a,b) (fabs((a)-(b)) < EPS)

static int tests_run = 0;
static int tests_passed = 0;
#define TEST(name) do { tests_run++; printf("  %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)

int main(void) {
    printf("=== Measurement Error Core Tests ===\n");

    /* Test: Absolute error */
    TEST("absolute_error");
    assert(NEAR(meas_absolute_error(10.5, 10.0), 0.5));
    assert(NEAR(meas_absolute_error(9.8, 10.0), -0.2));
    PASS();

    /* Test: Relative error */
    TEST("relative_error");
    assert(NEAR(meas_relative_error(10.5, 10.0), 0.05));
    assert(NEAR(meas_relative_error(0.0, 1.0), -1.0));
    PASS();

    /* Test: Percentage error */
    TEST("percentage_error");
    assert(NEAR(meas_percentage_error(10.5, 10.0), 5.0));
    assert(NEAR(meas_percentage_error(5.0, 10.0), 50.0));
    PASS();

    /* Test: FS error */
    TEST("fs_error");
    assert(NEAR(meas_fs_error(0.5, 10.0), 5.0));
    PASS();

    /* Test: Accuracy class */
    TEST("accuracy_class");
    assert(NEAR(meas_accuracy_class(0.03), 0.05));
    assert(NEAR(meas_accuracy_class(0.3), 0.5));
    assert(NEAR(meas_accuracy_class(3.0), 5.0));
    PASS();

    /* Test: Repeatability */
    TEST("repeatability");
    double rep_data[] = {10.0, 10.1, 9.9, 10.0, 10.2, 9.8, 10.0, 10.1, 9.9, 10.0};
    double r = meas_repeatability(rep_data, 10);
    assert(r > 0.0);
    PASS();

    /* Test: Hysteresis error */
    TEST("hysteresis_error");
    double hyst = meas_hysteresis_error(5.0, 4.8, 10.0);
    assert(NEAR(hyst, 0.02));
    PASS();

    /* Test: Endpoint nonlinearity */
    TEST("endpoint_nonlinearity");
    double in[]  = {0.0, 2.5, 5.0, 7.5, 10.0};
    double out[] = {0.0, 2.6, 5.1, 7.4, 10.0};
    double max_dev;
    double nl = meas_endpoint_nonlinearity(in, out, 5, &max_dev);
    assert(nl >= 0.0);
    assert(max_dev >= 0.0);
    PASS();

    /* Test: Quantization error RMS */
    TEST("quantization_error");
    double q_rms = meas_quantization_error_rms(5.0, 12);
    assert(q_rms > 0.0);
    /* 12-bit ADC over 5V: LSB=5/4096=0.0012207, Q/sqrt(12)=0.0003524 */
    double expected = (5.0/4096.0)/sqrt(12.0);
    assert(NEAR(q_rms, expected));
    PASS();

    /* Test: Loading error */
    TEST("loading_error");
    double le = meas_loading_error_pct(1000.0, 1000000.0);
    double le_expected = 100.0 * 1000.0 / (1000.0 + 1000000.0);
    assert(NEAR(le, le_expected));
    PASS();

    /* Test: CMRR error */
    TEST("cmrr_error");
    double cm = meas_cmrr_error(5.0, 80.0);  /* 80dB CMRR */
    assert(cm > 0.0 && cm < 0.001);
    PASS();

    /* Test: Thermal EMF error */
    TEST("thermal_emf");
    double te = meas_thermal_emf_error(40.0, 10.0);
    assert(NEAR(te, 400.0));  /* 40 uV/C * 10C = 400 uV */
    PASS();

    /* Test: Error budget construction */
    TEST("error_budget");
    ErrorBudget budget;
    error_budget_init(&budget, 2.0);
    assert(budget.num_components == 0);
    int rc = error_budget_add(&budget, ERR_SYSTEMATIC_BIAS, 0.1, 1.0, 1.0, "Bias");
    assert(rc == 0);
    rc = error_budget_add(&budget, ERR_RANDOM_NOISE, 0.05, 1.0, 1.0, "Noise");
    assert(rc == 0);
    error_budget_finalize(&budget);
    assert(budget.total_systematic > 0.0);
    assert(budget.total_random_rss > 0.0);
    assert(budget.combined_uncertainty > 0.0);
    PASS();

    /* Test: Measurement statistics */
    TEST("meas_stats");
    MeasurementStats stats;
    meas_stats_init(&stats);
    double sdata[] = {9.8, 10.0, 10.2, 9.9, 10.1};
    meas_stats_update(&stats, sdata, 5, 10.0);
    assert(NEAR(stats.mean, 10.0));
    assert(stats.std_dev > 0.0);
    assert(NEAR(stats.bias, 0.0));
    PASS();

    /* Test: Category names */
    TEST("category_names");
    assert(strcmp(error_category_name(ERR_NONE), "None") == 0);
    assert(strcmp(error_category_name(ERR_QUANTIZATION), "Quantization") == 0);
    PASS();

    /* Test: Measurement type names */
    TEST("type_names");
    assert(strcmp(measurement_type_name(MEAS_DIRECT), "Direct") == 0);
    assert(strcmp(measurement_type_name(MEAS_RATIOMETRIC), "Ratiometric") == 0);
    PASS();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
