#include "measurement_calibration.h"
#include "measurement_filtering.h"
#include "measurement_statistics.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define EPS 1e-10
#define NEAR(a,b) (fabs((a)-(b)) < EPS)

static int tests_run = 0, tests_passed = 0;
#define TEST(n) do { tests_run++; printf("  %s... ", n); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)

int main(void) {
    printf("=== Calibration, Filtering & Statistics Tests ===\n");

    /* Test: Linear calibration */
    TEST("linear_fit");
    CalibrationPoint pts[] = {
        {0.0, 0.1}, {2.5, 2.6}, {5.0, 5.0}, {7.5, 7.6}, {10.0, 10.1}
    };
    CalibrationResult result;
    int rc = cal_linear_fit(pts, 5, &result);
    assert(rc == 0);
    assert(NEAR(result.coefficients[1], 1.0));
    assert(result.r_squared > 0.99);
    PASS();

    /* Test: Polynomial fit */
    TEST("polynomial_fit");
    CalibrationPoint quad[] = {
        {0.0, 0.0}, {1.0, 1.0}, {2.0, 4.0}, {3.0, 9.0}, {4.0, 16.0}
    };
    rc = cal_polynomial_fit(quad, 5, 2, &result);
    assert(rc == 0);
    assert(result.r_squared > 0.999);
    PASS();

    /* Test: Exponential fit */
    TEST("exponential_fit");
    CalibrationPoint exp_pts[] = {
        {0.0, 1.0}, {1.0, 2.718}, {2.0, 7.389}, {3.0, 20.085}
    };
    rc = cal_exponential_fit(exp_pts, 4, &result);
    assert(rc == 0);
    double y_pred = cal_evaluate(&result, 1.0);
    assert(fabs(y_pred - exp(1.0)) < 0.1);
    PASS();

    /* Test: Model evaluation */
    TEST("cal_evaluate");
    CalibrationPoint lin[] = {{0.0, 1.0}, {10.0, 21.0}};
    rc = cal_linear_fit(lin, 2, &result);
    assert(rc == 0);
    assert(NEAR(cal_evaluate(&result, 5.0), 11.0));
    PASS();

    /* Test: Inverse evaluation */
    TEST("cal_inverse");
    double x_est;
    rc = cal_inverse_evaluate(&result, 11.0, 5.0, &x_est);
    assert(rc == 0);
    assert(NEAR(x_est, 5.0));
    PASS();

    /* Test: R-squared computation */
    TEST("r_squared");
    double r2 = cal_compute_r_squared(lin, 2, &result);
    assert(NEAR(r2, 1.0));
    PASS();

    /* Test: Moving average filter */
    TEST("moving_avg");
    MovingAvgFilter mavg;
    rc = mavg_init(&mavg, 5);
    assert(rc == 0);
    double ma_out = 0.0;
    for (int i = 0; i < 10; i++) {
        ma_out = mavg_update(&mavg, 10.0);
    }
    assert(NEAR(ma_out, 10.0));
    mavg_free(&mavg);
    PASS();

    /* Test: EMA filter */
    TEST("ema_filter");
    EMAFilter ema;
    ema_init(&ema, 0.3);
    double ema_out = ema_update(&ema, 10.0);
    assert(NEAR(ema_out, 10.0));
    ema_out = ema_update(&ema, 10.0);
    assert(NEAR(ema_out, 10.0));
    PASS();

    /* Test: Median filter */
    TEST("median_filter");
    double md[] = {1.0, 5.0, 3.0, 2.0, 4.0};
    double med = median_filter(md, 5);
    assert(NEAR(med, 3.0));
    PASS();

    /* Test: Kalman filter 1D */
    TEST("kalman_1d");
    KalmanFilter1D kf;
    kf1d_init(&kf, 0.0, 1.0, 0.01, 0.1, 1.0);
    double kf_out = kf1d_update(&kf, 10.0);
    assert(kf_out > 0.0);
    kf1d_predict(&kf);
    kf_out = kf1d_update(&kf, 10.0);
    assert(fabs(kf_out - 10.0) < 1.0);
    PASS();

    /* Test: Sensor fusion */
    TEST("sensor_fusion");
    double readings[] = {10.0, 10.2, 9.8};
    double weights[] = {0.5, 0.3, 0.2};
    double fused = sensor_fuse_weighted_average(readings, weights, 3);
    double expected_fused = 5.0 + 3.06 + 1.96;
    assert(NEAR(fused, expected_fused));
    PASS();

    /* Test: Inverse-variance fusion */
    TEST("inv_variance_fusion");
    double vars[] = {0.01, 0.04, 0.09};
    double fused_var;
    double fused_iv = sensor_fuse_inverse_variance(readings, vars, 3, &fused_var);
    assert(fused_var < vars[0]);
    assert(fused_iv > 0.0);
    PASS();

    /* Test: Gaussian MLE */
    TEST("gaussian_fit");
    double gdata[] = {9.8, 10.0, 10.2, 9.9, 10.1, 10.0, 9.9, 10.1, 10.0, 10.2};
    double mu_hat, sigma_hat;
    stats_fit_gaussian(gdata, 10, &mu_hat, &sigma_hat);
    assert(NEAR(mu_hat, 10.02));
    assert(sigma_hat > 0.0);
    PASS();

    /* Test: Grubbs test */
    TEST("grubbs_test");
    double gdata2[] = {10.0, 10.1, 9.9, 10.0, 15.0, 10.2};
    int outlier_idx;
    double G = stats_grubbs_test(gdata2, 6, &outlier_idx);
    assert(G > 0.0);
    assert(outlier_idx == 4);
    PASS();

    /* Test: Sample size */
    TEST("sample_size");
    size_t n_req = stats_required_sample_size(1.0, 0.1, 0.05);
    assert(n_req > 100 && n_req < 1000);
    PASS();

    /* Test: Allan variance */
    TEST("allan_variance");
    double adata[] = {1.0, 1.1, 0.9, 1.0, 1.2, 0.8, 1.0, 1.1, 0.9, 1.0,
                      1.0, 1.1, 0.9, 1.0, 1.2, 0.8};
    double avar = stats_allan_variance(adata, 16, 2);
    assert(avar >= 0.0);
    PASS();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
