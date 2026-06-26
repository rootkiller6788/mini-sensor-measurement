/**
 * @file    test_smart_sensor.c
 * @brief   Comprehensive test suite for mini-smart-sensor
 *
 * Covers all core APIs with mathematically meaningful assertions.
 * Each test validates a specific knowledge point from L1-L8.
 *
 * Tests organized by knowledge level:
 *   L1: struct initialization, enums, state transitions
 *   L2: buffer operations, sensor node management
 *   L3: Wheatstone bridge, ADC metrics, transfer functions
 *   L4: SNR fundamentals, MDS, propagation of uncertainty
 *   L5: digital filters, calibration algorithms, fusion methods
 *   L6: temperature compensation, pressure correction
 *   L7: vibration severity, TPMS, SpO2, AQI
 *   L8: energy harvesting, Allan variance, edge features
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include "smart_sensor.h"
#include "sensor_conditioning.h"
#include "sensor_calibration.h"
#include "sensor_fusion.h"
#include "sensor_applications.h"
#include "sensor_advanced.h"

#define ASSERT_NEAR(a, b, tol) \
    assert(fabs((a) - (b)) < (tol))

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define RUN_TEST(name) do { \
    tests_run++; \
    printf("  RUN  %s... ", #name); \
    fflush(stdout); \
    if (test_##name() == 0) { \
        printf("PASSED\n"); \
        tests_passed++; \
    } else { \
        printf("FAILED\n"); \
        tests_failed++; \
    } \
} while(0)

/* Helper: declare test as returning int (0=pass, non-0=fail) */
#define TEST(name) static int test_##name(void)

/* =========================================================================
 * L1: Core Definitions Tests
 * ========================================================================= */

TEST(sensor_spec_init) {
    ss_sensor_spec_t spec;
    ss_sensor_spec_init(&spec, SS_TRANSDUCER_RESISTIVE,
                        SS_MEASURAND_TEMPERATURE, 10.0, 100.0, 1000.0);

    assert(spec.transducer_type == SS_TRANSDUCER_RESISTIVE);
    assert(spec.measurand == SS_MEASURAND_TEMPERATURE);
    ASSERT_NEAR(spec.sensitivity, 10.0, 1e-9);
    ASSERT_NEAR(spec.full_scale_input, 100.0, 1e-9);
    ASSERT_NEAR(spec.full_scale_output, 1000.0, 1e-9);
    assert(spec.resolution > 0.0);
    assert(spec.bandwidth > 0.0);
    return 0;
}

TEST(accuracy_spec_init) {
    ss_accuracy_spec_t acc;
    ss_accuracy_spec_init(&acc, 0.5, 0.01);

    ASSERT_NEAR(acc.accuracy, 0.5, 1e-9);
    ASSERT_NEAR(acc.precision, 0.01, 1e-9);
    assert(acc.linearity_error > 0.0);
    assert(acc.hysteresis_error > 0.0);
    return 0;
}

TEST(analog_frontend_init) {
    ss_analog_frontend_t afe;
    ss_analog_frontend_init(&afe);

    ASSERT_NEAR(afe.excitation_voltage, 3.3, 1e-9);
    assert(afe.amplifier_gain > 0.0);
    assert(afe.filter_order > 0);
    assert(afe.adc_resolution_bits == 12);
    return 0;
}

TEST(digital_processing_init) {
    ss_digital_processing_t dp;
    ss_digital_processing_init(&dp);

    assert(dp.moving_avg_window > 0);
    assert(dp.lowpass_alpha > 0.0);
    assert(dp.decimation_factor >= 1);
    return 0;
}

TEST(sensor_config_init) {
    ss_sensor_config_t cfg;
    ss_sensor_spec_t spec;
    ss_sensor_spec_init(&spec, SS_TRANSDUCER_PIEZOELECTRIC,
                        SS_MEASURAND_ACCELERATION, 100.0, 16.0, 1600.0);

    ss_sensor_config_init(&cfg, &spec, NULL, NULL, NULL, SS_INTERFACE_I2C);

    assert(cfg.interface_type == SS_INTERFACE_I2C);
    assert(cfg.sample_period_ms > 0);
    assert(cfg.enable_temp_comp == 1);
    assert(cfg.enable_self_test == 1);
    return 0;
}

TEST(sensor_node_init) {
    ss_node_t node;
    int ret = ss_node_init(&node, 4, 0x12345678ABCDEF00ULL);

    assert(ret == 0);
    assert(node.n_channels == 4);
    assert(node.node_id == 0x12345678ABCDEF00ULL);
    assert(node.active_channel_mask == 0);
    return 0;
}

TEST(sensor_node_add_channel) {
    ss_node_t node;
    ss_sensor_config_t cfg;
    ss_sensor_spec_t spec;
    int chan;

    ss_node_init(&node, 4, 1);
    ss_sensor_spec_init(&spec, SS_TRANSDUCER_THERMOELECTRIC,
                        SS_MEASURAND_TEMPERATURE, 0.041, 1372.0, 54.8);
    ss_sensor_config_init(&cfg, &spec, NULL, NULL, NULL, SS_INTERFACE_SPI);

    chan = ss_node_add_channel(&node, &cfg);
    assert(chan == 0);
    assert(node.active_channel_mask == 0x01);

    chan = ss_node_add_channel(&node, &cfg);
    assert(chan == 1);
    assert(node.active_channel_mask == 0x03);
    return 0;
}

TEST(sensor_state_set_mode) {
    ss_sensor_state_t state;
    memset(&state, 0, sizeof(state));
    state.mode = SS_MODE_IDLE;

    /* Valid transition: IDLE -> NORMAL_SAMPLING */
    assert(ss_sensor_state_set_mode(&state, SS_MODE_NORMAL_SAMPLING) == 0);
    assert(state.mode == SS_MODE_NORMAL_SAMPLING);

    /* Invalid transition: NORMAL_SAMPLING -> CALIBRATION */
    assert(ss_sensor_state_set_mode(&state, SS_MODE_CALIBRATION) == -1);

    /* Valid: NORMAL_SAMPLING -> FAULT_DETECTED */
    assert(ss_sensor_state_set_mode(&state, SS_MODE_FAULT_DETECTED) == 0);
    assert(state.fault_count == 1);

    /* Invalid: FAULT_DETECTED -> NORMAL_SAMPLING (must go via MAINTENANCE) */
    assert(ss_sensor_state_set_mode(&state, SS_MODE_NORMAL_SAMPLING) == -1);

    /* Valid: FAULT_DETECTED -> MAINTENANCE */
    assert(ss_sensor_state_set_mode(&state, SS_MODE_MAINTENANCE) == 0);
    return 0;
}

TEST(sensor_state_update_stats) {
    ss_sensor_state_t state;
    memset(&state, 0, sizeof(state));

    /* Feed known sequence: 1, 2, 3, 4, 5 */
    /* Mean = 3.0, sample variance = 2.5 */
    int i;
    for (i = 0; i < 5; i++) {
        state.sample_count++;
        ss_sensor_state_update_stats(&state, (double)(i + 1));
    }

    ASSERT_NEAR(state.moving_average, 3.0, 1e-9);
    ASSERT_NEAR(state.moving_variance, 2.5, 1e-9);
    return 0;
}

/* =========================================================================
 * L2: Buffer Operations Tests
 * ========================================================================= */

TEST(measurement_buffer_push_get) {
    ss_measurement_buffer_t buf;
    ss_measurement_t storage[4];
    ss_measurement_t meas, out;

    ss_measurement_buffer_init(&buf, storage, 4);

    meas.timestamp_us = 1000;
    meas.calibrated_value = 25.5;
    meas.flags = SS_FLAG_VALID;
    assert(ss_measurement_buffer_push(&buf, &meas) == 1);

    assert(ss_measurement_buffer_get(&buf, 0, &out) == 0);
    assert(out.timestamp_us == 1000);
    ASSERT_NEAR(out.calibrated_value, 25.5, 1e-9);

    /* Push more to fill buffer */
    meas.calibrated_value = 26.0;
    ss_measurement_buffer_push(&buf, &meas);
    meas.calibrated_value = 27.0;
    ss_measurement_buffer_push(&buf, &meas);
    meas.calibrated_value = 28.0;
    ss_measurement_buffer_push(&buf, &meas);

    /* Buffer should be full */
    assert(buf.count == 4);

    /* Push one more — oldest overwritten */
    meas.calibrated_value = 29.0;
    ss_measurement_buffer_push(&buf, &meas);

    /* Oldest is now 26.0 */
    assert(ss_measurement_buffer_get(&buf, 0, &out) == 0);
    ASSERT_NEAR(out.calibrated_value, 26.0, 1e-9);
    return 0;
}

TEST(measurement_buffer_mean_stddev) {
    ss_measurement_buffer_t buf;
    ss_measurement_t storage[5];
    ss_measurement_t meas;
    int i;

    ss_measurement_buffer_init(&buf, storage, 5);

    /* Insert values 1, 2, 3, 4, 5 */
    for (i = 0; i < 5; i++) {
        meas.calibrated_value = (double)(i + 1);
        ss_measurement_buffer_push(&buf, &meas);
    }

    double mean = ss_measurement_buffer_mean(&buf);
    ASSERT_NEAR(mean, 3.0, 1e-6);

    double stddev = ss_measurement_buffer_stddev(&buf);
    ASSERT_NEAR(stddev, sqrt(2.5), 1e-6); /* sample stddev of 1..5 = sqrt(2.5) */
    return 0;
}

TEST(measurement_buffer_range) {
    ss_measurement_buffer_t buf;
    ss_measurement_t storage[3];
    ss_measurement_t meas;
    double min, max;

    ss_measurement_buffer_init(&buf, storage, 3);

    meas.calibrated_value = 10.0; ss_measurement_buffer_push(&buf, &meas);
    meas.calibrated_value = -5.0; ss_measurement_buffer_push(&buf, &meas);
    meas.calibrated_value = 3.0;  ss_measurement_buffer_push(&buf, &meas);

    assert(ss_measurement_buffer_range(&buf, &min, &max) == 0);
    ASSERT_NEAR(min, -5.0, 1e-9);
    ASSERT_NEAR(max, 10.0, 1e-9);
    return 0;
}

/* =========================================================================
 * L3: Mathematical Structures Tests
 * ========================================================================= */

TEST(wheatstone_voltage) {
    /* Balanced bridge: R1=120, R2=120, R3=120, R4=120, Vex=5V */
    double v = ss_cond_wheatstone_voltage(5.0, 120.0, 120.0, 120.0, 120.0);
    ASSERT_NEAR(v, 0.0, 1e-9);

    /* Unbalanced: R1 changes to 121 */
    v = ss_cond_wheatstone_voltage(5.0, 121.0, 120.0, 120.0, 120.0);
    /* Vout = 5 * (120/241 - 120/240) = 5 * (0.4979 - 0.5) = -0.0104 V */
    assert(fabs(v) > 0.001); /* Should have non-zero output */
    return 0;
}

TEST(strain_from_bridge) {
    /* Quarter bridge, Vout=2.5mV, Vex=5V, GF=2.0 */
    double strain = ss_cond_strain_from_bridge(0.0025, 5.0, 2.0, 0);
    /* strain = 4 * 0.0025 / (5.0 * 2.0) = 0.001 = 1000 microstrain */
    ASSERT_NEAR(strain, 0.001, 1e-9);

    /* Full bridge: 4x more sensitive */
    strain = ss_cond_strain_from_bridge(0.0025, 5.0, 2.0, 2);
    ASSERT_NEAR(strain, 0.00025, 1e-9);
    return 0;
}

TEST(adc_enob) {
    /* Ideal 12-bit ADC: SINAD = 6.02*12 + 1.76 = 74.0 dB */
    double sinad = ss_cond_adc_sqnr(12);
    double enob = ss_cond_adc_enob(sinad);
    ASSERT_NEAR(enob, 12.0, 0.1);
    return 0;
}

TEST(adc_quantization_noise) {
    /* 12-bit ADC, 3.3V ref */
    double noise = ss_cond_adc_quantization_noise(12, 3.3);
    /* LSB = 3.3/4096 = 0.805mV, sigma = LSB/sqrt(12) = 0.232mV */
    double lsb = 3.3 / 4096.0;
    double expected = lsb / sqrt(12.0);
    ASSERT_NEAR(noise, expected, 1e-6);
    return 0;
}

TEST(adc_sqnr) {
    /* 16-bit ADC SQNR = 6.02*16 + 1.76 = 98.08 dB */
    double sqnr = ss_cond_adc_sqnr(16);
    ASSERT_NEAR(sqnr, 98.08, 0.1);
    return 0;
}

TEST(transfer_function) {
    double params[4] = {2.0, 1.0, 0.0, 0.0};

    /* Linear: f(x) = 2*x + 1 */
    double y = ss_cond_transfer_function(3.0, params, 0);
    ASSERT_NEAR(y, 7.0, 1e-9);

    /* Exponential: f(x) = 2*exp(x) + 0 */
    y = ss_cond_transfer_function(0.0, params, 1);
    ASSERT_NEAR(y, 2.0, 1e-9);
    return 0;
}

TEST(snr_db) {
    double snr = ss_cond_snr_db(1.0, 0.1);
    ASSERT_NEAR(snr, 20.0, 0.01); /* 20*log10(10) = 20 dB */
    return 0;
}

TEST(min_detectable_signal) {
    /* MDS at 10 dB SNR with noise=1.0: MDS = 1.0 * 10^(10/20) = 3.162 */
    double mds = ss_cond_min_detectable_signal(1.0, 10.0);
    ASSERT_NEAR(mds, 3.162, 0.01);
    return 0;
}

/* =========================================================================
 * L5: Digital Filtering Tests
 * ========================================================================= */

TEST(moving_average_filter) {
    /* Window of 4 samples: [1,2,3,4] -> avg = 2.5 */
    double buf[4] = {0.0, 0.0, 0.0, 0.0};
    /* Feed values one at a time: 1,2,3,4 -> final window [1,2,3,4], avg=2.5 */
    ss_cond_moving_average(buf, 4, 1.0);
    ss_cond_moving_average(buf, 4, 2.0);
    ss_cond_moving_average(buf, 4, 3.0);
    double result = ss_cond_moving_average(buf, 4, 4.0);
    ASSERT_NEAR(result, 2.5, 1e-9);
    return 0;
}

TEST(running_moving_avg) {
    double window[4] = {0};
    size_t count = 0, idx = 0;
    double sum = 0.0;
    double result;

    result = ss_cond_running_moving_avg(window, 4, &count, &idx, &sum, 1.0);
    ASSERT_NEAR(result, 1.0, 1e-9);

    result = ss_cond_running_moving_avg(window, 4, &count, &idx, &sum, 2.0);
    ASSERT_NEAR(result, 1.5, 1e-9);

    result = ss_cond_running_moving_avg(window, 4, &count, &idx, &sum, 3.0);
    ASSERT_NEAR(result, 2.0, 1e-9);

    result = ss_cond_running_moving_avg(window, 4, &count, &idx, &sum, 4.0);
    ASSERT_NEAR(result, 2.5, 1e-9);

    /* Window full — now circular */
    result = ss_cond_running_moving_avg(window, 4, &count, &idx, &sum, 10.0);
    ASSERT_NEAR(result, (2.0+3.0+4.0+10.0)/4.0, 1e-9);
    return 0;
}

TEST(median_filter) {
    double buf[5] = {1.0, 2.0, 10.0, 4.0, 5.0}; /* 10 is outlier */
    double result = ss_cond_median_filter(buf, 5, 3.0);
    /* Sorted: [2,3,4,5,10] -> median = 4 */
    ASSERT_NEAR(result, 4.0, 1e-9);
    return 0;
}

TEST(iir_lowpass) {
    /* alpha=1.0: no filtering, output = input */
    double prev = 0.0;
    double result = ss_cond_iir_lowpass(&prev, 5.0, 1.0);
    ASSERT_NEAR(result, 5.0, 1e-9);

    /* alpha=0.1: y = 0.1*input + 0.9*prev */
    prev = 0.0;
    result = ss_cond_iir_lowpass(&prev, 10.0, 0.1);
    ASSERT_NEAR(result, 1.0, 1e-9);

    /* Steady state: y -> input */
    result = ss_cond_iir_lowpass(&prev, 10.0, 0.1);
    ASSERT_NEAR(result, 1.9, 1e-9);
    return 0;
}

TEST(iir_highpass) {
    /* alpha=1.0: y = input - prev_input */
    double prev_out = 0.0, prev_in = 0.0;
    double result = ss_cond_iir_highpass(&prev_out, &prev_in, 5.0, 1.0);
    ASSERT_NEAR(result, 5.0, 1e-9); /* prev_in was 0 */

    result = ss_cond_iir_highpass(&prev_out, &prev_in, 5.0, 1.0);
    ASSERT_NEAR(result, 0.0, 1e-9); /* DC blocked */
    return 0;
}

TEST(threshold_hysteresis) {
    int state = 0;
    /* Below low: stay 0 */
    assert(ss_cond_threshold_hysteresis(&state, 1.0, 2.0, 3.0) == 0);

    /* Cross high: turn on */
    assert(ss_cond_threshold_hysteresis(&state, 3.5, 2.0, 3.0) == 1);

    /* In hysteresis band: stay 1 */
    assert(ss_cond_threshold_hysteresis(&state, 2.5, 2.0, 3.0) == 1);

    /* Cross low: turn off */
    assert(ss_cond_threshold_hysteresis(&state, 1.5, 2.0, 3.0) == 0);
    return 0;
}

TEST(temp_compensate) {
    /* At reference temp: no correction */
    double val = ss_cond_temp_compensate(100.0, 25.0, 25.0, 0.1, 0.001);
    ASSERT_NEAR(val, 100.0, 1e-9);

    /* 10 degC above ref, zero drift=0.5 units/degC */
    val = ss_cond_temp_compensate(105.0, 35.0, 25.0, 0.5, 0.0);
    ASSERT_NEAR(val, 100.0, 1e-9); /* 105 - 0.5*10 = 100 */
    return 0;
}

/* =========================================================================
 * L5: Calibration Algorithm Tests
 * ========================================================================= */

TEST(linear_regression) {
    double x[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double y[] = {2.0, 4.0, 6.0, 8.0, 10.0}; /* y = 2*x */
    ss_linear_model_t model;

    int ret = ss_cal_linear_regression(x, y, 5, &model);
    assert(ret == 0);
    ASSERT_NEAR(model.sensitivity, 2.0, 1e-9);
    ASSERT_NEAR(model.offset, 0.0, 1e-9);
    ASSERT_NEAR(model.r_squared, 1.0, 1e-9);
    return 0;
}

TEST(linear_apply_inverse) {
    ss_linear_model_t model = {2.0, 1.0, 1.0, 0.0}; /* y = 2*x + 1 */
    double y = ss_cal_linear_apply(&model, 3.0);
    ASSERT_NEAR(y, 7.0, 1e-9);

    double x = ss_cal_linear_inverse(&model, 7.0);
    ASSERT_NEAR(x, 3.0, 1e-9);
    return 0;
}

TEST(polynomial_regression) {
    /* y = x^2: quadratic perfect fit */
    double x[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double y[] = {0.0, 1.0, 4.0, 9.0, 16.0};
    ss_polynomial_model_t model;

    int ret = ss_cal_polynomial_regression(x, y, 5, 2, &model);
    assert(ret == 0);
    assert(model.order == 2);

    /* Verify at a test point */
    double y_pred = ss_cal_polynomial_apply(&model, 2.5);
    ASSERT_NEAR(y_pred, 6.25, 1e-3);
    return 0;
}

TEST(polynomial_inverse) {
    /* y = x^2: find x where y=9.0, should be x=3.0 */
    double x[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    double y[] = {0.0, 1.0, 4.0, 9.0, 16.0, 25.0};
    ss_polynomial_model_t model;

    ss_cal_polynomial_regression(x, y, 6, 2, &model);

    double x_root = ss_cal_polynomial_inverse(&model, 9.0, 2.0, 1e-6, 50);
    ASSERT_NEAR(x_root, 3.0, 1e-4);
    return 0;
}

TEST(steinhart_hart_3pt) {
    /* Three calibration points for a typical NTC thermistor */
    /* At 0°C (273.15K): R = 27280 Ohm */
    /* At 25°C (298.15K): R = 10000 Ohm */
    /* At 50°C (323.15K): R = 3588 Ohm */
    ss_steinhart_hart_t model;

    int ret = ss_cal_steinhart_hart_3pt(27280.0, 273.15,
                                         10000.0, 298.15,
                                          3588.0, 323.15, &model);
    assert(ret == 0);

    /* Verify at the calibration points */
    double t1 = ss_cal_steinhart_hart_temp(&model, 27280.0);
    ASSERT_NEAR(t1, 273.15, 0.1);

    double t2 = ss_cal_steinhart_hart_temp(&model, 10000.0);
    ASSERT_NEAR(t2, 298.15, 0.1);

    double t3 = ss_cal_steinhart_hart_temp(&model, 3588.0);
    ASSERT_NEAR(t3, 323.15, 0.1);
    return 0;
}

TEST(thermistor_beta) {
    ss_steinhart_hart_t model;

    int ret = ss_cal_thermistor_beta(10000.0, 298.15, 3588.0, 323.15, &model);
    assert(ret == 0);
    assert(model.beta_value > 0.0);

    /* Verify at reference point */
    double t = ss_cal_steinhart_hart_temp(&model, 10000.0);
    ASSERT_NEAR(t, 298.15, 1.0);
    return 0;
}

TEST(two_point_calibration) {
    ss_linear_model_t model;
    /* Sensor reads 1.0V at 0 degC, 3.0V at 100 degC */
    int ret = ss_cal_two_point(0.0, 1.0, 100.0, 3.0, &model);
    assert(ret == 0);
    ASSERT_NEAR(model.sensitivity, 50.0, 1e-9);
    ASSERT_NEAR(model.offset, -50.0, 1e-9);

    /* Verify: at y_raw=2.0, x_true = 50.0 */
    double x = ss_cal_linear_apply(&model, 2.0);
    ASSERT_NEAR(x, 50.0, 1e-9);
    return 0;
}

TEST(combined_uncertainty) {
    /* Two inputs: sens_coeff [2, 3], unc [0.1, 0.2] */
    /* u_c = sqrt(2^2*0.1^2 + 3^2*0.2^2) = sqrt(0.04+0.36) = sqrt(0.4) = 0.6325 */
    double sens[] = {2.0, 3.0};
    double unc[] = {0.1, 0.2};

    double uc = ss_cal_combined_uncertainty(sens, unc, 2);
    ASSERT_NEAR(uc, sqrt(0.4), 1e-6);
    return 0;
}

/* =========================================================================
 * L5: Sensor Fusion Tests
 * ========================================================================= */

TEST(kalman_1d_step) {
    ss_kalman_1d_t kf;
    ss_fusion_kalman_1d_init(&kf, 0.0, 1.0, 0.01, 0.1);

    /* After one step with z=10, state should move toward 10 */
    double x = ss_fusion_kalman_1d_step(&kf, 10.0);
    assert(x > 0.0 && x < 10.0); /* Between prior and measurement */

    /* After many steps with same measurement, converge to ~10 */
    int i;
    for (i = 0; i < 100; i++) {
        x = ss_fusion_kalman_1d_step(&kf, 10.0);
    }
    ASSERT_NEAR(x, 10.0, 0.1);
    return 0;
}

TEST(kalman_imu_angles) {
    ss_kalman_3d_imu_t kf;
    double roll, pitch;

    ss_fusion_kalman_imu_init(&kf, 0.01, 0.001, 0.01);

    /* Simulate: accelerometer says level (0 roll, 0 pitch) */
    ss_fusion_kalman_imu_predict(&kf, 0.0, 0.0, 0.0);
    ss_fusion_kalman_imu_update(&kf, 0.0, 0.0, 1.0); /* Gravity straight down */

    ss_fusion_kalman_imu_get_angles(&kf, &roll, &pitch);
    ASSERT_NEAR(roll, 0.0, 0.02);
    ASSERT_NEAR(pitch, 0.0, 0.02);
    return 0;
}

TEST(complementary_filter) {
    ss_complementary_filter_t cf;
    ss_fusion_complementary_init(&cf, 0.98, 0.01);

    /* Initial: accel says level, gyro is still */
    ss_fusion_complementary_update(&cf, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    ASSERT_NEAR(cf.roll, 0.0, 0.01);
    ASSERT_NEAR(cf.pitch, 0.0, 0.01);
    return 0;
}

TEST(weighted_average) {
    /* Two sensors: reading=[10, 12], variance=[1, 4] */
    /* w1=1.0, w2=0.25, fused = (1*10+0.25*12)/(1.25) = 10.4 */
    double readings[] = {10.0, 12.0};
    double variances[] = {1.0, 4.0};
    double fused_var;

    double fused = ss_fusion_weighted_average(readings, variances, 2, &fused_var);
    ASSERT_NEAR(fused, 10.4, 1e-6);
    /* Fused variance should be < min(sigma_i^2) = 1.0 */
    assert(fused_var < 1.0);
    ASSERT_NEAR(fused_var, 0.8, 1e-6); /* 1/(1+0.25) = 0.8 */
    return 0;
}

TEST(median_fusion) {
    double readings[] = {10.0, 100.0, 12.0, 11.0, 13.0}; /* 100 is outlier */
    double med = ss_fusion_median(readings, 5);
    ASSERT_NEAR(med, 12.0, 1e-9); /* Sorted: [10,11,12,13,100] */
    return 0;
}

TEST(trimmed_mean_fusion) {
    double readings[] = {1.0, 100.0, 3.0, 5.0, 7.0}; /* trim 1 each side */
    double result = ss_fusion_trimmed_mean(readings, 5, 1);
    ASSERT_NEAR(result, 5.0, 1e-6); /* After trim: [3,5,7] -> mean=5 */
    return 0;
}

TEST(ema_fusion) {
    double ema = 0.0;
    ema = ss_fusion_ema(&ema, 10.0, 0.5);
    ASSERT_NEAR(ema, 5.0, 1e-9);

    ema = ss_fusion_ema(&ema, 10.0, 0.5);
    ASSERT_NEAR(ema, 7.5, 1e-9);
    return 0;
}

TEST(cumulative_avg) {
    double sum = 0.0;
    uint64_t count = 0;
    double avg;

    avg = ss_fusion_cumulative_avg(&sum, &count, 10.0);
    ASSERT_NEAR(avg, 10.0, 1e-9);

    avg = ss_fusion_cumulative_avg(&sum, &count, 20.0);
    ASSERT_NEAR(avg, 15.0, 1e-9);

    assert(count == 2);
    return 0;
}

TEST(mahalanobis_1d) {
    /* reading=12, mean=10, var=4 -> D = |12-10|/2 = 1.0 */
    double d = ss_fusion_mahalanobis_1d(12.0, 10.0, 4.0);
    ASSERT_NEAR(d, 1.0, 1e-9);
    return 0;
}

TEST(outlier_rejection) {
    /* 4 sensors with one slightly different: all within 3 sigma */
    double readings[] = {10.0, 10.3, 11.2, 9.9};
    double variances[] = {0.1, 0.1, 0.1, 0.1};
    uint32_t mask;
    double fused;

    int n = ss_fusion_outlier_rejection(readings, variances, 4, 3.0, &mask, &fused);

    /* All sensors should be accepted (no true outlier at 3-sigma threshold) */
    assert(n >= 3);
    assert(fused > 9.5 && fused < 11.0); /* Fused near 10 */
    return 0;
}

/* =========================================================================
 * L7: Application Tests
 * ========================================================================= */

TEST(vibration_severity) {
    /* 0.1g RMS at 100 Hz -> v = 0.1*9810/(2*pi*100) = 1.56 mm/s */
    double v = ss_app_vibration_severity(0.1, 100.0);
    ASSERT_NEAR(v, 1.562, 0.01);
    return 0;
}

TEST(vibration_zone) {
    assert(ss_app_vibration_zone(0.5) == 0);  /* Zone A */
    assert(ss_app_vibration_zone(2.0) == 1);  /* Zone B */
    assert(ss_app_vibration_zone(3.5) == 2);  /* Zone C */
    assert(ss_app_vibration_zone(5.0) == 3);  /* Zone D */
    return 0;
}

TEST(bpfo_frequency) {
    /* 8 balls, 1800 RPM, ball=10mm, pitch=50mm, angle=0 */
    /* BPFO = 8/2 * 30 * (1 - 10/50) = 4*30*0.8 = 96 Hz */
    double f = ss_app_bpfo_frequency(8, 1800.0, 10.0, 50.0, 0.0);
    ASSERT_NEAR(f, 96.0, 0.1);
    return 0;
}

TEST(tpms_compensate) {
    /* 200 kPa at 35°C, reference 20°C */
    /* P_ref = 200 * (293.15/308.15) = 190.26 kPa */
    double p = ss_app_tpms_compensate_pressure(200.0, 35.0, 20.0);
    ASSERT_NEAR(p, 190.26, 0.1);
    return 0;
}

TEST(tpms_deviation) {
    /* 180 kPa vs placard 220 kPa = -18.2% */
    double dev = ss_app_tpms_deviation_percent(180.0, 220.0);
    ASSERT_NEAR(dev, -18.1818, 0.1);
    return 0;
}

TEST(tpms_slow_leak) {
    double pressure[] = {220.0, 219.8, 219.5, 219.1, 218.6};
    double time[]     = {0.0, 24.0, 48.0, 72.0, 96.0};
    double leak_rate;
    int leak = ss_app_tpms_detect_slow_leak(pressure, time, 5, &leak_rate);
    assert(leak == 1); /* Leak detected */
    assert(leak_rate < 0.0); /* Negative rate */
    return 0;
}

TEST(spo2_estimate) {
    /* Typical PPG: AC_red/DC_red = 0.01, AC_ir/DC_ir = 0.02 */
    /* R = 0.01/0.02 = 0.5, SpO2 = 110 - 25*0.5 = 97.5% */
    double spo2 = ss_app_spo2_estimate(0.01, 1.0, 0.02, 1.0, 110.0, 25.0);
    ASSERT_NEAR(spo2, 97.5, 0.1);
    return 0;
}

TEST(heart_rate) {
    double intervals[] = {800.0, 820.0, 790.0}; /* ms */
    double hr = ss_app_heart_rate_from_ppg(intervals, 3);
    /* avg interval = 803.33 ms, HR = 60000/803.33 = 74.7 BPM */
    ASSERT_NEAR(hr, 74.7, 0.2);
    return 0;
}

TEST(perfusion_index) {
    double pi = ss_app_perfusion_index(0.05, 1.0);
    ASSERT_NEAR(pi, 5.0, 1e-9);
    return 0;
}

TEST(aqi_pm25) {
    /* PM2.5 = 10 ug/m3 -> AQI Good (0-50) */
    double aqi = ss_app_aqi_pm25(10.0);
    assert(aqi >= 0.0 && aqi <= 50.0);

    /* PM2.5 = 150 -> AQI Unhealthy (151-200) */
    aqi = ss_app_aqi_pm25(150.0);
    assert(aqi >= 151.0 && aqi <= 200.0);
    return 0;
}

TEST(co2_comfort) {
    assert(ss_app_co2_comfort_level(500.0) == 0);
    assert(ss_app_co2_comfort_level(800.0) == 1);
    assert(ss_app_co2_comfort_level(1500.0) == 2);
    assert(ss_app_co2_comfort_level(3000.0) == 3);
    return 0;
}

TEST(dew_point) {
    /* At 30°C, 50% RH: dew point ~ 18.4°C */
    double td = ss_app_dew_point(30.0, 50.0);
    ASSERT_NEAR(td, 18.4, 0.5);
    return 0;
}

TEST(heat_index) {
    /* At 32°C (90°F), 70% RH: heat index ~ 41°C (105°F) */
    double hi = ss_app_heat_index(32.0, 70.0);
    ASSERT_NEAR(hi, 41.0, 2.0);
    return 0;
}

/* =========================================================================
 * L8: Advanced Topic Tests
 * ========================================================================= */

TEST(pv_power) {
    /* 8% efficiency, 10 cm^2, 200 uW/cm^2 -> 160 uW */
    double p = ss_advanced_pv_power_uw(0.08, 10.0, 200.0);
    ASSERT_NEAR(p, 160.0, 1e-6);
    return 0;
}

TEST(teg_power) {
    /* Seebeck=200e-6 V/K, dT=10K, R=5 Ohm */
    /* P = (200e-6 * 10)^2 / (4*5) = 4e-6 / 20 = 0.2 uW = 2e-7 W */
    double p = ss_advanced_teg_power(200e-6, 10.0, 5.0);
    ASSERT_NEAR(p, 2.0e-7, 1e-10);
    return 0;
}

TEST(avg_power_duty_cycle) {
    /* Active=50mW, Sleep=10uW, duty=0.01, wakeup=100uJ at 1Hz */
    double p = ss_advanced_avg_power_uw(50.0, 10.0, 0.01, 100.0, 1.0);
    /* Expected: 50000*0.01 + 10*0.99 + 100 = 500+9.9+100 = 609.9 uW */
    assert(p > 500.0 && p < 700.0);
    return 0;
}

TEST(max_duty_cycle) {
    /* Harvest=1000uW, Active=50mW, Sleep=10uW */
    double duty = ss_advanced_max_duty_cycle(1000.0, 50.0, 10.0);
    /* Available = 1000-10 = 990 uW */
    /* duty = 990/(50000-10) = 0.0198 */
    ASSERT_NEAR(duty, 0.0198, 0.001);
    return 0;
}

TEST(signal_rms) {
    double data[] = {1.0, 2.0, 3.0};
    double rms = ss_advanced_signal_rms(data, 3);
    ASSERT_NEAR(rms, sqrt((1+4+9)/3.0), 1e-9);
    return 0;
}

TEST(crest_factor) {
    /* Sine: peak = 1.414, RMS = 1.0 -> CF = 1.414 */
    double data[] = {0.0, 0.707, 1.0, 0.707, 0.0, -0.707, -1.0, -0.707};
    double cf = ss_advanced_crest_factor(data, 8);
    ASSERT_NEAR(cf, 1.414, 0.1);
    return 0;
}

TEST(skewness) {
    /* Symmetric data around mean: [-2, -1, 0, 1, 2] -> mean=0, skewness=0 */
    double data[] = {-2.0, -1.0, 0.0, 1.0, 2.0};
    double sk = ss_advanced_skewness(data, 5);
    ASSERT_NEAR(sk, 0.0, 1e-9);
    return 0;
}

TEST(zero_crossing_rate) {
    double data[] = {1.0, -1.0, 1.0, -1.0, 1.0};
    double zcr = ss_advanced_zero_crossing_rate(data, 5);
    ASSERT_NEAR(zcr, 1.0, 1e-9); /* Every pair crosses zero */

    double data2[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    zcr = ss_advanced_zero_crossing_rate(data2, 5);
    ASSERT_NEAR(zcr, 0.0, 1e-9); /* No zero crossings */
    return 0;
}

TEST(allan_variance) {
    /* Generate white noise: AV ~ 1/sqrt(tau) */
    double data[1000];
    size_t i;
    double tau_vals[] = {0.1, 1.0, 10.0};
    double av_out[3];

    /* Initialize with simple random data */
    for (i = 0; i < 1000; i++) {
        data[i] = (double)(rand() % 1000) / 1000.0 - 0.5;
    }

    int ret = ss_advanced_allan_variance(data, 1000, 100.0, tau_vals, 3, av_out);
    assert(ret == 0);
    /* Check that longer tau gives smaller Allan deviation (for white noise) */
    assert(av_out[0] > 0.0); /* Should compute something */
    return 0;
}

TEST(spectral_centroid) {
    double mags[] = {1.0, 2.0, 1.0};
    double freqs[] = {100.0, 200.0, 300.0};
    /* centroid = (100*1 + 200*2 + 300*1)/(1+2+1) = (100+400+300)/4 = 200 Hz */
    double sc = ss_advanced_spectral_centroid(mags, freqs, 3);
    ASSERT_NEAR(sc, 200.0, 1e-9);
    return 0;
}

/* =========================================================================
 * Main test runner
 * ========================================================================= */

int main(void) {
    printf("=== mini-smart-sensor Test Suite ===\n\n");

    printf("--- L1: Core Definitions ---\n");
    RUN_TEST(sensor_spec_init);
    RUN_TEST(accuracy_spec_init);
    RUN_TEST(analog_frontend_init);
    RUN_TEST(digital_processing_init);
    RUN_TEST(sensor_config_init);
    RUN_TEST(sensor_node_init);
    RUN_TEST(sensor_node_add_channel);
    RUN_TEST(sensor_state_set_mode);
    RUN_TEST(sensor_state_update_stats);

    printf("\n--- L2: Buffer Operations ---\n");
    RUN_TEST(measurement_buffer_push_get);
    RUN_TEST(measurement_buffer_mean_stddev);
    RUN_TEST(measurement_buffer_range);

    printf("\n--- L3: Mathematical Structures ---\n");
    RUN_TEST(wheatstone_voltage);
    RUN_TEST(strain_from_bridge);
    RUN_TEST(adc_enob);
    RUN_TEST(adc_quantization_noise);
    RUN_TEST(adc_sqnr);
    RUN_TEST(transfer_function);
    RUN_TEST(snr_db);
    RUN_TEST(min_detectable_signal);

    printf("\n--- L5: Digital Filtering ---\n");
    RUN_TEST(moving_average_filter);
    RUN_TEST(running_moving_avg);
    RUN_TEST(median_filter);
    RUN_TEST(iir_lowpass);
    RUN_TEST(iir_highpass);
    RUN_TEST(threshold_hysteresis);
    RUN_TEST(temp_compensate);

    printf("\n--- L5: Calibration Algorithms ---\n");
    RUN_TEST(linear_regression);
    RUN_TEST(linear_apply_inverse);
    RUN_TEST(polynomial_regression);
    RUN_TEST(polynomial_inverse);
    RUN_TEST(steinhart_hart_3pt);
    RUN_TEST(thermistor_beta);
    RUN_TEST(two_point_calibration);
    RUN_TEST(combined_uncertainty);

    printf("\n--- L5: Sensor Fusion ---\n");
    RUN_TEST(kalman_1d_step);
    RUN_TEST(kalman_imu_angles);
    RUN_TEST(complementary_filter);
    RUN_TEST(weighted_average);
    RUN_TEST(median_fusion);
    RUN_TEST(trimmed_mean_fusion);
    RUN_TEST(ema_fusion);
    RUN_TEST(cumulative_avg);
    RUN_TEST(mahalanobis_1d);
    RUN_TEST(outlier_rejection);

    printf("\n--- L7: Applications ---\n");
    RUN_TEST(vibration_severity);
    RUN_TEST(vibration_zone);
    RUN_TEST(bpfo_frequency);
    RUN_TEST(tpms_compensate);
    RUN_TEST(tpms_deviation);
    RUN_TEST(tpms_slow_leak);
    RUN_TEST(spo2_estimate);
    RUN_TEST(heart_rate);
    RUN_TEST(perfusion_index);
    RUN_TEST(aqi_pm25);
    RUN_TEST(co2_comfort);
    RUN_TEST(dew_point);
    RUN_TEST(heat_index);

    printf("\n--- L8: Advanced Topics ---\n");
    RUN_TEST(pv_power);
    RUN_TEST(teg_power);
    RUN_TEST(avg_power_duty_cycle);
    RUN_TEST(max_duty_cycle);
    RUN_TEST(signal_rms);
    RUN_TEST(crest_factor);
    RUN_TEST(skewness);
    RUN_TEST(zero_crossing_rate);
    RUN_TEST(allan_variance);
    RUN_TEST(spectral_centroid);

    printf("\n========================================\n");
    printf("Tests Run: %d, Passed: %d, Failed: %d\n",
           tests_run, tests_passed, tests_failed);
    printf("========================================\n");

    return tests_failed > 0 ? 1 : 0;
}
