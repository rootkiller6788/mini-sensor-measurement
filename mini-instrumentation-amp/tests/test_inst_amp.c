/**
 * test_inst_amp.c - Comprehensive test suite for instrumentation amplifier library
 *
 * Tests all core functions with assert-based verification.
 * Run: make test
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "inst_amp_defs.h"
#include "inst_amp_topology.h"
#include "inst_amp_analysis.h"
#include "inst_amp_calibration.h"
#include "inst_amp_sensorif.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { printf("  %-50s", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define CHECK(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)
#define CHECK_CLOSE(a, b, tol, msg) do { if (fabs((a)-(b)) > (tol)) { printf("FAIL: %s (%.6f vs %.6f)\n", msg, a, b); tests_failed++; return; } } while(0)

/* ---- L2: 3-op-amp topology ---- */
static void test_three_opamp_gain(void) {
    TEST("3-op-amp gain at G=10");
    ia_three_opamp_config_t cfg = {25000.0, 10000.0, 10000.0, 5555.56, 0.1};
    double G = ia_three_opamp_gain(&cfg);
    CHECK_CLOSE(G, 10.0, 0.1, "G should be ~10");
    PASS();
}

static void test_three_opamp_rg_selection(void) {
    TEST("R_G selection for target gain");
    ia_three_opamp_config_t cfg = {25000.0, 10000.0, 10000.0, 1000.0, 0.1};
    double rg = ia_three_opamp_rg_for_gain(&cfg, 51.0);
    CHECK_CLOSE(rg, 1000.0, 10.0, "R_G should be ~1000 for G=51");
    PASS();
}

static void test_three_opamp_cmrr_resistor(void) {
    TEST("CMRR resistor mismatch limit");
    ia_three_opamp_config_t cfg = {25000.0, 10000.0, 10000.0, 1000.0, 0.1};
    double cmrr = ia_three_opamp_cmrr_resistor_limit(&cfg);
    CHECK(cmrr > 40.0 && cmrr < 80.0, "CMRR should be 40-80 dB for 0.1%");
    PASS();
}

static void test_three_opamp_output(void) {
    TEST("3-op-amp output voltage");
    ia_three_opamp_config_t cfg = {25000.0, 10000.0, 10000.0, 5555.56, 0.1};
    double vout = ia_three_opamp_output(&cfg, 0.051, -0.049, 0.0, 90.0);
    CHECK_CLOSE(vout, 1.0, 0.02, "Vout should be ~1V for 100mV diff at G=10");
    PASS();
}

/* ---- L2: 2-op-amp ---- */
static void test_two_opamp_gain(void) {
    TEST("2-op-amp gain");
    ia_two_opamp_config_t cfg = {10000.0, 10000.0, 10000.0, 10000.0, 5000.0};
    double G = ia_two_opamp_gain(&cfg);
    CHECK(G > 1.0, "Gain should be > 1");
    PASS();
}

static void test_two_opamp_cmrr_zero(void) {
    TEST("2-op-amp CMRR zero frequency");
    ia_two_opamp_config_t cfg = {10000.0, 10000.0, 10000.0, 10000.0, 5000.0};
    double fz = ia_two_opamp_cmrr_zero_freq(&cfg, 1.0);
    CHECK(fz > 0.0, "CMRR zero freq should be positive");
    PASS();
}

/* ---- L2: Current-mode ---- */
static void test_current_mode_gain(void) {
    TEST("Current-mode IA gain");
    ia_current_mode_config_t cfg = {1000.0, 10000.0, 1.0};
    double G = ia_current_mode_gain(&cfg);
    CHECK_CLOSE(G, 10.0, 0.01, "Gain should be 10");
    PASS();
}

/* ---- L2: Flying-capacitor ---- */
static void test_flying_cap_gain(void) {
    TEST("Flying-capacitor IA gain");
    ia_flying_cap_config_t cfg = {10.0, 1.0, 100.0, 5.0};
    double G = ia_flying_cap_gain(&cfg);
    CHECK_CLOSE(G, 10.0, 0.01, "Gain should be Cin/Cfb = 10");
    PASS();
}

/* ---- L2: Indirect current ---- */
static void test_indirect_current_gain(void) {
    TEST("Indirect current feedback gain");
    ia_indirect_current_config_t cfg = {1.0, 10000.0};
    double G = ia_indirect_current_gain(&cfg);
    CHECK_CLOSE(G, 10.0, 0.1, "Gain should be Gm*Rf = 10");
    PASS();
}

/* ---- L3: Noise analysis ---- */
static void test_thermal_noise(void) {
    TEST("Johnson-Nyquist thermal noise");
    double en = noise_thermal_nv_per_rhz(1000.0, 27.0);
    CHECK(en > 3.0 && en < 5.0, "1000 Ohm at 27C should be ~4.07 nV/rtHz");
    PASS();
}

static void test_shot_noise(void) {
    TEST("Shot noise");
    double in = noise_shot_fa_per_rhz(1e-9);
    CHECK(in > 0.0, "Shot noise should be positive");
    PASS();
}

static void test_noise_figure(void) {
    TEST("Noise figure");
    double nf = noise_figure_db(5.0, 1000.0, 27.0);
    CHECK(nf > 0.0, "NF should be > 0 dB for noisy amp");
    PASS();
}

/* ---- L3: Error budget ---- */
static void test_error_budget(void) {
    TEST("DC error budget RSS");
    ia_spec_t spec;
    ia_spec_init(&spec);
    double err = error_budget_rti_uv(&spec, 1000.0, 5.0, 0.01);
    CHECK(err > 0.0, "Error budget should be positive");
    PASS();
}

static void test_gain_error_budget(void) {
    TEST("Gain error budget");
    double err = gain_error_budget_pct(1.0, 100.0, 100000.0);
    CHECK(err > 0.0, "Gain error should be positive");
    PASS();
}

/* ---- L3: CMRR analysis ---- */
static void test_cmrr_combined(void) {
    TEST("Combined CMRR");
    double cmrr = cmrr_combined_db(100.0, 60.0, 120.0);
    CHECK(cmrr < 60.0, "Combined CMRR should be limited by worst contributor");
    PASS();
}

/* ---- L1: Spec validation ---- */
static void test_spec_validate(void) {
    TEST("IA spec validation");
    ia_spec_t spec;
    ia_spec_init(&spec);
    ia_error_code_t ec = ia_spec_validate(&spec);
    CHECK(ec == IA_OK, "Default spec should be valid");
    PASS();
}

/* ---- L1: Topology recommendation ---- */
static void test_topology_recommend(void) {
    TEST("Topology recommendation engine");
    ia_topology_t topo = ia_recommend_topology(SENSOR_STRAIN_GAUGE, 5.0);
    CHECK(topo == IA_TOPO_THREE_OPAMP, "Strain gauge should use 3-op-amp");
    PASS();
}

/* ---- L5: Calibration ---- */
static void test_two_point_cal(void) {
    TEST("Two-point calibration");
    double gain, offset;
    ia_cal_two_point(0.0, 0.01, 1.0, 10.01, &gain, &offset);
    CHECK_CLOSE(gain, 10.0, 0.01, "Gain should be ~10");
    CHECK_CLOSE(offset, 0.01, 0.01, "Offset should be ~0.01");
    PASS();
}

static void test_apply_calibration(void) {
    TEST("Apply calibration correction");
    double corrected = ia_cal_apply_correction(10.01, 10.0, 0.01);
    CHECK_CLOSE(corrected, 1.0, 0.001, "Corrected should be 1.0");
    PASS();
}

static void test_linear_regression(void) {
    TEST("Multi-point linear regression");
    double x[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double y[] = {0.05, 10.05, 20.05, 30.05, 40.05};
    double gain, offset;
    ia_cal_linear_regression(x, y, 5, &gain, &offset);
    CHECK_CLOSE(gain, 10.0, 0.01, "Gain should be ~10");
    CHECK_CLOSE(offset, 0.05, 0.01, "Offset should be ~0.05");
    PASS();
}

/* ---- L6: Strain gauge ---- */
static void test_strain_gauge_output(void) {
    TEST("Strain gauge bridge output");
    strain_gauge_config_t cfg = {2.0, 1000.0, 5.0, 5.0, 350.0};
    double v = strain_gauge_bridge_output_uv(&cfg, 500.0);
    CHECK(v > 0.0, "Bridge output should be positive");
    PASS();
}

static void test_strain_gauge_to_ue(void) {
    TEST("Strain ADC to microstrain");
    strain_gauge_config_t cfg = {2.0, 1000.0, 5.0, 5.0, 350.0};
    double G = strain_gauge_required_gain(&cfg);
    double v_bridge = strain_gauge_bridge_output_uv(&cfg, 500.0) * 1e-6;
    double v_adc = v_bridge * G;
    double ue = strain_gauge_adc_to_ue(v_adc, G, &cfg);
    CHECK_CLOSE(ue, 500.0, 10.0, "Should recover ~500 uStrain");
    PASS();
}

/* ---- L6: Thermocouple ---- */
static void test_thermocouple_voltage(void) {
    TEST("Thermocouple voltage (K-type)");
    thermocouple_config_t cfg = {TC_TYPE_K, 25.0, 125.0};
    double v = thermocouple_voltage_uv(&cfg);
    CHECK(v > 3000.0 && v < 5000.0, "K-type 100C delta should be ~4mV");
    PASS();
}

/* ---- L6: RTD ---- */
static void test_rtd_resistance(void) {
    TEST("RTD Pt100 resistance");
    rtd_config_t cfg = {100.0, 1.0, 0.00385};
    double R = rtd_resistance_at_temp(&cfg, 100.0);
    CHECK_CLOSE(R, 138.5, 0.5, "Pt100 at 100C should be ~138.5 Ohm");
    PASS();
}

/* ---- L6: Load cell ---- */
static void test_load_cell(void) {
    TEST("Load cell voltage");
    load_cell_config_t cfg = {2.0, 10.0, 5.0, 0.0};
    double vfs = load_cell_fullscale_voltage_mv(&cfg);
    CHECK_CLOSE(vfs, 10.0, 0.1, "Full-scale should be 10 mV");
    PASS();
}

/* ---- L6: ECG ---- */
static void test_ecg_cm_noise(void) {
    TEST("ECG common-mode noise");
    double cm_noise = ecg_cm_noise_rti_uv(1.0, 90.0);
    CHECK(cm_noise < 100.0, "CM noise should be < 100 uV RTI at 90dB CMRR");
    PASS();
}

/* ---- L6: Current shunt ---- */
static void test_current_shunt(void) {
    TEST("Current shunt measurement");
    current_shunt_config_t cfg = {0.001, 10.0, 48.0, true};
    double v = current_shunt_vout_at_current(&cfg, 5.0);
    CHECK_CLOSE(v, 0.005, 1e-6, "5A through 1mOhm should be 5mV");
    PASS();
}

/* ---- L3: Bridge ---- */
static void test_bridge_output(void) {
    TEST("Wheatstone bridge output");
    bridge_sensor_t sensor = {350.0, 0.35, 5.0, 1, 10.0, 0.0};
    double v_exact = bridge_output_voltage(&sensor, true);
    double v_linear = bridge_output_voltage(&sensor, false);
    CHECK(v_exact > 0.0 && v_linear > 0.0, "Bridge outputs should be positive");
    CHECK(v_linear >= v_exact, "Linear approximation overestimates");
    PASS();
}

/* ---- L1: E96 resistor ---- */
static void test_e96(void) {
    TEST("E96 standard resistor value");
    double r = ia_e96_nearest(5555.56);
    CHECK(r >= 5400.0 && r <= 5700.0, "~5555 should round to nearest E96 value (5.49k or 5.62k)");
    PASS();
}

/* ---- L5: Polynomial calibration ---- */
static void test_polynomial_cal(void) {
    TEST("Quadratic polynomial calibration");
    double x[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double y[] = {1.0, 4.0, 9.0, 16.0, 25.0};
    double coeffs[3];
    ia_cal_polynomial_fit(x, y, 5, 2, coeffs);
    double y_pred = ia_cal_eval_polynomial(coeffs, 2, 6.0);
    CHECK_CLOSE(y_pred, 36.0, 1.0, "y=x^2 should predict 36 at x=6");
    PASS();
}

/* ---- L3: SNR ---- */
static void test_snr(void) {
    TEST("Signal-to-noise ratio");
    double snr = signal_to_noise_ratio_db(1.0, 0.001);
    CHECK_CLOSE(snr, 60.0, 0.1, "SNR should be 60dB for 1000:1 ratio");
    PASS();
}

/* ---- L4: Gain sensitivity ---- */
static void test_gain_sensitivity(void) {
    TEST("Gain sensitivity to R_G");
    double s = gain_sensitivity_to_rg(100.0);
    CHECK(fabs(s + 0.99) < 0.01, "Sensitivity should be ~-0.99 at G=100");
    PASS();
}

/* ---- L4: Gain temperature coefficient ---- */
static void test_gain_tempco(void) {
    TEST("Gain temperature coefficient");
    double tc = gain_tempco_ppm_per_c(100.0, 50.0, 100.0);
    CHECK(tc > 0.0, "Gain tempco should be positive");
    PASS();
}

/* ---- L1: E24 ---- */
static void test_e24(void) {
    TEST("E24 standard resistor value");
    double r = ia_e24_nearest(4300.0);
    CHECK_CLOSE(r, 4300.0, 200.0, "4300 should round to standard E24 value");
    PASS();
}

/* ---- L1: State init ---- */
static void test_state_init(void) {
    TEST("IA state initialization");
    ia_state_t state;
    ia_state_init(&state);
    CHECK(!state.powered, "State should be off");
    CHECK(!state.calibrated, "State should be uncalibrated");
    PASS();
}

/* ---- L1: Default spec ---- */
static void test_spec_defaults(void) {
    TEST("Default IA specification");
    ia_spec_t spec;
    ia_spec_init(&spec);
    CHECK(spec.gain == 1.0, "Default gain should be 1");
    CHECK(spec.cmrr_db > 80.0, "Default CMRR should be > 80dB");
    PASS();
}

int main(void) {
    printf("=== mini-instrumentation-amp Test Suite ===\n\n");

    test_three_opamp_gain();
    test_three_opamp_rg_selection();
    test_three_opamp_cmrr_resistor();
    test_three_opamp_output();
    test_two_opamp_gain();
    test_two_opamp_cmrr_zero();
    test_current_mode_gain();
    test_flying_cap_gain();
    test_indirect_current_gain();
    test_thermal_noise();
    test_shot_noise();
    test_noise_figure();
    test_error_budget();
    test_gain_error_budget();
    test_cmrr_combined();
    test_spec_validate();
    test_topology_recommend();
    test_two_point_cal();
    test_apply_calibration();
    test_linear_regression();
    test_strain_gauge_output();
    test_strain_gauge_to_ue();
    test_thermocouple_voltage();
    test_rtd_resistance();
    test_load_cell();
    test_ecg_cm_noise();
    test_current_shunt();
    test_bridge_output();
    test_e96();
    test_polynomial_cal();
    test_snr();
    test_gain_sensitivity();
    test_gain_tempco();
    test_e24();
    test_state_init();
    test_spec_defaults();

    printf("\n=== Results: %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}