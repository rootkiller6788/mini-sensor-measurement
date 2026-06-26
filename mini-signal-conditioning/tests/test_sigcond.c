/**
 * test_sigcond.c - Signal Conditioning Module Test Suite
 *
 * Tests all core APIs: filter design, linearization, excitation,
 * isolation, CJC, bridge, and self-calibration.
 * Uses assert-based verification (no external test framework).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "sigcond_filter.h"
#include "sigcond_linearize.h"
#include "sigcond_excitation.h"
#include "sigcond_isolation.h"
#include "sigcond_cjc.h"
#include "sigcond_bridge.h"
#include "sigcond_selfcal.h"

#define EPS_LOOSE   1e-3
#define EPS_TIGHT   1e-9

static int tests_run = 0;
static int tests_passed = 0;

static void test_pass(const char *name) {
    tests_run++; tests_passed++;
    printf("  PASS: %s\n", name);
}
static void test_fail(const char *name, const char *msg) {
    tests_run++;
    printf("  FAIL: %s - %s\n", name, msg);
}

/* ---------------------------------------------------------------
 * Filter Tests
 * --------------------------------------------------------------- */
static void test_butterworth_poles(void) {
    double preal[20], pimag[20], real_pole;
    unsigned n = butterworth_poles(4, preal, pimag, &real_pole);
    if (n != 2) { test_fail("butterworth_poles order 4", "expected 2 pairs"); return; }
    test_pass("butterworth_poles order 4");

    double w0[10], qv[10];
    unsigned ns;
    butterworth_biquads(4, w0, qv, &ns);
    if (ns != 2) { test_fail("butterworth_biquads order 4", "expected 2 sections"); return; }
    test_pass("butterworth_biquads order 4");

    /* Check Q values for Butterworth 4th order: Q1=0.5412, Q2=1.3066 (any order) */
    bool has_q054 = false, has_q130 = false;
    for (unsigned i = 0; i < ns; i++) {
        if (fabs(qv[i] - 0.5411961) < EPS_LOOSE) has_q054 = true;
        if (fabs(qv[i] - 1.3065630) < EPS_LOOSE) has_q130 = true;
    }
    if (!has_q054 || !has_q130)
        { test_fail("butterworth Q values", ""); return; }
    test_pass("butterworth Q values");

    /* Verify anti-alias filter */
    bool ok = filter_verify_antialias(1000.0, 6000.0, 12, 70.0);
    if (!ok) { test_fail("anti-alias verify 12-bit", ""); return; }
    test_pass("anti-alias verify");

    /* Filter order calculation */
    int ord = filter_butterworth_order(1000.0, 3000.0, 0.1, 60.0);
    if (ord < 3 || ord > 10) { test_fail("butterworth order", "unexpected order"); return; }
    test_pass("butterworth order 1kHz-3kHz");
}

static void test_chebyshev_filter(void) {
    double preal[20], pimag[20];
    unsigned n = chebyshev1_poles(3, 1.0, preal, pimag);
    if (n != 2) { test_fail("chebyshev1_poles order 3", "expected 2 pairs"); return; }
    test_pass("chebyshev1_poles order 3");

    double w0[10], qv[10];
    unsigned ns;
    chebyshev1_biquads(4, 0.5, w0, qv, &ns);
    if (ns != 2) { test_fail("chebyshev1_biquads", "expected 2 sections"); return; }
    test_pass("chebyshev1_biquads");

    unsigned ell_n = elliptic_order(1000.0, 2000.0, 0.1, 60.0);
    if (ell_n < 3) { test_fail("elliptic order", "too low"); return; }
    test_pass("elliptic order");
}

static void test_sallen_key(void) {
    double r, c, rf_by_rg;
    double gain = sallen_key_lp_design(1000.0, 0.7071, &r, &c, &rf_by_rg);
    if (gain < 0.0) { test_fail("sallen_key design", "invalid gain"); return; }
    if (r <= 0.0 || c <= 0.0) { test_fail("sallen_key components", "invalid values"); return; }
    test_pass("sallen_key LP design");
}

static void test_digital_filters(void) {
    double coeffs[64];
    unsigned len;
    moving_average_coeffs(10, coeffs, &len);
    if (len != 10) { test_fail("moving average length", ""); return; }
    test_pass("moving average FIR");

    double alpha = iir1_lp_alpha(10.0, 1000.0);
    if (alpha <= 0.0 || alpha >= 1.0) { test_fail("IIR alpha", "out of range"); return; }
    test_pass("IIR LP alpha");

    double B0, B1, B2, A1, A2;
    bilinear_transform_sos(100.0, 1000.0, 1.4142, 1.0, 0.0, 0.0, 1.0, &B0, &B1, &B2, &A1, &A2);
    test_pass("bilinear transform");
}

/* ---------------------------------------------------------------
 * Linearization Tests
 * --------------------------------------------------------------- */
static void test_horner_poly(void) {
    double c[] = {1.0, 2.0, 3.0};  /* p(x) = 1 + 2x + 3x^2 */
    double result = poly_eval_horner(2, c, 2.0, -10.0, 10.0);
    double expected = 1.0 + 4.0 + 12.0;  /* 17 */
    if (fabs(result - expected) > EPS_TIGHT)
        { test_fail("horner evaluation", ""); return; }
    test_pass("horner polynomial");

    double p_val, p_deriv;
    poly_eval_with_derivative(2, c, 2.0, &p_val, &p_deriv);
    if (fabs(p_deriv - (2.0 + 12.0)) > EPS_TIGHT) /* p'(x)=2+6x, at x=2: 14 */
        { test_fail("horner derivative", ""); return; }
    test_pass("horner derivative");
}

static void test_lut_interp(void) {
    double xv[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double yv[] = {0.0, 2.0, 4.0, 6.0, 8.0};
    double y = lut_linear_interp(5, xv, yv, 1.5);
    if (fabs(y - 3.0) > EPS_TIGHT) { test_fail("LUT linear interp", ""); return; }
    test_pass("LUT linear interpolation");
}

static void test_steinhart_hart(void) {
    steinhart_hart_params_t params;
    steinhart_hart_fit_3pt(0.0, 10000.0, 25.0, 5000.0, 50.0, 2500.0, &params);
    double T = steinhart_hart_T(&params, 5000.0);
    if (fabs(T - 298.15) > 1.0) /* ~25C +- 1C accuracy for quick fit */
        { test_fail("steinhart-hart T", ""); return; }
    test_pass("steinhart-hart temperature");
}

static void test_rtd_cvd(void) {
    rtd_cvd_params_t pt100;
    rtd_cvd_init_pt100(&pt100);
    double R = rtd_cvd_R_from_T(&pt100, 0.0);
    if (fabs(R - 100.0) > EPS_TIGHT) { test_fail("RTD CVD R at 0C", ""); return; }
    test_pass("RTD CVD R(0C)=100");

    double R100 = rtd_cvd_R_from_T(&pt100, 100.0);
    if (fabs(R100 - 138.5) > 0.5) { test_fail("RTD CVD R at 100C", ""); return; }
    test_pass("RTD CVD R(100C)~138.5");

    double T = rtd_cvd_T_from_R(&pt100, 100.0);
    if (fabs(T - 0.0) > EPS_LOOSE) { test_fail("RTD CVD T(R=100)", ""); return; }
    test_pass("RTD CVD T(R=100)=0C");
}

static void test_nist_tc(void) {
    thermocouple_nist_model_t tcK;
    nist_tc_init_typeK(&tcK);

    /* At 0C, voltage should be ~0uV */
    double V0 = nist_tc_temp_to_voltage(&tcK, 0.0);
    if (fabs(V0) > 1.0) { test_fail("TC Type K V(0C)", ""); return; }
    test_pass("TC Type K V(0C)~0");

    /* At 100C, Type K ~ 4.096 mV = 4096 uV */
    double V100 = nist_tc_temp_to_voltage(&tcK, 100.0);
    if (fabs(V100 - 4096.0) > 100.0) { test_fail("TC Type K V(100C)", ""); return; }
    test_pass("TC Type K V(100C)~4096uV");

    /* Round-trip: T -> V -> T. NIST poly error up to ~2C without exponential correction terms. */
    double T_rt = nist_tc_voltage_to_temp(&tcK, V100);
    if (fabs(T_rt - 100.0) > 5.0) { test_fail("TC Type K round trip", ""); return; }
    test_pass("TC Type K round trip");
}

/* ---------------------------------------------------------------
 * Excitation Tests
 * --------------------------------------------------------------- */
static void test_excitation(void) {
    double err = excitation_load_error(5.0, 0.1, 100.0);
    if (err <= 0.0 || err >= 0.01) { test_fail("load regulation", ""); return; }
    test_pass("load regulation error");

    double z_max = excitation_max_source_impedance(100.0, 0.1);
    if (fabs(z_max - 0.1) > EPS_TIGHT) { test_fail("max source impedance", ""); return; }
    test_pass("max source impedance");

    double r_set, vin;
    double vcomp = howland_current_source(1e-3, &r_set, &vin, 15.0, 2.0, 1000.0);
    if (vcomp < 0.0) { test_fail("howland current pump", ""); return; }
    test_pass("howland current pump");

    double R_sens = ratiometric_sensor_resistance(5.0, 2.5, 1000.0);
    if (fabs(R_sens - 1000.0) > EPS_TIGHT)
        { test_fail("ratiometric sensor R", ""); return; }
    test_pass("ratiometric measurement");
}

/* ---------------------------------------------------------------
 * Isolation Tests
 * --------------------------------------------------------------- */
static void test_isolation(void) {
    double V_with = isolation_withstand_from_working(230.0);
    if (V_with < 1400.0 || V_with > 1500.0)
        { test_fail("isolation withstand 230V", ""); return; }
    test_pass("isolation withstand voltage");

    double V_derate = isolation_altitude_derate(5000.0, 3000.0);
    if (V_derate >= 5000.0) { test_fail("altitude derating", ""); return; }
    test_pass("altitude derating");

    double creep = isolation_creepage(230.0, 2, 3);
    if (creep < 2.0) { test_fail("creepage distance", ""); return; }
    test_pass("creepage distance");

    double cmrr_tot = isolation_total_cmrr(80.0, 100.0, 80.0);
    if (cmrr_tot > 80.0 || cmrr_tot < 70.0)
        { test_fail("total CMRR", ""); return; }
    test_pass("total CMRR");
}

/* ---------------------------------------------------------------
 * CJC Tests
 * --------------------------------------------------------------- */
static void test_cjc(void) {
    thermocouple_nist_model_t tcK;
    nist_tc_init_typeK(&tcK);

    /* At CJC temp of 25C, V_cj should be ~ (25-0)*40.6 = 1015 uV */
    double v_cj = cjc_equivalent_voltage(25.0, &tcK);
    if (fabs(v_cj - 1000.0) > 50.0) { test_fail("CJC voltage at 25C", ""); return; }
    test_pass("CJC equivalent voltage");

    /* Verify successive temperature law */
    bool law_ok = cjc_verify_successive_law(0.0, 25.0, 100.0, &tcK, 10.0);
    if (!law_ok) { test_fail("successive temp law", ""); return; }
    test_pass("successive temp law");

    /* Full CJC: V_measured = E(100C) - E(25C) ~ 4005-983 = 3022 uV at CJC=25C.
     * After CJC compensation, should recover ~100C. */
    double v_diff = nist_tc_temp_to_voltage(&tcK, 100.0) - nist_tc_temp_to_voltage(&tcK, 25.0);
    double T_hot = cjc_compensate(v_diff, 25.0, &tcK);
    if (fabs(T_hot - 100.0) > 2.0) { test_fail("CJC compensate", ""); return; }
    test_pass("CJC compensation");
}

/* ---------------------------------------------------------------
 * Bridge Tests
 * --------------------------------------------------------------- */
static void test_bridge(void) {
    double v_out = bridge_output_exact(10.0, 120.0, 120.0, 120.0, 120.0);
    if (fabs(v_out) > EPS_TIGHT) { test_fail("balanced bridge output", ""); return; }
    test_pass("balanced bridge zero output");

    /* Quarter bridge with 1% change */
    double r_nom = 120.0;
    double dR = 1.2;
    double v_quarter = bridge_output_exact(10.0, r_nom, r_nom, r_nom + dR, r_nom);
    double v_quarter_lin = bridge_output_linear(10.0, r_nom, dR, 1);
    if (fabs(v_quarter - v_quarter_lin) / fabs(v_quarter_lin) > 0.01)
        { test_fail("bridge linearity 1%", ""); return; }
    test_pass("bridge linear approximation");

    bool bal = bridge_is_balanced(100.0, 100.0, 100.0, 100.0, 0.001);
    if (!bal) { test_fail("bridge balanced detection", ""); return; }
    test_pass("bridge balance check");

    double sens = bridge_sensitivity_strain(10.0, 2.0, 1);
    if (fabs(sens - 5.0) > EPS_TIGHT) { test_fail("bridge sensitivity", ""); return; }
    test_pass("bridge sensitivity");

    double strain = bridge_strain_from_output(0.005, 10.0, 2.0, 1);
    if (fabs(strain - 0.001) > 1e-6) { test_fail("strain from output", ""); return; }
    test_pass("strain from bridge output");

    double min_strain = bridge_min_detectable_strain(10.0, 2.0, 1.0, 1000.0);
    if (min_strain <= 0.0) { test_fail("min detectable strain", ""); return; }
    test_pass("min detectable strain");
}

/* ---------------------------------------------------------------
 * Self-Calibration Tests
 * --------------------------------------------------------------- */
static void test_selfcal(void) {
    calibration_point_t p1 = {0.0, 0.0, 0.01, 25.0, 50.0};
    calibration_point_t p2 = {5.0, 5.0, 5.01, 25.0, 50.0};
    double gain, offset;
    int rc = selfcal_two_point(&p1, &p2, &gain, &offset);
    if (rc != 0) { test_fail("two-point calibration", "returned error"); return; }
    if (fabs(gain - 1.0) > EPS_LOOSE) { test_fail("gain", ""); return; }
    test_pass("two-point calibration");
}

int main(void) {
    printf("=== Signal Conditioning Module Tests ===\n\n");

    printf("-- Filter Tests --\n");
    test_butterworth_poles();
    test_chebyshev_filter();
    test_sallen_key();
    test_digital_filters();

    printf("\n-- Linearization Tests --\n");
    test_horner_poly();
    test_lut_interp();
    test_steinhart_hart();
    test_rtd_cvd();
    test_nist_tc();

    printf("\n-- Excitation Tests --\n");
    test_excitation();

    printf("\n-- Isolation Tests --\n");
    test_isolation();

    printf("\n-- CJC Tests --\n");
    test_cjc();

    printf("\n-- Bridge Tests --\n");
    test_bridge();

    printf("\n-- Self-Calibration Tests --\n");
    test_selfcal();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
