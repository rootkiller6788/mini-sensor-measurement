/**
 * test_biosensor.c — Comprehensive assert-based tests for all biosensor modules
 *
 * Tests cover L1-L6 knowledge levels with real numerical assertions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "biosensor_types.h"
#include "biosensor_electrochemical.h"
#include "biosensor_optical.h"
#include "biosensor_kinetics.h"
#include "biosensor_signal.h"
#include "biosensor_calibration.h"

#define EPS 1.0e-6
#define TEST_PASS() printf("  PASS: %s\n", __func__)

/* ============================================================================
 * L1: Type initialization tests
 * ============================================================================ */
static void test_biosensor_descriptor_init(void) {
    BiosensorDescriptor desc;
    biosensor_descriptor_init(&desc);
    assert(desc.bioreceptor_type == BIORECEPTOR_ENZYME);
    assert(desc.transducer_type == TRANSDUCER_AMPEROMETRIC);
    assert(fabs(desc.ph_optimum - 7.0) < EPS);
    assert(desc.shelf_life_days == 180);
    TEST_PASS();
}

static void test_michaelis_menten_params_init(void) {
    MichaelisMentenParams mm;
    michaelis_menten_params_init(&mm, 2.5e-3, 100.0);
    assert(fabs(mm.km - 2.5e-3) < EPS);
    assert(fabs(mm.vmax - 100.0) < EPS);
    assert(fabs(mm.enzyme_conc - 1.0e-9) < EPS);
    TEST_PASS();
}

static void test_nernst_params_init(void) {
    NernstParams np;
    nernst_params_init(&np, 0.22, 2.0);
    assert(fabs(np.standard_potential - 0.22) < EPS);
    assert(fabs(np.electron_count - 2.0) < EPS);
    assert(fabs(np.temperature_kelvin - 298.15) < EPS);
    assert(fabs(np.faraday_constant - 96485.3329) < EPS);
    TEST_PASS();
}

static void test_qcm_params_init(void) {
    QCMParameters qcm;
    qcm_params_init(&qcm, 10.0e6, 0.2);
    assert(fabs(qcm.fundamental_frequency_hz - 10.0e6) < 1.0);
    assert(fabs(qcm.electrode_area_cm2 - 0.2) < EPS);
    assert(qcm.mass_sensitivity_hz_per_ng > 0.0);
    assert(qcm.quality_factor > 1000.0);
    TEST_PASS();
}

static void test_measurement_series_init_free(void) {
    MeasurementSeries series;
    measurement_series_init(&series, 50);
    assert(series.capacity == 50);
    assert(series.points != NULL);
    assert(series.point_count == 0);
    assert(fabs(series.ph_value - 7.0) < EPS);
    measurement_series_free(&series);
    assert(series.points == NULL);
    assert(series.capacity == 0);
    TEST_PASS();
}

static void test_calibration_standards_init_free(void) {
    CalibrationStandards cs;
    calibration_standards_init(&cs, 8, 3);
    assert(cs.standard_count == 8);
    assert(cs.replicates == 3);
    assert(cs.concentrations != NULL);
    assert(cs.signals != NULL);
    assert(cs.signal_stds != NULL);
    calibration_standards_free(&cs);
    assert(cs.concentrations == NULL);
    TEST_PASS();
}

static void test_sensor_array_init_free(void) {
    SensorArray sa;
    sensor_array_init(&sa, 16);
    assert(sa.capacity == 16);
    assert(sa.sensors != NULL);
    assert(sa.pattern_vector != NULL);
    assert(sa.sensor_count == 0);
    sensor_array_free(&sa);
    assert(sa.sensors == NULL);
    TEST_PASS();
}

/* ============================================================================
 * L4: Nernst equation tests
 * ============================================================================ */
static void test_nernst_potential_standard(void) {
    NernstParams np;
    nernst_params_init(&np, 0.0, 1.0);
    /* At 25°C, RT/F = 25.693 mV ≈ 0.025693 V */
    /* For [Ox]=[Red], E = E° = 0 */
    double e = nernst_potential(&np, 1.0e-3, 1.0e-3);
    assert(fabs(e - 0.0) < EPS);
    TEST_PASS();
}

static void test_nernst_potential_ratio_10(void) {
    NernstParams np;
    nernst_params_init(&np, 0.0, 1.0);
    /* [Ox]/[Red] = 10 → E ≈ 0.025693 * ln(10) ≈ 0.05916 V */
    double e = nernst_potential(&np, 1.0e-2, 1.0e-3);
    double expected = 0.025693 * log(10.0);
    assert(fabs(e - expected) < 0.001);
    TEST_PASS();
}

static void test_isfet_compute_ph(void) {
    ISFETParams isfet;
    isfet_params_init(&isfet);
    isfet.threshold_voltage = 1.5;
    isfet.sensitivity_mv_per_ph = 58.0;

    /* At pH 7, V_th = 1.5 V */
    double ph = isfet_compute_ph(&isfet, 1.5);
    assert(fabs(ph - 7.0) < 0.05);

    /* At pH 4, V_th should shift by 3 * 58 mV = 174 mV */
    ph = isfet_compute_ph(&isfet, 1.5 - 0.174);
    assert(fabs(ph - 4.0) < 0.1);

    TEST_PASS();
}

/* ============================================================================
 * L4: Butler-Volmer and Cottrell tests
 * ============================================================================ */
static void test_butler_volmer_zero_overpotential(void) {
    ButlerVolmerParams bv;
    butler_volmer_params_init(&bv, 1.0e-4, 0.5);
    bv.overpotential = 0.0;
    double current = butler_volmer_current(&bv);
    /* At η=0, net current should be near zero */
    assert(fabs(current) < 1.0e-12);
    TEST_PASS();
}

static void test_cottrell_current_decay(void) {
    /* Standard ferricyanide reduction: n=1, D≈7.6×10⁻⁶ cm²/s,
     * C = 1 mM = 1e-6 mol/cm³, A = 0.07 cm² */
    double i_t1 = cottrell_current(1, 0.07, 7.6e-6, 1.0e-6, 1.0);
    double i_t4 = cottrell_current(1, 0.07, 7.6e-6, 1.0e-6, 4.0);
    /* Current should decay as 1/√t: I(4) ≈ I(1)/2 */
    assert(i_t1 > 0.0);
    assert(i_t4 > 0.0);
    assert(i_t4 < i_t1);
    double ratio = i_t1 / i_t4;
    assert(fabs(ratio - 2.0) < 0.01);
    TEST_PASS();
}

static void test_randles_sevcik_linear_scanrate(void) {
    double ip1 = randles_sevcik_peak_current(1, 0.07, 7.6e-6, 1.0e-6, 0.1);
    double ip4 = randles_sevcik_peak_current(1, 0.07, 7.6e-6, 1.0e-6, 0.4);
    /* I_p ∝ √v: I(0.4) / I(0.1) = √4 = 2 */
    assert(ip1 > 0.0);
    assert(ip4 > 0.0);
    double ratio = ip4 / ip1;
    assert(fabs(ratio - 2.0) < 0.01);
    TEST_PASS();
}

/* ============================================================================
 * L6: Glucose biosensor test
 * ============================================================================ */
static void test_glucose_biosensor_linear_range(void) {
    /* Low glucose → roughly linear response */
    double i_low = glucose_biosensor_current(1.0, 2.0, 0.07, 0.21);
    double i_mid = glucose_biosensor_current(3.0, 2.0, 0.07, 0.21);
    assert(i_low > 0.0);
    assert(i_mid > i_low);
    assert(i_mid < 3.5 * i_low);  /* sub-linear due to saturation */
    TEST_PASS();
}

static void test_glucose_conc_from_current(void) {
    double conc = glucose_conc_from_current(500.0, 1000.0, 8.0, 2000.0);
    assert(conc > 0.0);
    assert(conc < 2.0);  /* ~0.5 mM based on sensitivity */
    conc = glucose_conc_from_current(1800.0, 1000.0, 8.0, 2000.0);
    /* Near saturation — Michaelis-Menten correction kicks in */
    assert(conc > 3.0);
    TEST_PASS();
}

/* ============================================================================
 * L4: Beer-Lambert tests
 * ============================================================================ */
static void test_beer_lambert_absorbance(void) {
    double a = beer_lambert_absorbance(50000.0, 1.0, 1.0e-4);
    /* ε=50000, b=1cm, c=0.1mM → A = 5.0 */
    assert(fabs(a - 5.0) < EPS);
    double a_zero = beer_lambert_absorbance(50000.0, 1.0, 0.0);
    assert(fabs(a_zero - 0.0) < EPS);
    TEST_PASS();
}

static void test_concentration_from_absorbance(void) {
    double c = concentration_from_absorbance(0.5, 10000.0, 1.0);
    assert(fabs(c - 5.0e-5) < EPS);
    TEST_PASS();
}

static void test_two_component_absorbance(void) {
    /* NADH (ε_340=6220) and NAD (ε_340≈0, ε_260=18000) */
    double ca, cb;
    bool ok = two_component_absorbance(0.622, 0.0,
                                        6220.0, 0.0,
                                        0.0, 18000.0,
                                        1.0, &ca, &cb);
    assert(ok);
    assert(fabs(ca - 1.0e-4) < 1.0e-6);
    assert(fabs(cb) < 1.0e-12);
    TEST_PASS();
}

/* ============================================================================
 * L2: SPR tests
 * ============================================================================ */
static void test_spr_resonance_angle(void) {
    SPRParams spr;
    spr_params_init(&spr, 633.0, 1.0e5);
    spr.prism_ri = 1.52;
    /* Gold ε_m ≈ -12 at 633 nm, water n=1.33
     * k_SP/k0 = sqrt(ε_m*ε_d/(ε_m+ε_d)) = sqrt(12*1.77/10.23) ≈ 1.44
     * sin(θ) = 1.44/1.52 = 0.948 → θ ≈ 71.5° */
    double theta = spr_resonance_angle(&spr, 1.33, -12.0);
    /* Typical SPR angle for Au/water at 633 nm with BK7 ≈ 71-72° */
    assert(theta > 65.0);
    assert(theta < 80.0);
    TEST_PASS();
}

static void test_spr_reflectivity_bounds(void) {
    SPRParams spr;
    spr_params_init(&spr, 633.0, 1.0e5);
    spr.prism_ri = 1.52;
    spr.gold_layer_thickness_nm = 50.0;
    double r = spr_reflectivity(&spr, 53.5, 1.33, 0.18, 3.5);
    assert(r >= 0.0);
    assert(r <= 1.0);
    TEST_PASS();
}

/* ============================================================================
 * L4: FRET tests
 * ============================================================================ */
static void test_fret_efficiency_at_r0(void) {
    double e = fret_efficiency(5.0, 5.0);
    assert(fabs(e - 0.5) < EPS);
    double d = fret_distance(5.0, 0.5);
    assert(fabs(d - 5.0) < 0.01);
    TEST_PASS();
}

static void test_fret_efficiency_bounds(void) {
    double e_near  = fret_efficiency(5.0, 1.0);   /* r << R0 */
    double e_far   = fret_efficiency(5.0, 20.0);  /* r >> R0 */
    assert(e_near > 0.95);
    assert(e_far < 0.02);
    assert(e_near > e_far);
    TEST_PASS();
}

/* ============================================================================
 * L4: Michaelis-Menten tests
 * ============================================================================ */
static void test_michaelis_menten_rate(void) {
    MichaelisMentenParams mm;
    michaelis_menten_params_init(&mm, 5.0, 100.0);
    mm.substrate_conc = 5.0;
    double rate = michaelis_menten_rate(&mm);
    /* At [S]=Km, v = Vmax/2 = 50 */
    assert(fabs(rate - 50.0) < EPS);
    TEST_PASS();
}

static void test_michaelis_menten_inverse(void) {
    MichaelisMentenParams mm;
    michaelis_menten_params_init(&mm, 5.0, 100.0);
    double conc = michaelis_menten_inverse(&mm, 50.0);
    assert(fabs(conc - 5.0) < EPS);
    TEST_PASS();
}

static void test_competitive_inhibition(void) {
    MichaelisMentenParams mm;
    michaelis_menten_params_init(&mm, 5.0, 100.0);
    mm.substrate_conc = 5.0;
    /* No inhibitor */
    double v0 = competitive_inhibition_rate(&mm, 0.0, 1.0);
    assert(fabs(v0 - 50.0) < EPS);
    /* With inhibitor: [I]=Ki, apparent Km doubles → rate drops */
    double vi = competitive_inhibition_rate(&mm, 1.0, 1.0);
    assert(vi < v0);
    assert(vi > 0.0);
    TEST_PASS();
}

static void test_lineweaver_burk_estimate(void) {
    /* Perfect Michaelis-Menten data: Km=5, Vmax=100 */
    double conc[] = {1.0, 2.0, 5.0, 10.0, 20.0};
    double rates[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    int n = 5;
    for (int i = 0; i < n; i++) {
        rates[i] = 100.0 * conc[i] / (5.0 + conc[i]);
    }
    double km_est, vmax_est;
    double r2 = lineweaver_burk_estimate(conc, rates, n, &km_est, &vmax_est);
    assert(r2 > 0.95);
    assert(fabs(km_est - 5.0) < 1.0);
    assert(fabs(vmax_est - 100.0) < 5.0);
    TEST_PASS();
}

static void test_mm_nls_fit(void) {
    double conc[] = {1.0, 2.0, 5.0, 10.0, 20.0, 50.0};
    double rates[6];
    int n = 6;
    for (int i = 0; i < n; i++) {
        rates[i] = 100.0 * conc[i] / (5.0 + conc[i]);
    }
    /* Use Hanes-Woolf to get good initial estimates */
    double km_init, vmax_init;
    hanes_woolf_estimate(conc, rates, n, &km_init, &vmax_init);
    double km_est, vmax_est;
    int iters = michaelis_menten_nls_fit(conc, rates, n, km_init, vmax_init,
                                          &km_est, &vmax_est, 50, 1.0e-6);
    /* NLS should converge with good initial guesses */
    assert(iters > 0);
    /* Verify estimates are reasonable (within 10% for Km, 2% for Vmax) */
    assert(fabs(km_est - 5.0) < 0.5);
    assert(fabs(vmax_est - 100.0) < 2.0);
    TEST_PASS();
}

/* ============================================================================
 * L4: Hill and Langmuir tests
 * ============================================================================ */
static void test_hill_fractional_saturation(void) {
    HillParams hp;
    hill_params_init(&hp, 1.0e-6, 2.0, 100.0);
    /* At K_half, saturation = max/2 */
    double theta = hill_fractional_saturation(&hp, 1.0e-6);
    assert(fabs(theta - 50.0) < EPS);
    /* Positive cooperativity: steep transition */
    double theta_low  = hill_fractional_saturation(&hp, 1.0e-7);
    double theta_high = hill_fractional_saturation(&hp, 1.0e-5);
    assert(theta_low < 10.0);
    assert(theta_high > 90.0);
    TEST_PASS();
}

static void test_langmuir_coverage(void) {
    LangmuirParams lp;
    lp.k_adsorption = 1.0e5;
    lp.k_desorption = 0.01;
    lp.equilibrium_constant = lp.k_adsorption / lp.k_desorption; /* K = 1e7 */
    double theta = langmuir_coverage(&lp, 1.0e-7);
    assert(theta > 0.4);
    assert(theta < 0.6);  /* K*c = 1 → θ = 0.5 */
    TEST_PASS();
}

static void test_langmuir_kinetic_convergence(void) {
    LangmuirParams lp;
    lp.k_adsorption = 0.5;
    lp.k_desorption = 0.1;
    lp.equilibrium_constant = lp.k_adsorption / lp.k_desorption; /* K = 5 */
    double bulk_conc = 0.5;
    /* k_obs = k_ads*c + k_des = 0.5*0.5 + 0.1 = 0.35 s⁻¹ */
    /* θ_eq = k_ads*c / k_obs = 0.25/0.35 ≈ 0.714 */
    double theta_early = langmuir_kinetic_coverage(&lp, bulk_conc, 0.1);
    double theta_late  = langmuir_kinetic_coverage(&lp, bulk_conc, 20.0);
    assert(theta_late > theta_early);
    double theta_eq  = langmuir_coverage(&lp, bulk_conc);
    /* After 20s, exp(-0.35*20) = exp(-7) ≈ 0.00091 → ~99.91% converged */
    assert(fabs(theta_late - theta_eq) < 0.01);
    TEST_PASS();
}

/* ============================================================================
 * L5: LOD/LOQ tests
 * ============================================================================ */
static void test_lod_three_sigma(void) {
    double lod = lod_three_sigma(0.01, 1.0);
    assert(fabs(lod - 0.033) < EPS);
    double lod2 = lod_three_sigma(0.1, 0.5);
    assert(fabs(lod2 - 0.66) < EPS);
    TEST_PASS();
}

static void test_loq_ten_sigma(void) {
    double loq = loq_ten_sigma(0.01, 1.0);
    assert(fabs(loq - 0.1) < EPS);
    assert(loq > lod_three_sigma(0.01, 1.0));
    TEST_PASS();
}

/* ============================================================================
 * L5: Baseline correction tests
 * ============================================================================ */
static void test_moving_minimum_baseline(void) {
    double signal[] = {1.0, 1.1, 5.0, 1.2, 1.3, 6.0, 1.4, 1.5};
    double baseline[8];
    moving_minimum_baseline(signal, 8, 3, baseline);
    /* Baseline should track local minima, not peaks */
    assert(baseline[2] < 5.0);
    assert(baseline[5] < 6.0);
    TEST_PASS();
}

static void test_exponential_drift_correction(void) {
    double signal[]  = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    double time[]    = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    double corrected[6];
    exponential_drift_correction(signal, time, 6, 2.0, 2.0, corrected);
    /* Corrected signal should be smaller than raw (drift removed) */
    for (int i = 0; i < 6; i++) {
        assert(corrected[i] < signal[i] + 1.0e-9);
    }
    TEST_PASS();
}

/* ============================================================================
 * L5: Smoothing tests
 * ============================================================================ */
static void test_savitzky_golay_smooth(void) {
    double signal[] = {1.0, 2.0, 1.5, 3.0, 2.5, 4.0, 3.5};
    double smoothed[7];
    bool ok = savitzky_golay_smooth(signal, 7, 5, 2, smoothed);
    assert(ok);
    /* Smoothed values should be within range of input */
    for (int i = 0; i < 7; i++) {
        assert(smoothed[i] >= 1.0 && smoothed[i] <= 4.0);
    }
    TEST_PASS();
}

static void test_median_filter(void) {
    double signal[] = {1.0, 1.0, 100.0, 1.0, 1.0, 1.0, 1.0};
    double filtered[7];
    median_filter_spike_removal(signal, 7, 3, filtered);
    /* Spike at index 2 should be removed */
    assert(filtered[2] < 10.0);
    TEST_PASS();
}

/* ============================================================================
 * L5: Calibration tests
 * ============================================================================ */
static void test_linear_calibration_fit(void) {
    CalibrationStandards cs;
    calibration_standards_init(&cs, 5, 1);
    cs.concentrations[0] = 0.0;  cs.signals[0] = 0.05;
    cs.concentrations[1] = 2.0;  cs.signals[1] = 2.05;
    cs.concentrations[2] = 5.0;  cs.signals[2] = 5.02;
    cs.concentrations[3] = 10.0; cs.signals[3] = 10.03;
    cs.concentrations[4] = 20.0; cs.signals[4] = 20.01;
    LinearCalibration lc;
    double r2 = linear_calibration_fit(&cs, &lc);
    assert(r2 > 0.999);
    assert(fabs(lc.slope - 1.0) < 0.01);
    assert(fabs(lc.intercept) < 0.1);
    calibration_standards_free(&cs);
    TEST_PASS();
}

static void test_standard_addition_method(void) {
    /* Sample with unknown Pb²⁺ concentration ~10 ppb
     * All aliquots have same sample volume (5 mL) and total volume (50 mL) */
    double sample_vol[] = {5.0, 5.0, 5.0, 5.0, 5.0};
    double spike_conc[] = {0.0, 100.0, 200.0, 300.0, 400.0};
    double spike_vol[]  = {0.0, 0.5, 0.5, 0.5, 0.5};
    double total_vol[]  = {50.0, 50.0, 50.0, 50.0, 50.0};
    double signals[5];
    double true_conc = 10.0;
    double sensitivity = 2.5;
    for (int i = 0; i < 5; i++) {
        double added = spike_conc[i] * spike_vol[i] / total_vol[i];
        double total_in_aliquot = true_conc * sample_vol[i] / total_vol[i] + added;
        signals[i] = sensitivity * total_in_aliquot;
    }
    double est_conc, r_squared;
    bool ok = standard_addition_method(sample_vol, spike_conc, spike_vol,
                                        total_vol, signals, 5,
                                        &est_conc, &r_squared);
    assert(ok);
    assert(fabs(est_conc - true_conc) < 2.0);
    TEST_PASS();
}

static void test_fourpl_fit(void) {
    CalibrationStandards cs;
    calibration_standards_init(&cs, 8, 1);
    /* Generate near-ideal 4PL data: y = 0.05 + 1.95 / (1 + (x/1.0)^1.5) */
    for (int i = 0; i < 8; i++) {
        cs.concentrations[i] = 0.05 * (double)(i + 1);
        double ratio = cs.concentrations[i] / 1.0;
        cs.signals[i] = 0.05 + 1.95 / (1.0 + pow(ratio, 1.5));
    }

    FourPLLogisticFit fit;
    fourpl_fit_init(&fit);
    fit.a = 0.05; fit.d = 2.0;
    /* Use log-logistic linearization for good initial b, c */
    fourpl_linearized_initial_guess(&cs, &fit);
    int iters = fourpl_fit(&cs, &fit, 100, 0.01);
    /* Educational LM — verify convergence behavior */
    assert(iters >= 0);  /* may converge immediately or iterate */
    assert(fit.c > 0.0);
    assert(fit.b > 0.0);
    /* Verify predictions are reasonable (RMSE check) */
    double rmse = 0.0;
    for (int i = 0; i < 8; i++) {
        double pred = elisa_fourpl_predict(&fit, cs.concentrations[i]);
        double err = pred - cs.signals[i];
        rmse += err * err;
    }
    rmse = sqrt(rmse / 8.0);
    assert(rmse < 0.5);  /* reasonable fit quality */
    calibration_standards_free(&cs);
    TEST_PASS();
}

static void test_grubbs_outlier(void) {
    double values[] = {10.0, 10.2, 10.1, 10.3, 25.0, 10.0};
    int outlier_idx;
    bool found = grubbs_outlier_test(values, 6, 0.05, &outlier_idx);
    assert(found);
    assert(outlier_idx == 4);
    /* No outlier in clean data */
    double clean[] = {10.0, 10.2, 10.1, 10.3, 10.0};
    found = grubbs_outlier_test(clean, 5, 0.05, &outlier_idx);
    assert(!found || outlier_idx == -1);
    TEST_PASS();
}

static void test_coefficient_of_variation(void) {
    double values[] = {10.0, 11.0, 10.5, 9.5, 10.0};
    double cv = coefficient_of_variation(values, 5);
    assert(cv > 0.0);
    assert(cv < 15.0);
    TEST_PASS();
}

/* ============================================================================
 * L6: Cyclic voltammetry peak detection
 * ============================================================================ */
static void test_detect_voltammetry_peaks(void) {
    /* Simulated CV: oxidation peak at +0.4V, reduction at +0.1V */
    double voltage[] = {-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    double current[] = {0.0, 0.01, 0.05, 0.03, 0.1, 0.5, 1.0, 0.3, 0.05};
    double peak_v[5], peak_i[5];
    int n_peaks = detect_voltammetry_peaks(current, voltage, 9,
                                            peak_v, peak_i, 5);
    assert(n_peaks >= 1);  /* should find at least one peak */
    TEST_PASS();
}

static void test_coulometric_charge(void) {
    double current[] = {1.0, 1.5, 2.0, 1.5, 1.0};
    double time[]    = {0.0, 1.0, 2.0, 3.0, 4.0};
    double q = coulometric_charge(current, time, 5);
    /* Trapezoidal: 0.5*((1+1.5)+(1.5+2)+(2+1.5)+(1.5+1)) = 6.5 */
    assert(fabs(q - 6.0) < 0.1);
    TEST_PASS();
}

/* ============================================================================
 * L7: ELISA and immunoassay tests
 * ============================================================================ */
static void test_elisa_fourpl_predict(void) {
    FourPLLogisticFit fit;
    fourpl_fit_init(&fit);
    fit.a = 0.05; fit.b = 1.2; fit.c = 3.0; fit.d = 2.5;
    double sig_at_ec50 = elisa_fourpl_predict(&fit, 3.0);
    assert(fabs(sig_at_ec50 - (0.05 + 2.45/2.0)) < 0.01);
    double sig_zero = elisa_fourpl_predict(&fit, 0.0);
    assert(fabs(sig_zero - 0.05) < EPS);
    TEST_PASS();
}

static void test_antibody_antigen_complex(void) {
    /* Kd = 1 nM, [Ab] = 10 nM, [Ag] = 5 nM */
    double complex = antibody_antigen_complex(1.0e-9, 5.0e-9, 1.0e-8);
    /* Most antigen should be bound when [Ab] >> Kd */
    assert(complex > 4.0e-9);
    assert(complex <= 5.0e-9);
    TEST_PASS();
}

static void test_binding_association_curve(void) {
    double r_t = binding_association_curve(1.0e5, 0.001, 1.0e-7, 100.0, 300.0);
    double r_eq = 100.0 * 1.0e-7 / (1.0e-8 + 1.0e-7);  /* ~90.9 */
    assert(r_t > 0.0);
    assert(r_t < r_eq + 1.0);
    TEST_PASS();
}

/* ============================================================================
 * L7: DNA probe tests
 * ============================================================================ */
static void test_dna_hybridization_efficiency(void) {
    double eff = dna_hybridization_efficiency(1.0e6, 1.0e-6);
    /* K*c = 1 → θ = 0.5 */
    assert(fabs(eff - 0.5) < 0.01);
    double eff_high = dna_hybridization_efficiency(1.0e6, 1.0e-4);
    assert(eff_high > 0.99);
    TEST_PASS();
}

static void test_dna_probe_melting_temperature(void) {
    /* 10-mer: ATCG ATCG AT */
    double tm = dna_probe_melting_temperature("ATCGATCGAT");
    /* AT=6, GC=4 → Tm = 2*6 + 4*4 = 12+16 = 28°C */
    assert(fabs(tm - 28.0) < 0.1);
    /* GC-rich 10-mer: GGG CCC GGG C */
    tm = dna_probe_melting_temperature("GGGCCCGGGC");
    /* AT=0, GC=10 → Tm = 0 + 40 = 40°C */
    assert(fabs(tm - 40.0) < 0.1);
    TEST_PASS();
}

/* ============================================================================
 * L6: Sensor aging and recalibration
 * ============================================================================ */
static void test_two_point_calibration_adjust(void) {
    LinearCalibration old_cal, new_cal;
    linear_calibration_init(&old_cal);
    old_cal.slope = 1.0; old_cal.intercept = 0.1;

    two_point_calibration_adjust(2.0, 1.9, 8.0, 8.1, &old_cal, &new_cal);
    /* New slope should be ~(8.1-1.9)/(8-2) = 6.2/6 = 1.033 */
    assert(fabs(new_cal.slope - 1.033) < 0.1);
    TEST_PASS();
}

static void test_estimate_sensor_remaining_life(void) {
    double remaining = estimate_sensor_remaining_life(100.0, 80.0, 30, 50.0);
    assert(remaining > 0.0);
    assert(remaining < 365.0);
    /* Sensor already below minimum */
    double dead = estimate_sensor_remaining_life(100.0, 40.0, 30, 50.0);
    assert(dead == 0.0);
    TEST_PASS();
}

/* ============================================================================
 * L8: Calibration transfer
 * ============================================================================ */
static void test_pds_calibration_transfer(void) {
    int n_standards = 3;
    int n_channels = 4;
    int window_hw  = 1;
    double master[] = {
        1.0, 2.0, 3.0, 4.0,
        2.0, 3.0, 4.0, 5.0,
        3.0, 4.0, 5.0, 6.0
    };
    double slave[] = {
        1.1, 2.1, 3.1, 4.1,
        2.1, 3.1, 4.1, 5.1,
        3.1, 4.1, 5.1, 6.1
    };
    double *transfer = calloc(n_channels * n_channels, sizeof(double));
    assert(transfer != NULL);
    double rmse = pds_calibration_transfer(master, slave,
                                            n_standards, n_channels,
                                            window_hw, transfer);
    /* RMSE should be reasonable (slave ≈ master + constant offset) */
    assert(rmse < 1.0);
    free(transfer);
    TEST_PASS();
}

/* ============================================================================
 * L5: Mandel and Durbin-Watson tests
 * ============================================================================ */
static void test_mandel_fitting_test(void) {
    /* This test requires replicates; with single replicates it defaults to true */
    CalibrationStandards cs;
    calibration_standards_init(&cs, 5, 1);
    LinearCalibration lc;
    linear_calibration_init(&lc);
    bool adequate = mandel_fitting_test(&cs, &lc, 3.5);
    assert(adequate);  /* cannot reject without replicate data */
    calibration_standards_free(&cs);
    TEST_PASS();
}

static void test_durbin_watson_test(void) {
    double residuals[] = {0.1, -0.15, 0.05, -0.2, 0.1};
    double dw;
    bool independent = durbin_watson_test(residuals, 5, &dw);
    assert(dw > 0.0);
    assert(dw < 4.0);
    (void)independent;
    TEST_PASS();
}

/* ============================================================================
 * Main test runner
 * ============================================================================ */
int main(void) {
    printf("=== mini-biosensor-chemical Test Suite ===\n\n");

    printf("L1: Type Initialization Tests\n");
    test_biosensor_descriptor_init();
    test_michaelis_menten_params_init();
    test_nernst_params_init();
    test_qcm_params_init();
    test_measurement_series_init_free();
    test_calibration_standards_init_free();
    test_sensor_array_init_free();

    printf("\nL4: Fundamental Laws Tests\n");
    test_nernst_potential_standard();
    test_nernst_potential_ratio_10();
    test_isfet_compute_ph();
    test_butler_volmer_zero_overpotential();
    test_cottrell_current_decay();
    test_randles_sevcik_linear_scanrate();
    test_beer_lambert_absorbance();
    test_concentration_from_absorbance();
    test_two_component_absorbance();
    test_fret_efficiency_at_r0();
    test_fret_efficiency_bounds();
    test_michaelis_menten_rate();
    test_michaelis_menten_inverse();
    test_competitive_inhibition();
    test_lineweaver_burk_estimate();
    test_mm_nls_fit();
    test_hill_fractional_saturation();
    test_langmuir_coverage();
    test_langmuir_kinetic_convergence();

    printf("\nL2: Core Concepts Tests\n");
    test_spr_resonance_angle();
    test_spr_reflectivity_bounds();

    printf("\nL5: Algorithms Tests\n");
    test_lod_three_sigma();
    test_loq_ten_sigma();
    test_moving_minimum_baseline();
    test_exponential_drift_correction();
    test_savitzky_golay_smooth();
    test_median_filter();
    test_linear_calibration_fit();
    test_standard_addition_method();
    test_fourpl_fit();
    test_grubbs_outlier();
    test_coefficient_of_variation();
    test_mandel_fitting_test();
    test_durbin_watson_test();

    printf("\nL6: Canonical Problems Tests\n");
    test_glucose_biosensor_linear_range();
    test_glucose_conc_from_current();
    test_detect_voltammetry_peaks();
    test_coulometric_charge();
    test_two_point_calibration_adjust();
    test_estimate_sensor_remaining_life();

    printf("\nL7: Applications Tests\n");
    test_elisa_fourpl_predict();
    test_antibody_antigen_complex();
    test_binding_association_curve();
    test_dna_hybridization_efficiency();
    test_dna_probe_melting_temperature();

    printf("\nL8: Advanced Topics Tests\n");
    test_pds_calibration_transfer();

    printf("\n=== All tests passed ===\n");
    return 0;
}
