#include "measurement_applications.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>

/* ─── L7: Thermocouple Measurement Error Budget ──────────────────────────
 *
 * A complete Type K thermocouple measurement chain includes:
 *   - Thermocouple itself (Seebeck coefficient uncertainty)
 *   - Cold junction compensation sensor (accuracy)
 *   - Amplifier (offset, gain error, noise)
 *   - ADC (quantization, INL, DNL)
 *   - Voltage reference (accuracy, tempco)
 *   - EMI pickup
 *
 * Each component contributes to the total measurement uncertainty.
 */

void thermocouple_spec_init(ThermocoupleSystemSpec *spec) {
    if (!spec) return;
    memset(spec, 0, sizeof(*spec));
    /* Default Type K characteristics */
    spec->seebeck_coeff_mv_per_c = 0.041;     /* ~41 uV/C at 25C */
    spec->seebeck_uncertainty_pct = 0.75;      /* Standard tolerance */
    spec->cjc_sensor_accuracy_c = 0.5;         /* e.g., LM35 or thermistor */
    spec->amplifier_offset_uv = 10.0;          /* e.g., AD8495: 5 uV typ */
    spec->amplifier_gain_error_pct = 0.2;      /* 0.2% typical */
    spec->amplifier_noise_uv_rms = 1.5;        /* 0.1-10 Hz noise */
    spec->adc_resolution_bits = 16;
    spec->adc_inl_lsb = 1.0;
    spec->adc_dnl_lsb = 0.5;
    spec->adc_reference_error_pct = 0.05;      /* precision ref */
    spec->emi_induced_error_uv = 1.0;         /* estimated */
}

void thermocouple_error_budget(const ThermocoupleSystemSpec *spec,
                                double measured_temp_c,
                                UncertaintyBudget *budget) {
    if (!spec || !budget) return;
    unc_budget_init(budget, 0.95);

    /* 1. Thermocouple Seebeck coefficient uncertainty */
    double V_measured = measured_temp_c * spec->seebeck_coeff_mv_per_c;
    double u_seebeck = V_measured * spec->seebeck_uncertainty_pct / 100.0;
    /* Convert voltage uncertainty to temperature: 1 mV = 24.4 C for Type K */
    double sensitivity_temp = 1.0 / spec->seebeck_coeff_mv_per_c;
    unc_add_type_b(budget, "TC Seebeck Coeff", DIST_NORMAL,
                    u_seebeck, sensitivity_temp);

    /* 2. Cold junction compensation sensor */
    unc_add_type_b(budget, "CJC Sensor", DIST_NORMAL,
                    spec->cjc_sensor_accuracy_c, 1.0);

    /* 3. Amplifier offset (referred to input) */
    double u_amp_offset = spec->amplifier_offset_uv * 1e-3; /* uV -> mV */
    unc_add_type_b(budget, "Amp Offset", DIST_UNIFORM,
                    u_amp_offset, sensitivity_temp);

    /* 4. Amplifier gain error */
    double u_amp_gain = V_measured * spec->amplifier_gain_error_pct / 100.0;
    unc_add_type_b(budget, "Amp Gain Error", DIST_NORMAL,
                    u_amp_gain, sensitivity_temp);

    /* 5. Amplifier noise */
    double u_amp_noise = spec->amplifier_noise_uv_rms * 1e-3;
    unc_add_type_b(budget, "Amp Noise", DIST_NORMAL,
                    u_amp_noise, sensitivity_temp);

    /* 6. ADC quantization */
    double V_fsr = 5.0; /* typical ADC full-scale */
    double lsb = V_fsr / pow(2.0, spec->adc_resolution_bits);
    double u_adc_quant = lsb / sqrt(12.0); /* uniform distribution RMS */
    unc_add_type_b(budget, "ADC Quantization", DIST_UNIFORM,
                    u_adc_quant * 1000.0, sensitivity_temp); /* V -> mV */

    /* 7. ADC INL */
    double u_adc_inl = spec->adc_inl_lsb * lsb;
    unc_add_type_b(budget, "ADC INL", DIST_UNIFORM,
                    u_adc_inl * 1000.0, sensitivity_temp);

    /* 8. Voltage reference error */
    double u_ref = V_fsr * spec->adc_reference_error_pct / 100.0
                   / sqrt(3.0); /* uniform distribution */
    unc_add_type_b(budget, "Vref Error", DIST_UNIFORM,
                    u_ref * 1000.0, sensitivity_temp);

    /* 9. EMI */
    double u_emi = spec->emi_induced_error_uv * 1e-3;
    unc_add_type_b(budget, "EMI Pickup", DIST_UNIFORM,
                    u_emi, sensitivity_temp);

    unc_budget_finalize(budget);
}

/* ─── L7: Strain Gauge Measurement Uncertainty ─────────────────────────── */

void straingauge_spec_init(StrainGaugeSystemSpec *spec) {
    if (!spec) return;
    memset(spec, 0, sizeof(*spec));
    spec->excitation_voltage = 5.0;
    spec->excitation_accuracy_pct = 0.05;
    spec->gauge_factor = 2.07;
    spec->gf_uncertainty_pct = 0.5;
    spec->bridge_imbalance_uv_per_v = 2.0;
    spec->amplifier_gain = 100.0;
    spec->amplifier_offset_uv = 50.0;
    spec->amplifier_noise_uv_rms = 2.0;
    spec->adc_bits = 24;
    spec->adc_fsr_v = 5.0;
}

/* Convert bridge output voltage to strain for quarter-bridge:
 * epsilon = 4 * V_out / (V_ex * GF)   (linear approximation)
 */
double straingauge_strain_from_voltage(double V_out, double V_ex,
                                        double gauge_factor) {
    if (fabs(V_ex) < DBL_EPSILON || fabs(gauge_factor) < DBL_EPSILON)
        return 0.0;
    return 4.0 * V_out / (V_ex * gauge_factor);
}

void straingauge_uncertainty_budget(const StrainGaugeSystemSpec *spec,
                                     double measured_strain,
                                     UncertaintyBudget *budget) {
    if (!spec || !budget) return;
    unc_budget_init(budget, 0.95);

    /* 1. Excitation voltage uncertainty */
    double u_vex = spec->excitation_voltage
                   * spec->excitation_accuracy_pct / 100.0 / sqrt(3.0);
    /* Sensitivity of strain to V_ex: d(strain)/d(V_ex) = -strain / V_ex */
    double sens_vex = -measured_strain / spec->excitation_voltage;
    unc_add_type_b(budget, "Excitation V", DIST_UNIFORM,
                    u_vex, sens_vex);

    /* 2. Gauge factor uncertainty */
    double u_gf = spec->gauge_factor * spec->gf_uncertainty_pct / 100.0;
    double sens_gf = -measured_strain / spec->gauge_factor;
    unc_add_type_b(budget, "Gauge Factor", DIST_NORMAL,
                    u_gf, sens_gf);

    /* 3. Bridge imbalance (offset) */
    double V_imbalance = spec->bridge_imbalance_uv_per_v * 1e-6
                         * spec->excitation_voltage;
    double strain_equiv = straingauge_strain_from_voltage(
        V_imbalance, spec->excitation_voltage, spec->gauge_factor);
    unc_add_type_b(budget, "Bridge Imbalance", DIST_UNIFORM,
                    strain_equiv, 1.0);

    /* 4. Amplifier offset */
    double V_amp_offset = spec->amplifier_offset_uv * 1e-6;
    double strain_amp_off = straingauge_strain_from_voltage(
        V_amp_offset / spec->amplifier_gain,
        spec->excitation_voltage, spec->gauge_factor);
    unc_add_type_b(budget, "Amp Offset", DIST_UNIFORM,
                    strain_amp_off, 1.0);

    /* 5. Amplifier noise */
    double V_amp_noise = spec->amplifier_noise_uv_rms * 1e-6;
    double strain_amp_noise = straingauge_strain_from_voltage(
        V_amp_noise / spec->amplifier_gain,
        spec->excitation_voltage, spec->gauge_factor);
    unc_add_type_b(budget, "Amp Noise", DIST_NORMAL,
                    strain_amp_noise, 1.0);

    /* 6. ADC quantization */
    double lsb = spec->adc_fsr_v / pow(2.0, spec->adc_bits);
    double u_adc = lsb / sqrt(12.0);
    double strain_adc = straingauge_strain_from_voltage(
        u_adc / spec->amplifier_gain,
        spec->excitation_voltage, spec->gauge_factor);
    unc_add_type_b(budget, "ADC Quantization", DIST_UNIFORM,
                    strain_adc, 1.0);

    unc_budget_finalize(budget);
}

/* ─── L7: ADC Error Analysis ───────────────────────────────────────────── */

void adc_spec_init(ADCSpec *spec) {
    if (!spec) return;
    memset(spec, 0, sizeof(*spec));
    spec->bits = 12;
    spec->fsr_voltage = 3.3;
    spec->inl_lsb = 1.0;
    spec->dnl_lsb = 0.5;
    spec->offset_error_lsb = 2.0;
    spec->gain_error_lsb = 4.0;
    spec->noise_rms_lsb = 0.3;
    spec->reference_accuracy_pct = 0.1;
    spec->reference_tempco_ppm_per_c = 20.0;
    spec->jitter_ps_rms = 1.0;
    spec->input_frequency_hz = 1000.0;
}

/* Total RMS noise from all sources, referred to LSB */
double adc_total_noise_lsb(const ADCSpec *spec) {
    if (!spec) return 0.0;

    /* Quantization noise: uniform -> RMS = Q / sqrt(12) */
    double q_rms = 1.0 / sqrt(12.0);

    /* INL contribution treated as additional noise */
    double inl_rms = spec->inl_lsb / sqrt(12.0);

    /* DNL contribution */
    double dnl_rms = spec->dnl_lsb / sqrt(12.0);

    /* Thermal noise */
    double thermal_rms = spec->noise_rms_lsb;

    /* Reference noise (proportional to signal, assume at FSR/2) */
    double ref_noise = 0.5 * pow(2.0, spec->bits)
                       * spec->reference_accuracy_pct / 100.0 / sqrt(3.0);

    /* RSS combination */
    double total_sq = q_rms * q_rms + inl_rms * inl_rms
                      + dnl_rms * dnl_rms + thermal_rms * thermal_rms
                      + ref_noise * ref_noise;
    return sqrt(total_sq);
}

/* Effective Number of Bits: ENOB = (SINAD_dB - 1.76) / 6.02
 * SINAD = Signal-to-Noise And Distortion ratio
 */
double adc_enob(const ADCSpec *spec) {
    double sinad = adc_sinad_db(spec);
    return (sinad - 1.76) / 6.02;
}

/* SINAD including quantization, thermal noise, and jitter */
double adc_sinad_db(const ADCSpec *spec) {
    if (!spec) return 0.0;

    /* Total noise power in LSB^2 */
    double noise_lsb2 = adc_total_noise_lsb(spec);
    noise_lsb2 *= noise_lsb2;

    /* Full-scale sinusoidal signal: amplitude = 2^(N-1) LSB
     * RMS = 2^(N-1) / sqrt(2) LSB
     * Signal power = (2^(N-1) / sqrt(2))^2 = 2^(2N-3) LSB^2 */
    double signal_rms = pow(2.0, spec->bits - 1) / M_SQRT2;
    double signal_power = signal_rms * signal_rms;

    if (noise_lsb2 < DBL_EPSILON) return INFINITY;

    double sinad_linear = signal_power / noise_lsb2;

    /* Jitter degradation */
    if (spec->jitter_ps_rms > 0.0 && spec->input_frequency_hz > 0.0) {
        double snr_jitter = adc_jitter_snr_db(spec->input_frequency_hz,
                                               spec->jitter_ps_rms);
        double jitter_linear = pow(10.0, -snr_jitter / 10.0);
        sinad_linear = 1.0 / (1.0 / sinad_linear + jitter_linear);
    }

    return 10.0 * log10(sinad_linear);
}

/* Aperture jitter SNR: SNR_j = -20 * log10(2 * pi * f_in * t_jitter_rms) */
double adc_jitter_snr_db(double input_freq_hz, double jitter_ps_rms) {
    if (input_freq_hz <= 0.0 || jitter_ps_rms <= 0.0) return INFINITY;
    double tj = jitter_ps_rms * 1e-12; /* ps -> s */
    return -20.0 * log10(2.0 * M_PI * input_freq_hz * tj);
}

void adc_uncertainty_budget(const ADCSpec *spec, double input_voltage,
                             UncertaintyBudget *budget) {
    if (!spec || !budget) return;
    unc_budget_init(budget, 0.95);

    double lsb_voltage = spec->fsr_voltage / pow(2.0, spec->bits);

    /* Quantization */
    unc_add_type_b(budget, "Quantization", DIST_UNIFORM,
                    lsb_voltage / sqrt(3.0), 1.0);

    /* INL */
    double inl_v = spec->inl_lsb * lsb_voltage;
    unc_add_type_b(budget, "INL", DIST_UNIFORM, inl_v, 1.0);

    /* DNL */
    double dnl_v = spec->dnl_lsb * lsb_voltage;
    unc_add_type_b(budget, "DNL", DIST_UNIFORM, dnl_v, 1.0);

    /* Offset */
    double offset_v = spec->offset_error_lsb * lsb_voltage;
    unc_add_type_b(budget, "Offset Error", DIST_UNIFORM, offset_v, 1.0);

    /* Gain error */
    double gain_v = spec->gain_error_lsb * lsb_voltage;
    unc_add_type_b(budget, "Gain Error", DIST_UNIFORM, gain_v, 1.0);

    /* Reference */
    double ref_v = input_voltage * spec->reference_accuracy_pct / 100.0
                   / sqrt(3.0);
    unc_add_type_b(budget, "Reference", DIST_UNIFORM, ref_v, 1.0);

    /* Thermal noise */
    double thermal_v = spec->noise_rms_lsb * lsb_voltage;
    unc_add_type_b(budget, "Thermal Noise", DIST_NORMAL, thermal_v, 1.0);

    unc_budget_finalize(budget);
}

/* ─── L7: pH Sensor ──────────────────────────────────────────────────────
 *
 * pH electrode follows the Nernst equation:
 *   E = E0 + (RT / nF) * ln(10) * (pH - 7)
 * At 25C: slope = 59.16 mV/pH
 *
 * Temperature dependence:
 *   slope(T) = slope(25C) * (T + 273.15) / 298.15
 */

void ph_sensor_spec_init(PHSensorSpec *spec) {
    if (!spec) return;
    memset(spec, 0, sizeof(*spec));
    spec->zero_point_mv = 0.0;           /* ideal: 0 mV at pH 7 */
    spec->slope_mv_per_ph = 59.16;       /* Nernst slope at 25C */
    spec->slope_tempco_pct_per_c = 0.336; /* 100/T(K): ~0.336%/C at 25C */
    spec->drift_mv_per_day = 0.5;        /* typical aging drift */
    spec->noise_mv_rms = 0.1;            /* pick-up and thermal noise */
    spec->input_impedance_ohm = 1e12;    /* high impedance for glass electrode */
}

double ph_to_mv(double ph, const PHSensorSpec *spec) {
    if (!spec) return 0.0;
    return spec->zero_point_mv + spec->slope_mv_per_ph * (ph - 7.0);
}

double mv_to_ph(double mv, const PHSensorSpec *spec) {
    if (!spec || fabs(spec->slope_mv_per_ph) < DBL_EPSILON) return 7.0;
    return 7.0 + (mv - spec->zero_point_mv) / spec->slope_mv_per_ph;
}

void ph_error_budget(const PHSensorSpec *spec, double ph, double temp_c,
                      double days_since_cal, UncertaintyBudget *budget) {
    if (!spec || !budget) return;
    unc_budget_init(budget, 0.95);

    /* Temperature-corrected slope */
    double slope_t = spec->slope_mv_per_ph
                     * (1.0 + spec->slope_tempco_pct_per_c / 100.0
                        * (temp_c - 25.0));

    /* 1. Zero-point uncertainty (typically 0.1 pH equivalent) */
    double u_zero_ph = 0.1;
    unc_add_type_b(budget, "Zero Point", DIST_NORMAL, u_zero_ph, 1.0);

    /* 2. Slope uncertainty (typically 1% after calibration) */
    double u_slope_ph = fabs(ph - 7.0) * 0.01; /* 1% of reading span */
    unc_add_type_b(budget, "Slope Uncertainty", DIST_NORMAL,
                    u_slope_ph, 1.0);

    /* 3. Temperature effect on slope */
    double u_temp_ph = fabs(ph - 7.0)
                       * spec->slope_tempco_pct_per_c / 100.0 * 0.2;
    unc_add_type_b(budget, "Temperature Effect", DIST_UNIFORM,
                    u_temp_ph, 1.0);

    /* 4. Drift since calibration */
    double drift_mv = spec->drift_mv_per_day * days_since_cal;
    double drift_ph = drift_mv / slope_t;
    unc_add_type_b(budget, "Drift", DIST_UNIFORM, drift_ph, 1.0);

    /* 5. Electrical noise */
    double noise_ph = spec->noise_mv_rms / slope_t;
    unc_add_type_b(budget, "Electrical Noise", DIST_NORMAL, noise_ph, 1.0);

    /* 6. Loading error from finite input impedance */
    /* For pH electrodes: input impedance > 10^12 ohm needed.
     * Assume source impedance ~100 Mohm for glass electrode. */
    double source_z = 100e6; /* 100 MOhm typical glass electrode */
    double loading_error = source_z / (source_z + spec->input_impedance_ohm);
    unc_add_type_b(budget, "Loading Error", DIST_UNIFORM,
                    loading_error * fabs(ph - 7.0), 1.0);

    unc_budget_finalize(budget);
}

/* ─── L8: Monte Carlo Uncertainty Propagation (GUM Supplement 1) ─────────
 *
 * Monte Carlo Method (MCM) for uncertainty evaluation:
 * 1. Define the measurement model Y = f(X1, ..., Xn)
 * 2. Assign PDFs to each input quantity
 * 3. Draw M random samples from each input PDF
 * 4. Evaluate the model for each sample
 * 5. Compute summary statistics from the output distribution
 *
 * This handles nonlinear models and non-Gaussian distributions
 * where the GUM uncertainty framework (first-order Taylor) breaks down.
 */

/* Simple LCG pseudo-random number generator */
static unsigned long lcg_state = 123456789UL;

static void lcg_seed(unsigned long seed) {
    lcg_state = seed;
}

static double lcg_rand(void) {
    lcg_state = lcg_state * 1103515245UL + 12345UL;
    return (double)(lcg_state & 0x7FFFFFFFUL) / 2147483648.0;
}

/* Box-Muller transform: generate two independent N(0,1) samples */
static double box_muller_normal(double *second) {
    double u1, u2, r;
    do {
        u1 = 2.0 * lcg_rand() - 1.0;
        u2 = 2.0 * lcg_rand() - 1.0;
        r = u1 * u1 + u2 * u2;
    } while (r >= 1.0 || r < 1e-15);
    double factor = sqrt(-2.0 * log(r) / r);
    if (second) *second = u2 * factor;
    return u1 * factor;
}

/* Compare doubles for sorting */
static int cmp_double(const void *a, const void *b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

void mc_uncertainty_propagation(
    double (*model_func)(const double *inputs, void *context),
    const double *means, const double *stddevs,
    const UncertaintyDistribution *dists, size_t n_inputs,
    size_t n_samples, void *context,
    double *y_mean, double *y_std, double *y_ci_low, double *y_ci_high) {

    if (!model_func || !means || !stddevs || !dists || n_inputs == 0
        || n_samples < 100) {
        if (y_mean)     *y_mean     = 0.0;
        if (y_std)      *y_std      = 0.0;
        if (y_ci_low)   *y_ci_low   = 0.0;
        if (y_ci_high)  *y_ci_high  = 0.0;
        return;
    }

    double *outputs = (double*)malloc(n_samples * sizeof(double));
    double *inputs  = (double*)malloc(n_inputs * sizeof(double));
    if (!outputs || !inputs) {
        free(outputs); free(inputs);
        if (y_mean) *y_mean = 0.0;
        return;
    }

    lcg_seed(42UL);

    double sum_y = 0.0, sum_y2 = 0.0;
    for (size_t s = 0; s < n_samples; s++) {
        /* Generate random input vector */
        for (size_t i = 0; i < n_inputs; i++) {
            double z;
            switch (dists[i]) {
                case DIST_NORMAL: {
                    z = box_muller_normal(NULL);
                    inputs[i] = means[i] + z * stddevs[i];
                    break;
                }
                case DIST_UNIFORM: {
                    double u = lcg_rand();
                    inputs[i] = means[i] + (2.0 * u - 1.0) * sqrt(3.0)
                                * stddevs[i];
                    break;
                }
                case DIST_TRIANGULAR: {
                    double u = lcg_rand();
                    if (u < 0.5) {
                        inputs[i] = means[i] + (sqrt(2.0 * u) - 1.0)
                                    * sqrt(6.0) * stddevs[i];
                    } else {
                        inputs[i] = means[i] + (1.0 - sqrt(2.0 * (1.0 - u)))
                                    * sqrt(6.0) * stddevs[i];
                    }
                    break;
                }
                default: {
                    z = box_muller_normal(NULL);
                    inputs[i] = means[i] + z * stddevs[i];
                }
            }
        }

        /* Evaluate model */
        double y = model_func(inputs, context);
        outputs[s] = y;
        sum_y  += y;
        sum_y2 += y * y;
    }

    double N = (double)n_samples;
    double mean_y = sum_y / N;
    double var_y  = (sum_y2 - sum_y * sum_y / N) / (N - 1.0);
    double std_y  = sqrt(var_y > 0.0 ? var_y : 0.0);

    /* Compute percentiles by sorting */
    qsort(outputs, n_samples, sizeof(double), cmp_double);
    size_t idx_lo = (size_t)(0.025 * (double)n_samples);
    size_t idx_hi = (size_t)(0.975 * (double)n_samples);
    if (idx_lo >= n_samples) idx_lo = 0;
    if (idx_hi >= n_samples) idx_hi = n_samples - 1;

    if (y_mean)    *y_mean    = mean_y;
    if (y_std)     *y_std     = std_y;
    if (y_ci_low)  *y_ci_low  = outputs[idx_lo];
    if (y_ci_high) *y_ci_high = outputs[idx_hi];

    free(outputs);
    free(inputs);
}

/* ─── L8: Bayesian Sensor Fault Detection ─────────────────────────────────
 *
 * Implements a simple Bayesian hypothesis test for sensor faults.
 * H0: sensor is operating normally (measurement ~ N(mu_normal, sigma_normal^2))
 * H1: sensor is faulty            (measurement ~ N(mu_fault, sigma_fault^2))
 *
 * Using Bayes rule:
 *   P(H1 | data) = P(data | H1) * P(H1) / P(data)
 *   where P(data) = P(data|H0)*P(H0) + P(data|H1)*P(H1)
 *
 * The prior P(H1) reflects the expected sensor reliability.
 */

void bayes_fault_init(BayesianFaultDetector *det, double prior_fault_prob,
                       double normal_std, double fault_std, double threshold) {
    if (!det) return;
    det->prob_normal = 1.0 - prior_fault_prob;
    det->prob_fault  = prior_fault_prob;
    det->likelihood_normal_mean = 0.0;   /* zero-mean under normal operation */
    det->likelihood_normal_std  = normal_std;
    det->likelihood_fault_mean  = 0.0;   /* zero-mean fault */
    det->likelihood_fault_std   = fault_std;
    det->threshold = threshold;
}

double bayes_fault_update(BayesianFaultDetector *det, double measurement) {
    if (!det) return 0.0;

    /* Likelihood under H0 */
    double z0 = measurement / det->likelihood_normal_std;
    double lh0 = exp(-0.5 * z0 * z0) / (det->likelihood_normal_std
                                         * sqrt(2.0 * M_PI));

    /* Likelihood under H1 */
    double z1 = measurement / det->likelihood_fault_std;
    double lh1 = exp(-0.5 * z1 * z1) / (det->likelihood_fault_std
                                         * sqrt(2.0 * M_PI));

    /* Bayes update */
    double evidence = lh0 * det->prob_normal + lh1 * det->prob_fault;
    if (evidence < DBL_EPSILON) return det->prob_fault;

    double posterior = lh1 * det->prob_fault / evidence;
    det->prob_fault  = posterior;
    det->prob_normal = 1.0 - posterior;

    return posterior;
}

/* ─── L8: Time-Varying Drift Model (Wiener Process) ───────────────────────
 *
 * Models sensor drift as a Wiener process (Brownian motion with drift):
 *
 *   dX_t = mu * dt + sigma * dW_t
 *
 * where:
 *   mu      = drift rate (deterministic trend)
 *   sigma   = diffusion coefficient (random component)
 *   dW_t    = Wiener increment ~ N(0, dt)
 *
 * After time dt:
 *   X_{t+dt} ~ N(X_t + mu*dt, sigma^2 * dt)
 *
 * This model is used for predicting sensor aging, frequency drift in
 * oscillators, and long-term stability of references.
 */

void drift_model_init(DriftModel *model, double drift_rate, double diffusion) {
    if (!model) return;
    model->drift_rate = drift_rate;
    model->diffusion = diffusion;
    model->current_drift = 0.0;
    model->time_elapsed = 0.0;
}

double drift_model_evolve(DriftModel *model, double dt) {
    if (!model || dt < 0.0) return model ? model->current_drift : 0.0;

    /* Generate Wiener increment */
    double z = box_muller_normal(NULL);
    double dW = z * sqrt(dt);

    model->current_drift += model->drift_rate * dt
                            + model->diffusion * dW;
    model->time_elapsed += dt;

    return model->current_drift;
}

/* Probability that drift exceeds limit within time horizon.
 * For a Wiener process with drift mu and diffusion sigma:
 *   P( max_{0<=s<=T} X_s >= L ) = 1 - Phi((L - mu*T)/(sigma*sqrt(T)))
 *                                   + exp(2*mu*L/sigma^2)
 *                                     * Phi((-L - mu*T)/(sigma*sqrt(T)))
 *
 * This is the first-passage probability for Brownian motion with drift.
 */
double drift_exceedance_probability(const DriftModel *model,
                                     double limit, double horizon) {
    if (!model || horizon <= 0.0 || limit <= 0.0) return 0.0;

    double mu = model->drift_rate;
    double sigma = model->diffusion;
    double x0 = model->current_drift;

    /* Remaining distance to limit */
    double d = limit - x0;

    /* Compute using reflection principle / Bachelier-Levy formula */
    double sigma_sqrt_T = sigma * sqrt(horizon);
    if (sigma_sqrt_T < DBL_EPSILON) {
        /* Pure drift case */
        return (mu * horizon >= d) ? 1.0 : 0.0;
    }

    double z1 = (d - mu * horizon) / sigma_sqrt_T;
    double z2 = (-d - mu * horizon) / sigma_sqrt_T;

    /* Use error function for Phi */
    double Phi1 = 0.5 * (1.0 + erf(z1 / M_SQRT2));
    double Phi2 = 0.5 * (1.0 + erf(z2 / M_SQRT2));

    double exceed_prob = 1.0 - Phi1
                         + exp(2.0 * mu * d / (sigma * sigma)) * Phi2;

    if (exceed_prob < 0.0) exceed_prob = 0.0;
    if (exceed_prob > 1.0) exceed_prob = 1.0;

    return exceed_prob;
}
