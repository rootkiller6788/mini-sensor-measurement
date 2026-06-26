#ifndef MEASUREMENT_APPLICATIONS_H
#define MEASUREMENT_APPLICATIONS_H

#include <stddef.h>
#include "measurement_error.h"
#include "measurement_uncertainty.h"

/**
 * @file measurement_applications.h
 * @brief L7-L8: Applied measurement error analysis for real-world sensors
 *
 * L7 Applications:
 *   - Thermocouple measurement system error budget
 *   - Strain gauge bridge measurement uncertainty
 *   - ADC data acquisition system error analysis
 *   - pH sensor calibration and drift analysis
 *   - Industrial process temperature measurement
 *
 * L8 Advanced Topics:
 *   - Monte Carlo uncertainty propagation (GUM Supplement 1)
 *   - Bayesian sensor fault detection
 *   - Time-varying drift modeling (Wiener process)
 *
 * Courses: MIT 6.630 (sensors lab), Stanford EE359, ETH 227-0455
 */

/* --- L7: Thermocouple Measurement Error Budget --- */

typedef struct {
    double seebeck_coeff_mv_per_c;     /* Type K: ~0.041 mV/C */
    double seebeck_uncertainty_pct;    /* typical 0.75% */
    double cjc_sensor_accuracy_c;      /* cold-junction sensor accuracy */
    double amplifier_offset_uv;        /* input offset voltage */
    double amplifier_gain_error_pct;   /* gain error */
    double amplifier_noise_uv_rms;     /* input-referred noise */
    double adc_resolution_bits;        /* ADC bits */
    double adc_inl_lsb;                /* ADC INL in LSB */
    double adc_dnl_lsb;                /* ADC DNL in LSB */
    double adc_reference_error_pct;    /* voltage reference accuracy */
    double emi_induced_error_uv;       /* estimated EMI pickup */
} ThermocoupleSystemSpec;

void thermocouple_spec_init(ThermocoupleSystemSpec *spec);
void thermocouple_error_budget(const ThermocoupleSystemSpec *spec,
                                double measured_temp_c,
                                UncertaintyBudget *budget);

/* --- L7: Strain Gauge Measurement Uncertainty --- */

typedef struct {
    double excitation_voltage;
    double excitation_accuracy_pct;
    double gauge_factor;
    double gf_uncertainty_pct;
    double bridge_imbalance_uv_per_v;
    double amplifier_gain;
    double amplifier_offset_uv;
    double amplifier_noise_uv_rms;
    double adc_bits;
    double adc_fsr_v;
} StrainGaugeSystemSpec;

void straingauge_spec_init(StrainGaugeSystemSpec *spec);
void straingauge_uncertainty_budget(const StrainGaugeSystemSpec *spec,
                                     double measured_strain,
                                     UncertaintyBudget *budget);
double straingauge_strain_from_voltage(double V_out, double V_ex,
                                        double gauge_factor);

/* --- L7: ADC/DAQ System Error Analysis --- */

typedef struct {
    unsigned int bits;
    double       fsr_voltage;
    double       inl_lsb;
    double       dnl_lsb;
    double       offset_error_lsb;
    double       gain_error_lsb;
    double       noise_rms_lsb;
    double       reference_accuracy_pct;
    double       reference_tempco_ppm_per_c;
    double       jitter_ps_rms;
    double       input_frequency_hz;
} ADCSpec;

void adc_spec_init(ADCSpec *spec);

/** Total ADC RMS noise in LSB including quantization, INL, thermal noise */
double adc_total_noise_lsb(const ADCSpec *spec);

/** ADC effective number of bits: ENOB = (SINAD - 1.76) / 6.02 */
double adc_enob(const ADCSpec *spec);

/** ADC signal-to-noise-and-distortion ratio (dB) */
double adc_sinad_db(const ADCSpec *spec);

/** ADC aperture jitter error: SNR_jitter = -20*log10(2*pi*f_in*t_jitter) */
double adc_jitter_snr_db(double input_freq_hz, double jitter_ps_rms);

/** ADC error budget: fill uncertainty components from ADC spec */
void adc_uncertainty_budget(const ADCSpec *spec, double input_voltage,
                             UncertaintyBudget *budget);

/* --- L7: pH Sensor Calibration and Drift --- */

typedef struct {
    double zero_point_mv;        /* mV at pH 7 (theoretical: 0 mV) */
    double slope_mv_per_ph;      /* Nernst slope (~59.16 mV/pH at 25C) */
    double slope_tempco_pct_per_c;
    double drift_mv_per_day;
    double noise_mv_rms;
    double input_impedance_ohm;  /* typically > 10^12 ohm */
} PHSensorSpec;

void ph_sensor_spec_init(PHSensorSpec *spec);

/** Convert pH to mV using Nernst equation: E = E0 + slope * (pH - 7) */
double ph_to_mv(double ph, const PHSensorSpec *spec);

/** Convert mV to pH */
double mv_to_ph(double mv, const PHSensorSpec *spec);

/** pH measurement error budget */
void ph_error_budget(const PHSensorSpec *spec, double ph, double temp_c,
                      double days_since_cal, UncertaintyBudget *budget);

/* --- L8: Monte Carlo Uncertainty Propagation (GUM Supplement 1) --- */

/**
 * Propagate uncertainties through an arbitrary measurement model
 * using Monte Carlo simulation.
 *
 * GUM Supplement 1 describes MCM for uncertainty evaluation when:
 *  - The measurement model is strongly nonlinear
 *  - Input distributions are non-Gaussian
 *  - The Welch-Satterthwaite approximation is invalid
 *
 * @param model_func  y = f(x1, ..., xn): takes input vector and returns output
 * @param means       mean values of n input quantities
 * @param stddevs     standard uncertainties of n input quantities
 * @param dists       distributions for each input
 * @param n_inputs    number of input quantities
 * @param n_samples   number of MC trials (recommended: 100000 to 1000000)
 * @param y_mean      [out] mean of output distribution
 * @param y_std       [out] standard uncertainty of output
 * @param y_ci_low    [out] 2.5th percentile
 * @param y_ci_high   [out] 97.5th percentile
 */
void mc_uncertainty_propagation(
    double (*model_func)(const double *inputs, void *context),
    const double *means, const double *stddevs,
    const UncertaintyDistribution *dists, size_t n_inputs,
    size_t n_samples, void *context,
    double *y_mean, double *y_std, double *y_ci_low, double *y_ci_high);

/* --- L8: Bayesian Sensor Fault Detection --- */

typedef struct {
    double prob_normal;       /* P(H0): probability sensor is normal */
    double prob_fault;        /* P(H1): probability sensor is faulty */
    double likelihood_normal_mean;
    double likelihood_normal_std;
    double likelihood_fault_mean;
    double likelihood_fault_std;
    double threshold;         /* posterior fault probability threshold */
} BayesianFaultDetector;

void bayes_fault_init(BayesianFaultDetector *det, double prior_fault_prob,
                       double normal_std, double fault_std, double threshold);

/**
 * Update fault probability using Bayes rule on new measurement.
 * Returns posterior P(H1 | measurement).
 */
double bayes_fault_update(BayesianFaultDetector *det, double measurement);

/* --- L8: Time-Varying Drift Model (Wiener Process) --- */

typedef struct {
    double drift_rate;        /* mu: drift rate per unit time */
    double diffusion;         /* sigma: volatility per sqrt(time) */
    double current_drift;     /* current drift state */
    double time_elapsed;      /* time since last update */
} DriftModel;

void drift_model_init(DriftModel *model, double drift_rate, double diffusion);

/** Evolve drift over dt seconds. Returns predicted total drift. */
double drift_model_evolve(DriftModel *model, double dt);

/** Probability that drift exceeds limit within time horizon */
double drift_exceedance_probability(const DriftModel *model,
                                     double limit, double horizon);

#endif /* MEASUREMENT_APPLICATIONS_H */
