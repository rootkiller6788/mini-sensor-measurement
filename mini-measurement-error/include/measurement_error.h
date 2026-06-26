/**
 * @file measurement_error.h
 * @brief L1 Definitions — Core measurement error types, categories, and fundamental structures
 *
 * Covers:
 *   L1 - Accuracy, precision, bias, random error, hysteresis, linearity,
 *        quantization error, zero drift, span drift, repeatability, reproducibility
 *   L2 - Systematic vs random errors, error classification, measurement chains
 *
 * Reference: ISO 5725 (Accuracy), JCGM 200:2012 (VIM), GUM JCGM 100:2008
 * Courses: MIT 6.630, ETH 227-0455, Berkeley EE117, Tsinghua 信号与系统
 */

#ifndef MEASUREMENT_ERROR_H
#define MEASUREMENT_ERROR_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

/* Math constants for strict C11 compliance (not in standard <math.h>) */
#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif
#ifndef M_SQRT2
#define M_SQRT2     1.41421356237309504880
#endif

/* ─── L1: Error Category Enumeration ────────────────────────────────────────
 * Each enum value represents a distinct error mechanism in measurement systems.
 * ISO 5725-1:1994 defines accuracy as the combination of trueness and precision.
 */

typedef enum {
    ERR_NONE            = 0,
    ERR_SYSTEMATIC_BIAS = 1,
    ERR_RANDOM_NOISE    = 2,
    ERR_HYSTERESIS      = 3,
    ERR_NONLINEARITY    = 4,
    ERR_QUANTIZATION    = 5,
    ERR_ZERO_DRIFT      = 6,
    ERR_SPAN_DRIFT      = 7,
    ERR_RESOLUTION      = 8,
    ERR_LOADING         = 9,
    ERR_ENVIRONMENTAL   = 10,
    ERR_AGING           = 11,
    ERR_CROSSTALK       = 12,
    ERR_COMMON_MODE     = 13,
    ERR_GROUND_LOOP     = 14,
    ERR_THERMAL_EMF     = 15
} ErrorCategory;

const char *error_category_name(ErrorCategory cat);

/* ─── L1: Measurement Type Enumeration ────────────────────────────────────── */
typedef enum {
    MEAS_DIRECT       = 0,
    MEAS_INDIRECT     = 1,
    MEAS_DIFFERENTIAL = 2,
    MEAS_RATIOMETRIC  = 3,
    MEAS_NULL         = 4,
    MEAS_SUBSTITUTION = 5
} MeasurementType;

const char *measurement_type_name(MeasurementType type);

/* ─── L1: Core Measurement Point ──────────────────────────────────────────── */
typedef struct {
    double  value;
    double  true_value;
    double  systematic_error;
    double  random_error;
    double  uncertainty;
    int64_t timestamp_us;
    uint8_t category_mask;
} MeasurementPoint;

/* ─── L1: Measurement Statistics (Descriptive) ────────────────────────────── */
typedef struct {
    double  mean;
    double  std_dev;
    double  variance;
    size_t  count;
    double  min_value;
    double  max_value;
    double  range;
    double  bias;
    double  precision;
    double  rms_error;
} MeasurementStats;

/* ─── L1: Error Budget Structure ──────────────────────────────────────────── */
#define MAX_ERROR_SOURCES 32

typedef struct {
    ErrorCategory category;
    double        magnitude;
    double        sensitivity;
    double        divisor;
    char          source_name[64];
} ErrorComponent;

typedef struct {
    ErrorComponent components[MAX_ERROR_SOURCES];
    size_t         num_components;
    double         total_systematic;
    double         total_random_rss;
    double         combined_uncertainty;
    double         expanded_uncertainty;
    double         coverage_factor;
} ErrorBudget;

/* ─── L1: Sensor Specification ────────────────────────────────────────────── */
typedef struct {
    double  full_scale_range;
    double  sensitivity;
    double  offset;
    double  nonlinearity_pct;
    double  hysteresis_pct;
    double  repeatability_pct;
    double  resolution;
    double  noise_rms;
    double  temp_coeff_zero;
    double  temp_coeff_span;
    double  long_term_stability;
    double  bandwidth_hz;
    double  input_impedance;
} SensorSpec;

/* ─── L2: Core Error Computation Functions ────────────────────────────────── */

/** Compute absolute error: e = x_measured - x_true */
double meas_absolute_error(double measured, double true_value);

/** Compute relative error: epsilon = (x_measured - x_true) / x_true */
double meas_relative_error(double measured, double true_value);

/** Compute percentage error: delta = 100 * |x_measured - x_true| / x_true */
double meas_percentage_error(double measured, double true_value);

/** Full-scale error (% of FSR) */
double meas_fs_error(double absolute_error, double full_scale_range);

/** Accuracy class per IEC 60051 */
double meas_accuracy_class(double max_percentage_error);

/** Repeatability per ISO 5725-2:2019 section 3.21 */
double meas_repeatability(const double *measurements, size_t n);

/** Reproducibility — between-laboratory precision */
double meas_reproducibility(const double *measurements, size_t n,
                             double inter_lab_std);

/** Hysteresis error from up-scale and down-scale readings */
double meas_hysteresis_error(double x_up, double x_down, double fsr);

/** Endpoint nonlinearity: max deviation from endpoint line as % FSR */
double meas_endpoint_nonlinearity(const double *input, const double *output,
                                   size_t n, double *max_dev);

/** ADC quantization error RMS: Q / sqrt(12) */
double meas_quantization_error_rms(double full_scale_range, unsigned int bits);

/** Loading error percentage due to finite instrument input impedance */
double meas_loading_error_pct(double source_impedance, double instrument_impedance);

/** Common-mode rejection error */
double meas_cmrr_error(double common_mode_voltage, double cmrr_db);

/** Thermal EMF error from Seebeck effect */
double meas_thermal_emf_error(double seebeck_coeff_uv_per_c,
                               double temp_difference_c);

/* ─── Error Budget ────────────────────────────────────────────────────────── */
void error_budget_init(ErrorBudget *budget, double coverage_factor);
int  error_budget_add(ErrorBudget *budget, ErrorCategory cat, double magnitude,
                      double sensitivity, double divisor, const char *name);
void error_budget_finalize(ErrorBudget *budget);

/* ─── Measurement Statistics ──────────────────────────────────────────────── */
void meas_stats_init(MeasurementStats *stats);
void meas_stats_update(MeasurementStats *stats, const double *data, size_t n,
                        double true_value);

#endif /* MEASUREMENT_ERROR_H */
