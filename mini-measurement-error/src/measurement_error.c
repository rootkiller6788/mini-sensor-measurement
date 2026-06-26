#include "measurement_error.h"
#include <math.h>
#include <string.h>
#include <float.h>

/* ─── L1: Error Category Names ─────────────────────────────────────────── */

const char *error_category_name(ErrorCategory cat) {
    static const char *names[] = {
        [ERR_NONE]            = "None",
        [ERR_SYSTEMATIC_BIAS] = "Systematic Bias",
        [ERR_RANDOM_NOISE]    = "Random Noise",
        [ERR_HYSTERESIS]      = "Hysteresis",
        [ERR_NONLINEARITY]    = "Nonlinearity",
        [ERR_QUANTIZATION]    = "Quantization",
        [ERR_ZERO_DRIFT]      = "Zero Drift",
        [ERR_SPAN_DRIFT]      = "Span Drift",
        [ERR_RESOLUTION]      = "Resolution",
        [ERR_LOADING]         = "Loading",
        [ERR_ENVIRONMENTAL]   = "Environmental",
        [ERR_AGING]           = "Aging",
        [ERR_CROSSTALK]       = "Crosstalk",
        [ERR_COMMON_MODE]     = "Common Mode",
        [ERR_GROUND_LOOP]     = "Ground Loop",
        [ERR_THERMAL_EMF]     = "Thermal EMF"
    };
    if (cat <= ERR_THERMAL_EMF) return names[cat];
    return "Unknown";
}

const char *measurement_type_name(MeasurementType type) {
    static const char *names[] = {
        [MEAS_DIRECT]       = "Direct",
        [MEAS_INDIRECT]     = "Indirect",
        [MEAS_DIFFERENTIAL] = "Differential",
        [MEAS_RATIOMETRIC]  = "Ratiometric",
        [MEAS_NULL]         = "Null-Balance",
        [MEAS_SUBSTITUTION] = "Substitution"
    };
    if (type <= MEAS_SUBSTITUTION) return names[type];
    return "Unknown";
}

/* ─── L2: Core Error Functions ────────────────────────────────────────── */

double meas_absolute_error(double measured, double true_value) {
    return measured - true_value;
}

double meas_relative_error(double measured, double true_value) {
    if (fabs(true_value) < DBL_EPSILON) {
        return (fabs(measured) < DBL_EPSILON) ? 0.0 : INFINITY;
    }
    return (measured - true_value) / true_value;
}

double meas_percentage_error(double measured, double true_value) {
    if (fabs(true_value) < DBL_EPSILON) {
        return (fabs(measured) < DBL_EPSILON) ? 0.0 : INFINITY;
    }
    return 100.0 * fabs(measured - true_value) / fabs(true_value);
}

double meas_fs_error(double absolute_error, double full_scale_range) {
    if (fabs(full_scale_range) < DBL_EPSILON) return INFINITY;
    return 100.0 * fabs(absolute_error) / fabs(full_scale_range);
}

double meas_accuracy_class(double max_percentage_error) {
    static const double classes[] = {0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.5, 5.0};
    int n = sizeof(classes) / sizeof(classes[0]);
    for (int i = 0; i < n; i++) {
        if (max_percentage_error <= classes[i]) return classes[i];
    }
    return 5.0;
}

/* ─── L1: Internal Statistics ──────────────────────────────────────────── */

static double sample_mean(const double *x, size_t n) {
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += x[i];
    return (n > 0) ? sum / (double)n : 0.0;
}

static double sample_stddev(const double *x, size_t n, double mean) {
    if (n < 2) return 0.0;
    double sum_sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        double d = x[i] - mean;
        sum_sq += d * d;
    }
    return sqrt(sum_sq / (double)(n - 1));
}

/* ─── L1: Repeatability (ISO 5725-2:2019) ────────────────────────────────
 *
 * r = t_{95%, n-1} * sqrt(2) * s_r
 * where s_r is the within-laboratory repeatability standard deviation.
 * The t-value is the two-tailed Student t for 95% confidence with n-1 df.
 */

double meas_repeatability(const double *measurements, size_t n) {
    if (n < 2) return 0.0;
    double mean = sample_mean(measurements, n);
    double s_r = sample_stddev(measurements, n, mean);
    /* Small-sample Student-t approximation */
    double t_val = 1.96;  /* Large-sample limit */
    if (n <= 30) {
        t_val = 2.0 + 2.0 / sqrt((double)n);
    }
    return t_val * M_SQRT2 * s_r;
}

/* ─── L1: Reproducibility ────────────────────────────────────────────────
 *
 * R = t_{95%} * sqrt(2) * s_R
 * where s_R^2 = s_r^2 + s_L^2
 * s_L = between-laboratory standard deviation
 */

double meas_reproducibility(const double *measurements, size_t n,
                             double inter_lab_std) {
    if (n < 2) return 0.0;
    double mean = sample_mean(measurements, n);
    double s_r = sample_stddev(measurements, n, mean);
    double s_R = sqrt(s_r * s_r + inter_lab_std * inter_lab_std);
    double t_val = (n <= 30) ? (2.0 + 2.0 / sqrt((double)n)) : 1.96;
    return t_val * M_SQRT2 * s_R;
}

/* ─── L1: Hysteresis Error ───────────────────────────────────────────────
 *
 * hyst_frac = |x_down - x_up| / FSR
 * Returns fractional hysteresis (multiply by 100 for percentage).
 */

double meas_hysteresis_error(double x_up, double x_down, double fsr) {
    if (fabs(fsr) < DBL_EPSILON) return INFINITY;
    return fabs(x_down - x_up) / fabs(fsr);
}

/* ─── L2: Endpoint Nonlinearity ──────────────────────────────────────────
 *
 * Fits a straight line through (input[0], output[0]) and
 * (input[n-1], output[n-1]), then finds maximum vertical deviation
 * from that line. Returns nonlinearity as % of output FSR.
 */

double meas_endpoint_nonlinearity(const double *input, const double *output,
                                   size_t n, double *max_dev) {
    if (n < 2) {
        if (max_dev) *max_dev = 0.0;
        return 0.0;
    }
    double x0 = input[0], y0 = output[0];
    double x1 = input[n-1], y1 = output[n-1];
    double dx = x1 - x0;
    if (fabs(dx) < DBL_EPSILON) {
        if (max_dev) *max_dev = 0.0;
        return 0.0;
    }
    double slope = (y1 - y0) / dx;
    double intercept = y0 - slope * x0;
    double fsr = output[n-1] - output[0];
    if (fabs(fsr) < DBL_EPSILON) {
        if (max_dev) *max_dev = 0.0;
        return 0.0;
    }
    double max_deviation = 0.0;
    for (size_t i = 0; i < n; i++) {
        double y_line = slope * input[i] + intercept;
        double dev = fabs(output[i] - y_line);
        if (dev > max_deviation) max_deviation = dev;
    }
    if (max_dev) *max_dev = max_deviation;
    return 100.0 * max_deviation / fabs(fsr);
}

/* ─── L2: ADC Quantization Error ─────────────────────────────────────────
 *
 * Ideal N-bit ADC: Q (LSB) = FSR / 2^N
 * Quantization error is uniformly distributed over [-Q/2, Q/2].
 * RMS of uniform distribution: sigma = Q / sqrt(12)
 */

double meas_quantization_error_rms(double full_scale_range, unsigned int bits) {
    if (bits == 0) return INFINITY;
    double lsb = full_scale_range / (double)(1ULL << bits);
    return lsb / sqrt(12.0);
}

/* ─── L2: Loading Error ─────────────────────────────────────────────────
 *
 * When a measuring instrument with finite input impedance Z_in is
 * connected to a source with impedance Z_source, the voltage divider
 * effect introduces loading error.
 *
 * V_measured = V_actual * Z_in / (Z_source + Z_in)
 * error_pct = Z_source / (Z_source + Z_in) * 100
 */

double meas_loading_error_pct(double source_impedance, double instrument_impedance) {
    if (fabs(source_impedance + instrument_impedance) < DBL_EPSILON)
        return 0.0;
    return 100.0 * source_impedance / (source_impedance + instrument_impedance);
}

/* ─── L2: Common-Mode Rejection Error ────────────────────────────────────
 *
 * CMRR(dB) = 20 * log10(A_dm / A_cm)
 * Input-referred error = V_cm / 10^(CMRR_dB / 20)
 */

double meas_cmrr_error(double common_mode_voltage, double cmrr_db) {
    double cmrr_linear = pow(10.0, cmrr_db / 20.0);
    if (cmrr_linear < DBL_EPSILON) return INFINITY;
    return fabs(common_mode_voltage) / cmrr_linear;
}

/* ─── L2: Thermal EMF Error ──────────────────────────────────────────────
 *
 * Seebeck effect: when two dissimilar metals form a junction with a
 * temperature gradient, a thermoelectric voltage is generated.
 *
 * V_thermal (uV) = alpha * delta_T
 * alpha = Seebeck coefficient in uV/degree_C
 */

double meas_thermal_emf_error(double seebeck_coeff_uv_per_c,
                               double temp_difference_c) {
    return seebeck_coeff_uv_per_c * temp_difference_c;
}

/* ─── L1: Error Budget Construction ────────────────────────────────────── */

void error_budget_init(ErrorBudget *budget, double coverage_factor) {
    if (!budget) return;
    memset(budget, 0, sizeof(*budget));
    budget->coverage_factor = coverage_factor;
}

int error_budget_add(ErrorBudget *budget, ErrorCategory cat, double magnitude,
                     double sensitivity, double divisor, const char *name) {
    if (!budget || budget->num_components >= MAX_ERROR_SOURCES) return -1;
    size_t i = budget->num_components;
    budget->components[i].category = cat;
    budget->components[i].magnitude = fabs(magnitude);
    budget->components[i].sensitivity = sensitivity;
    budget->components[i].divisor = (divisor > 0.0) ? divisor : 1.0;
    if (name) {
        strncpy(budget->components[i].source_name, name, 63);
        budget->components[i].source_name[63] = '\0';
    } else {
        budget->components[i].source_name[0] = '\0';
    }
    budget->num_components++;
    return 0;
}

/* Finalize: systematic = algebraic sum, random = RSS.
 * u_c = sqrt( u_sys^2 + u_rand_rss^2 )
 * U = k * u_c
 */

void error_budget_finalize(ErrorBudget *budget) {
    if (!budget) return;
    double sum_sys = 0.0;
    double sum_rss_sq = 0.0;
    for (size_t i = 0; i < budget->num_components; i++) {
        double u_i = budget->components[i].magnitude
                     / budget->components[i].divisor;
        double c_i = budget->components[i].sensitivity;
        double contrib = c_i * u_i;
        switch (budget->components[i].category) {
            case ERR_SYSTEMATIC_BIAS:
            case ERR_ZERO_DRIFT:
            case ERR_SPAN_DRIFT:
            case ERR_LOADING:
                sum_sys += contrib;
                break;
            default:
                sum_rss_sq += contrib * contrib;
                break;
        }
    }
    budget->total_systematic = sum_sys;
    budget->total_random_rss = sqrt(sum_rss_sq);
    budget->combined_uncertainty = sqrt(sum_sys * sum_sys + sum_rss_sq);
    budget->expanded_uncertainty = budget->coverage_factor
                                    * budget->combined_uncertainty;
}

/* ─── L1: Measurement Statistics ────────────────────────────────────────── */

void meas_stats_init(MeasurementStats *stats) {
    if (!stats) return;
    memset(stats, 0, sizeof(*stats));
    stats->min_value = INFINITY;
    stats->max_value = -INFINITY;
}

void meas_stats_update(MeasurementStats *stats, const double *data, size_t n,
                        double true_value) {
    if (!stats || !data || n == 0) return;
    double sum = 0.0, sum_sq = 0.0;
    double min_val = data[0], max_val = data[0];
    for (size_t i = 0; i < n; i++) {
        sum += data[i];
        sum_sq += data[i] * data[i];
        if (data[i] < min_val) min_val = data[i];
        if (data[i] > max_val) max_val = data[i];
    }
    stats->count = n;
    stats->mean = sum / (double)n;
    if (n > 1) {
        stats->variance = (sum_sq - sum * sum / (double)n) / (double)(n - 1);
    } else {
        stats->variance = 0.0;
    }
    stats->std_dev = sqrt(stats->variance);
    stats->min_value = min_val;
    stats->max_value = max_val;
    stats->range = max_val - min_val;
    stats->bias = stats->mean - true_value;
    stats->precision = 2.0 * stats->std_dev;
    stats->rms_error = sqrt(stats->bias * stats->bias
                             + stats->std_dev * stats->std_dev);
}
