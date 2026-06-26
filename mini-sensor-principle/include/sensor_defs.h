/**
 * @file sensor_defs.h
 * @brief Core sensor type definitions, measurement terms, and classification
 *
 * Level: L1 — Definitions
 *
 * Key Definitions:
 *   - Sensitivity: S = dY/dX — output change per unit input change
 *   - Resolution: smallest detectable input change (limited by noise floor)
 *   - Accuracy: closeness of measured value to true value
 *   - Precision: repeatability of measurements under unchanged conditions
 *   - Hysteresis: maximum difference in output when input is approached
 *         from opposite directions (ascending vs descending)
 *   - Linearity: deviation from best-fit straight line, often as %FS
 *   - Drift: slow change in output over time with constant input
 *   - Span: difference between maximum and minimum measurable input
 *   - Offset/Zero error: output when input is zero
 *   - Response time: time to reach 63.2% (1 - 1/e) of final value
 *   - Noise Equivalent Power (NEP): input signal power yielding SNR=1
 *   - Specific Detectivity D*: normalized detectivity per unit area & bandwidth
 *
 * References:
 *   - Fraden, "Handbook of Modern Sensors", 5th Ed, Ch. 1-2
 *   - Webster, "Medical Instrumentation", 4th Ed, Ch. 2
 *   - Bentley, "Principles of Measurement Systems", 4th Ed, Ch. 1-3
 */

#ifndef SENSOR_DEFS_H
#define SENSOR_DEFS_H

#include <stdint.h>
#include <stddef.h>

/* ──── L1: Sensor Classification & Type Hierarchy ──── */

/** Sensor energy conversion classification */
typedef enum {
    SENSOR_ENERGY_PASSIVE = 0,   /* Modulating: requires external excitation
                                    (strain gauge, RTD, thermistor, LVDT) */
    SENSOR_ENERGY_ACTIVE  = 1,   /* Self-generating: produces own signal
                                    (thermocouple, piezoelectric, photovoltaic) */
} sensor_energy_type_t;

/** Physical domain of measurement (SI base + common derived) */
typedef enum {
    SENSOR_DOMAIN_TEMPERATURE = 0,
    SENSOR_DOMAIN_MECHANICAL  = 1,  /* force, pressure, strain, displacement */
    SENSOR_DOMAIN_ELECTRICAL  = 2,  /* voltage, current, charge, resistance */
    SENSOR_DOMAIN_MAGNETIC    = 3,
    SENSOR_DOMAIN_OPTICAL     = 4,  /* IR, visible, UV */
    SENSOR_DOMAIN_CHEMICAL    = 5,
    SENSOR_DOMAIN_ACOUSTIC    = 6,
    SENSOR_DOMAIN_RADIATION   = 7,  /* ionizing: alpha, beta, gamma, neutron */
    SENSOR_DOMAIN_BIOLOGICAL  = 8,
    SENSOR_DOMAIN_FLUIDIC     = 9,  /* flow, level, viscosity */
    SENSOR_DOMAIN_INERTIAL    = 10, /* acceleration, angular velocity */
} sensor_domain_t;

/** Sensor output signal type */
typedef enum {
    SENSOR_OUTPUT_ANALOG_VOLTAGE  = 0,
    SENSOR_OUTPUT_ANALOG_CURRENT  = 1,  /* 4-20 mA industrial standard */
    SENSOR_OUTPUT_RESISTANCE      = 2,
    SENSOR_OUTPUT_CAPACITANCE     = 3,
    SENSOR_OUTPUT_INDUCTANCE      = 4,
    SENSOR_OUTPUT_FREQUENCY       = 5,
    SENSOR_OUTPUT_DIGITAL         = 6,  /* I2C, SPI, UART, 1-Wire */
    SENSOR_OUTPUT_CHARGE          = 7,  /* piezoelectric */
    SENSOR_OUTPUT_PULSE           = 8,  /* frequency-encoded, PWM */
} sensor_output_type_t;

/** Measurement performance classification */
typedef enum {
    SENSOR_GRADE_INDUSTRIAL = 0,  /* 0.5% – 2% FS accuracy */
    SENSOR_GRADE_LABORATORY = 1,  /* 0.05% – 0.5% FS accuracy */
    SENSOR_GRADE_PRECISION  = 2,  /* < 0.05% FS accuracy */
    SENSOR_GRADE_CONSUMER   = 3,  /* 2% – 10% FS accuracy, low cost */
    SENSOR_GRADE_MILITARY   = 4,  /* rugged, wide temp range, moderate accuracy */
} sensor_grade_t;

/* ──── L1: Core Measurement Performance Structure ──── */

/**
 * @brief Complete static measurement performance specification.
 *
 * Theorem: Measurement Chain Error Budget (GUM, ISO/IEC Guide 98-3)
 *   Total expanded uncertainty U = k * sqrt(Σ u_i²)
 *   where u_i are Type A (statistical) and Type B (systematic) contributions.
 *   Coverage factor k=2 gives ~95% confidence interval for normal distribution.
 */
typedef struct {
    /* ── Transfer characteristic parameters ── */
    double   sensitivity;        /* [output_unit / input_unit] dY/dX at operating point */
    double   full_scale_input;   /* [input_unit] maximum rated input */
    double   full_scale_output;  /* [output_unit] corresponding output */

    /* ── Error parameters ── */
    double   offset_error;       /* [output_unit] output at zero input */
    double   offset_drift;       /* [output_unit / °C] thermal zero drift coefficient */
    double   gain_error;         /* [% FS] sensitivity deviation from nominal */
    double   gain_drift;         /* [ppm FS / °C] thermal gain drift coefficient */

    double   nonlinearity;       /* [% FS] max deviation from best-fit line
                                    (endpoint or least-squares fit) */
    double   hysteresis;         /* [% FS] max difference between up-scale
                                    and down-scale readings at same input */
    double   repeatability;      /* [% FS] max deviation under same conditions
                                    in short time interval (per ISO 5725) */

    double   resolution;         /* [input_unit] smallest detectable change,
                                    typically 3σ of noise floor */
    double   noise_rms;          /* [output_unit RMS] broadband output noise */

    /* ── Dynamic parameters ── */
    double   bandwidth;          /* [Hz] -3 dB cutoff frequency */
    double   rise_time;          /* [s] 10%-90% step response rise time
                                    ≈ 0.35/bandwidth for 1st-order system */
    double   settling_time;      /* [s] time to settle within ±1% of final value */
    double   sampling_rate_max;  /* [Hz] max sampling rate (digital output sensors) */

    /* ── Environmental ── */
    double   temp_range_min;     /* [°C] minimum operating temperature */
    double   temp_range_max;     /* [°C] maximum operating temperature */
    double   temp_coefficient;   /* [%/°C] overall temperature sensitivity coefficient */

    /* ── Long-term stability ── */
    double   long_term_drift;    /* [% FS / year] annual stability */
    double   mtbf_hours;         /* [hours] mean time between failures */

} sensor_specs_t;

/* ──── L1: Sensor Identification & Metadata ──── */

typedef struct {
    char     manufacturer[32];
    char     model[32];
    char     serial_number[32];
    uint32_t manufacturing_date; /* YYYYMMDD format */

    sensor_domain_t      domain;
    sensor_energy_type_t energy_type;
    sensor_output_type_t output_type;
    sensor_grade_t       grade;

    /* Interface details */
    double   excitation_required;  /* V or mA, 0 if self-generating */
    double   output_impedance;     /* [Ω] source impedance */
    double   input_impedance;      /* [Ω] load impedance seen by measurand */
    double   power_consumption;    /* [W] */

    double   mass_grams;           /* [g] physical mass */
    double   volume_cm3;           /* [cm³] */

} sensor_identity_t;

/* ──── L1: Measurement Data Structure ──── */

/** Single calibrated sensor reading with uncertainty */
typedef struct {
    double   value;           /* calibrated measurement value [input_unit] */
    double   uncertainty;     /* expanded uncertainty [input_unit], k=2 */
    double   raw_output;      /* raw sensor output [output_unit] before calibration */
    uint64_t timestamp_us;    /* measurement timestamp [µs] since epoch */
    uint8_t  valid;           /* 1 = measurement valid, 0 = error */
    uint8_t  saturation;      /* 1 = sensor saturated (input beyond range) */
    uint16_t sample_index;    /* sequential sample number */
} sensor_reading_t;

/**
 * @brief Generic sensor abstraction.
 *
 * This struct abstracts the common interface of all sensors enabling
 * polymorphic sensor handling in measurement systems.
 * Concrete sensors embed this and provide specific implementations.
 */
typedef struct sensor_t {
    sensor_identity_t identity;
    sensor_specs_t    specs;

    /* Calibration polynomial coefficients (up to 5th order):
     *   y_cal = c[0] + c[1]*x + c[2]*x² + c[3]*x³ + c[4]*x⁴ + c[5]*x⁵ */
    double cal_coeffs[6];
    uint8_t cal_order;  /* polynomial order (1=linear, 2=quadratic, etc.) */
    uint8_t cal_valid;  /* 1 = calibration coefficients are valid */

    /* Function pointers for polymorphic behavior */
    int   (*init)(struct sensor_t *s);
    int   (*read)(struct sensor_t *s, sensor_reading_t *out);
    int   (*calibrate)(struct sensor_t *s, const double *ref, const double *raw, size_t n);
    int   (*self_test)(struct sensor_t *s);
    void  (*sleep)(struct sensor_t *s);
    void  (*wake)(struct sensor_t *s);

    /* Private implementation data */
    void *priv;

} sensor_t;

/* ──── L1: Sensor Error Budget Structure ──── */

/**
 * @brief Measurement uncertainty budget per GUM (ISO/IEC Guide 98-3).
 *
 * Each uncertainty component characterized by:
 *   - Type A: evaluated by statistical analysis of repeated measurements
 *   - Type B: evaluated by other means (datasheet, calibration cert, experience)
 *
 * Combined standard uncertainty: u_c = sqrt(Σ u_i²)
 * Expanded uncertainty: U = k * u_c (k=2 for ~95% coverage)
 */
typedef struct {
    double type_a_repeatability;   /* [input_unit] std dev of repeated readings */
    double type_b_reference;       /* [input_unit] reference standard uncertainty */
    double type_b_resolution;      /* [input_unit] DMM / ADC quantization */
    double type_b_thermal;         /* [input_unit] thermal EMF, drift effects */
    double type_b_loading;         /* [input_unit] sensor loading of measurand */
    double type_b_noise;           /* [input_unit] electrical noise contribution */
    double type_b_calibration;     /* [input_unit] calibration curve fit residual */
    double combined_standard;      /* [input_unit] sqrt(Σ u_i²) */
    double expanded_uncertainty;   /* [input_unit] k * u_c, k=2 */
    size_t degrees_of_freedom;     /* effective degrees of freedom (Welch-Satterthwaite) */
} sensor_uncertainty_budget_t;

/* ──── L1: Sensor Noise Model Types ──── */

/** Fundamental noise mechanisms in sensors */
typedef enum {
    NOISE_THERMAL_JOHNSON  = 0,  /* Johnson-Nyquist: Vn = sqrt(4kTRB) */
    NOISE_SHOT             = 1,  /* Schottky: In = sqrt(2q I_dc B) */
    NOISE_FLICKER_1F       = 2,  /* 1/f noise: S(f) ∝ 1/f^α, α≈0.8-1.3 */
    NOISE_POPCORN_BURST    = 3,  /* Random telegraph noise (RTS) */
    NOISE_QUANTIZATION     = 4,  /* ADC: σ = q/√12 where q = FS/2^N */
    NOISE_AVALANCHE        = 5,  /* Excess noise in reverse-biased junctions */
    NOISE_ENVIRONMENTAL    = 6,  /* EMI, vibration, thermal fluctuations */
} noise_mechanism_t;

/* ──── L1: Allan Variance Structure ──── */

/**
 * @brief Allan variance / Allan deviation for sensor noise characterization.
 *
 * Theorem (Allan, 1966):
 *   σ²(τ) = ½ ⟨(ȳ_{k+1} - ȳ_k)²⟩
 *   where ȳ_k is the average over interval τ starting at time kτ.
 *
 * Allan variance reveals noise sources by their τ-dependence:
 *   σ²(τ) ∝ τ⁻¹  →  white noise / angle random walk
 *   σ²(τ) ∝ τ⁰   →  flicker (1/f) noise / bias instability
 *   σ²(τ) ∝ τ¹   →  random walk / rate random walk
 *   σ²(τ) ∝ τ²   →  drift / ramp
 *   σ²(τ) ∝ τ³   →  quantization noise (for frequency data)
 */
typedef struct {
    double  *tau_values;       /* [s] cluster times (observation intervals) */
    double  *allan_var;        /* [input_unit²] Allan variance at each τ */
    double  *allan_dev;        /* [input_unit] Allan deviation = sqrt(allan_var) */
    size_t   num_clusters;     /* number of τ values computed */

    /* Extracted noise coefficients from the Allan variance log-log slope */
    double   angle_random_walk;    /* white noise coefficient [unit/√Hz] */
    double   bias_instability;     /* flicker floor minimum [unit] */
    double   rate_random_walk;     /* random walk coefficient [unit/√s] */
    double   quantization_noise;   /* quantization [unit] */
    double   drift_ramp;           /* rate ramp [unit/s] */

    double   min_allan_dev;        /* minimum Allan deviation (bias instability floor) */
    double   tau_at_min;           /* τ at which minimum occurs */

} allan_variance_t;

/* ──── L1: Core Sensor Function Declarations ──── */

int    sensor_specs_validate(const sensor_specs_t *s);
double sensor_total_error_percent(const sensor_specs_t *s);
double sensor_dynamic_error_percent(double f_input, double f_bandwidth);
double sensor_simulate_reading(double x_true, const sensor_specs_t *s,
                                double s_noise, int ascending);
int    sensor_uncertainty_budget_compute(sensor_uncertainty_budget_t *b);
int    sensor_identity_same_model(const sensor_identity_t *a,
                                   const sensor_identity_t *b);
double sensor_estimate_power_passive(double excitation, double resistance,
                                      int is_current_excitation);
double sensor_thermal_noise_limit(double temperature_K, double resistance,
                                   double bandwidth, double sensitivity);
int    sensor_reading_check_valid(const sensor_reading_t *r,
                                   double x_min, double x_max);
double sensor_measurement_rate(const sensor_reading_t *r0,
                                const sensor_reading_t *r1);
double sensor_noise_exponent_estimate(double S1, double f1,
                                       double S2, double f2);
double sensor_noise_equivalent_bandwidth_1st_order(double tau);
double sensor_noise_factor(double total_noise_Vrms, double gain,
                            double bandwidth, double source_resistance);

#endif /* SENSOR_DEFS_H */
