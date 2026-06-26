/**
 * @file sensor_transfer.h
 * @brief Sensor transfer functions, sensitivity, linearity, and dynamic response
 *
 * Level: L2 — Core Concepts · L3 — Mathematical Structures · L4 — Fundamental Laws
 *
 * Every sensor is characterized by its transfer function H(s):
 *   Y(s) = H(s) * X(s) + N(s)
 *   where X is the measurand, Y is the output, N is additive noise.
 *
 * Static:  y = f(x)  (polynomial, exponential, rational)
 * Dynamic: H(s) = K / (τ s + 1)              [1st order]
 *          H(s) = K ω_n² / (s² + 2ζ ω_n s + ω_n²)  [2nd order]
 *
 * Key Concepts:
 *   - Sensitivity: S = dy/dx at operating point
 *   - Linearity error: max |f(x) - L(x)|/FS where L is best-fit line
 *   - Hysteresis: energy dissipation in sensing element
 *   - Bandwidth: -3 dB frequency from Bode plot
 *   - Resonance: peak in 2nd order systems near ω_n when ζ < 0.707
 *
 * References:
 *   - Doebelin, "Measurement Systems", 6th Ed, Ch. 3-4
 *   - Fraden, "Handbook of Modern Sensors", 5th Ed, Ch. 2
 *   - Ogata, "System Dynamics", 4th Ed, Ch. 5
 */

#ifndef SENSOR_TRANSFER_H
#define SENSOR_TRANSFER_H

#include "sensor_defs.h"
#include <stddef.h>

/* ──── L2: Polynomial Transfer Function ──── */

/**
 * @brief Polynomial sensor model (most common static model).
 *
 * General form: y = a₀ + a₁·x + a₂·x² + ... + aₙ·xⁿ
 *
 * Theorem (Weierstrass Approximation):
 *   Any continuous function f(x) on [a,b] can be uniformly approximated
 *   by a polynomial to any desired accuracy.
 *   This justifies polynomial calibration curves for nonlinear sensors.
 *
 * Common special forms:
 *   Order 0: y = a₀          → offset-only (ideal switch)
 *   Order 1: y = a₀ + a₁·x   → linear (ideal sensor)
 *   Order 2: y = a₀ + a₁·x + a₂·x²  → quadratic (thermocouple EMF)
 *   Order 3: y = a₀ + a₁·x + a₂·x² + a₃·x³ → cubic (RTD, thermistor approx.)
 *   Order 4: y = a₀ + a₁·x + ... + a₄·x⁴  → Callendar-Van Dusen (PRT)
 */
typedef struct {
    double  *coeff;   /* [a₀, a₁, ..., aₙ], length = order+1 */
    uint8_t  order;   /* polynomial degree (≥ 0) */
    double   x_min;   /* valid input range minimum */
    double   x_max;   /* valid input range maximum */
} poly_transfer_t;

int  poly_transfer_init(poly_transfer_t *pt, const double *coeff, uint8_t order,
                        double x_min, double x_max);
void poly_transfer_free(poly_transfer_t *pt);
double poly_transfer_evaluate(const poly_transfer_t *pt, double x);
double poly_transfer_inverse(const poly_transfer_t *pt, double y, double tol);
double poly_transfer_sensitivity(const poly_transfer_t *pt, double x);
double poly_transfer_nonlinearity(const poly_transfer_t *pt, size_t n_samples);

/* ──── L2: Thermistor Transfer Function (Steinhart-Hart) ──── */

/**
 * @brief Steinhart-Hart equation for NTC/PTC thermistors.
 *
 * Theorem (Steinhart & Hart, 1968):
 *   1/T = A + B·ln(R) + C·[ln(R)]³
 *   where T [K], R [Ω], and A, B, C are thermistor-specific constants.
 *
 * Inverse (R from T) requires solving the cubic in ln(R):
 *   C·[ln(R)]³ + 0·[ln(R)]² + B·[ln(R)] + (A - 1/T) = 0
 *
 * Simplified β-parameter model (2-coefficient):
 *   R(T) = R₀ · exp[β · (1/T - 1/T₀)]
 *   valid over ~50°C range with ~0.1°C accuracy
 *
 * References:
 *   - Steinhart & Hart, "Calibration curves for thermistors",
 *     Deep-Sea Research, 15(4):497-503, 1968
 */
typedef struct {
    double A;  /* [1/K] Steinhart-Hart coefficient */
    double B;  /* [1/K] Steinhart-Hart coefficient */
    double C;  /* [1/K] Steinhart-Hart coefficient */
    double R_at_T0; /* [Ω] reference resistance for β-model */
    double T0;      /* [K] reference temperature for β-model */
    double B_param; /* [K] β parameter for simplified model */
    uint8_t use_beta; /* 1 = use β-model, 0 = use full S-H model */
} thermistor_model_t;

int    thermistor_model_init_full(thermistor_model_t *tm,
                                  double A, double B, double C);
int    thermistor_model_init_beta(thermistor_model_t *tm,
                                  double R0, double T0_K, double beta);
double thermistor_temp_from_R(const thermistor_model_t *tm, double R);
double thermistor_R_from_temp(const thermistor_model_t *tm, double T_K);
double thermistor_sensitivity(const thermistor_model_t *tm, double T_K);
double thermistor_self_heating(double current, double R_th,
                                double dissipation_constant_mW_K);
double thermistor_linearization_Rp(const thermistor_model_t *tm,
                                    double T_center, double *Rp_out);

/* ──── L2: RTD Transfer Function (Callendar-Van Dusen) ──── */

/**
 * @brief Callendar-Van Dusen equation for platinum RTDs (IEC 60751).
 *
 * For T ≥ 0°C:
 *   R(T) = R₀ · [1 + A·T + B·T²]
 *
 * For T < 0°C (full Callendar-Van Dusen):
 *   R(T) = R₀ · [1 + A·T + B·T² + C·(T-100)·T³]
 *
 * IEC 60751 standard coefficients for α=0.00385055 (Pt100):
 *   A = 3.9083×10⁻³ °C⁻¹
 *   B = -5.775×10⁻⁷ °C⁻²
 *   C = -4.183×10⁻¹² °C⁻⁴
 *
 * R₀ = 100 Ω (Pt100), 500 Ω (Pt500), 1000 Ω (Pt1000)
 *
 * References:
 *   - IEC 60751:2008 "Industrial platinum resistance thermometers"
 *   - Callendar, "On the practical measurement of temperature",
 *     Phil. Trans. Roy. Soc. A, 178:161-230, 1887
 */
typedef struct {
    double R0;      /* [Ω] resistance at 0°C */
    double A;       /* [°C⁻¹] linear coefficient */
    double B;       /* [°C⁻²] quadratic coefficient */
    double C;       /* [°C⁻⁴] sub-zero correction coefficient */
    double alpha;   /* temperature coefficient: α = (R₁₀₀-R₀)/(100·R₀) */
} rtd_cvd_model_t;

int    rtd_cvd_init_iec751(rtd_cvd_model_t *rtd, double R0);
int    rtd_cvd_init_custom(rtd_cvd_model_t *rtd,
                            double R0, double alpha, double delta, double beta);
double rtd_R_from_T(const rtd_cvd_model_t *rtd, double T_C);
double rtd_T_from_R(const rtd_cvd_model_t *rtd, double R, double tol);
double rtd_sensitivity(const rtd_cvd_model_t *rtd, double T_C);
double rtd_self_heating(double I_meas, double R_element,
                         double thermal_resistance_K_W);
double rtd_lead_resistance_2wire(double R_lead, double R_measured, double *R_sensor);

/* ──── L2: Thermocouple Transfer Function ──── */

/**
 * @brief Thermocouple EMF (Seebeck effect) modeled by ITS-90 polynomials.
 *
 * Theorem (Seebeck, 1821):
 *   The open-circuit voltage between two dissimilar conductors
 *   with junctions at different temperatures is:
 *   V = ∫_{T_ref}^{T_meas} [S_A(T) - S_B(T)] dT
 *   where S(T) is the absolute Seebeck coefficient of each material.
 *
 * NIST ITS-90 polynomial (Type K, as example):
 *   E(μV) = Σ c_i · t^i + a₀·exp[a₁·(t - a₂)²]
 *   where t = T[°C], and coefficients depend on temperature range.
 *
 * Types defined: B, E, J, K, N, R, S, T (IEC 60584-1)
 *
 * Cold Junction Compensation (CJC):
 *   V_actual(T_meas, T_ref) = V_table(T_meas, 0°C) - V_table(T_ref, 0°C)
 *   where T_ref is measured at the cold junction (typically via thermistor/RTD).
 *
 * Law of Intermediate Metals:
 *   A third metal inserted in series does not affect the net EMF
 *   if both ends of the third metal are at the same temperature.
 *
 * References:
 *   - NIST Monograph 175, "Temperature-Electromotive Force Reference Functions"
 *   - IEC 60584-1:2013 "Thermocouples — Part 1: EMF specifications and tolerances"
 *   - Pollock, "Thermocouples: Theory and Properties", CRC Press, 1991
 */
typedef enum {
    TC_TYPE_B = 0,  /* Pt-30%Rh / Pt-6%Rh, 0 to 1820°C */
    TC_TYPE_E = 1,  /* Ni-Cr / Cu-Ni, -270 to 1000°C, highest sensitivity */
    TC_TYPE_J = 2,  /* Fe / Cu-Ni, -210 to 1200°C */
    TC_TYPE_K = 3,  /* Ni-Cr / Ni-Al, -270 to 1372°C, most common */
    TC_TYPE_N = 4,  /* Ni-Cr-Si / Ni-Si-Mg, -270 to 1300°C, improved stability */
    TC_TYPE_R = 5,  /* Pt-13%Rh / Pt, -50 to 1768°C */
    TC_TYPE_S = 6,  /* Pt-10%Rh / Pt, -50 to 1768°C, defining ITS-90 */
    TC_TYPE_T = 7,  /* Cu / Cu-Ni, -270 to 400°C, best for sub-zero */
} tc_type_t;

typedef struct {
    tc_type_t type;
    /* NIST ITS-90 polynomial coefficients (per sub-range) */
    double   c[16];          /* up to 15th order polynomial + exponential */
    double   a0_exp, a1_exp, a2_exp;  /* exponential correction coefficients */
    uint8_t  range_index;    /* which temperature sub-range is active */
    double   t_min_C;        /* valid temperature range minimum */
    double   t_max_C;        /* valid temperature range maximum */
    double   sensitivity_typical;  /* [µV/°C] at typical operating temp */
    double   tolerance_class1;     /* [°C] Class 1 tolerance at ref temp */
    double   tolerance_class2;     /* [°C] Class 2 tolerance at ref temp */
} tc_model_t;

int    tc_model_init(tc_model_t *tc, tc_type_t type);
double tc_emf_from_T(const tc_model_t *tc, double T_meas_C, double T_ref_C);
double tc_T_from_emf(const tc_model_t *tc, double emf_uV, double T_ref_C, double tol);
double tc_sensitivity(const tc_model_t *tc, double T_C);
double tc_cjc_compensate(const tc_model_t *tc, double V_measured,
                          double T_cj_measured, double T_cj_ref);
int    tc_law_intermediate_metals_check(double V_emf, double T1, double T2, double T3);

/* ──── L2: 1st-Order Dynamic Sensor Model ──── */

/**
 * @brief First-order dynamic system model for sensors.
 *
 * Transfer function: H(s) = K / (τ·s + 1)
 *
 * Step response: y(t) = K·A·[1 - exp(-t/τ)]
 *   - Time constant τ: time to reach 63.2% of final value
 *   - Rise time (10%-90%): t_r = τ·ln(9) ≈ 2.197τ
 *   - Settling time (1%): t_s = τ·ln(100) ≈ 4.605τ
 *
 * Frequency response: H(jω) = K / (1 + jωτ)
 *   - Magnitude: |H| = K / √(1 + ω²τ²)
 *   - Phase: φ = -arctan(ωτ)
 *   - Bandwidth: f_{-3dB} = 1/(2πτ)
 *
 * Common examples: thermocouple (thermal mass), pressure sensor (fluid fill),
 *                  gas sensor (diffusion), humidity sensor (sorption)
 *
 * Theorem (1st-order ODE solution):
 *   τ·dy/dt + y = K·x(t)  has complementary solution y_c = y₀·exp(-t/τ)
 *   and particular solution depends on x(t).
 */
typedef struct {
    double K;    /* static gain [output/input] = DC sensitivity */
    double tau;  /* time constant [s] */
    double y0;   /* initial condition */
    double t0;   /* initial time */
} sensor_1st_order_t;

int    sensor_1st_order_init(sensor_1st_order_t *s1, double K, double tau,
                              double y0, double t0);
double sensor_1st_order_step_response(const sensor_1st_order_t *s1,
                                       double A_step, double t);
double sensor_1st_order_freq_response_mag(const sensor_1st_order_t *s1,
                                            double f_hz);
double sensor_1st_order_freq_response_phase(const sensor_1st_order_t *s1,
                                              double f_hz);
double sensor_1st_order_bandwidth(const sensor_1st_order_t *s1);
double sensor_1st_order_rise_time(const sensor_1st_order_t *s1);
double sensor_1st_order_settling_time(const sensor_1st_order_t *s1,
                                       double percent);
int    sensor_1st_order_simulate(const sensor_1st_order_t *s1,
                                  const double *input, double *output,
                                  size_t n, double dt);

/* ──── L2: 2nd-Order Dynamic Sensor Model ──── */

/**
 * @brief Second-order dynamic system model for sensors.
 *
 * Transfer function: H(s) = K·ω_n² / (s² + 2ζ·ω_n·s + ω_n²)
 *
 * Parameters:
 *   ω_n: natural frequency [rad/s]
 *   ζ: damping ratio (dimensionless)
 *   K: static gain
 *
 * Damping regimes:
 *   ζ = 0: undamped, sustained oscillation at ω_n
 *   0 < ζ < 1: underdamped, oscillatory with overshoot
 *   ζ = 1: critically damped, fastest non-overshooting response
 *   ζ > 1: overdamped, sluggish, no overshoot
 *
 * Step response (ζ < 1):
 *   y(t) = K·A·[1 - (1/√(1-ζ²))·exp(-ζ·ω_n·t)·sin(ω_d·t + φ)]
 *   where ω_d = ω_n·√(1-ζ²), φ = arctan(√(1-ζ²)/ζ)
 *
 * Peak time: t_p = π / ω_d
 * Overshoot: M_p = exp(-π·ζ/√(1-ζ²)) · 100%
 * Settling time (2%): t_s ≈ 4/(ζ·ω_n)
 *
 * Examples: accelerometer (mass-spring-damper), pressure transducer (diaphragm),
 *           geophone, seismometer, MEMS gyroscope
 *
 * References:
 *   - Doebelin, "Measurement Systems", 6th Ed, Ch. 4
 *   - Thomson, "Theory of Vibration with Applications", 5th Ed
 */
typedef struct {
    double K;       /* static gain [output/input] */
    double omega_n; /* natural frequency [rad/s] */
    double zeta;    /* damping ratio [dimensionless] */
    double y0;      /* initial output */
    double dy0;     /* initial derivative of output */
} sensor_2nd_order_t;

int    sensor_2nd_order_init(sensor_2nd_order_t *s2, double K,
                              double omega_n, double zeta,
                              double y0, double dy0);
double sensor_2nd_order_step_response(const sensor_2nd_order_t *s2,
                                       double A_step, double t);
double sensor_2nd_order_peak_time(const sensor_2nd_order_t *s2);
double sensor_2nd_order_overshoot_percent(const sensor_2nd_order_t *s2);
double sensor_2nd_order_settling_time(const sensor_2nd_order_t *s2,
                                       double band_percent);
double sensor_2nd_order_freq_response_mag(const sensor_2nd_order_t *s2,
                                            double f_hz);
double sensor_2nd_order_freq_response_phase(const sensor_2nd_order_t *s2,
                                              double f_hz);
double sensor_2nd_order_resonant_frequency(const sensor_2nd_order_t *s2);
double sensor_2nd_order_resonant_magnification(const sensor_2nd_order_t *s2);
double sensor_2nd_order_bandwidth(const sensor_2nd_order_t *s2);
int    sensor_2nd_order_simulate(const sensor_2nd_order_t *s2,
                                  const double *input, double *output,
                                  double *doutput, size_t n, double dt);

/* ──── L2: Nonlinear Transfer Function via Lookup Table ──── */

/**
 * @brief Piecewise linear interpolation lookup table for arbitrary nonlinear sensors.
 *
 * Given n calibration points (x_i, y_i) in ascending x order,
 * evaluates y = f(x) by linear interpolation between nearest points.
 *
 * Error bound for linear interpolation:
 *   |f(x) - L(x)| ≤ (1/8) · max|f''(ξ)| · (Δx)²
 *   where Δx is the maximum interval width.
 */
typedef struct {
    double *x_table;    /* input (measurand) values, length = n_points */
    double *y_table;    /* output (sensor reading) values, length = n_points */
    size_t  n_points;
    double  x_min, x_max;
} sensor_lut_t;

int    sensor_lut_init(sensor_lut_t *lut, const double *x, const double *y, size_t n);
void   sensor_lut_free(sensor_lut_t *lut);
double sensor_lut_evaluate(const sensor_lut_t *lut, double x);
double sensor_lut_inverse(const sensor_lut_t *lut, double y);
double sensor_lut_sensitivity(const sensor_lut_t *lut, double x);
int    sensor_lut_get_interval(const sensor_lut_t *lut, double x,
                                size_t *idx_low, size_t *idx_high);

/* ──── L3: Sensor Linearity Metrics ──── */

/**
 * @brief Compute independent linearity (best-fit straight line method).
 *
 * Finds coefficients m (slope) and b (intercept) minimizing:
 *   Σ (y_i - (m·x_i + b))²
 *
 * Closed-form solution (normal equations):
 *   m = [n·Σ(x_i·y_i) - Σx_i·Σy_i] / [n·Σ(x_i²) - (Σx_i)²]
 *   b = [Σy_i - m·Σx_i] / n
 *
 * Independent linearity error = max|y_i - (m·x_i + b)| / (y_max - y_min) · 100%
 */
int    sensor_linearity_independent(const double *x, const double *y,
                                     size_t n, double *slope, double *intercept,
                                     double *max_error_percent);

/**
 * @brief Compute endpoint linearity (straight line through endpoints).
 *
 * Endpoint line connects (x_0, y_0) to (x_{n-1}, y_{n-1}).
 * Simpler than independent linearity but typically yields larger error.
 */
int    sensor_linearity_endpoint(const double *x, const double *y,
                                  size_t n, double *max_error_percent);

/**
 * @brief Compute terminal-based linearity (zero-based endpoint).
 *
 * Line passes through origin: y = m·x where m is best-fit slope through origin.
 * Used when zero input must produce zero output by physical law.
 */
int    sensor_linearity_terminal(const double *x, const double *y,
                                  size_t n, double *slope,
                                  double *max_error_percent);

/* ──── L3: Hysteresis Computation ──── */

/**
 * @brief Compute hysteresis error from ascending/descending calibration data.
 *
 * Hysteresis = max|y_up(x_i) - y_down(x_i)| / FS · 100%
 *
 * Physical origins:
 *   - Magnetic domain wall pinning (magnetic sensors)
 *   - Elastic after-effect (strain gauges)
 *   - Molecular friction (piezoelectric)
 *   - Capillary condensation (humidity sensors)
 *   - Charge trapping (semiconductor sensors)
 */
double sensor_hysteresis_percent(const double *x, const double *y_up,
                                  const double *y_down, size_t n,
                                  double full_scale);

/* ──── L3: Dynamic Error Estimation ──── */

/**
 * @brief Estimate dynamic error for a 1st-order sensor tracking a ramp input.
 *
 * For ramp input x(t) = r·t:
 *   Steady-state error: e_ss = r·τ
 *   (Sensor lags behind by τ seconds, producing a constant offset)
 *
 * For sinusoidal input x(t) = A·sin(ωt):
 *   Magnitude error: ε_mag = |K/√(1+ω²τ²) - K|/K · 100%
 *   Phase lag: no amplitude error for DC-referenced, but phase causes
 *              timing error Δt = φ/ω
 */
double sensor_dynamic_ramp_error_1st_order(double ramp_rate, double tau);
double sensor_dynamic_sinusoidal_mag_error_1st_order(double freq_hz, double tau);
double sensor_dynamic_phase_lag_1st_order(double freq_hz, double tau);

/**
 * @brief Estimate dynamic error for a 2nd-order sensor tracking a ramp input.
 *
 * For ramp input x(t) = r·t:
 *   Steady-state error: e_ss = 2·ζ·r / ω_n
 */
double sensor_dynamic_ramp_error_2nd_order(double ramp_rate,
                                            double omega_n, double zeta);

#endif /* SENSOR_TRANSFER_H */
