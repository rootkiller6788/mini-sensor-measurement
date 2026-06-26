/**
 * sigcond_linearize.h - Sensor Nonlinearity Correction
 *
 * Polynomial fitting, LUT interpolation, Steinhart-Hart (thermistor),
 * Callendar-Van Dusen (RTD), NIST ITS-90 (thermocouple).
 *
 * L2-L5: Horner evaluation, interpolation, inverse functions
 * Reference: Pallas-Areny & Webster (2001), NIST Monograph 175
 */
#ifndef SIGCOND_LINEARIZE_H
#define SIGCOND_LINEARIZE_H
#include "sigcond_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

/** Horner polynomial evaluation: p(x)=c[0]+x*(c[1]+x*(c[2]+...+x*c[n])). O(n) complexity. */
double poly_eval_horner(unsigned degree, const double c[], double x, double x_min, double x_max);

/** Simultaneous evaluation of p(x) and p'(x) for Newton-Raphson inversion */
void poly_eval_with_derivative(unsigned degree, const double c[], double x, double *p_val, double *p_deriv);

/** Polynomial inverse via Newton-Raphson: solve p(x)=y for x */
double poly_inverse_newton(unsigned degree, const double c[], double y_target, double x_guess, double tol, unsigned max_iter);

/** Linear interpolation from sorted LUT: y=y_k+(y_{k+1}-y_k)*(x-x_k)/(x_{k+1}-x_k) */
double lut_linear_interp(unsigned n_points, const double x_vals[], const double y_vals[], double x_query);

/** Quadratic Lagrange interpolation (3 nearest points) */
double lut_quadratic_interp(unsigned n_points, const double x_vals[], const double y_vals[], double x_query);

/** Cubic Catmull-Rom spline interpolation */
double lut_cubic_spline(unsigned n_points, const double x_vals[], const double y_vals[], double x_query);

/** Binary search for interval containing x. Returns i where x_vals[i] <= x < x_vals[i+1]. */
int binary_search_interval(unsigned n, const double x_vals[], double x, unsigned *idx);

/** Steinhart-Hart 3-param: 1/T = A + B*ln(R) + C*(ln(R))^3. Returns T in Kelvin. */
double steinhart_hart_T(const steinhart_hart_params_t *params, double resistance_ohm);

/** Steinhart-Hart inverse: compute R from T using cubic formula */
double steinhart_hart_R(const steinhart_hart_params_t *params, double temperature_kelvin);

/** Fit Steinhart-Hart coefficients from 3 (T,R) calibration points */
void steinhart_hart_fit_3pt(double t1_c, double r1_ohm, double t2_c, double r2_ohm, double t3_c, double r3_ohm, steinhart_hart_params_t *params);

/** NTC thermistor beta: B = ln(R1/R2) / (1/T1 - 1/T2) */
double thermistor_beta(double t1_kelvin, double r1_ohm, double t2_kelvin, double r2_ohm);

/** Initialize standard Pt100 CVD parameters per IEC 60751 */
void rtd_cvd_init_pt100(rtd_cvd_params_t *params);

/** Initialize standard Pt1000 CVD parameters */
void rtd_cvd_init_pt1000(rtd_cvd_params_t *params);

/** RTD resistance from temperature: R(t)=R0*(1+A*t+B*t^2+...), IEC 60751 */
double rtd_cvd_R_from_T(const rtd_cvd_params_t *params, double temperature_c);

/** Temperature from RTD resistance: analytical for t>=0, iterative for t<0 */
double rtd_cvd_T_from_R(const rtd_cvd_params_t *params, double resistance_ohm);

/** RTD alpha coefficient at reference */
double rtd_alpha(const rtd_cvd_params_t *params);

/** NIST ITS-90 forward: thermocouple T to V */
double nist_tc_temp_to_voltage(const thermocouple_nist_model_t *model, double temperature_c);

/** NIST ITS-90 inverse: thermocouple V to T */
double nist_tc_voltage_to_temp(const thermocouple_nist_model_t *model, double voltage_uv);

/** Initialize NIST models for Types K, J, T, E */
void nist_tc_init_typeK(thermocouple_nist_model_t *model);
void nist_tc_init_typeJ(thermocouple_nist_model_t *model);
void nist_tc_init_typeT(thermocouple_nist_model_t *model);
void nist_tc_init_typeE(thermocouple_nist_model_t *model);

/** Piecewise linearization with hysteresis model for sensors with different fwd/rev characteristics */
double piecewise_linear_hysteresis(unsigned n_segments, const double breakpoints[], const double fwd_slopes[], const double fwd_intercepts[], const double rev_slopes[], const double rev_intercepts[], double x, bool direction_increasing);

/** Build piecewise linear segments from calibration points */
void build_piecewise_linear(unsigned n_points, const calibration_point_t points[], double breakpoints[], double slopes[], double intercepts[], unsigned *n_segments);

/** Least-squares polynomial fit to calibration data (degree <= 5) */
void poly_fit_least_squares(unsigned n_points, const calibration_point_t points[], unsigned target_degree, polynomial_coeffs_t *result);

/** Maximum residual error of polynomial fit */
double poly_max_residual(unsigned n_points, const calibration_point_t points[], const polynomial_coeffs_t *poly);

#ifdef __cplusplus
}
#endif
#endif /* SIGCOND_LINEARIZE_H */
