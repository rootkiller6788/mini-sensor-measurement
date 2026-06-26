#ifndef MEASUREMENT_COMPENSATION_H
#define MEASUREMENT_COMPENSATION_H

#include <stddef.h>

/**
 * @file measurement_compensation.h
 * @brief L6: Canonical error compensation problems in measurement systems
 *
 * L6 Canonical Problems:
 *   - Temperature compensation (thermocouple cold-junction, strain gauge TCS)
 *   - Nonlinearity correction (polynomial inverse, piecewise linear LUT)
 *   - Zero/span drift compensation
 *   - Bridge nonlinearity correction (Wheatstone bridge)
 *   - Lead wire resistance compensation (3-wire and 4-wire RTD)
 *
 * Courses: MIT 6.630 (sensors), Berkeley EE117 (instrumentation),
 *          Michigan EECS 411 (microwave measurements)
 */

/* --- L6: Temperature Compensation --- */

typedef struct {
    double ref_temperature;     /* reference temperature (typically 25 C) */
    double tc_zero;             /* zero drift coefficient (units / C) */
    double tc_span;             /* span/gain drift coefficient (1 / C) */
    double tc_quadratic;        /* quadratic term (1 / C^2, optional) */
    int    has_quadratic;       /* nonzero if quadratic term is used */
} TempCompensationModel;

/**
 * Apply temperature compensation to a measurement.
 *
 * x_corrected = (x_raw - tc_zero * dT) / (1 + tc_span * dT + tc_quad * dT^2)
 * where dT = T_actual - T_ref
 *
 * @param x_raw       uncorrected measurement
 * @param T_actual    actual temperature (deg C)
 * @param model       temperature compensation parameters
 * @return temperature-compensated measurement
 */
double tempcomp_apply(double x_raw, double T_actual,
                       const TempCompensationModel *model);

/**
 * Thermocouple cold-junction compensation.
 * V_actual = V_measured + S * (T_cj - T_ref)
 * where S is the Seebeck coefficient of the thermocouple at T_cj.
 *
 * @param V_measured       measured thermocouple voltage (mV)
 * @param T_cold_junction  cold junction temperature (C)
 * @param T_reference      reference temperature (typically 0 C)
 * @param seebeck_coeff    Seebeck coefficient (mV/C)
 * @return compensated voltage
 */
double tempcomp_thermocouple_cjc(double V_measured, double T_cold_junction,
                                  double T_reference, double seebeck_coeff);

/**
 * Strain gauge temperature compensation.
 * Apparent strain = (alpha_gauge - alpha_substrate) * dT / GF
 * where GF is gauge factor.
 *
 * @param measured_strain    raw strain reading (microstrain)
 * @param alpha_gauge        CTE of gauge material (ppm/C)
 * @param alpha_substrate    CTE of substrate material (ppm/C)
 * @param dT                 temperature change from reference
 * @param gauge_factor       gauge factor (typically 2.0)
 * @return compensated strain
 */
double tempcomp_strain_gauge(double measured_strain,
                              double alpha_gauge, double alpha_substrate,
                              double dT, double gauge_factor);

/* --- L6: Nonlinearity Correction --- */

/**
 * Polynomial inverse correction using pre-computed coefficients.
 * If sensor has y = f(x) but we need x = g(y), precompute g coefficients.
 *
 * x_corrected = c0 + c1*y + c2*y^2 + c3*y^3
 *
 * @param y_measured  sensor output
 * @param coeffs      correction polynomial coefficients
 * @param order       polynomial order (1-5)
 * @return corrected value
 */
double nlin_polynomial_correction(double y_measured, const double *coeffs,
                                   size_t order);

/**
 * Piecewise linear LUT correction for arbitrary nonlinearity.
 *
 * For y in [y_i, y_{i+1}]:
 *   x = x_i + (x_{i+1} - x_i) * (y - y_i) / (y_{i+1} - y_i)
 *
 * @param y_measured    sensor output
 * @param y_breakpoints LUT input breakpoints (sensor output)
 * @param x_breakpoints LUT output breakpoints (corrected value)
 * @param n_points      number of breakpoints
 * @return corrected value
 */
double nlin_lut_correction(double y_measured,
                            const double *y_breakpoints,
                            const double *x_breakpoints,
                            size_t n_points);

/* --- L6: Wheatstone Bridge Nonlinearity Correction --- */

/**
 * Correct Wheatstone bridge nonlinearity for quarter-bridge configuration.
 *
 * Exact relationship for quarter bridge:
 *   V_out / V_ex = (GF * strain) / (4 + 2 * GF * strain)
 *
 * Nonlinearity error = GF * strain / 2  (fractional)
 *
 * @param measured_strain  apparent strain from linear approximation
 * @param gauge_factor     gauge factor
 * @return corrected strain
 */
double bridge_quarter_linearity_correction(double measured_strain,
                                            double gauge_factor);

/**
 * Half-bridge correction (push-pull, opposite arms).
 * Nonlinearity is reduced compared to quarter bridge but not zero.
 */
double bridge_half_linearity_correction(double measured_strain,
                                         double gauge_factor);

/* --- L6: Lead Wire Resistance Compensation --- */

/**
 * 3-wire RTD lead resistance compensation.
 * R_rtd = R_measured - R_lead
 * where R_lead is measured via the third wire.
 *
 * @param R_measured    total resistance including lead wires
 * @param R_lead        lead wire resistance per wire
 * @param num_leads     number of lead wires in measurement path (2 or 3)
 * @return compensated RTD resistance
 */
double leadwire_3wire_compensation(double R_measured, double R_lead,
                                    int num_leads);

/**
 * 4-wire (Kelvin) measurement: ideal compensation.
 * Forces current through outer leads, senses voltage on inner leads.
 * R_rtd = V_sense / I_force (no lead resistance error in theory).
 *
 * @param V_sense      voltage sensed across RTD
 * @param I_force      forced current
 * @param contact_R    contact resistance per lead (for error estimate)
 * @param error_est    [out] estimated residual error
 * @return RTD resistance
 */
double kelvin_4wire_measurement(double V_sense, double I_force,
                                 double contact_R, double *error_est);

/* --- L6: Combined Compensation Pipeline --- */

/**
 * Apply a full compensation pipeline to a raw measurement:
 * 1. Zero offset removal
 * 2. Nonlinearity correction
 * 3. Temperature compensation
 * 4. Span/gain correction
 *
 * @param raw             raw sensor output
 * @param zero_offset     additive offset
 * @param nlin_coeffs     nonlinearity correction (NULL = skip)
 * @param nlin_order      polynomial order for nlin
 * @param tempcomp        temperature compensation model (NULL = skip)
 * @param T_actual        actual temperature
 * @param span_gain       multiplicative gain correction
 * @return fully corrected measurement
 */
double compensation_pipeline(double raw, double zero_offset,
                              const double *nlin_coeffs, size_t nlin_order,
                              const TempCompensationModel *tempcomp,
                              double T_actual, double span_gain);

#endif /* MEASUREMENT_COMPENSATION_H */
