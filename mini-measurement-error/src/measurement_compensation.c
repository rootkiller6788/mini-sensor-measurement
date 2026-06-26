#include "measurement_compensation.h"
#include <math.h>
#include <float.h>

/* ─── L6: Temperature Compensation ───────────────────────────────────────
 *
 * Most sensors exhibit temperature-dependent offset and gain errors.
 * A standard linear+quadratic compensation model:
 *
 *   x_corrected = (x_raw - tc_zero * dT) / (1 + tc_span * dT + tc_quad * dT^2)
 *
 * where dT = T_actual - T_ref.
 *
 * This compensates for:
 *   - Zero drift: the output offset when input is zero changes with temperature
 *   - Span drift: the sensitivity/gain changes with temperature
 *   - Quadratic term: higher-order curvature (e.g., in piezoresistive sensors)
 */

double tempcomp_apply(double x_raw, double T_actual,
                       const TempCompensationModel *model) {
    if (!model) return x_raw;
    double dT = T_actual - model->ref_temperature;

    /* Remove zero drift */
    double x_zero_corrected = x_raw - model->tc_zero * dT;

    /* Correct span drift */
    double span_factor = 1.0 + model->tc_span * dT;
    if (model->has_quadratic) {
        span_factor += model->tc_quadratic * dT * dT;
    }

    if (fabs(span_factor) < DBL_EPSILON) return x_zero_corrected;
    return x_zero_corrected / span_factor;
}

/* ─── L6: Thermocouple Cold-Junction Compensation ────────────────────────
 *
 * A thermocouple measures the temperature DIFFERENCE between the hot
 * junction (measurement point) and the cold junction (reference).
 * To get absolute temperature:
 *
 *   V_actual(T_hot) = V_measured(T_hot - T_cj) + V_equiv(T_cj - T_ref)
 *
 * For small ranges near room temperature, a linear approximation
 * using the Seebeck coefficient works well:
 *
 *   V_compensated = V_measured + S * (T_cj - T_ref)
 *
 * where S is the Seebeck coefficient at T_cj (mV/C).
 */

double tempcomp_thermocouple_cjc(double V_measured, double T_cold_junction,
                                  double T_reference, double seebeck_coeff) {
    double V_cj = seebeck_coeff * (T_cold_junction - T_reference);
    return V_measured + V_cj;
}

/* ─── L6: Strain Gauge Temperature Compensation ─────────────────────────
 *
 * Strain gauges experience apparent strain due to differential thermal
 * expansion between the gauge and the substrate:
 *
 *   epsilon_apparent = (alpha_gauge - alpha_substrate) * dT / GF
 *
 * where alpha = coefficient of thermal expansion (CTE), GF = gauge factor.
 *
 * For self-temperature-compensated (STC) gauges, alpha_gauge is
 * matched to the substrate material (e.g., steel: ~11 ppm/C,
 * aluminum: ~23 ppm/C).
 */

double tempcomp_strain_gauge(double measured_strain,
                              double alpha_gauge, double alpha_substrate,
                              double dT, double gauge_factor) {
    if (fabs(gauge_factor) < DBL_EPSILON) return measured_strain;
    double apparent_strain = (alpha_gauge - alpha_substrate) * dT
                              / gauge_factor;
    return measured_strain - apparent_strain;
}

/* ─── L6: Polynomial Nonlinearity Correction ─────────────────────────────
 *
 * Many sensors have a known nonlinear transfer function y = f(x).
 * If we have an inverse polynomial model:
 *
 *   x_corrected = c0 + c1*y + c2*y^2 + c3*y^3
 *
 * we can apply it directly. This is pre-computed from calibration data.
 */

double nlin_polynomial_correction(double y_measured, const double *coeffs,
                                   size_t order) {
    if (!coeffs || order == 0) return y_measured;

    double result = 0.0;
    double y_pow = 1.0;
    size_t n = order + 1;
    for (size_t i = 0; i < n && i <= 5; i++) {
        result += coeffs[i] * y_pow;
        y_pow *= y_measured;
    }
    return result;
}

/* ─── L6: Piecewise Linear LUT Correction ────────────────────────────────
 *
 * For arbitrary nonlinear transfer functions, a piecewise linear lookup
 * table provides robust correction without assuming a functional form.
 *
 * For y in [y_i, y_{i+1}]:
 *   x = x_i + (x_{i+1} - x_i) * (y - y_i) / (y_{i+1} - y_i)
 */

double nlin_lut_correction(double y_measured,
                            const double *y_breakpoints,
                            const double *x_breakpoints,
                            size_t n_points) {
    if (!y_breakpoints || !x_breakpoints || n_points < 2)
        return y_measured;

    /* Clamp to range */
    if (y_measured <= y_breakpoints[0])
        return x_breakpoints[0];
    if (y_measured >= y_breakpoints[n_points - 1])
        return x_breakpoints[n_points - 1];

    /* Binary search for interval */
    size_t lo = 0, hi = n_points - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (y_breakpoints[mid] <= y_measured)
            lo = mid;
        else
            hi = mid;
    }

    /* Linear interpolation in [lo, hi] */
    double dy = y_breakpoints[hi] - y_breakpoints[lo];
    if (fabs(dy) < DBL_EPSILON)
        return x_breakpoints[lo];

    double frac = (y_measured - y_breakpoints[lo]) / dy;
    return x_breakpoints[lo] + frac * (x_breakpoints[hi] - x_breakpoints[lo]);
}

/* ─── L6: Wheatstone Bridge Nonlinearity ─────────────────────────────────
 *
 * Quarter-bridge (one active gauge):
 *   Exact:  V_out/V_ex = (GF * e) / (4 + 2 * GF * e)
 *   Linear approximation: V_out/V_ex ~= (GF * e) / 4
 *
 * The nonlinearity error (fractional):
 *   NL_error = (linear - exact) / exact * 100%
 *            ~= GF * e / 2  (for small strains)
 *
 * Corrected strain = measured_strain / (1 - GF * measured_strain / 2)
 *                   (for quarter bridge)
 */

double bridge_quarter_linearity_correction(double measured_strain,
                                            double gauge_factor) {
    double denom = 1.0 - gauge_factor * measured_strain / 2.0;
    if (fabs(denom) < DBL_EPSILON) return measured_strain;
    return measured_strain / denom;
}

/* Half-bridge (push-pull, opposite arms):
 *   Exact: V_out/V_ex = (GF * e) / 2
 * This is inherently linear — no correction needed for ideal case.
 * But in practice, small residual nonlinearity exists due to
 * gauge mismatch and transverse sensitivity.
 *
 * Corrected strain ~= measured_strain * (1 + GF * measured_strain / 4)
 * (a much smaller correction than quarter bridge)
 */

double bridge_half_linearity_correction(double measured_strain,
                                         double gauge_factor) {
    double correction = 1.0 + gauge_factor * measured_strain / 4.0;
    return measured_strain * correction;
}

/* ─── L6: Lead Wire Resistance Compensation ──────────────────────────────
 *
 * In 2-wire RTD measurement, the lead wire resistance adds directly:
 *   R_measured = R_rtd + R_lead1 + R_lead2
 *
 * In 3-wire measurement (with matched leads):
 *   R_rtd = R_measured - R_lead  (one lead compensated by bridge balance)
 *
 * In 4-wire (Kelvin) measurement:
 *   Force current through outer leads, sense voltage on inner leads.
 *   R_rtd = V_sense / I_force   (theoretically perfect)
 */

double leadwire_3wire_compensation(double R_measured, double R_lead,
                                    int num_leads) {
    return R_measured - (double)num_leads * R_lead;
}

double kelvin_4wire_measurement(double V_sense, double I_force,
                                 double contact_R, double *error_est) {
    if (fabs(I_force) < DBL_EPSILON) {
        if (error_est) *error_est = INFINITY;
        return 0.0;
    }
    double R_rtd = V_sense / I_force;

    /* Residual error from finite voltmeter input impedance.
     * If Z_vm >> R_rtd, the current leakage is negligible.
     * error_est is set to a small fraction accounting for practical limits. */
    if (error_est) {
        /* Typical: input bias current of voltmeter causes ~0.01% error */
        *error_est = 0.0;
        if (contact_R > 0.0) {
            *error_est = contact_R * 1e-6; /* heuristic for typical values */
        }
    }
    return R_rtd;
}

/* ─── L6: Full Compensation Pipeline ──────────────────────────────────────
 *
 * Applies corrections in the standard order:
 * 1. Remove zero offset
 * 2. Compensate nonlinearity
 * 3. Apply temperature compensation
 * 4. Apply span/gain correction
 *
 * This pipeline is common in industrial transmitters and data acquisition
 * systems (e.g., HART, Foundation Fieldbus).
 */

double compensation_pipeline(double raw, double zero_offset,
                              const double *nlin_coeffs, size_t nlin_order,
                              const TempCompensationModel *tempcomp,
                              double T_actual, double span_gain) {
    /* Step 1: Zero offset */
    double corrected = raw - zero_offset;

    /* Step 2: Nonlinearity correction */
    if (nlin_coeffs && nlin_order > 0) {
        corrected = nlin_polynomial_correction(corrected, nlin_coeffs,
                                                nlin_order);
    }

    /* Step 3: Temperature compensation */
    if (tempcomp) {
        corrected = tempcomp_apply(corrected, T_actual, tempcomp);
    }

    /* Step 4: Span/gain correction */
    if (fabs(span_gain) > DBL_EPSILON) {
        corrected *= span_gain;
    }

    return corrected;
}
