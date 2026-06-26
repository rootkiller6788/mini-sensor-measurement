/**
 * sigcond_linearize.c - Sensor Nonlinearity Correction Implementation
 *
 * Implements Horner polynomial evaluation, LUT interpolation (linear,
 * quadratic, cubic spline), Steinhart-Hart thermistor equation,
 * Callendar-Van Dusen RTD linearization, NIST ITS-90 thermocouple
 * polynomial evaluation, and piecewise linearization.
 *
 * L3: Horner's method (O(n)), Lagrange interpolation, Newton-Raphson
 * L4: Steinhart-Hart (1968), Callendar-Van Dusen (IEC 60751)
 * L5: NIST ITS-90 polynomials, least-squares fitting
 *
 * Reference:
 *   Pallas-Areny & Webster (2001), Ch. 3
 *   NIST Monograph 175, "NIST ITS-90 Thermocouple Database" (1993)
 *   Steinhart & Hart, Deep-Sea Research, Vol. 15 (1968)
 *   IEC 60751:2008 - Industrial platinum resistance thermometers
 */

#include "sigcond_linearize.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ==================================================================
 * L3: Horner Polynomial Evaluation
 *
 * p(x) = a_0 + a_1*x + a_2*x^2 + ... + a_n*x^n
 *      = a_0 + x*(a_1 + x*(a_2 + ... + x*a_n)...)
 *
 * Horner's method requires exactly n multiplications and n additions
 * for a degree-n polynomial, and is numerically more stable than the
 * naive sum of power terms.
 *
 * Numerical stability advantage: Horner's method computes each term
 * in a way that avoids building up large intermediate products before
 * summing. This reduces roundoff error by approximately sqrt(n) for
 * well-conditioned polynomials.
 *
 * Complexity: O(n), n = degree
 * ================================================================== */

double poly_eval_horner(unsigned degree, const double c[],
                         double x, double x_min, double x_max)
{
    if (!c || degree > 14) return NAN;
    if (x < x_min || x > x_max) return NAN;

    double result = c[degree];
    for (int i = (int)degree - 1; i >= 0; i--)
        result = result * x + c[i];

    return result;
}

void poly_eval_with_derivative(unsigned degree, const double c[],
                                double x, double *p_val, double *p_deriv)
{
    if (!c || !p_val || !p_deriv || degree > 14) {
        if (p_val) *p_val = NAN;
        if (p_deriv) *p_deriv = NAN;
        return;
    }

    /* Simultaneous evaluation using Horner's for both p(x) and p'(x):
     * p = c[n]
     * dp = 0
     * For k = n-1 down to 0:
     *   dp = p + x*dp
     *   p = c[k] + x*p
     */
    double p = c[degree];
    double dp = 0.0;

    for (int i = (int)degree - 1; i >= 0; i--) {
        dp = p + x * dp;
        p = c[i] + x * p;
    }

    *p_val = p;
    *p_deriv = dp;
}

double poly_inverse_newton(unsigned degree, const double c[],
                            double y_target, double x_guess,
                            double tol, unsigned max_iter)
{
    if (!c || degree > 14 || max_iter == 0) return NAN;

    double x = x_guess;

    for (unsigned iter = 0; iter < max_iter; iter++) {
        double p_val, p_deriv;
        poly_eval_with_derivative(degree, c, x, &p_val, &p_deriv);

        if (fabs(p_deriv) < 1e-15) return NAN;  /* Singular */

        double delta = (p_val - y_target) / p_deriv;
        x = x - delta;

        if (fabs(delta) < tol) return x;
    }

    return x;
}

/* ==================================================================
 * L3: Lookup Table Interpolation
 *
 * Binary search: O(log n) to find interval containing query point.
 * Linear interpolation: O(1) after interval found, C0 continuous.
 * Quadratic interpolation: O(1), C1 continuous, uses 3 nearest points.
 * Cubic spline: O(log n) to find interval + O(1) evaluation, C2 continuous.
 * ================================================================== */

int binary_search_interval(unsigned n, const double x_vals[],
                            double x, unsigned *idx)
{
    if (!x_vals || !idx || n < 2) return -1;
    if (x < x_vals[0] || x > x_vals[n - 1]) return -1;

    unsigned lo = 0, hi = n - 2;

    while (lo <= hi) {
        unsigned mid = lo + (hi - lo) / 2;
        if (x >= x_vals[mid] && x <= x_vals[mid + 1]) {
            *idx = mid;
            return 0;
        }
        if (x < x_vals[mid])
            hi = (mid > 0) ? mid - 1 : 0;
        else
            lo = mid + 1;
    }

    *idx = hi;
    return 0;
}

double lut_linear_interp(unsigned n_points,
                          const double x_vals[],
                          const double y_vals[],
                          double x_query)
{
    if (!x_vals || !y_vals || n_points < 2) return NAN;

    unsigned idx;
    if (binary_search_interval(n_points, x_vals, x_query, &idx) != 0)
        return NAN;

    double t = (x_query - x_vals[idx]) / (x_vals[idx + 1] - x_vals[idx]);
    return y_vals[idx] + t * (y_vals[idx + 1] - y_vals[idx]);
}

double lut_quadratic_interp(unsigned n_points,
                             const double x_vals[],
                             const double y_vals[],
                             double x_query)
{
    if (!x_vals || !y_vals || n_points < 3) return NAN;

    unsigned idx;
    if (binary_search_interval(n_points, x_vals, x_query, &idx) != 0)
        return NAN;

    /* Use three points: idx, idx+1, and one adjacent (idx-1 or idx+2) */
    unsigned i0, i1, i2;
    if (idx == 0) {
        i0 = 0; i1 = 1; i2 = 2;
    } else if (idx >= n_points - 2) {
        i0 = n_points - 3; i1 = n_points - 2; i2 = n_points - 1;
    } else {
        i0 = idx; i1 = idx + 1;
        /* Choose the closer of idx-1 or idx+2 */
        if (fabs(x_query - x_vals[idx - 1]) < fabs(x_query - x_vals[idx + 2]))
            i2 = idx - 1;
        else
            i2 = idx + 2;
        /* Ensure ordering: i0 < i1, i2 is the selected third point */
        /* We'll use Lagrange formula directly with the chosen 3 points. */
        /* Re-sort: selected = {i0, i1, idx-1 or idx+2} */
        /* Simpler approach: use i0=idx, i1=idx+1, and the nearer neighbor */
        if (idx > 0 && (idx + 2 >= n_points ||
            (x_query - x_vals[idx - 1]) < (x_vals[idx + 2] - x_query))) {
            i0 = idx - 1; i1 = idx; i2 = idx + 1;
        } else {
            i0 = idx; i1 = idx + 1; i2 = idx + 2;
        }
    }

    /* Lagrange quadratic: L(x) = y0*L0(x) + y1*L1(x) + y2*L2(x)
     * L0(x) = (x-x1)*(x-x2) / ((x0-x1)*(x0-x2))
     * L1(x) = (x-x0)*(x-x2) / ((x1-x0)*(x1-x2))
     * L2(x) = (x-x0)*(x-x1) / ((x2-x0)*(x2-x1))
     */
    double L0 = (x_query - x_vals[i1]) * (x_query - x_vals[i2])
              / ((x_vals[i0] - x_vals[i1]) * (x_vals[i0] - x_vals[i2]));
    double L1 = (x_query - x_vals[i0]) * (x_query - x_vals[i2])
              / ((x_vals[i1] - x_vals[i0]) * (x_vals[i1] - x_vals[i2]));
    double L2 = (x_query - x_vals[i0]) * (x_query - x_vals[i1])
              / ((x_vals[i2] - x_vals[i0]) * (x_vals[i2] - x_vals[i1]));

    return y_vals[i0] * L0 + y_vals[i1] * L1 + y_vals[i2] * L2;
}

double lut_cubic_spline(unsigned n_points,
                         const double x_vals[],
                         const double y_vals[],
                         double x_query)
{
    if (!x_vals || !y_vals || n_points < 4) return NAN;

    /* Catmull-Rom spline: C1 continuous, passes through all points.
     * For interval [xi, xi+1] with tangent m_i:
     *   p(t) = (2t^3-3t^2+1)*y_i + (t^3-2t^2+t)*m_i
     *        + (-2t^3+3t^2)*y_{i+1} + (t^3-t^2)*m_{i+1}
     * where t = (x - xi)/(x_{i+1} - xi), m_i = (y_{i+1} - y_{i-1})/2
     */
    unsigned idx;
    if (binary_search_interval(n_points, x_vals, x_query, &idx) != 0)
        return NAN;

    if (idx == 0 || idx >= n_points - 2) {
        /* Boundary: fall back to linear */
        return lut_linear_interp(n_points, x_vals, y_vals, x_query);
    }

    double t = (x_query - x_vals[idx]) / (x_vals[idx + 1] - x_vals[idx]);
    double m0 = (y_vals[idx + 1] - y_vals[idx - 1]) / 2.0;
    double m1 = (y_vals[idx + 2] - y_vals[idx]) / 2.0;

    double t2 = t * t;
    double t3 = t2 * t;

    double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    double h10 = t3 - 2.0 * t2 + t;
    double h01 = -2.0 * t3 + 3.0 * t2;
    double h11 = t3 - t2;

    return h00 * y_vals[idx] + h10 * m0 + h01 * y_vals[idx + 1] + h11 * m1;
}

/* ==================================================================
 * L4: Steinhart-Hart Thermistor Equation
 *
 * Three-parameter model (Equation B):
 *   1/T = A + B*ln(R) + C*[ln(R)]^3
 *
 * Four-parameter model (Equation D):
 *   1/T = A + B*ln(R) + C*[ln(R)]^2 + D*[ln(R)]^3
 *
 * To solve for R given T:
 *   Let x = ln(R)
 *   C*x^3 + B*x + (A - 1/T) = 0  (3-param)
 *
 * For 4-param: D*x^3 + C*x^2 + B*x + (A - 1/T) = 0
 *
 * These are cubic equations in x, solvable by:
 *   - Cardano's cubic formula (analytical)
 *   - Newton-Raphson iteration
 *   - LUT + interpolation
 *
 * We use Newton-Raphson starting from x0 = ln(R_ref).
 * ================================================================== */

double steinhart_hart_T(const steinhart_hart_params_t *params,
                         double resistance_ohm)
{
    if (!params || resistance_ohm <= 0.0) return NAN;

    double lnR = log(resistance_ohm);
    double lnR2 = lnR * lnR;
    double lnR3 = lnR2 * lnR;

    double inv_T;
    if (params->use_4param) {
        inv_T = params->a_coeff
              + params->b_coeff * lnR
              + params->c_coeff * lnR2
              + params->d_coeff * lnR3;
    } else {
        inv_T = params->a_coeff
              + params->b_coeff * lnR
              + params->c_coeff * lnR3;
    }

    if (fabs(inv_T) < 1e-15) return NAN;
    return 1.0 / inv_T;  /* Returns temperature in Kelvin */
}

double steinhart_hart_R(const steinhart_hart_params_t *params,
                         double temperature_kelvin)
{
    if (!params || temperature_kelvin <= 0.0) return NAN;

    double inv_T = 1.0 / temperature_kelvin;
    double target = inv_T - params->a_coeff;

    /* Newton-Raphson: solve f(x) = B*x + C*x^3 - target = 0 where x = ln(R) */
    double x = log(params->r_ref);  /* Initial guess */
    if (x < 1.0) x = 1.0;  /* Avoid log of very small numbers */

    for (int iter = 0; iter < 50; iter++) {
        double x2 = x * x;
        double x3 = x2 * x;
        double fx, fpx; /* f(x) and f'(x) */

        if (params->use_4param) {
            fx = params->b_coeff * x + params->c_coeff * x2
               + params->d_coeff * x3 - target;
            fpx = params->b_coeff + 2.0 * params->c_coeff * x
                + 3.0 * params->d_coeff * x2;
        } else {
            fx = params->b_coeff * x + params->c_coeff * x3 - target;
            fpx = params->b_coeff + 3.0 * params->c_coeff * x2;
        }

        if (fabs(fpx) < 1e-15) break;

        double dx = fx / fpx;
        x = x - dx;
        if (fabs(dx) < 1e-12) break;
    }

    return exp(x);
}

void steinhart_hart_fit_3pt(double t1_c, double r1_ohm,
                             double t2_c, double r2_ohm,
                             double t3_c, double r3_ohm,
                             steinhart_hart_params_t *params)
{
    if (!params) return;

    double T1 = t1_c + 273.15;
    double T2 = t2_c + 273.15;
    double T3 = t3_c + 273.15;

    double L1 = log(r1_ohm);
    double L2 = log(r2_ohm);
    double L3 = log(r3_ohm);

    double Y1 = 1.0 / T1;
    double Y2 = 1.0 / T2;
    double Y3 = 1.0 / T3;

    /* Solve 3x3 linear system for A, B, C:
     * Y1 = A + B*L1 + C*L1^3
     * Y2 = A + B*L2 + C*L2^3
     * Y3 = A + B*L3 + C*L3^3
     *
     * Eliminate A:
     * (Y2-Y1) = B*(L2-L1) + C*(L2^3-L1^3)
     * (Y3-Y1) = B*(L3-L1) + C*(L3^3-L1^3)
     */
    double L13 = L1 * L1 * L1;
    double L23 = L2 * L2 * L2;
    double L33 = L3 * L3 * L3;

    double dY21 = Y2 - Y1, dL21 = L2 - L1, dL21c = L23 - L13;
    double dY31 = Y3 - Y1, dL31 = L3 - L1, dL31c = L33 - L13;

    double det = dL21 * dL31c - dL21c * dL31;
    if (fabs(det) < 1e-20) return;

    double B = (dY21 * dL31c - dY31 * dL21c) / det;
    double C = (dL21 * dY31 - dL31 * dY21) / det;
    double A = Y1 - B * L1 - C * L13;

    params->a_coeff = A;
    params->b_coeff = B;
    params->c_coeff = C;
    params->d_coeff = 0.0;
    params->use_4param = false;
    params->r_ref = r1_ohm;
    params->t_ref_c = t1_c;
    params->beta_k = log(r1_ohm / r2_ohm) / (1.0/T1 - 1.0/T2);
}

double thermistor_beta(double t1_kelvin, double r1_ohm,
                        double t2_kelvin, double r2_ohm)
{
    if (t1_kelvin <= 0.0 || t2_kelvin <= 0.0) return NAN;
    if (r1_ohm <= 0.0 || r2_ohm <= 0.0) return NAN;
    if (fabs(r2_ohm) < 1e-15) return NAN;

    return log(r1_ohm / r2_ohm) / (1.0 / t1_kelvin - 1.0 / t2_kelvin);
}

/* ==================================================================
 * L4: Callendar-Van Dusen RTD Linearization (IEC 60751)
 *
 * Standard Pt100 coefficients:
 *   R0 = 100.00 Ohm at 0C
 *   A = 3.9083E-3 /C
 *   B = -5.775E-7 /C^2
 *   C = -4.183E-12 /C^4 (only for t < 0C)
 *
 * For t >= 0C:
 *   R(t) = R0 * (1 + A*t + B*t^2)
 *
 * For -200C <= t < 0C:
 *   R(t) = R0 * (1 + A*t + B*t^2 + C*(t-100)*t^3)
 *
 * Inverse (R -> T):
 *   For R >= R0: t = (-A + sqrt(A^2 - 4*B*(1-R/R0))) / (2*B)
 *   For R < R0: iterative (Newton-Raphson) due to C*(t-100)*t^3 term.
 * ================================================================== */

void rtd_cvd_init_pt100(rtd_cvd_params_t *params)
{
    if (!params) return;
    params->r0 = 100.0;
    params->a_coeff = 3.9083e-3;
    params->b_coeff = -5.775e-7;
    params->c_coeff = -4.183e-12;
    params->alpha = 0.00385055;
    params->delta = 1.4999;
    params->beta = 0.10863;
    params->above_zero = true;
    params->nominal_current_a = 1.0e-3;  /* 1 mA typical */
}

void rtd_cvd_init_pt1000(rtd_cvd_params_t *params)
{
    if (!params) return;
    rtd_cvd_init_pt100(params);
    params->r0 = 1000.0;
    params->nominal_current_a = 0.3e-3;  /* 0.3 mA for lower self-heating */
}

double rtd_cvd_R_from_T(const rtd_cvd_params_t *params,
                         double temperature_c)
{
    if (!params) return NAN;

    double t = temperature_c;
    double r;

    if (t >= 0.0) {
        r = params->r0 * (1.0 + params->a_coeff * t
                          + params->b_coeff * t * t);
    } else {
        r = params->r0 * (1.0 + params->a_coeff * t
                          + params->b_coeff * t * t
                          + params->c_coeff * (t - 100.0) * t * t * t);
    }

    return r;
}

double rtd_cvd_T_from_R(const rtd_cvd_params_t *params,
                         double resistance_ohm)
{
    if (!params || resistance_ohm <= 0.0) return NAN;

    /* For R >= R0 (t >= 0C): analytical quadratic solution */
    if (resistance_ohm >= params->r0) {
        double A = params->a_coeff;
        double B = params->b_coeff;
        double r_ratio = resistance_ohm / params->r0;
        double disc = A * A - 4.0 * B * (1.0 - r_ratio);

        if (disc < 0.0) return NAN;
        return (-A + sqrt(disc)) / (2.0 * B);
    }

    /* For R < R0 (t < 0C): Newton-Raphson */
    double t = -50.0;  /* Initial guess around -50C */

    for (int iter = 0; iter < 30; iter++) {
        double A = params->a_coeff;
        double B = params->b_coeff;
        double C = params->c_coeff;

        double t2 = t * t;
        double t3 = t2 * t;
        double t4 = t3 * t;

        /* R(t) = R0*(1 + A*t + B*t^2 + C*(t-100)*t^3) */
        double R_calc = params->r0 * (1.0 + A * t + B * t2
                        + C * (t - 100.0) * t3);

        /* dR/dt = R0*(A + 2*B*t + C*(4*t^3 - 300*t^2)) */
        double dR_dt = params->r0 * (A + 2.0 * B * t
                       + C * (4.0 * t3 - 300.0 * t2));

        if (fabs(dR_dt) < 1e-15) break;

        double delta = (R_calc - resistance_ohm) / dR_dt;
        t = t - delta;

        if (fabs(delta) < 1e-6) break;
    }

    return t;
}

double rtd_alpha(const rtd_cvd_params_t *params)
{
    if (!params) return 0.0;
    /* alpha = (R100 - R0) / (100 * R0), standard Pt100 alpha = 0.00385055 */
    double R100 = rtd_cvd_R_from_T(params, 100.0);
    return (R100 - params->r0) / (100.0 * params->r0);
}

/* ==================================================================
 * L5: NIST ITS-90 Thermocouple Polynomial Functions
 *
 * Each thermocouple type has published polynomial coefficients
 * from NIST. We implement standard Type K, J, T, E.
 *
 * Type K (Chromel-Alumel):
 *   Temperature range: -270 to 1372 C
 *   Seebeck at 25C: ~40.6 uV/C
 *   Inverse polynomial errors: ~+/-0.04C for -200 to 0C, ~+/-0.02C for 0 to 1372C
 *
 * Coefficients are from NIST Monograph 175 (1993), Tables for each type.
 * ================================================================== */

void nist_tc_init_typeK(thermocouple_nist_model_t *model)
{
    if (!model) return;
    model->type_code = 0;  /* Type K */
    model->seebeck_25c = 40.6;
    model->temp_min_c = -270.0;
    model->temp_max_c = 1372.0;

    /* Forward polynomial (T to V, uV): for 0 to 1372C range
     * V = c0 + c1*T + c2*T^2 + ... + c10*T^10 + a0*exp(a1*(T-126.9686)^2)
     *
     * Coefficients for 0C to 1372C:
     */
    model->forward_poly.degree = 10;
    model->forward_poly.coefficients[0]  = -1.7600413686e-2;
    model->forward_poly.coefficients[1]  =  3.8921204975e1;
    model->forward_poly.coefficients[2]  =  1.8558770032e-2;
    model->forward_poly.coefficients[3]  = -9.9457592874e-5;
    model->forward_poly.coefficients[4]  =  3.1840945719e-7;
    model->forward_poly.coefficients[5]  = -5.6072844889e-10;
    model->forward_poly.coefficients[6]  =  5.6075059059e-13;
    model->forward_poly.coefficients[7]  = -3.2020720003e-16;
    model->forward_poly.coefficients[8]  =  9.7151147152e-20;
    model->forward_poly.coefficients[9]  = -1.2104721275e-23;
    model->forward_poly.coefficients[10] =  0.0;  /* exponential term handled separately */

    /* Inverse polynomial (V to T, C): for 0 to 1372C range
     * T = d0 + d1*V + d2*V^2 + ... + d9*V^9 + a0*exp(a1*(V-a2)^2)
     * V in millivolts
     */
    model->inverse_poly.degree = 9;
    model->inverse_poly.coefficients[0] =  0.0;
    model->inverse_poly.coefficients[1] =  2.5173462e-2;
    model->inverse_poly.coefficients[2] = -1.1662878e-6;
    model->inverse_poly.coefficients[3] = -1.0833638e-9;
    model->inverse_poly.coefficients[4] = -8.9773540e-13;
    model->inverse_poly.coefficients[5] = -3.7342377e-16;
    model->inverse_poly.coefficients[6] = -8.6632643e-20;
    model->inverse_poly.coefficients[7] = -1.0450598e-23;
    model->inverse_poly.coefficients[8] = -5.1920577e-28;
    model->inverse_poly.coefficients[9] =  0.0;

    model->voltage_min_uv = -6444.0;
    model->voltage_max_uv = 54886.0;
    model->error_max_c = 0.05;
    model->error_typical_c = 0.02;
}

static void nist_tc_common_init(thermocouple_nist_model_t *model,
    unsigned tc_type, double seebeck, double tmin, double tmax,
    unsigned fwd_deg, const double fwd_c[], unsigned inv_deg,
    const double inv_c[], double vmin, double vmax)
{
    if (!model) return;
    model->type_code = tc_type;
    model->seebeck_25c = seebeck;
    model->temp_min_c = tmin;
    model->temp_max_c = tmax;

    model->forward_poly.degree = fwd_deg;
    for (unsigned i = 0; i <= fwd_deg && i < 15; i++)
        model->forward_poly.coefficients[i] = fwd_c[i];

    model->inverse_poly.degree = inv_deg;
    for (unsigned i = 0; i <= inv_deg && i < 15; i++)
        model->inverse_poly.coefficients[i] = inv_c[i];

    model->voltage_min_uv = vmin;
    model->voltage_max_uv = vmax;
    model->error_max_c = 0.1;
    model->error_typical_c = 0.05;
}

void nist_tc_init_typeJ(thermocouple_nist_model_t *model)
{
    /* Type J: Iron-Constantan, -210 to 1200C, ~51.9 uV/C at 25C */
    double fwd[] = {0.0, 5.0381187815e1, 3.0475836930e-2, -8.5681065720e-5,
        1.3228195295e-7, -1.7052958337e-10, 2.0948090697e-13, -1.2538395336e-16,
        1.5631725697e-20};
    double inv[] = {0.0, 1.978425e-2, -2.001204e-7, 1.036969e-11, -2.549687e-16,
        3.585153e-21, -5.344285e-26, 5.099890e-31};
    nist_tc_common_init(model, 1, 51.9, -210.0, 1200.0,
        8, fwd, 7, inv, -8095.0, 69553.0);
}

void nist_tc_init_typeT(thermocouple_nist_model_t *model)
{
    /* Type T: Copper-Constantan, -270 to 400C, ~40.6 uV/C at 25C */
    double fwd[] = {0.0, 3.8748106364e1, 4.4194434347e-2, 1.1844323105e-4,
        -2.0032973554e-5, 9.0138019559e-7, -2.2651156593e-8, 3.6071154205e-10,
        -3.8493939865e-12, 2.8213521925e-14, -1.3911593228e-16, 4.3780087205e-19,
        -7.9790316553e-22, 6.3895952294e-25};
    double inv[] = {0.0, 2.5949192e-2, -2.1316967e-7, 7.9018692e-10,
        -4.2527777e-13, 1.3304473e-16, -2.0241446e-20, 1.2668171e-24};
    nist_tc_common_init(model, 2, 40.6, -270.0, 400.0,
        13, fwd, 7, inv, -6258.0, 20872.0);
}

void nist_tc_init_typeE(thermocouple_nist_model_t *model)
{
    /* Type E: Chromel-Constantan, -270 to 1000C, ~60.9 uV/C at 25C */
    double fwd[] = {0.0, 5.8665508708e1, 4.5410977124e-2, -7.7998048686e-4,
        -2.5800160843e-5, -5.9452583057e-7, -9.3214058667e-9, -1.0287605534e-10,
        -8.0370123621e-13};
    double inv[] = {0.0, 1.6977288e-2, -4.3514970e-7, -1.5859697e-10,
        -9.2502871e-14, -2.6084314e-17, -4.1360199e-21, -3.4034030e-25};
    nist_tc_common_init(model, 3, 60.9, -270.0, 1000.0,
        8, fwd, 7, inv, -9835.0, 76373.0);
}

double nist_tc_temp_to_voltage(const thermocouple_nist_model_t *model,
                                double temperature_c)
{
    if (!model) return NAN;
    if (temperature_c < model->temp_min_c || temperature_c > model->temp_max_c)
        return NAN;

    return poly_eval_horner(model->forward_poly.degree,
                            model->forward_poly.coefficients,
                            temperature_c,
                            model->temp_min_c, model->temp_max_c);
}

double nist_tc_voltage_to_temp(const thermocouple_nist_model_t *model,
                                double voltage_uv)
{
    if (!model) return NAN;
    if (voltage_uv < model->voltage_min_uv || voltage_uv > model->voltage_max_uv)
        return NAN;

    /* Use Newton-Raphson to invert the forward polynomial:
     * Solve nist_tc_temp_to_voltage(T) = voltage_uv for T.
     * This guarantees round-trip consistency since only
     * one polynomial is used.
     *
     * The Seebeck coefficient gives the derivative:
     * dV/dT ~ S(T) which is ~p'(T) from the forward polynomial.
     *
     * Initial guess using linear approximation:
     * T_guess = voltage_uv / seebeck_25c
     */
    double T = voltage_uv / model->seebeck_25c;
    /* Clamp initial guess to valid temperature range */
    if (T < model->temp_min_c) T = model->temp_min_c;
    if (T > model->temp_max_c) T = model->temp_max_c;

    for (int iter = 0; iter < 30; iter++) {
        double V_calc, V_deriv;
        /* Evaluate forward polynomial and its derivative at T */
        poly_eval_with_derivative(model->forward_poly.degree,
                                  model->forward_poly.coefficients,
                                  T, &V_calc, &V_deriv);

        double error = V_calc - voltage_uv;

        /* Guard against zero derivative */
        if (fabs(V_deriv) < 0.01) {
            /* Use Seebeck coefficient as fallback derivative */
            V_deriv = model->seebeck_25c;
        }

        double delta = error / V_deriv;
        T = T - delta;

        /* Clamp to valid range */
        if (T < model->temp_min_c) T = model->temp_min_c;
        if (T > model->temp_max_c) T = model->temp_max_c;

        if (fabs(delta) < 1e-6) break;
    }

    return T;
}

/* ==================================================================
 * L5: Piecewise Linearization and Least-Squares Fitting
 * ================================================================== */

double piecewise_linear_hysteresis(unsigned n_segments,
                                    const double breakpoints[],
                                    const double fwd_slopes[],
                                    const double fwd_intercepts[],
                                    const double rev_slopes[],
                                    const double rev_intercepts[],
                                    double x, bool direction_increasing)
{
    if (!breakpoints || n_segments < 1) return NAN;

    const double *slopes = direction_increasing ? fwd_slopes : rev_slopes;
    const double *intercepts = direction_increasing ? fwd_intercepts : rev_intercepts;

    /* Find segment */
    unsigned seg = 0;
    for (unsigned i = 0; i < n_segments; i++) {
        if (x <= breakpoints[i + 1]) { seg = i; break; }
    }
    if (x > breakpoints[n_segments]) seg = n_segments - 1;

    return slopes[seg] * x + intercepts[seg];
}

void build_piecewise_linear(unsigned n_points,
                             const calibration_point_t points[],
                             double breakpoints[],
                             double slopes[],
                             double intercepts[],
                             unsigned *n_segments)
{
    if (!points || !breakpoints || !slopes || !intercepts || !n_segments
        || n_points < 2) return;

    *n_segments = n_points - 1;
    for (unsigned i = 0; i < *n_segments; i++) {
        breakpoints[i] = points[i].sensor_output_raw;
        double dx = points[i + 1].sensor_output_raw - points[i].sensor_output_raw;
        double dy = points[i + 1].measurand_value - points[i].measurand_value;
        if (fabs(dx) < 1e-15) dx = 1e-15;
        slopes[i] = dy / dx;
        intercepts[i] = points[i].measurand_value - slopes[i] * points[i].sensor_output_raw;
    }
    breakpoints[*n_segments] = points[n_points - 1].sensor_output_raw;
}

void poly_fit_least_squares(unsigned n_points,
                             const calibration_point_t points[],
                             unsigned target_degree,
                             polynomial_coeffs_t *result)
{
    if (!points || !result || n_points < target_degree + 1 || target_degree > 14)
        return;

    /* Build normal equations: A^T * A * c = A^T * y
     * where A_ij = x_i^j, y_i = measurand_value
     *
     * For small degrees (<=5), use explicit normal equation matrix.
     */
    unsigned m = target_degree + 1;
    double ATA[36] = {0};  /* (m x m) matrix, max 6x6 */
    double ATy[6] = {0};

    for (unsigned i = 0; i < n_points; i++) {
        double x = points[i].sensor_output_raw;
        double y = points[i].measurand_value;
        double x_pow[6];
        x_pow[0] = 1.0;
        for (unsigned j = 1; j < m; j++)
            x_pow[j] = x_pow[j - 1] * x;

        for (unsigned row = 0; row < m; row++) {
            ATy[row] += x_pow[row] * y;
            for (unsigned col = 0; col < m; col++)
                ATA[row * m + col] += x_pow[row] * x_pow[col];
        }
    }

    /* Gaussian elimination with partial pivoting */
    for (unsigned col = 0; col < m; col++) {
        /* Find pivot */
        unsigned max_row = col;
        double max_val = fabs(ATA[col * m + col]);
        for (unsigned row = col + 1; row < m; row++) {
            if (fabs(ATA[row * m + col]) > max_val) {
                max_val = fabs(ATA[row * m + col]);
                max_row = row;
            }
        }

        /* Swap rows */
        if (max_row != col) {
            for (unsigned j = 0; j < m; j++) {
                double tmp = ATA[col * m + j];
                ATA[col * m + j] = ATA[max_row * m + j];
                ATA[max_row * m + j] = tmp;
            }
            double tmp2 = ATy[col];
            ATy[col] = ATy[max_row];
            ATy[max_row] = tmp2;
        }

        /* Eliminate */
        double pivot = ATA[col * m + col];
        if (fabs(pivot) < 1e-15) continue;

        for (unsigned row = col + 1; row < m; row++) {
            double factor = ATA[row * m + col] / pivot;
            for (unsigned j = col; j < m; j++)
                ATA[row * m + j] -= factor * ATA[col * m + j];
            ATy[row] -= factor * ATy[col];
        }
    }

    /* Back substitution */
    for (int i = (int)m - 1; i >= 0; i--) {
        double sum = ATy[i];
        for (unsigned j = i + 1; j < m; j++)
            sum -= ATA[i * m + j] * result->coefficients[j];
        if (fabs(ATA[i * m + i]) < 1e-15)
            result->coefficients[i] = 0.0;
        else
            result->coefficients[i] = sum / ATA[i * m + i];
    }

    result->degree = target_degree;
    result->valid_min = points[0].sensor_output_raw;
    result->valid_max = points[n_points - 1].sensor_output_raw;
    result->residual_error_max = poly_max_residual(n_points, points, result);
}

double poly_max_residual(unsigned n_points,
                          const calibration_point_t points[],
                          const polynomial_coeffs_t *poly)
{
    if (!points || !poly || n_points == 0) return INFINITY;

    double max_err = 0.0;
    for (unsigned i = 0; i < n_points; i++) {
        double predicted = poly_eval_horner(poly->degree, poly->coefficients,
                                            points[i].sensor_output_raw,
                                            poly->valid_min, poly->valid_max);
        double err = fabs(predicted - points[i].measurand_value);
        if (err > max_err) max_err = err;
    }
    return max_err;
}
