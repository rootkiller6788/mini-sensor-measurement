/**
 * @file sensor_transfer.c
 * @brief Sensor transfer function implementations
 *
 * Levels: L2 Core Concepts — polynomial, thermistor, RTD, thermocouple transfer
 *         L3 Mathematical Structures — dynamic model ODE solutions
 *         L4 Fundamental Laws — Seebeck, Steinhart-Hart, Callendar-Van Dusen
 */

#include "sensor_transfer.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ═══════════════════════════════════════════════════════════════════════
 * L2: Polynomial Transfer Function
 * ═══════════════════════════════════════════════════════════════════════ */

int poly_transfer_init(poly_transfer_t *pt, const double *coeff, uint8_t order,
                        double x_min, double x_max)
{
    if (!pt || !coeff || order > 15 || x_min >= x_max) return -1;

    pt->order = order;
    pt->x_min = x_min;
    pt->x_max = x_max;
    pt->coeff = (double *)calloc(order + 1, sizeof(double));
    if (!pt->coeff) return -1;
    for (int i = 0; i <= order; i++) {
        pt->coeff[i] = coeff[i];
    }
    return 0;
}

void poly_transfer_free(poly_transfer_t *pt)
{
    if (pt && pt->coeff) {
        free(pt->coeff);
        pt->coeff = NULL;
        pt->order = 0;
    }
}

/**
 * @brief Horner's method for efficient polynomial evaluation.
 *
 * Complexity: O(n) multiplications, O(n) additions
 *
 * y = a₀ + x·(a₁ + x·(a₂ + x·(... + x·aₙ)))
 */
double poly_transfer_evaluate(const poly_transfer_t *pt, double x)
{
    if (!pt || !pt->coeff) return NAN;
    double y = 0.0;
    for (int i = pt->order; i >= 0; i--) {
        y = y * x + pt->coeff[i];
    }
    return y;
}

/**
 * @brief Inverse polynomial evaluation using Newton-Raphson iteration.
 *
 * Solve f(x) - y = 0 for x using:
 *   x_{k+1} = x_k - [f(x_k) - y] / f'(x_k)
 *
 * Convergence is quadratic near the root for well-behaved monotonic functions.
 * Initial guess uses linear approximation: x₀ = (y - a₀) / a₁.
 */
double poly_transfer_inverse(const poly_transfer_t *pt, double y, double tol)
{
    if (!pt || !pt->coeff || pt->order < 1 || tol <= 0.0) return NAN;

    /* Initial guess: linear approximation */
    double x = (y - pt->coeff[0]) / pt->coeff[1];
    if (x < pt->x_min) x = pt->x_min;
    if (x > pt->x_max) x = pt->x_max;

    for (int iter = 0; iter < 50; iter++) {
        /* Evaluate f(x) and f'(x) simultaneously using Horner */
        double f = 0.0, df = 0.0;
        for (int i = pt->order; i >= 0; i--) {
            df = df * x + f;
            f = f * x + pt->coeff[i];
        }
        /* f is now f(x), but we used Horner with accumulation.
         * Correct approach: separate passes for f and f' */
        /* Recompute properly */
        f = pt->coeff[pt->order];
        df = 0.0;
        for (int i = pt->order - 1; i >= 0; i--) {
            df = df * x + f;
            f = f * x + pt->coeff[i];
        }
        double residual = f - y;
        if (fabs(residual) < tol) return x;
        if (fabs(df) < 1e-30) break;
        x = x - residual / df;
        if (x < pt->x_min) x = pt->x_min;
        if (x > pt->x_max) x = pt->x_max;
    }
    return x;
}

double poly_transfer_sensitivity(const poly_transfer_t *pt, double x)
{
    if (!pt || !pt->coeff || pt->order < 1) return NAN;

    /* f'(x) = Σ_{i=1}^{n} i·a_i·x^{i-1} */
    double df = 0.0;
    double x_pow = 1.0;
    for (int i = 1; i <= pt->order; i++) {
        df += i * pt->coeff[i] * x_pow;
        x_pow *= x;
    }
    return df;
}

double poly_transfer_nonlinearity(const poly_transfer_t *pt, size_t n_samples)
{
    if (!pt || !pt->coeff || n_samples < 2) return 0.0;

    /* Best-fit line through endpoints
     * m = (f(x_max) - f(x_min)) / (x_max - x_min)
     * b = f(x_min) - m·x_min */
    double f_min = poly_transfer_evaluate(pt, pt->x_min);
    double f_max = poly_transfer_evaluate(pt, pt->x_max);
    double m = (f_max - f_min) / (pt->x_max - pt->x_min);
    double b = f_min - m * pt->x_min;
    double full_scale = f_max - f_min;

    double max_dev = 0.0;
    for (size_t i = 0; i <= n_samples; i++) {
        double x = pt->x_min + (pt->x_max - pt->x_min) * i / n_samples;
        double f_x = poly_transfer_evaluate(pt, x);
        double line = m * x + b;
        double dev = fabs(f_x - line);
        if (dev > max_dev) max_dev = dev;
    }

    if (full_scale < 1e-30) return 0.0;
    return max_dev / full_scale * 100.0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L2+L4: Thermistor Transfer Function (Steinhart-Hart, 1968)
 * ═══════════════════════════════════════════════════════════════════════ */

int thermistor_model_init_full(thermistor_model_t *tm,
                                double A, double B, double C)
{
    if (!tm) return -1;
    tm->A = A;
    tm->B = B;
    tm->C = C;
    tm->R_at_T0 = 0.0;
    tm->T0 = 0.0;
    tm->B_param = 0.0;
    tm->use_beta = 0;
    return 0;
}

int thermistor_model_init_beta(thermistor_model_t *tm,
                                double R0, double T0_K, double beta)
{
    if (!tm || R0 <= 0.0 || T0_K <= 0.0 || beta <= 0.0) return -1;
    tm->R_at_T0 = R0;
    tm->T0 = T0_K;
    tm->B_param = beta;
    tm->use_beta = 1;
    tm->A = tm->B = tm->C = 0.0;
    return 0;
}

/**
 * @brief Compute temperature from thermistor resistance.
 *
 * Full Steinhart-Hart: 1/T = A + B·ln(R) + C·[ln(R)]³
 *
 * β-model: T = 1 / [1/T₀ + (1/β)·ln(R/R₀)]
 */
double thermistor_temp_from_R(const thermistor_model_t *tm, double R)
{
    if (!tm || R <= 0.0) return NAN;

    if (tm->use_beta) {
        double inv_T = 1.0 / tm->T0 + (1.0 / tm->B_param) * log(R / tm->R_at_T0);
        if (inv_T <= 0.0) return NAN;
        return 1.0 / inv_T;
    } else {
        double lnR = log(R);
        double inv_T = tm->A + tm->B * lnR + tm->C * lnR * lnR * lnR;
        if (inv_T <= 0.0) return NAN;
        return 1.0 / inv_T;
    }
}

/**
 * @brief Compute thermistor resistance from temperature.
 *
 * Full S-H: Solve cubic C·[ln(R)]³ + B·[ln(R)] + (A - 1/T) = 0
 *
 * This cubic in ln(R) has the special form with no quadratic term
 * (depressed cubic). Use Cardano's formula or Newton's method.
 *
 * For the β-model: R = R₀ · exp[β·(1/T - 1/T₀)]
 */
double thermistor_R_from_temp(const thermistor_model_t *tm, double T_K)
{
    if (!tm || T_K <= 0.0) return NAN;

    if (tm->use_beta) {
        return tm->R_at_T0 * exp(tm->B_param * (1.0 / T_K - 1.0 / tm->T0));
    } else {
        /* Solve C·X³ + B·X + (A - 1/T) = 0 for X = ln(R) using Newton */
        double target = 1.0 / T_K - tm->A;
        /* Initial guess: linear term dominates */
        double X = target / tm->B;
        for (int iter = 0; iter < 30; iter++) {
            double X2 = X * X;
            double f = tm->C * X2 * X + tm->B * X - target;
            double df = 3.0 * tm->C * X2 + tm->B;
            if (fabs(df) < 1e-30) break;
            double dX = f / df;
            X -= dX;
            if (fabs(dX) < 1e-12) break;
        }
        return exp(X);
    }
}

/**
 * @brief Sensitivity of thermistor: dR/dT.
 *
 * For NTC thermistors, dR/dT is negative and typically -2% to -6% per °C.
 */
double thermistor_sensitivity(const thermistor_model_t *tm, double T_K)
{
    if (!tm || T_K <= 0.0) return NAN;

    if (tm->use_beta) {
        double R = thermistor_R_from_temp(tm, T_K);
        return -R * tm->B_param / (T_K * T_K);
    } else {
        double R = thermistor_R_from_temp(tm, T_K);
        double lnR = log(R);
        double d_invT_dR = (tm->B + 3.0 * tm->C * lnR * lnR) / R;
        if (fabs(d_invT_dR) < 1e-30) return NAN;
        return -1.0 / (T_K * T_K * d_invT_dR);
    }
}

/**
 * @brief Estimate self-heating temperature rise.
 *
 * ΔT = I² · R_th · R_thermistor / dissipation_constant
 *
 * The dissipation constant δ [mW/K] is the power required to
 * raise the thermistor temperature by 1 K above ambient.
 * Typical values: 0.5–10 mW/K (depends on size, encapsulation, medium).
 */
double thermistor_self_heating(double current, double R_th,
                                double dissipation_constant_mW_K)
{
    if (dissipation_constant_mW_K <= 0.0) return 0.0;
    double power_mW = current * current * R_th * 1000.0;
    return power_mW / dissipation_constant_mW_K;
}

/**
 * @brief Compute parallel resistor Rp for linearization at temperature T_center.
 *
 * A fixed resistor Rp in parallel with the thermistor creates a combined
 * resistance with reduced nonlinearity near T_center.
 *
 * For Rp chosen such that d(R_comb)/dT has zero curvature at T_center:
 *   Rp = R(T_center) · (β - 2·T_center) / (β + 2·T_center)
 *
 * This extends the linear range from ~10°C to ~50°C for NTC thermistors.
 */
double thermistor_linearization_Rp(const thermistor_model_t *tm,
                                    double T_center, double *Rp_out)
{
    if (!tm || T_center <= 0.0) return -1.0;
    double R_center = thermistor_R_from_temp(tm, T_center);

    if (tm->use_beta) {
        double beta = tm->B_param;
        double Rp = R_center * (beta - 2.0 * T_center) / (beta + 2.0 * T_center);
        if (Rp_out) *Rp_out = Rp;
        return Rp;
    }
    /* For full S-H, approximate using local β at T_center */
    double T2 = T_center + 10.0;
    double R2 = thermistor_R_from_temp(tm, T2);
    double beta_local = log(R2 / R_center) / (1.0 / T2 - 1.0 / T_center);
    double Rp = R_center * (beta_local - 2.0 * T_center) / (beta_local + 2.0 * T_center);
    if (Rp_out) *Rp_out = Rp;
    return Rp;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L2+L4: RTD Callendar-Van Dusen Equation (IEC 60751)
 * ═══════════════════════════════════════════════════════════════════════ */

int rtd_cvd_init_iec751(rtd_cvd_model_t *rtd, double R0)
{
    if (!rtd || R0 <= 0.0) return -1;
    rtd->R0 = R0;
    /* IEC 60751 standard coefficients for α = 0.00385055 */
    rtd->A = 3.9083e-3;
    rtd->B = -5.775e-7;
    rtd->C = -4.183e-12;
    rtd->alpha = 0.00385055;
    return 0;
}

int rtd_cvd_init_custom(rtd_cvd_model_t *rtd,
                         double R0, double alpha, double delta, double beta)
{
    if (!rtd || R0 <= 0.0 || alpha <= 0.0) return -1;

    /* Callendar-Van Dusen coefficients from α, δ, β:
     * A = α + α·δ/100
     * B = -α·δ/10000
     * C = -α·β/100000000  (for T < 0°C) */
    rtd->R0 = R0;
    rtd->alpha = alpha;
    rtd->A = alpha + alpha * delta / 100.0;
    rtd->B = -alpha * delta / 10000.0;
    rtd->C = -alpha * beta / 100000000.0;
    return 0;
}

/**
 * @brief Compute RTD resistance from temperature.
 *
 * For T ≥ 0°C: R(T) = R₀·[1 + A·T + B·T²]
 * For T < 0°C: R(T) = R₀·[1 + A·T + B·T² + C·(T-100)·T³]
 *
 * The (T-100)·T³ term is zero at T=0°C and T=100°C, providing continuity.
 */
double rtd_R_from_T(const rtd_cvd_model_t *rtd, double T_C)
{
    if (!rtd) return NAN;

    double R;
    if (T_C >= 0.0) {
        R = rtd->R0 * (1.0 + rtd->A * T_C + rtd->B * T_C * T_C);
    } else {
        double T3 = T_C * T_C * T_C;
        R = rtd->R0 * (1.0 + rtd->A * T_C + rtd->B * T_C * T_C
                        + rtd->C * (T_C - 100.0) * T3);
    }
    return R;
}

/**
 * @brief Compute temperature from RTD resistance using iterative inversion.
 *
 * For T ≥ 0°C: solve quadratic B·T² + A·T + (1 - R/R₀) = 0
 *   T = [-A + √(A² - 4B·(1 - R/R₀))] / (2B)   [positive root]
 *
 * For T < 0°C: use Newton-Raphson on the full quartic.
 */
double rtd_T_from_R(const rtd_cvd_model_t *rtd, double R, double tol)
{
    if (!rtd || R <= 0.0 || tol <= 0.0) return NAN;

    /* First, check if T ≥ 0 using the quadratic formula */
    double ratio = R / rtd->R0 - 1.0;
    double discriminant = rtd->A * rtd->A - 4.0 * rtd->B * ratio;

    if (discriminant >= 0.0) {
        /* Quadratic solution for T ≥ 0 */
        double T_pos = (-rtd->A + sqrt(discriminant)) / (2.0 * rtd->B);
        if (T_pos >= 0.0) return T_pos;
    }

    /* T < 0: Newton-Raphson on quartic */
    double T = -10.0;  /* initial guess */
    for (int iter = 0; iter < 50; iter++) {
        double T2 = T * T;
        double T3 = T2 * T;
        double f = rtd->R0 * (1.0 + rtd->A * T + rtd->B * T2
                               + rtd->C * (T - 100.0) * T3) - R;
        /* f'(T) = R₀·[A + 2B·T + C·(4T³ - 300T²)] */
        double df = rtd->R0 * (rtd->A + 2.0 * rtd->B * T
                                + rtd->C * (4.0 * T3 - 300.0 * T2));
        if (fabs(df) < 1e-30) break;
        double dT = f / df;
        T -= dT;
        if (fabs(dT) < tol) break;
    }
    return T;
}

double rtd_sensitivity(const rtd_cvd_model_t *rtd, double T_C)
{
    if (!rtd) return NAN;

    /* dR/dT = R₀·[A + 2B·T + C·(4T³ - 300T²)] for T < 0
     *        = R₀·[A + 2B·T] for T ≥ 0
     */
    if (T_C >= 0.0) {
        return rtd->R0 * (rtd->A + 2.0 * rtd->B * T_C);
    } else {
        double T2 = T_C * T_C;
        double T3 = T2 * T_C;
        return rtd->R0 * (rtd->A + 2.0 * rtd->B * T_C
                           + rtd->C * (4.0 * T3 - 300.0 * T2));
    }
}

double rtd_self_heating(double I_meas, double R_element,
                         double thermal_resistance_K_W)
{
    double P = I_meas * I_meas * R_element;
    return P * thermal_resistance_K_W;
}

/**
 * @brief Correct 2-wire RTD measurement for lead resistance.
 *
 * In 2-wire configuration, R_measured = R_sensor + 2·R_lead.
 * This function extracts R_sensor given known R_lead.
 *
 * The error can be significant: 1 Ω lead → ~2.6°C error for Pt100.
 * This is why 3-wire and 4-wire are preferred for precision applications.
 */
double rtd_lead_resistance_2wire(double R_lead, double R_measured,
                                  double *R_sensor)
{
    double R_s = R_measured - 2.0 * R_lead;
    if (R_sensor) *R_sensor = R_s;
    return R_s;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L2+L4: Thermocouple EMF (Seebeck Effect)
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Initialize thermocouple model with ITS-90 coefficients.
 *
 * Implements NIST ITS-90 polynomial reference functions for each type.
 * Coefficients from NIST Monograph 175 (1993), ranges from IEC 60584-1.
 */
int tc_model_init(tc_model_t *tc, tc_type_t type)
{
    if (!tc) return -1;
    memset(tc, 0, sizeof(tc_model_t));
    tc->type = type;

    switch (type) {
    case TC_TYPE_K:  /* Ni-Cr / Ni-Al */
        /* Range: -270 to 0°C: 10th-order polynomial
         * Range: 0 to 1372°C: 10th-order + exponential */
        tc->t_min_C = -270.0;
        tc->t_max_C = 1372.0;
        tc->sensitivity_typical = 41.0;  /* µV/°C near 25°C */
        tc->tolerance_class1 = 1.5;       /* ±1.5°C or 0.004·|t| */
        tc->tolerance_class2 = 2.5;
        /* Coefficients for 0 to 1372°C range */
        tc->c[0] = 0.0;
        tc->c[1] = 39.450145139403;       /* µV/°C linear term ≈ sensitivity */
        tc->c[2] = 2.361201631e-2;
        tc->c[3] = -3.2876385e-4;
        tc->c[4] = -4.817756e-6;
        tc->c[5] = -1.3959e-8;
        tc->c[6] = 1.06981e-10;
        tc->c[7] = -5.53e-14;
        tc->c[8] = 0.0;
        tc->c[9] = 0.0;
        tc->a0_exp = 1.185976e-1;
        tc->a1_exp = -1.183432e-4;
        tc->a2_exp = 126.9686;
        tc->range_index = 0;
        break;

    case TC_TYPE_J:  /* Fe / Cu-Ni */
        tc->t_min_C = -210.0;
        tc->t_max_C = 1200.0;
        tc->sensitivity_typical = 52.0;  /* µV/°C, highest of base-metal types */
        tc->tolerance_class1 = 1.5;
        tc->tolerance_class2 = 2.5;
        tc->c[0] = 0.0;
        tc->c[1] = 50.381187921231;
        tc->c[2] = 3.0475837e-2;
        tc->c[3] = -8.568106e-5;
        tc->c[4] = 1.32282e-7;
        tc->c[5] = -1.7053e-10;
        tc->c[6] = 2.095e-13;
        tc->c[7] = -1.25e-16;
        tc->c[8] = 0.0;
        tc->range_index = 0;
        break;

    case TC_TYPE_T:  /* Cu / Cu-Ni */
        tc->t_min_C = -270.0;
        tc->t_max_C = 400.0;
        tc->sensitivity_typical = 43.0;  /* µV/°C */
        tc->tolerance_class1 = 0.5;
        tc->tolerance_class2 = 1.0;
        tc->c[0] = 0.0;
        tc->c[1] = 38.748106403642;
        tc->c[2] = 3.329222e-2;
        tc->c[3] = 2.06182e-4;
        tc->c[4] = -2.1882e-6;
        tc->c[5] = 1.100e-8;
        tc->c[6] = -3.08e-11;
        tc->c[7] = 4.6e-14;
        tc->c[8] = -2.9e-17;
        tc->range_index = 0;
        break;

    case TC_TYPE_E:  /* Ni-Cr / Cu-Ni, highest sensitivity */
        tc->t_min_C = -270.0;
        tc->t_max_C = 1000.0;
        tc->sensitivity_typical = 61.0;  /* µV/°C, highest common TC */
        tc->tolerance_class1 = 1.5;
        tc->tolerance_class2 = 2.5;
        tc->c[0] = 0.0;
        tc->c[1] = 58.665508702355;
        tc->c[2] = 4.503275e-2;
        tc->c[3] = 2.87156e-5;
        tc->c[4] = -3.3259e-7;
        tc->c[5] = 3.40e-10;
        tc->c[6] = -3.0e-13;
        tc->c[7] = 0.0;
        tc->range_index = 0;
        break;

    case TC_TYPE_S:  /* Pt-10%Rh / Pt */
        tc->t_min_C = -50.0;
        tc->t_max_C = 1768.0;
        tc->sensitivity_typical = 6.0;   /* µV/°C at 0°C, very low */
        tc->tolerance_class1 = 1.0;
        tc->tolerance_class2 = 1.5;
        tc->c[0] = 0.0;
        tc->c[1] = 5.40313308631;
        tc->c[2] = 1.2593429e-2;
        tc->c[3] = -2.32478e-5;
        tc->c[4] = 3.2207e-8;
        tc->c[5] = -3.315e-11;
        tc->c[6] = 2.56e-14;
        tc->c[7] = -1.2e-17;
        tc->a0_exp = 3.43e-2;
        tc->a1_exp = -9.84e-5;
        tc->a2_exp = 55.0;
        tc->range_index = 0;
        break;

    default:
        return -1;
    }
    return 0;
}

/**
 * @brief Compute thermocouple EMF for temperature difference.
 *
 * E(T_meas, T_ref) = f(T_meas) - f(T_ref)
 * where f(T) is the ITS-90 reference function (EMF relative to 0°C).
 *
 * f(T) = Σ_{i=0}^{n} c_i · T^i + a₀·exp[a₁·(T - a₂)²]
 */
static double tc_emf_relative_to_0C(const tc_model_t *tc, double T_C)
{
    if (!tc) return NAN;

    double E = 0.0;
    double T_pow = 1.0;
    for (int i = 0; i < 16; i++) {
        if (fabs(tc->c[i]) > 1e-30) {
            E += tc->c[i] * T_pow;
        }
        T_pow *= T_C;
    }

    /* Exponential correction (if applicable) */
    if (fabs(tc->a0_exp) > 1e-30) {
        double arg = tc->a1_exp * (T_C - tc->a2_exp);
        /* Prevent overflow: clamp argument */
        if (arg > 100.0) arg = 100.0;
        if (arg < -100.0) arg = -100.0;
        double exp_term = tc->a0_exp * exp(arg * arg);
        E += exp_term;
    }

    return E;
}

double tc_emf_from_T(const tc_model_t *tc, double T_meas_C, double T_ref_C)
{
    if (!tc) return NAN;
    double E_meas = tc_emf_relative_to_0C(tc, T_meas_C);
    double E_ref  = tc_emf_relative_to_0C(tc, T_ref_C);
    return E_meas - E_ref;
}

/**
 * @brief Convert thermocouple EMF back to temperature.
 *
 * Uses the measured voltage and cold junction temperature.
 * First computes the EMF relative to 0°C, then inverts the reference function.
 */
double tc_T_from_emf(const tc_model_t *tc, double emf_uV, double T_ref_C,
                      double tol)
{
    if (!tc || tol <= 0.0) return NAN;

    /* E(T, 0) = V_measured + E(T_ref, 0) */
    double E_ref = tc_emf_relative_to_0C(tc, T_ref_C);
    double E_target = emf_uV + E_ref;

    /* Newton-Raphson inversion of E(T, 0) = E_target */
    double T = T_ref_C + emf_uV / tc->sensitivity_typical;
    if (T < tc->t_min_C) T = tc->t_min_C + 1.0;
    if (T > tc->t_max_C) T = tc->t_max_C - 1.0;

    for (int iter = 0; iter < 50; iter++) {
        double E_val = tc_emf_relative_to_0C(tc, T);
        double residual = E_val - E_target;

        /* Numerical derivative */
        double dT = 0.01;
        double E_plus = tc_emf_relative_to_0C(tc, T + dT);
        double dE = (E_plus - E_val) / dT;

        if (fabs(dE) < 1e-30) break;
        double delta = residual / dE;
        T -= delta;
        if (fabs(delta) < tol) break;
    }
    return T;
}

double tc_sensitivity(const tc_model_t *tc, double T_C)
{
    if (!tc) return NAN;
    double dT = 0.01;
    double E1 = tc_emf_relative_to_0C(tc, T_C);
    double E2 = tc_emf_relative_to_0C(tc, T_C + dT);
    return (E2 - E1) / dT;
}

/**
 * @brief Apply cold junction compensation.
 *
 * The thermocouple measures V = E(T_meas, T_cj).
 * To find T_meas, we need T_cj (measured independently, e.g., by thermistor).
 *
 * V_compensated = V_measured + E(T_cj, 0°C)  →  gives E(T_meas, 0°C)
 */
double tc_cjc_compensate(const tc_model_t *tc, double V_measured,
                          double T_cj_measured, double T_cj_ref)
{
    if (!tc) return NAN;
    /* EMF corresponding to cold junction relative to 0°C */
    double E_cj = tc_emf_relative_to_0C(tc, T_cj_measured);
    double E_cj_ref = tc_emf_relative_to_0C(tc, T_cj_ref);

    /* Compensated EMF = measured + correction */
    return V_measured + E_cj - E_cj_ref;
}

/**
 * @brief Check Law of Intermediate Metals.
 *
 * If three metals A, B, C form junctions at T1, T2, T3:
 *   E_total = E_AB(T1) + E_BC(T2) + E_CA(T3)
 *
 * If T2 = T3, then E_total = E_AB(T1) regardless of metal C.
 * This allows using copper PCB traces and solder joints without
 * affecting the thermocouple reading.
 *
 * Returns: 0 if law is satisfied within tolerance, 1 if violation detected.
 */
int tc_law_intermediate_metals_check(double V_emf, double T1, double T2, double T3)
{
    /* If the two intermediate junctions are at the same temperature,
     * their EMF contributions cancel regardless of the third metal.
     * This is always valid for any material as long as T2 == T3.
     * T1 is the hot junction temperature (used for reference context). */
    (void)T1;  /* T1 is contextual: the EMF was generated at this temperature */
    double delta_T = fabs(T2 - T3);
    /* Allow 0.1°C tolerance */
    if (delta_T < 0.1) return 0;

    /* If T2 ≠ T3, a correction is needed.
     * The magnitude depends on the Seebeck coefficients of the metals.
     * For typical Cu-PbSn (solder) junctions: ~1-3 µV/°C.
     * Here we flag if the expected error exceeds 1% of the measured EMF. */
    double expected_error = delta_T * 3.0;  /* worst-case ~3 µV/°C */
    if (fabs(V_emf) > 1.0 && expected_error > 0.01 * fabs(V_emf)) {
        return 1;  /* significant violation */
    }
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L3: 1st-Order Dynamic Sensor Model
 * ═══════════════════════════════════════════════════════════════════════ */

int sensor_1st_order_init(sensor_1st_order_t *s1, double K, double tau,
                           double y0, double t0)
{
    if (!s1 || tau <= 0.0) return -1;
    s1->K = K;
    s1->tau = tau;
    s1->y0 = y0;
    s1->t0 = t0;
    return 0;
}

double sensor_1st_order_step_response(const sensor_1st_order_t *s1,
                                       double A_step, double t)
{
    if (!s1) return NAN;
    double dt = t - s1->t0;
    if (dt < 0.0) return s1->y0;
    return s1->K * A_step * (1.0 - exp(-dt / s1->tau)) + s1->y0;
}

double sensor_1st_order_freq_response_mag(const sensor_1st_order_t *s1,
                                            double f_hz)
{
    if (!s1) return NAN;
    double omega = 2.0 * M_PI * f_hz;
    return s1->K / sqrt(1.0 + omega * omega * s1->tau * s1->tau);
}

double sensor_1st_order_freq_response_phase(const sensor_1st_order_t *s1,
                                              double f_hz)
{
    if (!s1) return NAN;
    double omega = 2.0 * M_PI * f_hz;
    return -atan(omega * s1->tau);
}

double sensor_1st_order_bandwidth(const sensor_1st_order_t *s1)
{
    if (!s1) return NAN;
    return 1.0 / (2.0 * M_PI * s1->tau);
}

double sensor_1st_order_rise_time(const sensor_1st_order_t *s1)
{
    if (!s1) return NAN;
    /* 10%-90% rise time */
    return s1->tau * log(9.0);
}

double sensor_1st_order_settling_time(const sensor_1st_order_t *s1,
                                       double percent)
{
    if (!s1 || percent <= 0.0 || percent >= 100.0) return NAN;
    /* Time to settle within ±percent% of final value.
     * percent = 1 means ±1% → t_s = τ·ln(100) */
    return s1->tau * log(100.0 / percent);
}

int sensor_1st_order_simulate(const sensor_1st_order_t *s1,
                               const double *input, double *output,
                               size_t n, double dt)
{
    if (!s1 || !input || !output || n < 2 || dt <= 0.0) return -1;

    double alpha = dt / s1->tau;
    /* Exponential discretization: y[k] = exp(-α)·y[k-1] + (1-exp(-α))·K·x[k] */
    double decay = exp(-alpha);
    double gain_factor = s1->K * (1.0 - decay);

    output[0] = s1->y0;
    for (size_t k = 1; k < n; k++) {
        output[k] = decay * output[k-1] + gain_factor * input[k];
    }
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L3: 2nd-Order Dynamic Sensor Model
 * ═══════════════════════════════════════════════════════════════════════ */

int sensor_2nd_order_init(sensor_2nd_order_t *s2, double K,
                           double omega_n, double zeta,
                           double y0, double dy0)
{
    if (!s2 || omega_n <= 0.0 || zeta < 0.0) return -1;
    s2->K = K;
    s2->omega_n = omega_n;
    s2->zeta = zeta;
    s2->y0 = y0;
    s2->dy0 = dy0;
    return 0;
}

double sensor_2nd_order_step_response(const sensor_2nd_order_t *s2,
                                       double A_step, double t)
{
    if (!s2) return NAN;
    if (t < 0.0) return s2->y0;

    double wn = s2->omega_n;
    double z = s2->zeta;
    double steady_state = s2->K * A_step;

    if (z < 1.0) {
        /* Underdamped */
        double wd = wn * sqrt(1.0 - z * z);
        double phi = atan(sqrt(1.0 - z * z) / z);
        double exp_term = exp(-z * wn * t);
        double decay = exp_term / sqrt(1.0 - z * z);
        return steady_state * (1.0 - decay * sin(wd * t + phi));
    } else if (fabs(z - 1.0) < 1e-6) {
        /* Critically damped */
        return steady_state * (1.0 - (1.0 + wn * t) * exp(-wn * t));
    } else {
        /* Overdamped */
        double r1 = wn * (z + sqrt(z * z - 1.0));
        double r2 = wn * (z - sqrt(z * z - 1.0));
        return steady_state * (1.0 - (r2 * exp(-r1 * t) - r1 * exp(-r2 * t))
                                / (r2 - r1));
    }
}

double sensor_2nd_order_peak_time(const sensor_2nd_order_t *s2)
{
    if (!s2 || s2->zeta >= 1.0) return INFINITY;  /* no overshoot */
    double wd = s2->omega_n * sqrt(1.0 - s2->zeta * s2->zeta);
    return M_PI / wd;
}

double sensor_2nd_order_overshoot_percent(const sensor_2nd_order_t *s2)
{
    if (!s2 || s2->zeta >= 1.0) return 0.0;
    return 100.0 * exp(-M_PI * s2->zeta / sqrt(1.0 - s2->zeta * s2->zeta));
}

double sensor_2nd_order_settling_time(const sensor_2nd_order_t *s2,
                                       double band_percent)
{
    if (!s2) return NAN;
    /* t_s ≈ -ln(band/100) / (ζ·ω_n) */
    if (s2->zeta <= 0.0) return INFINITY;
    return -log(band_percent / 100.0) / (s2->zeta * s2->omega_n);
}

double sensor_2nd_order_freq_response_mag(const sensor_2nd_order_t *s2,
                                            double f_hz)
{
    if (!s2) return NAN;
    double w = 2.0 * M_PI * f_hz;
    double wn = s2->omega_n;
    double r = w / wn;
    double z = s2->zeta;
    double denom = sqrt((1.0 - r * r) * (1.0 - r * r)
                         + 4.0 * z * z * r * r);
    if (denom < 1e-30) return INFINITY;
    return s2->K / denom;
}

double sensor_2nd_order_freq_response_phase(const sensor_2nd_order_t *s2,
                                              double f_hz)
{
    if (!s2) return NAN;
    double r = 2.0 * M_PI * f_hz / s2->omega_n;
    double denom = 1.0 - r * r;
    if (fabs(denom) < 1e-30) return -M_PI / 2.0;
    return -atan2(2.0 * s2->zeta * r, denom);
}

double sensor_2nd_order_resonant_frequency(const sensor_2nd_order_t *s2)
{
    if (!s2) return NAN;
    /* Resonant frequency: ω_r = ω_n·√(1 - 2ζ²)  [valid for ζ < 1/√2] */
    double z2 = s2->zeta * s2->zeta;
    if (z2 >= 0.5) return 0.0;  /* no resonant peak */
    return s2->omega_n * sqrt(1.0 - 2.0 * z2) / (2.0 * M_PI);
}

double sensor_2nd_order_resonant_magnification(const sensor_2nd_order_t *s2)
{
    if (!s2) return NAN;
    /* Q = 1/(2ζ·√(1-ζ²))  [for ζ < 1/√2] */
    double z2 = s2->zeta * s2->zeta;
    if (z2 >= 0.5) return 1.0;
    return 1.0 / (2.0 * s2->zeta * sqrt(1.0 - z2));
}

double sensor_2nd_order_bandwidth(const sensor_2nd_order_t *s2)
{
    if (!s2) return NAN;
    /* -3 dB bandwidth for 2nd-order system:
     * ω_BW = ω_n · √[(1-2ζ²) + √(4ζ⁴ - 4ζ² + 2)] */
    double z2 = s2->zeta * s2->zeta;
    double inner = sqrt(4.0 * z2 * z2 - 4.0 * z2 + 2.0);
    double omega_bw = s2->omega_n * sqrt((1.0 - 2.0 * z2) + inner);
    return omega_bw / (2.0 * M_PI);
}

int sensor_2nd_order_simulate(const sensor_2nd_order_t *s2,
                               const double *input, double *output,
                               double *doutput, size_t n, double dt)
{
    if (!s2 || !input || !output || n < 2 || dt <= 0.0) return -1;

    /* State-space discretization using forward Euler:
     * State: [y, dy/dt]^T
     * d²y/dt² + 2ζω_n·dy/dt + ω_n²·y = K·ω_n²·x(t)
     *
     * Forward Euler:
     *   y[k+1] = y[k] + dt·dy[k]
     *   dy[k+1] = dy[k] + dt·(K·ω_n²·x[k] - 2ζω_n·dy[k] - ω_n²·y[k])
     */
    double y = s2->y0;
    double dy = s2->dy0;
    double wn2 = s2->omega_n * s2->omega_n;
    double damping = 2.0 * s2->zeta * s2->omega_n;

    output[0] = y;
    if (doutput) doutput[0] = dy;

    for (size_t k = 0; k < n - 1; k++) {
        double ddy = s2->K * wn2 * input[k] - damping * dy - wn2 * y;
        y = y + dt * dy;
        dy = dy + dt * ddy;
        output[k+1] = y;
        if (doutput) doutput[k+1] = dy;
    }
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L2: Lookup Table Transfer Function
 * ═══════════════════════════════════════════════════════════════════════ */

int sensor_lut_init(sensor_lut_t *lut, const double *x, const double *y, size_t n)
{
    if (!lut || !x || !y || n < 2) return -1;

    lut->n_points = n;
    lut->x_table = (double *)calloc(n, sizeof(double));
    lut->y_table = (double *)calloc(n, sizeof(double));
    if (!lut->x_table || !lut->y_table) {
        free(lut->x_table);
        free(lut->y_table);
        return -1;
    }
    for (size_t i = 0; i < n; i++) {
        lut->x_table[i] = x[i];
        lut->y_table[i] = y[i];
    }
    lut->x_min = x[0];
    lut->x_max = x[n-1];
    return 0;
}

void sensor_lut_free(sensor_lut_t *lut)
{
    if (lut) {
        free(lut->x_table);
        free(lut->y_table);
        lut->x_table = NULL;
        lut->y_table = NULL;
        lut->n_points = 0;
    }
}

int sensor_lut_get_interval(const sensor_lut_t *lut, double x,
                             size_t *idx_low, size_t *idx_high)
{
    if (!lut || !lut->x_table || lut->n_points < 2) return -1;

    /* Binary search for interval */
    if (x <= lut->x_table[0]) {
        *idx_low = 0;
        *idx_high = 1;
        return 0;
    }
    if (x >= lut->x_table[lut->n_points - 1]) {
        *idx_low = lut->n_points - 2;
        *idx_high = lut->n_points - 1;
        return 0;
    }

    size_t lo = 0, hi = lut->n_points - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (lut->x_table[mid] <= x) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    *idx_low = lo;
    *idx_high = hi;
    return 0;
}

double sensor_lut_evaluate(const sensor_lut_t *lut, double x)
{
    if (!lut) return NAN;
    size_t lo, hi;
    if (sensor_lut_get_interval(lut, x, &lo, &hi) < 0) return NAN;

    double x_lo = lut->x_table[lo];
    double x_hi = lut->x_table[hi];

    /* Linear interpolation */
    if (fabs(x_hi - x_lo) < 1e-30) return lut->y_table[lo];
    double alpha = (x - x_lo) / (x_hi - x_lo);
    return lut->y_table[lo] * (1.0 - alpha) + lut->y_table[hi] * alpha;
}

double sensor_lut_inverse(const sensor_lut_t *lut, double y)
{
    if (!lut) return NAN;
    if (y <= lut->y_table[0]) return lut->x_table[0];
    if (y >= lut->y_table[lut->n_points - 1]) return lut->x_table[lut->n_points - 1];

    /* Linear search for interval in y (assumes monotonic) */
    for (size_t i = 0; i < lut->n_points - 1; i++) {
        double y_lo = lut->y_table[i];
        double y_hi = lut->y_table[i+1];
        if ((y >= y_lo && y <= y_hi) || (y <= y_lo && y >= y_hi)) {
            if (fabs(y_hi - y_lo) < 1e-30)
                return lut->x_table[i];
            double alpha = (y - y_lo) / (y_hi - y_lo);
            return lut->x_table[i] * (1.0 - alpha) + lut->x_table[i+1] * alpha;
        }
    }
    return NAN;
}

double sensor_lut_sensitivity(const sensor_lut_t *lut, double x)
{
    if (!lut) return NAN;
    size_t lo, hi;
    if (sensor_lut_get_interval(lut, x, &lo, &hi) < 0) return NAN;

    double dx = lut->x_table[hi] - lut->x_table[lo];
    if (fabs(dx) < 1e-30) return 0.0;
    return (lut->y_table[hi] - lut->y_table[lo]) / dx;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L3: Linearity & Hysteresis Analysis
 * ═══════════════════════════════════════════════════════════════════════ */

int sensor_linearity_independent(const double *x, const double *y,
                                  size_t n, double *slope, double *intercept,
                                  double *max_error_percent)
{
    if (!x || !y || n < 3) return -1;

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_x2 += x[i] * x[i];
    }

    double denom = n * sum_x2 - sum_x * sum_x;
    if (fabs(denom) < 1e-30) return -1;

    double m = (n * sum_xy - sum_x * sum_y) / denom;
    double b = (sum_y - m * sum_x) / n;

    if (slope) *slope = m;
    if (intercept) *intercept = b;

    /* Find max deviation */
    double y_min = y[0], y_max = y[0];
    double max_dev = 0.0;
    for (size_t i = 0; i < n; i++) {
        double y_fit = m * x[i] + b;
        double dev = fabs(y[i] - y_fit);
        if (dev > max_dev) max_dev = dev;
        if (y[i] < y_min) y_min = y[i];
        if (y[i] > y_max) y_max = y[i];
    }

    double fs = y_max - y_min;
    if (max_error_percent) {
        *max_error_percent = (fs > 1e-30) ? (max_dev / fs * 100.0) : 0.0;
    }
    return 0;
}

int sensor_linearity_endpoint(const double *x, const double *y,
                               size_t n, double *max_error_percent)
{
    if (!x || !y || n < 2) return -1;

    double x0 = x[0], y0 = y[0];
    double xn = x[n-1], yn = y[n-1];
    double dx = xn - x0;

    if (fabs(dx) < 1e-30) return -1;
    double m = (yn - y0) / dx;
    double b = y0 - m * x0;

    double y_min = y[0], y_max = y[0];
    double max_dev = 0.0;
    for (size_t i = 0; i < n; i++) {
        double y_line = m * x[i] + b;
        double dev = fabs(y[i] - y_line);
        if (dev > max_dev) max_dev = dev;
        if (y[i] < y_min) y_min = y[i];
        if (y[i] > y_max) y_max = y[i];
    }

    double fs = y_max - y_min;
    if (max_error_percent) {
        *max_error_percent = (fs > 1e-30) ? (max_dev / fs * 100.0) : 0.0;
    }
    return 0;
}

int sensor_linearity_terminal(const double *x, const double *y,
                               size_t n, double *slope,
                               double *max_error_percent)
{
    if (!x || !y || n < 2) return -1;

    /* Best-fit line through origin: y = m·x
     * m = Σ(x_i·y_i) / Σ(x_i²) */
    double sum_xy = 0.0, sum_x2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum_xy += x[i] * y[i];
        sum_x2 += x[i] * x[i];
    }
    if (sum_x2 < 1e-30) return -1;
    double m = sum_xy / sum_x2;

    if (slope) *slope = m;

    double y_min = y[0], y_max = y[0];
    double max_dev = 0.0;
    for (size_t i = 0; i < n; i++) {
        double dev = fabs(y[i] - m * x[i]);
        if (dev > max_dev) max_dev = dev;
        if (y[i] < y_min) y_min = y[i];
        if (y[i] > y_max) y_max = y[i];
    }

    double fs = y_max - y_min;
    if (max_error_percent) {
        *max_error_percent = (fs > 1e-30) ? (max_dev / fs * 100.0) : 0.0;
    }
    return 0;
}

double sensor_hysteresis_percent(const double *x, const double *y_up,
                                  const double *y_down, size_t n,
                                  double full_scale)
{
    if (!x || !y_up || !y_down || n < 2 || full_scale <= 0.0) return 0.0;

    double max_hyst = 0.0;
    for (size_t i = 0; i < n; i++) {
        double diff = fabs(y_up[i] - y_down[i]);
        if (diff > max_hyst) max_hyst = diff;
    }
    return max_hyst / full_scale * 100.0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L3: Dynamic Error Estimation
 * ═══════════════════════════════════════════════════════════════════════ */

double sensor_dynamic_ramp_error_1st_order(double ramp_rate, double tau)
{
    return ramp_rate * tau;
}

double sensor_dynamic_sinusoidal_mag_error_1st_order(double freq_hz, double tau)
{
    double omega = 2.0 * M_PI * freq_hz;
    double mag = 1.0 / sqrt(1.0 + omega * omega * tau * tau);
    return fabs(1.0 - mag) * 100.0;
}

double sensor_dynamic_phase_lag_1st_order(double freq_hz, double tau)
{
    double omega = 2.0 * M_PI * freq_hz;
    return -atan(omega * tau) * 180.0 / M_PI;
}

double sensor_dynamic_ramp_error_2nd_order(double ramp_rate,
                                            double omega_n, double zeta)
{
    if (omega_n <= 0.0) return INFINITY;
    return 2.0 * zeta * ramp_rate / omega_n;
}
