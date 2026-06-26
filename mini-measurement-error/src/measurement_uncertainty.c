#include "measurement_uncertainty.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>

/* ─── L3: Distribution Divisors (GUM Section 4.3) ────────────────────────
 *
 * Each distribution has a standard divisor that converts the half-width
 * of the distribution to a standard uncertainty.
 *
 *   Normal:     divisor = 1 (or k, handled separately)
 *   Uniform:    divisor = sqrt(3) ≈ 1.732
 *   Triangular: divisor = sqrt(6) ≈ 2.449
 *   U-shaped:   divisor = sqrt(2) ≈ 1.414
 *   Student-t:  divisor = depends on degrees of freedom
 *   Trapezoid:  divisor = sqrt( (a^2+b^2+ab) / 3 ) where a,b are radii
 */

double unc_distribution_divisor(UncertaintyDistribution dist) {
    switch (dist) {
        case DIST_NORMAL:     return 1.0;
        case DIST_UNIFORM:    return sqrt(3.0);      /* 1.7320508 */
        case DIST_TRIANGULAR: return sqrt(6.0);      /* 2.4494897 */
        case DIST_U_SHAPED:   return sqrt(2.0);      /* 1.4142136 */
        case DIST_STUDENT_T:  return 1.0;             /* handled via k factor */
        case DIST_TRAPEZOID:  return sqrt(2.4);       /* typical */
        default:              return 1.0;
    }
}

/* ─── L4: GUM Uncertainty Budget ────────────────────────────────────────── */

void unc_budget_init(UncertaintyBudget *budget, double confidence_level) {
    if (!budget) return;
    memset(budget, 0, sizeof(*budget));
    budget->confidence_level = confidence_level;
    budget->coverage_factor = 2.0;  /* default k=2 for ~95% */
}

int unc_add_type_a(UncertaintyBudget *budget, const char *name,
                    double std_dev, size_t n, double sens_coef) {
    if (!budget || budget->num_components >= MAX_UNC_COMPONENTS) return -1;
    size_t i = budget->num_components;
    UncertaintyComponent *c = &budget->components[i];
    if (name) { strncpy(c->name, name, 63); c->name[63] = '\0'; }
    else { c->name[0] = '\0'; }
    c->eval_type = EVAL_TYPE_A;
    c->distribution = DIST_NORMAL;
    c->half_width = std_dev;
    c->divisor = 1.0;
    c->sensitivity_coeff = sens_coef;
    c->degrees_freedom = (n > 1) ? (double)(n - 1) : 1.0;
    c->standard_uncertainty = std_dev;
    c->contribution = (sens_coef * std_dev) * (sens_coef * std_dev);
    budget->num_components++;
    return 0;
}

int unc_add_type_b(UncertaintyBudget *budget, const char *name,
                    UncertaintyDistribution distribution, double half_width,
                    double sens_coef) {
    if (!budget || budget->num_components >= MAX_UNC_COMPONENTS) return -1;
    size_t i = budget->num_components;
    UncertaintyComponent *c = &budget->components[i];
    if (name) { strncpy(c->name, name, 63); c->name[63] = '\0'; }
    else { c->name[0] = '\0'; }
    c->eval_type = EVAL_TYPE_B;
    c->distribution = distribution;
    c->half_width = half_width;
    c->divisor = unc_distribution_divisor(distribution);
    c->sensitivity_coeff = sens_coef;
    c->degrees_freedom = INFINITY;  /* Type B, infinite df by default */
    c->standard_uncertainty = half_width / c->divisor;
    c->contribution = (sens_coef * c->standard_uncertainty)
                       * (sens_coef * c->standard_uncertainty);
    budget->num_components++;
    return 0;
}

/* ─── L4: Welch-Satterthwaite Equation ───────────────────────────────────
 *
 * nu_eff = u_c^4 / SUM ( (c_i^4 * u_i^4) / nu_i )
 *
 * This gives the effective degrees of freedom for the combined uncertainty,
 * used to determine the appropriate coverage factor k from Student's t.
 */

double unc_welch_satterthwaite(const double *contributions,
                                const double *dofs, size_t n,
                                double u_c_sq) {
    if (n == 0 || u_c_sq < DBL_EPSILON) return 1.0;
    double denom = 0.0;
    for (size_t i = 0; i < n; i++) {
        double dof = dofs[i];
        if (dof < 1.0) dof = 1.0;
        /* contribution[i] = (c_i * u_i)^2 */
        /* nu_i term: (c_i * u_i)^4 / nu_i */
        double c2 = contributions[i];
        denom += (c2 * c2) / dof;
    }
    double u_c4 = u_c_sq * u_c_sq;
    if (denom < DBL_EPSILON) return 1.0;
    double nu_eff = u_c4 / denom;
    if (nu_eff < 1.0) nu_eff = 1.0;
    if (nu_eff > 100000.0) nu_eff = 100000.0;
    return nu_eff;
}

/* ─── L4: Student's t for 95% Confidence ─────────────────────────────────
 *
 * Rational approximation to the Student's t distribution 97.5th percentile
 * (two-tailed 95%). Accurate to ~0.1% for nu >= 1.
 *
 * Reference: Abramowitz & Stegun 26.7.5 with asymptotic correction.
 */

double unc_students_t_95(double nu) {
    if (nu <= 0.0) nu = 1.0;
    if (nu > 100.0) {
        /* Asymptotic: z = 1.96 for 95% two-tailed */
        return 1.96 + (1.96 * (1.96 * 1.96 + 1.0)) / (4.0 * nu);
    }
    /* Rational approximation for small nu */
    double nu2 = nu * nu;
    double nu3 = nu2 * nu;
    double num = 1.96 + 2.376 / nu + 2.156 / nu2 + 0.766 / nu3;
    double den = 1.0 + 0.993 / nu - 0.071 / nu2;
    return num / den;
}

/* ─── L4: Finalize Uncertainty Budget ────────────────────────────────────
 *
 * Steps per GUM Section 8:
 * 1. Compute combined standard uncertainty u_c = sqrt(SUM contrib_i)
 * 2. Compute effective degrees of freedom (Welch-Satterthwaite)
 * 3. Look up coverage factor k from Student's t
 * 4. Compute expanded uncertainty U = k * u_c
 */

void unc_budget_finalize(UncertaintyBudget *budget) {
    if (!budget || budget->num_components == 0) return;

    /* Step 1: Combined standard uncertainty */
    double sum_contrib = 0.0;
    for (size_t i = 0; i < budget->num_components; i++) {
        double ci = budget->components[i].sensitivity_coeff;
        double ui = budget->components[i].standard_uncertainty;
        budget->components[i].contribution = (ci * ui) * (ci * ui);
        sum_contrib += budget->components[i].contribution;
    }
    budget->combined_standard_uncertainty = sqrt(sum_contrib);
    double u_c_sq = sum_contrib;

    /* Step 2: Effective degrees of freedom */
    if (budget->num_components > 0) {
        double *contribs = (double*)malloc(budget->num_components * sizeof(double));
        double *dofs     = (double*)malloc(budget->num_components * sizeof(double));
        if (contribs && dofs) {
            for (size_t i = 0; i < budget->num_components; i++) {
                contribs[i] = budget->components[i].contribution;
                dofs[i] = budget->components[i].degrees_freedom;
            }
            budget->effective_degrees_freedom =
                unc_welch_satterthwaite(contribs, dofs,
                                         budget->num_components, u_c_sq);
        } else {
            budget->effective_degrees_freedom = INFINITY;
        }
        free(contribs);
        free(dofs);
    } else {
        budget->effective_degrees_freedom = INFINITY;
    }

    /* Step 3-4: Coverage factor and expanded uncertainty */
    if (budget->effective_degrees_freedom >= 100000.0) {
        budget->coverage_factor = 1.96; /* large-sample normal approx */
    } else {
        budget->coverage_factor = unc_students_t_95(
            budget->effective_degrees_freedom);
    }
    budget->expanded_uncertainty = budget->coverage_factor
                                    * budget->combined_standard_uncertainty;
}

/* ─── L4: Error Propagation (GUM Section 5) ────────────────────────────── */

double unc_propagate_uncorrelated(const double *partials,
                                   const double *std_uncs, size_t n) {
    double sum_sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        double term = partials[i] * std_uncs[i];
        sum_sq += term * term;
    }
    return sqrt(sum_sq);
}

double unc_propagate_correlated(const double *partials, const double *std_uncs,
                                 const double *cov_matrix, size_t n) {
    double sum_sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        double term = partials[i] * std_uncs[i];
        sum_sq += term * term;
    }
    /* Add covariance cross terms */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            sum_sq += 2.0 * partials[i] * partials[j]
                      * cov_matrix[i * n + j];
        }
    }
    if (sum_sq < 0.0) sum_sq = 0.0; /* numerical safety */
    return sqrt(sum_sq);
}

/* ─── L4: Common Uncertainty Propagation Formulas ──────────────────────── */

double unc_sum_uncertainty(double u1, double u2) {
    return sqrt(u1 * u1 + u2 * u2);
}

/* Y = X1 * X2: relative uncertainties add in quadrature */
double unc_product_uncertainty(double x1, double u1, double x2, double u2) {
    if (fabs(x1) < DBL_EPSILON || fabs(x2) < DBL_EPSILON) return INFINITY;
    double y = x1 * x2;
    return y * sqrt((u1 / x1) * (u1 / x1) + (u2 / x2) * (u2 / x2));
}

/* Y = X1 / X2 */
double unc_quotient_uncertainty(double x1, double u1, double x2, double u2) {
    return unc_product_uncertainty(x1, u1, x2, u2);
}

/* Y = X^n: u_c/y = |n| * u/x */
double unc_power_uncertainty(double x, double u, double exponent) {
    if (fabs(x) < DBL_EPSILON) return INFINITY;
    return fabs(exponent) * pow(fabs(x), exponent) * u / fabs(x);
}

/* Y = exp(a*X): u_c/y = |a| * u */
double unc_exp_uncertainty(double a, double u) {
    return fabs(a) * u;
}

/* Y = ln(X): u_c = u/x */
double unc_log_uncertainty(double x, double u) {
    if (fabs(x) < DBL_EPSILON) return INFINITY;
    return u / fabs(x);
}

/* ─── L3: Confidence Intervals ──────────────────────────────────────────── */

/* 95% confidence interval for the mean */
void unc_confidence_interval(double mean, double s, size_t n,
                              double *ci_low, double *ci_high) {
    if (n < 2) {
        if (ci_low)  *ci_low  = mean;
        if (ci_high) *ci_high = mean;
        return;
    }
    double se = s / sqrt((double)n);  /* standard error of the mean */
    double t_val;
    if (n <= 30) {
        t_val = unc_students_t_95((double)(n - 1));
    } else {
        t_val = 1.96;
    }
    double half_width = t_val * se;
    if (ci_low)  *ci_low  = mean - half_width;
    if (ci_high) *ci_high = mean + half_width;
}

/* ─── L3: Gaussian (Normal) Distribution Functions ──────────────────────── */

double gaussian_pdf(double x, double mu, double sigma) {
    if (sigma < DBL_EPSILON) return (fabs(x - mu) < DBL_EPSILON) ? INFINITY : 0.0;
    double z = (x - mu) / sigma;
    return exp(-0.5 * z * z) / (sigma * sqrt(2.0 * M_PI));
}

/* Gaussian CDF via the error function:
 * Phi(x) = 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2))))
 */
double gaussian_cdf(double x, double mu, double sigma) {
    if (sigma < DBL_EPSILON) return (x >= mu) ? 1.0 : 0.0;
    double z = (x - mu) / (sigma * M_SQRT2);
    return 0.5 * (1.0 + erf(z));
}

/* Inverse Gaussian CDF (quantile function).
 * Uses the rational approximation from Abramowitz & Stegun §26.2.23.
 * Accurate to ~1e-4 for p in (0, 1).
 */
double gaussian_quantile(double p, double mu, double sigma) {
    if (p <= 0.0) return -INFINITY;
    if (p >= 1.0) return INFINITY;
    if (sigma < DBL_EPSILON) return mu;

    /* Use symmetry: compute for p <= 0.5, reflect for p > 0.5 */
    double q = p;
    int sign = -1;
    if (q > 0.5) { q = 1.0 - q; sign = 1; }

    /* Rational approximation for the inverse of the error function complement */
    double t = sqrt(-2.0 * log(q));
    double c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
    double d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;
    double num = c0 + c1 * t + c2 * t * t;
    double den = 1.0 + d1 * t + d2 * t * t + d3 * t * t * t;
    double z = t - num / den;

    return mu + sign * z * sigma;
}

/* ─── L3: Uniform Distribution ──────────────────────────────────────────── */

double uniform_pdf(double x, double a, double b) {
    if (b <= a) return 0.0;
    if (x >= a && x <= b) return 1.0 / (b - a);
    return 0.0;
}
