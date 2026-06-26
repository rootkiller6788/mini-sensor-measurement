
#ifndef MEASUREMENT_UNCERTAINTY_H
#define MEASUREMENT_UNCERTAINTY_H

#include <stddef.h>
#include "measurement_error.h"

/* ─── L3: Distribution Enumeration ────────────────────────────────────────── */
typedef enum {
    DIST_NORMAL     = 0,
    DIST_UNIFORM    = 1,
    DIST_TRIANGULAR = 2,
    DIST_U_SHAPED   = 3,
    DIST_STUDENT_T  = 4,
    DIST_TRAPEZOID  = 5
} UncertaintyDistribution;

double unc_distribution_divisor(UncertaintyDistribution dist);

/* ─── L3: Uncertainty Component ───────────────────────────────────────────── */
typedef enum {
    EVAL_TYPE_A = 0,
    EVAL_TYPE_B = 1
} UncertaintyEvalType;

typedef struct {
    char                    name[64];
    UncertaintyEvalType     eval_type;
    UncertaintyDistribution distribution;
    double                  half_width;
    double                  divisor;
    double                  sensitivity_coeff;
    double                  degrees_freedom;
    double                  standard_uncertainty;
    double                  contribution;
} UncertaintyComponent;

/* ─── L3: Uncertainty Budget (GUM section 8) ──────────────────────────────── */
#define MAX_UNC_COMPONENTS 32

typedef struct {
    UncertaintyComponent components[MAX_UNC_COMPONENTS];
    size_t               num_components;
    double               combined_standard_uncertainty;
    double               effective_degrees_freedom;
    double               coverage_factor;
    double               expanded_uncertainty;
    double               confidence_level;
} UncertaintyBudget;

/* ─── L4: Core GUM Functions ───────────────────────────────────────────────── */
void unc_budget_init(UncertaintyBudget *budget, double confidence_level);
int  unc_add_type_a(UncertaintyBudget *budget, const char *name,
                    double std_dev, size_t n, double sens_coef);
int  unc_add_type_b(UncertaintyBudget *budget, const char *name,
                    UncertaintyDistribution distribution, double half_width,
                    double sens_coef);
void unc_budget_finalize(UncertaintyBudget *budget);

/* ─── L4: Error Propagation (GUM section 5) ───────────────────────────────── */

/** First-order Taylor propagation, uncorrelated inputs */
double unc_propagate_uncorrelated(const double *partials,
                                   const double *std_uncs, size_t n);

/** First-order Taylor propagation, correlated inputs (includes covariance) */
double unc_propagate_correlated(const double *partials, const double *std_uncs,
                                 const double *cov_matrix, size_t n);

/* ─── L4: Common Uncertainty Propagation Formulas ─────────────────────────── */
double unc_sum_uncertainty(double u1, double u2);
double unc_product_uncertainty(double x1, double u1, double x2, double u2);
double unc_quotient_uncertainty(double x1, double u1, double x2, double u2);
double unc_power_uncertainty(double x, double u, double exponent);
double unc_exp_uncertainty(double a, double u);
double unc_log_uncertainty(double x, double u);

/* ─── L3: Confidence Intervals ────────────────────────────────────────────── */
void unc_confidence_interval(double mean, double s, size_t n,
                              double *ci_low, double *ci_high);

/* ─── L4: Welch-Satterthwaite equation ────────────────────────────────────── */
double unc_welch_satterthwaite(const double *contributions,
                                const double *dofs, size_t n,
                                double u_c_sq);

/* ─── L4: Student's t lookup (95% two-tailed approximation) ───────────────── */
double unc_students_t_95(double nu);

/* ─── L3: Gaussian Distribution Functions ─────────────────────────────────── */
double gaussian_pdf(double x, double mu, double sigma);
double gaussian_cdf(double x, double mu, double sigma);
double gaussian_quantile(double p, double mu, double sigma);
double uniform_pdf(double x, double a, double b);

#endif /* MEASUREMENT_UNCERTAINTY_H */
