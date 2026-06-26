#include "measurement_uncertainty.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define EPS 1e-10
#define NEAR(a,b) (fabs((a)-(b)) < EPS)

static int tests_run = 0, tests_passed = 0;
#define TEST(n) do { tests_run++; printf("  %s... ", n); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)

int main(void) {
    printf("=== Uncertainty Analysis Tests ===\n");

    /* Test: Distribution divisors */
    TEST("distribution_divisors");
    assert(NEAR(unc_distribution_divisor(DIST_UNIFORM), sqrt(3.0)));
    assert(NEAR(unc_distribution_divisor(DIST_TRIANGULAR), sqrt(6.0)));
    assert(NEAR(unc_distribution_divisor(DIST_U_SHAPED), sqrt(2.0)));
    PASS();

    /* Test: Uncertainty budget (GUM) */
    TEST("unc_budget");
    UncertaintyBudget budget;
    unc_budget_init(&budget, 0.95);
    unc_add_type_a(&budget, "Repeatability", 0.02, 10, 1.0);
    unc_add_type_b(&budget, "Calibration", DIST_NORMAL, 0.05, 1.0);
    unc_add_type_b(&budget, "Resolution", DIST_UNIFORM, 0.01, 1.0);
    unc_budget_finalize(&budget);
    assert(budget.combined_standard_uncertainty > 0.0);
    assert(budget.expanded_uncertainty
           > budget.combined_standard_uncertainty);
    PASS();

    /* Test: Propagation uncorrelated */
    TEST("propagate_uncorrelated");
    double partials[] = {1.0, 1.0};
    double std_uncs[] = {0.1, 0.2};
    double uc = unc_propagate_uncorrelated(partials, std_uncs, 2);
    assert(NEAR(uc, sqrt(0.01 + 0.04)));
    PASS();

    /* Test: Common propagation formulas */
    TEST("sum_uncertainty");
    assert(NEAR(unc_sum_uncertainty(3.0, 4.0), 5.0));
    PASS();

    TEST("product_uncertainty");
    double up = unc_product_uncertainty(10.0, 0.1, 20.0, 0.2);
    assert(up > 0.0);
    PASS();

    /* Test: Confidence interval */
    TEST("confidence_interval");
    double ci_low, ci_high;
    unc_confidence_interval(50.0, 1.0, 30, &ci_low, &ci_high);
    assert(ci_low < 50.0 && ci_high > 50.0);
    PASS();

    /* Test: Welch-Satterthwaite */
    TEST("welch_satterthwaite");
    double contribs[] = {0.01, 0.04};
    double dofs[] = {9.0, 100.0};
    double nu_eff = unc_welch_satterthwaite(contribs, dofs, 2, 0.05);
    assert(nu_eff > 0.0);
    PASS();

    /* Test: Student's t */
    TEST("students_t_95");
    double t_inf = unc_students_t_95(1e6);
    assert(fabs(t_inf - 1.96) < 0.01);
    double t_10 = unc_students_t_95(10.0);
    assert(t_10 > 1.96);
    PASS();

    /* Test: Gaussian PDF */
    TEST("gaussian_pdf");
    double pdf0 = gaussian_pdf(0.0, 0.0, 1.0);
    assert(NEAR(pdf0, 1.0 / sqrt(2.0 * M_PI)));
    double pdf_sym = gaussian_pdf(1.0, 0.0, 1.0);
    assert(NEAR(pdf_sym, gaussian_pdf(-1.0, 0.0, 1.0)));
    PASS();

    /* Test: Gaussian CDF */
    TEST("gaussian_cdf");
    double cdf0 = gaussian_cdf(0.0, 0.0, 1.0);
    assert(NEAR(cdf0, 0.5));
    double cdf_hi = gaussian_cdf(100.0, 0.0, 1.0);
    assert(NEAR(cdf_hi, 1.0));
    PASS();

    /* Test: Gaussian quantile */
    TEST("gaussian_quantile");
    double q50 = gaussian_quantile(0.5, 0.0, 1.0);
    assert(fabs(q50) < 0.001);
    double q975 = gaussian_quantile(0.975, 0.0, 1.0);
    assert(fabs(q975 - 1.96) < 0.05);
    PASS();

    /* Test: Uniform PDF */
    TEST("uniform_pdf");
    assert(NEAR(uniform_pdf(0.5, 0.0, 1.0), 1.0));
    assert(NEAR(uniform_pdf(1.5, 0.0, 1.0), 0.0));
    PASS();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
