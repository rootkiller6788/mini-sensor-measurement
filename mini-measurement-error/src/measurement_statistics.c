#include "measurement_statistics.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>

/* ─── L3: MLE Gaussian Fit ────────────────────────────────────────────────
 *
 * Maximum Likelihood Estimation for Normal distribution:
 *   mu_hat = (1/n) * SUM x_i
 *   sigma_unbiased = sqrt( (1/(n-1)) * SUM (x_i - mu)^2 )
 *
 * These are the standard sample mean and unbiased sample std dev.
 */

void stats_fit_gaussian(const double *data, size_t n,
                         double *mu_hat, double *sigma_hat) {
    if (!data || n < 2 || !mu_hat || !sigma_hat) {
        if (mu_hat) *mu_hat = 0.0;
        if (sigma_hat) *sigma_hat = 0.0;
        return;
    }
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += data[i];
    *mu_hat = sum / (double)n;

    double sum_sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        double d = data[i] - *mu_hat;
        sum_sq += d * d;
    }
    *sigma_hat = sqrt(sum_sq / (double)(n - 1));
}

/* ─── L3: Kolmogorov-Smirnov Test ─────────────────────────────────────────
 *
 * KS statistic: D_n = sup |F_n(x) - F(x)|
 * where F_n is the empirical CDF and F is the theoretical (normal) CDF.
 *
 * The empirical CDF at point x: F_n(x) = (number of samples <= x) / n
 */

/* Forward declaration for Gaussian CDF (implemented in measurement_uncertainty.c) */
static double gs_cdf(double x, double mu, double sigma) {
    if (sigma < DBL_EPSILON) return (x >= mu) ? 1.0 : 0.0;
    double z = (x - mu) / (sigma * M_SQRT2);
    return 0.5 * (1.0 + erf(z));
}

/* Comparator for qsort */
static int cmp_double(const void *a, const void *b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

double stats_ks_test_normal(const double *data, size_t n,
                             double mu, double sigma) {
    if (!data || n < 3 || sigma < DBL_EPSILON) return 0.0;

    /* Sort data copy */
    double *sorted = (double*)malloc(n * sizeof(double));
    if (!sorted) return 0.0;
    memcpy(sorted, data, n * sizeof(double));
    qsort(sorted, n, sizeof(double), cmp_double);

    double d_max = 0.0;
    for (size_t i = 0; i < n; i++) {
        double ecdf = (double)(i + 1) / (double)n;
        double tcdf = gs_cdf(sorted[i], mu, sigma);
        double d_plus = fabs(ecdf - tcdf);
        if (d_plus > d_max) d_max = d_plus;

        /* Also check the step just before each observation */
        double ecdf_prev = (double)i / (double)n;
        double d_minus = fabs(ecdf_prev - tcdf);
        if (d_minus > d_max) d_max = d_minus;
    }
    free(sorted);
    return d_max;
}

/* ─── L3: Shapiro-Wilk W Statistic (Simplified) ───────────────────────────
 *
 * Tests the null hypothesis that the sample comes from a normal distribution.
 * W = (SUM a_i * x_(i))^2 / SUM (x_i - xbar)^2
 * where a_i are coefficients and x_(i) are order statistics.
 *
 * This is a simplified approximation valid for 3 <= n <= 50.
 * Reference: Shapiro & Wilk (1965), Biometrika 52, 591-611.
 */

double stats_shapiro_wilk(const double *data, size_t n) {
    if (!data || n < 3 || n > 50) return 0.0;

    double *sorted = (double*)malloc(n * sizeof(double));
    if (!sorted) return 0.0;
    memcpy(sorted, data, n * sizeof(double));
    qsort(sorted, n, sizeof(double), cmp_double);

    /* Compute mean and SS */
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += sorted[i];
    double mean = sum / (double)n;
    double ss = 0.0;
    for (size_t i = 0; i < n; i++) {
        double d = sorted[i] - mean;
        ss += d * d;
    }
    if (ss < DBL_EPSILON) { free(sorted); return 0.0; }

    /* Blom's plotting positions for approximate Shapiro-Wilk coefficients.
     * The exact coefficients require solving for expected normal order stats.
     * We use the simplified Blom score: a_i ~= inv_norm((i - 3/8) / (n + 1/4))
     */
    double b = 0.0;
    double a2_sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        double p = ((double)(i + 1) - 0.375) / ((double)n + 0.25);
        /* Approximate inverse normal via rational approximation */
        double q = (p <= 0.5) ? p : 1.0 - p;
        int sign = (p > 0.5) ? 1 : -1;
        double t_val = sqrt(-2.0 * log(q > 1e-15 ? q : 1e-15));
        double c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
        double d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;
        double num = c0 + c1 * t_val + c2 * t_val * t_val;
        double den = 1.0 + d1 * t_val + d2 * t_val * t_val
                     + d3 * t_val * t_val * t_val;
        double mi = sign * (t_val - num / den);
        b += mi * sorted[i];
        a2_sum += mi * mi;
    }
    free(sorted);

    double w = (b * b) / ss;
    return w;
}

/* ─── L3: One-Sample t-Test ──────────────────────────────────────────────
 *
 * H0: population mean = hypothesized_mu
 * t = (xbar - mu0) / (s / sqrt(n))
 * Under H0, t follows Student's t with n-1 degrees of freedom.
 */

double stats_one_sample_t_test(const double *data, size_t n,
                                double hypothesized_mu) {
    if (!data || n < 2) return 0.0;
    double sum = 0.0, sum_sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += data[i];
        sum_sq += data[i] * data[i];
    }
    double xbar = sum / (double)n;
    double s2 = (sum_sq - sum * sum / (double)n) / (double)(n - 1);
    if (s2 < DBL_EPSILON) return 0.0;
    double s = sqrt(s2);
    return (xbar - hypothesized_mu) / (s / sqrt((double)n));
}

/* ─── L3: Two-Sample t-Test (Welch's) ───────────────────────────────────
 *
 * H0: mu1 = mu2 (no difference between two independent samples)
 * Welch's t-test does NOT assume equal variances.
 *
 * t = (xbar1 - xbar2) / sqrt(s1^2/n1 + s2^2/n2)
 */

double stats_two_sample_t_test(const double *data1, size_t n1,
                                const double *data2, size_t n2) {
    if (!data1 || !data2 || n1 < 2 || n2 < 2) return 0.0;

    double sum1 = 0.0, sum_sq1 = 0.0;
    for (size_t i = 0; i < n1; i++) {
        sum1 += data1[i];
        sum_sq1 += data1[i] * data1[i];
    }
    double xbar1 = sum1 / (double)n1;
    double s1_2 = (sum_sq1 - sum1 * sum1 / (double)n1) / (double)(n1 - 1);

    double sum2 = 0.0, sum_sq2 = 0.0;
    for (size_t i = 0; i < n2; i++) {
        sum2 += data2[i];
        sum_sq2 += data2[i] * data2[i];
    }
    double xbar2 = sum2 / (double)n2;
    double s2_2 = (sum_sq2 - sum2 * sum2 / (double)n2) / (double)(n2 - 1);

    double se = sqrt(s1_2 / (double)n1 + s2_2 / (double)n2);
    if (se < DBL_EPSILON) return 0.0;
    return (xbar1 - xbar2) / se;
}

/* ─── L3: F-Test for Equality of Variances ───────────────────────────────
 *
 * H0: sigma1^2 = sigma2^2
 * F = s_larger^2 / s_smaller^2
 * Under H0, F ~ F(n_larger-1, n_smaller-1)
 */

double stats_f_test_variances(const double *data1, size_t n1,
                               const double *data2, size_t n2) {
    if (!data1 || !data2 || n1 < 2 || n2 < 2) return 0.0;

    double sum1 = 0.0, sum_sq1 = 0.0;
    for (size_t i = 0; i < n1; i++) {
        sum1 += data1[i];
        sum_sq1 += data1[i] * data1[i];
    }
    double var1 = (sum_sq1 - sum1 * sum1 / (double)n1) / (double)(n1 - 1);

    double sum2 = 0.0, sum_sq2 = 0.0;
    for (size_t i = 0; i < n2; i++) {
        sum2 += data2[i];
        sum_sq2 += data2[i] * data2[i];
    }
    double var2 = (sum_sq2 - sum2 * sum2 / (double)n2) / (double)(n2 - 1);

    if (var1 < DBL_EPSILON && var2 < DBL_EPSILON) return 1.0;
    if (var2 < DBL_EPSILON) return INFINITY;
    if (var1 < DBL_EPSILON) return 0.0;

    /* Put larger variance in numerator */
    if (var1 >= var2) return var1 / var2;
    else return var2 / var1;
}

/* ─── L3: Chi-Squared Variance Test ─────────────────────────────────────
 *
 * H0: sigma^2 = sigma0^2
 * chi^2 = (n - 1) * s^2 / sigma0^2
 * Under H0, chi^2 ~ ChiSq(n-1)
 */

double stats_chi_sq_variance_test(const double *data, size_t n,
                                   double sigma0_sq) {
    if (!data || n < 2 || sigma0_sq < DBL_EPSILON) return 0.0;

    double sum = 0.0, sum_sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += data[i];
        sum_sq += data[i] * data[i];
    }
    double s2 = (sum_sq - sum * sum / (double)n) / (double)(n - 1);
    return (double)(n - 1) * s2 / sigma0_sq;
}

/* ─── L3: Prediction Interval ─────────────────────────────────────────────
 *
 * PI = xbar +/- t_{alpha/2, n-1} * s * sqrt(1 + 1/n)
 *
 * This is wider than the confidence interval for the mean because
 * it accounts for the variability of a single new observation.
 */

void stats_prediction_interval(const double *data, size_t n, double alpha,
                                double *pi_low, double *pi_high) {
    if (!data || n < 2) {
        if (pi_low) *pi_low = 0.0;
        if (pi_high) *pi_high = 0.0;
        return;
    }

    double sum = 0.0, sum_sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += data[i];
        sum_sq += data[i] * data[i];
    }
    double xbar = sum / (double)n;
    double s = sqrt((sum_sq - sum * sum / (double)n) / (double)(n - 1));

    /* Approximate t-value for given alpha */
    double t_val;
    if (alpha <= 0.01) t_val = 2.58;
    else if (alpha <= 0.05) t_val = 1.96;
    else if (alpha <= 0.10) t_val = 1.645;
    else t_val = 1.0;

    if (n <= 30) {
        t_val = t_val + 1.0 / sqrt((double)n);
    }

    double half_width = t_val * s * sqrt(1.0 + 1.0 / (double)n);
    if (pi_low)  *pi_low  = xbar - half_width;
    if (pi_high) *pi_high = xbar + half_width;
}

/* ─── L3: Required Sample Size ───────────────────────────────────────────
 *
 * n_required = (z_{alpha/2} * sigma / E)^2
 * where E is the desired margin of error.
 */

size_t stats_required_sample_size(double sigma, double margin, double alpha) {
    if (margin < DBL_EPSILON || sigma < DBL_EPSILON) return 0;
    double z_val;
    if (alpha <= 0.01) z_val = 2.576;
    else if (alpha <= 0.05) z_val = 1.96;
    else if (alpha <= 0.10) z_val = 1.645;
    else z_val = 1.0;

    double n_float = (z_val * sigma / margin) * (z_val * sigma / margin);
    size_t n = (size_t)ceil(n_float);
    return (n < 2) ? 2 : n;
}

/* ─── L8: Allan Variance (IEEE Std 1139-2008) ────────────────────────────
 *
 * Allan variance characterizes frequency/time stability as a function
 * of averaging time. It reveals different noise processes.
 *
 * sigma_y^2(tau) = 1 / (2*(M-1)) * SUM_{i=1}^{M-1} (ybar_{i+1} - ybar_i)^2
 * where ybar_i is the average over the i-th interval of length tau,
 * and M = floor(N/tau) is the number of intervals.
 */

double stats_allan_variance(const double *data, size_t n, size_t tau) {
    if (!data || n < 2 || tau == 0) return 0.0;

    size_t M = n / tau;  /* number of intervals */
    if (M < 2) return 0.0;

    /* Compute cluster averages */
    double *clusters = (double*)malloc(M * sizeof(double));
    if (!clusters) return 0.0;

    for (size_t i = 0; i < M; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < tau; j++) {
            sum += data[i * tau + j];
        }
        clusters[i] = sum / (double)tau;
    }

    /* Overlapping Allan variance for better confidence */
    double sum_sq = 0.0;
    size_t count = 0;
    for (size_t i = 0; i + 1 < M; i++) {
        double diff = clusters[i + 1] - clusters[i];
        sum_sq += diff * diff;
        count++;
    }

    free(clusters);
    if (count == 0) return 0.0;
    return sum_sq / (2.0 * (double)count);
}

/* Allan deviation at multiple tau values */
void stats_allan_deviation_multi(const double *data, size_t n,
                                  const size_t *taus, size_t num_taus,
                                  double *adev) {
    if (!data || !taus || !adev) return;
    for (size_t i = 0; i < num_taus; i++) {
        double avar = stats_allan_variance(data, n, taus[i]);
        adev[i] = sqrt(avar);
    }
}

/* ─── L8: Grubbs Test for Outliers (ISO 5725-2) ──────────────────────────
 *
 * G = max_i |x_i - xbar| / s
 * Tests whether the most extreme value is a statistical outlier.
 */

double stats_grubbs_test(const double *data, size_t n, int *outlier_idx) {
    if (!data || n < 3) {
        if (outlier_idx) *outlier_idx = -1;
        return 0.0;
    }

    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += data[i];
    double xbar = sum / (double)n;

    double sum_sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        double d = data[i] - xbar;
        sum_sq += d * d;
    }
    double s = sqrt(sum_sq / (double)(n - 1));
    if (s < DBL_EPSILON) {
        if (outlier_idx) *outlier_idx = -1;
        return 0.0;
    }

    double g_max = 0.0;
    int max_idx = -1;
    for (size_t i = 0; i < n; i++) {
        double g = fabs(data[i] - xbar) / s;
        if (g > g_max) { g_max = g; max_idx = (int)i; }
    }

    if (outlier_idx) *outlier_idx = max_idx;
    return g_max;
}

/* ─── L8: Modified Z-Score Outlier Detection ─────────────────────────────
 *
 * Uses median (robust to outliers) and MAD (median absolute deviation):
 *   M_i = 0.6745 * (x_i - median) / MAD
 *   |M_i| > 3.5 => outlier
 *
 * This is more robust than standard Z-scores when data may contain
 * multiple outliers or come from heavy-tailed distributions.
 */

static int cmp_double_asc(const void *a, const void *b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

static double compute_median(double *sorted, size_t n) {
    if (n % 2 == 1) return sorted[n / 2];
    else return 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]);
}

void stats_modified_zscore_outliers(const double *data, size_t n, int *flags) {
    if (!data || !flags || n < 3) {
        if (flags) {
            for (size_t i = 0; i < n; i++) flags[i] = 0;
        }
        return;
    }

    double *sorted = (double*)malloc(n * sizeof(double));
    if (!sorted) { for (size_t i = 0; i < n; i++) flags[i] = 0; return; }
    memcpy(sorted, data, n * sizeof(double));
    qsort(sorted, n, sizeof(double), cmp_double_asc);

    double median = compute_median(sorted, n);

    /* Compute MAD: median of absolute deviations from median */
    double *abs_devs = (double*)malloc(n * sizeof(double));
    if (!abs_devs) { free(sorted); return; }
    for (size_t i = 0; i < n; i++) {
        abs_devs[i] = fabs(data[i] - median);
    }
    qsort(abs_devs, n, sizeof(double), cmp_double_asc);
    double mad = compute_median(abs_devs, n);
    free(abs_devs);
    free(sorted);

    if (mad < DBL_EPSILON) {
        for (size_t i = 0; i < n; i++) flags[i] = 0;
        return;
    }

    for (size_t i = 0; i < n; i++) {
        double mz = 0.6745 * fabs(data[i] - median) / mad;
        flags[i] = (mz > 3.5) ? 1 : 0;
    }
}
