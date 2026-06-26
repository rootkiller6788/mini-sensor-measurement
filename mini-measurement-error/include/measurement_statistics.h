#ifndef MEASUREMENT_STATISTICS_H
#define MEASUREMENT_STATISTICS_H

#include <stddef.h>
#include "measurement_error.h"

/**
 * @file measurement_statistics.h
 * @brief L3: Statistical distributions, hypothesis testing for measurements
 *
 * L3 Mathematical Structures:
 *   - Normal/Gaussian distribution fitting and testing
 *   - Chi-squared distribution for variance testing
 *   - Student t-distribution for mean comparisons
 *   - F-distribution for variance ratio tests
 *
 * L8 Advanced Topics:
 *   - Allan variance for noise characterization (IEEE Std 1139)
 *   - Grubbs test for outlier detection (ISO 5725-2)
 *
 * Courses: MIT 6.431 (Probability), Stanford STATS 200, ETH 401-2604
 */

/* --- L3: Distribution Fitting --- */

/** MLE of Gaussian parameters from data */
void stats_fit_gaussian(const double *data, size_t n,
                         double *mu_hat, double *sigma_hat);

/** Kolmogorov-Smirnov test statistic for normality */
double stats_ks_test_normal(const double *data, size_t n,
                             double mu, double sigma);

/** Shapiro-Wilk W statistic for normality (simplified, n <= 50) */
double stats_shapiro_wilk(const double *data, size_t n);

/* --- L3: Hypothesis Testing --- */

/** One-sample t-test: H0: mu = hypothesized_mu */
double stats_one_sample_t_test(const double *data, size_t n,
                                double hypothesized_mu);

/** Two-sample independent t-test (Welch, unequal variances) */
double stats_two_sample_t_test(const double *data1, size_t n1,
                                const double *data2, size_t n2);

/** F-test for equality of variances: H0: sigma1^2 = sigma2^2 */
double stats_f_test_variances(const double *data1, size_t n1,
                               const double *data2, size_t n2);

/** Chi-squared test for variance: H0: sigma^2 = sigma0^2 */
double stats_chi_sq_variance_test(const double *data, size_t n,
                                   double sigma0_sq);

/* --- L3: Confidence and Prediction Intervals --- */

/** Prediction interval for a single future measurement */
void stats_prediction_interval(const double *data, size_t n, double alpha,
                                double *pi_low, double *pi_high);

/** Sample size for desired margin of error */
size_t stats_required_sample_size(double sigma, double margin, double alpha);

/* --- L8: Allan Variance (IEEE Std 1139-2008) --- */

/** Allan variance at timescale tau (number of samples) */
double stats_allan_variance(const double *data, size_t n, size_t tau);

/** Overlapping Allan deviation at multiple timescales */
void stats_allan_deviation_multi(const double *data, size_t n,
                                  const size_t *taus, size_t num_taus,
                                  double *adev);

/* --- L8: Outlier Detection --- */

/** Grubbs test for a single outlier (ISO 5725-2) */
double stats_grubbs_test(const double *data, size_t n, int *outlier_idx);

/** Modified Z-score outlier detection (robust, uses median/MAD) */
void stats_modified_zscore_outliers(const double *data, size_t n, int *flags);

#endif /* MEASUREMENT_STATISTICS_H */
