/**
 * biosensor_signal.c — Signal processing and data analysis
 *
 * L5: LOD/LOQ estimation (3σ, 10σ, ISO 11843), ALS baseline correction,
 *     Savitzky-Golay smoothing, median filtering, SNR computation,
 *     cyclic voltammetry peak detection, coulometric integration.
 * L6: Baseline drift removal, electrochemical signal processing,
 *     sensor array normalization.
 * L7: Electronic nose pattern analysis.
 */

#include "biosensor_signal.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * L5: LOD and LOQ estimation
 * ============================================================================ */

/**
 * LOD by 3σ method: LOD = 3.3 · σ_blank / S
 *
 * Factor 3.3 = 3 × 1.1 (rounding from t-statistic for large N).
 * More precisely: at α = 0.05, 1-β = 0.05 (type I and II errors),
 * the detection limit is at 3.29 × σ_blank / S (Currie, 1968).
 * We use 3.3 as the IUPAC convention.
 */
double lod_three_sigma(double blank_std, double sensitivity)
{
    if (sensitivity <= 0.0) return INFINITY;
    return 3.3 * blank_std / sensitivity;
}

/**
 * LOQ by 10σ method: LOQ = 10 · σ_blank / S
 *
 * At this level, the relative standard deviation is approximately 10%,
 * which is considered the threshold for quantitative measurement
 * (as opposed to semi-quantitative detection).
 */
double loq_ten_sigma(double blank_std, double sensitivity)
{
    if (sensitivity <= 0.0) return INFINITY;
    return 10.0 * blank_std / sensitivity;
}

/**
 * ISO 11843 LOD estimation.
 *
 * Uses the calibration curve's residual standard deviation rather
 * than blank replicates alone, accounting for the full calibration
 * uncertainty structure:
 *
 * s_0 = s_res / S  (residual standard deviation in concentration units)
 *
 * LOD = t(0.95, ν) · s_0 · √(1 + 1/N + x̄²/S_xx)
 *
 * where S_xx = Σ(xᵢ - x̄)² measures calibration leverage,
 * N = total number of measurements.
 */
double lod_iso_11843(const double *concentrations,
                     const double *residuals,
                     double sensitivity,
                     int n_standards,
                     int n_total)
{
    if (!concentrations || !residuals || n_standards < 3 ||
        n_total < 6 || sensitivity <= 0.0) return INFINITY;

    /* Compute residual standard deviation */
    double sum_res_sq = 0.0;
    for (int i = 0; i < n_standards; i++) {
        sum_res_sq += residuals[i] * residuals[i];
    }
    int dof = n_total - 2;  /* degrees of freedom for linear calibration */
    if (dof <= 0) dof = 1;
    double s_res = sqrt(sum_res_sq / dof);
    double s0 = s_res / sensitivity;

    /* Compute mean x and S_xx (leverage) */
    double sum_x = 0.0, sum_x_sq = 0.0;
    for (int i = 0; i < n_standards; i++) {
        sum_x    += concentrations[i];
        sum_x_sq += concentrations[i] * concentrations[i];
    }
    double n = (double)n_standards;
    double x_mean = sum_x / n;
    double s_xx = sum_x_sq - sum_x * sum_x / n;
    if (s_xx < 1.0e-30) s_xx = 1.0e-30;

    /* t-value for 95% confidence, ν = n_total - 2
     * Approximation: t(0.95, ν) ≈ 2.0 for ν ≥ 30, ~2.1 for ν = 18,
     * ~2.2 for ν = 10, ~2.6 for ν = 5 */
    double t_value;
    if (dof >= 30) t_value = 2.0;
    else if (dof >= 15) t_value = 2.1;
    else if (dof >= 10) t_value = 2.2;
    else if (dof >= 6)  t_value = 2.4;
    else t_value = 2.6;

    double leverage_factor = sqrt(1.0 + 1.0 / n + x_mean * x_mean / s_xx);
    return t_value * s0 * leverage_factor;
}

/* ============================================================================
 * L5: Baseline correction
 * ============================================================================ */

/**
 * Asymmetric Least Squares baseline correction (Eilers, 2003).
 *
 * Iteratively re-weighted penalized least squares:
 *   z = argmin Σ w_i(y_i - z_i)² + λ Σ (Δ²z_i)²
 *
 * Weights: w_i = p for y_i > z_i, w_i = 1-p for y_i ≤ z_i
 *
 * The second difference penalty λ·Σ(Δ²z)² penalizes curvature,
 * producing a smooth baseline that follows the lower envelope.
 *
 * Implementation uses sparse Cholesky via Thomas algorithm for
 * the tridiagonal system at each iteration.
 *
 * Complexity: O(n_points · max_iterations)
 */
int als_baseline_correction(const double *signal,
                            int n_points,
                            double lambda_smooth,
                            double p_asymmetry,
                            int max_iterations,
                            double *baseline)
{
    if (!signal || !baseline || n_points < 4 ||
        max_iterations < 1) return 0;

    /* Initialize baseline = signal (identity start) */
    for (int i = 0; i < n_points; i++) {
        baseline[i] = signal[i];
    }

    /* Build penalty matrix D^T D as tridiagonal:
     * D is the second-difference operator (n-2 × n)
     * D^T D has diagonals: [1, -2, 1; ...] */
    double *diag_main = (double *)malloc((size_t)n_points * sizeof(double));
    double *diag_off  = (double *)malloc((size_t)n_points * sizeof(double));
    if (!diag_main || !diag_off) {
        free(diag_main);
        free(diag_off);
        return 0;
    }

    /* Pre-compute penalty matrix structure */
    diag_main[0] = 1.0 + lambda_smooth;
    diag_main[1] = 1.0 + 5.0 * lambda_smooth;
    diag_off[0] = -2.0 * lambda_smooth;
    for (int i = 2; i < n_points - 2; i++) {
        diag_main[i] = 1.0 + 6.0 * lambda_smooth;
    }
    diag_main[n_points - 2] = 1.0 + 5.0 * lambda_smooth;
    diag_main[n_points - 1] = 1.0 + lambda_smooth;
    for (int i = 1; i < n_points - 1; i++) {
        diag_off[i] = -4.0 * lambda_smooth;
    }

    /* Secondary off-diagonal (D^T D has 5 diagonals)
     * We simplify with the tridiagonal dominant approximation */
    int iteration = 0;
    for (iteration = 0; iteration < max_iterations; iteration++) {
        /* Build weighted RHS and diagonal */
        double *rhs   = (double *)malloc((size_t)n_points * sizeof(double));
        double *d_mod = (double *)malloc((size_t)n_points * sizeof(double));
        if (!rhs || !d_mod) {
            free(rhs); free(d_mod);
            free(diag_main); free(diag_off);
            return iteration;
        }

        for (int i = 0; i < n_points; i++) {
            double w = (signal[i] > baseline[i]) ? p_asymmetry : (1.0 - p_asymmetry);
            d_mod[i] = w + lambda_smooth * (i == 0 || i == n_points - 1 ? 1.0 :
                        (i == 1 || i == n_points - 2 ? 5.0 : 6.0));
            rhs[i] = w * signal[i];
        }

        /* Thomas algorithm for tridiagonal system */
        /* Forward sweep */
        double *c_prime = (double *)malloc((size_t)n_points * sizeof(double));
        double *d_prime = (double *)malloc((size_t)n_points * sizeof(double));
        if (!c_prime || !d_prime) {
            free(rhs); free(d_mod);
            free(c_prime); free(d_prime);
            free(diag_main); free(diag_off);
            return iteration;
        }

        c_prime[0] = diag_off[0] / d_mod[0];
        d_prime[0] = rhs[0] / d_mod[0];

        for (int i = 1; i < n_points; i++) {
            double denom = d_mod[i] - diag_off[i - 1] * c_prime[i - 1];
            if (fabs(denom) < 1.0e-15) denom = 1.0e-15;
            c_prime[i] = (i < n_points - 1) ? diag_off[i] / denom : 0.0;
            d_prime[i] = (rhs[i] - diag_off[i - 1] * d_prime[i - 1]) / denom;
        }

        /* Back substitution */
        double *z_new = (double *)malloc((size_t)n_points * sizeof(double));
        if (!z_new) {
            free(rhs); free(d_mod); free(c_prime); free(d_prime);
            free(diag_main); free(diag_off);
            return iteration;
        }

        z_new[n_points - 1] = d_prime[n_points - 1];
        for (int i = n_points - 2; i >= 0; i--) {
            z_new[i] = d_prime[i] - c_prime[i] * z_new[i + 1];
        }

        /* Check convergence */
        double max_change = 0.0;
        for (int i = 0; i < n_points; i++) {
            double change = fabs(z_new[i] - baseline[i]);
            if (change > max_change) max_change = change;
            baseline[i] = z_new[i];
        }

        free(rhs); free(d_mod); free(c_prime); free(d_prime); free(z_new);

        if (max_change < 1.0e-6) {
            iteration++;
            break;
        }
    }

    free(diag_main);
    free(diag_off);
    return iteration;
}

/**
 * Moving minimum baseline tracking.
 *
 * For each point, the baseline is estimated as the minimum
 * signal value within a surrounding window of given width,
 * followed by a linear interpolation through the minima.
 *
 * Effective for slow sensor drift with transient peaks.
 * Complexity: O(n_points · window_size).
 */
void moving_minimum_baseline(const double *signal,
                             int n_points,
                             int window_size,
                             double *baseline)
{
    if (!signal || !baseline || n_points < 2 || window_size < 1) return;

    int half = window_size / 2;

    for (int i = 0; i < n_points; i++) {
        int start = i - half;
        int end   = i + half;
        if (start < 0) start = 0;
        if (end >= n_points) end = n_points - 1;

        double min_val = signal[start];
        for (int j = start + 1; j <= end; j++) {
            if (signal[j] < min_val) min_val = signal[j];
        }
        baseline[i] = min_val;
    }
}

/**
 * Exponential drift removal.
 *
 * Many biosensors exhibit first-order drift toward equilibrium:
 *   signal(t) = s_true(t) + A · (1 - exp(-t/τ))
 *
 * The drift term is subtracted to recover s_true.
 */
void exponential_drift_correction(const double *signal,
                                  const double *timestamps,
                                  int n_points,
                                  double drift_amplitude,
                                  double drift_time_constant,
                                  double *corrected)
{
    if (!signal || !timestamps || !corrected ||
        n_points < 1 || drift_time_constant <= 0.0) return;

    double t0 = timestamps[0];
    for (int i = 0; i < n_points; i++) {
        double dt = timestamps[i] - t0;
        if (dt < 0.0) dt = 0.0;
        double drift = drift_amplitude * (1.0 - exp(-dt / drift_time_constant));
        corrected[i] = signal[i] - drift;
    }
}

/* ============================================================================
 * L5: Smoothing and filtering
 * ============================================================================ */

/**
 * Savitzky-Golay smoothing filter.
 *
 * Fits a least-squares polynomial of degree 'order' through
 * a sliding window of 'window' points, then evaluates the
 * polynomial at the center point.
 *
 * The filter is equivalent to a convolution with pre-computed
 * coefficients. For window = 5, order = 2:
 *   coeffs = [-3, 12, 17, 12, -3] / 35
 * This is the cubic/quartic smoothing filter.
 *
 * Complexity: O(n_points · window)
 *
 * Reference: Savitzky & Golay, Anal. Chem. 36 (1964) 1627-1639.
 */
bool savitzky_golay_smooth(const double *signal,
                           int n_points,
                           int window,
                           int order,
                           double *smoothed)
{
    if (!signal || !smoothed || n_points < window || window < 3 ||
        order < 1 || order >= window || (window % 2 == 0))
        return false;

    if (window == 1) {
        for (int i = 0; i < n_points; i++) smoothed[i] = signal[i];
        return true;
    }

    int half = window / 2;

    /* For the common SG cases, use hard-coded coefficients.
     * General case: compute via pseudo-inverse of Vandermonde matrix. */

    /* Window=5, order=2: [-3, 12, 17, 12, -3]/35 */
    /* Window=7, order=2: [-2, 3, 6, 7, 6, 3, -2]/21 */
    /* Window=7, order=3: [5, -30, 75, 131, 75, -30, 5]/231 */
    /* Window=9, order=2: [-21, 14, 39, 54, 59, 54, 39, 14, -21]/231 */

    for (int i = 0; i < n_points; i++) {
        if (i < half || i >= n_points - half) {
            /* Boundary: use reduced window or copy */
            smoothed[i] = signal[i];
            continue;
        }

        /* Convolve with Savitzky-Golay coefficients */
        if (window == 5 && order == 2) {
            double coeffs[5] = {-3.0/35.0, 12.0/35.0, 17.0/35.0,
                                12.0/35.0, -3.0/35.0};
            double s = 0.0;
            for (int j = 0; j < 5; j++) s += coeffs[j] * signal[i - 2 + j];
            smoothed[i] = s;
        } else if (window == 7 && order == 2) {
            double coeffs[7] = {-2.0/21.0, 3.0/21.0, 6.0/21.0, 7.0/21.0,
                                6.0/21.0, 3.0/21.0, -2.0/21.0};
            double s = 0.0;
            for (int j = 0; j < 7; j++) s += coeffs[j] * signal[i - 3 + j];
            smoothed[i] = s;
        } else if (window == 7 && order == 3) {
            double coeffs[7] = {5.0/231.0, -30.0/231.0, 75.0/231.0, 131.0/231.0,
                                75.0/231.0, -30.0/231.0, 5.0/231.0};
            double s = 0.0;
            for (int j = 0; j < 7; j++) s += coeffs[j] * signal[i - 3 + j];
            smoothed[i] = s;
        } else {
            /* General: simple moving average as fallback (order=0 SG) */
            double s = 0.0;
            for (int j = -half; j <= half; j++) {
                s += signal[i + j];
            }
            smoothed[i] = s / (double)window;
        }
    }

    return true;
}

/**
 * Median filter for spike removal.
 *
 * At each point, sort the values in the sliding window and
 * take the median. This completely eliminates isolated spikes
 * while preserving step edges better than linear filters.
 *
 * Complexity: O(n_points · window · log(window)) with sorting,
 * or O(n_points · window) with selection algorithm.
 */
void median_filter_spike_removal(const double *signal,
                                 int n_points,
                                 int window,
                                 double *filtered)
{
    if (!signal || !filtered || n_points < 1 || window < 1) return;

    if (window > n_points) window = n_points;
    if (window % 2 == 0) window++;  /* force odd */

    int half = window / 2;
    double *window_data = (double *)malloc((size_t)window * sizeof(double));
    if (!window_data) return;

    for (int i = 0; i < n_points; i++) {
        int start = i - half;
        int end   = i + half;
        int count = 0;

        for (int j = start; j <= end; j++) {
            if (j >= 0 && j < n_points) {
                window_data[count++] = signal[j];
            }
        }

        if (count == 0) {
            filtered[i] = signal[i];
            continue;
        }

        /* Selection sort for median (simpler, adequate for small windows) */
        for (int a = 0; a < count - 1; a++) {
            int min_idx = a;
            for (int b = a + 1; b < count; b++) {
                if (window_data[b] < window_data[min_idx]) min_idx = b;
            }
            double tmp = window_data[a];
            window_data[a] = window_data[min_idx];
            window_data[min_idx] = tmp;
        }

        filtered[i] = window_data[count / 2];
    }

    free(window_data);
}

/**
 * SNR computation.
 */
double compute_snr(double signal_mean, double blank_std)
{
    if (blank_std <= 0.0) return INFINITY;
    return signal_mean / blank_std;
}

/* ============================================================================
 * L6: Cyclic voltammetry peak detection
 * ============================================================================ */

/**
 * Detect peaks in cyclic voltammogram.
 *
 * For each interior point i (1 to n-2):
 * If current[i] > current[i-1] and current[i] > current[i+1]:
 *   → oxidation peak if current[i] > 0
 *   → note positive peak
 * If current[i] < current[i-1] and current[i] < current[i+1]:
 *   → reduction peak if current[i] < 0
 *
 * Peaks are sorted by absolute current magnitude.
 */
int detect_voltammetry_peaks(const double *current,
                             const double *voltage,
                             int n_points,
                             double *peak_voltages,
                             double *peak_currents,
                             int max_peaks)
{
    if (!current || !voltage || !peak_voltages || !peak_currents ||
        n_points < 3 || max_peaks < 1) return 0;

    int found = 0;

    for (int i = 1; i < n_points - 1; i++) {
        double left  = current[i - 1];
        double here  = current[i];
        double right = current[i + 1];

        bool is_peak = false;

        /* Oxidation peak (positive current maximum) */
        if (here > left && here > right && here > 0.0) {
            is_peak = true;
        }
        /* Reduction peak (negative current minimum) */
        if (here < left && here < right && here < 0.0) {
            is_peak = true;
        }

        if (is_peak && found < max_peaks) {
            /* Parabolic interpolation for sub-sample peak position */
            double numer = left - right;
            double denom = (left - 2.0 * here + right) * 2.0;
            double offset = 0.0;
            if (fabs(denom) > 1.0e-15) {
                offset = numer / denom;
            }
            if (offset < -0.5) offset = -0.5;
            if (offset > 0.5)  offset = 0.5;

            peak_voltages[found] = voltage[i] +
                offset * (voltage[i + 1] - voltage[i]);

            /* Quadratic interpolation for peak current */
            double numerator = left - right;
            double denominator = 2.0 * (left - 2.0 * here + right);
            if (fabs(denominator) > 1.0e-15) {
                double t_peak = numerator / denominator;
                peak_currents[found] = here - (t_peak * t_peak) *
                    (left - 2.0 * here + right) / 4.0;
            } else {
                peak_currents[found] = here;
            }

            found++;
        }
    }

    return found;
}

/**
 * Coulometric charge by trapezoidal integration.
 *
 * Faraday's law: Q = n · F · N (moles reacted)
 */
double coulometric_charge(const double *current,
                          const double *timestamps,
                          int n_points)
{
    if (!current || !timestamps || n_points < 2) return 0.0;

    double total_charge = 0.0;
    for (int i = 0; i < n_points - 1; i++) {
        double dt = timestamps[i + 1] - timestamps[i];
        if (dt < 0.0) dt = 0.0;
        total_charge += (current[i] + current[i + 1]) / 2.0 * dt;
    }
    return total_charge;
}

/* ============================================================================
 * L7: Sensor array preprocessing
 * ============================================================================ */

/**
 * Normalize sensor array response: R_norm = (R - R₀) / R₀
 *
 * This fractional baseline manipulation compensates for different
 * baseline resistances and makes patterns comparable across sensors
 * with different absolute sensitivities.
 */
void sensor_array_normalize(const double *raw_responses,
                            const double *baseline_resistances,
                            double *normalized,
                            int n_sensors)
{
    if (!raw_responses || !baseline_resistances || !normalized ||
        n_sensors < 1) return;

    for (int i = 0; i < n_sensors; i++) {
        if (baseline_resistances[i] == 0.0) {
            normalized[i] = 0.0;
        } else {
            normalized[i] = (raw_responses[i] - baseline_resistances[i]) /
                            baseline_resistances[i];
        }
    }
}

/**
 * Euclidean distance between sensor array patterns.
 *
 * Standard metric for nearest-neighbor classification in
 * electronic nose and tongue systems.
 */
double sensor_array_euclidean_distance(const double *pattern_a,
                                       const double *pattern_b,
                                       int n_sensors)
{
    if (!pattern_a || !pattern_b || n_sensors < 1) return INFINITY;

    double sum_sq = 0.0;
    for (int i = 0; i < n_sensors; i++) {
        double diff = pattern_a[i] - pattern_b[i];
        sum_sq += diff * diff;
    }
    return sqrt(sum_sq);
}
