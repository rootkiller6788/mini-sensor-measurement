#include "measurement_filtering.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ─── L5: Moving Average Filter (FIR) ────────────────────────────────────
 *
 * Simple moving average (SMA) with a sliding window.
 * y[n] = (1/N) * SUM_{k=0}^{N-1} x[n-k]
 *
 * Maintains a ring buffer and running sum for O(1) per-sample update.
 * Ideal for reducing white noise in steady-state measurements.
 * Noise reduction: sigma_out = sigma_in / sqrt(N)
 */

int mavg_init(MovingAvgFilter *filt, size_t window_size) {
    if (!filt || window_size == 0) return -1;
    filt->buffer = (double*)calloc(window_size, sizeof(double));
    if (!filt->buffer) return -1;
    filt->window_size = window_size;
    filt->index = 0;
    filt->count = 0;
    filt->sum = 0.0;
    filt->initialized = 1;
    return 0;
}

void mavg_free(MovingAvgFilter *filt) {
    if (!filt) return;
    free(filt->buffer);
    filt->buffer = NULL;
    filt->initialized = 0;
}

void mavg_reset(MovingAvgFilter *filt) {
    if (!filt || !filt->buffer) return;
    memset(filt->buffer, 0, filt->window_size * sizeof(double));
    filt->index = 0;
    filt->count = 0;
    filt->sum = 0.0;
}

double mavg_update(MovingAvgFilter *filt, double sample) {
    if (!filt || !filt->buffer || filt->window_size == 0) return sample;

    if (filt->count < filt->window_size) {
        /* Still filling the window */
        filt->buffer[filt->index] = sample;
        filt->sum += sample;
        filt->count++;
        filt->index = (filt->index + 1) % filt->window_size;
        return filt->sum / (double)filt->count;
    } else {
        /* Window full: subtract oldest, add newest */
        double oldest = filt->buffer[filt->index];
        filt->buffer[filt->index] = sample;
        filt->sum = filt->sum - oldest + sample;
        filt->index = (filt->index + 1) % filt->window_size;
        return filt->sum / (double)filt->window_size;
    }
}

/* ─── L5: Exponential Moving Average (1st-order IIR) ────────────────────
 *
 * y[k] = alpha * x[k] + (1 - alpha) * y[k-1]
 *
 * Time constant tau relates to alpha:
 *   alpha = 1 - exp(-Ts / tau)  where Ts = sample period
 *   For Ts << tau:  alpha ~= Ts / tau
 *
 * 3dB cutoff: f_c = alpha / (2 * pi * Ts) (approximately for alpha << 1)
 * Equivalent window: N_eq = (2 - alpha) / alpha
 */

void ema_init(EMAFilter *filt, double alpha) {
    if (!filt) return;
    if (alpha <= 0.0) alpha = 0.5;
    if (alpha > 1.0)  alpha = 1.0;
    filt->alpha = alpha;
    filt->output = 0.0;
    filt->initialized = 0;
}

double ema_update(EMAFilter *filt, double sample) {
    if (!filt) return sample;

    if (!filt->initialized) {
        filt->output = sample;
        filt->initialized = 1;
        return sample;
    }

    filt->output = filt->alpha * sample
                   + (1.0 - filt->alpha) * filt->output;
    return filt->output;
}

double ema_alpha_from_tau(double sample_period, double time_constant) {
    if (time_constant < DBL_EPSILON) return 1.0;
    if (sample_period < DBL_EPSILON) return 1.0;
    return 1.0 - exp(-sample_period / time_constant);
}

/* ─── L5: Median Filter (Nonlinear) ──────────────────────────────────────
 *
 * Median filter is optimal for removing impulsive (salt-and-pepper) noise
 * while preserving edges. Uses quickselect algorithm.
 */

static void swap_double(double *a, double *b) {
    double tmp = *a; *a = *b; *b = tmp;
}

static double quickselect(double *arr, size_t n, size_t k) {
    if (n <= 1) return arr[0];

    /* Use middle element as pivot */
    double pivot = arr[n / 2];

    /* Three-way partition */
    size_t lt = 0, eq = 0, gt = n;
    while (eq < gt) {
        if (arr[eq] < pivot) {
            swap_double(&arr[lt], &arr[eq]);
            lt++; eq++;
        } else if (arr[eq] > pivot) {
            gt--;
            swap_double(&arr[eq], &arr[gt]);
        } else {
            eq++;
        }
    }

    if (k < lt) return quickselect(arr, lt, k);
    if (k < eq) return pivot;
    return quickselect(arr + eq, n - eq, k - eq);
}

double median_filter(const double *data, size_t n) {
    if (!data || n == 0) return 0.0;
    if (n == 1) return data[0];

    double *copy = (double*)malloc(n * sizeof(double));
    if (!copy) return 0.0;
    memcpy(copy, data, n * sizeof(double));

    double result = quickselect(copy, n, n / 2);
    if (n % 2 == 0) {
        double second = quickselect(copy, n / 2, n / 2 - 1);
        result = 0.5 * (result + second);
    }
    free(copy);
    return result;
}

void running_median_filter(const double *input, double *output,
                           size_t length, size_t window_size) {
    if (!input || !output || length == 0 || window_size < 3) return;

    /* Make window_size odd for true median */
    if (window_size % 2 == 0) window_size++;

    for (size_t i = 0; i < length; i++) {
        size_t start, end;
        if (i < window_size / 2) {
            start = 0;
            end = window_size;
        } else if (i + window_size / 2 >= length) {
            start = length - window_size;
            end = length;
        } else {
            start = i - window_size / 2;
            end = i + window_size / 2 + 1;
        }
        if (end > length) end = length;
        output[i] = median_filter(&input[start], end - start);
    }
}

/* ─── L5: 1D Kalman Filter for Sensor Measurements ──────────────────────
 *
 * State model: x_k = F * x_{k-1} + w_k    (w ~ N(0, Q))
 * Measurement: z_k = H * x_k + v_k        (v ~ N(0, R))
 *
 * Predict:
 *   x_hat^- = F * x_hat
 *   P^- = F * P * F + Q
 *
 * Update:
 *   K = P^- * H / (H * P^- * H + R)
 *   x_hat = x_hat^- + K * (z - H * x_hat^-)
 *   P = (1 - K * H) * P^-
 *
 * For simple 1D case with F=1, H=1:
 *   K = P^- / (P^- + R)
 *   x_hat = x_hat^- + K * (z - x_hat^-)
 *   P = (1 - K) * P^-
 */

void kf1d_init(KalmanFilter1D *kf, double x0, double P0,
               double Q, double R, double F) {
    if (!kf) return;
    kf->x_hat = x0;
    kf->P = P0;
    kf->Q = Q;
    kf->R = R;
    kf->F = F;
    kf->H = 1.0;
    kf->initialized = 1;
}

void kf1d_predict(KalmanFilter1D *kf) {
    if (!kf || !kf->initialized) return;
    kf->x_hat = kf->F * kf->x_hat;
    kf->P = kf->F * kf->P * kf->F + kf->Q;
}

double kf1d_update(KalmanFilter1D *kf, double measurement) {
    if (!kf || !kf->initialized) return measurement;

    /* Innovation (measurement residual) */
    double y = measurement - kf->H * kf->x_hat;

    /* Innovation covariance */
    double S = kf->H * kf->P * kf->H + kf->R;

    /* Kalman gain */
    double K = kf->P * kf->H / S;

    /* State update */
    kf->x_hat = kf->x_hat + K * y;

    /* Covariance update (Joseph form for numerical stability) */
    double I_KH = 1.0 - K * kf->H;
    kf->P = I_KH * kf->P * I_KH + K * kf->R * K;

    return kf->x_hat;
}

/* ─── L5: Sensor Fusion ────────────────────────────────────────────────── */

double sensor_fuse_weighted_average(const double *readings,
                                     const double *weights, size_t n) {
    if (!readings || n == 0) return 0.0;
    if (!weights) {
        /* Equal weights if none provided */
        double sum = 0.0;
        for (size_t i = 0; i < n; i++) sum += readings[i];
        return sum / (double)n;
    }

    double sum_w = 0.0, sum_wx = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum_w  += weights[i];
        sum_wx += weights[i] * readings[i];
    }
    return (sum_w > DBL_EPSILON) ? (sum_wx / sum_w) : 0.0;
}

/* Inverse-variance weighted fusion.
 * The optimal linear combination for independent Gaussian sensors.
 * w_i = (1/sigma_i^2) / SUM_j (1/sigma_j^2)
 * var_fused = 1 / SUM_j (1/sigma_j^2)
 */

double sensor_fuse_inverse_variance(const double *readings,
                                     const double *variances, size_t n,
                                     double *fused_var) {
    if (!readings || !variances || n == 0) {
        if (fused_var) *fused_var = INFINITY;
        return 0.0;
    }

    double *inv_vars = (double*)malloc(n * sizeof(double));
    if (!inv_vars) {
        if (fused_var) *fused_var = INFINITY;
        return 0.0;
    }

    double sum_inv = 0.0;
    double sum_inv_x = 0.0;
    for (size_t i = 0; i < n; i++) {
        double iv = (variances[i] > DBL_EPSILON) ? 1.0 / variances[i] : 0.0;
        inv_vars[i] = iv;
        sum_inv += iv;
        sum_inv_x += iv * readings[i];
    }

    double fused = (sum_inv > DBL_EPSILON) ? (sum_inv_x / sum_inv) : 0.0;
    if (fused_var) {
        *fused_var = (sum_inv > DBL_EPSILON) ? (1.0 / sum_inv) : INFINITY;
    }
    free(inv_vars);
    return fused;
}

double sensor_fuse_variance_reduction(size_t n_sensors) {
    return (n_sensors > 0) ? (double)n_sensors : 1.0;
}
