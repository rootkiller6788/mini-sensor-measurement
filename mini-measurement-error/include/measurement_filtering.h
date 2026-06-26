#ifndef MEASUREMENT_FILTERING_H
#define MEASUREMENT_FILTERING_H

#include <stddef.h>

/**
 * @file measurement_filtering.h
 * @brief L5: Noise reduction and sensor fusion algorithms for measurements
 *
 * L5 Algorithms/Methods:
 *   - Moving Average (FIR) filtering for steady-state measurements
 *   - Exponential Moving Average (EMA / 1st-order IIR)
 *   - Median filter for impulsive noise removal
 *   - 1D Kalman filter for dynamic measurements
 *   - Sensor fusion: weighted averaging of redundant sensors
 *
 * Courses: MIT 6.003 (filtering), Stanford EE359 (Kalman), Berkeley EE123 (DSP)
 */

/* --- L5: Moving Average Filter --- */

typedef struct {
    double *buffer;
    size_t  window_size;
    size_t  index;
    size_t  count;
    double  sum;
    int     initialized;
} MovingAvgFilter;

/** Initialize moving average filter; allocates window buffer */
int mavg_init(MovingAvgFilter *filt, size_t window_size);

/** Free filter resources */
void mavg_free(MovingAvgFilter *filt);

/** Reset filter state (reuse buffer) */
void mavg_reset(MovingAvgFilter *filt);

/** Insert new sample and return filtered output */
double mavg_update(MovingAvgFilter *filt, double sample);

/* --- L5: Exponential Moving Average (1st-order IIR lowpass) --- */

typedef struct {
    double alpha;       /* smoothing factor in (0, 1] */
    double output;      /* current EMA value */
    int    initialized;
} EMAFilter;

/** Initialize EMA filter with smoothing factor alpha */
void ema_init(EMAFilter *filt, double alpha);

/** Update EMA: y[k] = alpha * x[k] + (1-alpha) * y[k-1] */
double ema_update(EMAFilter *filt, double sample);

/** Compute alpha from time constant: alpha = 1 - exp(-Ts/tau) */
double ema_alpha_from_tau(double sample_period, double time_constant);

/* --- L5: Median Filter (Nonlinear, for impulsive noise) --- */

/**
 * Compute median of a data window using quickselect algorithm.
 * O(n) average time complexity.
 *
 * @param data     array of samples
 * @param n        number of samples
 * @return median value
 */
double median_filter(const double *data, size_t n);

/**
 * Running median filter: apply median filter to a sliding window
 * over the input signal.
 *
 * @param input        input signal
 * @param output       filtered output (must be pre-allocated, same length)
 * @param length       signal length
 * @param window_size  odd number, typically 3, 5, or 7
 */
void running_median_filter(const double *input, double *output,
                           size_t length, size_t window_size);

/* --- L5: 1-D Kalman Filter for Sensor Measurements --- */

typedef struct {
    double x_hat;          /* state estimate */
    double P;              /* error covariance */
    double Q;              /* process noise covariance */
    double R;              /* measurement noise covariance */
    double F;              /* state transition (default 1.0 for constant model) */
    double H;              /* measurement matrix (default 1.0) */
    int    initialized;
} KalmanFilter1D;

/**
 * Initialize 1D Kalman filter.
 *
 * @param kf     filter struct
 * @param x0     initial state estimate
 * @param P0     initial covariance (larger = less confident in x0)
 * @param Q      process noise (larger = trust model less, follow measurements more)
 * @param R      measurement noise (larger = trust sensor less, smoother output)
 * @param F      state transition (1.0 for constant, dt for velocity)
 */
void kf1d_init(KalmanFilter1D *kf, double x0, double P0,
               double Q, double R, double F);

/**
 * Kalman predict step (called between measurements).
 * Projects state and covariance forward in time.
 */
void kf1d_predict(KalmanFilter1D *kf);

/**
 * Kalman update step (called when new measurement arrives).
 * Corrects state estimate using the measurement.
 *
 * @param kf            filter
 * @param measurement   new sensor reading z_k
 * @return updated state estimate x_hat
 */
double kf1d_update(KalmanFilter1D *kf, double measurement);

/* --- L5: Sensor Fusion --- */

/**
 * Simple weighted average fusion for redundant sensors.
 * weights should sum to 1.0 (or will be normalized).
 *
 * x_fused = SUM w_i * x_i
 *
 * @param readings  readings from N sensors
 * @param weights   weights (typically inverse of variance)
 * @param n         number of sensors
 * @return fused estimate
 */
double sensor_fuse_weighted_average(const double *readings,
                                     const double *weights, size_t n);

/**
 * Inverse-variance weighted fusion (optimal for independent Gaussian sensors).
 * w_i = (1 / sigma_i^2) / SUM (1 / sigma_j^2)
 *
 * @param readings  sensor readings
 * @param variances measurement variances for each sensor
 * @param n         number of sensors
 * @param fused_var [out] variance of fused estimate
 * @return fused estimate
 */
double sensor_fuse_inverse_variance(const double *readings,
                                     const double *variances, size_t n,
                                     double *fused_var);

/**
 * Compute the variance reduction factor from sensor fusion.
 * For N independent sensors with equal variance:
 *   var_fused = var_single / N
 *   reduction_factor = N
 */
double sensor_fuse_variance_reduction(size_t n_sensors);

#endif /* MEASUREMENT_FILTERING_H */
