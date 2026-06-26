#ifndef SENSOR_FUSION_H
#define SENSOR_FUSION_H
#include "smart_sensor.h"
#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * Sensor Fusion Module
 *
 * Implements multi-sensor data fusion algorithms:
 *   L5: Kalman filter (1D/3D), complementary filter, weighted
 *       average, exponential smoothing, median fusion
 *   L4: Kalman-Bucy theory, Bayesian sensor fusion
 *   L6: IMU orientation estimation, multi-sensor temperature fusion
 *   L8: Mahalanobis distance outlier rejection, Dempster-Shafer
 *
 * Reference:
 *   Kalman, R.E. (1960), "A New Approach to Linear Filtering and
 *     Prediction Problems", Trans. ASME J. Basic Eng., 82:35-45
 *   Maybeck, P.S. (1979), "Stochastic Models, Estimation and Control",
 *     Academic Press, Vol. 1
 *   Mahony, R. et al. (2008), "Nonlinear Complementary Filters on
 *     the Special Orthogonal Group", IEEE TAC 53(5):1203-1218
 *   Hall, D.L. & Llinas, J. (2001), "Handbook of Multisensor Data Fusion",
 *     CRC Press
 *
 * Course Mapping:
 *   MIT 6.450 — Digital Communications (estimation theory)
 *   Stanford EE363 — Linear Dynamical Systems (Kalman filtering)
 *   Berkeley EE221A — Linear Systems (observer theory)
 *   ETH 227-0436 — Communications (stochastic estimation)
 * ========================================================================= */

/* =========================================================================
 * L5: 1-Dimensional Kalman Filter for Scalar Sensor Fusion
 *
 * State-space model:
 *   x_k = A * x_{k-1} + B * u_k + w_k    (process, w_k ~ N(0, Q))
 *   z_k = H * x_k + v_k                   (measurement, v_k ~ N(0, R))
 *
 * Kalman gain:      K_k = P_k^- * H / (H * P_k^- * H + R)
 * State update:     x_k = x_k^- + K_k * (z_k - H * x_k^-)
 * Covariance update: P_k = (1 - K_k * H) * P_k^-
 * ========================================================================= */

/**
 * @brief 1D Kalman filter state
 */
typedef struct {
    double x;          /* State estimate */
    double P;          /* Error covariance */
    double A;          /* State transition (default 1.0 for constant model) */
    double H;          /* Observation matrix (default 1.0 for direct meas.) */
    double Q;          /* Process noise covariance */
    double R;          /* Measurement noise covariance */
    uint64_t iteration;/* Filter iteration count */
} ss_kalman_1d_t;

/**
 * @brief Initialize 1D Kalman filter
 *
 * @param kf     Filter state
 * @param x0     Initial state estimate
 * @param P0     Initial error covariance (large = unsure of initial state)
 * @param Q      Process noise covariance (models system uncertainty)
 * @param R      Measurement noise covariance (models sensor noise)
 */
void ss_fusion_kalman_1d_init(ss_kalman_1d_t *kf, double x0,
                              double P0, double Q, double R);

/**
 * @brief Kalman filter prediction step (time update)
 *
 * Projects state and covariance forward in time:
 *   x_k^- = A * x_{k-1}
 *   P_k^- = A * P_{k-1} * A + Q
 *
 * @param kf  Filter state (updated in place)
 */
void ss_fusion_kalman_1d_predict(ss_kalman_1d_t *kf);

/**
 * @brief Kalman filter update step (measurement update)
 *
 * Incorporates new measurement to correct the state estimate:
 *   K = P_k^- * H / (H * P_k^- * H + R)
 *   x_k = x_k^- + K * (z - H * x_k^-)
 *   P_k = (1 - K * H) * P_k^-
 *
 * @param kf  Filter state (updated in place)
 * @param z   New measurement
 * @return    Updated state estimate
 */
double ss_fusion_kalman_1d_update(ss_kalman_1d_t *kf, double z);

/**
 * @brief Combined predict + update in one step
 * @param kf  Filter state
 * @param z   New measurement
 * @return    Updated state estimate
 */
double ss_fusion_kalman_1d_step(ss_kalman_1d_t *kf, double z);

/* =========================================================================
 * L5: 3-Axis Kalman Filter for Orientation Estimation
 * ========================================================================= */

/**
 * @brief 3-axis Kalman filter for IMU sensor fusion
 *
 * Estimates roll, pitch angles from gyroscope and accelerometer.
 * State: [roll, pitch]^T
 */
typedef struct {
    double x[2];        /* State vector [roll, pitch] in radians */
    double P[4];        /* 2x2 covariance matrix (row-major) */
    double Q_gyro[4];   /* Gyroscope process noise (2x2) */
    double R_accel[4];  /* Accelerometer measurement noise (2x2) */
    double dt;          /* Time step in seconds */
    uint64_t iteration;
} ss_kalman_3d_imu_t;

/**
 * @brief Initialize 3-axis IMU Kalman filter
 *
 * @param kf       Filter state
 * @param dt       Sampling period in seconds
 * @param q_gyro   Gyroscope noise variance (rad^2/s^2)
 * @param r_accel  Accelerometer noise variance (g^2)
 */
void ss_fusion_kalman_imu_init(ss_kalman_3d_imu_t *kf, double dt,
                               double q_gyro, double r_accel);

/**
 * @brief IMU Kalman predict step using gyroscope rates
 *
 * State prediction using gyroscope integration:
 *   roll_pred  = roll + dt * (gx + gy*sin(roll)*tan(pitch) + gz*cos(roll)*tan(pitch))
 *   pitch_pred = pitch + dt * (gy*cos(roll) - gz*sin(roll))
 *
 * @param kf   Filter state
 * @param gx   Gyroscope x-axis angular rate (rad/s)
 * @param gy   Gyroscope y-axis angular rate (rad/s)
 * @param gz   Gyroscope z-axis angular rate (rad/s)
 */
void ss_fusion_kalman_imu_predict(ss_kalman_3d_imu_t *kf,
                                  double gx, double gy, double gz);

/**
 * @brief IMU Kalman update using accelerometer measurements
 *
 * Gravity vector measurement model:
 *   roll_meas  = atan2(ay, az)
 *   pitch_meas = atan2(-ax, sqrt(ay^2 + az^2))
 *
 * @param kf  Filter state
 * @param ax  Accelerometer x-axis (g)
 * @param ay  Accelerometer y-axis (g)
 * @param az  Accelerometer z-axis (g)
 */
void ss_fusion_kalman_imu_update(ss_kalman_3d_imu_t *kf,
                                 double ax, double ay, double az);

/**
 * @brief Read current roll/pitch estimates
 * @param kf    Filter state
 * @param roll  Output: roll angle in radians
 * @param pitch Output: pitch angle in radians
 */
void ss_fusion_kalman_imu_get_angles(const ss_kalman_3d_imu_t *kf,
                                     double *roll, double *pitch);

/* =========================================================================
 * L5: Complementary Filter for IMU Orientation
 *
 * Advantages over Kalman: lower computational cost, no covariance
 * matrices, single tuning parameter alpha.
 *
 * angle = alpha * (angle + gyro*dt) + (1-alpha) * accel_angle
 *
 * Low-pass filters accelerometer (noise), high-pass filters gyroscope (drift).
 *
 * Reference: Mahony et al. (2008), IEEE TAC 53(5):1203-1218
 * ========================================================================= */

/**
 * @brief Complementary filter state for IMU
 */
typedef struct {
    double roll;         /* Estimated roll angle (rad) */
    double pitch;        /* Estimated pitch angle (rad) */
    double alpha;        /* Filter coefficient (0-1, typically 0.98) */
    double dt;           /* Time step (s) */
    uint64_t iteration;
} ss_complementary_filter_t;

/**
 * @brief Initialize complementary filter
 *
 * @param cf     Filter state
 * @param alpha  Filter coefficient (0-1)
 *               alpha=1.0: trust gyroscope entirely
 *               alpha=0.0: trust accelerometer entirely
 *               Typical: 0.98 for most MEMS IMUs
 * @param dt     Sampling period in seconds
 */
void ss_fusion_complementary_init(ss_complementary_filter_t *cf,
                                  double alpha, double dt);

/**
 * @brief Complementary filter update step
 *
 * Combines gyroscope integration (short-term) with accelerometer
 * measurement (long-term reference).
 *
 * @param cf    Filter state
 * @param gx    Gyroscope x rate (rad/s)
 * @param gy    Gyroscope y rate (rad/s)
 * @param gz    Gyroscope z rate (rad/s)
 * @param ax    Accelerometer x (g)
 * @param ay    Accelerometer y (g)
 * @param az    Accelerometer z (g)
 */
void ss_fusion_complementary_update(ss_complementary_filter_t *cf,
                                    double gx, double gy, double gz,
                                    double ax, double ay, double az);

/* =========================================================================
 * L5: Weighted Average and Statistical Sensor Fusion
 * ========================================================================= */

/**
 * @brief Inverse-variance weighted average of multiple sensor readings
 *
 * Optimal fusion under the assumption of independent Gaussian errors.
 * Weight for sensor i:
 *   w_i = (1/sigma_i^2) / sum_j (1/sigma_j^2)
 *
 * The fused estimate has minimum variance:
 *   sigma_fused^2 = 1 / sum_i (1/sigma_i^2)
 *
 * Reference: Cochran, W.G. (1937), "Problems arising in the analysis
 *            of a series of similar experiments", JRSS Suppl. 4:102-118
 *
 * @param readings          Array of sensor readings [n_sensors]
 * @param variances         Array of sensor variances [n_sensors]
 * @param n_sensors         Number of sensors
 * @param fused_variance    Output: variance of fused estimate (can be NULL)
 * @return                  Optimal fused estimate
 */
double ss_fusion_weighted_average(const double *readings,
                                  const double *variances,
                                  size_t n_sensors,
                                  double *fused_variance);

/**
 * @brief Median-based fusion for outlier rejection
 *
 * Uses the median as a robust estimate of central tendency.
 * Median has breakdown point of 50% (resists up to 50% outliers).
 *
 * For comparison: mean breakdown point = 0% (single outlier corrupts it).
 *
 * @param readings   Array of sensor readings [n_sensors]
 * @param n_sensors  Number of sensors
 * @return           Median value
 */
double ss_fusion_median(const double *readings, size_t n_sensors);

/**
 * @brief Robust fusion: trim extremes, then average the rest
 *
 * Sorts readings, removes k lowest and k highest, averages the middle.
 * Alpha-trimmed mean with alpha = k/n_sensors.
 *
 * @param readings   Array of sensor readings [n_sensors]
 * @param n_sensors  Number of sensors
 * @param trim_n     Number of extremes to trim from each end
 * @return           Trimmed mean
 */
double ss_fusion_trimmed_mean(const double *readings, size_t n_sensors,
                              size_t trim_n);

/**
 * @brief Exponential moving average (EMA) sensor fusion
 *
 * Recursive low-pass filter for sensor readings:
 *   y_k = alpha * x_k + (1 - alpha) * y_{k-1}
 *
 * Equivalent to a first-order IIR filter with time constant:
 *   tau = -dt / ln(1 - alpha) ≈ dt / alpha  (for small alpha)
 *
 * @param current_ema   Current EMA value (updated in place)
 * @param new_reading   New sensor reading
 * @param alpha         Smoothing factor (0-1, typically 0.05-0.3)
 * @return              Updated EMA value
 */
double ss_fusion_ema(double *current_ema, double new_reading, double alpha);

/**
 * @brief Cumulative moving average for stable sensors
 *
 * Maintains running sum and count for exact average without storing
 * all samples. Suitable for sensor burn-in and initial calibration.
 *
 *   mean_n = mean_{n-1} + (x_n - mean_{n-1}) / n
 *
 * @param sum    Running sum accumulator (updated in place)
 * @param count  Running count accumulator (updated in place)
 * @param value  New sensor reading
 * @return       Current cumulative average
 */
double ss_fusion_cumulative_avg(double *sum, uint64_t *count, double value);

/* =========================================================================
 * L8: Advanced Sensor Fusion — Mahalanobis Distance Outlier Detection
 *
 * Mahalanobis distance: D_M(x) = sqrt((x - mu)^T * S^{-1} * (x - mu))
 *
 * Accounts for correlation in multi-sensor readings. A reading is
 * considered an outlier if D_M exceeds a chi-squared threshold.
 * ========================================================================= */

/**
 * @brief Compute Mahalanobis distance for outlier detection
 *
 * @param reading   Observed sensor reading
 * @param mean      Expected mean of the sensor
 * @param variance  Expected variance of the sensor
 * @return          Mahalanobis distance (1D: |reading - mean| / sigma)
 */
double ss_fusion_mahalanobis_1d(double reading, double mean, double variance);

/**
 * @brief Multi-sensor outlier rejection using Mahalanobis distance
 *
 * For each sensor reading, computes D_M from the fused mean.
 * Readings with D_M > threshold are rejected.
 *
 * @param readings      Array of sensor readings [n_sensors]
 * @param variances     Array of sensor variances [n_sensors]
 * @param n_sensors     Number of sensors
 * @param threshold     Mahalanobis distance threshold (typically 2.5-3.0)
 * @param accepted_mask Output: bitmask of accepted sensors (bit i set = accepted)
 * @param fused_value   Output: fused estimate from accepted readings only
 * @return              Number of accepted readings
 */
int ss_fusion_outlier_rejection(const double *readings,
                                const double *variances,
                                size_t n_sensors,
                                double threshold,
                                uint32_t *accepted_mask,
                                double *fused_value);

#ifdef __cplusplus
}
#endif
#endif /* SENSOR_FUSION_H */
