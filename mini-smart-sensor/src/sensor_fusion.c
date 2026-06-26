/**
 * @file    sensor_fusion.c
 * @brief   Multi-sensor data fusion algorithms
 *
 * Knowledge Coverage:
 *   L5 (Algorithms): 1D Kalman filter, 3-axis IMU Kalman filter,
 *                    complementary filter, weighted average fusion,
 *                    median fusion, trimmed mean, EMA, cumulative avg
 *   L4 (Fundamental Laws): Kalman-Bucy filtering, Bayesian fusion,
 *                          inverse-variance weighting optimality
 *   L6 (Canonical Problems): IMU orientation estimation, multi-sensor
 *                            temperature fusion with outlier rejection
 *   L8 (Advanced): Mahalanobis distance outlier detection
 *
 * Reference:
 *   Kalman, R.E. (1960), Trans. ASME J. Basic Eng., 82:35-45
 *   Maybeck, P.S. (1979), "Stochastic Models, Estimation and Control", Vol. 1
 *   Mahony, R. et al. (2008), IEEE TAC 53(5):1203-1218
 *   Hall, D.L. & Llinas, J. (2001), "Handbook of Multisensor Data Fusion", CRC
 *   Welch, G. & Bishop, G. (2006), "An Introduction to the Kalman Filter", UNC TR 95-041
 *
 * Course Mapping:
 *   MIT 6.450 — Digital Communications (estimation, Kalman filtering)
 *   Stanford EE363 — Linear Dynamical Systems (Kalman filter)
 *   Berkeley EE221A — Linear Systems (state estimation, observers)
 *   Michigan EECS 455 — Communications (stochastic estimation)
 */

#include "smart_sensor.h"
#include "sensor_fusion.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* =========================================================================
 * L5: 1-Dimensional Kalman Filter
 *
 * The Kalman filter is the optimal recursive state estimator for
 * linear Gaussian systems. It minimizes the mean-squared error of
 * the state estimate given noisy measurements.
 *
 * State-space model:
 *   x_k = A * x_{k-1} + w_k      w_k ~ N(0, Q)   [process]
 *   z_k = H * x_k + v_k           v_k ~ N(0, R)   [measurement]
 *
 * The filter recursively produces:
 *   x_hat_k|k:    state estimate given measurements up to time k
 *   P_k|k:        error covariance of the state estimate
 *
 * Optimality: Under Gaussian assumptions, the Kalman filter is the
 * minimum-variance unbiased estimator (MVUE). If noise is non-Gaussian,
 * it is the best linear unbiased estimator (BLUE).
 *
 * Reference: Kalman (1960); Welch & Bishop (2006)
 * ========================================================================= */

/**
 * @brief Initialize 1D Kalman filter for scalar sensor fusion
 *
 * Selects reasonable default parameters for a constant (non-evolving)
 * state model: A=1.0, H=1.0.
 *
 * Tuning guidance:
 *   - P0 large (e.g., 1.0): filter is uncertain about initial state
 *     → rapid convergence toward measurements
 *   - P0 small (e.g., 0.01): filter trusts initial guess more
 *     → slower adaptation
 *   - Q large: model is uncertain → trusts measurements more
 *     → faster response, noisier output
 *   - Q small: model is trusted → smoother output, slower response
 *   - R = sensor measurement variance (estimate from noise analysis)
 *
 * @param kf  Filter state
 * @param x0  Initial state estimate
 * @param P0  Initial error covariance
 * @param Q   Process noise covariance
 * @param R   Measurement noise covariance
 */
void ss_fusion_kalman_1d_init(ss_kalman_1d_t *kf, double x0,
                              double P0, double Q, double R) {
    if (!kf) return;
    memset(kf, 0, sizeof(*kf));
    kf->x = x0;
    kf->P = P0;
    kf->A = 1.0;
    kf->H = 1.0;
    kf->Q = Q;
    kf->R = R;
    kf->iteration = 0;
}

/**
 * @brief Kalman predict step (time update / a priori estimate)
 *
 * Projects the state and covariance forward one time step:
 *   x_k^- = A * x_{k-1}
 *   P_k^- = A * P_{k-1} * A + Q
 *
 * For the default constant-state model (A=1):
 *   x_k^- = x_{k-1}       (state unchanged — best guess is current estimate)
 *   P_k^- = P_{k-1} + Q   (uncertainty grows by process noise)
 *
 * This step must be called before each measurement update.
 *
 * @param kf  Filter state (updated in place)
 */
void ss_fusion_kalman_1d_predict(ss_kalman_1d_t *kf) {
    if (!kf) return;
    kf->x = kf->A * kf->x;
    kf->P = kf->A * kf->P * kf->A + kf->Q;
}

/**
 * @brief Kalman update step (measurement update / a posteriori estimate)
 *
 * Incorporates a new measurement z to correct the predicted state:
 *   K = P_k^- * H / (H * P_k^- * H + R)       [Kalman gain]
 *   x_k = x_k^- + K * (z - H * x_k^-)          [state update]
 *   P_k = (I - K * H) * P_k^-                  [covariance update]
 *
 * The Kalman gain K balances the prediction vs. measurement:
 *   - K → 0 when R is large (noisy measurement → trust prediction)
 *   - K → 1 when P is large (uncertain prediction → trust measurement)
 *
 * @param kf  Filter state (updated in place)
 * @param z   New measurement
 * @return    Updated (a posteriori) state estimate
 */
double ss_fusion_kalman_1d_update(ss_kalman_1d_t *kf, double z) {
    double K, innovation;

    if (!kf) return z;

    /* Kalman gain: K = P^- * H^T / (H * P^- * H^T + R) */
    K = kf->P * kf->H / (kf->H * kf->P * kf->H + kf->R);

    /* Innovation (measurement residual): z - H*x^- */
    innovation = z - kf->H * kf->x;

    /* State update: x = x^- + K * innovation */
    kf->x = kf->x + K * innovation;

    /* Covariance update: P = (I - K*H) * P^- */
    kf->P = (1.0 - K * kf->H) * kf->P;

    kf->iteration++;
    return kf->x;
}

/**
 * @brief Combined predict + update in one call
 *
 * Convenience wrapper for the common case where prediction and
 * measurement occur at the same time step (most sensor fusion
 * applications).
 *
 * @param kf  Filter state
 * @param z   New measurement
 * @return    Updated state estimate
 */
double ss_fusion_kalman_1d_step(ss_kalman_1d_t *kf, double z) {
    ss_fusion_kalman_1d_predict(kf);
    return ss_fusion_kalman_1d_update(kf, z);
}

/* =========================================================================
 * L5: 3-Axis IMU Kalman Filter for Orientation
 *
 * State vector: [roll, pitch]^T  (Euler angles in radians)
 *
 * Process model (gyroscope integration):
 *   roll_k  = roll_{k-1}  + dt * gx  (simplified — full kinematics
 *   pitch_k = pitch_{k-1} + dt * gy   requires Euler rate equations)
 *
 * Measurement model (accelerometer):
 *   roll_meas  = atan2(ay, az)
 *   pitch_meas = atan2(-ax, sqrt(ay^2 + az^2))
 *
 * Note: This is a simplified 2-state model. A full attitude estimator
 * would use quaternions and a 4-state or 7-state EKF. This model is
 * suitable for low-dynamics applications where pitch < 60 degrees
 * (small-angle approximation is valid).
 *
 * Reference: Valenti, R.G. et al. (2015), Sensors 15(8):19302-19330
 * ========================================================================= */

/**
 * @brief Initialize 3-axis IMU Kalman filter
 *
 * @param kf       Filter state
 * @param dt       Sampling period (seconds)
 * @param q_gyro   Gyroscope process noise variance (rad^2/s^2)
 * @param r_accel  Accelerometer measurement noise variance (g^2)
 */
void ss_fusion_kalman_imu_init(ss_kalman_3d_imu_t *kf, double dt,
                               double q_gyro, double r_accel) {
    if (!kf) return;
    memset(kf, 0, sizeof(*kf));

    /* Initial state: level (zero roll, zero pitch) */
    kf->x[0] = 0.0;  /* Roll */
    kf->x[1] = 0.0;  /* Pitch */
    kf->dt = dt;

    /* Initial covariance: moderately uncertain */
    kf->P[0] = 0.1; kf->P[1] = 0.0;
    kf->P[2] = 0.0; kf->P[3] = 0.1;

    /* Process noise (gyroscope): Q = G * q_gyro * G^T * dt^2 */
    {
        double q = q_gyro * dt * dt;
        kf->Q_gyro[0] = q;   kf->Q_gyro[1] = 0.0;
        kf->Q_gyro[2] = 0.0; kf->Q_gyro[3] = q;
    }

    /* Measurement noise (accelerometer) */
    kf->R_accel[0] = r_accel; kf->R_accel[1] = 0.0;
    kf->R_accel[2] = 0.0;     kf->R_accel[3] = r_accel;

    kf->iteration = 0;
}

/**
 * @brief IMU Kalman predict step using gyroscope measurements
 *
 * Simplified Euler angle kinematics (valid for small angles):
 *   roll_pred  = roll + dt * gx
 *   pitch_pred = pitch + dt * gy
 *
 * For a full implementation, the Euler rate equations should be used:
 *   roll_dot  = gx + gy*sin(roll)*tan(pitch) + gz*cos(roll)*tan(pitch)
 *   pitch_dot = gy*cos(roll) - gz*sin(roll)
 *
 * @param kf  Filter state
 * @param gx  Gyro x-axis rate (rad/s)
 * @param gy  Gyro y-axis rate (rad/s)
 * @param gz  Gyro z-axis rate (rad/s)
 */
void ss_fusion_kalman_imu_predict(ss_kalman_3d_imu_t *kf,
                                  double gx, double gy, double gz) {
    (void)gz; /* z-axis gyro not used in simplified 2-state model */

    if (!kf) return;

    /* State prediction: Euler integration */
    kf->x[0] += kf->dt * gx;   /* Roll update */
    kf->x[1] += kf->dt * gy;   /* Pitch update */

    /* Covariance prediction: P = P + Q */
    kf->P[0] += kf->Q_gyro[0];
    kf->P[1] += kf->Q_gyro[1];
    kf->P[2] += kf->Q_gyro[2];
    kf->P[3] += kf->Q_gyro[3];
}

/**
 * @brief IMU Kalman update using accelerometer (gravity vector) measurement
 *
 * Converts accelerometer readings to roll and pitch angles using
 * the gravity vector geometry:
 *   roll_meas  = atan2(ay, az)
 *   pitch_meas = atan2(-ax, sqrt(ay*ay + az*az))
 *
 * The measurement Jacobian H = I (2x2 identity) since we are
 * directly measuring roll and pitch.
 *
 * Note: This measurement model assumes the only acceleration is gravity.
 * During dynamic motion, accelerometer readings include linear acceleration
 * and will produce erroneous angle estimates. The Kalman gain automatically
 * reduces trust in accelerometer when innovation is large.
 *
 * @param kf  Filter state
 * @param ax  Accelerometer x-axis (g)
 * @param ay  Accelerometer y-axis (g)
 * @param az  Accelerometer z-axis (g)
 */
void ss_fusion_kalman_imu_update(ss_kalman_3d_imu_t *kf,
                                 double ax, double ay, double az) {
    double z_roll, z_pitch;       /* Measurements */
    double innov_roll, innov_pitch; /* Innovations */
    double S_roll, S_pitch;       /* Innovation covariances */
    double K_roll, K_pitch;       /* Kalman gains */

    if (!kf) return;

    /* Convert accelerometer to roll and pitch (gravity reference) */
    z_roll  = atan2(ay, az);
    z_pitch = atan2(-ax, sqrt(ay * ay + az * az));

    /* Innovation: z - H*x (H = I, so innovation = z - x) */
    innov_roll  = z_roll  - kf->x[0];
    innov_pitch = z_pitch - kf->x[1];

    /* Innovation covariance: S = H*P*H^T + R = P + R (since H=I) */
    S_roll  = kf->P[0] + kf->R_accel[0];
    S_pitch = kf->P[3] + kf->R_accel[3];

    /* Kalman gains: K = P*H^T * S^{-1} (scalar division since diagonal) */
    if (fabs(S_roll) > 1e-30) K_roll = kf->P[0] / S_roll;
    else K_roll = 0.0;

    if (fabs(S_pitch) > 1e-30) K_pitch = kf->P[3] / S_pitch;
    else K_pitch = 0.0;

    /* State update: x = x + K*innovation */
    kf->x[0] += K_roll  * innov_roll;
    kf->x[1] += K_pitch * innov_pitch;

    /* Covariance update: P = (I - K*H)*P (H=I, diagonal) */
    kf->P[0] = (1.0 - K_roll)  * kf->P[0];
    kf->P[3] = (1.0 - K_pitch) * kf->P[3];
    /* Off-diagonal terms also need updating: P01 = (1-K_roll)*P01 */
    kf->P[1] = (1.0 - K_roll) * kf->P[1];
    kf->P[2] = (1.0 - K_pitch) * kf->P[2];

    kf->iteration++;
}

/**
 * @brief Get current roll and pitch estimates from IMU Kalman filter
 *
 * @param kf     Filter state
 * @param roll   Output: roll angle (rad)
 * @param pitch  Output: pitch angle (rad)
 */
void ss_fusion_kalman_imu_get_angles(const ss_kalman_3d_imu_t *kf,
                                     double *roll, double *pitch) {
    if (kf) {
        if (roll)  *roll  = kf->x[0];
        if (pitch) *pitch = kf->x[1];
    }
}

/* =========================================================================
 * L5: Complementary Filter for IMU Orientation
 *
 * The complementary filter exploits the frequency-domain complementarity
 * of gyroscope and accelerometer measurements:
 *   - Gyroscope: accurate at high frequencies (low drift in short term)
 *   - Accelerometer: accurate at low frequencies (gravity reference)
 *
 * Filter equation:
 *   angle = alpha * (angle + gyro_rate * dt) + (1 - alpha) * accel_angle
 *
 * This is equivalent to:
 *   - High-pass filtering the gyroscope integration (removes drift)
 *   - Low-pass filtering the accelerometer (removes vibration noise)
 *   - Summing the two filtered signals
 *
 * Transfer function analysis:
 *   H_gyro(s)  = alpha * s / (s + (1-alpha)/dt)     [highpass]
 *   H_accel(s) = (1-alpha)/dt / (s + (1-alpha)/dt)  [lowpass]
 *   H_gyro(s) + H_accel(s) = 1                     [all-pass — complementary!]
 *
 * The parameter alpha determines the crossover frequency:
 *   f_crossover = (1-alpha) / (2*pi*dt*alpha)
 *
 * Typical values:
 *   alpha = 0.98, dt = 0.01s → f_crossover ≈ 0.3 Hz
 *   (trust gyro above 0.3 Hz, trust accel below 0.3 Hz)
 *
 * Reference: Mahony et al. (2008), IEEE TAC 53(5):1203-1218
 * ========================================================================= */

void ss_fusion_complementary_init(ss_complementary_filter_t *cf,
                                  double alpha, double dt) {
    if (!cf) return;
    memset(cf, 0, sizeof(*cf));
    cf->roll  = 0.0;
    cf->pitch = 0.0;
    cf->alpha = alpha;
    cf->dt    = dt;
    cf->iteration = 0;
}

/**
 * @brief Complementary filter update
 *
 * Combines gyroscope angular rates and accelerometer gravity vector
 * to estimate roll and pitch angles.
 *
 * Process:
 *   1. Integrate gyroscope: angle_gyro = previous_angle + gyro_rate * dt
 *   2. Compute accel angle from gravity vector
 *   3. Blend: angle_new = alpha * angle_gyro + (1-alpha) * angle_accel
 *
 * @param cf  Filter state
 * @param gx  Gyro x rate (rad/s)
 * @param gy  Gyro y rate (rad/s)
 * @param gz  Gyro z rate (rad/s)
 * @param ax  Accel x (g)
 * @param ay  Accel y (g)
 * @param az  Accel z (g)
 */
void ss_fusion_complementary_update(ss_complementary_filter_t *cf,
                                    double gx, double gy, double gz,
                                    double ax, double ay, double az) {
    double roll_accel, pitch_accel;
    double roll_gyro, pitch_gyro;
    (void)gz;

    if (!cf) return;

    /* Accel-based angles: gravity vector decomposition */
    roll_accel  = atan2(ay, az);
    pitch_accel = atan2(-ax, sqrt(ay * ay + az * az));

    /* Gyro-based angles: integrate rate */
    roll_gyro  = cf->roll  + gx * cf->dt;
    pitch_gyro = cf->pitch + gy * cf->dt;

    /* Complementary blend */
    cf->roll  = cf->alpha * roll_gyro  + (1.0 - cf->alpha) * roll_accel;
    cf->pitch = cf->alpha * pitch_gyro + (1.0 - cf->alpha) * pitch_accel;

    cf->iteration++;
}

/* =========================================================================
 * L5: Weighted Average and Statistical Fusion
 *
 * Inverse-variance weighting is optimal for combining independent
 * measurements with known (or estimated) variances. It produces
 * the minimum-variance unbiased estimator of the true value.
 *
 * Derivation (Lagrange multiplier optimization):
 *   Minimize: Var(theta_hat) = sum(w_i^2 * sigma_i^2)
 *   Subject to: sum(w_i) = 1 (unbiased constraint)
 *   Solution: w_i = (1/sigma_i^2) / sum_j(1/sigma_j^2)
 *
 * The fused variance is always less than or equal to the smallest
 * individual sensor variance:
 *   1/sigma_fused^2 = sum_i(1/sigma_i^2)
 *   → sigma_fused = 1/sqrt(sum_i(1/sigma_i^2))
 *
 * Example: Three sensors with sigma = [2, 3, 5]
 *   sigma_fused = 1/sqrt(1/4 + 1/9 + 1/25) = 1/sqrt(0.4) ≈ 1.58
 *   (better than the best individual sensor at sigma=2)
 * ========================================================================= */

/**
 * @brief Inverse-variance weighted average (optimal fusion for Gaussian sensors)
 *
 * @param readings        Sensor readings [n_sensors]
 * @param variances       Sensor noise variances [n_sensors]
 * @param n_sensors       Number of sensors
 * @param fused_variance  Output: variance of the fused estimate (may be NULL)
 * @return                Optimal weighted average
 */
double ss_fusion_weighted_average(const double *readings,
                                  const double *variances,
                                  size_t n_sensors,
                                  double *fused_variance) {
    size_t i;
    double sum_weighted = 0.0;
    double sum_weights  = 0.0;

    if (!readings || !variances || n_sensors == 0) return 0.0;
    if (n_sensors == 1) {
        if (fused_variance) *fused_variance = variances[0];
        return readings[0];
    }

    for (i = 0; i < n_sensors; i++) {
        if (variances[i] <= 0.0) continue; /* Skip zero-variance sensors? Try small epsilon */
        {
            double w = 1.0 / variances[i];
            sum_weighted += w * readings[i];
            sum_weights  += w;
        }
    }

    if (sum_weights <= 0.0) {
        /* All variances invalid — fall back to simple average */
        double avg = 0.0;
        for (i = 0; i < n_sensors; i++) avg += readings[i];
        avg /= (double)n_sensors;
        if (fused_variance) *fused_variance = 1.0;
        return avg;
    }

    if (fused_variance) {
        *fused_variance = 1.0 / sum_weights;
    }

    return sum_weighted / sum_weights;
}

/**
 * @brief Median: robust estimator of central tendency
 *
 * The median has a breakdown point of 50% — meaning that up to 50%
 * of the data can be arbitrary outliers before the median is
 * corrupted. This makes it ideal for sensor fusion in the presence
 * of intermittent faults.
 *
 * In contrast, the mean has a breakdown point of 0% — a single
 * outlier can arbitrarily corrupt the result.
 *
 * Complexity: O(n log n) due to sorting (uses qsort)
 *
 * @param readings    Sensor readings [n_sensors]
 * @param n_sensors   Number of sensors
 * @return            Median value
 */

/* Comparison function for qsort */
static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

double ss_fusion_median(const double *readings, size_t n_sensors) {
    double *sorted;
    double median_val;

    if (!readings || n_sensors == 0) return 0.0;
    if (n_sensors == 1) return readings[0];

    /* Copy to temp array for sorting */
    sorted = (double *)malloc(n_sensors * sizeof(double));
    if (!sorted) {
        /* Fallback: simple average if allocation fails */
        size_t i;
        double avg = 0.0;
        for (i = 0; i < n_sensors; i++) avg += readings[i];
        return avg / (double)n_sensors;
    }

    memcpy(sorted, readings, n_sensors * sizeof(double));
    qsort(sorted, n_sensors, sizeof(double), cmp_double);

    if (n_sensors % 2 == 1) {
        /* Odd: middle element */
        median_val = sorted[n_sensors / 2];
    } else {
        /* Even: average of two middle elements */
        median_val = (sorted[n_sensors / 2 - 1] + sorted[n_sensors / 2]) / 2.0;
    }

    free(sorted);
    return median_val;
}

/**
 * @brief Trimmed mean: remove extremes, average the rest
 *
 * Sorts readings, discards trim_n smallest and trim_n largest,
 * then computes the mean of the remaining values.
 *
 * The trimming fraction alpha = trim_n / n_sensors. Common choices:
 *   - alpha = 0: ordinary mean (no trimming)
 *   - alpha = 0.1-0.2: mild trimming — good balance for typical data
 *   - alpha = 0.5: approaches the median (when 2*trim_n ≈ n_sensors - 1)
 *
 * @param readings    Array of sensor readings [n_sensors]
 * @param n_sensors   Number of sensors
 * @param trim_n      Number of extremes to trim from each side
 * @return            Trimmed mean
 */
double ss_fusion_trimmed_mean(const double *readings, size_t n_sensors,
                              size_t trim_n) {
    double *sorted;
    size_t i, start, end;
    double sum = 0.0;

    if (!readings || n_sensors == 0) return 0.0;
    if (2 * trim_n >= n_sensors) {
        /* Trimming too much — fall back to median */
        return ss_fusion_median(readings, n_sensors);
    }

    sorted = (double *)malloc(n_sensors * sizeof(double));
    if (!sorted) {
        double avg = 0.0;
        for (i = 0; i < n_sensors; i++) avg += readings[i];
        return avg / (double)n_sensors;
    }

    memcpy(sorted, readings, n_sensors * sizeof(double));
    qsort(sorted, n_sensors, sizeof(double), cmp_double);

    start = trim_n;
    end   = n_sensors - trim_n;

    for (i = start; i < end; i++) {
        sum += sorted[i];
    }

    free(sorted);
    return sum / (double)(end - start);
}

/**
 * @brief Exponential moving average (EMA) — recursive lowpass filter
 *
 * y_k = alpha * x_k + (1 - alpha) * y_{k-1}
 *
 * This is identical to the first-order IIR lowpass filter.
 * In the fusion context, it acts as a temporal fusion method,
 * combining the current reading with past history.
 *
 * The half-life (number of samples to reduce a step to half amplitude):
 *   N_half = ln(0.5) / ln(1 - alpha) ≈ 0.693 / alpha (for small alpha)
 *
 * @param current_ema  Current EMA value (updated in place)
 * @param new_reading  New sensor reading
 * @param alpha        Smoothing factor (0-1)
 * @return             Updated EMA value
 */
double ss_fusion_ema(double *current_ema, double new_reading, double alpha) {
    if (!current_ema) return new_reading;
    if (alpha <= 0.0) return *current_ema;
    if (alpha >= 1.0) {
        *current_ema = new_reading;
        return new_reading;
    }
    *current_ema = alpha * new_reading + (1.0 - alpha) * (*current_ema);
    return *current_ema;
}

/**
 * @brief Cumulative moving average (exact, no storage of past samples)
 *
 * Maintains a running sum and count. Each new reading updates:
 *   sum_new = sum_old + reading
 *   count_new = count_old + 1
 *   mean_new = sum_new / count_new
 *
 * Advantages:
 *   - No storage of past samples needed
 *   - Exact (no filter lag like EMA)
 *   - Each sample has equal weight regardless of age
 *
 * Disadvantages:
 *   - Does not adapt to changing conditions (all samples weighted equally)
 *   - Numerical precision may degrade for very large counts (use double for sum)
 *
 * @param sum    Running sum accumulator (updated in place)
 * @param count  Running count (updated in place)
 * @param value  New sensor reading
 * @return       Current cumulative average
 */
double ss_fusion_cumulative_avg(double *sum, uint64_t *count, double value) {
    if (!sum || !count) return value;

    *sum += value;
    (*count)++;

    if (*count == 0) return 0.0;
    return *sum / (double)(*count);
}

/* =========================================================================
 * L8: Mahalanobis Distance for Outlier Detection
 *
 * The Mahalanobis distance measures how many standard deviations a
 * point is from the mean, accounting for the covariance structure:
 *
 *   D_M(x) = sqrt((x - mu)^T * S^{-1} * (x - mu))
 *
 * In 1D: D_M(x) = |x - mu| / sigma (the absolute z-score)
 *
 * A reading with D_M > threshold (typically 2.5-3.0, corresponding to
 * p < 0.01 under normality) is flagged as an outlier.
 *
 * Reference: Mahalanobis, P.C. (1936), Proc. Nat. Inst. Sci. India 2:49-55
 * ========================================================================= */

double ss_fusion_mahalanobis_1d(double reading, double mean, double variance) {
    if (variance <= 0.0) return (reading != mean) ? INFINITY : 0.0;
    return fabs(reading - mean) / sqrt(variance);
}

/**
 * @brief Multi-sensor outlier rejection via Mahalanobis distance
 *
 * Algorithm:
 *   1. Compute inverse-variance weighted mean of all readings
 *   2. For each sensor: compute D_M from the fused mean
 *   3. Reject readings with D_M > threshold
 *   4. Re-compute fused mean using only accepted readings
 *
 * This implements a simple 1-iteration hard-decision outlier rejection.
 * For more sophisticated approaches, see RANSAC or iteratively
 * reweighted least squares.
 *
 * @param readings       Sensor readings [n_sensors]
 * @param variances      Sensor variances [n_sensors]
 * @param n_sensors      Number of sensors
 * @param threshold       Mahalanobis distance threshold (2.5-3.0 recommended)
 * @param accepted_mask  Output: bitmask of accepted sensors
 * @param fused_value    Output: fused mean from accepted sensors
 * @return               Number of accepted sensors
 */
int ss_fusion_outlier_rejection(const double *readings,
                                const double *variances,
                                size_t n_sensors,
                                double threshold,
                                uint32_t *accepted_mask,
                                double *fused_value) {
    size_t i;
    double init_mean;
    double fused_var;
    uint32_t mask = 0;
    int n_accepted = 0;

    if (!readings || !variances || n_sensors == 0) {
        if (accepted_mask) *accepted_mask = 0;
        if (fused_value) *fused_value = 0.0;
        return 0;
    }

    /* Step 1: Compute initial fused mean (all sensors) */
    init_mean = ss_fusion_weighted_average(readings, variances, n_sensors,
                                           &fused_var);

    /* Step 2: Test each sensor */
    for (i = 0; i < n_sensors; i++) {
        double d = ss_fusion_mahalanobis_1d(readings[i], init_mean,
                                            variances[i] + fused_var);
        if (d <= threshold) {
            mask |= (1u << i);
            n_accepted++;
        }
    }

    /* Step 3: Re-fuse using only accepted sensors */
    {
        double accepted_readings[SS_MAX_CHANNELS];
        double accepted_variances[SS_MAX_CHANNELS];
        int j = 0;

        for (i = 0; i < n_sensors && j < SS_MAX_CHANNELS; i++) {
            if (mask & (1u << i)) {
                accepted_readings[j] = readings[i];
                accepted_variances[j] = variances[i];
                j++;
            }
        }

        if (j > 0) {
            *fused_value = ss_fusion_weighted_average(accepted_readings,
                                                      accepted_variances,
                                                      (size_t)j, NULL);
        } else {
            /* No sensor accepted — fall back to all-sensor mean */
            *fused_value = init_mean;
            n_accepted = 1; /* Force at least one accepted for safety */
            mask = 1u;      /* Accept first sensor */
        }
    }

    if (accepted_mask) *accepted_mask = mask;
    return n_accepted;
}
