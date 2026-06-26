/**
 * @file mems_fusion.h
 * @brief L5-L6 MEMS Sensor Fusion — Complementary filter, Mahony, Madgwick,
 *        Kalman filter for IMU orientation estimation
 *
 * Covers:
 *   L5 - Complementary filter, Mahony explicit complementary filter,
 *        Madgwick gradient descent AHRS, Kalman filter for IMU
 *   L6 - Attitude and heading reference system (AHRS), roll/pitch from
 *        accelerometer, yaw from magnetometer, quaternion kinematics,
 *        dead reckoning
 *
 * Reference: Mahony et al., "Nonlinear Complementary Filters on SO(3)" (2008)
 *            Madgwick, "An Efficient Orientation Filter for IMU and MARG" (2010)
 *            Markley & Crassidis, "Fundamentals of Spacecraft Attitude" (2014)
 *            Crassidis, Markley, Cheng, "Nonlinear Attitude Filtering Methods" (2007)
 * Courses: MIT 16.485 (GNSS/Inertial), Stanford EE267 (VR/AR),
 *          ETH 227-0436, Tsinghua 导航原理
 */

#ifndef MEMS_FUSION_H
#define MEMS_FUSION_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include "mems_basics.h"

/* ─── L6: Attitude Representation ──────────────────────────────────────────── */

/**
 * @brief Euler angles (intrinsic ZYX: yaw-pitch-roll)
 *        Roll φ: rotation about x-axis
 *        Pitch θ: rotation about y-axis
 *        Yaw ψ: rotation about z-axis
 */
typedef struct {
    double roll;    /* [rad] */
    double pitch;   /* [rad] */
    double yaw;     /* [rad] */
} MemsEuler;

/* ─── L6: Quaternion Operations ─────────────────────────────────────────────── */

/**
 * @brief Quaternion identity: q = (1, 0, 0, 0)
 */
MemsQuaternion mems_quat_identity(void);

/**
 * @brief Quaternion from axis-angle: q = (cos(θ/2), sin(θ/2)·u)
 */
MemsQuaternion mems_quat_from_axis_angle(double angle_rad,
                                          double ux, double uy, double uz);

/**
 * @brief Quaternion multiplication: q_out = q1 ⊗ q2
 *        Hamilton product. Reference: Markley & Crassidis (2014) §2.2.1
 */
MemsQuaternion mems_quat_multiply(MemsQuaternion q1, MemsQuaternion q2);

/**
 * @brief Quaternion conjugate: q* = (w, -x, -y, -z)
 */
MemsQuaternion mems_quat_conjugate(MemsQuaternion q);

/**
 * @brief Quaternion normalization
 */
MemsQuaternion mems_quat_normalize(MemsQuaternion q);

/**
 * @brief Rotate vector by quaternion: v' = q ⊗ v ⊗ q*
 *        where v is treated as pure quaternion (0, vx, vy, vz)
 */
MemsVector3 mems_quat_rotate(MemsQuaternion q, MemsVector3 v);

/**
 * @brief Convert quaternion to Euler angles (ZYX intrinsic)
 *        Reference: Markley & Crassidis (2014) §2.4
 */
MemsEuler mems_quat_to_euler(MemsQuaternion q);

/**
 * @brief Convert Euler angles to quaternion
 */
MemsQuaternion mems_euler_to_quat(MemsEuler e);

/* ─── L5: Complementary Filter ─────────────────────────────────────────────── */

/**
 * @brief Complementary filter parameters
 */
typedef struct {
    double alpha;         /* filter weight [0,1] — higher = trust gyro more */
    double dt;            /* time step [s] */
    MemsQuaternion q_est;  /* current orientation estimate */
} MemsCompFilter;

/**
 * @brief Initialize complementary filter
 * @param alpha Filter coefficient (0 = trust accel, 1 = trust gyro)
 * @param dt    Sample period [s]
 */
MemsCompFilter mems_comp_filter_init(double alpha, double dt);

/**
 * @brief Update complementary filter with gyro + accel measurements
 *        q_gyro = q_est + (dt/2) · q_est ⊗ ω
 *        q_accel = from gravity vector
 *        q_est = α · q_gyro + (1-α) · q_accel  (then normalize)
 *        Reference: Euston et al., "A Complementary Filter for Attitude" (2008)
 *
 * @param cf    Filter state
 * @param gyro  Gyro measurement [rad/s] including Earth rate
 * @param accel Accelerometer measurement [m/s²] or [g]
 * @return Updated orientation quaternion
 */
MemsQuaternion mems_comp_filter_update(MemsCompFilter *cf,
                                        MemsVector3 gyro,
                                        MemsVector3 accel);

/* ─── L5: Mahony Explicit Complementary Filter ─────────────────────────────── */

/**
 * @brief Mahony filter parameters (SO(3) complementary filter)
 *        Reference: Mahony et al., IEEE Trans. Auto. Control (2008)
 */
typedef struct {
    double kp;            /* proportional gain */
    double ki;            /* integral gain */
    double dt;            /* time step [s] */
    MemsQuaternion q_est;  /* current estimate */
    MemsVector3 bias_int;  /* integral of gyro bias error */
    double beta;          /* alternative: explicit complementary gain */
} MemsMahonyFilter;

/**
 * @brief Initialize Mahony filter
 */
MemsMahonyFilter mems_mahony_init(double kp, double ki, double dt);

/**
 * @brief Mahony filter update step (AHRS with accel + mag)
 *        ω_corrected = ω_gyro + ω_correction
 *        where ω_correction = Kp·e + Ki·∫e dt
 *        and e = v_accel × v_estimated or similar
 *
 * @param mf    Filter state
 * @param gyro  Gyro [rad/s]
 * @param accel Accelerometer [m/s²]
 * @param mag   Magnetometer [μT] (optional, pass zero vector if unused)
 */
MemsQuaternion mems_mahony_update(MemsMahonyFilter *mf,
                                   MemsVector3 gyro,
                                   MemsVector3 accel,
                                   MemsVector3 mag);

/* ─── L5: Madgwick Gradient Descent AHRS ───────────────────────────────────── */

/**
 * @brief Madgwick filter parameters
 *        Reference: Madgwick, "An Efficient Orientation Filter" (2010)
 *                   Madgwick, Harrison, Vaidyanathan (2011)
 */
typedef struct {
    double beta;           /* algorithm gain (tuned to gyro noise) */
    double dt;             /* time step [s] */
    MemsQuaternion q_est;   /* current estimate */
    double gyro_trust;     /* effective trust in gyro (derived from beta) */
} MemsMadgwickFilter;

/**
 * @brief Initialize Madgwick filter
 * @param beta Algorithm gain — higher values increase accel/mag correction
 * @param dt   Sample period [s]
 */
MemsMadgwickFilter mems_madgwick_init(double beta, double dt);

/**
 * @brief Madgwick filter update (6-axis: gyro + accel, no mag)
 *        q_est = q_gyro - β · ∇f/||∇f||
 *        where ∇f is the gradient of the orientation error function.
 *        Complexity: O(1), ~40 scalar multiplications per update
 *
 * @param mf    Filter state
 * @param gyro  Gyro [rad/s]
 * @param accel Accelerometer [m/s²]
 */
MemsQuaternion mems_madgwick_update_imu(MemsMadgwickFilter *mf,
                                         MemsVector3 gyro,
                                         MemsVector3 accel);

/**
 * @brief Madgwick filter update (9-axis: gyro + accel + mag)
 *        Full MARG (Magnetic, Angular Rate, Gravity) AHRS
 *
 * @param mf    Filter state
 * @param gyro  Gyro [rad/s]
 * @param accel Accelerometer [m/s²]
 * @param mag   Magnetometer [μT]
 */
MemsQuaternion mems_madgwick_update_marg(MemsMadgwickFilter *mf,
                                          MemsVector3 gyro,
                                          MemsVector3 accel,
                                          MemsVector3 mag);

/* ─── L5: 1D Kalman Filter for Sensor Fusion ───────────────────────────────── */

/**
 * @brief 1-D Kalman filter for fusing two noisy measurements
 *        (e.g., gyro integration + accelerometer tilt)
 *
 *        State: angle θ
 *        Predict: θ̂ₖ = θ̂ₖ₋₁ + ω·dt, P = P + Q
 *        Update:  K = P / (P + R), θ̂ += K·(z - θ̂), P *= (1 - K)
 */
typedef struct {
    double state;         /* current state estimate [rad] */
    double covariance;    /* estimation error covariance */
    double process_noise; /* Q: process noise spectral density */
    double meas_noise;    /* R: measurement noise variance */
    double dt;            /* time step [s] */
} MemsKalman1D;

/**
 * @brief Initialize 1D Kalman filter
 */
MemsKalman1D mems_kalman1d_init(double initial_state, double initial_cov,
                                 double process_noise, double meas_noise,
                                 double dt);

/**
 * @brief Kalman predict step: integrate gyro rate
 */
void mems_kalman1d_predict(MemsKalman1D *kf, double rate);

/**
 * @brief Kalman update step: fuse measurement
 */
void mems_kalman1d_update(MemsKalman1D *kf, double measurement);

/* ─── L6: Roll/Pitch from Accelerometer ────────────────────────────────────── */

/**
 * @brief Compute roll and pitch from accelerometer in static conditions
 *        roll  = atan2(a_y, a_z)
 *        pitch = atan2(-a_x, sqrt(a_y² + a_z²))
 *        Reference: Pedley (2013), "Tilt Sensing Using Accelerometers"
 *
 * @param accel Accelerometer reading [m/s²] in body frame
 * @return Roll and pitch in radians
 */
MemsEuler mems_accel_to_roll_pitch(MemsVector3 accel);

/* ─── L6: Yaw from Magnetometer (with tilt compensation) ─────────────────── */

/**
 * @brief Compute yaw (heading) from magnetometer with tilt compensation
 *        Reference: ST Microsystems AN3192, "Tilt-Compensated eCompass"
 *
 * @param mag    Magnetometer reading [μT] in body frame
 * @param roll   Roll angle [rad]
 * @param pitch  Pitch angle [rad]
 * @param declination Magnetic declination [rad] for the location
 * @return Yaw angle [rad], 0 = magnetic north
 */
double mems_mag_to_yaw(MemsVector3 mag, double roll, double pitch,
                        double declination);

/* ─── L6: Dead Reckoning ───────────────────────────────────────────────────── */

/**
 * @brief Dead reckoning with IMU data
 *        Integrate acceleration to velocity, velocity to position.
 *        Rotation compensated by orientation quaternion.
 *        Reference: Groves (2013) §5.4
 *
 * @param accel_body   Acceleration in body frame [m/s²]
 * @param orientation  Current orientation quaternion
 * @param dt           Time step [s]
 * @param position     Position [m], updated in-place
 * @param velocity     Velocity [m/s], updated in-place
 * @param gravity      Gravity vector in world frame [m/s²] (typ. (0,0,9.81))
 */
void mems_dead_reckon(MemsVector3 accel_body, MemsQuaternion orientation,
                       double dt, MemsVector3 *position, MemsVector3 *velocity,
                       MemsVector3 gravity);

#endif /* MEMS_FUSION_H */
