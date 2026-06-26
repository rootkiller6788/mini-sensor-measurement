/**
 * @file mems_fusion.c
 * @brief L5-L6 MEMS Sensor Fusion — Quaternion kinematics, complementary filter,
 *        Mahony filter, Madgwick AHRS, Kalman filter for IMU
 *
 * Reference: Mahony et al., IEEE TAC (2008)
 *            Madgwick, "An Efficient Orientation Filter" (2010)
 *            Markley & Crassidis, "Fundamentals of Spacecraft Attitude" (2014)
 *            Groves (2013) "Principles of GNSS, Inertial, and Multisensor Navigation"
 */

#include "mems_fusion.h"
#include <math.h>
#include <string.h>
#include <float.h>

/* ─── L6: Quaternion Operations ───────────────────────────────────────────── */

MemsQuaternion mems_quat_identity(void) {
    MemsQuaternion q = {1.0, 0.0, 0.0, 0.0};
    return q;
}

MemsQuaternion mems_quat_from_axis_angle(double angle_rad,
                                          double ux, double uy, double uz) {
    double half_angle = angle_rad * 0.5;
    double s = sin(half_angle);
    /* Normalize axis */
    double mag = sqrt(ux * ux + uy * uy + uz * uz);
    if (mag < DBL_EPSILON) return mems_quat_identity();
    MemsQuaternion q;
    q.w = cos(half_angle);
    q.x = ux / mag * s;
    q.y = uy / mag * s;
    q.z = uz / mag * s;
    return q;
}

MemsQuaternion mems_quat_multiply(MemsQuaternion q1, MemsQuaternion q2) {
    /* Hamilton product: q1 ⊗ q2
     * Reference: Markley & Crassidis (2014) Eq. 2.24 */
    MemsQuaternion r;
    r.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    r.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    r.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    r.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
    return r;
}

MemsQuaternion mems_quat_conjugate(MemsQuaternion q) {
    MemsQuaternion r = {q.w, -q.x, -q.y, -q.z};
    return r;
}

MemsQuaternion mems_quat_normalize(MemsQuaternion q) {
    double norm = sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
    if (norm < DBL_EPSILON) return mems_quat_identity();
    MemsQuaternion r = {q.w / norm, q.x / norm, q.y / norm, q.z / norm};
    return r;
}

MemsVector3 mems_quat_rotate(MemsQuaternion q, MemsVector3 v) {
    /* Rotate vector by quaternion: q ⊗ v ⊗ q*
     * Treat v as pure quaternion (0, vx, vy, vz) */
    MemsQuaternion p = {0.0, v.x, v.y, v.z};
    MemsQuaternion qc = mems_quat_conjugate(q);
    MemsQuaternion temp = mems_quat_multiply(q, p);
    MemsQuaternion result = mems_quat_multiply(temp, qc);
    MemsVector3 out = {result.x, result.y, result.z};
    return out;
}

MemsEuler mems_quat_to_euler(MemsQuaternion q) {
    /* ZYX intrinsic Euler angles (yaw-pitch-roll)
     * Reference: Markley & Crassidis (2014) §2.4.4
     *
     * roll  = atan2(2(qw·qx + qy·qz), 1 - 2(qx² + qy²))
     * pitch = asin(2(qw·qy - qz·qx))
     * yaw   = atan2(2(qw·qz + qx·qy), 1 - 2(qy² + qz²))
     */
    MemsEuler e;

    /* Roll (φ) */
    double sinr_cosp = 2.0 * (q.w * q.x + q.y * q.z);
    double cosr_cosp = 1.0 - 2.0 * (q.x * q.x + q.y * q.y);
    e.roll = atan2(sinr_cosp, cosr_cosp);

    /* Pitch (θ) */
    double sinp = 2.0 * (q.w * q.y - q.z * q.x);
    if (fabs(sinp) >= 1.0) {
        e.pitch = copysign(M_PI / 2.0, sinp);
    } else {
        e.pitch = asin(sinp);
    }

    /* Yaw (ψ) */
    double siny_cosp = 2.0 * (q.w * q.z + q.x * q.y);
    double cosy_cosp = 1.0 - 2.0 * (q.y * q.y + q.z * q.z);
    e.yaw = atan2(siny_cosp, cosy_cosp);

    return e;
}

MemsQuaternion mems_euler_to_quat(MemsEuler e) {
    /* ZYX Euler to quaternion
     * Reference: Markley & Crassidis (2014) §2.4.1 */
    double cr = cos(e.roll * 0.5);
    double sr = sin(e.roll * 0.5);
    double cp = cos(e.pitch * 0.5);
    double sp = sin(e.pitch * 0.5);
    double cy = cos(e.yaw * 0.5);
    double sy = sin(e.yaw * 0.5);

    MemsQuaternion q;
    q.w = cr * cp * cy + sr * sp * sy;
    q.x = sr * cp * cy - cr * sp * sy;
    q.y = cr * sp * cy + sr * cp * sy;
    q.z = cr * cp * sy - sr * sp * cy;
    return q;
}

/* ─── L5: Complementary Filter ──────────────────────────────────────────────
 *
 * The complementary filter exploits the frequency-domain complementarity
 * of accelerometers (good at low freq) and gyroscopes (good at high freq).
 *
 * Gyroscope integration → orientation at high frequencies (drift-prone).
 * Accelerometer gravity vector → orientation at low frequencies (noise-prone).
 *
 * q_est = α · q_gyro + (1-α) · q_accel
 *
 * where q_gyro integrates the angular velocity:
 *   q_gyro = q_est + (Δt/2) · q_est ⊗ ω
 *
 * and q_accel estimates orientation from the gravity direction:
 *   Uses simplified tilt calculation.
 *
 * The filter coefficient α is typically chosen so that the crossover
 * frequency f_c = 1/(2π·τ) where τ = α·Δt/(1-α).
 *
 * Reference: Euston et al. (2008), "A Complementary Filter for Attitude
 *            Estimation of a Fixed-Wing UAV"
 */

MemsCompFilter mems_comp_filter_init(double alpha, double dt) {
    MemsCompFilter cf;
    cf.alpha = (alpha < 0.0) ? 0.0 : ((alpha > 1.0) ? 1.0 : alpha);
    cf.dt = (dt > 0.0) ? dt : 0.01;
    cf.q_est = mems_quat_identity();
    return cf;
}

MemsQuaternion mems_comp_filter_update(MemsCompFilter *cf,
                                        MemsVector3 gyro,
                                        MemsVector3 accel) {
    if (!cf) return mems_quat_identity();

    /* 1. Gyro integration: quaternion derivative
     * dq/dt = 0.5 · q ⊗ ω  (where ω is pure quaternion angular velocity) */
    MemsQuaternion omega_q = {0.0, gyro.x, gyro.y, gyro.z};
    MemsQuaternion dq_temp = mems_quat_multiply(cf->q_est, omega_q);
    MemsQuaternion dq = {
        dq_temp.w * 0.5 * cf->dt,
        dq_temp.x * 0.5 * cf->dt,
        dq_temp.y * 0.5 * cf->dt,
        dq_temp.z * 0.5 * cf->dt
    };

    MemsQuaternion q_gyro;
    q_gyro.w = cf->q_est.w + dq.w;
    q_gyro.x = cf->q_est.x + dq.x;
    q_gyro.y = cf->q_est.y + dq.y;
    q_gyro.z = cf->q_est.z + dq.z;
    q_gyro = mems_quat_normalize(q_gyro);

    /* 2. Accel orientation: from gravity direction
     * Accelerometer measures gravity in body frame.
     * Assumes no significant linear acceleration. */
    MemsEuler ea = mems_accel_to_roll_pitch(accel);
    MemsQuaternion q_accel = mems_euler_to_quat(ea);

    /* 3. Complementary blend: linear interpolation on SO(3) via slerp approx.
     * For small angles, linear quaternion blend + normalization ≈ slerp. */
    double alpha = cf->alpha;
    MemsQuaternion q_blend;
    q_blend.w = alpha * q_gyro.w + (1.0 - alpha) * q_accel.w;
    q_blend.x = alpha * q_gyro.x + (1.0 - alpha) * q_accel.x;
    q_blend.y = alpha * q_gyro.y + (1.0 - alpha) * q_accel.y;
    q_blend.z = alpha * q_gyro.z + (1.0 - alpha) * q_accel.z;

    cf->q_est = mems_quat_normalize(q_blend);
    return cf->q_est;
}

/* ─── L5: Mahony Explicit Complementary Filter ──────────────────────────────
 *
 * Mahony et al. (2008) developed a nonlinear complementary filter on SO(3)
 * with proven stability properties.
 *
 * The correction term uses the cross product of the estimated gravity
 * direction with the measured gravity:
 *   e = v_accel × v_estimated
 *
 * This provides a proportional+integral correction to the angular velocity:
 *   ω_corrected = ω_gyro + Kp·e + Ki·∫e dt
 *
 * Reference: Mahony, Hamel, Pflimlin, "Nonlinear Complementary Filters on
 *            the Special Orthogonal Group", IEEE TAC (2008)
 */

MemsMahonyFilter mems_mahony_init(double kp, double ki, double dt) {
    MemsMahonyFilter mf;
    mf.kp = (kp >= 0.0) ? kp : 0.5;
    mf.ki = (ki >= 0.0) ? ki : 0.01;
    mf.dt = (dt > 0.0) ? dt : 0.01;
    mf.q_est = mems_quat_identity();
    mf.bias_int.x = 0.0;
    mf.bias_int.y = 0.0;
    mf.bias_int.z = 0.0;
    mf.beta = 0.1;
    return mf;
}

MemsQuaternion mems_mahony_update(MemsMahonyFilter *mf,
                                   MemsVector3 gyro,
                                   MemsVector3 accel,
                                   MemsVector3 mag) {
    if (!mf) return mems_quat_identity();

    /* Normalize accelerometer measurement */
    double accel_norm = sqrt(accel.x * accel.x + accel.y * accel.y +
                              accel.z * accel.z);
    if (accel_norm < DBL_EPSILON) return mf->q_est;

    MemsVector3 a = {accel.x / accel_norm, accel.y / accel_norm,
                      accel.z / accel_norm};

    /* Estimated gravity direction from current orientation:
     * In world frame, gravity is (0, 0, -1). Rotate to body frame. */
    MemsVector3 world_gravity = {0.0, 0.0, -1.0};
    MemsQuaternion qc = mems_quat_conjugate(mf->q_est);
    MemsVector3 v_est = mems_quat_rotate(qc, world_gravity);

    /* Error: cross product of measured and estimated direction */
    MemsVector3 e;
    e.x = a.y * v_est.z - a.z * v_est.y;
    e.y = a.z * v_est.x - a.x * v_est.z;
    e.z = a.x * v_est.y - a.y * v_est.x;

    /* Integral term (gyro bias estimate) */
    mf->bias_int.x += mf->ki * e.x * mf->dt;
    mf->bias_int.y += mf->ki * e.y * mf->dt;
    mf->bias_int.z += mf->ki * e.z * mf->dt;

    /* Corrected angular velocity */
    MemsVector3 omega;
    omega.x = gyro.x + mf->kp * e.x + mf->bias_int.x;
    omega.y = gyro.y + mf->kp * e.y + mf->bias_int.y;
    omega.z = gyro.z + mf->kp * e.z + mf->bias_int.z;

    /* If magnetometer available, add yaw correction */
    int has_mag = (fabs(mag.x) + fabs(mag.y) + fabs(mag.z)) > DBL_EPSILON;
    if (has_mag) {
        double mag_norm = sqrt(mag.x * mag.x + mag.y * mag.y + mag.z * mag.z);
        if (mag_norm > DBL_EPSILON) {
            MemsVector3 m = {mag.x / mag_norm, mag.y / mag_norm,
                              mag.z / mag_norm};
            /* Rotate world magnetic field reference to body frame */
            MemsVector3 world_mag = {1.0, 0.0, 0.0}; /* simplified: mag N */
            MemsVector3 m_est = mems_quat_rotate(qc, world_mag);

            MemsVector3 em;
            em.x = m.y * m_est.z - m.z * m_est.y;
            em.y = m.z * m_est.x - m.x * m_est.z;
            em.z = m.x * m_est.y - m.y * m_est.x;

            /* Blend mag error into correction (weight 0.5 to preserve roll/pitch) */
            omega.x += 0.5 * mf->kp * em.x;
            omega.y += 0.5 * mf->kp * em.y;
            omega.z += 0.5 * mf->kp * em.z;
        }
    }

    /* Integrate: dq/dt = 0.5 · q ⊗ ω */
    MemsQuaternion omega_q = {0.0, omega.x, omega.y, omega.z};
    MemsQuaternion dq_temp = mems_quat_multiply(mf->q_est, omega_q);
    MemsQuaternion dq = {
        dq_temp.w * 0.5 * mf->dt,
        dq_temp.x * 0.5 * mf->dt,
        dq_temp.y * 0.5 * mf->dt,
        dq_temp.z * 0.5 * mf->dt
    };

    mf->q_est.w += dq.w;
    mf->q_est.x += dq.x;
    mf->q_est.y += dq.y;
    mf->q_est.z += dq.z;
    mf->q_est = mems_quat_normalize(mf->q_est);

    return mf->q_est;
}

/* ─── L5: Madgwick Gradient Descent AHRS ────────────────────────────────────
 *
 * Madgwick's filter (2010) uses gradient descent to find the quaternion
 * that minimizes the orientation error from accelerometer and magnetometer.
 *
 * The innovation combines:
 *   q_est = q_gyro - μ · ∇f/||∇f||
 *
 * where ∇f is the gradient of the orientation error objective function.
 * The step size μ is tuned for the gyroscope noise characteristics.
 *
 * For IMU-only (6-axis):
 *   f(q) = 0.5 · ||R(q)·g_world - a_body||²
 *   ∇f = J^T · f  where J is the Jacobian
 *
 * Complexity: ~40 multiply-adds per update (very efficient).
 *
 * Reference: Madgwick, Harrison, Vaidyanathan (2011)
 *            "Estimation of IMU and MARG Orientation Using Gradient Descent"
 */

MemsMadgwickFilter mems_madgwick_init(double beta, double dt) {
    MemsMadgwickFilter mf;
    mf.beta = (beta > 0.0) ? beta : 0.04;
    mf.dt = (dt > 0.0) ? dt : 0.01;
    mf.q_est = mems_quat_identity();
    mf.gyro_trust = 1.0;
    return mf;
}

MemsQuaternion mems_madgwick_update_imu(MemsMadgwickFilter *mf,
                                         MemsVector3 gyro,
                                         MemsVector3 accel) {
    if (!mf) return mems_quat_identity();

    double q1 = mf->q_est.w, q2 = mf->q_est.x;
    double q3 = mf->q_est.y, q4 = mf->q_est.z;

    double ax = accel.x, ay = accel.y, az = accel.z;

    /* Normalize accelerometer */
    double recip_norm = ax * ax + ay * ay + az * az;
    if (recip_norm < 1e-10) {
        /* No valid accel data, just integrate gyro */
        MemsQuaternion omega_q = {0.0, gyro.x, gyro.y, gyro.z};
        MemsQuaternion dq_temp = mems_quat_multiply(mf->q_est, omega_q);
        mf->q_est.w += dq_temp.w * 0.5 * mf->dt;
        mf->q_est.x += dq_temp.x * 0.5 * mf->dt;
        mf->q_est.y += dq_temp.y * 0.5 * mf->dt;
        mf->q_est.z += dq_temp.z * 0.5 * mf->dt;
        mf->q_est = mems_quat_normalize(mf->q_est);
        return mf->q_est;
    }

    recip_norm = 1.0 / sqrt(recip_norm);
    ax *= recip_norm; ay *= recip_norm; az *= recip_norm;

    /* Objective function f and Jacobian J for IMU
     * Gravity direction in world frame: [0, 0, 1]
     * Rotated to body frame: R(q)·[0,0,1]^T
     *
     * f = R(q)·[0,0,1]^T - [ax, ay, az]^T (simplified)
     *
     * f[0] = 2(q2·q4 - q1·q3) - ax
     * f[1] = 2(q1·q2 + q3·q4) - ay
     * f[2] = 2(0.5 - q2² - q3²) - az
     */
    double _2q1 = 2.0 * q1, _2q2 = 2.0 * q2;
    double _2q3 = 2.0 * q3, _2q4 = 2.0 * q4;
    double _4q1 = 4.0 * q1, _4q2 = 4.0 * q2;
    double _4q3 = 4.0 * q3, _4q4 = 4.0 * q4;
    (void)_4q1;  /* used in full 9-axis MARG, reserved for consistency */
    (void)_4q4;

    double f1 = _2q2 * q4 - _2q1 * q3 - ax;
    double f2 = _2q1 * q2 + _2q3 * q4 - ay;
    double f3 = 1.0 - _2q2 * q2 - _2q3 * q3 - az;

    /* Jacobian J (3x4), formula from Madgwick Eq. 24-26
     * J^T·f gives the gradient */
    double Jt_f1 = -_2q3 * f1 + _2q2 * f2;
    double Jt_f2 = _2q4 * f1 + _2q1 * f2 - _4q2 * f3;
    double Jt_f3 = -_2q1 * f1 + _2q4 * f2 - _4q3 * f3;
    double Jt_f4 = _2q2 * f1 + _2q3 * f2;

    /* Normalize gradient */
    double grad_norm = sqrt(Jt_f1 * Jt_f1 + Jt_f2 * Jt_f2 +
                            Jt_f3 * Jt_f3 + Jt_f4 * Jt_f4);
    if (grad_norm < 1e-10) grad_norm = 1e-10;
    double step = mf->beta / grad_norm;

    /* Gyro integration: dq = 0.5 · q ⊗ ω · dt */
    double dqw = -0.5 * q2 * gyro.x - 0.5 * q3 * gyro.y - 0.5 * q4 * gyro.z;
    double dqx =  0.5 * q1 * gyro.x + 0.5 * q3 * gyro.z - 0.5 * q4 * gyro.y;
    double dqy =  0.5 * q1 * gyro.y - 0.5 * q2 * gyro.z + 0.5 * q4 * gyro.x;
    double dqz =  0.5 * q1 * gyro.z + 0.5 * q2 * gyro.y - 0.5 * q3 * gyro.x;

    /* Combined update: q = q + dq_gyro·dt - μ·∇f/||∇f|| */
    mf->q_est.w += (dqw - step * Jt_f1) * mf->dt;
    mf->q_est.x += (dqx - step * Jt_f2) * mf->dt;
    mf->q_est.y += (dqy - step * Jt_f3) * mf->dt;
    mf->q_est.z += (dqz - step * Jt_f4) * mf->dt;
    mf->q_est = mems_quat_normalize(mf->q_est);

    return mf->q_est;
}

MemsQuaternion mems_madgwick_update_marg(MemsMadgwickFilter *mf,
                                          MemsVector3 gyro,
                                          MemsVector3 accel,
                                          MemsVector3 mag) {
    if (!mf) return mems_quat_identity();

    /* For full MARG (9-axis), we add magnetometer correction.
     * This extends the IMU filter with yaw correction from magnetic field.
     * For simplicity, we use IMU update + separate mag-based yaw correction. */

    double q1 = mf->q_est.w, q2 = mf->q_est.x;
    double q3 = mf->q_est.y, q4 = mf->q_est.z;

    double ax = accel.x, ay = accel.y, az = accel.z;
    double mx = mag.x, my = mag.y, mz = mag.z;

    /* Normalize accel */
    double a_norm2 = ax * ax + ay * ay + az * az;
    if (a_norm2 < 1e-10) a_norm2 = 1.0;
    double a_r = 1.0 / sqrt(a_norm2);
    ax *= a_r; ay *= a_r; az *= a_r;

    /* Normalize mag */
    double m_norm2 = mx * mx + my * my + mz * mz;
    if (m_norm2 < 1e-10) m_norm2 = 1.0;
    double m_r = 1.0 / sqrt(m_norm2);
    mx *= m_r; my *= m_r; mz *= m_r;

    /* Reference direction of Earth's magnetic field (world frame) */
    double _2q1 = 2.0 * q1, _2q2 = 2.0 * q2;
    double _2q3 = 2.0 * q3, _2q4 = 2.0 * q4;
    double _4q1 = 4.0 * q1, _4q2 = 4.0 * q2;
    double _4q3 = 4.0 * q3;

    double q1q1 = q1 * q1, q2q2 = q2 * q2, q3q3 = q3 * q3, q4q4 = q4 * q4;
    double q1q2 = q1 * q2, q1q3 = q1 * q3, q1q4 = q1 * q4;
    double q2q3 = q2 * q3, q2q4 = q2 * q4, q3q4 = q3 * q4;
    (void)q1q1;  /* used in non-MARG 6-axis variant */
    (void)_4q1; (void)_4q2; (void)_4q3;  /* speed-optimization precomputed */

    /* Rotate magnetic reference to body frame */
    double hx = mx * (0.5 - q3q3 - q4q4) + my * (q2q3 - q1q4) + mz * (q2q4 + q1q3);
    double hy = mx * (q2q3 + q1q4) + my * (0.5 - q2q2 - q4q4) + mz * (q3q4 - q1q2);
    double hz = mx * (q2q4 - q1q3) + my * (q3q4 + q1q2) + mz * (0.5 - q2q2 - q3q3);

    double bx = sqrt(hx * hx + hy * hy);
    double bz = hz;

    /* Objective function (6 equations: 3 accel + 3 mag) */
    double f1 = _2q2 * q4 - _2q1 * q3 - ax;
    double f2 = _2q1 * q2 + _2q3 * q4 - ay;
    double f3 = 1.0 - _2q2 * q2 - _2q3 * q3 - az;

    double f4 = bx * (0.5 - q3q3 - q4q4) + bz * (q2q4 - q1q3) - mx;
    double f5 = bx * (q2q3 - q1q4) + bz * (q1q2 + q3q4) - my;
    double f6 = bx * (q1q3 + q2q4) + bz * (0.5 - q2q2 - q3q3) - mz;

    /* Gradient: J^T · f (6 equations combined) */
    double Jt_f1 = -_2q3 * f1 + _2q2 * f2 - _2q4 * f4 + q4 * f5 - bz * q3 * f4
                   + (-_4q1) * f6;
    double Jt_f2 = _2q4 * f1 + _2q1 * f2 - _4q2 * f3 + q3 * f4 + q2 * f5
                   + (-_4q3) * f6;
    double Jt_f3 = -_2q1 * f1 + _2q4 * f2 - _4q3 * f3 + q4 * f4 + q1 * f5
                   + (-_4q2) * f6;
    double Jt_f4 = _2q2 * f1 + _2q3 * f2 + q1 * f4 + q2 * f5 + q3 * f6;

    /* Normalize gradient */
    double grad_norm = sqrt(Jt_f1 * Jt_f1 + Jt_f2 * Jt_f2 +
                            Jt_f3 * Jt_f3 + Jt_f4 * Jt_f4);
    if (grad_norm < 1e-10) grad_norm = 1e-10;
    double step = mf->beta / grad_norm;

    /* Gyro integration */
    double dqw = -0.5 * q2 * gyro.x - 0.5 * q3 * gyro.y - 0.5 * q4 * gyro.z;
    double dqx =  0.5 * q1 * gyro.x + 0.5 * q3 * gyro.z - 0.5 * q4 * gyro.y;
    double dqy =  0.5 * q1 * gyro.y - 0.5 * q2 * gyro.z + 0.5 * q4 * gyro.x;
    double dqz =  0.5 * q1 * gyro.z + 0.5 * q2 * gyro.y - 0.5 * q3 * gyro.x;

    mf->q_est.w += (dqw - step * Jt_f1) * mf->dt;
    mf->q_est.x += (dqx - step * Jt_f2) * mf->dt;
    mf->q_est.y += (dqy - step * Jt_f3) * mf->dt;
    mf->q_est.z += (dqz - step * Jt_f4) * mf->dt;
    mf->q_est = mems_quat_normalize(mf->q_est);

    return mf->q_est;
}

/* ─── L5: 1D Kalman Filter ─────────────────────────────────────────────── */

MemsKalman1D mems_kalman1d_init(double initial_state, double initial_cov,
                                 double process_noise, double meas_noise,
                                 double dt) {
    MemsKalman1D kf;
    kf.state = initial_state;
    kf.covariance = (initial_cov > 0.0) ? initial_cov : 1.0;
    kf.process_noise = (process_noise >= 0.0) ? process_noise : 0.001;
    kf.meas_noise = (meas_noise > 0.0) ? meas_noise : 0.01;
    kf.dt = (dt > 0.0) ? dt : 0.01;
    return kf;
}

void mems_kalman1d_predict(MemsKalman1D *kf, double rate) {
    if (!kf) return;
    /* State update: θ̂ = θ̂ + ω·dt */
    kf->state += rate * kf->dt;
    /* Covariance update: P = P + Q */
    kf->covariance += kf->process_noise;
}

void mems_kalman1d_update(MemsKalman1D *kf, double measurement) {
    if (!kf) return;
    /* Kalman gain */
    double K = kf->covariance / (kf->covariance + kf->meas_noise);
    /* State update */
    kf->state += K * (measurement - kf->state);
    /* Covariance update (Joseph form) */
    kf->covariance = (1.0 - K) * kf->covariance;
    if (kf->covariance < 1e-15) kf->covariance = 1e-15;
}

/* ─── L6: Roll/Pitch from Accelerometer ──────────────────────────────────── */

MemsEuler mems_accel_to_roll_pitch(MemsVector3 accel) {
    MemsEuler e;
    /* Roll: rotation about x-axis
     * roll = atan2(ay, az) */
    e.roll = atan2(accel.y, accel.z);

    /* Pitch: rotation about y-axis
     * pitch = atan2(-ax, sqrt(ay² + az²)) */
    double ayz = sqrt(accel.y * accel.y + accel.z * accel.z);
    e.pitch = atan2(-accel.x, ayz);

    /* Yaw cannot be determined from accelerometer alone */
    e.yaw = 0.0;

    return e;
}

/* ─── L6: Tilt-Compensated Yaw from Magnetometer ──────────────────────────
 *
 * To get true heading, the magnetometer must be tilt-compensated.
 * Rotate magnetic field to the horizontal plane using roll and pitch.
 */

double mems_mag_to_yaw(MemsVector3 mag, double roll, double pitch,
                        double declination) {
    /* De-rotate magnetometer reading using roll and pitch */
    double cr = cos(roll), sr = sin(roll);
    double cp = cos(pitch), sp = sin(pitch);

    /* Tilt-compensated x (north) and y (east) */
    double mx_h = mag.x * cp + mag.y * sr * sp + mag.z * cr * sp;
    double my_h = mag.y * cr - mag.z * sr;

    /* Yaw (heading = 0 is magnetic north) */
    double yaw = atan2(-my_h, mx_h);

    /* Apply magnetic declination */
    yaw += declination;

    /* Normalize to [-π, π] */
    while (yaw > M_PI) yaw -= 2.0 * M_PI;
    while (yaw < -M_PI) yaw += 2.0 * M_PI;

    return yaw;
}

/* ─── L6: Dead Reckoning ────────────────────────────────────────────────── */

void mems_dead_reckon(MemsVector3 accel_body, MemsQuaternion orientation,
                       double dt, MemsVector3 *position,
                       MemsVector3 *velocity, MemsVector3 gravity) {
    if (!position || !velocity || dt <= 0.0) return;

    /* Rotate acceleration to world frame */
    MemsVector3 accel_world = mems_quat_rotate(orientation, accel_body);

    /* Remove gravity */
    MemsVector3 accel_linear;
    accel_linear.x = accel_world.x - gravity.x;
    accel_linear.y = accel_world.y - gravity.y;
    accel_linear.z = accel_world.z - gravity.z;

    /* Integrate velocity */
    velocity->x += accel_linear.x * dt;
    velocity->y += accel_linear.y * dt;
    velocity->z += accel_linear.z * dt;

    /* Integrate position */
    position->x += velocity->x * dt;
    position->y += velocity->y * dt;
    position->z += velocity->z * dt;
}
