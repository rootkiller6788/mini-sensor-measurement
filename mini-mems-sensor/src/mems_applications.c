/**
 * @file mems_applications.c
 * @brief L7-L8 MEMS Applications — Smartphone IMU, drone AHRS, automotive
 *        airbag/ESC, pedometer, vibration monitoring, GPS-denied INS
 *
 * Reference: Niu et al. (2012) "Using Inertial Sensors in Smartphones"
 *            Beard & McLain (2012) "Small Unmanned Aircraft"
 *            Fleming (2001) "Overview of Automotive Sensors"
 *            Foxlin (2005) "Pedestrian Tracking with Shoe-Mounted IMU"
 *            ISO 10816-1:1995 (Vibration Severity)
 *            Randall & Antoni (2011) "Bearing Fault Detection"
 */

#include "mems_applications.h"
#include <math.h>
#include <string.h>
#include <float.h>

/* Standard gravity [m/s²] */
#define G_STANDARD 9.80665
/* Standard sea-level pressure [Pa] */
#define P0_SEA_LEVEL 101325.0
/* Standard temperature at sea level [K] */
#define T0_SEA_LEVEL 288.15
/* Temperature lapse rate [K/m] */
#define LAPSE_RATE -0.0065
/* Gas constant for dry air [J/(kg·K)] */
#define R_AIR 287.058
/* Ratio R·L/g */
#define BARO_EXPONENT (R_AIR * (-LAPSE_RATE) / G_STANDARD)

/* ─── L7: Smartphone Orientation ─────────────────────────────────────────── */

void mems_phone_orientation_update(const MemsVector3 *accel,
                                    const MemsVector3 *gyro,
                                    const MemsVector3 *mag,
                                    double dt,
                                    MemsPhoneOrientation *state,
                                    double declination) {
    if (!accel || !gyro || !mag || !state || dt <= 0.0) return;

    /* Low-pass filter accelerometer to separate gravity */
    static const double alpha_lp = 0.05;  /* low-pass coefficient */
    MemsVector3 gravity_lp;
    gravity_lp.x = (1.0 - alpha_lp) * state->gravity.x + alpha_lp * accel->x;
    gravity_lp.y = (1.0 - alpha_lp) * state->gravity.y + alpha_lp * accel->y;
    gravity_lp.z = (1.0 - alpha_lp) * state->gravity.z + alpha_lp * accel->z;

    /* Normalize gravity estimate */
    double g_norm = sqrt(gravity_lp.x * gravity_lp.x +
                          gravity_lp.y * gravity_lp.y +
                          gravity_lp.z * gravity_lp.z);
    if (g_norm > 0.1) {
        state->gravity.x = gravity_lp.x / g_norm;
        state->gravity.y = gravity_lp.y / g_norm;
        state->gravity.z = gravity_lp.z / g_norm;
    }

    /* User acceleration = raw - low-pass */
    state->user_accel.x = accel->x - gravity_lp.x;
    state->user_accel.y = accel->y - gravity_lp.y;
    state->user_accel.z = accel->z - gravity_lp.z;

    /* Store raw rotation rate */
    state->rotation_rate = *gyro;

    /* Magnetic field */
    state->magnetic_field = *mag;

    /* Compute orientation using simple complementary filter */
    MemsCompFilter cf;
    cf.alpha = 0.98;
    cf.dt = dt;
    cf.q_est = state->device_attitude;
    MemsVector3 unit_gravity = state->gravity;
    MemsVector3 unit_grav_neg = {-unit_gravity.x, -unit_gravity.y,
                                   -unit_gravity.z};
    state->device_attitude = mems_comp_filter_update(&cf, *gyro,
                                                      unit_grav_neg);
    state->device_euler = mems_quat_to_euler(state->device_attitude);

    /* Tilt-compensated heading */
    state->heading = mems_mag_to_yaw(*mag, state->device_euler.roll,
                                      state->device_euler.pitch,
                                      declination);

    /* Motion detection: accelerometer magnitude variance */
    double accel_mag = sqrt(accel->x * accel->x +
                             accel->y * accel->y +
                             accel->z * accel->z);
    double accel_diff = fabs(accel_mag - G_STANDARD);
    if (accel_diff < 0.3) state->device_motion = 0;           /* static */
    else if (accel_diff < 2.0) state->device_motion = 1;      /* moving */
    else state->device_motion = 2;                            /* shaking */
}

/* ─── L7: Barometric Altitude ──────────────────────────────────────────────
 *
 * The barometric formula relates pressure to altitude in the ISA
 * (International Standard Atmosphere):
 *
 *   h = (T₀/L) · (1 - (P/P₀)^(R·L/g))
 *
 * For the troposphere (h < 11000 m):
 *   T₀ = 288.15 K, P₀ = 101325 Pa, L = -0.0065 K/m, R = 287.058 J/(kg·K)
 */

double mems_baro_to_altitude(double pressure_pa, double temperature_k,
                              double sea_level_pressure_pa) {
    (void)temperature_k;  /* used in non-ISA implementations */
    if (pressure_pa <= 0.0 || sea_level_pressure_pa <= 0.0) return 0.0;

    double ratio = pressure_pa / sea_level_pressure_pa;
    if (ratio <= 0.0) return INFINITY;

    /* ISA barometric formula */
    double exponent = R_AIR * (-LAPSE_RATE) / G_STANDARD;
    double altitude = (T0_SEA_LEVEL / (-LAPSE_RATE)) *
                       (1.0 - pow(ratio, exponent));

    return altitude;
}

/* ─── L7: Drone Attitude ─────────────────────────────────────────────────── */

void mems_drone_attitude_update(const MemsVector3 *accel,
                                 const MemsVector3 *gyro,
                                 const MemsVector3 *mag,
                                 double pressure_pa,
                                 double temperature_c,
                                 double dt,
                                 MemsDroneAttitude *drone) {
    if (!accel || !gyro || !mag || !drone) return;

    /* Use Mahony filter for robust attitude (better in dynamic conditions) */
    MemsMahonyFilter mf;
    mf.kp = 1.0;
    mf.ki = 0.05;
    mf.dt = dt;
    mf.q_est = drone->attitude;
    mf.bias_int = drone->angular_rate;  /* reuse angular_rate as bias integral */

    drone->attitude = mems_mahony_update(&mf, *gyro, *accel, *mag);
    drone->euler = mems_quat_to_euler(drone->attitude);

    /* Filtered angular rate */
    static const double beta_gyro = 0.15;
    drone->angular_rate.x = (1.0 - beta_gyro) * drone->angular_rate.x +
                             beta_gyro * gyro->x;
    drone->angular_rate.y = (1.0 - beta_gyro) * drone->angular_rate.y +
                             beta_gyro * gyro->y;
    drone->angular_rate.z = (1.0 - beta_gyro) * drone->angular_rate.z +
                             beta_gyro * gyro->z;

    /* Body-frame linear acceleration */
    drone->linear_accel = *accel;

    /* Barometric altitude */
    double temp_k = temperature_c + 273.15;
    drone->altitude = mems_baro_to_altitude(pressure_pa, temp_k,
                                             P0_SEA_LEVEL);

    /* Vertical speed from baro delta (simple differentiation with low-pass) */
    static double prev_altitude = 0.0;
    if (prev_altitude > 0.0) {
        double raw_vs = (drone->altitude - prev_altitude) / dt;
        drone->vertical_speed = 0.7 * drone->vertical_speed + 0.3 * raw_vs;
    }
    prev_altitude = drone->altitude;
}

/* ─── L7: Automotive Airbag Crash Detection ─────────────────────────────────
 *
 * Airbag deployment algorithms must discriminate between real crashes
 * and non-crash events (potholes, door slams, etc.).
 *
 * Key criteria (per Bosch/Continental algorithms):
 *   - Δv (velocity change) exceeds threshold
 *   - Deceleration magnitude exceeds threshold for minimum time
 *   - Jerk profile matches crash pulse
 *
 * The algorithm monitors a sliding window of acceleration samples.
 * Reference: Fleming (2001), ISO 12097-2
 */

void mems_airbag_crash_check(const MemsVector3 *accel, double dt,
                              MemsAirbagState *state) {
    if (!accel || !state || dt <= 0.0) return;

    state->accel = *accel;

    /* Resultant acceleration magnitude */
    state->abs_accel = sqrt(accel->x * accel->x +
                             accel->y * accel->y +
                             accel->z * accel->z);

    /* Integrate acceleration to get velocity change */
    MemsVector3 accel_increment = {accel->x * dt, accel->y * dt,
                                    accel->z * dt};
    state->velocity.x += accel_increment.x;
    state->velocity.y += accel_increment.y;
    state->velocity.z += accel_increment.z;

    state->velocity_change = sqrt(state->velocity.x * state->velocity.x +
                                   state->velocity.y * state->velocity.y +
                                   state->velocity.z * state->velocity.z);

    /* Jerk: derivative of acceleration (sliding window difference)
     * Simple single-sample jerk: jerk = (a_k - a_{k-1})/dt */
    static double prev_accel_mag = 0.0;
    state->jerk = (state->abs_accel - prev_accel_mag) / dt;
    prev_accel_mag = state->abs_accel;

    /* Track maximum deceleration in g's */
    double decel_g = state->abs_accel / G_STANDARD;
    if (decel_g > state->max_deceleration)
        state->max_deceleration = decel_g;

    /* Crash detection logic:
     * - Fire pretensioner: Δv > 5 m/s OR decel > 2g sustained
     * - Fire airbag: Δv > 12 m/s AND decel > 5g for > 10ms
     *
     * Simplified for demonstration. Real algorithms use much more sophisticated
     * pattern recognition (wavelet, neural network, etc.).
     */
    if (state->velocity_change > 5.0 || decel_g > 2.0) {
        state->fire_decision = 1;  /* pretensioner */
        if (state->crash_time_ms < 1.0)
            state->crash_time_ms = 0.0;
    }

    if (state->velocity_change > 12.0 && decel_g > 5.0) {
        state->fire_decision = 2;  /* airbag */
    }

    /* Crash severity */
    if (decel_g < 3.0) state->crash_severity = 1;
    else if (decel_g < 8.0) state->crash_severity = 2;
    else state->crash_severity = 3;

    state->crash_time_ms += dt * 1000.0;
}

/* ─── L7: Wearable Step Counter ───────────────────────────────────────────── */

MemsPedometer mems_pedometer_init(double sample_rate) {
    MemsPedometer ped;
    memset(&ped, 0, sizeof(ped));
    ped.sample_rate = (sample_rate > 0.0) ? sample_rate : 50.0;
    ped.threshold_high = 10.5;   /* m/s², above gravity */
    ped.threshold_low  = 9.5;    /* m/s², below gravity */
    ped.step_length = 0.75;      /* m, typical adult */
    return ped;
}

int mems_pedometer_update(MemsPedometer *ped, double accel_mag,
                           double timestamp) {
    if (!ped) return 0;

    int step_detected = 0;

    /* Simple peak detection with hysteresis */
    if (!ped->walking && accel_mag > ped->threshold_high) {
        /* Rising edge: potential step start */
        ped->walking = 1;
    } else if (ped->walking && accel_mag < ped->threshold_low) {
        /* Falling edge: step completed */
        double time_since_last = timestamp - ped->last_step_time;
        /* Enforce minimum step time (max cadence ~ 200 steps/min) */
        if (time_since_last > 0.3) {
            ped->step_count++;
            ped->distance_walked += ped->step_length;
            ped->cadence = 60.0 / time_since_last;
            ped->last_step_time = timestamp;
            step_detected = 1;
        }
        ped->walking = 0;
    }

    /* Dynamic threshold adaptation */
    double step_avg = 0.5 * (ped->threshold_high + ped->threshold_low);
    if (step_avg > 0.0) {
        ped->threshold_high = step_avg + 0.3;
        ped->threshold_low  = step_avg - 0.3;
    }

    return step_detected;
}

/* ─── L7: Industrial Vibration Monitoring (ISO 10816-1) ──────────────────── */

int mems_iso_10816_zone(double vrms, int machine_type) {
    /* ISO 10816-1:1995 vibration severity zones
     * Class I:  Small machines (< 15 kW)
     * Class II: Medium machines (15-75 kW)
     * Class III: Large rigid foundation (> 75 kW)
     * Class IV: Large flexible foundation
     *
     * Returns zone: 0=A (good), 1=B (acceptable), 2=C (alert), 3=D (danger)
     *
     * Zone boundaries [mm/s rms]:
     *         | A/B  | B/C  | C/D
     * Class I | 1.4  | 2.8  | 4.5
     * Class II| 2.8  | 4.5  | 7.1
     * ClassIII| 4.5  | 7.1  | 11.0
     * Class IV| 7.1  | 11.0 | 18.0
     */
    static const double boundaries[4][3] = {
        {1.4,  2.8,  4.5 },   /* Class I */
        {2.8,  4.5,  7.1 },   /* Class II */
        {4.5,  7.1,  11.0},   /* Class III */
        {7.1,  11.0, 18.0}    /* Class IV */
    };

    if (machine_type < 1 || machine_type > 4) machine_type = 2;
    int idx = machine_type - 1;

    if (vrms <= boundaries[idx][0]) return 0;  /* A: Good */
    if (vrms <= boundaries[idx][1]) return 1;  /* B: Acceptable */
    if (vrms <= boundaries[idx][2]) return 2;  /* C: Alert */
    return 3;                                   /* D: Danger */
}

void mems_bearing_fault_freqs(double shaft_hz, int n_balls,
                               double ball_dia, double pitch_dia,
                               double contact_deg,
                               double *bpfo, double *bpfi,
                               double *bsf, double *ftf) {
    if (!bpfo || !bpfi || !bsf || !ftf) return;
    if (shaft_hz <= 0.0 || n_balls < 1 ||
        ball_dia <= 0.0 || pitch_dia <= 0.0) {
        *bpfo = *bpfi = *bsf = *ftf = 0.0;
        return;
    }

    double contact_rad = contact_deg * M_PI / 180.0;
    double cos_phi = cos(contact_rad);
    double d_div_D = ball_dia / pitch_dia;

    /* Ball Pass Frequency Outer Race (BPFO) */
    *bpfo = (double)n_balls / 2.0 * shaft_hz * (1.0 - d_div_D * cos_phi);

    /* Ball Pass Frequency Inner Race (BPFI) */
    *bpfi = (double)n_balls / 2.0 * shaft_hz * (1.0 + d_div_D * cos_phi);

    /* Ball Spin Frequency (BSF) */
    *bsf = pitch_dia / (2.0 * ball_dia) * shaft_hz *
           (1.0 - d_div_D * d_div_D * cos_phi * cos_phi);

    /* Fundamental Train Frequency (FTF) = Cage Frequency */
    *ftf = 0.5 * shaft_hz * (1.0 - d_div_D * cos_phi);
}

/* ─── L7: Automotive ESC (Electronic Stability Control) ────────────────────
 *
 * ESC compares driver-intended yaw rate (from steering wheel angle) with
 * measured yaw rate (from gyroscope). Deviation indicates understeer or
 * oversteer.
 *
 * Reference yaw rate from bicycle model:
 *   ψ̇_ref = (v / L) · δ / (1 + K_US · v²)
 *
 * where L = wheelbase, δ = steer angle, K_US = understeer gradient,
 * v = vehicle speed.
 *
 * Reference: van Zanten (2000), "Bosch ESP Systems"
 *            Rajamani (2012), "Vehicle Dynamics and Control"
 */

void mems_esc_check(double yaw_rate, double steering_angle,
                     double vehicle_speed, double wheelbase,
                     double understeer_grad, MemsEscState *state) {
    if (!state) return;

    state->yaw_rate = yaw_rate;
    state->steering_angle = steering_angle;
    state->vehicle_speed = vehicle_speed;

    /* Reference yaw rate from bicycle model */
    if (wheelbase > 0.0 && fabs(vehicle_speed) > 0.1) {
        state->yaw_rate_ref = (vehicle_speed / wheelbase) *
                               steering_angle /
                               (1.0 + understeer_grad *
                                vehicle_speed * vehicle_speed);
    } else {
        state->yaw_rate_ref = 0.0;
    }

    /* Sideslip angle estimate from lateral acceleration + yaw rate
     * β = ∫(a_y/v - ψ̇) dt  (simplified kinematic)
     * Use static approximation: β ≈ (a_y/v - ψ̇) · τ */
    double lat_accel = state->lateral_accel;
    if (fabs(vehicle_speed) > 0.1) {
        double tau = 0.3;  /* sideslip time constant */
        state->sideslip_angle = (lat_accel / vehicle_speed - yaw_rate) * tau;
    }

    /* Stability assessment */
    double yaw_error = yaw_rate - state->yaw_rate_ref;
    double threshold = 0.05;  /* rad/s */

    if (fabs(yaw_error) < threshold) {
        state->stability_flag = 0;  /* stable */
        state->esc_intervention = 0;
    } else if (yaw_error < -threshold) {
        /* Understeer: vehicle turns less than intended */
        state->stability_flag = 1;
        state->esc_intervention = 1;  /* brake inner rear wheel */
    } else {
        /* Oversteer: vehicle turns more than intended */
        state->stability_flag = 2;
        state->esc_intervention = 2;  /* brake outer front wheel */
    }
}

/* ─── L8: GPS-Denied Navigation ──────────────────────────────────────────── */

MemsInsState mems_ins_init(void) {
    MemsInsState ins;
    memset(&ins, 0, sizeof(ins));
    ins.attitude = mems_quat_identity();
    ins.stationary_thresh = 0.15;  /* m²/s⁴ acceleration variance */
    ins.stationary = 1;
    ins.zupt_count = 0;
    return ins;
}

void mems_ins_zupt_update(const MemsVector3 *accel, const MemsVector3 *gyro,
                           double dt, MemsInsState *ins) {
    if (!accel || !gyro || !ins || dt <= 0.0) return;

    /* Remove estimated bias */
    MemsVector3 accel_corrected, gyro_corrected;
    accel_corrected.x = accel->x - ins->accel_bias[0];
    accel_corrected.y = accel->y - ins->accel_bias[1];
    accel_corrected.z = accel->z - ins->accel_bias[2];
    gyro_corrected.x = gyro->x - ins->gyro_bias[0];
    gyro_corrected.y = gyro->y - ins->gyro_bias[1];
    gyro_corrected.z = gyro->z - ins->gyro_bias[2];

    /* Gyro integration: update quaternion */
    MemsQuaternion omega_q = {0.0, gyro_corrected.x, gyro_corrected.y,
                               gyro_corrected.z};
    MemsQuaternion dq_temp = mems_quat_multiply(ins->attitude, omega_q);
    ins->attitude.w += dq_temp.w * 0.5 * dt;
    ins->attitude.x += dq_temp.x * 0.5 * dt;
    ins->attitude.y += dq_temp.y * 0.5 * dt;
    ins->attitude.z += dq_temp.z * 0.5 * dt;
    ins->attitude = mems_quat_normalize(ins->attitude);

    /* Rotate acceleration to world frame and remove gravity */
    MemsVector3 accel_world = mems_quat_rotate(ins->attitude,
                                                accel_corrected);
    MemsVector3 gravity = {0.0, 0.0, G_STANDARD};

    /* Dead reckoning */
    MemsVector3 accel_linear;
    accel_linear.x = accel_world.x - gravity.x;
    accel_linear.y = accel_world.y - gravity.y;
    accel_linear.z = accel_world.z - gravity.z;

    ins->velocity.x += accel_linear.x * dt;
    ins->velocity.y += accel_linear.y * dt;
    ins->velocity.z += accel_linear.z * dt;

    ins->position.x += ins->velocity.x * dt;
    ins->position.y += ins->velocity.y * dt;
    ins->position.z += ins->velocity.z * dt;

    /* Zero-velocity update if stationary */
    if (ins->stationary) {
        /* ZUPT: zero out velocity and bound position drift */
        ins->velocity.x *= 0.1;
        ins->velocity.y *= 0.1;
        ins->velocity.z *= 0.1;
        /* Estimate bias from gyro when stationary */
        for (int i = 0; i < 3; i++) {
            double alpha_bias = 0.001;
            if (i == 0) {
                ins->gyro_bias[0] += alpha_bias * gyro->x;
                ins->accel_bias[0] += alpha_bias *
                    (accel->x - ins->attitude.x * G_STANDARD);
            } else if (i == 1) {
                ins->gyro_bias[1] += alpha_bias * gyro->y;
                ins->accel_bias[1] += alpha_bias *
                    (accel->y - ins->attitude.y * G_STANDARD);
            } else {
                ins->gyro_bias[2] += alpha_bias * gyro->z;
                ins->accel_bias[2] += alpha_bias *
                    (accel->z - ins->attitude.z * G_STANDARD);
            }
        }
        ins->zupt_count++;
    }
}

int mems_detect_stationary(const MemsVector3 *accel_buf,
                            const MemsVector3 *gyro_buf,
                            size_t buf_len,
                            double accel_thresh,
                            double gyro_thresh) {
    if (!accel_buf || !gyro_buf || buf_len < 3) return 1;

    /* Compute acceleration magnitude variance */
    double mean_a = 0.0, mean_a_sq = 0.0;
    for (size_t i = 0; i < buf_len; i++) {
        double a_mag = sqrt(accel_buf[i].x * accel_buf[i].x +
                             accel_buf[i].y * accel_buf[i].y +
                             accel_buf[i].z * accel_buf[i].z);
        mean_a += a_mag;
        mean_a_sq += a_mag * a_mag;
    }
    mean_a /= (double)buf_len;
    mean_a_sq /= (double)buf_len;
    double var_a = mean_a_sq - mean_a * mean_a;
    if (var_a < 0.0) var_a = 0.0;

    /* Compute gyro magnitude mean */
    double mean_g = 0.0;
    for (size_t i = 0; i < buf_len; i++) {
        double g_mag = sqrt(gyro_buf[i].x * gyro_buf[i].x +
                             gyro_buf[i].y * gyro_buf[i].y +
                             gyro_buf[i].z * gyro_buf[i].z);
        mean_g += g_mag;
    }
    mean_g /= (double)buf_len;

    /* Stationary if both variance and angular rate are below threshold */
    return (var_a < accel_thresh && mean_g < gyro_thresh) ? 1 : 0;
}

/* ─── L8: MEMS Micromirror for Structured Light ──────────────────────────── */

double mems_mirror_optical_angle(double mechanical_angle_rad) {
    /* Optical scan angle = 2 × mechanical angle (law of reflection) */
    return 2.0 * mechanical_angle_rad;
}

double mems_mirror_resonant_freq(double beam_length, double beam_width,
                                  double beam_thickness, double shear_modulus,
                                  double mirror_mass, double mirror_radius) {
    /* MEMS torsional mirror resonance:
     *
     * Torsional spring constant for a rectangular beam:
     *   k_torsion = (G · β · w · t³) / L
     *   where β ≈ 1/3 for w >> t (thin beam approximation)
     *         β ≈ 1/3 - 0.21·(t/w)·(1 - t⁴/(12·w⁴)) for all aspect ratios
     *
     * Mirror moment of inertia (disk):
     *   I_mirror = 0.5 · m · r²
     *
     * Resonant frequency:
     *   f_n = (1/2π) · sqrt(k_torsion / I_mirror)
     */
    if (beam_length <= 0.0 || beam_width <= 0.0 || beam_thickness <= 0.0 ||
        mirror_mass <= 0.0 || mirror_radius <= 0.0) return 0.0;

    /* Torsion constant β */
    double t_over_w = beam_thickness / beam_width;
    double beta;
    if (beam_width > beam_thickness * 10.0) {
        beta = 1.0 / 3.0;  /* thin beam approximation */
    } else {
        beta = 1.0 / 3.0 - 0.21 * t_over_w *
               (1.0 - pow(t_over_w, 4) / 12.0);
    }
    if (beta < 0.05) beta = 0.05;

    /* Torsional spring constant */
    double k_torsion = shear_modulus * beta * beam_width *
                        beam_thickness * beam_thickness * beam_thickness /
                        beam_length;

    /* Mirror moment of inertia (disk about diameter) */
    double I_mirror = 0.25 * mirror_mass * mirror_radius * mirror_radius;

    if (I_mirror <= 0.0) return 0.0;
    return sqrt(k_torsion / I_mirror) / (2.0 * M_PI);
}
