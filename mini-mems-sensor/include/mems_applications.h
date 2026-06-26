/**
 * @file mems_applications.h
 * @brief L7-L8 MEMS Sensor Applications — Smartphone IMU, drone stabilization,
 *        automotive safety, wearable fitness, industrial condition monitoring
 *
 * Covers:
 *   L7 - Smartphone IMU (iPhone-like orientation), drone AHRS (Quadrotor),
 *        automotive airbag crash detection, wearable step counter,
 *        industrial vibration monitoring, automotive ESC (Electronic Stability)
 *   L8 - GPS-denied navigation, MEMS-based magnetic anomaly detection,
 *        sensor array processing, structured-light MEMS mirror
 *
 * Reference: Niu et al., "Using Inertial Sensors in Smartphones" (2012)
 *            Beard & McLain, "Small Unmanned Aircraft: Theory and Practice" (2012)
 *            Fleming, "Overview of Automotive Sensors" (2001)
 *            Chopra, "MEMS for Automotive Applications" (2006)
 * Courses: MIT 16.485 (GNSS), Stanford EE267 (VR/AR), Michigan EECS 411,
 *          Tsinghua 汽车电子
 */

#ifndef MEMS_APPLICATIONS_H
#define MEMS_APPLICATIONS_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include "mems_basics.h"
#include "mems_fusion.h"

/* ─── L7: Smartphone IMU Application ───────────────────────────────────────── */

/**
 * @brief Smartphone orientation state (iPhone CoreMotion equivalent)
 *        Reference: Apple CoreMotion Framework, Android SensorManager
 *        Real-world keywords: iPhone, Samsung, Huawei
 */
typedef struct {
    MemsQuaternion device_attitude;  /* device orientation */
    MemsEuler      device_euler;     /* roll/pitch/yaw */
    MemsVector3    user_accel;       /* acceleration without gravity */
    MemsVector3    gravity;          /* gravity vector */
    MemsVector3    rotation_rate;    /* gyro rate */
    MemsVector3    magnetic_field;   /* calibrated mag field */
    double         heading;          /* magnetic heading [rad] */
    double         heading_accuracy;  /* [rad] */
    uint8_t        device_motion;    /* 0=static, 1=moving, 2=shaking */
} MemsPhoneOrientation;

/**
 * @brief Process raw IMU data into phone orientation
 *        Modeled on iPhone 15 Pro / Galaxy S25 sensor pipeline:
 *        - Low-pass accel for gravity direction
 *        - High-pass accel for user acceleration
 *        - Complementary filter for orientation
 *        - Tilt-compensated compass for heading
 *
 * @param accel  Raw accelerometer [m/s²]
 * @param gyro   Raw gyroscope [rad/s]
 * @param mag    Raw magnetometer [μT]
 * @param dt     Time step [s]
 * @param state  Previous state (updated in-place)
 * @param declination Local magnetic declination [rad]
 */
void mems_phone_orientation_update(const MemsVector3 *accel,
                                    const MemsVector3 *gyro,
                                    const MemsVector3 *mag,
                                    double dt,
                                    MemsPhoneOrientation *state,
                                    double declination);

/* ─── L7: Drone Stabilization AHRS ──────────────────────────────────────────── */

/**
 * @brief Quadrotor/drone attitude state
 *        Keywords: Quadrotor, DJI, ArduPilot, PX4, drone
 */
typedef struct {
    MemsQuaternion attitude;       /* current attitude */
    MemsEuler      euler;          /* roll/pitch/yaw */
    MemsVector3    angular_rate;   /* filtered gyro [rad/s] */
    MemsVector3    linear_accel;   /* body-frame acceleration */
    double         altitude;       /* baro-derived altitude [m] */
    double         vertical_speed; /* baro-derived rate [m/s] */
    double         motor_throttle; /* normalized [0,1] */
    uint8_t        armed;          /* flight status */
} MemsDroneAttitude;

/**
 * @brief Drone stabilization filter using IMU + optional barometer
 *        Implements Madgwick AHRS for attitude + altitude hold filter.
 *        Reference: Beard & McLain (2012) §8, ArduPilot EKF implementation
 */
void mems_drone_attitude_update(const MemsVector3 *accel,
                                 const MemsVector3 *gyro,
                                 const MemsVector3 *mag,
                                 double pressure_pa,
                                 double temperature_c,
                                 double dt,
                                 MemsDroneAttitude *drone);

/**
 * @brief Convert pressure to altitude using barometric formula
 *        h = (T₀/L) · (1 - (P/P₀)^(R·L/g))
 *        ISA standard: T₀ = 288.15 K, P₀ = 101325 Pa, L = -0.0065 K/m
 */
double mems_baro_to_altitude(double pressure_pa, double temperature_k,
                              double sea_level_pressure_pa);

/* ─── L7: Automotive Airbag Crash Detection ───────────────────────────────── */

/**
 * @brief Airbag crash sensor state
 *        Keywords: Airbag, Bosch, Continental, ESC, NCAP, crash
 *        Reference: Fleming (2001), ISO 12097-2 (Airbag testing)
 */
typedef struct {
    MemsVector3 accel;            /* current acceleration [m/s²] */
    MemsVector3 velocity;         /* integrated velocity [m/s] */
    double      jerk;             /* derivative of acceleration [m/s³] */
    double      abs_accel;        /* resultant acceleration magnitude */
    double      velocity_change;  /* Δv [m/s] */
    double      max_deceleration; /* peak deceleration [g] */
    double      crash_time_ms;    /* time since crash detection [ms] */
    uint8_t     fire_decision;    /* 0=normal, 1=pretensioner, 2=airbag */
    uint8_t     crash_severity;   /* 0-3 */
} MemsAirbagState;

/**
 * @brief Evaluate crash condition from accelerometer data
 *        Criteria: Δv threshold, deceleration magnitude, jerk profile.
 *        Based on Bosch/Continental airbag algorithms (Fleming 2001).
 */
void mems_airbag_crash_check(const MemsVector3 *accel, double dt,
                              MemsAirbagState *state);

/* ─── L7: Wearable Step Counter ────────────────────────────────────────────── */

/**
 * @brief Pedometer step counter state
 *        Keywords: Fitbit, Apple Watch, Garmin, Nike, step
 *        Reference: Zhao, "Pedestrian Dead Reckoning" (2016)
 */
typedef struct {
    uint32_t step_count;         /* total steps */
    double   step_length;        /* estimated stride length [m] */
    double   distance_walked;    /* total distance [m] */
    double   cadence;            /* steps per minute */
    double   last_step_time;     /* timestamp of last step [s] */
    double   accel_magnitude;    /* current accel mag */
    double   threshold_high;     /* upper threshold for step detection */
    double   threshold_low;      /* lower threshold (dynamic) */
    double   sample_rate;        /* [Hz] */
    uint8_t  walking;           /* currently walking flag */
} MemsPedometer;

/**
 * @brief Initialize pedometer with sample rate
 */
MemsPedometer mems_pedometer_init(double sample_rate);

/**
 * @brief Update pedometer with accelerometer magnitude
 *        Returns 1 if a step was detected this update
 */
int mems_pedometer_update(MemsPedometer *ped, double accel_mag, double timestamp);

/* ─── L7: Industrial Condition Monitoring ─────────────────────────────────── */

/**
 * @brief Vibration analysis state for rotating machinery
 *        Keywords: ISO 10816, Bently Nevada, SKF, bearing, pump
 *        Reference: ISO 10816-1:1995 (Vibration Severity), Mobley (2002)
 */
typedef struct {
    double vibration_rms;        /* RMS vibration [mm/s rms] */
    double vibration_peak;       /* peak vibration [mm/s] */
    double crest_factor;         /* peak/RMS ratio */
    double iso_10816_zone;       /* A=good, B=acceptable, C=alert, D=danger */
    double dominant_freq;        /* dominant frequency [Hz] */
    double running_speed;        /* RPM */
    double bearing_freq[4];      /* BPFO, BPFI, BSF, FTF [Hz] */
    uint8_t severity_level;      /* 0=normal, 1=advisory, 2=warning, 3=alarm */
} MemsVibrationMonitor;

/**
 * @brief Evaluate vibration severity per ISO 10816-1
 *        Classifies into zones A-D based on machine type and power rating.
 *
 * @param vrms       RMS vibration velocity [mm/s]
 * @param machine_type 1=small, 2=medium, 3=large rigid, 4=large flexible
 * @return Zone: 0=A, 1=B, 2=C, 3=D
 */
int mems_iso_10816_zone(double vrms, int machine_type);

/**
 * @brief Compute bearing fault frequencies
 *        BPFO = (N/2)·f_r·(1 - (d/D)·cos(φ))
 *        BPFI = (N/2)·f_r·(1 + (d/D)·cos(φ))
 *        BSF = (D/(2d))·f_r·(1 - (d/D)²·cos²(φ))
 *        FTF = (1/2)·f_r·(1 - (d/D)·cos(φ))
 *        where N=balls, f_r=shaft freq, d=ball dia, D=pitch dia, φ=contact angle
 *        Reference: Randall & Antoni, MSSP (2011)
 *
 * @param shaft_hz    Shaft rotation frequency [Hz] (RPM/60)
 * @param n_balls     Number of rolling elements
 * @param ball_dia    Ball diameter [mm]
 * @param pitch_dia   Pitch diameter [mm]
 * @param contact_deg Contact angle [degrees]
 * @param bpfo        Output: Ball Pass Frequency Outer Race [Hz]
 * @param bpfi        Output: Ball Pass Frequency Inner Race [Hz]
 * @param bsf         Output: Ball Spin Frequency [Hz]
 * @param ftf         Output: Fundamental Train Frequency [Hz]
 */
void mems_bearing_fault_freqs(double shaft_hz, int n_balls,
                               double ball_dia, double pitch_dia,
                               double contact_deg,
                               double *bpfo, double *bpfi,
                               double *bsf, double *ftf);

/* ─── L7: Automotive ESC (Electronic Stability Control) ──────────────────── */

/**
 * @brief ESC vehicle dynamics state
 *        Keywords: ESC, Bosch ESP, Continental, Toyota, Tesla
 */
typedef struct {
    double yaw_rate;             /* measured yaw rate [rad/s] */
    double yaw_rate_ref;         /* reference yaw rate from steering [rad/s] */
    double sideslip_angle;       /* vehicle sideslip angle [rad] */
    double steering_angle;       /* from steering wheel [rad] */
    double vehicle_speed;        /* [m/s] */
    double wheel_speed_fl;       /* [m/s] */
    double wheel_speed_fr;
    double wheel_speed_rl;
    double wheel_speed_rr;
    double lateral_accel;        /* lateral acceleration [m/s²] */
    uint8_t esc_intervention;    /* 0=off, 1=brake, 2=brake+throttle */
    uint8_t stability_flag;      /* 0=stable, 1=understeer, 2=oversteer */
} MemsEscState;

/**
 * @brief Evaluate ESC intervention need
 *        Compares measured yaw rate to reference yaw rate (driver intent).
 *        Reference: van Zanten, Bosch ESP (2000), Rajamani (2012)
 *
 * @param yaw_rate         Measured yaw rate from gyro [rad/s]
 * @param steering_angle   Steering wheel angle [rad]
 * @param vehicle_speed    Vehicle speed [m/s]
 * @param wheelbase        Vehicle wheelbase [m]
 * @param understeer_grad  Understeer gradient [rad·s²/m]
 * @param state            ESC state (updated in-place)
 */
void mems_esc_check(double yaw_rate, double steering_angle,
                     double vehicle_speed, double wheelbase,
                     double understeer_grad, MemsEscState *state);

/* ─── L8: GPS-Denied Navigation ──────────────────────────────────────────────
 * Indoor/underground navigation using IMU dead reckoning with zero-velocity
 * updates (ZUPT). Reference: Foxlin (2005), "Pedestrian Tracking with Shoe-Mounted
 * Inertial Sensors". Keywords: GPS-denied, indoor, tunnel, subway.
 */

typedef struct {
    MemsVector3 position;
    MemsVector3 velocity;
    MemsQuaternion attitude;
    double      covariance[3][3];
    double      accel_bias[3];
    double      gyro_bias[3];
    uint8_t     stationary;
    double      stationary_thresh;
    uint32_t    zupt_count;
} MemsInsState;

/**
 * @brief Initialize INS state for GPS-denied navigation
 */
MemsInsState mems_ins_init(void);

/**
 * @brief ZUPT-aided INS update
 *        If stationary detected, apply zero-velocity update to bound drift.
 *        Reference: Foxlin (2005), Jiménez et al. (2010)
 *
 * @param accel  Accelerometer [m/s²] body frame
 * @param gyro   Gyroscope [rad/s] body frame
 * @param dt     Time step [s]
 * @param ins    INS state (updated in-place)
 */
void mems_ins_zupt_update(const MemsVector3 *accel, const MemsVector3 *gyro,
                           double dt, MemsInsState *ins);

/**
 * @brief Detect stationary condition from accelerometer + gyro
 *        Uses acceleration magnitude variance and gyro magnitude threshold.
 *
 * @param accel_buf  Ring buffer of recent accel measurements
 * @param gyro_buf   Ring buffer of recent gyro measurements  
 * @param buf_len    Number of samples in buffers
 * @param accel_thresh Acceleration magnitude variance threshold [m²/s⁴]
 * @param gyro_thresh  Gyro magnitude threshold [rad/s]
 * @return 1 if stationary, 0 otherwise
 */
int mems_detect_stationary(const MemsVector3 *accel_buf,
                            const MemsVector3 *gyro_buf,
                            size_t buf_len,
                            double accel_thresh,
                            double gyro_thresh);

/* ─── L8: MEMS Mirror for Structured Light ────────────────────────────────── */

/**
 * @brief 1D MEMS micromirror scan state
 *        Keywords: structured light, Texas Instruments DLP, STMicro LBS
 *        Reference: Yalcinkaya et al., "MEMS Scanner" (2006)
 */
typedef struct {
    double scan_angle;           /* current mechanical angle [rad] */
    double scan_amplitude;       /* peak-to-peak amplitude [rad] */
    double scan_frequency;       /* resonant scan frequency [Hz] */
    double drive_voltage;        /* drive signal voltage [V] */
    double mirror_diameter;      /* [mm] */
    double optical_angle;        /* optical angle = 2× mechanical [rad] */
    double q_factor;             /* mirror Q */
    uint8_t scan_axis;           /* 0=horizontal, 1=vertical */
} MemsMirrorState;

/**
 * @brief Compute MEMS mirror optical scan angle from mechanical
 *        θ_optical = 2 · θ_mechanical
 */
double mems_mirror_optical_angle(double mechanical_angle_rad);

/**
 * @brief Compute mirror resonant frequency from dimensions
 *        f_n = (1/2π) · sqrt(k_torsion / I_mirror)
 *        where k_torsion depends on torsion beam geometry
 *        and I_mirror is the mirror moment of inertia.
 *
 * @param beam_length   Torsion beam length [m]
 * @param beam_width    Torsion beam width [m]
 * @param beam_thickness Torsion beam thickness [m]
 * @param shear_modulus G [Pa], Si ≈ 79.9 GPa
 * @param mirror_mass   Mirror mass [kg]
 * @param mirror_radius Mirror radius [m] (assumes disk shape)
 */
double mems_mirror_resonant_freq(double beam_length, double beam_width,
                                  double beam_thickness, double shear_modulus,
                                  double mirror_mass, double mirror_radius);

#endif /* MEMS_APPLICATIONS_H */
