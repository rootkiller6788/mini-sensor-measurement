/**
 * @file example_sensor_fusion.c
 * @brief L6 Example: MEMS IMU sensor fusion for orientation estimation
 *
 * Demonstrates complementary filter, Mahony filter, and Madgwick AHRS
 * for fusing gyroscope + accelerometer + magnetometer data to estimate
 * 3D orientation (attitude and heading reference system).
 *
 * Real-world systems:
 *   - DJI Phantom/Matrice drone flight controller (Mahony/Madgwick)
 *   - Apple iPhone ARKit / CoreMotion (proprietary filter)
 *   - ArduPilot / PX4 autopilot EKF
 *   - Boston Dynamics Spot robot IMU
 *   - GoPro HyperSmooth stabilization
 *   - Tesla vehicle dynamics IMU
 */

#include <stdio.h>
#include <math.h>
#include "mems_basics.h"
#include "mems_fusion.h"
#include "mems_calibration.h"

#define G 9.80665

int main(void) {
    printf("=== MEMS Sensor Fusion for Orientation Estimation ===\n\n");

    /* ─── Initialize filters ─────────────────────────────────────── */
    printf("1. Sensor Fusion Filter Initialization\n\n");

    double dt = 0.01;  /* 100 Hz IMU sample rate */
    double beta = 0.04; /* Madgwick filter gain */

    MemsCompFilter cf = mems_comp_filter_init(0.98, dt);
    MemsMahonyFilter mf = mems_mahony_init(1.0, 0.05, dt);
    MemsMadgwickFilter mwf = mems_madgwick_init(beta, dt);

    printf("   Sample rate: %.0f Hz (dt = %.3f s)\n", 1.0 / dt, dt);
    printf("   Comp filter alpha: 0.98 (gyro trust)\n");
    printf("   Mahony Kp: 1.0, Ki: 0.05\n");
    printf("   Madgwick beta: 0.04 (gyro noise tuned)\n\n");

    /* ─── Simulated IMU data ─────────────────────────────────────── */
    printf("2. Processing Simulated IMU Data\n\n");

    /* Simulate a 90° rotation about the Y axis (pitch up from
     * flat to upright) over 1 second, then hold.

     * At t=0: accel = (0,0,g), gyro = (0,0,0), mag = (Bx,0,Bz)
     * At t=1: accel = (-g,0,0), gyro = (0,0,0)
     * During rotation: gyro_y ≈ π/2 rad/s */

    double rotation_rate = M_PI / 2.0;  /* 90 deg/s for 1 second */
    MemsQuaternion q_comp = mems_quat_identity();
    MemsQuaternion q_mahony = mems_quat_identity();
    MemsQuaternion q_madgwick = mems_quat_identity();

    printf("   Simulating 90° pitch rotation over 1 second...\n\n");
    printf("   Time(ms) | Comp Filter(R/P)  | Mahony(R/P)       | Madgwick(R/P)\n");
    printf("   ---------|--------------------|--------------------|-----------------\n");

    for (int step = 0; step <= 100; step++) {
        double t = step * dt;

        /* Generate synthetic sensor data */
        MemsVector3 gyro = {0.0, 0.0, 0.0};
        MemsVector3 accel = {0.0, 0.0, G};
        MemsVector3 mag = {30.0, 0.0, 40.0};  /* 50μT total, ~53° inclination */

        if (t < 1.0) {
            /* During rotation: constant angular velocity about Y */
            gyro.y = rotation_rate;

            /* Gravity vector rotates in XZ plane */
            double angle = rotation_rate * t;
            accel.x = -G * sin(angle);
            accel.z = G * cos(angle);

            /* Mag field rotates similarly */
            double bx = 30.0;
            double bz = 40.0;
            mag.x = bx * cos(angle) + bz * sin(angle);
            mag.z = -bx * sin(angle) + bz * cos(angle);
        } else {
            /* End state: total 90° rotation */
            accel.x = -G;
            accel.z = 0.0;
            mag.x = 40.0;
            mag.z = -30.0;
        }

        /* Update filters */
        cf.q_est = q_comp;
        q_comp = mems_comp_filter_update(&cf, gyro, accel);

        mf.q_est = q_mahony;
        q_mahony = mems_mahony_update(&mf, gyro, accel, mag);

        mwf.q_est = q_madgwick;
        q_madgwick = mems_madgwick_update_marg(&mwf, gyro, accel, mag);

        /* Print every 100ms */
        if (step % 10 == 0) {
            MemsEuler e_comp = mems_quat_to_euler(q_comp);
            MemsEuler e_mah = mems_quat_to_euler(q_mahony);
            MemsEuler e_mad = mems_quat_to_euler(q_madgwick);

            printf("   %7.0f | (%+5.1f,%+5.1f)        | (%+5.1f,%+5.1f)        | (%+5.1f,%+5.1f)\n",
                   t * 1000.0,
                   mems_rad_to_deg(e_comp.roll),
                   mems_rad_to_deg(e_comp.pitch),
                   mems_rad_to_deg(e_mah.roll),
                   mems_rad_to_deg(e_mah.pitch),
                   mems_rad_to_deg(e_mad.roll),
                   mems_rad_to_deg(e_mad.pitch));
        }
    }

    /* ─── Final Results ─────────────────────────────────────────── */
    printf("\n3. Final Orientation Results (expected: roll=0°, pitch=90°)\n\n");

    MemsEuler e_comp_final = mems_quat_to_euler(q_comp);
    MemsEuler e_mah_final = mems_quat_to_euler(q_mahony);
    MemsEuler e_mad_final = mems_quat_to_euler(q_madgwick);

    printf("   Filter        | Roll (deg) | Pitch (deg) | Yaw (deg)\n");
    printf("   --------------|------------|-------------|----------\n");
    printf("   Complementary |   %+7.2f  |   %+7.2f   | %+7.2f\n",
           mems_rad_to_deg(e_comp_final.roll),
           mems_rad_to_deg(e_comp_final.pitch),
           mems_rad_to_deg(e_comp_final.yaw));
    printf("   Mahony        |   %+7.2f  |   %+7.2f   | %+7.2f\n",
           mems_rad_to_deg(e_mah_final.roll),
           mems_rad_to_deg(e_mah_final.pitch),
           mems_rad_to_deg(e_mah_final.yaw));
    printf("   Madgwick MARG |   %+7.2f  |   %+7.2f   | %+7.2f\n",
           mems_rad_to_deg(e_mad_final.roll),
           mems_rad_to_deg(e_mad_final.pitch),
           mems_rad_to_deg(e_mad_final.yaw));

    printf("\n   Note: Complementary filter cannot estimate yaw from accel alone.\n");
    printf("         Mahony & Madgwick MARG use magnetometer for yaw.\n");

    /* ─── Drift Analysis ─────────────────────────────────────────── */
    printf("\n4. Gyro-Only Drift Estimate (no fusion)\n\n");

    /* If we only integrate gyro, the orientation drifts due to bias */
    double gyro_bias = 0.01;  /* rad/s bias */
    double drift_1s = gyro_bias * 1.0;  /* 1 second */
    double drift_60s = gyro_bias * 60.0;  /* 1 minute */
    double drift_3600s = gyro_bias * 3600.0;  /* 1 hour */

    printf("   Gyro bias: %.3f rad/s (%.2f deg/s)\n",
           gyro_bias, mems_rad_to_deg(gyro_bias));
    printf("   Drift in 1 second:  %.2f deg\n", mems_rad_to_deg(drift_1s));
    printf("   Drift in 1 minute:  %.1f deg\n", mems_rad_to_deg(drift_60s));
    printf("   Drift in 1 hour:    %.1f deg\n\n", mems_rad_to_deg(drift_3600s));

    printf("   → This is why sensor fusion is essential!\n");
    printf("     Gyro alone drifts unboundedly.\n");
    printf("     Accelerometer corrects roll/pitch (gravity reference).\n");
    printf("     Magnetometer corrects yaw (magnetic north reference).\n");

    /* ─── Application Note ────────────────────────────────────────── */
    printf("\n5. Application Recommendations\n\n");
    printf("   | Use Case                | Filter      | Reason\n");
    printf("   |------------------------|-------------|----------------------\n");
    printf("   | Smartphone AR           | Madgwick    | Fast, 6/9-axis\n");
    printf("   | Drone FC (ArduPilot)    | Mahony/EKF  | Stable, drift-reject\n");
    printf("   | VR headset (fast motion)| Comp Filter | Low latency, simple\n");
    printf("   | Automotive ESC          | Kalman      | Optimal, linearized\n");
    printf("   | Human motion tracking   | Madgwick    | Magnetic immunity\n");

    printf("\n=== Sensor Fusion Demo Complete ===\n");
    return 0;
}
