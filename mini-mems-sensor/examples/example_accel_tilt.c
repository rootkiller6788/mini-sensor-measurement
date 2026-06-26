/**
 * @file example_accel_tilt.c
 * @brief L6 Example: Tilt sensing with a 3-axis MEMS accelerometer
 *
 * Demonstrates how a smartphone/tablet/game controller uses a MEMS
 * accelerometer to determine roll and pitch (device orientation).
 *
 * Key applications:
 *   - Smartphone screen rotation (Apple iPhone, Samsung Galaxy)
 *   - Game controller tilt (Nintendo Switch Joy-Con)
 *   - Digital camera level/horizon (Sony Alpha, Canon EOS)
 *   - Laptop free-fall detection
 *   - Construction digital level tool
 */

#include <stdio.h>
#include <math.h>
#include "mems_basics.h"
#include "mems_fusion.h"
#include "mems_dynamics.h"

#define G 9.80665

static const char *orientation_name(double roll_deg, double pitch_deg) {
    if (fabs(pitch_deg) < 15.0 && fabs(roll_deg) < 15.0)
        return "LANDSCAPE (phone flat, screen up)";
    if (fabs(pitch_deg) > 75.0)
        return "PORTRAIT UPRIGHT";
    if (fabs(pitch_deg + 90.0) < 15.0)
        return "LANDSCAPE LEFT (rotated CCW)";
    if (fabs(pitch_deg - 90.0) < 15.0)
        return "LANDSCAPE RIGHT (rotated CW)";
    if (fabs(roll_deg - 90.0) < 15.0)
        return "PORTRAIT LEFT (phone on left edge)";
    if (fabs(roll_deg + 90.0) < 15.0)
        return "PORTRAIT RIGHT (phone on right edge)";
    if (fabs(pitch_deg + 180.0) < 15.0 || fabs(pitch_deg - 180.0) < 15.0)
        return "FACE DOWN";
    return "TILTED";
}

int main(void) {
    printf("=== MEMS Accelerometer Tilt Sensing Example ===\n\n");

    /* Print the physical principle */
    printf("Physical Principle:\n");
    printf("  A 3-axis MEMS accelerometer measures the gravity vector\n");
    printf("  projected onto its sensing axes. From this, roll (φ) and\n");
    printf("  pitch (θ) are computed:\n\n");
    printf("    roll  = atan2(a_y, a_z)\n");
    printf("    pitch = atan2(-a_x, sqrt(a_y^2 + a_z^2))\n\n");

    /* ─── Scenario 1: Phone flat on table ─────────────────────────── */
    printf("─── Scenario 1: Phone Flat on Table (screen up) ───\n");
    MemsVector3 accel1 = {0.0, 0.0, G};
    MemsEuler e1 = mems_accel_to_roll_pitch(accel1);
    printf("  Accel reading:  (0.00, 0.00, %.2f) m/s^2\n", G);
    printf("  Roll:  %+.1f deg, Pitch: %+.1f deg\n",
           mems_rad_to_deg(e1.roll), mems_rad_to_deg(e1.pitch));
    printf("  Orientation: %s\n\n",
           orientation_name(mems_rad_to_deg(e1.roll),
                           mems_rad_to_deg(e1.pitch)));

    /* ─── Scenario 2: Phone tilted 30° ────────────────────────────── */
    printf("─── Scenario 2: Phone Tilted 30° (reading in landscape) ───\n");
    double tilt = 30.0 * M_PI / 180.0;
    MemsVector3 accel2 = {
        G * sin(tilt),
        0.0,
        G * cos(tilt)
    };
    MemsEuler e2 = mems_accel_to_roll_pitch(accel2);
    printf("  Accel reading:  (%.2f, 0.00, %.2f) m/s^2\n",
           accel2.x, accel2.z);
    printf("  Roll:  %+.1f deg, Pitch: %+.1f deg\n",
           mems_rad_to_deg(e2.roll), mems_rad_to_deg(e2.pitch));
    printf("  Orientation: %s\n\n",
           orientation_name(mems_rad_to_deg(e2.roll),
                           mems_rad_to_deg(e2.pitch)));

    /* ─── Scenario 3: Phone portrait upright ──────────────────────── */
    printf("─── Scenario 3: Phone Portrait Upright (screen facing user) ───\n");
    MemsVector3 accel3 = {0.0, G, 0.0};
    MemsEuler e3 = mems_accel_to_roll_pitch(accel3);
    printf("  Accel reading:  (0.00, %.2f, 0.00) m/s^2\n", G);
    printf("  Roll:  %+.1f deg, Pitch: %+.1f deg\n",
           mems_rad_to_deg(e3.roll), mems_rad_to_deg(e3.pitch));
    printf("  Orientation: %s\n\n",
           orientation_name(mems_rad_to_deg(e3.roll),
                           mems_rad_to_deg(e3.pitch)));

    /* ─── Scenario 4: Free-fall detection ─────────────────────────── */
    printf("─── Scenario 4: Free-Fall Detection ───\n");
    MemsVector3 accel4 = {0.0, 0.0, 0.0};  /* free fall */
    double mag = sqrt(accel4.x * accel4.x + accel4.y * accel4.y +
                       accel4.z * accel4.z);
    printf("  Accel reading:  (0.00, 0.00, 0.00) m/s^2\n");
    printf("  Acceleration magnitude: %.2f m/s^2 (%.2f g)\n", mag, mag / G);
    printf("  FREE FALL DETECTED! (|a| < 0.5g)\n");
    if (mag < 0.35 * G)
        printf("  → Disk head parking / laptop HDD protection triggered\n\n");

    /* ─── Scenario 5: Temperature effects on MEMS tilt ────────────── */
    printf("─── Scenario 5: MEMS Accelerometer Temperature Drift ───\n");
    printf("  (Simulating Bosch BMI270 temperature behavior)\n\n");

    MemsSpec accel_spec;
    mems_spec_typical_accel(&accel_spec);

    printf("  Spec: Noise density = %.0f ug/√Hz\n", accel_spec.noise_density);
    printf("  Spec: Bias tempco = %.2f mg/°C\n", accel_spec.bias_tempco);
    printf("  Spec: Bandwidth = %.0f Hz\n", accel_spec.bandwidth);

    /* Noise in 200 Hz bandwidth */
    double noise_rms = accel_spec.noise_density * 1e-6 * G *
                        sqrt(accel_spec.bandwidth);
    printf("  RMS noise at %d Hz BW: %.3f mg\n",
           (int)accel_spec.bandwidth,
           noise_rms / G * 1000.0);

    /* Temperature error over 50°C range */
    double temp_drift = accel_spec.bias_tempco * 50.0; /* mg */
    double tilt_error_rad = asin(temp_drift / 1000.0);
    printf("  Bias drift over 50°C: %.1f mg → tilt error: %.2f deg\n",
           temp_drift, mems_rad_to_deg(tilt_error_rad));
    printf("  This is the MEMS temperature-dependent tilt error budget.\n\n");

    printf("=== MEMS Tilt Sensing Complete ===");

    return 0;
}
