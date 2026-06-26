/**
 * @file ex_accelerometer_sim.c
 * @brief Example: MEMS accelerometer dynamic response simulation
 *
 * This example demonstrates:
 *   1. 2nd-order mass-spring-damper model of a MEMS accelerometer
 *   2. Frequency response (Bode plot data)
 *   3. Step response with overshoot and settling time
 *   4. Bandwidth and resonant behavior
 *   5. Noise simulation (thermal + flicker)
 *
 * This addresses L6 canonical problem: inertial sensor dynamic characterization.
 * MEMS accelerometers are used in GPS/IMU (Boeing, Airbus), automotive airbag
 * (Detroit), smartphone orientation (Apple iPhone), and quadrotor stabilization.
 *
 * References:
 *   - Yazdi, "Micromachined Inertial Sensors", Proc. IEEE 86(8), 1998
 *   - Titterton & Weston, "Strapdown Inertial Navigation Technology", 2nd Ed
 */

#include <stdio.h>
#include <math.h>
#include "../include/sensor_transfer.h"
#include "../include/sensor_noise.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== MEMS Accelerometer Dynamic Simulation ===\n\n");

    /* ──── Accelerometer Model Parameters ────
     * Typical MEMS accelerometer:
     *   - Natural frequency: 1-10 kHz (limited by proof mass and spring)
     *   - Damping ratio: 0.3-0.7 (underdamped for fast response)
     *   - Sensitivity: 100-1000 mV/g
     *   - Noise floor: 10-100 µg/√Hz
     *   - Full scale: ±2g to ±200g
     */
    sensor_2nd_order_t accel;
    double K_sens = 0.400;       /* 400 mV/g sensitivity */
    double f_n = 2000.0;         /* 2 kHz natural frequency */
    double omega_n = 2.0 * M_PI * f_n;
    double zeta = 0.5;           /* damping ratio (slightly underdamped) */
    sensor_2nd_order_init(&accel, K_sens, omega_n, zeta, 0.0, 0.0);

    printf("Accelerometer Parameters:\n");
    printf("  Sensitivity:     %.1f mV/g\n", K_sens * 1000.0);
    printf("  Natural freq:    %.0f Hz (ω_n = %.0f rad/s)\n", f_n, omega_n);
    printf("  Damping ratio:   %.2f\n", zeta);
    printf("  Bandwidth (-3dB): %.0f Hz\n",
           sensor_2nd_order_bandwidth(&accel));
    if (zeta < 0.707) {
        printf("  Resonant freq:   %.0f Hz\n",
               sensor_2nd_order_resonant_frequency(&accel));
        printf("  Q-factor:        %.1f\n",
               sensor_2nd_order_resonant_magnification(&accel));
    }

    /* ──── Frequency Response ──── */
    printf("\n┌─── Frequency Response (Bode Magnitude) ──────────────────────┐\n");
    printf("│  Freq (Hz)   |H(f)| (mV/g)    |H| (dB rel DC)    Phase (°)  │\n");
    printf("│  ─────────   ──────────────    ───────────────    ─────────  │\n");

    double freqs[] = {10, 50, 100, 500, 1000, 2000, 3000, 5000};
    for (int i = 0; i < 8; i++) {
        double f = freqs[i];
        double mag = sensor_2nd_order_freq_response_mag(&accel, f);
        double phase = sensor_2nd_order_freq_response_phase(&accel, f);
        double mag_db = (K_sens > 0) ? 20.0 * log10(mag / K_sens) : -99.0;

        printf("│  %7.0f      %9.3f          %+7.2f          %+8.2f     │\n",
               f, mag * 1000.0, mag_db, phase * 180.0 / M_PI);
    }
    printf("└──────────────────────────────────────────────────────────────┘\n");

    /* ──── Step Response — 1g step ──── */
    printf("\n┌─── Step Response to 1g Acceleration ─────────────────────────┐\n");
    printf("│  Time (ms)    Output (mV)    %% of Final                    │\n");
    printf("│  ─────────    ───────────    ──────────                     │\n");

    double step_g = 1.0;
    double t_peak = sensor_2nd_order_peak_time(&accel);
    double Mp = sensor_2nd_order_overshoot_percent(&accel);
    double t_settle = sensor_2nd_order_settling_time(&accel, 2.0);

    double t_samples[] = {0.0, t_peak/3, t_peak*2/3, t_peak, t_peak*2,
                          t_settle/2, t_settle, t_settle*2};
    char *labels[] = {"0", "tp/3", "2tp/3", "tp (peak)", "2tp",
                      "ts/2", "ts (settle)", "2ts"};

    for (int i = 0; i < 8; i++) {
        double y = sensor_2nd_order_step_response(&accel, step_g, t_samples[i]);
        double pct = y / (K_sens * step_g) * 100.0;
        if (i >= 4) {
            printf("│  %-11s  %9.3f      %7.1f%%                    │\n",
                   labels[i], y * 1000.0, pct);
        } else {
            printf("│  %-11s  %9.3f      %7.1f%%                    │\n",
                   labels[i], y * 1000.0, pct);
        }
    }
    printf("│                                                            │\n");
    printf("│  Peak time: %.3f ms     Overshoot: %.1f%%                    │\n",
           t_peak * 1000.0, Mp);
    printf("│  Settling time (2%%): %.3f ms                               │\n",
           t_settle * 1000.0);
    printf("└──────────────────────────────────────────────────────────────┘\n");

    /* ──── Noise Floor Analysis ──── */
    printf("\n┌─── Noise Floor Analysis ─────────────────────────────────────┐\n");

    double temp_K = 300.0;
    double bw_hz = 100.0;  /* measurement bandwidth */
    /* Assume MEMS has equivalent resistance noise from readout electronics */
    double R_equiv = 100e3;  /* feedback resistor in transimpedance amp */
    double vn_thermal = johnson_noise_voltage(temp_K, R_equiv, bw_hz);
    double vn_per_root_hz = vn_thermal / sqrt(bw_hz);

    printf("│  Temperature: %.0f K (%.0f °C)                              │\n",
           temp_K, temp_K - 273.15);
    printf("│  Bandwidth: %.0f Hz                                          │\n", bw_hz);
    printf("│  Johnson noise: %.2f µVrms (%.1f nV/√Hz)                    │\n",
           vn_thermal * 1e6, vn_per_root_hz * 1e9);

    /* Noise-equivalent acceleration (in g) */
    double nea = vn_thermal / (K_sens / 1000.0);  /* in mg */
    printf("│  Noise-equivalent acceleration: %.3f mg (%.1f µg/√Hz)       │\n",
           nea, nea / sqrt(bw_hz) * 1000.0);

    double resolution = resolution_from_noise(vn_thermal, K_sens);
    printf("│  Resolution (3σ): %.3f mg                                    │\n",
           resolution * 1000.0);
    printf("└──────────────────────────────────────────────────────────────┘\n");

    /* ──── Ramp Response (constant acceleration test) ──── */
    printf("\n┌─── Dynamic Error for Ramp Input ─────────────────────────────┐\n");
    double ramp_rate = 1.0;  /* 1 g/s acceleration ramp */
    double e_ramp = sensor_dynamic_ramp_error_2nd_order(ramp_rate,
                                                         omega_n, zeta);
    printf("│  Ramp rate: %.1f g/s                                         │\n", ramp_rate);
    printf("│  Steady-state tracking error: %.4f g                         │\n", e_ramp);
    printf("│  → Faster sensor (higher ω_n) reduces tracking error        │\n");
    printf("└──────────────────────────────────────────────────────────────┘\n");

    return 0;
}
