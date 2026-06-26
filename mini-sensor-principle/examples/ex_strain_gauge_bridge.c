/**
 * @file ex_strain_gauge_bridge.c
 * @brief Example: Strain gauge measurement with Wheatstone bridge
 *
 * This example demonstrates:
 *   1. Quarter, half, and full bridge strain measurement
 *   2. Bridge nonlinearity and its correction
 *   3. Effect of lead wire resistance in 2-wire measurement
 *   4. Amplifier gain and noise considerations
 *   5. Strain-to-stress conversion (Hooke's Law)
 *
 * This is a canonical L6 problem: strain measurement for structural
 * health monitoring, as used in Boeing 747, Apollo, F-35, and
 * Detroit automotive testing.
 */

#include <stdio.h>
#include <math.h>
#include "../include/sensor_bridge.h"

int main(void) {
    printf("=== Strain Gauge Bridge Measurement ===\n\n");

    const double R0 = 350.0;       /* 350Ω strain gauge (common for steel) */
    const double V_ex = 5.0;       /* 5V bridge excitation */
    const double GF = 2.06;        /* Constantan gauge factor */
    const double E_mod = 200e9;    /* Young's modulus for steel [Pa] */
    const double R_lead = 0.5;     /* Lead wire resistance per conductor */

    /* ──── Quarter Bridge Example ──── */
    printf("┌─── Quarter Bridge (1 active gauge) ───────────────────────────┐\n");

    wheatstone_bridge_t qb;
    bridge_init_balanced(&qb, R0, V_ex);

    double strains[] = {100e-6, 500e-6, 1000e-6, 2000e-6, 5000e-6}; /* µε */
    printf("│  Strain (µε)   ΔR (Ω)    V_out (mV)   Linear ε    Error %%   │\n");
    printf("│  ───────────   ──────    ──────────   ─────────    ──────   │\n");

    for (int i = 0; i < 5; i++) {
        double strain = strains[i];
        double dR = R0 * GF * strain;
        bridge_set_quarter_active(&qb, R0, R0 + dR, 1);

        double V_out = bridge_output_voltage_exact(&qb);
        double strain_meas = strain_from_bridge_output(V_out, V_ex, GF,
                                                        BRIDGE_QUARTER);
        /* Nonlinearity correction would be applied in precision applications */
        double error_pct = (strain_meas - strain) / strain * 100.0;

        printf("│  %8.0f     %6.3f     %10.4f    %8.0f      %+6.2f    │\n",
               strain * 1e6, dR, V_out * 1000.0, strain_meas * 1e6, error_pct);
    }
    printf("└───────────────────────────────────────────────────────────────┘\n");

    /* ──── Half Bridge vs Quarter Bridge ──── */
    printf("\n┌─── Bridge Type Comparison (ε = 1000 µε) ──────────────────────┐\n");
    printf("│  Type          V_out (mV)    Sensitivity (mV/V/ε)             │\n");
    printf("│  ──────        ──────────    ────────────────────             │\n");

    /* Quarter bridge */
    bridge_set_quarter_active(&qb, R0, R0 + R0 * GF * 0.001, 1);
    double Vq = bridge_output_voltage_exact(&qb);
    printf("│  Quarter         %8.4f        %7.4f                      │\n",
           Vq * 1000.0, (Vq / V_ex / 0.001) * 1000.0);

    /* Half bridge (adjacent, temperature compensated) */
    wheatstone_bridge_t hb;
    bridge_init_balanced(&hb, R0, V_ex);
    bridge_set_half_adjacent(&hb, R0, R0 + R0 * GF * 0.001, R0 - R0 * GF * 0.001);
    double Vh = bridge_output_voltage_exact(&hb);
    printf("│  Half-Adj        %8.4f        %7.4f                      │\n",
           Vh * 1000.0, (Vh / V_ex / 0.001) * 1000.0);

    /* Full bridge */
    wheatstone_bridge_t fb;
    bridge_init_balanced(&fb, R0, V_ex);
    bridge_set_full_pushpull(&fb, R0, R0 * GF * 0.001);
    double Vf = bridge_output_voltage_exact(&fb);
    printf("│  Full            %8.4f        %7.4f                      │\n",
           Vf * 1000.0, (Vf / V_ex / 0.001) * 1000.0);

    printf("│  Full/Quarter sensitivity ratio: %.1fx                       │\n",
           Vf / Vq);
    printf("└───────────────────────────────────────────────────────────────┘\n");

    /* ──── Lead Wire Error ──── */
    printf("\n┌─── Lead Wire Error Analysis ──────────────────────────────────┐\n");
    double err_2w = bridge_lead_wire_error_2wire(R0, R_lead);
    printf("│  R_gauge = %.0f Ω, R_lead/wire = %.1f Ω                      │\n",
           R0, R_lead);
    printf("│  2-wire fractional error: %.3f%%                              │\n",
           err_2w * 100.0);
    printf("│  Equivalent strain error: %.0f µε                             │\n",
           err_2w / GF * 1e6);
    printf("│  → Use 3-wire or 4-wire (Kelvin) connection for precision     │\n");
    printf("└───────────────────────────────────────────────────────────────┘\n");

    /* ──── Signal Conditioning Requirements ──── */
    printf("\n┌─── Amplifier Requirements (for 10-bit ADC, Vref=5.0V) ───────┐\n");
    double V_min = Vq;  /* quarter bridge at 1000 µε */
    double V_adc_fs = 5.0;
    double gain = bridge_amplifier_gain_required(V_min, V_adc_fs);
    printf("│  Min differential signal: %.3f mV                             │\n",
           V_min * 1000.0);
    printf("│  Required gain: %.0f V/V (to use full ADC range)              │\n",
           gain);
    double cmrr = bridge_amplifier_cmrr_required(V_ex, V_min, 0.001);
    printf("│  Required CMRR: %.0f (%.0f dB) for 0.1%% error                 │\n",
           cmrr, 20.0 * log10(cmrr));
    printf("└───────────────────────────────────────────────────────────────┘\n");

    /* ──── Stress Calculation (Hooke's Law) ──── */
    printf("\n┌─── Strain-to-Stress (Steel, E=200 GPa) ───────────────────────┐\n");
    printf("│  ε (µε)    Stress (MPa)    Application                        │\n");
    printf("│  ──────    ────────────    ───────────                        │\n");
    printf("│    100       %6.1f        Normal service                     │\n",
           100e-6 * E_mod / 1e6);
    printf("│   1000       %6.1f        Design limit                       │\n",
           1000e-6 * E_mod / 1e6);
    printf("│   2000       %6.1f        Yield point (mild steel)           │\n",
           2000e-6 * E_mod / 1e6);
    printf("└───────────────────────────────────────────────────────────────┘\n");

    return 0;
}
