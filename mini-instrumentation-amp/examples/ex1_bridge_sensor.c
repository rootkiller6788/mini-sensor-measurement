/**
 * ex1_bridge_sensor.c - Complete Bridge Sensor Readout Example
 *
 * L6 Canonical Problem: Wheatstone bridge strain gauge readout
 * from sensor to engineering units using a 3-op-amp IA.
 *
 * Demonstrates: bridge output -> IA gain selection -> ADC conversion
 * -> calibration correction -> strain in microstrain
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "inst_amp_defs.h"
#include "inst_amp_topology.h"
#include "inst_amp_calibration.h"
#include "inst_amp_sensorif.h"

int main(void) {
    printf("=== Bridge Sensor Readout Example ===\n\n");

    /* 1. Configure the strain gauge bridge */
    strain_gauge_config_t sg = {
        .gauge_factor = 2.0,        /* Constantan foil gauge */
        .max_strain_ue = 1000.0,    /* +/-1000 microstrain range */
        .v_excitation = 5.0,        /* 5V bridge excitation */
        .v_adc_fullscale = 5.0,     /* 5V ADC range */
        .bridge_r_ohm = 350.0       /* 350 Ohm bridge */
    };

    /* 2. Calculate expected bridge output at full scale */
    double v_bridge_fs_uv = strain_gauge_bridge_output_uv(&sg, sg.max_strain_ue);
    printf("Full-scale bridge output: %.2f uV (%.2f mV)\n",
           v_bridge_fs_uv, v_bridge_fs_uv * 1e-3);

    /* 3. Select IA gain to map sensor FS to ADC FS */
    double G = strain_gauge_required_gain(&sg);
    printf("Required IA gain: %.0f V/V\n", G);

    /* 4. Configure 3-op-amp IA */
    ia_three_opamp_config_t ia = {
        .r1_ohm = 25000.0,
        .r2_ohm = 10000.0,
        .r3_ohm = 10000.0,
        .rg_ohm = 0.0,
        .resistor_tol_pct = 0.1
    };
    ia.rg_ohm = ia_three_opamp_rg_for_gain(&ia, G);
    printf("R_G = %.2f Ohm -> standard value: %.0f Ohm\n",
           ia.rg_ohm, ia_e96_nearest(ia.rg_ohm));

    /* 5. Simulate a measurement at 250 microstrain */
    double strain_true = 250.0;
    double v_bridge_v = strain_gauge_bridge_output_uv(&sg, strain_true) * 1e-6;
    double v_ia_out = ia_three_opamp_output(&ia, v_bridge_v/2, -v_bridge_v/2, 0.0, 90.0);

    printf("\nAt %.0f uStrain:\n", strain_true);
    printf("  Bridge output: %.6f mV\n", v_bridge_v * 1000.0);
    printf("  IA output:     %.4f V\n", v_ia_out);

    /* 6. Simulate ADC reading and calibration */
    double v_adc = v_ia_out;  /* Ideal ADC */

    /* 7. Two-point calibration simulation */
    double v_cal1 = 0.0, v_out1 = 0.001;  /* Zero input, small offset */
    double v_cal2 = strain_gauge_bridge_output_uv(&sg, sg.max_strain_ue) * 1e-6;
    double v_out2 = ia_three_opamp_output(&ia, v_cal2/2, -v_cal2/2, 0.0, 90.0);

    double cal_gain, cal_offset;
    ia_cal_two_point(v_cal1, v_out1, v_cal2, v_out2, &cal_gain, &cal_offset);
    printf("\nCalibration: G_actual=%.4f, V_os_output=%.4f V\n", cal_gain, cal_offset);

    /* 8. Correct the measurement */
    double v_corrected = ia_cal_apply_correction(v_adc, cal_gain, cal_offset);
    double strain_measured = strain_gauge_adc_to_ue(v_adc, cal_gain, &sg);
    printf("Corrected input: %.6f V -> %.1f uStrain\n", v_corrected, strain_measured);
    printf("Error: %.2f uStrain (%.2f%%)\n", strain_measured - strain_true,
           100.0 * (strain_measured - strain_true) / strain_true);

    printf("\n=== Bridge Readout Complete ===\n");
    return 0;
}