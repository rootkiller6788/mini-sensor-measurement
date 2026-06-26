/**
 * ex3_ecg_frontend.c - ECG Biopotential Front-End Design
 *
 * L6/L7 Canonical Problem: ECG instrumentation amplifier front-end
 * design per IEC 60601-2-25 requirements.
 *
 * Demonstrates: gain selection, CMRR analysis, electrode offset
 * handling, and noise budget for medical ECG acquisition.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "inst_amp_defs.h"
#include "inst_amp_topology.h"
#include "inst_amp_analysis.h"
#include "inst_amp_sensorif.h"

int main(void) {
    printf("=== ECG Front-End Design Example ===\n\n");

    /* 1. ECG signal characteristics */
    double ecg_amplitude_mv = 2.0;     /* R-wave peak amplitude */
    double ecg_bw_low = 0.05;          /* Diagnostic low cutoff (Hz) */
    double ecg_bw_high = 150.0;        /* Diagnostic high cutoff (Hz) */
    double adc_fs = 3.3;               /* ADC full-scale (V) */

    /* 2. Determine required gain */
    double G = ecg_required_gain(ecg_amplitude_mv, adc_fs);
    printf("ECG amplitude: %.1f mV, ADC FS: %.1f V\n", ecg_amplitude_mv, adc_fs);
    printf("Required gain: %.0f V/V\n\n", G);

    /* 3. Electrode offset handling */
    ecg_frontend_config_t ecg_cfg = {
        .gain = G,
        .cmrr_db = 100.0,
        .highpass_hz = 0.05,
        .lowpass_hz = 150.0,
        .electrode_offset_mv = 300.0
    };

    double hp_att = ecg_highpass_attenuation(ecg_cfg.lowpass_hz,
                                              ecg_cfg.highpass_hz);
    printf("Highpass filter: f_c=%.2f Hz, attenuation at %.0f Hz = %.4f\n",
           ecg_cfg.highpass_hz, ecg_cfg.lowpass_hz, hp_att);
    printf("Electrode offset (%.0f mV) blocked by AC coupling\n\n",
           ecg_cfg.electrode_offset_mv);

    /* 4. Common-mode interference analysis */
    double v_cm_50hz = 1.0;   /* 1Vrms mains interference */
    double cm_noise_rti = ecg_cm_noise_rti_uv(v_cm_50hz, ecg_cfg.cmrr_db);
    printf("Common-mode interference: %.1f Vrms at 50 Hz\n", v_cm_50hz);
    printf("  CMRR: %.0f dB\n", ecg_cfg.cmrr_db);
    printf("  Residue RTI: %.2f uVrms\n", cm_noise_rti);
    printf("  Residue at output: %.2f mVrms\n\n", cm_noise_rti * G * 1e-3);

    /* 5. Noise budget */
    double en_amp = 10.0;   /* nV/rtHz */
    double in_amp = 1.0;    /* pA/rtHz */
    double rs = 5000.0;     /* Electrode-source impedance (Ohm) */

    double en_total = noise_total_rti_nv_per_rhz(en_amp, in_amp, rs, 0, 0, G, 27.0);
    double vn_rms_rti = en_total * sqrt(ecg_bw_high - ecg_bw_low) * 1e-3; /* uVrms */

    printf("Noise budget (0.05-150 Hz):\n");
    printf("  Amplifier en: %.1f nV/rtHz\n", en_amp);
    printf("  Total RTI noise density: %.1f nV/rtHz\n", en_total);
    printf("  Total RTI noise (BW=%.0f Hz): %.2f uVrms\n",
           ecg_bw_high - ecg_bw_low, vn_rms_rti);
    printf("  SNR at %.0f mV input: %.1f dB\n",
           ecg_amplitude_mv, signal_to_noise_ratio_db(ecg_amplitude_mv * 1e-3,
                                                       vn_rms_rti * 1e-6));
    printf("  Effective bits: %.1f\n\n",
           ia_effective_bits(ecg_amplitude_mv * 1e-3, vn_rms_rti * 1e-6));

    /* 6. Design summary */
    printf("=== Design Summary ===\n");
    printf("  Topology: 3-op-amp IA with AC coupling\n");
    printf("  Gain: %.0f V/V\n", G);
    printf("  Bandwidth: %.2f - %.0f Hz\n", ecg_cfg.highpass_hz, ecg_cfg.lowpass_hz);
    printf("  CMRR: %.0f dB minimum\n", ecg_cfg.cmrr_db);
    printf("  Input impedance: > 100 MOhm\n");
    printf("  Compliance: IEC 60601-2-25 diagnostic ECG\n");

    printf("\n=== ECG Front-End Design Complete ===\n");
    return 0;
}