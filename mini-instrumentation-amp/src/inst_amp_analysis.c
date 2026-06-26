/**
 * inst_amp_analysis.c - IA Noise, CMRR, and Error Budget Analysis
 *
 * L3 Mathematical Structures: Complete noise analysis (thermal,
 * shot, 1/f), CMRR modeling, error budget synthesis, and
 * sensitivity analysis for instrumentation amplifier designs.
 *
 * L4 Fundamental Laws:
 *   - Johnson-Nyquist noise: Vn = sqrt(4kTR * BW)
 *   - Shot noise: In = sqrt(2qI * BW)
 *   - 1/f noise: Vn^2 = K/f * BW
 *   - RSS error summation: err_total = sqrt(sum err_i^2)
 *
 * Courses: Stanford EE247, MIT 6.630, Berkeley EE105
 * Reference:
 *   Motchenbacher & Connelly, "Low-Noise Electronic System Design" (1993)
 *   Gray, Hurst, Lewis, Meyer, "Analysis and Design of Analog ICs" (5th ed., 2009)
 *   Kitchin & Counts (2006), Ch. 6 - Error Budget Analysis
 */

#include "inst_amp_defs.h"
#include "inst_amp_topology.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ==================================================================
 * L3/L4: Johnson-Nyquist (Thermal) Noise
 *
 * Physical origin: random thermal motion of charge carriers in
 * any resistive material above absolute zero (Johnson, 1928;
 * Nyquist, 1928).
 *
 * Theorem (Nyquist, 1928):
 *   v_n_rms = sqrt(4 * k_B * T * R * BW)
 *   where k_B = 1.380649e-23 J/K (Boltzmann constant)
 *         T = absolute temperature (K)
 *         R = resistance (Ohm)
 *         BW = noise bandwidth (Hz)
 *
 * At room temperature (27C = 300K):
 *   v_n = sqrt(4 * 1.38e-23 * 300 * R * BW)
 *       ~ 0.129 * sqrt(R) nV/rtHz
 *   For R = 1k:   en = 4.07 nV/rtHz
 *   For R = 10k:  en = 12.9 nV/rtHz
 *   For R = 100k: en = 40.7 nV/rtHz
 *
 * This is the fundamental noise floor that no amplifier can beat.
 * For a 350-Ohm strain gauge bridge at 27C: en = 2.4 nV/rtHz.
 *
 * Complexity: O(1)
 * ================================================================== */
double noise_thermal_nv_per_rhz(double r_ohm, double temp_c) {
    if (r_ohm <= 0.0) return 0.0;

    double k_B = 1.380649e-23;
    double T_K = temp_c + 273.15;

    /* v_n = sqrt(4*k*T*R) in V/rtHz, convert to nV/rtHz */
    return sqrt(4.0 * k_B * T_K * r_ohm) * 1e9;
}

double noise_thermal_rms_uv(double r_ohm, double bw_hz, double temp_c) {
    return noise_thermal_nv_per_rhz(r_ohm, temp_c) * sqrt(bw_hz) * 1e-3;
}

/* ==================================================================
 * L3/L4: Shot Noise
 *
 * Physical origin: discrete nature of charge carriers crossing
 * a potential barrier (p-n junction). Each carrier arrival is
 * a Poisson process (Schottky, 1918).
 *
 * Theorem (Schottky, 1918):
 *   i_n_rms = sqrt(2 * q * I_DC * BW)
 *   where q = 1.602176634e-19 C (elementary charge)
 *         I_DC = DC bias current (A)
 *         BW = noise bandwidth (Hz)
 *
 * In bipolar transistor input IAs: I_B creates shot noise.
 * In CMOS input IAs: gate leakage < 1pA, shot noise negligible.
 *
 * For I_B = 500 pA (AD620): i_n = sqrt(2*1.6e-19*5e-10) = 0.4 fA/rtHz
 * For I_B = 1 nA (AD8421):    i_n = sqrt(2*1.6e-19*1e-9) = 0.57 fA/rtHz
 *
 * Complexity: O(1)
 * ================================================================== */
double noise_shot_fa_per_rhz(double i_dc_a) {
    if (i_dc_a <= 0.0) return 0.0;

    double q = 1.602176634e-19;
    /* i_n = sqrt(2*q*I_DC) in A/rtHz, convert to fA/rtHz */
    return sqrt(2.0 * q * i_dc_a) * 1e15;
}

/* ==================================================================
 * L3/L4: 1/f Noise (Flicker Noise)
 *
 * Physical origin: trap states in semiconductor-oxide interface
 * and in bulk material capture and release carriers with a
 * distribution of time constants, producing a 1/f spectrum
 * (Hooge, 1969; McWhorter, 1957).
 *
 * Empirical model:
 *   v_n_1f^2(f) = K_f / f^alpha  (alpha ~ 0.8-1.2, typically 1)
 *   v_n_1f^2(BW) = K_f * ln(f_high / f_low)
 *
 * Corner frequency f_c: frequency where 1/f noise power equals
 * thermal noise power:
 *   en_1f(f_c) = en_thermal
 *   i.e., sqrt(K_f / f_c) = en_thermal
 *   -> f_c = K_f / en_thermal^2
 *
 * Below f_c, 1/f noise dominates. Above f_c, thermal noise dominates.
 *
 * For a typical precision op-amp: f_c = 1-100 Hz
 * For chopper-stabilized amps: f_c = 0.01 Hz (virtually no 1/f!)
 *
 * Reference: Hooge, "1/f noise is no surface effect", Phys. Lett. A, 1969
 * Complexity: O(1)
 * ================================================================== */
double noise_one_over_f_nv_rms(double k_f_nv2, double f_low_hz, double f_high_hz) {
    if (f_low_hz <= 0.0 || f_high_hz <= f_low_hz) return 0.0;
    return sqrt(k_f_nv2 * log(f_high_hz / f_low_hz));
}

double noise_corner_frequency(double k_f_nv2, double en_thermal_nv2_per_hz) {
    if (en_thermal_nv2_per_hz <= 0.0) return INFINITY;
    return k_f_nv2 / en_thermal_nv2_per_hz;
}

/* ==================================================================
 * L3: Total Amplifier Input-Referred Noise (Complete Model)
 *
 * Combines all noise sources referred to input:
 *   1. Amplifier voltage noise (en_amp): from input stage
 *   2. Amplifier current noise (in_amp): flows through source R_s
 *   3. Source thermal noise: 4kTR_s
 *   4. Gain-setting resistor noise: for 3-op-amp, R_G contributes
 *      v_n_RG = sqrt(4kTR_G) * 2R1/R_G (referred to input)
 *
 * Total:
 *   en_total = sqrt(en_amp^2 + (in_amp*R_s)^2 + 4kTR_s +
 *                   4kTR_G*(2R1/R_G)^2/ G^2)
 *
 * Complexity: O(1)
 * ================================================================== */
double noise_total_rti_nv_per_rhz(double en_amp_nv_rhz, double in_amp_pa_rhz,
                                   double rs_ohm, double rg_ohm,
                                   double r1_ohm, double gain,
                                   double temp_c) {
    double k_B = 1.380649e-23;
    double T_K = temp_c + 273.15;

    /* Amplifier voltage noise */
    double en_amp2 = en_amp_nv_rhz * en_amp_nv_rhz;

    /* Amplifier current noise through source R */
    double in_amp_nv_rhz = (in_amp_pa_rhz * 1e-12) * rs_ohm * 1e9;
    double en_curr2 = in_amp_nv_rhz * in_amp_nv_rhz;

    /* Source thermal noise */
    double en_src = 4.0 * k_B * T_K * rs_ohm;
    double en_src_nv2 = en_src * 1e18;  /* Convert V^2/Hz to nV^2/Hz */

    /* R_G noise, referred to input */
    double en_rg2 = 0.0;
    if (rg_ohm > 0.0 && gain > 0.0) {
        double en_rg = 4.0 * k_B * T_K * rg_ohm;
        double en_rg_rti = en_rg * (2.0 * r1_ohm / rg_ohm) / gain;
        en_rg2 = en_rg_rti * 1e18;
    }

    return sqrt(en_amp2 + en_curr2 + en_src_nv2 + en_rg2);
}

/* ==================================================================
 * L3: Noise Figure (NF) Computation
 *
 * Noise Figure quantifies how much an amplifier degrades the SNR:
 *   NF_dB = 10 * log10(SNR_in / SNR_out)
 *         = 10 * log10(1 + en_amp_total^2 / en_source^2)
 *
 * Where en_source^2 = 4kTR_s (thermal noise of source).
 *
 * For an ideal (noiseless) amplifier: NF = 0 dB.
 * For AD620 with R_s = 1k: NF ~ 2 dB (G=100)
 *
 * Cascaded noise figure (Friis formula):
 *   NF_total = NF1 + (NF2-1)/G1 + (NF3-1)/(G1*G2) + ...
 *
 * Reference: Friis, "Noise Figures of Radio Receivers", Proc. IRE, 1944
 * Complexity: O(1)
 * ================================================================== */
double noise_figure_db(double en_amp_total_nv_rhz, double rs_ohm, double temp_c) {
    if (rs_ohm <= 0.0) return 0.0;

    double en_source = noise_thermal_nv_per_rhz(rs_ohm, temp_c);

    if (en_source <= 0.0) return 0.0;
    return 10.0 * log10(1.0 + (en_amp_total_nv_rhz * en_amp_total_nv_rhz)
                              / (en_source * en_source));
}

/* ==================================================================
 * L3: Cascaded Noise Figure (Friis Formula)
 *
 * NF_cascade_dB = NF1_dB + (NF2_linear - 1)/G1_linear + ...
 *
 * Complexity: O(n_stages)
 * ================================================================== */
double noise_figure_cascade_db(const double *nf_db, const double *gain_linear,
                                int n_stages) {
    if (!nf_db || !gain_linear || n_stages <= 0) return 0.0;

    double nf_linear = 0.0;
    double cum_gain = 1.0;

    for (int i = 0; i < n_stages; i++) {
        double nf_i_linear = pow(10.0, nf_db[i] / 10.0);
        nf_linear += (nf_i_linear - 1.0) / cum_gain;
        cum_gain *= gain_linear[i];
    }

    /* Add 1 for the noise contribution that would exist even without stages */
    /* The Friis formula sums excess noise normalized by preceding gain */
    double nf_total_linear = 1.0 + nf_linear;
    return 10.0 * log10(nf_total_linear);
}

/* ==================================================================
 * L3: Combined DC Error Budget (RSS Method)
 *
 * All independent error sources sum in root-sum-square (RSS):
 *   V_err_total_rti = sqrt(V_os^2 + (I_os*R_s)^2
 *                         + (V_cm/CMRR_linear)^2
 *                         + (V_psrr*V_supply_ripple)^2
 *                         + ...)
 *
 * This assumes errors are uncorrelated Gaussian distributions.
 * For worst-case analysis, use simple sum (non-RSS).
 *
 * Each error term is first "referred to input" (RTI), then the
 * total RTI error is multiplied by gain for output error.
 *
 * Complexity: O(1)
 * ================================================================== */
double error_budget_rti_uv(const ia_spec_t *spec, double rs_ohm,
                            double v_cm_v, double v_ripple_v) {
    if (!spec) return 0.0;

    double total_var = 0.0;

    /* 1. Input offset voltage */
    total_var += spec->vos_uv * spec->vos_uv;

    /* 2. Input offset current through source resistance */
    double vos_from_ios = (spec->ios_pa * 1e-12) * rs_ohm * 1e6;
    total_var += vos_from_ios * vos_from_ios;

    /* 3. CMRR error: V_cm rejected by CMRR, residue appears RTI */
    double cmrr_linear = pow(10.0, spec->cmrr_db / 20.0);
    if (cmrr_linear > 0.0) {
        double v_cm_error_rti = v_cm_v / cmrr_linear * 1e6;  /* uV */
        total_var += v_cm_error_rti * v_cm_error_rti;
    }

    /* 4. PSRR error: supply ripple attenuated by PSRR, residue RTI */
    double psrr_linear = pow(10.0, spec->psrr_db / 20.0);
    if (psrr_linear > 0.0) {
        double v_psrr_error_rti = v_ripple_v / psrr_linear * 1e6;
        total_var += v_psrr_error_rti * v_psrr_error_rti;
    }

    return sqrt(total_var);
}

/* ==================================================================
 * L3: Gain Error Budget
 *
 * Total gain error combines:
 *   - Resistor tolerance (dominant for 3-op-amp with R_G)
 *   - Op-amp finite open-loop gain: G_err = G/(A_OL) approx
 *   - CMRR effect on gain (small, second-order)
 *
 * Total gain error (RSS):
 *   G_err_total = sqrt(G_err_RG^2 + G_err_AOL^2)
 *
 * Complexity: O(1)
 * ================================================================== */
double gain_error_budget_pct(double rg_tol_pct, double G, double A_OL) {
    double err_rg = rg_tol_pct;  /* R_G tolerance directly affects gain */

    /* Finite A_OL creates gain error: G_actual = G_ideal / (1 + G_ideal/A_OL) */
    double err_aol = 0.0;
    if (A_OL > 0.0) {
        double g_actual = G / (1.0 + G / A_OL);
        err_aol = fabs(g_actual - G) / G * 100.0;
    }

    return sqrt(err_rg * err_rg + err_aol * err_aol);
}

/* ==================================================================
 * L3: Signal-to-Noise Ratio (SNR) at IA Output
 *
 * SNR_dB = 20 * log10(V_signal_rms / V_noise_rms)
 *
 * Minimum detectable signal (MDS): signal level where SNR = 0 dB
 *   MDS = V_noise_rms  (input-referred)
 *
 * Dynamic range: DR_dB = 20 * log10(V_max / V_noise_floor)
 *
 * For a 16-bit ADC with 10V FS: LSB = 10/65536 = 153 uV
 * For meaningful 1-LSB resolution, V_noise_RTI * G < 153 uV
 *
 * Complexity: O(1)
 * ================================================================== */
double signal_to_noise_ratio_db(double v_signal_rms, double v_noise_rms) {
    if (v_noise_rms <= 0.0) return INFINITY;
    return 20.0 * log10(v_signal_rms / v_noise_rms);
}

double minimum_detectable_signal(double v_noise_rms) {
    return v_noise_rms;  /* SNR = 0 dB => signal = noise */
}

double dynamic_range_db(double v_max, double v_noise_floor) {
    if (v_noise_floor <= 0.0) return INFINITY;
    return 20.0 * log10(v_max / v_noise_floor);
}

/* ==================================================================
 * L3: CMRR as Function of Source Impedance Imbalance
 *
 * Even a perfect IA has finite CMRR when source impedances are
 * unbalanced. The imbalance converts common-mode voltage to
 * differential:
 *
 *   V_diff_error = V_cm * (deltaZ / Z_in_cm)
 *
 * Effective CMRR from impedance imbalance:
 *   CMRR_Z = Z_in_cm / deltaZ  (linear)
 *
 * Example: Z_in_cm = 10 GOhm, deltaZ = 1 kOhm (from electrode mismatch)
 *   CMRR_Z = 10e9 / 1000 = 10e6 = 140 dB
 *
 * But if Z_in_cm = 10 MOhm (2-op-amp topology), deltaZ = 1 kOhm:
 *   CMRR_Z = 10e6 / 1000 = 1000 = 60 dB  (much worse!)
 *
 * This is why high input impedance is critical for CMRR.
 *
 * Reference: Metting van Rijn et al., "The Isolation Mode Rejection
 *   Ratio in Bioelectric Amplifiers", IEEE BME, 1991
 * Complexity: O(1)
 * ================================================================== */
double cmrr_from_impedance_imbalance(double zin_cm_ohm, double delta_z_ohm) {
    if (delta_z_ohm <= 0.0) return INFINITY;
    double cmrr_linear = zin_cm_ohm / delta_z_ohm;
    return 20.0 * log10(cmrr_linear);
}

/* ==================================================================
 * L3: Combined CMRR from Multiple Sources
 *
 * Total CMRR from op-amp CMRR, resistor mismatch, and impedance
 * imbalance:
 *   1/CMRR_total = 1/CMRR_amp + 1/CMRR_res + 1/CMRR_Z
 *
 * All CMRR values must be in linear units.
 *
 * Complexity: O(1)
 * ================================================================== */
double cmrr_combined_db(double cmrr_amp_db, double cmrr_res_db,
                         double cmrr_impedance_db) {
    double c1 = pow(10.0, -cmrr_amp_db / 20.0);
    double c2 = pow(10.0, -cmrr_res_db / 20.0);
    double c3 = pow(10.0, -cmrr_impedance_db / 20.0);

    double c_total_inv = c1 + c2 + c3;
    if (c_total_inv <= 0.0) return INFINITY;
    return -20.0 * log10(c_total_inv);
}

/* ==================================================================
 * L3: Sensitivity Analysis - Effect of R_G Tolerance on Gain
 *
 * dG/G = - (G-1)/G * dR_G/R_G  (for 3-op-amp with R2=R3)
 *
 * For G=100: dG/G = -0.99 * dR_G/R_G
 *   1% R_G error -> 0.99% gain error  (nearly 1:1)
 * For G=2:   dG/G = -0.5 * dR_G/R_G
 *   1% R_G error -> 0.5% gain error
 *
 * Sensitivity diminishes at low gains.
 *
 * Complexity: O(1)
 * ================================================================== */
double gain_sensitivity_to_rg(double G) {
    if (G <= 1.0) return 0.0;
    return -(G - 1.0) / G;  /* dG/G per dR_G/R_G */
}

/* ==================================================================
 * L3: Sensitivity Analysis - Effect of Internal Resistor Drift
 *
 * For 3-op-amp IA with internal thin-film resistors:
 *   dG/G = (2R1/R_G)/G * dR1/R1 + (2R1/R_G)/G * dR_G/R_G
 *
 * The internal R1 has tempco of +/-50 ppm/C (thin-film SiCr).
 * For delta_T = 50C: dR1/R1 = 2500 ppm = 0.25%
 *
 * Complexity: O(1)
 * ================================================================== */
double gain_tempco_ppm_per_c(double G, double r1_ppm, double rg_ppm) {
    if (G <= 1.0) return 0.0;

    /* Sensitivity coefficients */
    double s_r1 = 2.0 * (G - 1.0) / G;
    double s_rg = -2.0 * (G - 1.0) / G;

    return sqrt((s_r1 * r1_ppm) * (s_r1 * r1_ppm)
              + (s_rg * rg_ppm) * (s_rg * rg_ppm));
}

/* ==================================================================
 * L3: Offset Voltage Drift over Temperature
 *
 * V_os(T) = V_os(25C) + TC_Vos * (T - 25)
 *
 * Total drift over a temperature range:
 *   delta_Vos = TC_Vos * (T_max - T_min)
 *
 * For AD620: TC_Vos = 0.5 uV/C (G>=100), delta_T = 100C -> 50 uV drift
 * For chopper amp: TC_Vos = 0.005 uV/C -> negligible
 *
 * Complexity: O(1)
 * ================================================================== */
double offset_drift_uv(double vos_25c_uv, double tc_vos_nv_per_c,
                        double temp_c) {
    return vos_25c_uv + tc_vos_nv_per_c * (temp_c - 25.0) * 1e-3;
}

/* ==================================================================
 * L3: Input Bias Current Drift and Its Effect
 *
 * Input bias current flows through the source resistance, creating
 * a voltage drop that appears as an additional offset:
 *   V_os_ib = I_B * R_s
 *
 * For I_B = 500 pA and R_s = 10 kOhm: V_os_ib = 5 uV
 * For I_B = 500 pA and R_s = 1 MOhm: V_os_ib = 500 uV (significant!)
 *
 * For high source impedance sensors (pH electrodes, EEG),
 * ultra-low I_B amplifiers (AD549: I_B ~ 40 fA) are required.
 *
 * Complexity: O(1)
 * ================================================================== */
double bias_current_offset_uv(double ib_pa, double rs_ohm) {
    return (ib_pa * 1e-12) * rs_ohm * 1e6;
}

/* ==================================================================
 * L3: Total Unadjusted Error (TUE) Computation
 *
 * TUE combines all DC errors at the IA output:
 *   V_out_error = G * V_err_RTI + V_ref_error
 *
 * This is the error before calibration correction.
 *
 * Complexity: O(1)
 * ================================================================== */
double total_unadjusted_error_mv(const ia_spec_t *spec, double G,
                                   double v_rti_error_uv) {
    (void)spec;
    return G * v_rti_error_uv * 1e-3;
}

/* ==================================================================
 * L5: Optimal Source Resistance for Minimum Noise Figure
 *
 * Noise figure is minimized when the amplifier's voltage and
 * current noise contributions are balanced:
 *
 *   R_s_opt = en_amp / in_amp
 *
 * Where en_amp is the amplifier's input voltage noise density
 * and in_amp is its input current noise density.
 *
 * For AD620: en = 9 nV/rtHz, in = 0.1 pA/rtHz (at 1 kHz)
 *   R_s_opt = 9e-9 / 0.1e-12 = 90 kOhm
 *
 * For source R < R_s_opt: voltage noise dominates
 * For source R > R_s_opt: current noise dominates
 *
 * Complexity: O(1)
 * ================================================================== */
double optimal_source_resistance(double en_nv_rhz, double in_pa_rhz) {
    if (in_pa_rhz <= 0.0) return INFINITY;
    return (en_nv_rhz * 1e-9) / (in_pa_rhz * 1e-12);
}

/* ==================================================================
 * L3: Resistor Noise from Gain Network
 *
 * In a 3-op-amp IA, the R_G resistor contributes noise:
 *   v_n_RG_RTI = sqrt(4kTR_G) * dG/dR_G * R_G / G
 *
 * At G=100, R1=25k, R_G=500:
 *   v_n_RG = sqrt(4*1.38e-23*300*500) = 2.87 nV/rtHz
 *   Referred to input: v_n_RG_RTI = 2.87 * (2*25k/500) / 100 = 2.87 nV/rtHz
 *
 * At G=1000, R_G=50:
 *   v_n_RG = sqrt(4*1.38e-23*300*50) = 0.91 nV/rtHz
 *
 * Lower R_G values give lower resistor noise but higher current
 * and power dissipation.
 *
 * Complexity: O(1)
 * ================================================================== */
double rg_noise_contribution_rti_nv_rhz(double rg_ohm, double r1_ohm,
                                          double G, double temp_c) {
    if (rg_ohm <= 0.0 || G <= 1.0) return 0.0;
    double k_B = 1.380649e-23;
    double T_K = temp_c + 273.15;

    double vn_rg = sqrt(4.0 * k_B * T_K * rg_ohm);
    /* dVout/dR_G * R_G/G = voltage gain sensitivity * R_G/G */
    double sensitivity = (2.0 * r1_ohm / rg_ohm) / G;
    return vn_rg * sensitivity * 1e9;  /* nV/rtHz */
}