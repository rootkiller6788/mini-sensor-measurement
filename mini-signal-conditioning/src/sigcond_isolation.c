/**
 * sigcond_isolation.c - Galvanic Isolation Implementation
 *
 * Implements isolation safety margins, leakage current calculation,
 * CMRR analysis for capacitive barriers, CMTI requirements,
 * creepage/clearance per IEC standards, and power transfer estimation.
 *
 * L2: galvanic isolation, common-mode rejection
 * L3: barrier capacitance modeling
 * L4: IEC 60601-1/60664-1 safety standards
 */

#include "sigcond_isolation.h"
#include <stdlib.h>

/* IEC 60601-1: V_withstand = 2 * V_work + 1000V (for V_work <= 1000V)
 * Above 1000V: V_withstand = V_work + 2000V */
double isolation_withstand_from_working(double working_vrms)
{
    if (working_vrms <= 0.0) return 0.0;
    if (working_vrms <= 1000.0)
        return 2.0 * working_vrms + 1000.0;
    else
        return working_vrms + 2000.0;
}

/* Altitude derating: V_allowable = V_rated * (1 - h/30000) for h < 30000ft
 * In meters: V_allowable = V_rated * (1 - h/9144) for h < 9144m
 * Beyond 9144m: atmospheric pressure too low, consult Paschen's law minimum. */
double isolation_altitude_derate(double rated_voltage, double altitude_meters)
{
    if (altitude_meters <= 0.0) return rated_voltage;
    if (altitude_meters >= 9144.0) return rated_voltage * 0.3; /* Approximate */
    return rated_voltage * (1.0 - altitude_meters / 9144.0);
}

/* IEC 60664-1 creepage minimum:
 * Based on working voltage, pollution degree, and material group.
 * Simplified formula (approximate, for design guidance):
 *   Creepage ~ V_work / K where K depends on pollution degree.
 *   PD1 (clean): K ~ 1.0 kV/mm
 *   PD2 (typical): K ~ 0.5 kV/mm
 *   PD3 (conductive dust): K ~ 0.2 kV/mm
 *
 * Material Group I: CTI >= 600V (glass, ceramic)
 * Material Group II: CTI 400-600V
 * Material Group IIIa/IIIb: CTI < 400V (FR-4 PCB)
 */
double isolation_creepage(double working_vrms, int pollution_degree,
                           int material_group)
{
    if (working_vrms <= 0.0) return 0.0;

    /* Base creepage from IEC 60664-1 Table F.1 */
    double base_mm;

    if (working_vrms <= 50.0)       base_mm = 0.5;
    else if (working_vrms <= 100.0) base_mm = 0.8;
    else if (working_vrms <= 150.0) base_mm = 1.5;
    else if (working_vrms <= 300.0) base_mm = 3.0;
    else if (working_vrms <= 600.0) base_mm = 5.5;
    else if (working_vrms <= 1000.0) base_mm = 8.0;
    else                            base_mm = working_vrms * 0.008;

    /* Pollution degree multiplier */
    double pd_mult;
    switch (pollution_degree) {
        case 1: pd_mult = 0.8;  break;  /* Clean, controlled environment */
        case 2: pd_mult = 1.0;  break;  /* Normal, non-conductive pollution */
        case 3: pd_mult = 1.6;  break;  /* Conductive pollution possible */
        case 4: pd_mult = 2.5;  break;  /* Persistent conductivity */
        default: pd_mult = 1.0; break;
    }

    /* Material group multiplier (I-II-III) */
    double mg_mult;
    switch (material_group) {
        case 1: mg_mult = 1.0;  break;  /* CTI >= 600 */
        case 2: mg_mult = 1.4;  break;  /* CTI 400-600 */
        case 3: mg_mult = 1.8;  break;  /* CTI 175-400 */
        default: mg_mult = 1.4; break;
    }

    return base_mm * pd_mult * mg_mult;
}

double isolation_clearance(double working_vrms, int overvoltage_cat,
                            double altitude_m)
{
    if (working_vrms <= 0.0) return 0.0;

    /* Base impulse withstand from overvoltage category */
    double impulse_kv;
    switch (overvoltage_cat) {
        case 1: impulse_kv = 0.5;  break;  /* Protected electronic */
        case 2: impulse_kv = 2.5;  break;  /* Appliance */
        case 3: impulse_kv = 4.0;  break;  /* Distribution */
        case 4: impulse_kv = 6.0;  break;  /* Origin of installation */
        default: impulse_kv = 2.5; break;
    }

    /* Clearance ~ 0.01 mm/V for uniform field, 0.02 mm/V for non-uniform
     * Simplified: d(mm) = impulse_kV * 0.8 (per IEC 60664-1 Table F.2) */
    double base_clearance = impulse_kv * 0.8; /* mm */

    /* Altitude correction: multiply by 1.4 at 3000m, 2.2 at 5000m */
    double alt_factor = 1.0 + altitude_m * 0.00014;

    return base_clearance * alt_factor;
}

double isolation_leakage_current(double barrier_cap_pf,
                                  double cm_voltage_v,
                                  double cm_frequency_hz)
{
    /* I_leak = C * dV/dt = C * V_cm * 2*pi*f */
    return barrier_cap_pf * 1e-12 * cm_voltage_v * 2.0 * M_PI * cm_frequency_hz;
}

double isolation_cmti_required(double signal_bandwidth_hz,
                                double max_cm_voltage_v,
                                double max_bit_error)
{
    /* Required CMTI to keep jitter < 1/(2*pi*BW*SNR)
     * Simplified: CMTI >= dV_cm / (1/signal_BW)
     *            = dV_cm * signal_BW
     */
    if (signal_bandwidth_hz <= 0.0) return 0.0;
    return max_cm_voltage_v * signal_bandwidth_hz / 1e9; /* kV/us */
}

double isolation_mtbf_estimate(double working_v, double rated_v,
                                double base_mtbf_hours)
{
    /* Voltage stress derating:
     * MTBF = MTBF_base * exp(-alpha * V_work/V_rated)
     * Simplified: linear derating for low stress ratios.
     */
    double stress_ratio = working_v / rated_v;
    if (stress_ratio >= 1.0) return base_mtbf_hours * 0.01;
    if (stress_ratio <= 0.0) return base_mtbf_hours * 10.0;

    /* Exponential derating model (typical for dielectric breakdown) */
    return base_mtbf_hours * exp(-2.0 * stress_ratio);
}

double isolation_total_cmrr(double cmrr_pre_barrier_db,
                             double cmrr_barrier_db,
                             double cmrr_post_barrier_db)
{
    /* 1/CMRR_total = 1/CMRR_pre + 1/CMRR_barrier + 1/CMRR_post
     * Convert to linear, sum reciprocals, convert back.
     */
    double cmrr_pre_lin = pow(10.0, cmrr_pre_barrier_db / 20.0);
    double cmrr_bar_lin = pow(10.0, cmrr_barrier_db / 20.0);
    double cmrr_post_lin = pow(10.0, cmrr_post_barrier_db / 20.0);

    if (cmrr_pre_lin < 1.0) cmrr_pre_lin = 1.0;
    if (cmrr_bar_lin < 1.0) cmrr_bar_lin = 1.0;
    if (cmrr_post_lin < 1.0) cmrr_post_lin = 1.0;

    double inv_total = 1.0/cmrr_pre_lin + 1.0/cmrr_bar_lin + 1.0/cmrr_post_lin;
    double cmrr_total_lin = 1.0 / inv_total;

    return 20.0 * log10(cmrr_total_lin);
}

double isolation_cap_barrier_cmrr(double barrier_cap_pf,
                                   double diff_impedance_ohm,
                                   double frequency_hz)
{
    if (barrier_cap_pf <= 0.0 || diff_impedance_ohm <= 0.0 || frequency_hz <= 0.0)
        return 0.0;

    double z_barrier = 1.0 / (2.0 * M_PI * frequency_hz * barrier_cap_pf * 1e-12);
    return 20.0 * log10(z_barrier / diff_impedance_ohm);
}

double isolation_pcb_cmrr(double c_imbalance_pf,
                           double common_mode_z_ohm,
                           double frequency_hz)
{
    if (c_imbalance_pf <= 0.0 || common_mode_z_ohm <= 0.0 || frequency_hz <= 0.0)
        return 120.0; /* Ideal */

    double imbalance_z = 1.0 / (2.0 * M_PI * frequency_hz * c_imbalance_pf * 1e-12);
    return 20.0 * log10(imbalance_z / common_mode_z_ohm);
}

double isolation_power_transfer(double v_primary, double f_switching_hz,
                                 double b_max_t, double a_core_m2,
                                 unsigned n_primary)
{
    /* Transformer power transfer:
     * P = (V_primary * duty * A_core * B_max * f_sw) / N_primary
     * Duty cycle assumed 0.4 for push-pull.
     */
    if (n_primary == 0 || f_switching_hz <= 0.0) return 0.0;
    double duty = 0.4;
    return v_primary * duty * a_core_m2 * b_max_t * f_switching_hz / n_primary;
}

double isolation_power_budget(unsigned num_channels,
                               double power_per_channel_mw,
                               double adc_power_mw,
                               double mcu_power_mw)
{
    double total_mw = num_channels * power_per_channel_mw + adc_power_mw + mcu_power_mw;
    return total_mw * 1.3; /* 30% margin */
}

double isolation_gain_error(double barrier_insertion_loss_db)
{
    /* Gain_error = 1 - 10^(-IL_dB/20) */
    return 1.0 - pow(10.0, -barrier_insertion_loss_db / 20.0);
}
