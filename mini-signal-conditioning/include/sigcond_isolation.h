/**
 * sigcond_isolation.h - Galvanic Isolation Design
 *
 * Optical, transformer, capacitive, RF-coupled isolation barriers.
 * Safety margins, CMTI, creepage/clearance, CMRR analysis.
 *
 * L2-L4: galvanic isolation, barrier capacitance, safety standards
 * Reference: Kester (1999), IEC 60601-1, VDE 0884-10/11
 */
#ifndef SIGCOND_ISOLATION_H
#define SIGCOND_ISOLATION_H
#include "sigcond_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

/** Withstand voltage from working voltage per IEC 60601-1: V_withstand = 2*V_work + 1000 */
double isolation_withstand_from_working(double working_vrms);

/** Altitude derating: V_allowable(h) = V_rated * (1 - h/30000) for h < 30k ft */
double isolation_altitude_derate(double rated_voltage, double altitude_meters);

/** Minimum creepage per IEC 60664-1 based on working voltage and pollution degree */
double isolation_creepage(double working_vrms, int pollution_degree, int material_group);

/** Minimum clearance per IEC 60664-1 */
double isolation_clearance(double working_vrms, int overvoltage_cat, double altitude_m);

/** Leakage current through barrier: I_leak = C_barrier * dV_cm/dt */
double isolation_leakage_current(double barrier_cap_pf, double cm_voltage_v, double cm_frequency_hz);

/** Required CMTI for target error rate: CMTI >= dV_cm_max / (1/(2*BW_signal)) */
double isolation_cmti_required(double signal_bandwidth_hz, double max_cm_voltage_v, double max_bit_error);

/** Barrier MTBF estimate with voltage stress derating */
double isolation_mtbf_estimate(double working_v, double rated_v, double base_mtbf_hours);

/** Total CMRR of isolated chain: 1/CMRR_total = sum(1/CMRR_i) */
double isolation_total_cmrr(double cmrr_pre_barrier_db, double cmrr_barrier_db, double cmrr_post_barrier_db);

/** Capacitive barrier CMRR at frequency: CMRR(f)=20*log10(1/(2*pi*f*C_barrier*R_diff)) */
double isolation_cap_barrier_cmrr(double barrier_cap_pf, double diff_impedance_ohm, double frequency_hz);

/** PCB CMRR degradation from capacitive imbalance */
double isolation_pcb_cmrr(double c_imbalance_pf, double common_mode_z_ohm, double frequency_hz);

/** Power transfer across isolation transformer */
double isolation_power_transfer(double v_primary, double f_switching_hz, double b_max_t, double a_core_m2, unsigned n_primary);

/** Isolated side power budget estimate */
double isolation_power_budget(unsigned num_channels, double power_per_channel_mw, double adc_power_mw, double mcu_power_mw);

/** Isolation amplifier gain error from barrier insertion loss */
double isolation_gain_error(double barrier_insertion_loss_db);

#ifdef __cplusplus
}
#endif
#endif /* SIGCOND_ISOLATION_H */
