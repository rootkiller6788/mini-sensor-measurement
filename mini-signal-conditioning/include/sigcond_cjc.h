/**
 * sigcond_cjc.h - Cold Junction Compensation for Thermocouples
 *
 * CJC using thermistor, RTD, IC sensor, or diode. NIST ITS-90
 * polynomial evaluation for compensated temperature.
 *
 * L2-L4: Seebeck effect, intermediate metals law, CJC principle
 * Reference: NIST Monograph 175 (1993), ASTM E230/E230M
 */
#ifndef SIGCOND_CJC_H
#define SIGCOND_CJC_H
#include "sigcond_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

/** Thermistor-based CJC: measure T_cj from NTC thermistor resistance using Steinhart-Hart */
double cjc_measure_thermistor(double thermistor_r_ohm, const steinhart_hart_params_t *params);

/** RTD-based CJC: measure T_cj from Pt100/Pt1000 resistance using CVD equation */
double cjc_measure_rtd(double rtd_resistance_ohm, const rtd_cvd_params_t *params);

/** IC sensor CJC: V_out(mV) = offset + slope * T(C). LM35: 10mV/C, 0C=0mV. TMP36: 10mV/C, 25C=750mV. */
double cjc_measure_ic_sensor(double sensor_voltage_mv, double offset_mv, double slope_mv_per_c);

/** Diode junction CJC: V_f(T) = V_f0 - (T-25C)*2mV/C */
double cjc_measure_diode(double v_forward_mv, double v_f0_at_25c_mv, double temp_coefficient_mv_per_c);

/** Full CJC algorithm: 1) E_cj = fwd_poly(T_cj), 2) E_total = V_meas + E_cj, 3) T_hot = inv_poly(E_total) */
double cjc_compensate(double v_measured_uv, double cj_temp_c, const thermocouple_nist_model_t *tc_model);

/** Equivalent thermoelectric voltage for a given CJC temperature */
double cjc_equivalent_voltage(double cj_temp_c, const thermocouple_nist_model_t *tc_model);

/** Verify law of successive temperatures: |E(T1,T3) - E(T1,T2) - E(T2,T3)| < tol */
bool cjc_verify_successive_law(double t1_c, double t2_c, double t3_c, const thermocouple_nist_model_t *tc_model, double tolerance_uv);

/** Isothermal error from terminal block temperature gradient dT_CJ */
double cjc_isothermal_error(double delta_t_cj_c, const thermocouple_nist_model_t *tc_model);

/** Total CJC uncertainty (RSS of sensor, gradient, and PCB thermoelectric errors) */
double cjc_total_uncertainty(double sensor_accuracy_c, double gradient_error_c, double pcb_thermoelectric_error_c);

/** CJC sensor placement verification: thermal coupling adequacy */
bool cjc_verify_placement(double sensor_to_terminal_mm, double terminal_thermal_resistance_c_per_w, double max_power_dissipation_w, double max_allowable_gradient_c);

/** Self-heating correction for CJC thermistor sensor */
double cjc_thermistor_self_heat_correction(double measured_temp_c, double excitation_v, double thermistor_r_ohm);

#ifdef __cplusplus
}
#endif
#endif /* SIGCOND_CJC_H */
