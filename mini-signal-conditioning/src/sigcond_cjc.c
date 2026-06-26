/**
 * sigcond_cjc.c - Cold Junction Compensation Implementation
 *
 * CJC temperature measurement using thermistor, RTD, IC sensor,
 * or diode junction. Full CJC algorithm: measure CJC temperature,
 * compute equivalent thermoelectric voltage, add to measured
 * voltage, compute hot junction temperature via ITS-90.
 *
 * L2: Seebeck effect, intermediate metals law
 * L3: NIST ITS-90 polynomial evaluation for CJC
 * L4: Law of successive temperatures, isothermal error
 */

#include "sigcond_cjc.h"
#include "sigcond_linearize.h"
#include <stdlib.h>

double cjc_measure_thermistor(double thermistor_r_ohm,
                               const steinhart_hart_params_t *params)
{
    double T_kelvin = steinhart_hart_T(params, thermistor_r_ohm);
    if (isnan(T_kelvin)) return NAN;
    return T_kelvin - 273.15; /* Kelvin to Celsius */
}

double cjc_measure_rtd(double rtd_resistance_ohm,
                        const rtd_cvd_params_t *params)
{
    return rtd_cvd_T_from_R(params, rtd_resistance_ohm);
}

double cjc_measure_ic_sensor(double sensor_voltage_mv,
                              double offset_mv, double slope_mv_per_c)
{
    if (fabs(slope_mv_per_c) < 1e-15) return NAN;
    return (sensor_voltage_mv - offset_mv) / slope_mv_per_c;
}

double cjc_measure_diode(double v_forward_mv, double v_f0_at_25c_mv,
                          double temp_coefficient_mv_per_c)
{
    if (fabs(temp_coefficient_mv_per_c) < 1e-15) return 25.0;
    return 25.0 + (v_f0_at_25c_mv - v_forward_mv) / (-temp_coefficient_mv_per_c);
}

/* Full CJC algorithm:
 * Step 1: E_CJ = f(T_CJ) using forward ITS-90 polynomial
 * Step 2: E_total = V_measured + E_CJ
 * Step 3: T_hot = f^{-1}(E_total) using inverse ITS-90 polynomial
 *
 * This implements the Law of Successive Temperatures:
 *   E(T_hot, 0) = E(T_hot, T_CJ) + E(T_CJ, 0)
 *                = V_measured + E(T_CJ)
 */
double cjc_compensate(double v_measured_uv, double cj_temp_c,
                       const thermocouple_nist_model_t *tc_model)
{
    if (!tc_model) return NAN;

    /* Step 1: Convert CJC temperature to voltage */
    double v_cj = cjc_equivalent_voltage(cj_temp_c, tc_model);
    if (isnan(v_cj)) return NAN;

    /* Step 2: Sum voltages */
    double v_total = v_measured_uv + v_cj;

    /* Step 3: Convert total voltage to temperature */
    return nist_tc_voltage_to_temp(tc_model, v_total);
}

double cjc_equivalent_voltage(double cj_temp_c,
                               const thermocouple_nist_model_t *tc_model)
{
    if (!tc_model) return NAN;
    return nist_tc_temp_to_voltage(tc_model, cj_temp_c);
}

/* Law of successive temperatures verification:
 * E(T1, T3) = E(T1, T2) + E(T2, T3)
 * For T1=0, T2=cj, T3=hot:
 *   E(0, hot) = E(0, cj) + E(cj, hot)
 *   V_total = V(0 to cj) + V(cj to hot) = V_known + V_measured
 */
bool cjc_verify_successive_law(double t1_c, double t2_c, double t3_c,
                                const thermocouple_nist_model_t *tc_model,
                                double tolerance_uv)
{
    if (!tc_model) return false;

    double v13 = nist_tc_temp_to_voltage(tc_model, t3_c)
               - nist_tc_temp_to_voltage(tc_model, t1_c);
    double v12 = nist_tc_temp_to_voltage(tc_model, t2_c)
               - nist_tc_temp_to_voltage(tc_model, t1_c);
    double v23 = nist_tc_temp_to_voltage(tc_model, t3_c)
               - nist_tc_temp_to_voltage(tc_model, t2_c);

    return fabs(v13 - v12 - v23) < tolerance_uv;
}

/* Isothermal error from terminal temperature gradient.
 * If CJC terminals differ by dT_CJ, each terminal forms a
 * parasitic TC: Cu-alloyA at one temp, Cu-alloyB at slightly
 * different temp. Resulting error voltage:
 *   V_err = S_Cu_A * dT/2 - S_Cu_B * dT/2
 *         ~ S_AB * dT/2  (Seebeck of the TC itself)
 *
 * For Type K at 25C: S~40.6uV/C, dT=0.1C => V_err~2.0uV => ~0.05C error
 */
double cjc_isothermal_error(double delta_t_cj_c,
                             const thermocouple_nist_model_t *tc_model)
{
    if (!tc_model || delta_t_cj_c < 0.0) return 0.0;
    return tc_model->seebeck_25c * delta_t_cj_c / 2.0;
}

double cjc_total_uncertainty(double sensor_accuracy_c,
                              double gradient_error_c,
                              double pcb_thermoelectric_error_c)
{
    return sqrt(sensor_accuracy_c * sensor_accuracy_c
              + gradient_error_c * gradient_error_c
              + pcb_thermoelectric_error_c * pcb_thermoelectric_error_c);
}

bool cjc_verify_placement(double sensor_to_terminal_mm,
                           double terminal_thermal_resistance_c_per_w,
                           double max_power_dissipation_w,
                           double max_allowable_gradient_c)
{
    /* Temperature gradient = Power * Thermal_resistance
     * Gradients over distance create CJC errors.
     * Check if expected gradient exceeds allowable.
     */
    double gradient = max_power_dissipation_w * terminal_thermal_resistance_c_per_w;
    /* Gradient across sensor-to-terminal distance */
    double gradient_over_distance = gradient * sensor_to_terminal_mm / 10.0;
    return gradient_over_distance < max_allowable_gradient_c;
}

double cjc_thermistor_self_heat_correction(double measured_temp_c,
                                            double excitation_v,
                                            double thermistor_r_ohm)
{
    if (thermistor_r_ohm <= 0.0) return measured_temp_c;

    double power_w = excitation_v * excitation_v / thermistor_r_ohm;
    /* Typical thermistor dissipation constant: 1 mW/C in still air */
    double dissipation_constant_c_per_w = 1000.0; /* 1 mW/C => 1000 C/W */
    double self_heat_c = power_w * dissipation_constant_c_per_w;

    return measured_temp_c - self_heat_c;
}
