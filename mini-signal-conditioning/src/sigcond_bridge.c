/**
 * sigcond_bridge.c - Wheatstone Bridge Implementation
 *
 * Bridge transfer functions, completion network design, balance,
 * lead wire compensation (3-wire, 4-wire Kelvin), sensitivity
 * optimization, and nonlinearity correction.
 *
 * L2: Wheatstone bridge principle, lead compensation
 * L3: bridge transfer functions, nonlinearity correction
 * L4: Kirchhoff laws, bridge balance condition
 */

#include "sigcond_bridge.h"
#include <stdlib.h>

/* Exact bridge output:
 * V_out = V_exc * (R3/(R2+R3) - R4/(R1+R4))
 *
 * This is the most general form. No assumptions about resistor
 * equality or small changes.
 */
double bridge_output_exact(double v_exc, double r1, double r2,
                            double r3, double r4)
{
    if (r2 + r3 <= 0.0 || r1 + r4 <= 0.0) return 0.0;
    double va = v_exc * r3 / (r2 + r3);  /* Voltage at node a */
    double vb = v_exc * r4 / (r1 + r4);  /* Voltage at node b */
    return va - vb;
}

/* Linear approximation:
 * For R1=R2=R3=R4=R and one active arm:
 *   V_out ~ V_exc * dR / (4R)   quarter bridge
 *   V_out ~ V_exc * dR / (2R)   half bridge
 *   V_out ~ V_exc * dR / R      full bridge
 */
double bridge_output_linear(double v_exc, double r_nominal,
                             double delta_r, unsigned active_arms)
{
    if (r_nominal <= 0.0 || active_arms == 0) return 0.0;

    double factor;
    switch (active_arms) {
        case 1: factor = 4.0; break;  /* Quarter bridge */
        case 2: factor = 2.0; break;  /* Half bridge */
        case 4: factor = 1.0; break;  /* Full bridge */
        default: factor = 4.0; break;
    }
    return v_exc * delta_r / (factor * r_nominal);
}

/* Nonlinearity error for quarter bridge with single active arm R+delta_R:
 * Exact: V_out = V_exc * delta_R / (4R + 2*delta_R)
 * Linear: V_out_lin = V_exc * delta_R / (4R)
 * NL = (V_linear/V_exact - 1) * 100%
 *    = delta_R / (2R) * 100%  (first-order)
 */
double bridge_nonlinearity_error(double r_nominal, double delta_r,
                                  unsigned active_arms)
{
    if (r_nominal <= 0.0) return 0.0;

    double strain = delta_r / r_nominal;
    /* For quarter bridge: NL ~ strain/2
     * For half bridge:   NL ~ strain/4 (adjacent arms cancel)
     * For full bridge:   NL ~ 0 (symmetric cancellation)
     */
    switch (active_arms) {
        case 1: return fabs(strain) / 2.0 * 100.0;
        case 2: return fabs(strain) / 4.0 * 100.0;
        case 4: return 0.0;  /* Ideally zero */
        default: return fabs(strain) / 2.0 * 100.0;
    }
}

/* Nonlinearity correction for full bridge:
 * V_out = V_exc * (delta_R/R)  (exact for full bridge with equal arms)
 * For GF gauge factor: delta_R/R = GF * strain
 * strain = V_out / (V_exc * GF)
 *
 * For quarter bridge:
 * V_out = V_exc * delta_R / (4R + 2*delta_R)
 * Solving for delta_R:
 *   delta_R = 4R * V_out / (V_exc - 2*V_out)
 *   strain = delta_R / (R * GF)
 */
double bridge_nonlinearity_correct(double v_measured, double v_exc,
                                    double gauge_factor)
{
    if (fabs(gauge_factor) < 1e-15 || fabs(v_exc) < 1e-15) return 0.0;

    /* Assume quarter bridge correction */
    double denom = v_exc - 2.0 * v_measured;
    if (fabs(denom) < 1e-15) return 0.0;

    double delta_r_over_r = 4.0 * v_measured / denom;
    return delta_r_over_r / gauge_factor;
}

void bridge_completion_design(bridge_topology_t topology,
                               double r_sensor_nominal,
                               bridge_completion_t *completion)
{
    if (!completion || r_sensor_nominal <= 0.0) return;

    completion->r1 = r_sensor_nominal;
    completion->r2 = r_sensor_nominal;
    completion->r3 = r_sensor_nominal;

    switch (topology) {
        case BRIDGE_QUARTER_SINGLE:
            completion->num_wires = 2;
            break;
        case BRIDGE_QUARTER_3WIRE:
            completion->num_wires = 3;
            break;
        case BRIDGE_QUARTER_4WIRE:
            completion->num_wires = 4;
            break;
        case BRIDGE_HALF_ADJACENT:
        case BRIDGE_HALF_OPPOSITE:
            completion->num_wires = 3;
            break;
        case BRIDGE_FULL:
            completion->num_wires = 4;
            break;
        default:
            completion->num_wires = 2;
            break;
    }

    completion->r_balance_min = 0.0;
    completion->r_balance_max = r_sensor_nominal * 0.01; /* 1% balance pot */
    completion->r_lead_per_wire = 0.1;  /* ~1m of 24 AWG copper */
    completion->tcr_ppm_per_c = 10.0;   /* Low-TCR completion resistors */
}

double bridge_completion_error(double r_sensor_nominal,
                                double r_completion_actual)
{
    if (r_sensor_nominal <= 0.0) return INFINITY;
    return fabs(r_completion_actual - r_sensor_nominal) / r_sensor_nominal * 100.0;
}

bool bridge_is_balanced(double r1, double r2, double r3, double r4,
                         double tolerance_ratio)
{
    double left = r1 * r3;
    double right = r2 * r4;
    double r_nom_sq = r1 * r3; /* Use geometric mean */
    if (r_nom_sq <= 0.0) return true;

    return fabs(left - right) / r_nom_sq < tolerance_ratio;
}

double bridge_balance_pot_value(double r_nominal, double max_imbalance_pct)
{
    return r_nominal * max_imbalance_pct / 100.0;
}

double bridge_software_zero(double v_measured, double v_zero_offset)
{
    return v_measured - v_zero_offset;
}

/* Lead wire resistance: R = rho * L / A
 * Copper resistivity at 20C: rho = 1.68e-8 Ohm*m
 * AWG to area (mm^2): A = pi/4 * (0.127 * 92^((36-AWG)/39))^2
 * Temperature correction: rho(T) = rho_20 * (1 + alpha*(T-20))
 *   alpha_copper = 0.00393 /C (IACS)
 */
double bridge_lead_resistance(double length_m, double awg,
                               double temperature_c)
{
    /* Wire diameter in mm: d = 0.127 * 92^((36-AWG)/39) */
    double d_mm = 0.127 * pow(92.0, (36.0 - awg) / 39.0);
    double area_m2 = M_PI / 4.0 * d_mm * d_mm * 1e-6;

    if (area_m2 <= 0.0) return 0.0;

    double rho_20 = 1.68e-8;  /* Ohm*m */
    double alpha = 0.00393;    /* /C */
    double rho = rho_20 * (1.0 + alpha * (temperature_c - 20.0));

    return rho * length_m / area_m2;
}

double bridge_3wire_compensate(double v_exc, double v_sense_high,
                                double v_sense_low, double r_ref,
                                double r_lead_estimated)
{
    /* In 3-wire configuration:
     * V_high = I * (R_sensor + R_lead1)
     * V_low = I * R_lead3
     * If R_lead1 = R_lead3 (matched leads):
     *   R_sensor = (V_high - V_low) / I
     *            = (V_high - V_low) * R_ref / V_exc
     *
     * If leads not perfectly matched, residual error = (R_lead1-R_lead3).
     */
    if (fabs(v_exc) < 1e-15) return 0.0;

    double i_exc = v_exc / r_ref;
    double r_sensor = (v_sense_high - v_sense_low) / i_exc;
    return r_sensor;
}

double bridge_4wire_kelvin(double i_excitation_a,
                            double v_sensed_v)
{
    if (fabs(i_excitation_a) < 1e-15) return INFINITY;
    return v_sensed_v / i_excitation_a;
}

double bridge_lead_temp_error(double r_lead_20c, double temp_c,
                               double r_sensor_nominal)
{
    if (r_sensor_nominal <= 0.0) return 0.0;
    double alpha = 0.00393;  /* Copper TCR */
    double r_lead_at_T = r_lead_20c * (1.0 + alpha * (temp_c - 20.0));
    double r_lead_change = r_lead_at_T - r_lead_20c;
    return r_lead_change / r_sensor_nominal * 100.0;
}

double bridge_sensitivity_strain(double v_exc, double gauge_factor,
                                  unsigned active_arms)
{
    if (active_arms == 0) return 0.0;
    double factor;
    switch (active_arms) {
        case 1: factor = 4.0; break;
        case 2: factor = 2.0; break;
        case 4: factor = 1.0; break;
        default: factor = 4.0; break;
    }
    return v_exc * gauge_factor / factor;
}

double bridge_sensitivity_loadcell(double v_exc,
                                    double rated_output_mv_per_v,
                                    double full_scale_load)
{
    if (full_scale_load <= 0.0) return 0.0;
    return v_exc * rated_output_mv_per_v * 1e-3 / full_scale_load;
}

double bridge_min_detectable_strain(double v_exc, double gauge_factor,
                                     double noise_rti_uv, double gain)
{
    if (fabs(v_exc * gauge_factor * gain) < 1e-15) return INFINITY;
    return noise_rti_uv * 4.0 / (v_exc * gauge_factor * gain);
}

double bridge_output_strain(double v_exc, double gauge_factor,
                             double strain, unsigned active_arms)
{
    double factor;
    switch (active_arms) {
        case 1: factor = 4.0; break;
        case 2: factor = 2.0; break;
        case 4: factor = 1.0; break;
        default: factor = 4.0; break;
    }
    return v_exc * gauge_factor * strain / factor;
}

double bridge_strain_from_output(double v_out, double v_exc,
                                  double gauge_factor,
                                  unsigned active_arms)
{
    if (fabs(v_exc * gauge_factor) < 1e-15) return 0.0;

    double factor;
    switch (active_arms) {
        case 1: factor = 4.0; break;
        case 2: factor = 2.0; break;
        case 4: factor = 1.0; break;
        default: factor = 4.0; break;
    }
    return v_out * factor / (v_exc * gauge_factor);
}
