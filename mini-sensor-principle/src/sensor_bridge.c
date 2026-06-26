/**
 * @file sensor_bridge.c
 * @brief Wheatstone bridge circuit analysis and signal conditioning
 *
 * Level: L2 — Core Concepts · L3 — Mathematical Structures
 */

#include "sensor_bridge.h"
#include <math.h>

/* ═══════════════════════════════════════════════════════════════════════
 * L2: Bridge Initialization
 * ═══════════════════════════════════════════════════════════════════════ */

int bridge_init_balanced(wheatstone_bridge_t *b, double R0, double V_ex)
{
    if (!b || R0 <= 0.0) return -1;
    b->R1 = b->R2 = b->R3 = b->R4 = R0;
    b->V_excitation = V_ex;
    b->type = BRIDGE_QUARTER;
    return 0;
}

int bridge_set_quarter_active(wheatstone_bridge_t *b,
                               double R_nominal, double R_sensed,
                               int arm_index)
{
    if (!b || R_nominal <= 0.0 || R_sensed <= 0.0 || arm_index < 1 || arm_index > 4)
        return -1;

    /* Set all to nominal, then override the active arm */
    b->R1 = b->R2 = b->R3 = b->R4 = R_nominal;
    switch (arm_index) {
    case 1: b->R1 = R_sensed; break;
    case 2: b->R2 = R_sensed; break;
    case 3: b->R3 = R_sensed; break;
    case 4: b->R4 = R_sensed; break;
    default: return -1;
    }
    b->type = BRIDGE_QUARTER;
    return 0;
}

int bridge_set_half_adjacent(wheatstone_bridge_t *b,
                              double R_nominal,
                              double R_upper_active, double R_lower_active)
{
    if (!b || R_nominal <= 0.0) return -1;
    b->R1 = R_upper_active;
    b->R2 = R_nominal;
    b->R3 = R_nominal;
    b->R4 = R_lower_active;
    b->type = BRIDGE_HALF_ADJACENT;
    return 0;
}

int bridge_set_full_pushpull(wheatstone_bridge_t *b,
                              double R_nominal, double delta_R)
{
    if (!b || R_nominal <= 0.0) return -1;
    /* Push-pull: opposite arms change in opposite directions
     * R1 and R3: R_nominal + ΔR (increase under tension)
     * R2 and R4: R_nominal - ΔR (decrease under tension for Poisson effect)
     */
    b->R1 = R_nominal + delta_R;
    b->R2 = R_nominal - delta_R;
    b->R3 = R_nominal + delta_R;
    b->R4 = R_nominal - delta_R;
    b->type = BRIDGE_FULL;
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L3: Bridge Output Voltage & Analysis
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Exact bridge output voltage.
 *
 * V_out = V_ex · [R₃/(R₂+R₃) - R₄/(R₁+R₄)]
 *
 * This is the fundamental Wheatstone bridge equation (Wheatstone, 1843).
 * It is derived from voltage division on each leg.
 *
 * Balanced condition: R₁·R₃ = R₂·R₄ ↔ V_out = 0
 */
double bridge_output_voltage_exact(const wheatstone_bridge_t *b)
{
    if (!b) return NAN;

    double v_left  = b->R4 / (b->R1 + b->R4);  /* voltage divider, left leg */
    double v_right = b->R3 / (b->R2 + b->R3);   /* voltage divider, right leg */

    return b->V_excitation * (v_right - v_left);
}

/**
 * @brief Approximate bridge output for small ΔR.
 *
 * This linearized form is used for quick estimates:
 *
 * Quarter bridge: V_out ≈ V_ex · ΔR / (4·R₀)
 * Half bridge:    V_out ≈ V_ex · ΔR / (2·R₀)
 * Full bridge:    V_out ≈ V_ex · ΔR / R₀
 *
 * The approximation error for quarter bridge:
 *   |ε| ≈ (ΔR/R₀)² / (2 + ΔR/R₀)  →  < 0.5% for |ΔR/R₀| < 0.01
 */
double bridge_output_voltage_approx(const wheatstone_bridge_t *b, double R0,
                                     double delta_R)
{
    if (!b || R0 <= 0.0) return NAN;

    switch (b->type) {
    case BRIDGE_QUARTER:
        return b->V_excitation * delta_R / (4.0 * R0);
    case BRIDGE_HALF_ADJACENT:
        return b->V_excitation * delta_R / (2.0 * R0);
    case BRIDGE_HALF_OPPOSITE:
        return 0.0;  /* symmetric ΔR, no output */
    case BRIDGE_FULL:
        return b->V_excitation * delta_R / R0;
    default:
        return 0.0;
    }
}

double bridge_nonlinearity_error(const wheatstone_bridge_t *b, double R0,
                                  double delta_R)
{
    if (!b || R0 <= 0.0) return NAN;

    /* Compute the exact output with ΔR */
    wheatstone_bridge_t temp = *b;
    if (b->type == BRIDGE_QUARTER) {
        temp.R1 = R0 + delta_R;
    }

    double V_exact = bridge_output_voltage_exact(&temp);
    double V_approx = bridge_output_voltage_approx(b, R0, delta_R);

    if (fabs(V_exact) < 1e-30) return 0.0;
    return fabs(V_exact - V_approx) / fabs(V_exact) * 100.0;
}

/**
 * @brief Numerically compute bridge sensitivity dV_out/dR for a specific arm.
 *
 * Uses central finite difference with 0.01 Ω perturbation.
 */
double bridge_sensitivity(const wheatstone_bridge_t *b, int active_arm)
{
    if (!b || active_arm < 1 || active_arm > 4) return NAN;

    wheatstone_bridge_t b_hi = *b;
    wheatstone_bridge_t b_lo = *b;
    double R_nom;
    double dR_factor = 0.001;  /* 0.1% perturbation */

    switch (active_arm) {
    case 1: R_nom = b->R1; b_hi.R1 *= (1.0 + dR_factor); b_lo.R1 *= (1.0 - dR_factor); break;
    case 2: R_nom = b->R2; b_hi.R2 *= (1.0 + dR_factor); b_lo.R2 *= (1.0 - dR_factor); break;
    case 3: R_nom = b->R3; b_hi.R3 *= (1.0 + dR_factor); b_lo.R3 *= (1.0 - dR_factor); break;
    case 4: R_nom = b->R4; b_hi.R4 *= (1.0 + dR_factor); b_lo.R4 *= (1.0 - dR_factor); break;
    default: return NAN;
    }

    double V_hi = bridge_output_voltage_exact(&b_hi);
    double V_lo = bridge_output_voltage_exact(&b_lo);
    double dR = R_nom - R_nom * (1.0 - dR_factor);  /* effective ΔR */

    if (dR < 1e-30) return 0.0;
    return (V_hi - V_lo) / (2.0 * dR);
}

/**
 * @brief Thevenin equivalent resistance seen looking into the bridge output terminals.
 *
 * R_th = (R₁ ∥ R₄) + (R₂ ∥ R₃)
 *
 * For a balanced bridge (all R₀): R_th = R₀/2 + R₀/2 = R₀
 *
 * This is the source resistance that the amplifier sees, which determines
 * amplifier loading and noise contribution.
 * For low noise, R_th should be ≤ ~1 kΩ (Johnson noise ≤ 4 nV/√Hz).
 */
double bridge_thevenin_resistance(const wheatstone_bridge_t *b)
{
    if (!b) return NAN;

    double R14 = (b->R1 * b->R4) / (b->R1 + b->R4);
    double R23 = (b->R2 * b->R3) / (b->R2 + b->R3);

    return R14 + R23;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L2: Bridge Non-Ideal Effects
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Lead wire error in 2-wire quarter bridge.
 *
 * In a 2-wire measurement, the lead resistance adds directly to the
 * gauge resistance: R_measured = R_gauge + 2·R_lead.
 *
 * Fractional error: ε = 2·R_lead / (R_gauge + 2·R_lead) ≈ 2·R_lead / R_0
 *
 * For R₀ = 120 Ω, R_lead = 0.5 Ω (AWG 24, 6 m): ε ≈ 2·0.5/120 = 0.83%
 *
 * The 3-wire connection eliminates this by measuring lead resistance
 * and subtracting it. The 4-wire (Kelvin) connection excites through
 * one pair and senses through another, virtually eliminating lead error.
 */
double bridge_lead_wire_error_2wire(double R_gauge, double R_lead_per_wire)
{
    if (R_gauge <= 0.0) return 1.0;
    return 2.0 * R_lead_per_wire / (R_gauge + 2.0 * R_lead_per_wire);
}

/**
 * @brief Thermal EMF offset estimate.
 *
 * When dissimilar metals form a junction at a temperature gradient,
 * a parasitic thermocouple is created (Seebeck effect).
 *
 * Common sources in bridge circuits:
 *   - Copper wire to solder (Cu-SnPb): ~1-3 µV/°C
 *   - Relay contacts: ~20-50 µV/°C
 *   - Connector pins (various platings): ~5-40 µV/°C
 *
 * For precision measurements (µV-level), AC excitation or
 * current-reversal techniques are used to cancel thermal EMF.
 */
double bridge_thermal_emf_offset(double seebeck_coeff_uV_per_K,
                                  double delta_T_K)
{
    return seebeck_coeff_uV_per_K * delta_T_K;
}

/**
 * @brief Power dissipated in a single bridge arm.
 *
 * P_arm = I²·R_arm
 *
 * The total bridge current: I_total = V_ex / R_eq
 * where R_eq = (R₁+R₄) ∥ (R₂+R₃)
 *
 * Current through left leg: I_left = V_ex / (R₁+R₄)
 * Power in R₁: P₁ = [V_ex/(R₁+R₄)]² · R₁
 */
double bridge_arm_power_dissipation(const wheatstone_bridge_t *b, int arm_index)
{
    if (!b) return NAN;

    double R_sum, R_arm;
    switch (arm_index) {
    case 1: R_sum = b->R1 + b->R4; R_arm = b->R1; break;
    case 2: R_sum = b->R2 + b->R3; R_arm = b->R2; break;
    case 3: R_sum = b->R2 + b->R3; R_arm = b->R3; break;
    case 4: R_sum = b->R1 + b->R4; R_arm = b->R4; break;
    default: return NAN;
    }

    if (R_sum <= 0.0) return 0.0;
    double I_leg = b->V_excitation / R_sum;
    return I_leg * I_leg * R_arm;
}

double bridge_self_heating_temp_rise(double power_per_arm,
                                      double thermal_resistance)
{
    return power_per_arm * thermal_resistance;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L3: Bridge Signal Conditioning
 * ═══════════════════════════════════════════════════════════════════════ */

double bridge_amplifier_gain_required(double V_bridge_fs, double V_adc_fs)
{
    if (V_bridge_fs <= 0.0 || V_adc_fs <= 0.0) return 1.0;
    return V_adc_fs / V_bridge_fs;
}

/**
 * @brief Required CMRR for bridge amplifier.
 *
 * The bridge output contains a large common-mode voltage (≈V_ex/2)
 * and a small differential signal. The amplifier must reject the
 * common-mode while amplifying the differential component.
 *
 * CMRR = |A_diff| / |A_cm|
 *
 * CMRR_min = V_cm / (V_diff_min · ε_allowed)
 *
 * In dB: CMRR_dB = 20·log₁₀(CMRR)
 *
 * Example: V_ex=5V → V_cm=2.5V, V_diff=10mV, ε=0.1%
 *          CMRR = 2.5/(0.01·0.001) = 250,000 → 108 dB
 *          → Requires a high-quality instrumentation amplifier
 */
double bridge_amplifier_cmrr_required(double V_excitation,
                                       double V_diff_min,
                                       double error_fraction_allowed)
{
    double V_cm = V_excitation / 2.0;
    double V_error_max = V_diff_min * error_fraction_allowed;
    if (V_error_max <= 0.0) return INFINITY;
    return V_cm / V_error_max;
}

/* ═══════════════════════════════════════════════════════════════════════
 * L2: Strain Calculation
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Compute mechanical strain from bridge output voltage.
 *
 * For quarter bridge: ε = (4·V_out / V_ex) / GF
 * For half bridge:    ε = (2·V_out / V_ex) / GF
 * For full bridge:    ε = (V_out / V_ex) / GF
 *
 * GF (gauge factor) = (ΔR/R₀) / ε
 *
 * Typical strain range:
 *   - Metal structures: ±2000 µε (microstrain = 10⁻⁶)
 *   - Composite materials: ±10000 µε
 *   - Yield point of steel: ~2000 µε
 */
double strain_from_bridge_output(double V_out, double V_ex,
                                  double GF, bridge_type_t type)
{
    if (V_ex <= 0.0 || GF <= 0.0) return NAN;

    double ratio = V_out / V_ex;

    switch (type) {
    case BRIDGE_QUARTER:
        return 4.0 * ratio / GF;
    case BRIDGE_HALF_ADJACENT:
    case BRIDGE_HALF_OPPOSITE:
        return 2.0 * ratio / GF;
    case BRIDGE_FULL:
        return ratio / GF;
    default:
        return NAN;
    }
}

/**
 * @brief Correct quarter-bridge nonlinearity.
 *
 * The exact relationship for quarter bridge:
 *   V_out/V_ex = ΔR/(4R₀+2ΔR)
 *
 * Inverting: ΔR/R₀ = 4(V_out/V_ex) / (1 - 2(V_out/V_ex))
 *
 * True strain: ε_true = (ΔR/R₀) / GF
 * Apparent strain (linear approx): ε_app = 4(V_out/V_ex)/GF
 *
 * Correction: ε_true = ε_app / (1 - GF·ε_app/2)
 *
 * For GF=2, ε=10000 µε: ε_app≈10000, ε_true≈10000/(1-0.01)=10101 µε
 * → 1% error if uncorrected at moderate strains
 *
 * For GF=2, ε=50000 µε (plastic region): error ≈ 5.3%
 */
double bridge_strain_nonlinearity_corrected(double strain_apparent,
                                             double gauge_factor)
{
    double denom = 1.0 - gauge_factor * strain_apparent / 2.0;
    if (fabs(denom) < 1e-6) return strain_apparent * 2.0;  /* degenerate */
    return strain_apparent / denom;
}
