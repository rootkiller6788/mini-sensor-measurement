/**
 * @file sensor_bridge.h
 * @brief Wheatstone bridge circuits, signal conditioning, and bridge analysis
 *
 * Level: L2 — Core Concepts · L3 — Mathematical Structures
 *
 * The Wheatstone bridge is the most widely used circuit for precision
 * measurement of small resistance changes (strain gauges, RTDs, thermistors,
 * pressure sensors, load cells, magnetoresistors).
 *
 * Bridge Configurations:
 *
 *   Quarter Bridge (1 active arm):
 *     V_out = V_ex · ΔR / (4R₀ + 2ΔR) ≈ V_ex · ΔR / (4R₀)  [if ΔR << R₀]
 *
 *   Half Bridge (2 active arms, adjacent):
 *     V_out = V_ex · ΔR / (2R₀)   [double sensitivity, temperature compensated]
 *
 *   Full Bridge (4 active arms):
 *     V_out = V_ex · ΔR / R₀       [maximum sensitivity, temperature compensated]
 *
 *   Half Bridge (2 active arms, opposite):
 *     V_out = 0 for symmetric ΔR → measures bending, not axial strain
 *
 * Key Formulas:
 *   V_out = V_ex · (R₃/(R₂+R₃) - R₄/(R₁+R₄))
 *
 *   Balanced condition: R₁/R₂ = R₄/R₃  ↔  V_out = 0
 *
 * Key Effects:
 *   - Lead wire resistance: significant in quarter bridge with long leads
 *   - Thermal EMF: DC offset from dissimilar metal junctions
 *   - Self-heating: sensor excitation current causes temperature rise
 *   - Nonlinearity: exact bridge equation is nonlinear for large ΔR
 *   - Common-mode voltage: requires differential amplifier with high CMRR
 *
 * References:
 *   - Wheatstone, Phil. Trans. Roy. Soc., 133:303-327, 1843
 *   - Window, "Strain Gauge Technology", 2nd Ed, Ch. 4
 *   - Hoffmann, "Applying the Wheatstone Bridge Circuit", HBM, 2001
 */

#ifndef SENSOR_BRIDGE_H
#define SENSOR_BRIDGE_H

#include <stddef.h>

/* ──── L2: Wheatstone Bridge Types ──── */

typedef enum {
    BRIDGE_QUARTER = 0,     /* 1 active element, 3 fixed (completion) resistors */
    BRIDGE_HALF_ADJACENT = 1, /* 2 active elements, adjacent arms */
    BRIDGE_HALF_OPPOSITE = 2, /* 2 active elements, opposite arms */
    BRIDGE_FULL = 3,          /* 4 active elements */
} bridge_type_t;

/**
 * @brief Wheatstone bridge configuration and state.
 *
 *           R1 (top-left)     R2 (top-right)
 *         ┌───[ R1 ]───┬───[ R2 ]───┐
 *         │   + V_ex -  │            │
 *         │             ├── V_out ──┤
 *         │             │            │
 *         └───[ R4 ]───┴───[ R3 ]───┘
 *           R4 (bottom-left) R3 (bottom-right)
 *
 * V_out is measured between the two midpoints (right arm midpoint minus left).
 */
typedef struct {
    double R1;  /* [Ω] top-left arm */
    double R2;  /* [Ω] top-right arm */
    double R3;  /* [Ω] bottom-right arm */
    double R4;  /* [Ω] bottom-left arm */
    double V_excitation; /* [V] excitation voltage across bridge */
    bridge_type_t type;
} wheatstone_bridge_t;

/* ──── L2: Bridge Circuit Functions ──── */

/** Initialize a balanced bridge with nominal resistance R0 in all arms. */
int    bridge_init_balanced(wheatstone_bridge_t *b, double R0, double V_ex);

/** Set one arm to active resistance, other three to nominal R0 (quarter bridge). */
int    bridge_set_quarter_active(wheatstone_bridge_t *b,
                                  double R_nominal, double R_sensed,
                                  int arm_index);  /* 1-4, arm to modify */

/** Set two adjacent arms (half bridge, temperature-compensated). */
int    bridge_set_half_adjacent(wheatstone_bridge_t *b,
                                 double R_nominal,
                                 double R_upper_active, double R_lower_active);

/** Set four arms in push-pull configuration (full bridge). */
int    bridge_set_full_pushpull(wheatstone_bridge_t *b,
                                 double R_nominal, double delta_R);

/* ──── L3: Bridge Output Voltage ──── */

/**
 * @brief Compute exact bridge output voltage (no approximations).
 *
 * V_out = V_ex · [R₃/(R₂+R₃) - R₄/(R₁+R₄)]
 *
 * This is the exact nonlinear equation valid for all resistance values.
 */
double bridge_output_voltage_exact(const wheatstone_bridge_t *b);

/**
 * @brief Compute bridge output voltage with the small-ΔR approximation.
 *
 * For quarter bridge: V_out ≈ V_ex · ΔR / (4·R₀)
 * For half bridge:    V_out ≈ V_ex · ΔR / (2·R₀)
 * For full bridge:    V_out ≈ V_ex · ΔR / R₀
 *
 * Valid when |ΔR| / R₀ < 0.01 (error < 0.5%).
 */
double bridge_output_voltage_approx(const wheatstone_bridge_t *b, double R0,
                                     double delta_R);

/**
 * @brief Compute bridge nonlinearity error.
 *
 * Error = |V_exact - V_approx| / V_exact · 100%
 *
 * For quarter bridge with ΔR/R₀ = 0.1, nonlinearity ≈ 5%.
 */
double bridge_nonlinearity_error(const wheatstone_bridge_t *b, double R0,
                                  double delta_R);

/**
 * @brief Compute bridge sensitivity [V/Ω] at the operating point.
 *
 * S = dV_out/dR_active evaluated numerically at current resistances.
 */
double bridge_sensitivity(const wheatstone_bridge_t *b, int active_arm);

/**
 * @brief Compute bridge equivalent (Thevenin) resistance seen by the
 *        measurement amplifier input.
 *
 * R_th = (R₁ ∥ R₄) + (R₂ ∥ R₃)
 *      = R₁·R₄/(R₁+R₄) + R₂·R₃/(R₂+R₃)
 *
 * This determines the loading effect and noise contribution of the amplifier.
 */
double bridge_thevenin_resistance(const wheatstone_bridge_t *b);

/* ──── L2: Bridge Non-Ideal Effects ──── */

/**
 * @brief Estimate lead wire voltage drop effect in 2-wire measurement.
 *
 * For quarter bridge with lead resistance R_lead in each of the two measurement leads:
 *   Error ≈ 2·R_lead / (R_gauge + 2·R_lead)  [fractional error]
 *
 * 3-wire and 4-wire (Kelvin) connections eliminate this error.
 */
double bridge_lead_wire_error_2wire(double R_gauge, double R_lead_per_wire);

/**
 * @brief Estimate thermal EMF offset from parasitic thermocouple junctions.
 *
 * Cu-PbSn solder: ~1-3 µV/°C
 * Kovar-Cu: ~40 µV/°C
 *
 * Thermal EMF can be suppressed by AC excitation or current-reversal techniques.
 */
double bridge_thermal_emf_offset(double seebeck_coeff_uV_per_K,
                                  double delta_T_K);

/**
 * @brief Bridge self-heating power dissipation per active arm.
 *
 * P_arm = I²·R = (V_excitation / R_total)² · R_arm
 *
 * For quarter bridge: P ≈ V_ex² / (4·R₀)
 *
 * Self-heating causes measurement error proportional to the thermal
 * resistance of the sensor to ambient.
 */
double bridge_arm_power_dissipation(const wheatstone_bridge_t *b, int arm_index);

/**
 * @brief Estimate temperature rise from self-heating.
 *
 * ΔT = P_dissipated · R_th
 * where R_th [K/W] is the thermal resistance to ambient.
 *
 * For typical strain gauges on metal: R_th ≈ 10-50 K/W
 * For RTDs in still air: R_th ≈ 50-200 K/W
 */
double bridge_self_heating_temp_rise(double power_per_arm, double thermal_resistance);

/* ──── L3: Bridge Signal Conditioning & Amplifier ──── */

/**
 * @brief Instrumentation amplifier gain required for bridge signal.
 *
 * Required gain: G = V_adc_fullscale / V_bridge_fullscale
 *
 * Common-mode voltage at amplifier input:
 *   V_cm = V_ex/2  (for balanced bridge)
 *   V_cm varies as bridge unbalances (second-order effect)
 *
 * CMRR requirement: CMRR_dB > 20·log₁₀(V_cm / V_error_allowed)
 */
double bridge_amplifier_gain_required(double V_bridge_fs, double V_adc_fs);

/**
 * @brief Compute required CMRR for bridge amplifier.
 *
 * Common-mode voltage V_cm ≈ V_ex/2.
 * Differential signal V_diff << V_cm.
 *
 * Rejection needed: CMRR_min = V_cm / (V_diff · ε_allowed)
 * e.g., V_ex=5V, V_diff=25mV, ε=0.1% → CMRR > 86 dB
 */
double bridge_amplifier_cmrr_required(double V_excitation,
                                       double V_diff_min,
                                       double error_fraction_allowed);

/* ──── L2: Bridge Calibration Constants ──── */

/**
 * @brief Standard gauge factor definitions.
 *
 * Gauge Factor GF = (ΔR/R₀) / ε
 *   where ε = ΔL/L is mechanical strain.
 *
 * Typical GF values:
 *   - Constantan (Cu-Ni): GF ≈ 2.0 (most common)
 *   - Karma (Ni-Cr): GF ≈ 2.1
 *   - Isoelastic (Fe-Ni-Cr-Mo): GF ≈ 3.6
 *   - Platinum-tungsten: GF ≈ 4.0
 *   - Semiconductor (Si): GF ≈ ±50 to ±200 (but large temp coefficient)
 *
 * For semiconductor gauges, GF temperature coefficient can be 0.1-1%/°C,
 * so temperature compensation is essential.
 */
double strain_from_bridge_output(double V_out, double V_ex,
                                  double GF, bridge_type_t type);

/**
 * @brief Compute nonlinearity correction for quarter bridge.
 *
 * True strain ε_true = ε_approx / (1 - GF·ε_approx/2)
 *                     ≈ ε_approx · (1 + GF·ε_approx/2)
 *
 * For GF=2, ε=10000 µε: error ≈ 1% if uncorrected.
 */
double bridge_strain_nonlinearity_corrected(double strain_apparent,
                                             double gauge_factor);

#endif /* SENSOR_BRIDGE_H */
