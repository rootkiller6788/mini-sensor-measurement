/**
 * inst_amp_topology.c - IA Topology Implementations
 *
 * L2 Core Concepts: Complete implementations of all IA topologies
 * including 2-op-amp, current-mode, flying-capacitor, and indirect
 * current feedback architectures.
 *
 * Courses: Stanford EE247, Michigan EECS 411, Berkeley EE105
 * Reference:
 *   Van der Horn & Huijsing (1999) "Integrated Instrumentation Amplifiers"
 *   Toumazou, Lidgey, Haigh (1990) "Analogue IC Design: The Current-Mode Approach"
 *   Kitchin & Counts (2006) Ch. 4, 5
 */

#include "inst_amp_defs.h"
#include "inst_amp_topology.h"
#include <stdlib.h>
#include <math.h>

/* ==================================================================
 * L2: Two-Op-Amp IA Gain Computation
 *
 * Theorem (L4): For a 2-op-amp IA with R1=R4 and R2=R3:
 *   G = 1 + R4/R3 + 2*R4/R_G
 *
 * Derivation:
 *   A1 (non-inv): V_o1 = V_in- * (1 + R2/R1) - V_in+ * (R2/R1)
 *   A2 (diff):    V_out = V_in+ * (1 + R4/R3) - V_o1 * (R4/R3)
 *   Combined:     V_out = (V_in+ - V_in-) * (1 + R4/R3 + 2*R4/R_G)
 *
 * Note: This topology has asymmetric input impedance.
 *   Z_in+ = R3 + R4 || (R_G ...)  >> Z_in-
 * The input impedance mismatch can degrade CMRR when source
 * impedances are unbalanced.
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_two_opamp_gain(const ia_two_opamp_config_t *cfg) {
    if (!cfg) return 0.0;
    if (cfg->r3_ohm <= 0.0 || cfg->rg_ohm <= 0.0) return 0.0;

    double ratio = cfg->r4_ohm / cfg->r3_ohm;
    return 1.0 + ratio + 2.0 * cfg->r4_ohm / cfg->rg_ohm;
}

/* ==================================================================
 * L4: Two-Op-Amp CMRR Zero Frequency
 *
 * In the 2-op-amp topology, the common-mode signal takes different
 * paths through the two op-amp stages, creating a phase lag that
 * degrades CMRR at high frequencies.
 *
 * The CMRR reaches 0 dB (unity) at approximately:
 *   f_zero = GBW / (2 * pi * G)
 *
 * At this frequency, the CM gain equals the differential gain.
 * Above f_zero, the CM gain exceeds the diff gain (CMRR < 0 dB).
 *
 * For comparison, the 3-op-amp topology's CMRR zero is at:
 *   f_zero_3amp = GBW * (CMRR_res / G)
 * Which is typically 100-1000x higher.
 *
 * Reference: Kitchin & Counts (2006), Ch. 4, Fig 4-8
 * Complexity: O(1)
 * ================================================================== */
double ia_two_opamp_cmrr_zero_freq(const ia_two_opamp_config_t *cfg, double gbp_mhz) {
    if (!cfg) return 0.0;

    double G = ia_two_opamp_gain(cfg);
    if (G <= 0.0) return INFINITY;

    /* GBW in MHz -> Hz */
    return (gbp_mhz * 1e6) / (2.0 * M_PI * G);
}

/* ==================================================================
 * L2: Two-Op-Amp Input Impedance Calculation
 *
 * Asymmetric input impedance is a key limitation.
 *   Z_in+ = R3 + R4 (for typical values, ~20k-100k)
 *   Z_in- = R1 (typically 1k-10k)
 *
 * The lower Z_in- loads the sensor asymmetrically, which converts
 * common-mode to differential signal if source impedances differ.
 *
 * This is NOT a problem for the 3-op-amp topology where both
 * inputs see the op-amp's high input impedance.
 *
 * Complexity: O(1)
 * ================================================================== */
void ia_two_opamp_input_impedance(const ia_two_opamp_config_t *cfg,
                                   double *zin_plus, double *zin_minus) {
    if (!cfg || !zin_plus || !zin_minus) return;

    *zin_plus = cfg->r3_ohm + cfg->r4_ohm;  /* Approximate, ignoring R_G path */
    *zin_minus = cfg->r1_ohm;
}

/* ==================================================================
 * L2: Current-Mode IA Gain
 *
 * Principle (L3): Using a second-generation current conveyor (CCII+),
 * the input voltage is converted to current through R_G:
 *   I_diff = (V_in+ - V_in-) / R_G
 *
 * This current is then mirrored to the output and converted back
 * to voltage through R_L:
 *   V_out = I_diff * R_L
 *   G = R_L / R_G
 *
 * The gain is set by resistor ratio, independent of supply voltage.
 * This enables operation down to ~1.8V rails where voltage-mode
 * amplifiers cannot maintain sufficient headroom.
 *
 * The transconductance gm sets the bandwidth:
 *   BW = gm / (2*pi*C_comp)
 *
 * Reference: Toumazou, Lidgey, Haigh (1990), Ch. 3
 * Complexity: O(1)
 * ================================================================== */
double ia_current_mode_gain(const ia_current_mode_config_t *cfg) {
    if (!cfg) return 0.0;
    if (cfg->rg_ohm <= 0.0) return 0.0;

    return cfg->rl_ohm / cfg->rg_ohm;
}

/* ==================================================================
 * L2: Current-Mode IA Bandwidth Estimation
 *
 * BW is determined by the transconductance input stage:
 *   BW = gm / (2 * pi * C_comp)
 *
 * Where C_comp is the compensation capacitance (typically internal).
 * For gm = 1 mS and C_comp = 10 pF: BW ~ 16 MHz.
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_current_mode_bandwidth(const ia_current_mode_config_t *cfg,
                                  double c_comp_pf) {
    if (!cfg || c_comp_pf <= 0.0) return 0.0;
    return cfg->gm_ms * 1e-3 / (2.0 * M_PI * c_comp_pf * 1e-12);
}

/* ==================================================================
 * L2: Current-Mode IA Output Swing
 *
 * Unlike voltage-mode IAs, current-mode topologies can swing
 * very close to the rails (within 50-100mV) because the output
 * stage operates in the current domain.
 *
 * Output compliance voltage:
 *   V_out_max = V_cc - V_sat (typically V_cc - 0.1V)
 *   V_out_min = V_ee + V_sat (typically V_ee + 0.1V)
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_current_mode_output_max(const ia_current_mode_config_t *cfg,
                                   double v_cc, double v_sat) {
    (void)cfg;
    return v_cc - v_sat;
}

/* ==================================================================
 * L2: Flying-Capacitor IA Gain
 *
 * The gain is set by the capacitor ratio:
 *   G = C_in / C_fb
 *
 * Capacitor ratio accuracy can be 0.1% or better with on-chip
 * capacitors (much better than resistor matching).
 *
 * Sampling operation:
 *   phi1: C_in charged to V_in+ - V_in-
 *   phi2: Charge Q = C_in * V_diff transferred to C_fb
 *         V_out = Q / C_fb = (C_in/C_fb) * V_diff
 *
 * The switching frequency f_sw sets the effective bandwidth:
 *   BW_effective = f_sw / (2 * pi * G)
 *
 * Charge injection error from the switches:
 *   V_err = Q_injection / C_in
 *   where Q_injection is the switch charge injection (fC).
 *
 * Reference: Burt & Zhang, IEEE JSSC, 2006
 * Complexity: O(1)
 * ================================================================== */
double ia_flying_cap_gain(const ia_flying_cap_config_t *cfg) {
    if (!cfg) return 0.0;
    if (cfg->cfb_pf <= 0.0) return 0.0;

    return cfg->cin_pf / cfg->cfb_pf;
}

double ia_flying_cap_cmrr_ideal(void) {
    /* With perfectly matched capacitors on the same die,
     * CMRR can exceed 140 dB. Return linear value. */
    return 1.0e7;  /* 140 dB */
}

/* ==================================================================
 * L2: Flying-Capacitor Charge Injection Error
 *
 * When the MOS switches open, charge stored in the channel is
 * injected into the sampling capacitor, causing an offset error:
 *   V_offset = Q_injection / C_in
 *
 * For Q_inj = 10 fC and C_in = 1 pF: V_os = 10 uV
 *
 * This can be mitigated by:
 *   - Dummy switches (half-size, opposite phase)
 *   - Fully differential topology (CM rejection)
 *   - Correlated double sampling (CDS)
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_flying_cap_charge_injection_error(const ia_flying_cap_config_t *cfg) {
    if (!cfg || cfg->cin_pf <= 0.0) return 0.0;

    /* Charge injection in fC, capacitance in pF -> result in uV */
    return cfg->charge_injection_fc / cfg->cin_pf;
}

/* ==================================================================
 * L2: Indirect Current Feedback IA Gain
 *
 * This topology uses a transconductance input stage (Gm) and
 * a transimpedance feedback stage:
 *   G = Gm * R_feedback
 *
 * Key advantage: Bandwidth is approximately constant across
 * all gain settings because the feedback is in the current domain
 * rather than voltage domain.
 *
 * For Gm = 1 mS (1 mA/V) and R_f = 10k: G = 10 V/V
 * For Gm = 1 mS and R_f = 100k: G = 100 V/V
 *
 * Bandwidth ~ constant at ~GBW_Gm / G ... but since Gm stage
 * has very high bandwidth, the overall BW is limited by the
 * transimpedance stage pole: BW = 1 / (2*pi*R_f*C_parasitic)
 *
 * As R_f increases for higher gain, BW decreases proportionally,
 * but this is still better than VFB where BW ~ GBW/G.
 *
 * Reference: AD8421 datasheet, Analog Devices
 * Complexity: O(1)
 * ================================================================== */
double ia_indirect_current_gain(const ia_indirect_current_config_t *cfg) {
    if (!cfg) return 0.0;

    /* Gm in mS, convert to S; R_f in Ohm; G dimensionless */
    return cfg->gm_input_ms * 1e-3 * cfg->rf_ohm;
}

/* ==================================================================
 * L2: Indirect Current Feedback Bandwidth
 *
 * BW = 1 / (2 * pi * R_f * C_parasitic)
 * where C_parasitic is the total capacitance at the summing node.
 *
 * For R_f = 10k and C_par = 2 pF: BW = 8 MHz
 * For R_f = 100k and C_par = 2 pF: BW = 800 kHz
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_indirect_current_bandwidth(const ia_indirect_current_config_t *cfg,
                                      double c_parasitic_pf) {
    if (!cfg || c_parasitic_pf <= 0.0 || cfg->rf_ohm <= 0.0) return 0.0;

    return 1.0 / (2.0 * M_PI * cfg->rf_ohm * c_parasitic_pf * 1e-12);
}

/* ==================================================================
 * L3: Topology Comparison - CMRR vs Frequency
 *
 * Compares CMRR across topologies at a given frequency.
 *
 *   Topology     | CMRR_DC (dB) | CMRR Pole (Hz)
 *   -------------|-------------|----------------
 *   3-op-amp     | 80-120      | ~10-100 kHz (resistor limited)
 *   2-op-amp     | 70-90       | ~GBW/G
 *   Current-mode | 60-80       | ~gm/C (few MHz)
 *   Flying-cap   | 120-160     | ~f_sw/G
 *   Indirect CFB | 80-100      | ~1/(2*pi*R_f*C_par)
 *
 * Complexity: O(1)
 * ================================================================== */
double ia_topology_cmrr_at_freq(ia_topology_t topo, double cmrr_dc_db,
                                 double f_hz, double gbp_mhz, double gain) {
    double f_pole;

    switch (topo) {
        case IA_TOPO_THREE_OPAMP:
            /* Dominated by resistor mismatch, f_pole ~ 10-100kHz */
            f_pole = 50000.0;  /* 50 kHz typical */
            break;
        case IA_TOPO_TWO_OPAMP:
            f_pole = (gbp_mhz * 1e6) / gain;
            break;
        case IA_TOPO_CURRENT_MODE:
            f_pole = 1e6;  /* ~1 MHz for gm input stage */
            break;
        case IA_TOPO_FLYING_CAP:
            f_pole = 100e3;  /* Limited by switching frequency */
            break;
        case IA_TOPO_INDIRECT_CURRENT:
            f_pole = 500e3;  /* ~500 kHz */
            break;
        default:
            f_pole = 10000.0;
            break;
    }

    if (f_pole <= 0.0) return cmrr_dc_db;
    double rolloff = 10.0 * log10(1.0 + (f_hz / f_pole) * (f_hz / f_pole));
    return cmrr_dc_db - rolloff;
}

/* ==================================================================
 * L3: Gain-Bandwidth Trade-off Across Topologies
 *
 * For voltage-mode topologies (3-op-amp, 2-op-amp):
 *   BW_closed_loop = GBW / G
 *
 * For current-mode:
 *   BW ~ constant (gain-independent above some minimum)
 *
 * For indirect CFB:
 *   BW ~ constant for moderate gains, then drops at high R_f
 *
 * This function returns the estimated closed-loop bandwidth
 * for each topology at a given gain setting.
 * Complexity: O(1)
 * ================================================================== */
double ia_topology_bandwidth(ia_topology_t topo, double gain,
                              double gbp_mhz, double c_comp_pf) {
    switch (topo) {
        case IA_TOPO_THREE_OPAMP:
            return (gain > 0.0) ? (gbp_mhz * 1e6) / gain : INFINITY;
        case IA_TOPO_TWO_OPAMP:
            return (gain > 0.0) ? (gbp_mhz * 1e6) / (gain + 1.0) : INFINITY;
        case IA_TOPO_CURRENT_MODE:
            /* BW approximately independent of gain */
            return 1.0 / (2.0 * M_PI * c_comp_pf * 1e-12);
        case IA_TOPO_INDIRECT_CURRENT:
            /* Moderate BW, slightly gain-dependent */
            if (gain <= 0.0) return INFINITY;
            return 10e6 / sqrt(gain);  /* Empirical model */
        default:
            return (gain > 0.0) ? (gbp_mhz * 1e6) / gain : INFINITY;
    }
}

/* ==================================================================
 * L5: Auto-Range Gain Selection Algorithm
 *
 * For sensors with wide dynamic range (e.g., load cells from
 * 0 to full scale), auto-ranging selects the optimal gain
 * to maximize ADC resolution without saturation.
 *
 * Algorithm:
 *   if    V_in < V_fs / 1000 -> G = G_max
 *   elif  V_in < V_fs / 100  -> G = G_max / 10
 *   elif  V_in < V_fs / 10   -> G = G_max / 100
 *   else                      -> G = G_min
 *
 * Complexity: O(log(ratios)) = O(1)
 * ================================================================== */
double ia_autorange_gain(double v_in_rms, double v_fs_adc,
                          double g_min, double g_max, int num_ranges) {
    if (num_ranges < 1) num_ranges = 1;
    if (num_ranges > 8) num_ranges = 8;

    double g_step = pow(g_max / g_min, 1.0 / (num_ranges - 1));

    for (int r = 0; r < num_ranges; r++) {
        double g_test = g_min * pow(g_step, num_ranges - 1 - r);
        double v_max_input = v_fs_adc / g_test;
        if (v_in_rms <= v_max_input) return g_test;
    }

    return g_min;
}

/* ==================================================================
 * L5: Guard Band for Auto-Range Hysteresis
 *
 * Prevents gain hunting (rapid toggling) when input is near a
 * range boundary. Uses hysteresis: upward transition at a higher
 * threshold than downward transition.
 *
 *   G_up at V_in < V_thresh * 0.8  (reduce gain only when clearly needed)
 *   G_dn at V_in > V_thresh * 0.95
 *
 * Complexity: O(1)
 * ================================================================== */
bool ia_autorange_should_increase_gain(double v_current, double v_max_for_gain,
                                        double hysteresis) {
    return v_current < v_max_for_gain * 0.8 * hysteresis;
}

bool ia_autorange_should_decrease_gain(double v_current, double v_max_for_gain,
                                        double hysteresis) {
    return v_current > v_max_for_gain * 0.95 * hysteresis;
}