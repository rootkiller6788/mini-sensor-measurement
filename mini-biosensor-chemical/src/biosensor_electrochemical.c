/**
 * biosensor_electrochemical.c — Electrochemical biosensor implementations
 *
 * L2: Amperometry, potentiometry, impedimetry, voltammetry.
 * L4: Nernst equation, Butler-Volmer kinetics, Cottrell equation,
 *     Randles-Ševčík equation, Nikolskii-Eisenman equation.
 * L6: Glucose biosensor, pH-ISFET, ion-selective electrodes.
 * L7: Lactate biosensor for point-of-care.
 */

#include "biosensor_electrochemical.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Physical constants
 * ============================================================================ */

#define FARADAY_CONSTANT         96485.3329
#define GAS_CONSTANT             8.314462618
#define STANDARD_TEMPERATURE_K   298.15

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * L4: Nernst Equation
 * ============================================================================ */

/**
 * Nernst potential: E = E° + (RT/nF) · ln([Ox]/[Red])
 *
 * This is the thermodynamic foundation of all potentiometric biosensors.
 * At 25°C, (RT/F) = 25.693 mV, so for a 1-electron process the slope
 * is 59.16 mV per decade of concentration ratio change.
 */
double nernst_potential(const NernstParams *np,
                        double conc_oxidized,
                        double conc_reduced)
{
    if (!np) return 0.0;
    if (conc_reduced <= 0.0 || conc_oxidized <= 0.0) return np->standard_potential;

    double rt_over_nF = (np->gas_constant * np->temperature_kelvin) /
                        (np->electron_count * np->faraday_constant);
    double ratio = conc_oxidized / conc_reduced;
    if (ratio <= 0.0) return np->standard_potential;

    return np->standard_potential + rt_over_nF * log(ratio);
}

/* ============================================================================
 * L4/L6: ISFET pH sensor
 * ============================================================================ */

/**
 * ISFET pH computation from threshold voltage.
 *
 * The pH-sensitive insulator (Si₃N₄, Al₂O₃, Ta₂O₅) follows
 * site-binding theory, resulting in Nernstian or near-Nernstian
 * response. Ideal sensitivity = 59.16 mV/pH at 25°C.
 *
 * ΔV_th = sensitivity_mV/pH · (pH - pH_ref)
 *
 * pH = pH_ref + (V_th_measured - V_th_ref) / sensitivity
 */
double isfet_compute_ph(const ISFETParams *isfet, double v_th_measured)
{
    if (!isfet) return 7.0;

    /* Reference: V_th at pH 7 */
    double sensitivity_v_per_ph = isfet->sensitivity_mv_per_ph * 1.0e-3;
    if (fabs(sensitivity_v_per_ph) < 1.0e-12) return 7.0;

    double delta_v = v_th_measured - isfet->threshold_voltage;
    return 7.0 + delta_v / sensitivity_v_per_ph;
}

/* ============================================================================
 * L4: Butler-Volmer equation
 * ============================================================================ */

/**
 * Butler-Volmer: i = i₀ · [exp(α_a·n·F·η/(RT)) - exp(-α_c·n·F·η/(RT))]
 *
 * This is the fundamental equation describing electrode kinetics.
 * At small overpotential (|η| < 10 mV), linearizes to i ≈ i₀·n·F·η/(RT).
 * At large overpotential, one exponential term dominates (Tafel region).
 */
double butler_volmer_current(const ButlerVolmerParams *bv)
{
    if (!bv) return 0.0;

    double n = 1.0;  /* simplifying to 1-electron unless specified */
    double f_factor = n * FARADAY_CONSTANT /
                      (GAS_CONSTANT * STANDARD_TEMPERATURE_K);

    double eta = bv->overpotential;
    if (fabs(eta) < 1.0e-15) return 0.0;

    double anodic_term  = bv->anodic_transfer_coeff * f_factor * eta;
    double cathodic_term = bv->cathodic_transfer_coeff * f_factor * eta;

    double current_density = bv->exchange_current_density *
                             (exp(anodic_term) - exp(-cathodic_term));
    return current_density * bv->electrode_area;
}

/* ============================================================================
 * L4: Cottrell equation — chronoamperometry
 * ============================================================================ */

/**
 * Cottrell: I(t) = n·F·A·√D·C / √(π·t)
 *
 * Describes the diffusion-limited current following a potential step
 * to a value where the analyte is consumed at the maximum rate.
 *
 * The current decays with 1/√t because the diffusion layer thickness
 * grows as √(π·D·t), reducing the concentration gradient.
 *
 * Unit consistency:
 *   n: dimensionless, F: C/mol, A: cm², D: cm²/s,
 *   C: mol/cm³, t: s, result: C/s = A
 */
double cottrell_current(int n_electrons,
                        double electrode_area,
                        double diffusion_coeff,
                        double bulk_conc,
                        double time_sec)
{
    if (n_electrons <= 0 || electrode_area <= 0.0 ||
        diffusion_coeff <= 0.0 || bulk_conc <= 0.0 ||
        time_sec <= 0.0) return 0.0;

    double sqrt_pi_t = sqrt(M_PI * time_sec);
    return n_electrons * FARADAY_CONSTANT * electrode_area *
           sqrt(diffusion_coeff) * bulk_conc / sqrt_pi_t;
}

/* ============================================================================
 * L4: Randles-Ševčík equation — cyclic voltammetry
 * ============================================================================ */

/**
 * Randles-Ševčík: I_p = 2.69×10⁵ · n^(3/2) · A · √D · C · √v
 *
 * For a reversible electron transfer at 25°C.
 * The constant 2.69×10⁵ has units of C·mol⁻¹·V^(-1/2).
 *
 * I_p is the peak current in Amperes,
 * n is the number of electrons,
 * A is electrode area in cm²,
 * D is diffusion coefficient in cm²/s,
 * C is bulk concentration in mol/cm³,
 * v is scan rate in V/s.
 */
double randles_sevcik_peak_current(int n_electrons,
                                   double electrode_area,
                                   double diffusion_coeff,
                                   double bulk_conc,
                                   double scan_rate)
{
    if (n_electrons <= 0 || electrode_area <= 0.0 ||
        diffusion_coeff <= 0.0 || bulk_conc <= 0.0 ||
        scan_rate <= 0.0) return 0.0;

    double n_pow = pow((double)n_electrons, 1.5);
    return 2.69e5 * n_pow * electrode_area *
           sqrt(diffusion_coeff) * bulk_conc * sqrt(scan_rate);
}

/* ============================================================================
 * L2/L6: Glucose biosensor — amperometric detection
 * ============================================================================ */

/**
 * Glucose biosensor current model.
 *
 * GOx catalyzes: Glucose + O₂ → Gluconic acid + H₂O₂
 * H₂O₂ electro-oxidation at Pt electrode (+0.6V vs Ag/AgCl):
 *   H₂O₂ → O₂ + 2H⁺ + 2e⁻
 *
 * The response is modeled as a combination of:
 * 1. Linear regime: I = sensitivity · [Glucose]   (when [Glu] << K_m)
 * 2. Michaelis-Menten saturation at high concentrations
 * 3. O₂ co-substrate limitation: if [O₂] < [Glucose], O₂ is rate-limiting
 *
 * Reference: Clark & Lyons (1962), Ann. NY Acad. Sci.
 */
double glucose_biosensor_current(double glucose_conc_mM,
                                 double enzyme_activity_u,
                                 double electrode_area_cm2,
                                 double oxygen_conc_mM)
{
    if (glucose_conc_mM < 0.0 || enzyme_activity_u <= 0.0 ||
        electrode_area_cm2 <= 0.0) return 0.0;

    /* Typical GOx K_m for glucose ≈ 5-10 mM (immobilized) */
    double km_glucose = 8.0;   /* mM, apparent K_m for immobilized GOx */
    double km_oxygen  = 0.5;   /* mM, K_m for O₂ as co-substrate */

    /* Enzyme-limited maximum current = k_cat · [E_total] · electrode coverage */
    /* Typical: 1 U/cm² GOx produces ~500 nA/(cm²·mM) sensitivity */
    double sensitivity = enzyme_activity_u * electrode_area_cm2 * 500.0; /* nA/mM */

    /* Michaelis-Menten for glucose at given O₂ level */
    double oxygen_factor = oxygen_conc_mM / (km_oxygen + oxygen_conc_mM);
    double glucose_factor = glucose_conc_mM / (km_glucose + glucose_conc_mM);

    /* Combined rate limited by the slower of the two substrates */
    double vmax_effective = sensitivity * km_glucose * oxygen_factor;
    return vmax_effective * glucose_factor;
}

/**
 * Inverse model: estimate glucose concentration from measured current.
 */
double glucose_conc_from_current(double current_na,
                                 double sensitivity_na_per_mM,
                                 double km_apparent_mM,
                                 double imax_na)
{
    if (sensitivity_na_per_mM <= 0.0 || current_na < 0.0) return 0.0;

    /* Linear approximation for low currents */
    double conc_linear = current_na / sensitivity_na_per_mM;

    /* Michaelis-Menten correction for saturation */
    if (imax_na > 0.0 && current_na > 0.5 * imax_na) {
        double frac = current_na / imax_na;
        if (frac >= 1.0) frac = 0.999;
        conc_linear = km_apparent_mM * frac / (1.0 - frac);
    }

    return conc_linear;
}

/* ============================================================================
 * L2: EIS — Randles equivalent circuit
 * ============================================================================ */

/**
 * Randles circuit impedance: Z(ω) = R_s + (R_ct + Z_w) || C_dl
 *
 * Z_w = σ/√(ω) · (1 - j)  (Warburg impedance, semi-infinite diffusion)
 *
 * Real:  Z' = R_s + (R_ct + σ/√ω) / [(1 + σ·ω^{1/2}·C_dl)² + ω²·C_dl²·(R_ct + σ/√ω)²]
 *
 * We compute the parallel combination directly:
 *   Z = R_s + 1 / (1/(R_ct + Z_w) + jωC_dl)
 */
void eis_randles_impedance(const EISParams *eis,
                           double freq_hz,
                           double *z_real,
                           double *z_imag)
{
    if (!eis || !z_real || !z_imag) return;

    double omega = 2.0 * M_PI * freq_hz;
    double sqrt_omega = sqrt(omega);

    /* Warburg impedance */
    double sigma_over_sqrt_omega = eis->warburg_coefficient / sqrt_omega;
    double rct_plus_zw_real = eis->charge_transfer_resistance + sigma_over_sqrt_omega;
    double zw_imag = -sigma_over_sqrt_omega;

    /* Parallel combination: (R_ct + Z_w) || C_dl */
    /* admittance = 1/(R_ct + Z_w) + jωC_dl */
    double denom_complex = rct_plus_zw_real * rct_plus_zw_real +
                           zw_imag * zw_imag;
    if (denom_complex < 1.0e-30) {
        *z_real = eis->solution_resistance;
        *z_imag = 0.0;
        return;
    }

    double adm_real = rct_plus_zw_real / denom_complex;
    double adm_imag = -zw_imag / denom_complex + omega * eis->double_layer_capacitance;

    /* Z_parallel = 1 / Y_parallel */
    double adm_mag_sq = adm_real * adm_real + adm_imag * adm_imag;
    if (adm_mag_sq < 1.0e-30) {
        *z_real = eis->solution_resistance;
        *z_imag = 0.0;
        return;
    }

    double zp_real = adm_real / adm_mag_sq;
    double zp_imag = -adm_imag / adm_mag_sq;

    *z_real = eis->solution_resistance + zp_real;
    *z_imag = zp_imag;
}

/**
 * Estimate R_ct from Nyquist data by fitting the high-frequency semicircle.
 *
 * The Nyquist semicircle in the Randles model has center at (R_s + R_ct/2, 0)
 * and radius R_ct/2. We estimate R_s from the high-frequency x-intercept
 * and R_ct from the width of the arc.
 *
 * Methodology:
 * 1. Identify the high-frequency region where -Z″ is increasing
 * 2. Fit a circle to the Nyquist points
 * 3. R_s = x-intercept of the circle, R_ct = diameter
 */
double eis_fit_charge_transfer_resistance(const double *freq_array,
                                          const double *z_real_array,
                                          const double *z_imag_array,
                                          int n_points,
                                          double *rs_est,
                                          double *rct_est)
{
    if (!freq_array || !z_real_array || !z_imag_array ||
        !rs_est || !rct_est || n_points < 5) {
        if (rs_est)  *rs_est  = 0.0;
        if (rct_est) *rct_est = 0.0;
        return -1.0;
    }

    /* Find the high-frequency region: first few points with decreasing -Z″ */
    int hf_end = 0;
    double max_neg_zimag = 0.0;
    for (int i = 0; i < n_points; i++) {
        double neg_zimag = -z_imag_array[i];
        if (neg_zimag > max_neg_zimag) {
            max_neg_zimag = neg_zimag;
            hf_end = i;
        }
    }
    if (hf_end < 3) hf_end = n_points / 2;
    if (hf_end > n_points) hf_end = n_points;

    /* Simple chord method: the Nyquist semicircle's x-intercepts give
     * R_s (left) and R_s + R_ct (right).
     *
     * At high frequencies, Z″ → 0, Z′ → R_s.
     * We estimate R_s from the minimum Z′ at the highest frequencies. */
    double rs_sum = 0.0;
    int rs_count = 0;
    /* Take the 3 highest frequency points */
    int high_freq_count = (n_points >= 3) ? 3 : n_points;
    for (int i = 0; i < high_freq_count; i++) {
        rs_sum += z_real_array[i];
        rs_count++;
    }
    *rs_est = rs_sum / rs_count;

    /* Find the lowest frequency point in the semicircle region
     * (near the right x-intercept of the semicircle where Z″ returns to 0) */
    double max_zreal = 0.0;
    for (int i = hf_end / 2; i < n_points; i++) {
        if (z_real_array[i] > max_zreal) {
            max_zreal = z_real_array[i];
        }
    }

    *rct_est = max_zreal - *rs_est;
    if (*rct_est < 0.0) *rct_est = 0.0;

    /* RMS error: compute deviation from ideal semicircle
     * with center at (R_s + R_ct/2, 0), radius = R_ct/2 */
    double center_x = *rs_est + *rct_est / 2.0;
    double radius  = *rct_est / 2.0;
    double sum_sq_err = 0.0;
    int fit_count = 0;

    for (int i = 0; i < hf_end && i < n_points; i++) {
        /* Only consider points with negative Z″ (capacitive) */
        if (z_imag_array[i] >= 0.0) continue;
        double dx = z_real_array[i] - center_x;
        double dy = z_imag_array[i];  /* negative = below x-axis */
        double r_measured = sqrt(dx * dx + dy * dy);
        double err = r_measured - radius;
        sum_sq_err += err * err;
        fit_count++;
    }

    if (fit_count > 0 && radius > 1.0e-12) {
        return sqrt(sum_sq_err / fit_count);
    }
    return 0.0;
}

/* ============================================================================
 * L4/L6: Nikolskii-Eisenman equation — Ion-Selective Electrode
 * ============================================================================ */

/**
 * Nikolskii-Eisenman: E = E° + (RT/z_i·F)·ln(a_i + Σ K_ij·a_j^(z_i/z_j))
 *
 * This equation generalizes the Nernst equation for ion-selective
 * electrodes (ISEs) by accounting for interfering ions through
 * selectivity coefficients K_ij.
 *
 * For a perfectly selective electrode, K_ij = 0 for all j,
 * reducing to the simple Nernst equation.
 */
double nikolskii_eisenman_potential(const NernstParams *np,
                                    double primary_activity,
                                    const double *interferent_activities,
                                    const double *selectivity_coeffs,
                                    const int *interferent_charges,
                                    int n_interferents)
{
    if (!np) return 0.0;

    double interference_sum = 0.0;

    if (interferent_activities && selectivity_coeffs &&
        interferent_charges && n_interferents > 0) {
        for (int j = 0; j < n_interferents; j++) {
            if (interferent_charges[j] == 0) continue;
            double charge_ratio = (double)np->electron_count /
                                  (double)labs(interferent_charges[j]);
            double activity_term = pow(interferent_activities[j], charge_ratio);
            interference_sum += selectivity_coeffs[j] * activity_term;
        }
    }

    double total_activity = primary_activity + interference_sum;
    if (total_activity <= 0.0) return np->standard_potential;

    double rt_over_zF = (np->gas_constant * np->temperature_kelvin) /
                        (fabs((double)np->electron_count) * np->faraday_constant);

    return np->standard_potential + rt_over_zF * log(total_activity);
}

/* ============================================================================
 * L3: Debye-Hückel activity coefficient
 * ============================================================================ */

/**
 * Extended Debye-Hückel: log₁₀(γ) = -A·z²·√I/(1 + B·a·√I)
 *
 * Constants for water at 25°C: A = 0.509, B = 0.328 (mol/L basis, Å)
 *
 * Used to convert ionic concentrations to thermodynamic activities
 * for accurate ISE and potentiometric biosensor measurements.
 */
double debye_huckel_activity_coeff(int ion_charge,
                                   double ionic_strength,
                                   double ion_size_param)
{
    if (ionic_strength <= 0.0 || ion_size_param <= 0.0) return 1.0;

    /* Debye-Hückel constants for water at 25°C (molality scale) */
    const double A = 0.509;  /* (L/mol)^(1/2) */
    const double B = 0.328;  /* (L/mol)^(1/2) · Å⁻¹ */

    double z_sq = (double)(ion_charge * ion_charge);
    double sqrt_I = sqrt(ionic_strength);

    double numerator   = -A * z_sq * sqrt_I;
    double denominator = 1.0 + B * ion_size_param * sqrt_I;

    double log_gamma = numerator / denominator;
    return pow(10.0, log_gamma);
}

/* ============================================================================
 * L7: Lactate biosensor for point-of-care diagnostics
 * ============================================================================ */

/**
 * Lactate biosensor model.
 *
 * Blood lactate is a critical biomarker for:
 * - Sepsis severity (lactate > 2 mM = abnormal)
 * - Athletic performance (lactate threshold testing)
 * - Tissue hypoxia / shock
 *
 * The sensor uses lactate oxidase (LOx) immobilized on a Pt electrode.
 * Temperature and hematocrit corrections are essential for whole blood accuracy.
 */
double lactate_biosensor_current(double lactate_conc_mM,
                                 double temp_celsius,
                                 double hematocrit_percent)
{
    if (lactate_conc_mM < 0.0) return 0.0;

    /* Lactate oxidase K_m ≈ 0.5 mM for L-lactate */
    double km_lactate = 0.5;
    double vmax_25c   = 100.0;  /* nA at 25°C for typical sensor geometry */

    /* Temperature correction using Arrhenius: Q₁₀ ≈ 1.5 for LOx */
    double q10 = 1.5;
    double temp_factor = pow(q10, (temp_celsius - 25.0) / 10.0);

    /* Michaelis-Menten rate */
    double rate = vmax_25c * temp_factor *
                  lactate_conc_mM / (km_lactate + lactate_conc_mM);

    /* Hematocrit correction: higher Hct → lower plasma volume → lower current
     * Correction factor from i-STAT clinical validation data:
     * CF ≈ 1 + 0.008·(Hct - 42) for Hct 15-65% */
    double hct_correction = 1.0 + 0.008 * (hematocrit_percent - 42.0);
    if (hct_correction < 0.7)  hct_correction = 0.7;
    if (hct_correction > 1.5)  hct_correction = 1.5;

    return rate / hct_correction;
}
