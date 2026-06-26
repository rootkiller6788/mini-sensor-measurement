/**
 * biosensor_optical.c — Optical biosensor implementations
 *
 * L2: SPR, fluorescence, absorbance, evanescent wave sensing.
 * L4: Beer-Lambert law, FRET, Fresnel equations, Snell's law.
 * L6: ELISA quantification, DNA microarray, fiber-optic biosensor.
 * L7: DNA microarray fluorescence analysis for gene expression.
 */

#include "biosensor_optical.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * L4: Beer-Lambert Law
 * ============================================================================ */

/**
 * Beer-Lambert: A = ε · b · c
 *
 * The absorbance A (also called Optical Density, OD) is dimensionless.
 * ε is the molar extinction coefficient in M⁻¹·cm⁻¹,
 * b is the path length in cm,
 * c is the concentration in mol/L (M).
 *
 * This law is linear only for dilute solutions (typically A < 1.0)
 * and monochromatic light. Deviations occur due to:
 * - Stray light in the spectrophotometer
 * - Chemical equilibrium shifts
 * - High concentration (intermolecular interactions)
 */
double beer_lambert_absorbance(double extinction_coeff,
                               double path_length_cm,
                               double concentration_m)
{
    if (extinction_coeff < 0.0 || path_length_cm <= 0.0 ||
        concentration_m < 0.0) return 0.0;
    return extinction_coeff * path_length_cm * concentration_m;
}

/**
 * Inverse Beer-Lambert: c = A / (ε · b)
 */
double concentration_from_absorbance(double absorbance,
                                     double extinction_coeff,
                                     double path_length_cm)
{
    if (absorbance < 0.0 || extinction_coeff <= 0.0 ||
        path_length_cm <= 0.0) return 0.0;
    return absorbance / (extinction_coeff * path_length_cm);
}

/**
 * Two-component spectrophotometric analysis.
 *
 * For a binary mixture, the absorbances at two wavelengths
 * form a 2×2 linear system:
 *   [A₁] = [ε_a1  ε_b1] · [c_a]
 *   [A₂]   [ε_a2  ε_b2]   [c_b]  · b
 *
 * Solution via Cramer's rule:
 *   det = ε_a1·ε_b2 - ε_a2·ε_b1
 *   c_a = (A₁·ε_b2 - A₂·ε_b1) / (b · det)
 *   c_b = (A₂·ε_a1 - A₁·ε_a2) / (b · det)
 */
bool two_component_absorbance(double abs_lambda1, double abs_lambda2,
                              double eps_a1, double eps_a2,
                              double eps_b1, double eps_b2,
                              double path_length,
                              double *conc_a, double *conc_b)
{
    if (!conc_a || !conc_b) return false;
    if (path_length <= 0.0) return false;

    double det = eps_a1 * eps_b2 - eps_a2 * eps_b1;

    /* Check condition number: poorly conditioned if |det| is small */
    if (fabs(det) < 1.0e-12) {
        *conc_a = 0.0;
        *conc_b = 0.0;
        return false;
    }

    *conc_a = (abs_lambda1 * eps_b2 - abs_lambda2 * eps_b1) /
              (path_length * det);
    *conc_b = (abs_lambda2 * eps_a1 - abs_lambda1 * eps_a2) /
              (path_length * det);

    /* Physical constraint: concentrations cannot be negative */
    if (*conc_a < -1.0e-12 || *conc_b < -1.0e-12) {
        if (*conc_a < 0.0) *conc_a = 0.0;
        if (*conc_b < 0.0) *conc_b = 0.0;
        return false;
    }

    return true;
}

/* ============================================================================
 * L2: Surface Plasmon Resonance (SPR)
 * ============================================================================ */

/**
 * SPR resonance angle in the Kretschmann configuration.
 *
 * Surface plasmon wavevector:
 *   k_SP = (2π/λ) · √(ε_m · ε_d / (ε_m + ε_d))
 *
 * Where ε_m = ε_m_real (gold permittivity at λ, imaginary part small),
 * ε_d = n_sample² (sample dielectric constant).
 *
 * Resonance condition: k_SP = k_photon∥ = (2π/λ) · n_prism · sin(θ)
 *   → θ_SPR = arcsin( √(ε_m·ε_d/(ε_m+ε_d)) / n_prism )
 *
 * Valid when ε_m_real < 0 and |ε_m_real| > ε_d (gold at visible: ε_m ≈ -12 at 633 nm).
 */
double spr_resonance_angle(const SPRParams *spr,
                           double sample_refractive_index,
                           double gold_permittivity_real)
{
    if (!spr || sample_refractive_index <= 1.0 ||
        gold_permittivity_real >= 0.0) return 0.0;

    double eps_d = sample_refractive_index * sample_refractive_index;
    double eps_m = gold_permittivity_real;  /* negative for metals at visible */

    double eps_product = eps_m * eps_d;
    double eps_sum     = eps_m + eps_d;

    if (eps_sum >= 0.0) return 0.0;  /* total internal reflection regime */

    double k_sp_factor = sqrt(eps_product / eps_sum);
    double sin_theta = k_sp_factor / spr->prism_ri;

    if (sin_theta <= 0.0 || sin_theta >= 1.0) return 0.0;

    return asin(sin_theta) * 180.0 / M_PI;
}

/**
 * SPR reflectivity via Fresnel 3-layer model.
 *
 * Transfer matrix method for p-polarized light:
 * Each interface j has Fresnel coefficient:
 *   r_jk = (k_zj/ε_j - k_zk/ε_k) / (k_zj/ε_j + k_zk/ε_k)
 *
 * Full system reflectivity R = |r_total|².
 *
 * k_zj = √(ε_j·(2π/λ)² - k_∥²)  where k_∥ = n_prism·(2π/λ)·sin(θ)
 */
double spr_reflectivity(const SPRParams *spr,
                        double angle_deg,
                        double sample_ri,
                        double gold_ri_real,
                        double gold_ri_imag)
{
    if (!spr) return 1.0;

    double theta_rad = angle_deg * M_PI / 180.0;
    double k0 = 2.0 * M_PI / (spr->wavelength_nm * 1.0e-9);
    double k_parallel = spr->prism_ri * k0 * sin(theta_rad);

    /* Layer 0: prism, layer 1: gold, layer 2: sample */
    double eps_prism = spr->prism_ri * spr->prism_ri;
    double eps_gold_real = gold_ri_real * gold_ri_real - gold_ri_imag * gold_ri_imag;
    double eps_gold_imag = 2.0 * gold_ri_real * gold_ri_imag;

    double eps_sample = sample_ri * sample_ri;

    /* z-component of wavevector for each layer: k_zj = √(ε_j·k₀² - k_∥²) */
    double kz_prism_sq = eps_prism * k0 * k0 - k_parallel * k_parallel;
    double kz_gold_sq  = eps_gold_real * k0 * k0 - k_parallel * k_parallel;
    double kz_sample_sq = eps_sample * k0 * k0 - k_parallel * k_parallel;

    if (kz_prism_sq <= 0.0) return 1.0;  /* beyond critical angle? */

    double kz_prism = sqrt(kz_prism_sq);

    /* Gold layer: complex k_z */
    double kz_gold_real, kz_gold_imag;
    {
        double a = kz_gold_sq;
        double b = eps_gold_imag * k0 * k0;
        double mag = sqrt(a * a + b * b);
        kz_gold_real = sqrt((mag + a) / 2.0);
        kz_gold_imag = (b >= 0.0 ? 1.0 : -1.0) * sqrt((mag - a) / 2.0);
    }

    /* Sample layer */
    double kz_sample;
    if (kz_sample_sq > 0.0) {
        kz_sample = sqrt(kz_sample_sq);
    } else {
        /* Evanescent in sample */
        kz_sample = 0.0;  /* lossy sample approximation */
    }

    /* Fresnel reflection coefficients (p-polarization) */
    /* Interface prism→gold: r_01 = (k_z0/ε₀ - k_z1/ε₁) / (k_z0/ε₀ + k_z1/ε₁) */
    double denom_01_real = kz_prism / eps_prism + kz_gold_real / eps_gold_real;
    double numer_01_real = kz_prism / eps_prism - kz_gold_real / eps_gold_real;

    /* Simplified: treat gold as having real effective RI for amplitude */
    double r_01 = numer_01_real / denom_01_real;

    /* Phase shift through gold layer: φ = k_z1 · d */
    double phase = kz_gold_real * spr->gold_layer_thickness_nm * 1.0e-9;
    double exp_factor = exp(-2.0 * kz_gold_imag * spr->gold_layer_thickness_nm * 1.0e-9);

    /* Interface gold→sample */
    double denom_12 = kz_gold_real / eps_gold_real + kz_sample / eps_sample;
    double numer_12 = kz_gold_real / eps_gold_real - kz_sample / eps_sample;
    double r_12 = (fabs(denom_12) < 1.0e-30) ? 1.0 : numer_12 / denom_12;

    /* Total reflection coefficient: r = (r_01 + r_12·e^(2i·k_z1·d)) /
     *                                      (1 + r_01·r_12·e^(2i·k_z1·d)) */
    double cos_2phase = cos(2.0 * phase);
    double sin_2phase = sin(2.0 * phase);

    double num_real = r_01 + r_12 * exp_factor * cos_2phase;
    double num_imag = r_12 * exp_factor * sin_2phase;
    double den_real = 1.0 + r_01 * r_12 * exp_factor * cos_2phase;
    double den_imag = r_01 * r_12 * exp_factor * sin_2phase;

    double den_mag_sq = den_real * den_real + den_imag * den_imag;
    if (den_mag_sq < 1.0e-30) return 1.0;

    double r_mag_sq = (num_real * num_real + num_imag * num_imag) / den_mag_sq;
    return (r_mag_sq > 1.0) ? 1.0 : r_mag_sq;
}

/**
 * Adlayer thickness estimation from SPR angle shift.
 *
 * Δθ ∝ (dn/dc) · Γ · (1 - e^(-2d/l_d))
 *
 * Surface mass density Γ = d · ρ_adlayer
 * For protein films: ρ ≈ 1.37 g/cm³, dn/dc ≈ 0.182 mL/g
 */
double spr_adlayer_thickness(double angle_shift_mdeg,
                             double dn_dc,
                             double molecular_weight,
                             double sensitivity_mdeg_per_riu,
                             double decay_length_nm)
{
    if (sensitivity_mdeg_per_riu <= 0.0 || decay_length_nm <= 0.0) return 0.0;

    /* Convert angle shift to refractive index change */
    double delta_n = angle_shift_mdeg / sensitivity_mdeg_per_riu;

    /* delta_n = (dn/dc) · Γ where Γ is surface concentration [g/cm²]
     * Γ = delta_n / (dn/dc) */
    if (dn_dc <= 0.0) dn_dc = 0.182;  /* typical for proteins */
    double surface_conc = delta_n / dn_dc;  /* g/cm² → convert to g/m² */

    /* For thin films: delta_n ∝ d (linear regime)
     * For thick films: delta_n saturates as 1 - e^(-2d/l_d)
     *
     * In the thin-film approximation (d << l_d):
     *   d = delta_n · l_d / (2 · (n_adlayer - n_buffer))
     *
     * We estimate using a typical RI difference of 0.15 for protein vs buffer. */
    double dn_adlayer_buffer = 0.15;  /* n_protein ≈ 1.48, n_buffer ≈ 1.33 */
    double d_thin = delta_n * decay_length_nm / (2.0 * dn_adlayer_buffer);

    /* For thicker films, solve d = (l_d/2) · ln(Δn_sat / (Δn_sat - Δn))
     * where Δn_sat is the maximum RI shift for an infinitely thick layer. */
    double dn_sat = dn_adlayer_buffer;
    if (delta_n < 0.95 * dn_sat) {
        double ratio = dn_sat / (dn_sat - delta_n);
        if (ratio > 1.0) {
            d_thin = (decay_length_nm / 2.0) * log(ratio);
        }
    }

    return (d_thin > 0.0) ? d_thin : 0.0;
}

/* ============================================================================
 * L2/L4: Fluorescence
 * ============================================================================ */

/**
 * Fluorescence intensity: I_f = I₀ · Φ · (1 - 10^(-A)) · f_geo
 */
double fluorescence_intensity(const FluorescenceParams *fp,
                              double excitation_intensity,
                              double concentration_m,
                              double path_length_cm,
                              double geometric_factor)
{
    if (!fp || excitation_intensity < 0.0 || concentration_m < 0.0) return 0.0;

    double absorbance = fp->extinction_coefficient * path_length_cm * concentration_m;
    double absorption_fraction = 1.0 - pow(10.0, -absorbance);

    /* For very dilute solutions: 1 - 10^(-A) ≈ 2.303 · A as A → 0 */
    if (absorbance < 0.01) {
        absorption_fraction = 2.302585 * absorbance;  /* ln(10) · A */
    }

    double geo = (geometric_factor >= 0.0 && geometric_factor <= 1.0)
                 ? geometric_factor : 1.0;

    return excitation_intensity * fp->quantum_yield * absorption_fraction * geo;
}

/**
 * FRET efficiency: E = 1 / (1 + (r/R₀)⁶)
 *
 * The strong distance dependence (1/r⁶) makes FRET exquisitely
 * sensitive to nanoscale conformational changes, ideal for
 * biosensing molecular interactions.
 */
double fret_efficiency(double forster_radius_nm, double separation_nm)
{
    if (forster_radius_nm <= 0.0 || separation_nm <= 0.0) return 0.0;

    double ratio = separation_nm / forster_radius_nm;
    double ratio_pow6 = ratio * ratio * ratio;  /* r³/R₀³ */
    ratio_pow6 = ratio_pow6 * ratio_pow6;        /* r⁶/R₀⁶ */

    return 1.0 / (1.0 + ratio_pow6);
}

/**
 * FRET distance: r = R₀ · (1/E - 1)^(1/6)
 */
double fret_distance(double forster_radius_nm, double measured_efficiency)
{
    if (forster_radius_nm <= 0.0 || measured_efficiency <= 0.0) return INFINITY;
    if (measured_efficiency >= 1.0) return 0.0;

    double term = 1.0 / measured_efficiency - 1.0;
    return forster_radius_nm * pow(term, 1.0 / 6.0);
}

/* ============================================================================
 * L6: Fiber-optic evanescent wave biosensor
 * ============================================================================ */

/**
 * Evanescent wave absorbance in fiber-optic biosensor.
 *
 * The evanescent field extends beyond the fiber core into the
 * sensing medium. Only the fraction η of optical power in the
 * evanescent tail interacts with the analyte.
 *
 * A_eff = η · ε · c · L_eff
 * where η = V²/(V²+1) for weakly-guiding approximation.
 */
double fiber_evanescent_absorbance(const FiberOpticBiosensorParams *fiber_params,
                                   double extinction_coeff,
                                   double concentration_m)
{
    if (!fiber_params || extinction_coeff <= 0.0 || concentration_m < 0.0) return 0.0;

    /* Evanescent power fraction depends on V-parameter */
    double V = (2.0 * M_PI * fiber_params->fiber_core_diameter_um * 1.0e-6 / 2.0) /
               (fiber_params->evanescent_decay_length_nm * 1.0e-9) *
               fiber_params->numerical_aperture;

    /* Fraction of power in cladding (evanescent) for step-index fiber */
    double eta;
    if (V < 1.0) {
        eta = 0.5;  /* most power in cladding for V << 1 */
    } else if (V > 10.0) {
        eta = 4.0 / (3.0 * sqrt(V));  /* weakly guided, most power in core */
    } else {
        eta = 1.0 - (2.405 / V) * (2.405 / V);  /* approximation */
    }

    /* Effective path length: sensing region length × evanescent fraction */
    double effective_path = fiber_params->sensing_region_length_mm * 1.0e-1 * eta; /* cm */
    double absorbance = extinction_coeff * effective_path * concentration_m;

    /* Account for propagation loss */
    double loss_factor = pow(10.0,
        -fiber_params->propagation_loss_db_per_cm *
         fiber_params->sensing_region_length_mm * 1.0e-1 / 10.0);

    return absorbance * fiber_params->coupling_efficiency * loss_factor;
}

/**
 * Evanescent power fraction from V-parameter.
 *
 * V = (2π·a/λ) · √(n_core² - n_clad²) = (2π·a/λ) · NA
 *
 * For single-mode fiber (V < 2.405), the fraction of power
 * propagating in the cladding (evanescent tail) is approximately:
 *   η ≈ 4 / (3 · √V)  for V > 2.4
 */
double evanescent_power_fraction(const FiberOpticBiosensorParams *fiber_params,
                                 double wavelength_nm)
{
    if (!fiber_params || wavelength_nm <= 0.0) return 0.0;

    double a = fiber_params->fiber_core_diameter_um * 1.0e-6 / 2.0;  /* core radius in m */
    double V = (2.0 * M_PI * a) / (wavelength_nm * 1.0e-9) *
               fiber_params->numerical_aperture;

    if (V <= 0.0) return 0.0;

    if (V < 2.4) {
        /* Single-mode regime: use Marcuse approximation */
        double w = 1.1428 * V - 0.996;
        double gamma = sqrt(w * w);
        return 1.0 / (V * V * log(V));
    }

    /* Multi-mode: power fraction in evanescent field of highest-order mode */
    return 4.0 / (3.0 * sqrt(V));
}

/* ============================================================================
 * L6: ELISA — 4PL model
 * ============================================================================ */

/**
 * 4PL prediction: signal = a + (d - a) / (1 + (conc/c)^b)
 */
double elisa_fourpl_predict(const FourPLLogisticFit *fit, double concentration)
{
    if (!fit) return 0.0;
    if (concentration <= 0.0) return fit->a;

    double ratio = concentration / fit->c;
    if (ratio <= 0.0) return fit->a;

    double ratio_pow = pow(ratio, fit->b);
    return fit->a + (fit->d - fit->a) / (1.0 + ratio_pow);
}

/**
 * Inverse 4PL: concentration = c · ((d - a)/(signal - a) - 1)^(1/b)
 */
double elisa_fourpl_inverse(const FourPLLogisticFit *fit, double signal)
{
    if (!fit) return NAN;
    if (fit->b == 0.0) return NAN;

    double range = fit->d - fit->a;
    if (fabs(range) < 1.0e-15) return NAN;

    if (signal <= fit->a) return 0.0;
    if (signal >= fit->d) return INFINITY;

    double frac = range / (signal - fit->a) - 1.0;
    if (frac <= 0.0) return INFINITY;

    return fit->c * pow(frac, 1.0 / fit->b);
}

/**
 * Weighted mean concentration from replicate ELISA wells.
 *
 * Inverse-variance weighting gives more weight to precise measurements.
 * Variance increases with distance from EC50 (heteroscedasticity).
 */
double elisa_weighted_mean(const double *concentrations,
                           const double *variances,
                           int n_replicates)
{
    if (!concentrations || !variances || n_replicates <= 0) return 0.0;

    double sum_weighted = 0.0;
    double sum_weights  = 0.0;

    for (int i = 0; i < n_replicates; i++) {
        if (variances[i] <= 0.0) continue;  /* skip zero-variance (invalid) */
        double weight = 1.0 / variances[i];
        sum_weighted += weight * concentrations[i];
        sum_weights  += weight;
    }

    return (sum_weights > 0.0) ? (sum_weighted / sum_weights) : 0.0;
}

/* ============================================================================
 * L7: DNA microarray hybridization
 * ============================================================================ */

/**
 * DNA hybridization efficiency via Langmuir isotherm.
 *
 * At equilibrium, the fraction of surface probes hybridized is:
 *   θ = K · c / (1 + K · c)
 *
 * K = equilibrium constant, related to ΔG via K = exp(-ΔG/RT).
 * Typical K for 25-mer DNA at 45°C: 10⁶ to 10⁸ M⁻¹.
 */
double dna_hybridization_efficiency(double equilibrium_constant,
                                    double target_conc_m)
{
    if (equilibrium_constant <= 0.0 || target_conc_m < 0.0) return 0.0;

    double kc = equilibrium_constant * target_conc_m;
    return kc / (1.0 + kc);
}

/**
 * Wallace rule for short oligonucleotide melting temperature.
 *
 * T_m = 2°C × (A + T) + 4°C × (G + C)
 *
 * Valid for oligonucleotides with 4-14 bases at standard salt
 * (≈ 50 mM Na⁺). For longer probes or different salt conditions,
 * use the SantaLucia nearest-neighbor model.
 *
 * Formula is in °C. The double-stranded DNA denatures (melts)
 * at T_m where 50% of the strands are dissociated.
 */
double dna_probe_melting_temperature(const char *sequence)
{
    if (!sequence) return 0.0;

    int at_count = 0;
    int gc_count = 0;

    for (const char *p = sequence; *p != '\0'; p++) {
        switch (*p) {
            case 'A': case 'a': at_count++; break;
            case 'T': case 't': at_count++; break;
            case 'G': case 'g': gc_count++; break;
            case 'C': case 'c': gc_count++; break;
            default: break;  /* ignore non-standard characters */
        }
    }

    return 2.0 * at_count + 4.0 * gc_count;
}
