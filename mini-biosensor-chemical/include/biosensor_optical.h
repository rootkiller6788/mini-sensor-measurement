/**
 * biosensor_optical.h — Optical biosensor transduction mechanisms
 *
 * L2 Core Concepts: Surface Plasmon Resonance (SPR), fluorescence,
 *   absorbance, bioluminescence, evanescent wave sensing.
 * L4 Fundamental Laws: Beer-Lambert law, Fresnel equations, Förster resonance
 *   energy transfer (FRET), Snell's law.
 * L6 Canonical Problems: SPR immunoassay, DNA microarray fluorescence,
 *   fiber-optic glucose sensor.
 *
 * References:
 *   - Lakowicz, "Principles of Fluorescence Spectroscopy" (2006)
 *   - Homola, "Surface Plasmon Resonance Based Sensors" (2006)
 *   - Wolfbeis, "Fiber-Optic Chemical Sensors and Biosensors" (2008)
 *
 * L9 alignment: Stanford EE 247 (Optical), MIT 6.555, Cambridge Part II Biosensors
 */

#ifndef BIOSENSOR_OPTICAL_H
#define BIOSENSOR_OPTICAL_H

#include "biosensor_types.h"

/* ============================================================================
 * L4: Beer-Lambert Law — absorbance-based sensing
 * ============================================================================ */

/**
 * Compute absorbance from concentration (Beer-Lambert Law, L4).
 *
 * A = ε · b · c
 *
 * Where: A = absorbance (optical density, dimensionless)
 *        ε = molar extinction coefficient [M⁻¹·cm⁻¹]
 *        b = path length [cm]
 *        c = concentration [M]
 *
 * This is the foundational law for all spectrophotometric biosensors
 * including ELISA readers, oximeters, and colorimetric strips.
 *
 * @param extinction_coeff  ε [M⁻¹·cm⁻¹]
 * @param path_length_cm    b [cm]
 * @param concentration_m   c [M]
 * @return                  Absorbance A (dimensionless)
 */
double beer_lambert_absorbance(double extinction_coeff,
                               double path_length_cm,
                               double concentration_m);

/**
 * Inverse Beer-Lambert: concentration from measured absorbance.
 *
 * @param absorbance         Measured absorbance
 * @param extinction_coeff   ε [M⁻¹·cm⁻¹]
 * @param path_length_cm     b [cm]
 * @return                   Concentration [M]
 */
double concentration_from_absorbance(double absorbance,
                                     double extinction_coeff,
                                     double path_length_cm);

/**
 * Two-component absorbance analysis (simultaneous equations).
 *
 * For a mixture of two absorbing species with overlapping spectra,
 * solve the system:
 *   A_λ1 = ε_a1 · b · c_a + ε_b1 · b · c_b
 *   A_λ2 = ε_a2 · b · c_a + ε_b2 · b · c_b
 *
 * L6 Canonical Problem: Multicomponent spectrophotometric analysis.
 *
 * @param abs_lambda1    Absorbance at wavelength 1
 * @param abs_lambda2    Absorbance at wavelength 2
 * @param eps_a1,eps_a2  Extinction coefficients of species A at λ1,λ2
 * @param eps_b1,eps_b2  Extinction coefficients of species B at λ1,λ2
 * @param path_length    Path length [cm]
 * @param conc_a         [out] Concentration of species A [M]
 * @param conc_b         [out] Concentration of species B [M]
 * @return               true if solution is well-conditioned, false if singular
 */
bool two_component_absorbance(double abs_lambda1, double abs_lambda2,
                              double eps_a1, double eps_a2,
                              double eps_b1, double eps_b2,
                              double path_length,
                              double *conc_a, double *conc_b);

/* ============================================================================
 * L2: Surface Plasmon Resonance (SPR) — label-free optical biosensing
 * ============================================================================ */

/**
 * Compute SPR resonance angle (L2 + L3: Fresnel equations).
 *
 * For the Kretschmann configuration (prism / Au / sample):
 *   k_SP = (2π / λ) * √(ε_m * ε_d / (ε_m + ε_d))
 *   θ_SPR = arcsin( k_SP / (n_prism * 2π/λ) )
 *
 * Where k_SP is the surface plasmon wavevector, ε_m is the metal
 * permittivity (Au at λ), ε_d is the dielectric constant of the
 * sensing medium (sample).
 *
 * @param spr              SPR system parameters
 * @param sample_refractive_index RI of sample above gold
 * @param gold_permittivity_real Real part of gold permittivity at λ
 * @return                  SPR resonance angle [degrees]
 */
double spr_resonance_angle(const SPRParams *spr,
                           double sample_refractive_index,
                           double gold_permittivity_real);

/**
 * Compute SPR reflectivity curve (Fresnel multi-layer model).
 *
 * R(θ) for a 3-layer system (prism/Au/sample). Used to simulate
 * the full SPR dip for quantitative analysis of binding kinetics.
 *
 * L3 Mathematical Structure: Transfer matrix method for multi-layer optics.
 *
 * @param spr               SPR system parameters
 * @param angle_deg         Incident angle [degrees]
 * @param sample_ri         Sample refractive index
 * @param gold_ri_real      Real part of gold RI
 * @param gold_ri_imag      Imaginary part of gold RI
 * @return                  Reflectivity R (0 to 1)
 */
double spr_reflectivity(const SPRParams *spr,
                        double angle_deg,
                        double sample_ri,
                        double gold_ri_real,
                        double gold_ri_imag);

/**
 * Estimate biomolecular layer thickness from SPR angle shift (L5).
 *
 * Δθ_SPR ∝ Δn · d · exp(-d / l_d)
 *
 * Where d is the adlayer thickness, Δn is the RI difference between
 * adlayer and buffer, l_d is the evanescent decay length.
 *
 * Used in label-free kinetic analysis of protein-protein interactions.
 *
 * @param angle_shift_mdeg  Measured SPR angle shift [millidegrees]
 * @param dn_dc             Refractive index increment of analyte [mL/g]
 * @param molecular_weight  Analyte molecular weight [Da]
 * @param sensitivity_mdeg_per_riu Instrument sensitivity [mdeg/RIU]
 * @param decay_length_nm   Evanescent field decay length [nm]
 * @return                  Estimated adlayer thickness [nm]
 */
double spr_adlayer_thickness(double angle_shift_mdeg,
                             double dn_dc,
                             double molecular_weight,
                             double sensitivity_mdeg_per_riu,
                             double decay_length_nm);

/* ============================================================================
 * L2: Fluorescence-based biosensing
 * ============================================================================ */

/**
 * Compute fluorescence intensity (L4: Fundamental relationship).
 *
 * I_f = I₀ · Φ · (1 - 10^(-A)) · f_geo
 *
 * Where I₀ is excitation intensity, Φ is quantum yield,
 * A is absorbance at excitation wavelength, f_geo is geometric factor.
 *
 * For dilute solutions (A < 0.05): I_f ≈ 2.303 · I₀ · Φ · ε · b · c · f_geo
 *
 * @param fp              Fluorescence parameters
 * @param excitation_intensity  I₀ [W/cm² or arbitrary units]
 * @param concentration_m Analyte concentration [M]
 * @param path_length_cm  Optical path length [cm]
 * @param geometric_factor Detection geometry factor (0 to 1)
 * @return                Fluorescence intensity [a.u.]
 */
double fluorescence_intensity(const FluorescenceParams *fp,
                              double excitation_intensity,
                              double concentration_m,
                              double path_length_cm,
                              double geometric_factor);

/**
 * Compute FRET efficiency (Förster Resonance Energy Transfer).
 *
 * E = R₀⁶ / (R₀⁶ + r⁶) = 1 / (1 + (r/R₀)⁶)
 *
 * Where R₀ is the Förster radius (50% efficiency distance),
 * r is the donor-acceptor separation distance.
 *
 * L4 + L6: FRET is the dominant mechanism for homogeneous immunoassays
 * and DNA hybridization detection in biosensors.
 *
 * @param forster_radius_nm  R₀ [nm] — distance at 50% efficiency
 * @param separation_nm      r [nm] — donor-acceptor distance
 * @return                   FRET efficiency E (0 to 1)
 */
double fret_efficiency(double forster_radius_nm, double separation_nm);

/**
 * Compute donor-acceptor distance from measured FRET efficiency.
 *
 * r = R₀ · (1/E - 1)^(1/6)
 *
 * Used in molecular ruler applications — measuring nanoscale distances
 * in protein conformational changes, DNA hybridization, and
 * aptamer-target binding events.
 *
 * @param forster_radius_nm  R₀ [nm]
 * @param measured_efficiency E (0 to 1)
 * @return                   Estimated distance r [nm]
 */
double fret_distance(double forster_radius_nm, double measured_efficiency);

/* ============================================================================
 * L6: Fiber-optic biosensor — evanescent wave sensing
 * ============================================================================ */

/**
 * Compute evanescent wave absorbance in a fiber-optic biosensor.
 *
 * When the fiber cladding is removed and replaced by an absorbing
 * sensing medium, the evanescent field interacts with the analyte:
 *
 * A_ev = η · ε · c · L · (d_p / d_core)
 *
 * Where η is the fraction of power in the evanescent field,
 * d_p is the penetration depth, d_core is the core diameter.
 *
 * @param fiber_params        Fiber-optic biosensor parameters
 * @param extinction_coeff    Analyte ε [M⁻¹·cm⁻¹]
 * @param concentration_m     Analyte concentration [M]
 * @return                    Effective evanescent absorbance
 */
double fiber_evanescent_absorbance(const FiberOpticBiosensorParams *fiber_params,
                                   double extinction_coeff,
                                   double concentration_m);

/**
 * Compute the fraction of optical power in the evanescent field.
 *
 * For a step-index fiber with core RI n_co, cladding RI n_cl,
 * at wavelength λ, the V-parameter determines the mode distribution:
 *
 * V = (2π · a / λ) · √(n_co² - n_cl²)
 *
 * η ≈ 4 / (3 · √V) for weakly guiding fibers (V > 2.4, single-mode regime).
 *
 * @param fiber_params    Fiber biosensor parameters
 * @param wavelength_nm   Operating wavelength [nm]
 * @return                Evanescent power fraction η (0 to 1)
 */
double evanescent_power_fraction(const FiberOpticBiosensorParams *fiber_params,
                                 double wavelength_nm);

/* ============================================================================
 * L6: ELISA — Enzyme-Linked Immunosorbent Assay (optical canonical)
 * ============================================================================ */

/**
 * Model the ELISA calibration curve (4PL standard).
 *
 * The 4-parameter logistic model is the gold standard for immunoassay
 * calibration:
 *
 * signal = a + (d - a) / (1 + (conc / c)^b)
 *
 * Where a = background, d = saturation, c = EC50, b = slope.
 *
 * L6 Canonical Problem: ELISA quantification is the most widely used
 * biosensor format in clinical diagnostics.
 *
 * @param fit             Fitted 4PL parameters
 * @param concentration   Analyte concentration
 * @return                Predicted assay signal [OD, RFU, etc.]
 */
double elisa_fourpl_predict(const FourPLLogisticFit *fit, double concentration);

/**
 * Inverse 4PL: concentration from ELISA signal.
 *
 * conc = c · ((d - a) / (signal - a) - 1)^(1/b)
 *
 * @param fit             Fitted 4PL parameters
 * @param signal          Measured signal
 * @return                Estimated concentration (NaN if outside model range)
 */
double elisa_fourpl_inverse(const FourPLLogisticFit *fit, double signal);

/**
 * Compute the weighted mean concentration from replicate ELISA wells.
 *
 * Uses inverse-variance weighting:
 *   c̄ = Σ(c_i / σ_i²) / Σ(1 / σ_i²)
 *
 * For L5: Statistical methods for biosensor data reduction.
 *
 * @param concentrations  Array of concentration estimates
 * @param variances       Array of variances (σ²) for each estimate
 * @param n_replicates    Number of replicate measurements
 * @return                Weighted mean concentration
 */
double elisa_weighted_mean(const double *concentrations,
                           const double *variances,
                           int n_replicates);

/* ============================================================================
 * L7: Application — DNA microarray fluorescence analysis
 * ============================================================================ */

/**
 * Compute hybridization efficiency for DNA microarray spot.
 *
 * The Langmuir isotherm governs DNA hybridization at microarray surfaces:
 *
 * θ = K · c / (1 + K · c)
 *
 * Where θ is the fraction of probes hybridized (0 to 1),
 * K is the equilibrium constant, c is the target concentration.
 *
 * L7 Application: High-throughput gene expression profiling.
 *
 * @param equilibrium_constant  K [M⁻¹] — binding affinity
 * @param target_conc_m         Target DNA concentration [M]
 * @return                      Fraction of probes hybridized θ (0 to 1)
 */
double dna_hybridization_efficiency(double equilibrium_constant,
                                    double target_conc_m);

/**
 * Melting temperature for short DNA probes (Wallace rule, L4).
 *
 * T_m = 2°C × (A+T count) + 4°C × (G+C count)
 *
 * Valid for oligonucleotides ≤ 14 bases. For longer probes,
 * the SantaLucia nearest-neighbor model should be used.
 *
 * Fundamental to DNA biosensor design.
 *
 * @param sequence     DNA probe sequence (null-terminated, uppercase ACGT)
 * @return             Melting temperature T_m [°C]
 */
double dna_probe_melting_temperature(const char *sequence);

#endif /* BIOSENSOR_OPTICAL_H */
