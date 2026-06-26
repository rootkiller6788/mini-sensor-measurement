/**
 * biosensor_electrochemical.h — Electrochemical transduction mechanisms
 *
 * L2 Core Concepts: Amperometry, potentiometry, impedimetry, voltammetry.
 * L4 Fundamental Laws: Nernst equation, Butler-Volmer kinetics, Cottrell equation.
 * L6 Canonical Problems: Glucose biosensor amperometric response, pH-ISFET calibration.
 *
 * References:
 *   - Bard & Faulkner, "Electrochemical Methods" (2001)
 *   - Wang, "Analytical Electrochemistry" (2006)
 *   - Thévenot et al., Biosensors & Bioelectronics 16 (2001) 121-131
 *
 * L9 alignment: MIT 6.555 (Bioelectronics), UC Berkeley BioE 121
 */

#ifndef BIOSENSOR_ELECTROCHEMICAL_H
#define BIOSENSOR_ELECTROCHEMICAL_H

#include "biosensor_types.h"

/* ============================================================================
 * L2: Amperometric biosensor — fixed-potential current measurement
 * ============================================================================ */

/**
 * Compute the Nernst equilibrium potential.
 *
 * Theorem (Nernst, 1889): E = E° + (RT / nF) * ln([Ox]/[Red])
 *   Where E° is the standard reduction potential, R is the gas constant,
 *   T is temperature, n is electrons transferred, F is Faraday constant.
 *
 * This is L4 — a fundamental law of electrochemistry.
 *
 * @param np     Initialized NernstParams (E°, n, T, and constants)
 * @param conc_oxidized  Concentration of oxidized species [M]
 * @param conc_reduced   Concentration of reduced species [M]
 * @return       Half-cell potential [V]
 */
double nernst_potential(const NernstParams *np,
                        double conc_oxidized,
                        double conc_reduced);

/**
 * Compute pH from ISFET threshold voltage shift (L4: Nernst-derived).
 *
 * The ISFET uses a pH-sensitive gate insulator (Si₃N₄, Al₂O₃, Ta₂O₅).
 * ΔV_th = sensitivity * (pH - pH_ref)
 *
 * @param isfet    ISFET parameters with calibrated sensitivity
 * @param v_th_measured  Measured threshold voltage [V]
 * @return         Calculated pH value
 */
double isfet_compute_ph(const ISFETParams *isfet, double v_th_measured);

/**
 * Compute the Butler-Volmer current density (L4).
 *
 * i = i₀ * [exp(α_a * n * F * η / (R * T)) - exp(-α_c * n * F * η / (R * T))]
 *
 * Describes the relationship between overpotential η and net current density
 * at an electrode-electrolyte interface.
 *
 * @param bv     Butler-Volmer parameters (i₀, α_a, α_c, area)
 * @return       Net current [A]
 */
double butler_volmer_current(const ButlerVolmerParams *bv);

/**
 * Compute the Cottrell current decay (L4: diffusion-limited amperometry).
 *
 * I(t) = (n * F * A * √D * C) / √(π * t)
 *
 * This equation governs the transient current response after a potential step,
 * limited by semi-infinite linear diffusion of analyte to the electrode.
 *
 * @param n_electrons    Electrons per redox event
 * @param electrode_area Electrode surface area [cm²]
 * @param diffusion_coeff Diffusion coefficient D [cm²/s]
 * @param bulk_conc      Bulk analyte concentration [mol/cm³]
 * @param time_sec       Time after potential step [s]
 * @return               Diffusion-limited current [A]
 */
double cottrell_current(int n_electrons,
                        double electrode_area,
                        double diffusion_coeff,
                        double bulk_conc,
                        double time_sec);

/**
 * Randles-Ševčík equation for peak current in cyclic voltammetry (L4).
 *
 * I_p = 2.69e5 * n^(3/2) * A * √D * C * √v   (at 25°C, reversible)
 *
 * @param n_electrons    Electrons transferred
 * @param electrode_area Electrode area [cm²]
 * @param diffusion_coeff D [cm²/s]
 * @param bulk_conc      Bulk concentration [mol/cm³]
 * @param scan_rate      Scan rate v [V/s]
 * @return               Peak current [A]
 */
double randles_sevcik_peak_current(int n_electrons,
                                   double electrode_area,
                                   double diffusion_coeff,
                                   double bulk_conc,
                                   double scan_rate);

/* ============================================================================
 * L2: Glucose biosensor (canonical amperometric sensor)
 * ============================================================================ */

/**
 * Model glucose oxidase (GOx) based amperometric glucose biosensor.
 *
 * Reaction: Glucose + O₂ →(GOx)→ Gluconic acid + H₂O₂
 * Detection: H₂O₂ → O₂ + 2H⁺ + 2e⁻  (at +0.6 to +0.7 V vs Ag/AgCl)
 *
 * The current is proportional to glucose concentration within the
 * linear range, limited by oxygen co-substrate availability and
 * enzyme saturation (Michaelis-Menten).
 *
 * L2 + L6: Core concept + canonical problem.
 *
 * @param glucose_conc_mM   Glucose concentration [mM]
 * @param enzyme_activity_u  Enzyme activity on electrode [U/cm²]
 * @param electrode_area_cm2 Electrode area [cm²]
 * @param oxygen_conc_mM     Dissolved O₂ concentration [mM]
 * @return                   Expected biosensor current [nA]
 */
double glucose_biosensor_current(double glucose_conc_mM,
                                 double enzyme_activity_u,
                                 double electrode_area_cm2,
                                 double oxygen_conc_mM);

/**
 * Estimate glucose concentration from measured amperometric current (inverse model).
 *
 * Uses linear + Michaelis-Menten hybrid model for accuracy across range.
 *
 * @param current_na        Measured current [nA]
 * @param sensitivity_na_per_mM Sensor sensitivity [nA/mM]
 * @param km_apparent_mM    Apparent Km of immobilized GOx [mM]
 * @param imax_na           Maximum current at enzyme saturation [nA]
 * @return                  Estimated glucose concentration [mM]
 */
double glucose_conc_from_current(double current_na,
                                 double sensitivity_na_per_mM,
                                 double km_apparent_mM,
                                 double imax_na);

/* ============================================================================
 * L2: Electrochemical Impedance Spectroscopy (EIS) analysis
 * ============================================================================ */

/**
 * Compute the complex impedance of a Randles equivalent circuit.
 *
 * Z(ω) = R_s + [R_ct + σ/√(jω)] || C_dl
 *
 * Returns both real (Z') and imaginary (Z″) components.
 *
 * L2 Core Concept: Impedimetric biosensing — label-free detection
 * of binding events via interfacial impedance change.
 *
 * @param eis      EIS parameters (R_s, R_ct, C_dl, σ)
 * @param freq_hz  Measurement frequency [Hz]
 * @param z_real   [out] Real part Z' [Ω]
 * @param z_imag   [out] Imaginary part Z″ [Ω]
 */
void eis_randles_impedance(const EISParams *eis,
                           double freq_hz,
                           double *z_real,
                           double *z_imag);

/**
 * Fit R_ct from EIS Nyquist data using the Randles model.
 *
 * Under charge-transfer control (high frequency), the Nyquist plot is
 * a semicircle with diameter = R_ct. This function estimates R_ct from
 * the high-frequency arc via circle fitting.
 *
 * L5 Algorithm: Circle-fit method for Nyquist plot analysis.
 *
 * @param freq_array    Frequency array [Hz]
 * @param z_real_array  Z' array [Ω]
 * @param z_imag_array  Z″ array [Ω] (negative values as measured)
 * @param n_points      Number of frequency points
 * @param rs_est        [out] Estimated R_s from high-f intercept
 * @param rct_est       [out] Estimated R_ct from semicircle diameter
 * @return              RMS fit error of the circle model
 */
double eis_fit_charge_transfer_resistance(const double *freq_array,
                                          const double *z_real_array,
                                          const double *z_imag_array,
                                          int n_points,
                                          double *rs_est,
                                          double *rct_est);

/* ============================================================================
 * L6: Ion-Selective Electrode (ISE) — potentiometric biosensor
 * ============================================================================ */

/**
 * Nikolskii-Eisenman equation for ion-selective electrode (L4).
 *
 * E = E° + (RT / z_i F) * ln( a_i + Σ K_ij * a_j^(z_i/z_j) )
 *
 * Where K_ij is the selectivity coefficient for interfering ion j.
 *
 * @param np              Nernst equation parameters
 * @param primary_activity Activity of primary ion [mol/L]
 * @param interferent_activities Array of interferent ion activities
 * @param selectivity_coeffs      Array of selectivity coefficients K_ij
 * @param interferent_charges     Array of interferent ion charges
 * @param n_interferents   Number of interfering ions
 * @return                 Electrode potential [V]
 */
double nikolskii_eisenman_potential(const NernstParams *np,
                                    double primary_activity,
                                    const double *interferent_activities,
                                    const double *selectivity_coeffs,
                                    const int *interferent_charges,
                                    int n_interferents);

/**
 * Compute ion activity from Debye-Hückel theory (L3).
 *
 * log₁₀(γ) = -A * z² * √I / (1 + B * a * √I)
 *
 * Extended Debye-Hückel equation for ionic strength correction.
 *
 * @param ion_charge       z: charge number of ion
 * @param ionic_strength   I = 0.5 * Σ c_i * z_i² [mol/L]
 * @param ion_size_param   a: ion size parameter [Å], typical 3-9 Å
 * @return                 Activity coefficient γ (0 < γ ≤ 1)
 */
double debye_huckel_activity_coeff(int ion_charge,
                                   double ionic_strength,
                                   double ion_size_param);

/* ============================================================================
 * L7: Point-of-care diagnostic — lactate biosensor
 * ============================================================================ */

/**
 * Lactate biosensor model (L7 Application: sports medicine / critical care).
 *
 * Lactate + O₂ →(Lactate Oxidase)→ Pyruvate + H₂O₂
 *
 * Used in point-of-care devices (e.g., i-STAT, Lactate Scout) for
 * assessing tissue oxygenation, sepsis, and athletic performance.
 *
 * @param lactate_conc_mM   Blood lactate concentration [mM]
 * @param temp_celsius      Measurement temperature [°C]
 * @param hematocrit_percent Blood hematocrit [%] for correction
 * @return                  Sensor current [nA]
 */
double lactate_biosensor_current(double lactate_conc_mM,
                                 double temp_celsius,
                                 double hematocrit_percent);

#endif /* BIOSENSOR_ELECTROCHEMICAL_H */
