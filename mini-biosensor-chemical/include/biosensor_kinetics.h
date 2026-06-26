/**
 * biosensor_kinetics.h — Biorecognition kinetics and thermodynamic models
 *
 * L3 Mathematical Structures: Michaelis-Menten, Hill, Langmuir, competitive
 *   binding, diffusion-reaction, mass transport.
 * L4 Fundamental Laws: Michaelis-Menten equation, Fick's laws,
 *   Hill equation, Langmuir isotherm, van't Hoff equation.
 * L5 Algorithms: Parameter estimation by nonlinear regression,
 *   Lineweaver-Burk, Eadie-Hofstee, Scatchard plots.
 * L6 Canonical Problems: Enzyme inhibition analysis,
 *   antigen-antibody binding kinetics, DNA hybridization thermodynamics.
 *
 * References:
 *   - Michaelis & Menten, Biochemische Zeitschrift 49 (1913) 333-369
 *   - Briggs & Haldane, Biochem. J. 19 (1925) 338-339
 *   - Hill, J. Physiol. 40 (1910) iv-vii
 *   - Langmuir, JACS 40 (1918) 1361-1403
 *   - Fick, Annalen der Physik 170 (1855) 59-86
 *
 * L9 alignment: MIT 20.320 (Biomolecular Kinetics), Stanford BioE 301
 */

#ifndef BIOSENSOR_KINETICS_H
#define BIOSENSOR_KINETICS_H

#include "biosensor_types.h"

/* ============================================================================
 * L4: Michaelis-Menten enzyme kinetics
 * ============================================================================ */

/**
 * Compute the reaction rate v for a given substrate concentration (L4).
 *
 * v = V_max · [S] / (K_m + [S])
 *
 * This is the foundational equation of enzyme kinetics. It describes
 * the rate of an enzymatic reaction as a function of substrate
 * concentration, assuming steady-state conditions (Briggs-Haldane).
 *
 * The assumptions are:
 *   1. [S] >> [E] (excess substrate)
 *   2. Steady state: d[ES]/dt = 0
 *   3. Initial velocity measurement (negligible product inhibition)
 *   4. Single-substrate, single-intermediate mechanism
 *
 * @param params     Michaelis-Menten parameters (K_m, V_max)
 * @return           Reaction rate v [same units as V_max]
 */
double michaelis_menten_rate(const MichaelisMentenParams *params);

/**
 * Compute substrate concentration from measured rate (inverse MM).
 *
 * [S] = K_m · v / (V_max - v)
 *
 * @param params     Michaelis-Menten parameters
 * @param rate       Measured reaction rate v
 * @return           Substrate concentration [same units as K_m]
 */
double michaelis_menten_inverse(const MichaelisMentenParams *params,
                                double rate);

/**
 * Competitive inhibition rate (L6: canonical inhibition problem).
 *
 * v = V_max · [S] / (K_m · (1 + [I]/K_i) + [S])
 *
 * The inhibitor competes with substrate for the active site.
 * Apparent K_m increases by factor (1 + [I]/K_i), V_max unchanged.
 *
 * @param params          Michaelis-Menten parameters
 * @param inhibitor_conc  Inhibitor concentration [I]
 * @param inhibition_constant  K_i (dissociation constant of EI complex)
 * @return                Inhibited reaction rate
 */
double competitive_inhibition_rate(const MichaelisMentenParams *params,
                                   double inhibitor_conc,
                                   double inhibition_constant);

/**
 * Uncompetitive inhibition rate.
 *
 * v = V_max · [S] / (K_m + [S] · (1 + [I]/K_i))
 *
 * Inhibitor binds only to ES complex, not free enzyme.
 * Both apparent K_m and V_max decrease by factor (1 + [I]/K_i).
 *
 * @param params          Michaelis-Menten parameters
 * @param inhibitor_conc  [I]
 * @param inhibition_constant  K_i
 * @return                Inhibited reaction rate
 */
double uncompetitive_inhibition_rate(const MichaelisMentenParams *params,
                                     double inhibitor_conc,
                                     double inhibition_constant);

/**
 * Non-competitive inhibition rate.
 *
 * v = (V_max / (1 + [I]/K_i)) · [S] / (K_m + [S])
 *
 * Inhibitor binds to both E and ES with equal affinity.
 * V_max decreases, K_m unchanged.
 *
 * @param params          Michaelis-Menten parameters
 * @param inhibitor_conc  [I]
 * @param inhibition_constant  K_i
 * @return                Inhibited reaction rate
 */
double noncompetitive_inhibition_rate(const MichaelisMentenParams *params,
                                      double inhibitor_conc,
                                      double inhibition_constant);

/* ============================================================================
 * L5: Lineweaver-Burk linearization methods
 * ============================================================================ */

/**
 * Compute Lineweaver-Burk transformed coordinates (L5).
 *
 * The Lineweaver-Burk plot linearizes Michaelis-Menten kinetics:
 *
 *   1/v = (K_m / V_max) · (1/[S]) + 1/V_max
 *
 * y = slope · x + intercept
 *   where slope = K_m/V_max, intercept = 1/V_max
 *
 * This allows linear regression estimation of K_m and V_max,
 * though it is biased by data weighting at low [S].
 *
 * L5 Algorithm: Linearization for parameter estimation.
 *
 * @param params             Michaelis-Menten parameters
 * @param concentrations     Array of [S] values
 * @param rates              Array of measured v values
 * @param n_points           Number of data points
 * @param km_est             [out] Estimated K_m
 * @param vmax_est           [out] Estimated V_max
 * @return                   R² of the Lineweaver-Burk fit
 */
double lineweaver_burk_estimate(const double *concentrations,
                                const double *rates,
                                int n_points,
                                double *km_est,
                                double *vmax_est);

/**
 * Eadie-Hofstee linearization (L5).
 *
 * v = -K_m · (v/[S]) + V_max
 *
 * Less biased than Lineweaver-Burk because both axes contain v,
 * distributing errors more uniformly.
 *
 * @param concentrations     Array of [S] values
 * @param rates              Array of v values
 * @param n_points           Number of data points
 * @param km_est             [out] K_m = -slope
 * @param vmax_est           [out] V_max = intercept
 * @return                   R² of the Eadie-Hofstee fit
 */
double eadie_hofstee_estimate(const double *concentrations,
                              const double *rates,
                              int n_points,
                              double *km_est,
                              double *vmax_est);

/**
 * Hanes-Woolf linearization (L5).
 *
 * [S]/v = (1/V_max) · [S] + K_m/V_max
 *
 * Generally the most statistically robust of the three linearizations
 * because [S] is the independent variable and errors are symmetric in [S]/v.
 *
 * @param concentrations     Array of [S] values
 * @param rates              Array of v values
 * @param n_points           Number of data points
 * @param km_est             [out] K_m = intercept/slope
 * @param vmax_est           [out] V_max = 1/slope
 * @return                   R² of the Hanes-Woolf fit
 */
double hanes_woolf_estimate(const double *concentrations,
                            const double *rates,
                            int n_points,
                            double *km_est,
                            double *vmax_est);

/**
 * Direct nonlinear least-squares estimation of K_m and V_max (L5).
 *
 * Gauss-Newton iteration for minimizing Σ(v_i - V_max·[S_i]/(K_m+[S_i]))².
 * This is the gold standard — all linearization methods introduce bias.
 *
 * @param concentrations    [S] values
 * @param rates             Measured v values
 * @param n_points          Number of data points
 * @param km_init           Initial guess for K_m
 * @param vmax_init         Initial guess for V_max
 * @param km_est            [out] Optimal K_m estimate
 * @param vmax_est          [out] Optimal V_max estimate
 * @param max_iterations    Maximum Gauss-Newton iterations
 * @param tolerance         Convergence tolerance
 * @return                  Number of iterations used (0 if failed)
 */
int michaelis_menten_nls_fit(const double *concentrations,
                             const double *rates,
                             int n_points,
                             double km_init, double vmax_init,
                             double *km_est, double *vmax_est,
                             int max_iterations, double tolerance);

/* ============================================================================
 * L4: Hill equation — cooperative binding
 * ============================================================================ */

/**
 * Compute fractional saturation from the Hill equation (L4).
 *
 * θ = [L]^n / (K_d^n + [L]^n)
 *
 * Where θ is the fraction of binding sites occupied,
 * n is the Hill coefficient (n > 1 = positive cooperativity),
 * K_d is the dissociation constant (concentration at θ = 0.5).
 *
 * Alternative form: θ = 1 / (1 + (K_d/[L])^n)
 *
 * @param hp             Hill parameters (K_d equivalent to half_max_conc)
 * @param ligand_conc    Free ligand concentration [L]
 * @return               Fractional saturation θ (0 to 1)
 */
double hill_fractional_saturation(const HillParams *hp, double ligand_conc);

/**
 * Hill plot transformation for estimating cooperativity (L5).
 *
 * log(θ / (1 - θ)) = n · log[L] - n · log(K_d)
 *
 * The slope gives the Hill coefficient n, y-intercept gives -n·log(K_d).
 * This is the classic diagnostic for cooperative binding in hemoglobin,
 * transcription factors, and multi-subunit biosensors.
 *
 * L5 Algorithm: Graphical cooperativity analysis.
 *
 * @param concentrations   Array of ligand concentrations [L]
 * @param saturations      Array of measured θ values (0 to 1)
 * @param n_points         Number of data points
 * @param n_hill           [out] Estimated Hill coefficient
 * @param log_kd           [out] Estimated log₁₀(K_d)
 * @return                 R² of the Hill plot linear fit
 */
double hill_plot_estimate(const double *concentrations,
                          const double *saturations,
                          int n_points,
                          double *n_hill,
                          double *log_kd);

/* ============================================================================
 * L4: Langmuir adsorption isotherm + diffusion
 * ============================================================================ */

/**
 * Langmuir isotherm surface coverage (L4).
 *
 * θ = K · c / (1 + K · c)
 *
 * Where θ = q/q_max = fractional surface coverage,
 * K = equilibrium constant (k_ads/k_des),
 * c = bulk concentration.
 *
 * @param lp              Langmuir parameters
 * @param bulk_conc       Bulk concentration c
 * @return                Surface coverage θ (0 to 1)
 */
double langmuir_coverage(const LangmuirParams *lp, double bulk_conc);

/**
 * Langmuir kinetic model — time-dependent surface coverage (L6).
 *
 * dθ/dt = k_ads · c · (1 - θ) - k_des · θ
 *
 * Analytical solution:
 * θ(t) = θ_eq · (1 - exp(-k_obs · t))
 * where k_obs = k_ads · c + k_des, θ_eq = k_ads · c / k_obs
 *
 * This governs real-time SPR and QCM biosensor binding curves.
 *
 * @param lp              Langmuir parameters
 * @param bulk_conc       Constant bulk concentration c
 * @param time_sec        Time t [s]
 * @return                Surface coverage θ(t)
 */
double langmuir_kinetic_coverage(const LangmuirParams *lp,
                                 double bulk_conc,
                                 double time_sec);

/**
 * Langmuir dissociation model — pure desorption (L6: SPR wash step).
 *
 * θ(t) = θ₀ · exp(-k_des · t)   (after buffer wash, c ≈ 0)
 *
 * Used to extract k_des from the dissociation phase of SPR sensorgrams.
 *
 * @param lp              Langmuir parameters
 * @param initial_coverage θ₀ at dissociation start
 * @param time_sec        Time since start of dissociation [s]
 * @return                Remaining surface coverage θ(t)
 */
double langmuir_dissociation(const LangmuirParams *lp,
                             double initial_coverage,
                             double time_sec);

/**
 * One-dimensional Fickian diffusion model (semi-infinite, L4).
 *
 * C(x,t) = C₀ · erfc(x / √(4Dt))
 *
 * Where erfc is the complementary error function.
 * This describes analyte transport from bulk solution to the
 * sensor surface — the mass-transport limitation in biosensors.
 *
 * @param bulk_conc       Bulk concentration C₀ [M]
 * @param diffusion_coeff Diffusion coefficient D [cm²/s]
 * @param distance_cm     Distance from bulk x [cm]
 * @param time_sec        Time t [s]
 * @return                Local concentration C(x,t) [M]
 */
double fickian_concentration(double bulk_conc,
                             double diffusion_coeff,
                             double distance_cm,
                             double time_sec);

/**
 * Damköhler number — ratio of reaction rate to diffusion rate.
 *
 * Da = k_cat · Γ / (D · C₀ / δ)
 *
 * Where Γ is the surface enzyme density, δ is the diffusion layer thickness.
 * Da << 1: kinetically limited; Da >> 1: diffusion limited.
 *
 * L3: Dimensionless analysis for biosensor design optimization.
 *
 * @param kcat              Catalytic rate constant [1/s]
 * @param surface_enzyme_density Γ [mol/cm²]
 * @param diffusion_coeff   D [cm²/s]
 * @param bulk_conc         C₀ [mol/cm³]
 * @param diffusion_layer_thickness δ [cm]
 * @return                  Damköhler number Da
 */
double damkohler_number(double kcat,
                        double surface_enzyme_density,
                        double diffusion_coeff,
                        double bulk_conc,
                        double diffusion_layer_thickness);

/* ============================================================================
 * L7: Application — Immunoassay kinetics (antibody-antigen binding)
 * ============================================================================ */

/**
 * Antibody-antigen equilibrium binding.
 *
 * K_d = [Ab][Ag] / [Ab·Ag]
 *
 * The law of mass action applied to immunoassays:
 * Bound fraction = [Ag] / (K_d + [Ag])  (for [Ab]₀ >> [Ag·Ab])
 *
 * This is the theoretical basis for all immunoassay biosensors
 * (ELISA, SPR immunoassay, lateral flow tests).
 *
 * L7 Application: Clinical immunoassay design and QC.
 *
 * @param kd                 Dissociation constant [M], typically nM to pM
 * @param antigen_conc       Free antigen concentration [M]
 * @param antibody_conc      Total antibody binding sites [M]
 * @return                   Concentration of Ab·Ag complex [M]
 */
double antibody_antigen_complex(double kd,
                                double antigen_conc,
                                double antibody_conc);

/**
 * Binding kinetics curve (association phase).
 *
 * R(t) = R_eq · (1 - exp(-(k_on · C + k_off) · t))
 *
 * Used in SPR and BLI (Bio-Layer Interferometry) for measuring
 * k_on and k_off rate constants of biomolecular interactions.
 *
 * @param k_on              Association rate constant [M⁻¹s⁻¹]
 * @param k_off             Dissociation rate constant [s⁻¹]
 * @param analyte_conc      Analyte concentration C [M]
 * @param r_max             Maximum response at saturation
 * @param time_sec          Time [s]
 * @return                  Binding response R(t) [RU or nm shift]
 */
double binding_association_curve(double k_on, double k_off,
                                 double analyte_conc, double r_max,
                                 double time_sec);

#endif /* BIOSENSOR_KINETICS_H */
