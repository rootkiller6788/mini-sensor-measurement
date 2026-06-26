/**
 * @file mems_dynamics.h
 * @brief L2-L3 MEMS Mechanical Dynamics — Spring-mass-damper, transfer functions,
 *        Coriolis effect, structural analysis
 *
 * Covers:
 *   L2 - Spring-mass-damper systems, mechanical resonance, damping mechanisms
 *        (squeeze-film, Couette, thermoelastic), Coriolis effect in gyroscopes
 *   L3 - ODE solutions (under/critical/overdamped), transfer functions (Laplace),
 *        frequency response, Bode plots, Stoney's formula
 *
 * Reference: Bao, "Analysis and Design Principles of MEMS Devices" (2005)
 *            Senturia, "Microsystem Design" (2001)
 *            Thomson, "Theory of Vibration with Applications" (1996)
 * Courses: MIT 6.630 (EM Waves & Sensors), Stanford EE247 (MEMS & Sensors),
 *          ETH 227-0455, Tsinghua 传感器原理
 */

#ifndef MEMS_DYNAMICS_H
#define MEMS_DYNAMICS_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif

/* ─── L2: Spring-Mass-Damper State ───────────────────────────────────────────
 * Second-order mechanical system: m·ẍ + c·ẋ + k·x = F_ext(t)
 */
typedef struct {
    double mass;             /* m [kg] */
    double spring_k;         /* k [N/m] */
    double damping_c;        /* c [N·s/m] */
    double displacement;     /* x [m] */
    double velocity;         /* ẋ [m/s] */
    double acceleration;     /* ẍ [m/s²] */
    double force_external;   /* F_ext [N] */
    double time;             /* current simulation time [s] */
} MemsSpringMassDamper;

/**
 * @brief Step the mass-spring-damper by dt using semi-implicit Euler
 *        ẍ = (F_ext - c·ẋ - k·x) / m   (Newton's 2nd law)
 *        ẋ(t+dt) = ẋ(t) + ẍ·dt
 *        x(t+dt) = x(t) + ẋ(t+dt)·dt
 *
 * @param sys The mechanical system state
 * @param dt  Time step [s]
 * Complexity: O(1)
 */
void mems_smd_step(MemsSpringMassDamper *sys, double dt);

/**
 * @brief Compute steady-state amplitude for harmonic forcing
 *        X = F0 / sqrt((k - m·ω²)² + (c·ω)²)
 *        This is the solution to m·ẍ + c·ẋ + k·x = F0·cos(ωt)
 */
double mems_harmonic_amplitude(double mass, double spring_k,
                                double damping_c, double force_amplitude,
                                double angular_freq);

/**
 * @brief Compute phase lag for harmonic forcing
 *        φ = atan2(c·ω, k - m·ω²)
 */
double mems_harmonic_phase(double mass, double spring_k,
                            double damping_c, double angular_freq);

/**
 * @brief Compute the normalized frequency response magnitude
 *        |H(jω)| = 1 / sqrt((1 - r²)² + (2ζr)²)
 *        where r = ω/ω_n is normalized frequency
 */
double mems_freq_response_mag(double r_norm, double zeta);

/**
 * @brief Compute the normalized frequency response phase
 *        ∠H(jω) = -atan2(2ζr, 1 - r²)
 */
double mems_freq_response_phase(double r_norm, double zeta);

/**
 * @brief Compute natural frequency: ω_n = sqrt(k/m) [rad/s]
 */
double mems_natural_freq_rad(double mass, double spring_k);

/**
 * @brief Compute damped natural frequency: ω_d = ω_n·sqrt(1 - ζ²)
 *        Only defined for ζ < 1 (underdamped)
 */
double mems_damped_freq_rad(double wn, double zeta);

/**
 * @brief Compute displacement under DC acceleration: x = m·a / k
 *        This is the static sensitivity of an accelerometer.
 *        Reference: Bao (2005) §3.2
 */
double mems_static_displacement(double mass, double spring_k,
                                 double acceleration);

/* ─── L2: Damping Mechanisms ───────────────────────────────────────────────── */

/**
 * @brief Squeeze-film damping coefficient for parallel plates
 *        c_sf = μ·A / (g³) · β(w/L)
 *        where μ is dynamic viscosity, A is plate area, g is gap,
 *        and β is aspect-ratio correction.
 *        Reference: Bao (2005) §4.3, Blech (1983)
 *
 * @param viscosity Dynamic viscosity of fluid [Pa·s], air ≈ 1.85e-5
 * @param area      Plate area [m²]
 * @param gap       Gap distance [m]
 * @param width     Plate width [m]
 * @param length    Plate length [m]
 */
double mems_squeeze_film_damping(double viscosity, double area,
                                  double gap, double width, double length);

/**
 * @brief Couette (slide-film) damping coefficient
 *        c_couette = μ·A / g
 *        Reference: Senturia (2001) §7.3
 */
double mems_couette_damping(double viscosity, double area, double gap);

/**
 * @brief Thermoelastic damping (TED) — Zener's model
 *        Q_TED⁻¹ = (E·α²·T₀ / ρ·C_p) · (ω·τ / (1 + ω²·τ²))
 *        where τ = b²/(π²·χ) is thermal relaxation time.
 *        Reference: Zener (1937), Lifshitz & Roukes (2000)
 *
 * @param young_modulus   E [Pa]
 * @param thermal_exp     α [1/K]
 * @param temperature     T₀ [K]
 * @param density         ρ [kg/m³]
 * @param specific_heat   C_p [J/(kg·K)]
 * @param freq_rad        ω [rad/s]
 * @param beam_width      b [m], beam width in direction of heat flow
 * @param thermal_diff    χ [m²/s], thermal diffusivity
 */
double mems_thermoelastic_q_inv(double young_modulus, double thermal_exp,
                                 double temperature, double density,
                                 double specific_heat, double freq_rad,
                                 double beam_width, double thermal_diff);

/* ─── L3: Coriolis Effect (Gyroscope Physics) ──────────────────────────────── */

/**
 * @brief Coriolis force magnitude on a MEMS gyro proof mass
 *        F_coriolis = 2 · m · v_drive · Ω_input
 *        where v_drive is the drive velocity and Ω is the input rate.
 *        Reference: IEEE Std 1431-2004 (Coriolis vibratory gyroscopes)
 *
 * @param mass       Proof mass [kg]
 * @param v_drive    Drive-mode velocity [m/s]
 * @param omega_in   Input angular rate [rad/s]
 */
double mems_coriolis_force(double mass, double v_drive, double omega_in);

/**
 * @brief Sense displacement from Coriolis force
 *        x_sense = F_c / (k_sense · sqrt((1-r²)² + (2ζr)²))
 *        where r = ω_drive / ω_sense and ω_sense = ω_drive for mode-matched.
 *        Reference: Acar & Shkel, "MEMS Vibratory Gyroscopes" (2009)
 */
double mems_gyro_sense_displacement(double coriolis_force, double k_sense,
                                     double r_norm, double zeta);

/**
 * @brief Gyro scale factor: SF = x_sense / Ω_input
 *        Mode-matched case: SF = 2·m·v_drive·Q / k_sense
 */
double mems_gyro_scale_factor(double mass, double v_drive, double q_factor,
                               double k_sense);

/* ─── L3: Stoney's Formula (Thin-Film Stress) ──────────────────────────────── */

/**
 * @brief Stoney's formula for thin-film stress from curvature
 *        σ_f = (E_s · t_s²) / (6 · (1-ν_s) · t_f · R)
 *        where E_s is substrate Young's modulus, t_s is substrate thickness,
 *        ν_s is Poisson ratio, t_f is film thickness, R is curvature radius.
 *        Reference: Stoney (1909), Janssen et al. (2009) — extended formula
 *
 * @param E_s    Substrate Young's modulus [Pa]
 * @param t_s    Substrate thickness [m]
 * @param nu_s   Substrate Poisson's ratio [-]
 * @param t_f    Film thickness [m]
 * @param R      Radius of curvature [m], positive = tensile
 */
double mems_stoney_stress(double E_s, double t_s, double nu_s,
                           double t_f, double R);

/**
 * @brief Curvature from known stress (inverse Stoney)
 */
double mems_stoney_curvature(double E_s, double t_s, double nu_s,
                              double t_f, double stress_f);

/* ─── L2: MEMS Resonator Mass Sensing ────────────────────────────────────── */

double mems_mass_loading_sensitivity(double effective_mass);
double mems_mass_resolution(double f_resonant, double q_factor,
                             double effective_mass, double snr);

/* ─── L2: Capacitive Readout ──────────────────────────────────────────────── */

double mems_capacitive_displacement(double c1, double c2, double gap_nominal);
double mems_capacitance_parallel_plate(double epsilon, double area,
                                        double gap, double displacement);
double mems_pullin_voltage(double spring_k, double gap_nominal,
                            double epsilon, double area);

#endif /* MEMS_DYNAMICS_H */
