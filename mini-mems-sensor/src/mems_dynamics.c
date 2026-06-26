/**
 * @file mems_dynamics.c
 * @brief L2-L3 MEMS Mechanical Dynamics — Spring-mass-damper ODE,
 *        frequency response, damping mechanisms, Coriolis gyro physics,
 *        Stoney's thin-film stress formula
 *
 * Reference: Bao (2005) — "Analysis and Design Principles of MEMS Devices"
 *            Senturia (2001) — "Microsystem Design"
 *            Thomson (1996) — "Theory of Vibration with Applications"
 */

#include "mems_dynamics.h"
#include <math.h>
#include <float.h>
#include <string.h>

/* ─── L2: Spring-Mass-Damper Time Evolution ──────────────────────────────────
 * The MEMS inertial sensor is fundamentally a 2nd-order mass-spring-damper.
 *
 * Newton's 2nd law: m·ẍ + c·ẋ + k·x = F_ext(t)
 *
 * For an accelerometer, F_ext(t) = m·a(t) where a(t) is external acceleration.
 * The displacement x(t) is then sensed capacitively.
 *
 * For a gyroscope, Coriolis force couples drive and sense modes:
 *   F_coriolis = 2·m·v_drive × Ω_input
 */

void mems_smd_step(MemsSpringMassDamper *sys, double dt) {
    if (!sys || dt <= 0.0) return;
    if (sys->mass <= 0.0) return;

    /* Newton's 2nd law: ẍ = (F_ext - c·ẋ - k·x) / m */
    sys->acceleration = (sys->force_external -
                          sys->damping_c * sys->velocity -
                          sys->spring_k * sys->displacement) / sys->mass;

    /* Semi-implicit Euler integration */
    sys->velocity     += sys->acceleration * dt;
    sys->displacement += sys->velocity * dt;
    sys->time         += dt;
}

double mems_harmonic_amplitude(double mass, double spring_k,
                                double damping_c, double force_amplitude,
                                double angular_freq) {
    /* Steady-state solution to m·ẍ + c·ẋ + k·x = F₀·cos(ωt)
     * X = F₀ / sqrt((k - m·ω²)² + (c·ω)²) */
    if (mass <= 0.0 || spring_k < 0.0) return -1.0;
    double k_minus_mw2 = spring_k - mass * angular_freq * angular_freq;
    double denom = sqrt(k_minus_mw2 * k_minus_mw2 +
                         damping_c * angular_freq * damping_c * angular_freq);
    if (denom < DBL_EPSILON) return INFINITY; /* resonance with zero damping */
    return force_amplitude / denom;
}

double mems_harmonic_phase(double mass, double spring_k,
                            double damping_c, double angular_freq) {
    /* φ = atan2(c·ω, k - m·ω²) */
    if (mass <= 0.0) return 0.0;
    return atan2(damping_c * angular_freq,
                  spring_k - mass * angular_freq * angular_freq);
}

double mems_freq_response_mag(double r_norm, double zeta) {
    /* |H(jω)| = 1 / sqrt((1 - r²)² + (2ζr)²)
     * where r = ω/ω_n is the frequency ratio */
    if (r_norm < 0.0 || zeta < 0.0) return -1.0;
    double one_minus_r2 = 1.0 - r_norm * r_norm;
    double denom = sqrt(one_minus_r2 * one_minus_r2 +
                         4.0 * zeta * zeta * r_norm * r_norm);
    if (denom < DBL_EPSILON) return INFINITY;
    return 1.0 / denom;
}

double mems_freq_response_phase(double r_norm, double zeta) {
    /* ∠H(jω) = -atan2(2ζr, 1 - r²) */
    return -atan2(2.0 * zeta * r_norm, 1.0 - r_norm * r_norm);
}

double mems_natural_freq_rad(double mass, double spring_k) {
    if (mass <= 0.0 || spring_k <= 0.0) return 0.0;
    return sqrt(spring_k / mass);
}

double mems_damped_freq_rad(double wn, double zeta) {
    if (zeta >= 1.0 || zeta < 0.0) return 0.0;  /* critically/overdamped */
    return wn * sqrt(1.0 - zeta * zeta);
}

double mems_static_displacement(double mass, double spring_k,
                                 double acceleration) {
    if (spring_k <= 0.0) return 0.0;
    return mass * acceleration / spring_k;
}

/* ─── L2: Damping Mechanisms ─────────────────────────────────────────────────
 * MEMS devices use various energy dissipation mechanisms:
 *
 * 1. Squeeze-film damping: gas trapped between parallel plates
 *    Dominant for low-gap capacitive MEMS at atmospheric pressure.
 *    c_sf ∝ μ·A/g³ with aspect ratio correction factor β.
 *
 * 2. Slide-film (Couette) damping: gas sheared between parallel moving plates
 *    c_couette = μ·A/g
 *
 * 3. Thermoelastic damping (TED): internal friction from thermal gradients
 *    Due to thermoelastic coupling: compression→heat→diffusion→loss.
 *    Fundamental limit for vacuum-operated MEMS resonators.
 *    Q_TED ∝ (ρ·C_p) / (E·α²·T₀)
 */

double mems_squeeze_film_damping(double viscosity, double area,
                                  double gap, double width, double length) {
    if (gap <= 0.0 || area <= 0.0 || length <= 0.0) return 0.0;
    /* Aspect ratio factor: β ≈ 1 - 0.21·(w/L) for w/L < 1
     * Reference: Veijola (1995), modified Blech model */
    double ar = (width < length) ? width / length : length / width;
    double beta = 1.0 - 0.21 * ar + 0.09 * ar * ar * ar;
    if (beta < 0.1) beta = 0.1;  /* clamp for extreme aspect ratios */
    return viscosity * area * beta / (gap * gap * gap);
}

double mems_couette_damping(double viscosity, double area, double gap) {
    if (gap <= 0.0) return 0.0;
    return viscosity * area / gap;
}

double mems_thermoelastic_q_inv(double young_modulus, double thermal_exp,
                                 double temperature, double density,
                                 double specific_heat, double freq_rad,
                                 double beam_width, double thermal_diff) {
    /* Zener's thermoelastic damping model:
     * Q⁻¹ = (E·α²·T₀ / (ρ·C_p)) · (ω·τ / (1 + ω²·τ²))
     * where τ = b²/(π²·χ) is the thermal relaxation time.
     *
     * Reference: Zener (1937) "Internal Friction in Solids"
     *            Lifshitz & Roukes (2000) "Thermoelastic Damping in MEMS"
     */
    if (density <= 0.0 || specific_heat <= 0.0 || thermal_diff <= 0.0)
        return 0.0;

    /* Thermal relaxation constant: τ = b² / (π²·χ) */
    double tau = (beam_width * beam_width) / (M_PI * M_PI * thermal_diff);
    if (tau < DBL_EPSILON) return 0.0;

    /* Zener modulus: Δ_E = E·α²·T₀ / (ρ·C_p) */
    double delta_E = young_modulus * thermal_exp * thermal_exp * temperature /
                     (density * specific_heat);

    /* Q⁻¹ = Δ_E · (ω·τ) / (1 + (ω·τ)²) */
    double omega_tau = freq_rad * tau;
    return delta_E * omega_tau / (1.0 + omega_tau * omega_tau);
}

/* ─── L3: Coriolis Effect (Gyroscope Principle) ──────────────────────────────
 *
 * A MEMS vibratory gyroscope operates by detecting Coriolis-induced
 * displacement in the sense direction.
 *
 * The proof mass is driven at resonance in the drive axis (x) at amplitude A_d:
 *   x(t) = A_d · sin(ω_d·t)
 *   v_x(t) = ω_d·A_d · cos(ω_d·t)
 *
 * When the gyro rotates about the z-axis at rate Ω_z, a Coriolis force appears
 * in the y-direction:
 *   F_c = 2 · m · v_x · Ω_z
 *
 * This force excites the sense mode (y-axis), whose displacement is measured.
 * For mode-matched operation (ω_sense = ω_drive), the displacement is:
 *   y = (2 · m · ω_d · A_d / k_sense) · Q · Ω_z
 *
 * Scale factor: SF = y/Ω_z = 2·m·ω_d·A_d·Q / k_sense
 */

double mems_coriolis_force(double mass, double v_drive, double omega_in) {
    /* F_c = 2 · m · v_drive × Ω_input (magnitude for orthogonal axes) */
    return 2.0 * mass * v_drive * omega_in;
}

double mems_gyro_sense_displacement(double coriolis_force, double k_sense,
                                     double r_norm, double zeta) {
    /* Sense displacement amplified by mechanical Q:
     * x_sense = F_c / (k_sense · sqrt((1-r²)² + (2ζr)²)) */
    if (k_sense <= 0.0) return 0.0;
    double denom = k_sense * sqrt((1.0 - r_norm * r_norm) *
                                   (1.0 - r_norm * r_norm) +
                                   4.0 * zeta * zeta * r_norm * r_norm);
    if (denom < DBL_EPSILON) return INFINITY;
    return coriolis_force / denom;
}

double mems_gyro_scale_factor(double mass, double v_drive, double q_factor,
                               double k_sense) {
    /* SF = 2·m·v_drive·Q / k_sense  (mode-matched approximation) */
    if (k_sense <= 0.0) return 0.0;
    return 2.0 * mass * v_drive * q_factor / k_sense;
}

/* ─── L3: Stoney's Formula (Thin-Film Stress Measurement) ────────────────────
 *
 * Stoney's formula relates thin-film residual stress to substrate curvature.
 * This is fundamental for MEMS fabrication process control.
 *
 * Original (Stoney, 1909): σ_f = (E_s · t_s²) / (6 · t_f · R)
 *
 * Extended (Janssen et al., 2009): includes biaxial modulus E_s/(1-ν_s)
 *   σ_f = (E_s · t_s²) / (6 · (1-ν_s) · t_f · R)
 *
 * where:
 *   E_s  = substrate Young's modulus [Pa]
 *   t_s  = substrate thickness [m]
 *   ν_s  = substrate Poisson's ratio
 *   t_f  = film thickness [m]
 *   R    = radius of curvature [m] (positive = tensile stress)
 *
 * Assumptions: t_f ≪ t_s, isotropic elasticity, uniform stress,
 *              small deflection, no delamination.
 */

double mems_stoney_stress(double E_s, double t_s, double nu_s,
                           double t_f, double R) {
    /* σ_f = (E_s · t_s²) / (6 · (1-ν_s) · t_f · R)
     * Returns: stress in [Pa], positive = tensile */
    if (t_s <= 0.0 || t_f <= DBL_EPSILON || fabs(R) < DBL_EPSILON)
        return 0.0;
    if (nu_s >= 1.0 || nu_s <= -1.0) return 0.0;  /* physical Poisson range */
    double biaxial_modulus = E_s / (1.0 - nu_s);
    return biaxial_modulus * t_s * t_s / (6.0 * t_f * R);
}

double mems_stoney_curvature(double E_s, double t_s, double nu_s,
                              double t_f, double stress_f) {
    /* R = (E_s · t_s²) / (6 · (1-ν_s) · t_f · σ_f)
     * Inverse of Stoney's formula */
    if (t_s <= 0.0 || t_f <= DBL_EPSILON || fabs(stress_f) < DBL_EPSILON)
        return INFINITY;  /* zero stress = flat wafer */
    if (nu_s >= 1.0 || nu_s <= -1.0) return 0.0;
    double biaxial_modulus = E_s / (1.0 - nu_s);
    return biaxial_modulus * t_s * t_s / (6.0 * t_f * stress_f);
}

/* ─── L2: MEMS Resonator Analysis ────────────────────────────────────────────
 *
 * MEMS resonators are used as timing references (replacing quartz crystals),
 * RF filters (BAW/SAW), and mass sensors (gravimetric).
 *
 * The resonance frequency depends on effective mass and stiffness:
 *   f_n = (1/2π) · sqrt(k_eff / m_eff)
 *
 * Mass loading sensitivity (for gravimetric sensing):
 *   Δf/f = -Δm / (2·m_eff)
 *
 * Reference: Sauerbrey (1959), "Use of Quartz Crystal Vibrators for Thin Films"
 *            Ekinci et al. (2004), "Ultrasensitive NEMS Mass Sensors"
 */

double mems_mass_loading_sensitivity(double effective_mass) {
    /* Fractional frequency change per unit mass: Δf/f per kg
     * S_mass = -1 / (2·m_eff)  [1/kg] */
    if (effective_mass <= 0.0) return 0.0;
    return -1.0 / (2.0 * effective_mass);
}

double mems_mass_resolution(double f_resonant, double q_factor,
                             double effective_mass, double snr) {
    /* Minimum detectable mass from frequency noise:
     * δm_min = (2·m_eff / Q) · (1/SNR)
     * The Allan deviation minimum at f_resonant gives the frequency stability.
     * Reference: Ekinci et al., J. Appl. Phys. (2004) */
    (void)f_resonant;  /* used in full expression with Allan deviation */
    if (q_factor <= 0.0 || snr <= 0.0) return INFINITY;
    return 2.0 * effective_mass / (q_factor * snr);
}

/* ─── L2: Capacitive Readout Analysis ────────────────────────────────────────
 * Most MEMS inertial sensors use capacitive detection:
 * differential parallel-plate capacitors measuring proof mass displacement.
 *
 * For a single capacitor: C = ε₀·A/(d₀ - x)
 * For differential pair: C₁ = ε₀·A/(d₀ - x), C₂ = ε₀·A/(d₀ + x)
 *
 * The difference normalized by sum gives linear displacement:
 * (C₁ - C₂)/(C₁ + C₂) = x/d₀  (exact, not just small-displacement!)
 *
 * Reference: Boser & Howe (1996), "Surface Micromachined Accelerometers"
 */

double mems_capacitive_displacement(double c1, double c2, double gap_nominal) {
    /* x = d₀ · (C₁ - C₂) / (C₁ + C₂) */
    if (gap_nominal <= 0.0) return 0.0;
    double sum = c1 + c2;
    if (sum < DBL_EPSILON) return 0.0;
    return gap_nominal * (c1 - c2) / sum;
}

double mems_capacitance_parallel_plate(double epsilon, double area,
                                        double gap, double displacement) {
    /* C = ε·A / (d - x), with x toward decreasing gap */
    double effective_gap = gap - displacement;
    if (effective_gap <= DBL_EPSILON) return INFINITY;  /* pull-in */
    return epsilon * area / effective_gap;
}

double mems_pullin_voltage(double spring_k, double gap_nominal,
                            double epsilon, double area) {
    /* Electrostatic pull-in voltage for parallel-plate actuator:
     * V_pi = sqrt(8·k·d₀³ / (27·ε·A))
     * At V_pi, the plate snaps down to 2·d₀/3.
     * Reference: Senturia (2001) §6.3 */
    if (spring_k <= 0.0 || gap_nominal <= 0.0 || area <= 0.0) return INFINITY;
    double numerator = 8.0 * spring_k * gap_nominal * gap_nominal * gap_nominal;
    double denominator = 27.0 * epsilon * area;
    if (denominator <= 0.0) return INFINITY;
    return sqrt(numerator / denominator);
}
