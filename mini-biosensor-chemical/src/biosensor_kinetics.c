/**
 * biosensor_kinetics.c — Biorecognition kinetics and thermodynamic models
 *
 * L3: Michaelis-Menten, Hill, Langmuir, Fickian diffusion,
 *     antibody-antigen binding.
 * L4: Michaelis-Menten equation, Briggs-Haldane steady-state,
 *     Hill equation for cooperativity, Langmuir isotherm,
 *     Fick's laws of diffusion, mass-action law.
 * L5: Lineweaver-Burk, Eadie-Hofstee, Hanes-Woolf linearizations,
 *     nonlinear least-squares parameter estimation,
 *     Hill plot cooperativity analysis.
 * L6: Competitive/noncompetitive/uncompetitive inhibition,
 *     enzyme kinetics parameter extraction, SPR binding curves.
 * L7: Immunoassay binding kinetics for clinical diagnostics.
 */

#include "biosensor_kinetics.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * L4: Michaelis-Menten equation
 * ============================================================================ */

/**
 * Michaelis-Menten: v = V_max · [S] / (K_m + [S])
 *
 * Derived under the steady-state assumption (Briggs-Haldane, 1925):
 *   d[ES]/dt = k₁·[E][S] - (k₋₁ + k_cat)·[ES] ≈ 0
 *   K_m = (k₋₁ + k_cat) / k₁
 *   V_max = k_cat · [E]₀
 *
 * K_m has units of concentration; it equals the [S] at which v = V_max/2.
 */
double michaelis_menten_rate(const MichaelisMentenParams *params)
{
    if (!params) return 0.0;
    if (params->km <= 0.0 || params->substrate_conc < 0.0) return 0.0;

    return params->vmax * params->substrate_conc /
           (params->km + params->substrate_conc);
}

/**
 * Inverse Michaelis-Menten: [S] = K_m · v / (V_max - v)
 */
double michaelis_menten_inverse(const MichaelisMentenParams *params,
                                double rate)
{
    if (!params || params->km <= 0.0) return 0.0;
    if (rate < 0.0) return 0.0;
    if (rate >= params->vmax) return (params->km * rate) / 1.0e-10;  /* large */

    return params->km * rate / (params->vmax - rate);
}

/* ============================================================================
 * L4/L6: Enzyme inhibition types
 * ============================================================================ */

/**
 * Competitive inhibition:
 * Apparent K_m' = K_m · (1 + [I]/K_i), V_max unchanged.
 */
double competitive_inhibition_rate(const MichaelisMentenParams *params,
                                   double inhibitor_conc,
                                   double inhibition_constant)
{
    if (!params || params->km <= 0.0 || inhibition_constant <= 0.0)
        return michaelis_menten_rate(params);

    double km_apparent = params->km * (1.0 + inhibitor_conc / inhibition_constant);
    double s = params->substrate_conc;

    return params->vmax * s / (km_apparent + s);
}

/**
 * Uncompetitive inhibition:
 * Inhibitor binds only to ES complex → both K_m and V_max decrease.
 *
 * Apparent V_max' = V_max / (1 + [I]/K_i)
 * Apparent K_m'  = K_m  / (1 + [I]/K_i)
 */
double uncompetitive_inhibition_rate(const MichaelisMentenParams *params,
                                     double inhibitor_conc,
                                     double inhibition_constant)
{
    if (!params || params->km <= 0.0 || inhibition_constant <= 0.0)
        return michaelis_menten_rate(params);

    double factor = 1.0 + inhibitor_conc / inhibition_constant;
    double vmax_apparent = params->vmax / factor;
    double km_apparent   = params->km  / factor;
    double s = params->substrate_conc;

    return vmax_apparent * s / (km_apparent + s);
}

/**
 * Non-competitive inhibition:
 * Inhibitor binds to both E and ES with equal affinity.
 * V_max decreases: V_max' = V_max / (1 + [I]/K_i)
 * K_m unchanged.
 */
double noncompetitive_inhibition_rate(const MichaelisMentenParams *params,
                                      double inhibitor_conc,
                                      double inhibition_constant)
{
    if (!params || params->km <= 0.0 || inhibition_constant <= 0.0)
        return michaelis_menten_rate(params);

    double vmax_apparent = params->vmax / (1.0 + inhibitor_conc / inhibition_constant);
    double s = params->substrate_conc;

    return vmax_apparent * s / (params->km + s);
}

/* ============================================================================
 * L5: Linearization methods for K_m and V_max estimation
 * ============================================================================ */

/**
 * Lineweaver-Burk: 1/v = (K_m/V_max)·(1/[S]) + 1/V_max
 *
 * Linear regression: y = m·x + b
 * where m = K_m/V_max, b = 1/V_max
 * → V_max = 1/b, K_m = m/b
 *
 * Known issue: strongly weights low-[S] data points (1/v is large).
 */
double lineweaver_burk_estimate(const double *concentrations,
                                const double *rates,
                                int n_points,
                                double *km_est,
                                double *vmax_est)
{
    if (!concentrations || !rates || !km_est || !vmax_est || n_points < 3) {
        if (km_est)   *km_est   = 0.0;
        if (vmax_est) *vmax_est = 0.0;
        return 0.0;
    }

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0, sum_yy = 0.0;
    int valid = 0;

    for (int i = 0; i < n_points; i++) {
        if (concentrations[i] <= 0.0 || rates[i] <= 0.0) continue;
        double x = 1.0 / concentrations[i];  /* 1/[S] */
        double y = 1.0 / rates[i];            /* 1/v */
        sum_x  += x;
        sum_y  += y;
        sum_xy += x * y;
        sum_xx += x * x;
        sum_yy += y * y;
        valid++;
    }

    if (valid < 2) {
        *km_est = 0.0; *vmax_est = 0.0;
        return 0.0;
    }

    double n = (double)valid;
    double denom = n * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1.0e-30) {
        *km_est = 0.0; *vmax_est = 0.0;
        return 0.0;
    }

    double slope     = (n * sum_xy - sum_x * sum_y) / denom;
    double intercept = (sum_y - slope * sum_x) / n;

    if (intercept <= 0.0) {
        *km_est = 0.0; *vmax_est = 0.0;
        return 0.0;
    }

    *vmax_est = 1.0 / intercept;
    *km_est   = slope / intercept;

    /* R² computation */
    double ss_res = 0.0, ss_tot = 0.0;
    double y_mean = sum_y / n;
    for (int i = 0; i < n_points; i++) {
        if (concentrations[i] <= 0.0 || rates[i] <= 0.0) continue;
        double x = 1.0 / concentrations[i];
        double y = 1.0 / rates[i];
        double y_pred = slope * x + intercept;
        ss_res += (y - y_pred) * (y - y_pred);
        ss_tot += (y - y_mean) * (y - y_mean);
    }
    return (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;
}

/**
 * Eadie-Hofstee: v = -K_m · (v/[S]) + V_max
 *
 * Plot v vs v/[S]. Slope = -K_m, y-intercept = V_max.
 * Both axes contain the dependent variable v, reducing the
 * heteroscedasticity issue of Lineweaver-Burk.
 */
double eadie_hofstee_estimate(const double *concentrations,
                              const double *rates,
                              int n_points,
                              double *km_est,
                              double *vmax_est)
{
    if (!concentrations || !rates || !km_est || !vmax_est || n_points < 3) {
        if (km_est) *km_est = 0.0;
        if (vmax_est) *vmax_est = 0.0;
        return 0.0;
    }

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0, sum_yy = 0.0;
    int valid = 0;

    for (int i = 0; i < n_points; i++) {
        if (concentrations[i] <= 0.0) continue;
        double x = rates[i] / concentrations[i];  /* v/[S] */
        double y = rates[i];                       /* v */
        sum_x  += x;
        sum_y  += y;
        sum_xy += x * y;
        sum_xx += x * x;
        sum_yy += y * y;
        valid++;
    }

    if (valid < 2) return 0.0;
    double n = (double)valid;
    double denom = n * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1.0e-30) return 0.0;

    double slope     = (n * sum_xy - sum_x * sum_y) / denom;  /* = -K_m */
    double intercept = (sum_y - slope * sum_x) / n;            /* = V_max */

    *vmax_est = intercept;
    *km_est   = -slope;
    if (*km_est < 0.0) *km_est = 0.0;

    double y_mean = sum_y / n;
    double ss_res = 0.0, ss_tot = 0.0;
    for (int i = 0; i < n_points; i++) {
        if (concentrations[i] <= 0.0) continue;
        double y_pred = slope * (rates[i] / concentrations[i]) + intercept;
        ss_res += (rates[i] - y_pred) * (rates[i] - y_pred);
        ss_tot += (rates[i] - y_mean) * (rates[i] - y_mean);
    }
    return (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;
}

/**
 * Hanes-Woolf: [S]/v = (1/V_max)·[S] + K_m/V_max
 *
 * The independent variable [S] is error-free in the experimental
 * design, making this the most robust of the three linearizations.
 */
double hanes_woolf_estimate(const double *concentrations,
                            const double *rates,
                            int n_points,
                            double *km_est,
                            double *vmax_est)
{
    if (!concentrations || !rates || !km_est || !vmax_est || n_points < 3) {
        if (km_est) *km_est = 0.0;
        if (vmax_est) *vmax_est = 0.0;
        return 0.0;
    }

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0, sum_yy = 0.0;
    int valid = 0;

    for (int i = 0; i < n_points; i++) {
        if (rates[i] <= 0.0) continue;
        double x = concentrations[i];           /* [S] */
        double y = concentrations[i] / rates[i]; /* [S]/v */
        sum_x  += x;
        sum_y  += y;
        sum_xy += x * y;
        sum_xx += x * x;
        sum_yy += y * y;
        valid++;
    }

    if (valid < 2) return 0.0;
    double n = (double)valid;
    double denom = n * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1.0e-30) return 0.0;

    double slope     = (n * sum_xy - sum_x * sum_y) / denom;  /* = 1/V_max */
    double intercept = (sum_y - slope * sum_x) / n;            /* = K_m/V_max */

    if (slope <= 0.0) return 0.0;
    *vmax_est = 1.0 / slope;
    *km_est   = intercept / slope;

    double y_mean = sum_y / n;
    double ss_res = 0.0, ss_tot = 0.0;
    for (int i = 0; i < n_points; i++) {
        if (rates[i] <= 0.0) continue;
        double y_pred = slope * concentrations[i] + intercept;
        ss_res += (concentrations[i] / rates[i] - y_pred) *
                  (concentrations[i] / rates[i] - y_pred);
        ss_tot += (concentrations[i] / rates[i] - y_mean) *
                  (concentrations[i] / rates[i] - y_mean);
    }
    return (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;
}

/**
 * Direct nonlinear least-squares (Gauss-Newton) estimation of K_m, V_max.
 *
 * Minimizes: Σ (v_i - V_max·[S_i]/(K_m + [S_i]))²
 *
 * Jacobian:
 *   df/dK_m   = -V_max · [S] / (K_m + [S])²
 *   df/dV_max = [S] / (K_m + [S])
 *
 * Gauss-Newton update: θ_{k+1} = θ_k - (J^T J)⁻¹ J^T r
 */
int michaelis_menten_nls_fit(const double *concentrations,
                             const double *rates,
                             int n_points,
                             double km_init, double vmax_init,
                             double *km_est, double *vmax_est,
                             int max_iterations, double tolerance)
{
    if (!concentrations || !rates || !km_est || !vmax_est ||
        n_points < 3 || km_init <= 0.0 || vmax_init <= 0.0 ||
        max_iterations <= 0) {
        if (km_est)   *km_est   = km_init;
        if (vmax_est) *vmax_est = vmax_init;
        return 0;
    }

    double km   = km_init;
    double vmax = vmax_init;

    for (int iter = 0; iter < max_iterations; iter++) {
        double j11 = 0.0, j12 = 0.0, j22 = 0.0;  /* J^T J elements */
        double g1 = 0.0, g2 = 0.0;                /* J^T r (gradient) */
        double ss_res = 0.0;

        for (int i = 0; i < n_points; i++) {
            double s = concentrations[i];
            if (s <= 0.0) continue;

            double denom = km + s;
            double denom_sq = denom * denom;

            double pred = vmax * s / denom;
            double residual = rates[i] - pred;

            /* df/dK_m = -Vmax * s / (Km + s)^2 */
            double dK = -vmax * s / denom_sq;
            /* df/dV_max = s / (Km + s) */
            double dV = s / denom;

            j11 += dK * dK;
            j12 += dK * dV;
            j22 += dV * dV;

            g1 += dK * residual;
            g2 += dV * residual;

            ss_res += residual * residual;
        }

        /* Solve: [j11 j12; j12 j22] · Δθ = -[g1; g2] */
        double det = j11 * j22 - j12 * j12;
        if (fabs(det) < 1.0e-30) break;

        double delta_km   = -(j22 * g1 - j12 * g2) / det;
        double delta_vmax = -(-j12 * g1 + j11 * g2) / det;

        /* Levenberg-Marquardt damping for robustness */
        double damping = 1.0;
        double new_km   = km   + damping * delta_km;
        double new_vmax = vmax + damping * delta_vmax;

        /* Enforce positivity */
        if (new_km   < 1.0e-12) new_km   = km   * 0.5;
        if (new_vmax < 1.0e-12) new_vmax = vmax * 0.5;

        /* Check convergence */
        double rel_change = fabs(delta_km / km) + fabs(delta_vmax / vmax);
        km   = new_km;
        vmax = new_vmax;

        if (rel_change < tolerance) {
            *km_est   = km;
            *vmax_est = vmax;
            return iter + 1;
        }
    }

    *km_est   = km;
    *vmax_est = vmax;
    return max_iterations;
}

/* ============================================================================
 * L4: Hill equation — cooperative binding
 * ============================================================================ */

/**
 * Hill equation: θ = [L]^n / (K_d^n + [L]^n)
 *
 * When n > 1: positive cooperativity (binding of one ligand
 *   facilitates binding of subsequent ligands — hemoglobin, O₂).
 * When n = 1: non-cooperative (hyperbolic — equivalent to Langmuir/M-M).
 * When n < 1: negative cooperativity (binding hinders further binding).
 */
double hill_fractional_saturation(const HillParams *hp, double ligand_conc)
{
    if (!hp) return 0.0;
    if (ligand_conc <= 0.0) return hp->baseline;

    double kd = hp->half_max_conc;
    if (kd <= 0.0) return hp->max_response;

    double conc_pow = pow(ligand_conc, hp->hill_coefficient);
    double kd_pow   = pow(kd, hp->hill_coefficient);

    double theta = conc_pow / (kd_pow + conc_pow);
    return hp->baseline + (hp->max_response - hp->baseline) * theta;
}

/**
 * Hill plot: log(θ/(1-θ)) = n · log[L] - n · log(K_d)
 *
 * This transformation linearizes cooperative binding data.
 * The slope is the Hill coefficient n, and the x-intercept gives log(K_d).
 */
double hill_plot_estimate(const double *concentrations,
                          const double *saturations,
                          int n_points,
                          double *n_hill,
                          double *log_kd)
{
    if (!concentrations || !saturations || !n_hill || !log_kd || n_points < 3) {
        if (n_hill) *n_hill = 1.0;
        if (log_kd) *log_kd = 0.0;
        return 0.0;
    }

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0, sum_yy = 0.0;
    int valid = 0;

    for (int i = 0; i < n_points; i++) {
        if (concentrations[i] <= 0.0) continue;
        double theta = saturations[i];
        /* Clamp to avoid log(0) and division by zero */
        if (theta <= 0.01 || theta >= 0.99) continue;

        double x = log10(concentrations[i]);
        double y = log10(theta / (1.0 - theta));

        sum_x  += x;
        sum_y  += y;
        sum_xy += x * y;
        sum_xx += x * x;
        sum_yy += y * y;
        valid++;
    }

    if (valid < 2) return 0.0;
    double n = (double)valid;
    double denom = n * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1.0e-30) return 0.0;

    double slope     = (n * sum_xy - sum_x * sum_y) / denom;  /* n_Hill */
    double intercept = (sum_y - slope * sum_x) / n;           /* -n·log(K_d) */

    *n_hill = slope;
    *log_kd = -intercept / slope;  /* log10(K_d) */

    double y_mean = sum_y / n;
    double ss_res = 0.0, ss_tot = 0.0;
    for (int i = 0; i < n_points; i++) {
        if (concentrations[i] <= 0.0) continue;
        double theta = saturations[i];
        if (theta <= 0.01 || theta >= 0.99) continue;
        double y_pred = slope * log10(concentrations[i]) + intercept;
        double y_meas = log10(theta / (1.0 - theta));
        ss_res += (y_meas - y_pred) * (y_meas - y_pred);
        ss_tot += (y_meas - y_mean) * (y_meas - y_mean);
    }
    return (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;
}

/* ============================================================================
 * L4: Langmuir adsorption isotherm
 * ============================================================================ */

/**
 * Langmuir isotherm: θ = K·c / (1 + K·c)
 *
 * Assumptions:
 * 1. Adsorption occurs at specific, identical sites
 * 2. Each site holds exactly one molecule (monolayer)
 * 3. No lateral interactions between adsorbed molecules
 * 4. Reversible adsorption-desorption equilibrium
 */
double langmuir_coverage(const LangmuirParams *lp, double bulk_conc)
{
    if (!lp || bulk_conc < 0.0) return 0.0;
    if (lp->equilibrium_constant <= 0.0) return 0.0;

    double kc = lp->equilibrium_constant * bulk_conc;
    return kc / (1.0 + kc);
}

/**
 * Langmuir kinetics: θ(t) = θ_eq · (1 - exp(-k_obs · t))
 *
 * where k_obs = k_ads·c + k_des, θ_eq = k_ads·c / k_obs
 *
 * This is the analytical solution to the first-order linear ODE
 * describing time-dependent adsorption under constant bulk concentration.
 */
double langmuir_kinetic_coverage(const LangmuirParams *lp,
                                 double bulk_conc,
                                 double time_sec)
{
    if (!lp || bulk_conc < 0.0 || time_sec < 0.0) return 0.0;
    if (lp->k_adsorption <= 0.0) return 0.0;

    double k_obs = lp->k_adsorption * bulk_conc + lp->k_desorption;
    if (k_obs <= 0.0) return 0.0;

    double theta_eq = lp->k_adsorption * bulk_conc / k_obs;
    return theta_eq * (1.0 - exp(-k_obs * time_sec));
}

/**
 * Langmuir dissociation: θ(t) = θ₀ · exp(-k_des · t)
 *
 * After the bulk solution is replaced with pure buffer (c ≈ 0),
 * the bound analyte desorbs with first-order kinetics.
 */
double langmuir_dissociation(const LangmuirParams *lp,
                             double initial_coverage,
                             double time_sec)
{
    if (!lp || initial_coverage <= 0.0 || time_sec < 0.0) return 0.0;
    if (lp->k_desorption <= 0.0) return initial_coverage;

    return initial_coverage * exp(-lp->k_desorption * time_sec);
}

/* ============================================================================
 * L4: Fickian diffusion
 * ============================================================================ */

/**
 * Fick's 2nd law solution for semi-infinite medium:
 *   C(x,t) = C₀ · erfc(x / √(4Dt))
 *
 * erfc(z) = 1 - erf(z) = (2/√π) ∫_z^∞ exp(-u²) du
 *
 * Approximation using the Abramowitz-Stegun 7.1.26 formula
 * for the error function (fractional error < 1.5×10⁻⁷).
 */
double fickian_concentration(double bulk_conc,
                             double diffusion_coeff,
                             double distance_cm,
                             double time_sec)
{
    if (diffusion_coeff <= 0.0 || time_sec <= 0.0) return bulk_conc;
    if (distance_cm < 0.0) distance_cm = 0.0;

    double z = distance_cm / sqrt(4.0 * diffusion_coeff * time_sec);

    /* Abramowitz & Stegun erf approximation */
    double t = 1.0 / (1.0 + 0.3275911 * z);
    double erf_z = 1.0 - (0.254829592 * t -
                          0.284496736 * t * t +
                          1.421413741 * t * t * t -
                          1.453152027 * t * t * t * t +
                          1.061405429 * t * t * t * t * t) * exp(-z * z);

    if (z > 5.0) erf_z = 1.0;  /* erfc ≈ 0 for z > 5 */

    return bulk_conc * (1.0 - erf_z);  /* C₀ · erfc(z) */
}

/**
 * Damköhler number: Da = (k_cat · Γ) / (D · C₀ / δ)
 *
 * Interprets the balance between reaction and diffusion:
 * - Da << 1: reaction-limited (sensor responds to bulk concentration)
 * - Da >> 1: diffusion-limited (sensor depletes local analyte)
 * - Da ≈ 1: mixed control (most realistic for biosensors)
 */
double damkohler_number(double kcat,
                        double surface_enzyme_density,
                        double diffusion_coeff,
                        double bulk_conc,
                        double diffusion_layer_thickness)
{
    if (diffusion_coeff <= 0.0 || diffusion_layer_thickness <= 0.0 ||
        bulk_conc <= 0.0) return 0.0;

    double reaction_rate = kcat * surface_enzyme_density;
    double mass_transfer_rate = diffusion_coeff * bulk_conc /
                                diffusion_layer_thickness;

    if (mass_transfer_rate <= 0.0) return INFINITY;
    return reaction_rate / mass_transfer_rate;
}

/* ============================================================================
 * L7: Immunoassay binding kinetics
 * ============================================================================ */

/**
 * Antibody-antigen equilibrium binding (mass action law).
 *
 * K_d = [Ab][Ag] / [Ab·Ag]
 *
 * For [Ab]₀ total antibody binding sites and [Ag] total antigen:
 *   [Ab·Ag] = 0.5 · ([Ab]₀ + [Ag] + K_d -
 *              √(([Ab]₀ + [Ag] + K_d)² - 4[Ab]₀[Ag]))
 *
 * This accounts for antigen depletion when [Ab]₀ is not >> [Ag].
 *
 * References: Feldman (1972) Mathematical theory of complex
 *   ligand-binding systems.
 *
 * L7 Application: Clinical immunoassay quantification.
 */
double antibody_antigen_complex(double kd,
                                double antigen_conc,
                                double antibody_conc)
{
    if (kd < 0.0 || antigen_conc < 0.0 || antibody_conc <= 0.0) return 0.0;

    if (antigen_conc <= 0.0) return 0.0;

    /* Full quadratic solution for antigen-antibody equilibrium */
    double a = 1.0;
    double b = -(antibody_conc + antigen_conc + kd);
    double c = antibody_conc * antigen_conc;

    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0) discriminant = 0.0;

    /* The physically meaningful root is the smaller one
     * (complex concentration cannot exceed limiting reagent) */
    double sqrt_disc = sqrt(discriminant);
    double root1 = (-b - sqrt_disc) / (2.0 * a);
    double root2 = (-b + sqrt_disc) / (2.0 * a);

    double complex_conc = (root1 >= 0.0 && root1 <= antibody_conc &&
                           root1 <= antigen_conc) ? root1 : root2;

    if (complex_conc < 0.0) complex_conc = 0.0;
    if (complex_conc > antibody_conc) complex_conc = antibody_conc;
    if (complex_conc > antigen_conc)  complex_conc = antigen_conc;

    return complex_conc;
}

/**
 * Binding association curve: R(t) = R_eq · (1 - exp(-k_obs · t))
 *
 * with k_obs = k_on · [A] + k_off, R_eq = R_max · [A]/(K_d + [A])
 * where K_d = k_off / k_on.
 *
 * This is the standard model for SPR and BLI sensorgrams.
 */
double binding_association_curve(double k_on, double k_off,
                                 double analyte_conc, double r_max,
                                 double time_sec)
{
    if (k_on <= 0.0 || analyte_conc < 0.0 || time_sec < 0.0) return 0.0;

    double kd = (k_on > 0.0) ? (k_off / k_on) : INFINITY;
    double r_eq = r_max * analyte_conc / (kd + analyte_conc);

    double k_obs = k_on * analyte_conc + k_off;
    if (k_obs <= 0.0) return 0.0;

    return r_eq * (1.0 - exp(-k_obs * time_sec));
}
