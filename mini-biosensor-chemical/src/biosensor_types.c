/**
 * biosensor_types.c — Constructor/destructor and utility implementations
 *
 * L1: All type initialization functions for biosensor module.
 * Each init function sets sensible defaults matching real-world
 * biosensor specifications.
 */

#include "biosensor_types.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Physical constants */
#define FARADAY_CONSTANT        96485.3329   /* F [C/mol] */
#define GAS_CONSTANT            8.314462618  /* R [J/(mol·K)] */
#define STANDARD_TEMPERATURE_K  298.15       /* 25 °C */
#define QUARTZ_DENSITY_GCM3     2.648        /* AT-cut quartz */
#define QUARTZ_SHEAR_MODULUS    2.947e11     /* μ_q [g/(cm·s²)] */

void biosensor_descriptor_init(BiosensorDescriptor *desc)
{
    if (!desc) return;
    memset(desc, 0, sizeof(BiosensorDescriptor));
    desc->bioreceptor_type   = BIORECEPTOR_ENZYME;
    desc->transducer_type    = TRANSDUCER_AMPEROMETRIC;
    desc->measurement_mode   = MODE_KINETIC;
    desc->ph_optimum         = 7.0;
    desc->temp_optimum_celsius = 25.0;
    desc->response_time_sec  = 30.0;
    desc->recovery_time_sec  = 60.0;
    desc->shelf_life_days    = 180;
    desc->reuse_cycles       = 1;
    desc->selectivity_coefficient = 0.01; /* 1% cross-reactivity */
}

void michaelis_menten_params_init(MichaelisMentenParams *mm,
                                   double km, double vmax)
{
    if (!mm) return;
    memset(mm, 0, sizeof(MichaelisMentenParams));
    mm->km   = (km > 0.0) ? km : 1.0e-3;   /* default 1 mM */
    mm->vmax = (vmax > 0.0) ? vmax : 1.0;
    mm->kcat = 0.0;        /* must be set from Vmax/[E0] */
    mm->substrate_conc = 0.0;
    mm->enzyme_conc    = 1.0e-9;  /* 1 nM typical enzyme loading */
}

void hill_params_init(HillParams *hp, double k_half, double n_h, double rmax)
{
    if (!hp) return;
    memset(hp, 0, sizeof(HillParams));
    hp->half_max_conc  = (k_half > 0.0) ? k_half : 1.0e-6;   /* 1 μM */
    hp->hill_coefficient = (n_h > 0.0) ? n_h : 1.0;          /* non-cooperative default */
    hp->max_response   = (rmax > 0.0) ? rmax : 100.0;
    hp->baseline       = 0.0;
}

void nernst_params_init(NernstParams *np, double e0, double n)
{
    if (!np) return;
    memset(np, 0, sizeof(NernstParams));
    np->standard_potential = e0;
    np->electron_count     = (n > 0.0) ? n : 1.0;
    np->temperature_kelvin = STANDARD_TEMPERATURE_K;
    np->faraday_constant   = FARADAY_CONSTANT;
    np->gas_constant       = GAS_CONSTANT;
}

void butler_volmer_params_init(ButlerVolmerParams *bv, double i0, double alpha_a)
{
    if (!bv) return;
    memset(bv, 0, sizeof(ButlerVolmerParams));
    bv->exchange_current_density = (i0 > 0.0) ? i0 : 1.0e-4;  /* 0.1 mA/cm² */
    bv->anodic_transfer_coeff   = alpha_a;
    bv->cathodic_transfer_coeff = 1.0 - alpha_a;
    bv->electrode_area = 0.0707;   /* 3 mm dia = 0.0707 cm² */
    bv->overpotential  = 0.0;
}

void eis_params_init(EISParams *eis)
{
    if (!eis) return;
    memset(eis, 0, sizeof(EISParams));
    eis->solution_resistance       = 100.0;   /* 100 Ω typical PBS */
    eis->charge_transfer_resistance = 1000.0; /* 1 kΩ bare electrode */
    eis->double_layer_capacitance  = 1.0e-6;  /* 1 μF */
    eis->warburg_coefficient       = 500.0;   /* 500 Ω/√s */
    eis->constant_phase_exponent   = 0.8;     /* slightly non-ideal */
    eis->frequency_min_hz = 0.1;
    eis->frequency_max_hz = 1.0e5;
}

void isfet_params_init(ISFETParams *isfet)
{
    if (!isfet) return;
    memset(isfet, 0, sizeof(ISFETParams));
    isfet->threshold_voltage   = 1.5;    /* V_th at pH 7 */
    isfet->sensitivity_mv_per_ph = 58.0; /* Nernstian = 59.16 at 25°C */
    isfet->transconductance    = 100.0;  /* 100 μS */
    isfet->drain_current       = 100.0;  /* 100 μA */
    isfet->reference_voltage   = 2.5;    /* V */
}

void spr_params_init(SPRParams *spr, double wavelength_nm, double ri_sens)
{
    if (!spr) return;
    memset(spr, 0, sizeof(SPRParams));
    spr->wavelength_nm        = wavelength_nm;
    spr->refractive_index_sensitivity = (ri_sens > 0.0) ? ri_sens : 1.0e5;
    spr->penetration_depth_nm = 200.0;  /* typical Au at 633 nm */
    spr->gold_layer_thickness_nm = 50.0;
    spr->prism_ri             = 1.52;   /* BK7 glass at 633 nm */
    spr->resonance_angle_deg  = 0.0;    /* to be calculated */
}

void fluorescence_params_init(FluorescenceParams *fp, double ex_nm, double em_nm)
{
    if (!fp) return;
    memset(fp, 0, sizeof(FluorescenceParams));
    fp->excitation_wavelength_nm = ex_nm;
    fp->emission_wavelength_nm   = em_nm;
    fp->quantum_yield       = 0.5;
    fp->extinction_coefficient = 50000.0;  /* typical fluorophore ε */
    fp->stokes_shift_nm     = em_nm - ex_nm;
    fp->photobleaching_halflife_sec = 600.0;  /* 10 minutes */
    fp->forster_radius_nm   = 5.0;   /* typical FRET R₀ */
}

void qcm_params_init(QCMParameters *qcm, double f0_hz, double area_cm2)
{
    if (!qcm) return;
    memset(qcm, 0, sizeof(QCMParameters));
    qcm->fundamental_frequency_hz = (f0_hz > 0.0) ? f0_hz : 10.0e6; /* 10 MHz */
    qcm->electrode_area_cm2  = (area_cm2 > 0.0) ? area_cm2 : 0.2;   /* 5 mm dia */
    qcm->density_quartz_g_per_cm3 = QUARTZ_DENSITY_GCM3;
    qcm->shear_modulus_gpa  = QUARTZ_SHEAR_MODULUS;
    /* Sauerbrey mass sensitivity: Δf = -2f₀²/(A·√(ρ_q·μ_q)) · Δm */
    {
        double sqrt_rho_mu = sqrt(QUARTZ_DENSITY_GCM3 * QUARTZ_SHEAR_MODULUS);
        qcm->mass_sensitivity_hz_per_ng =
            2.0 * f0_hz * f0_hz / (area_cm2 * sqrt_rho_mu) * 1.0e-9; /* Hz/ng */
    }
    qcm->quality_factor       = 5000.0;
    qcm->motional_resistance_ohm = 10.0;
}

void measurement_series_init(MeasurementSeries *series, int32_t capacity)
{
    if (!series) return;
    memset(series, 0, sizeof(MeasurementSeries));
    series->capacity = (capacity > 0) ? capacity : 100;
    series->points = (MeasurementPoint *)calloc((size_t)series->capacity,
                                                 sizeof(MeasurementPoint));
    series->point_count = 0;
    series->temperature_celsius = 25.0;
    series->ph_value = 7.0;
}

void measurement_series_free(MeasurementSeries *series)
{
    if (!series) return;
    free(series->points);
    series->points = NULL;
    series->capacity = 0;
    series->point_count = 0;
}

void calibration_standards_init(CalibrationStandards *cs,
                                 int32_t n_standards, int32_t reps)
{
    if (!cs) return;
    memset(cs, 0, sizeof(CalibrationStandards));
    if (n_standards > 0) {
        cs->concentrations = (double *)calloc((size_t)n_standards, sizeof(double));
        cs->signals        = (double *)calloc((size_t)n_standards, sizeof(double));
        cs->signal_stds    = (double *)calloc((size_t)n_standards, sizeof(double));
    }
    cs->standard_count = n_standards;
    cs->replicates     = (reps > 0) ? reps : 1;
}

void calibration_standards_free(CalibrationStandards *cs)
{
    if (!cs) return;
    free(cs->concentrations);
    free(cs->signals);
    free(cs->signal_stds);
    cs->concentrations = NULL;
    cs->signals        = NULL;
    cs->signal_stds    = NULL;
    cs->standard_count = 0;
    cs->replicates     = 0;
}

void linear_calibration_init(LinearCalibration *lc)
{
    if (!lc) return;
    memset(lc, 0, sizeof(LinearCalibration));
    lc->slope     = 1.0;
    lc->intercept = 0.0;
}

void fourpl_fit_init(FourPLLogisticFit *fplf)
{
    if (!fplf) return;
    memset(fplf, 0, sizeof(FourPLLogisticFit));
    fplf->a = 0.01;  /* lower asymptote */
    fplf->b = 1.0;   /* slope */
    fplf->c = 1.0;   /* EC50 */
    fplf->d = 2.0;   /* upper asymptote */
}

void fivepl_fit_init(FivePLLogisticFit *fplf)
{
    if (!fplf) return;
    memset(fplf, 0, sizeof(FivePLLogisticFit));
    fplf->a = 0.01;
    fplf->b = 1.0;
    fplf->c = 1.0;
    fplf->d = 2.0;
    fplf->g = 1.0;   /* symmetric by default */
}

void sensor_array_init(SensorArray *sa, int32_t capacity)
{
    if (!sa) return;
    memset(sa, 0, sizeof(SensorArray));
    sa->capacity = (capacity > 0) ? capacity : 8;
    sa->sensors = (ChemiresistorElement *)calloc((size_t)sa->capacity,
                                                  sizeof(ChemiresistorElement));
    sa->sensor_count = 0;
    sa->pattern_dimension = sa->capacity;
    sa->pattern_vector = (double *)calloc((size_t)sa->pattern_dimension,
                                           sizeof(double));
    sa->sampling_rate_hz = 1.0;
}

void sensor_array_free(SensorArray *sa)
{
    if (!sa) return;
    free(sa->sensors);
    free(sa->pattern_vector);
    sa->sensors   = NULL;
    sa->pattern_vector = NULL;
    sa->capacity  = 0;
    sa->sensor_count = 0;
    sa->pattern_dimension = 0;
}
