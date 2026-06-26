/**
 * biosensor_types.h — Core type definitions for biosensor and chemical sensor systems
 *
 * L1 Definitions: All fundamental biosensor data structures.
 *
 * References:
 *   - IUPAC definition of biosensor (Thévenot et al., 1999)
 *   - Michaelis & Menten (1913) enzyme kinetics framework
 *   - Clark & Lyons (1962) first enzyme electrode
 *   - ISFET theory: Bergveld (1970)
 *
 * L9 alignment: MIT 6.555 (Bioelectronics), Stanford EE 293 (Biosensors),
 *   Berkeley BioE 121 (Bioinstrumentation)
 */

#ifndef BIOSENSOR_TYPES_H
#define BIOSENSOR_TYPES_H

#include <stdint.h>
#include <stdbool.h>

/* ============================================================================
 * L1: Sensor classification enums
 * ============================================================================ */

/** Biological recognition element types */
typedef enum {
    BIORECEPTOR_ENZYME       = 0,  /* Glucose oxidase, urease, etc. */
    BIORECEPTOR_ANTIBODY     = 1,  /* IgG, IgM, monoclonal, polyclonal */
    BIORECEPTOR_DNA          = 2,  /* ssDNA probes, aptamers */
    BIORECEPTOR_CELL         = 3,  /* Whole-cell microbial biosensors */
    BIORECEPTOR_TISSUE       = 4,  /* Plant/animal tissue slices */
    BIORECEPTOR_MIP          = 5,  /* Molecularly Imprinted Polymer (synthetic) */
    BIORECEPTOR_LECTIN       = 6,  /* Carbohydrate-binding proteins */
    BIORECEPTOR_BACTERIOPHAGE = 7, /* Phage-based recognition */
    BIORECEPTOR_NANOPORE     = 8   /* Alpha-hemolysin, solid-state */
} BioreceptorType;

/** Transduction mechanism types */
typedef enum {
    TRANSDUCER_AMPEROMETRIC   = 0,  /* Current measurement at fixed voltage */
    TRANSDUCER_POTENTIOMETRIC = 1,  /* Voltage measurement at zero current */
    TRANSDUCER_IMPEDIMETRIC   = 2,  /* AC impedance spectroscopy */
    TRANSDUCER_CONDUCTOMETRIC = 3,  /* Conductance change measurement */
    TRANSDUCER_OPTICAL        = 4,  /* Fluorescence, SPR, absorbance */
    TRANSDUCER_PIEZOELECTRIC  = 5,  /* QCM, SAW mass-sensitive */
    TRANSDUCER_CALORIMETRIC   = 6,  /* Thermistor-based enthalpy sensing */
    TRANSDUCER_ISFET          = 7,  /* Ion-Sensitive Field-Effect Transistor */
    TRANSDUCER_NANOWIRE       = 8   /* SiNW, CNT FET sensors */
} TransducerType;

/** Measurement modes */
typedef enum {
    MODE_END_POINT    = 0,  /* Single reading after incubation */
    MODE_KINETIC      = 1,  /* Continuous rate measurement */
    MODE_DIFFERENTIAL = 2,  /* Signal difference from baseline */
    MODE_RATIOMETRIC  = 3,  /* Ratio of two wavelengths/electrodes */
    MODE_PULSED       = 4,  /* Pulsed amperometric detection (PAD) */
    MODE_CYCLIC       = 5   /* Cyclic voltammetry */
} MeasurementMode;

/* ============================================================================
 * L1: Core biosensor descriptor
 * ============================================================================ */

/** A complete biosensor definition combining bioreceptor + transducer */
typedef struct {
    BioreceptorType bioreceptor_type;
    TransducerType   transducer_type;
    MeasurementMode  measurement_mode;

    const char *analyte_name;       /* e.g., "glucose", "cortisol", "Hg2+" */
    const char *bioreceptor_name;   /* e.g., "Glucose Oxidase", "anti-cortisol IgG" */
    const char *sample_matrix;      /* e.g., "blood", "urine", "river water" */

    double linear_range_min;        /* Lower limit of linear range [analyte units] */
    double linear_range_max;        /* Upper limit of linear range [analyte units] */

    double sensitivity;             /* Slope of calibration curve [signal/conc] */
    double limit_of_detection;      /* LOD = 3.3 * sigma_blank / sensitivity */
    double limit_of_quantification; /* LOQ = 10 * sigma_blank / sensitivity */

    double response_time_sec;       /* t90 time to 90% of final signal */
    double recovery_time_sec;       /* Time to return to baseline */

    double selectivity_coefficient; /* K_ij per IUPAC: response ratio */
    double ph_optimum;              /* pH for maximum activity */
    double temp_optimum_celsius;    /* Temperature for maximum activity */

    int32_t shelf_life_days;        /* Storage stability */
    int32_t reuse_cycles;           /* Number of measurements per sensor */
} BiosensorDescriptor;

/* ============================================================================
 * L1: Kinetic model parameters
 * ============================================================================ */

/** Michaelis-Menten enzyme kinetics parameters */
typedef struct {
    double km;                      /* Michaelis constant [M or mM] */
    double vmax;                    /* Maximum reaction rate [signal/s] */
    double kcat;                    /* Catalytic constant [1/s] = Vmax/[E0] */
    double substrate_conc;          /* Current substrate concentration for calc */
    double enzyme_conc;             /* Active enzyme concentration [M] */
} MichaelisMentenParams;

/** Hill equation parameters (cooperative binding) */
typedef struct {
    double half_max_conc;           /* K_d or K_0.5 [concentration units] */
    double hill_coefficient;        /* n_H: > 1 positive cooperativity,
                                     *      < 1 negative cooperativity,
                                     *      = 1 independent binding */
    double max_response;            /* R_max asymptotic response */
    double baseline;                /* R_0 baseline response at zero concentration */
} HillParams;

/** Langmuir adsorption isotherm */
typedef struct {
    double k_adsorption;            /* Adsorption rate constant */
    double k_desorption;            /* Desorption rate constant */
    double q_max;                   /* Maximum adsorption capacity */
    double equilibrium_constant;    /* K = k_ads / k_des */
} LangmuirParams;

/* ============================================================================
 * L1: Electrochemical sensor parameters
 * ============================================================================ */

/** Nernst equation parameters */
typedef struct {
    double standard_potential;      /* E° [V] */
    double electron_count;          /* n: number of electrons transferred */
    double temperature_kelvin;      /* T [K] */
    double faraday_constant;        /* F = 96485.3329 C/mol (default) */
    double gas_constant;            /* R = 8.314462618 J/(mol·K) (default) */
} NernstParams;

/** Butler-Volmer kinetics parameters */
typedef struct {
    double exchange_current_density; /* i0 [A/m²] */
    double anodic_transfer_coeff;    /* α_a or β */
    double cathodic_transfer_coeff;  /* α_c */
    double electrode_area;           /* A [m²] */
    double overpotential;            /* η = E - E_eq [V] */
} ButlerVolmerParams;

/** Electrochemical impedance spectroscopy (EIS) parameters */
typedef struct {
    double solution_resistance;     /* R_s [Ω] — ohmic resistance of electrolyte */
    double charge_transfer_resistance; /* R_ct [Ω] — Faradaic resistance */
    double double_layer_capacitance;   /* C_dl [F] — electrode-electrolyte interface */
    double warburg_coefficient;     /* σ [Ω/√s] — diffusion impedance coefficient */
    double constant_phase_exponent;  /* n: 0.5 (Warburg) to 1.0 (ideal capacitor) */
    double frequency_min_hz;        /* Lower bound of measurement */
    double frequency_max_hz;        /* Upper bound of measurement */
} EISParams;

/** ISFET (Ion-Sensitive Field-Effect Transistor) parameters */
typedef struct {
    double threshold_voltage;       /* V_th [V] at pH 7 */
    double sensitivity_mv_per_ph;   /* Sub-Nernstian sensitivity [mV/pH] */
    double transconductance;        /* g_m [μS] */
    double drain_current;           /* I_DS [μA] in saturation */
    double reference_voltage;       /* V_ref [V] applied to reference electrode */
} ISFETParams;

/* ============================================================================
 * L1: Optical biosensor parameters
 * ============================================================================ */

/** Surface Plasmon Resonance (SPR) parameters */
typedef struct {
    double resonance_angle_deg;     /* θ_SPR [degrees] */
    double refractive_index_sensitivity; /* [RIU^-1] angular shift per RIU change */
    double penetration_depth_nm;    /* Evanescent field decay length [nm] */
    double wavelength_nm;           /* Incident light wavelength [nm] */
    double gold_layer_thickness_nm; /* Au film thickness [nm] */
    double prism_ri;                /* Prism refractive index at λ */
} SPRParams;

/** Fluorescence-based biosensor parameters */
typedef struct {
    double excitation_wavelength_nm;   /* λ_ex */
    double emission_wavelength_nm;     /* λ_em */
    double quantum_yield;              /* Φ: 0 to 1 */
    double extinction_coefficient;     /* ε [M⁻¹·cm⁻¹] */
    double stokes_shift_nm;            /* λ_em - λ_ex */
    double photobleaching_halflife_sec; /* t_1/2 under illumination */
    double forster_radius_nm;          /* R_0 for FRET [nm] at 50% efficiency */
} FluorescenceParams;

/** Fiber-optic biosensor parameters */
typedef struct {
    double fiber_core_diameter_um;     /* Core diameter [μm] */
    double numerical_aperture;         /* NA = √(n_core² - n_clad²) */
    double evanescent_decay_length_nm; /* Penetration depth into sensing region */
    double sensing_region_length_mm;   /* Length of exposed core */
    double coupling_efficiency;        /* η: fraction of light coupled into fiber */
    double propagation_loss_db_per_cm; /* Attenuation in sensing region */
} FiberOpticBiosensorParams;

/* ============================================================================
 * L1: Piezoelectric mass-sensitive sensor parameters
 * ============================================================================ */

/** Quartz Crystal Microbalance (QCM) parameters */
typedef struct {
    double fundamental_frequency_hz;    /* f_0 [Hz], typically 5-10 MHz */
    double electrode_area_cm2;         /* Overlapping electrode area [cm²] */
    double density_quartz_g_per_cm3;   /* ρ_q = 2.648 g/cm³ for AT-cut quartz */
    double shear_modulus_gpa;          /* μ_q = 2.947×10^11 g/(cm·s²) for AT-cut */
    double mass_sensitivity_hz_per_ng; /* Sauerbrey sensitivity [Hz/ng] */
    double quality_factor;             /* Q-factor of crystal in liquid */
    double motional_resistance_ohm;    /* R_m from Butterworth-Van Dyke model */
} QCMParameters;

/* ============================================================================
 * L1: Measurement result container
 * ============================================================================ */

/** Single measurement data point */
typedef struct {
    double timestamp_sec;           /* Time of measurement */
    double signal_value;            /* Raw signal (current, voltage, absorbance...) */
    double analyte_concentration;   /* Calculated concentration [analyte units] */
    double signal_to_noise_ratio;   /* SNR = mean_signal / std_noise */
    bool   is_outlier;              /* Grubbs test flag */
} MeasurementPoint;

/** Complete measurement run */
typedef struct {
    MeasurementPoint *points;       /* Array of measurement points */
    int32_t           point_count;
    int32_t           capacity;
    double            temperature_celsius;   /* Measured temperature */
    double            ph_value;              /* Measured pH */
    double            background_signal;     /* Blank/buffer response */
    double            noise_std;             /* Standard deviation of blank */
    const char       *timestamp_iso8601;     /* ISO 8601 timestamp */
} MeasurementSeries;

/* ============================================================================
 * L1: Calibration model structures
 * ============================================================================ */

/** Calibration standards (known concentrations with measured signals) */
typedef struct {
    double *concentrations;         /* Known analyte concentrations */
    double *signals;                /* Corresponding measured signals */
    double *signal_stds;            /* Standard deviation at each concentration */
    int32_t standard_count;         /* Number of calibration standards */
    int32_t replicates;             /* Replicates per standard */
} CalibrationStandards;

/** 4-Parameter Logistic (4PL) fit result */
typedef struct {
    double a;                       /* Lower asymptote (baseline) */
    double b;                       /* Slope factor (Hill coefficient) */
    double c;                       /* Inflection point (IC50 / EC50) */
    double d;                       /* Upper asymptote (max signal) */
    double r_squared;               /* Goodness of fit R² */
    double rmse;                    /* Root Mean Square Error */
} FourPLLogisticFit;

/** 5-Parameter Logistic (5PL) fit result (adds asymmetry) */
typedef struct {
    double a;                       /* Lower asymptote */
    double b;                       /* Slope factor */
    double c;                       /* Inflection point */
    double d;                       /* Upper asymptote */
    double g;                       /* Asymmetry factor (g=1 → symmetric 4PL) */
    double r_squared;
    double rmse;
} FivePLLogisticFit;

/** Linear regression calibration */
typedef struct {
    double slope;                   /* Sensitivity [signal/concentration] */
    double intercept;               /* Background signal */
    double r_squared;               /* R² goodness of fit */
    double standard_error_slope;    /* SE of slope estimate */
    double standard_error_intercept; /* SE of intercept estimate */
    double prediction_interval_95;  /* 95% prediction interval half-width */
} LinearCalibration;

/* ============================================================================
 * L1: Sensor array / electronic nose/tongue
 * ============================================================================ */

/** Generic sensor element in an array */
typedef struct {
    int32_t           sensor_id;
    const char       *sensor_material;     /* e.g., "polyaniline", "SnO2", "Pd-doped CNT" */
    double            baseline_resistance_ohm; /* R_0 in clean air/reference */
    double            sensitivity_ppm;     /* ΔR/R_0 per ppm of target gas */
    double            selectivity_vector[8]; /* Response ratios to 8 interferents */
    double            operating_temp_celsius;
    double            heater_power_watt;
    bool              is_preheated;
} ChemiresistorElement;

/** Electronic nose / sensor array */
typedef struct {
    ChemiresistorElement *sensors;
    int32_t               sensor_count;
    int32_t               capacity;
    const char           *array_name;
    double                sampling_rate_hz;
    double               *pattern_vector;   /* Normalized response pattern */
    int32_t               pattern_dimension;
} SensorArray;

/* ============================================================================
 * Constructor/Destructor prototypes (implemented in src/biosensor_types.c)
 * ============================================================================ */

void biosensor_descriptor_init(BiosensorDescriptor *desc);
void michaelis_menten_params_init(MichaelisMentenParams *mm, double km, double vmax);
void hill_params_init(HillParams *hp, double k_half, double n_h, double rmax);
void nernst_params_init(NernstParams *np, double e0, double n);
void butler_volmer_params_init(ButlerVolmerParams *bv, double i0, double alpha_a);
void eis_params_init(EISParams *eis);
void isfet_params_init(ISFETParams *isfet);
void spr_params_init(SPRParams *spr, double wavelength_nm, double ri_sens);
void fluorescence_params_init(FluorescenceParams *fp, double ex_nm, double em_nm);
void qcm_params_init(QCMParameters *qcm, double f0_hz, double area_cm2);
void measurement_series_init(MeasurementSeries *series, int32_t capacity);
void measurement_series_free(MeasurementSeries *series);
void calibration_standards_init(CalibrationStandards *cs, int32_t n_standards, int32_t reps);
void calibration_standards_free(CalibrationStandards *cs);
void linear_calibration_init(LinearCalibration *lc);
void fourpl_fit_init(FourPLLogisticFit *fplf);
void fivepl_fit_init(FivePLLogisticFit *fplf);
void sensor_array_init(SensorArray *sa, int32_t capacity);
void sensor_array_free(SensorArray *sa);

#endif /* BIOSENSOR_TYPES_H */
