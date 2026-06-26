/**
 * @file mems_basics.h
 * @brief L1 Definitions — MEMS sensor types, core structures, and fundamental parameters
 *
 * Covers:
 *   L1 - MEMS sensor types (accelerometer, gyroscope, magnetometer, pressure,
 *        microphone, resonator), sensitivity, full-scale range, bandwidth,
 *        resolution, noise density, bias instability, scale factor error,
 *        cross-axis sensitivity, non-orthogonality, misalignment
 *   L2 - MEMS transduction principles (capacitive, piezoresistive, piezoelectric,
 *        thermal, optical), fabrication overview, mechanical-thermal noise
 *
 * Reference: IEEE Std 1293-2018 (Inertial Sensor Terminology),
 *            IEEE Std 2700-2017 (MEMS Performance Parameters),
 *            NIST MEMS Measurement Handbook
 * Courses: MIT 6.630 EM Waves & Sensors, Stanford EE247 MEMS & Sensors,
 *          Berkeley EE117 EM/Sensors, Michigan EECS 411 Microwave/MEMS,
 *          ETH 227-0455 EM & Sensors, Tsinghua 传感器原理
 */

#ifndef MEMS_BASICS_H
#define MEMS_BASICS_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif
#ifndef M_SQRT2
#define M_SQRT2     1.41421356237309504880
#endif

/* ─── L1: MEMS Sensor Type Enumeration ─────────────────────────────────────── */
typedef enum {
    MEMS_TYPE_ACCELEROMETER = 0,
    MEMS_TYPE_GYROSCOPE     = 1,
    MEMS_TYPE_MAGNETOMETER  = 2,
    MEMS_TYPE_PRESSURE      = 3,
    MEMS_TYPE_MICROPHONE    = 4,
    MEMS_TYPE_RESONATOR     = 5,
    MEMS_TYPE_TEMPERATURE   = 6,
    MEMS_TYPE_HUMIDITY      = 7,
    MEMS_TYPE_GAS           = 8,
    MEMS_TYPE_FLOW          = 9,
    MEMS_TYPE_FORCE         = 10,
    MEMS_TYPE_TILT          = 11
} MemsSensorType;

const char *mems_sensor_type_name(MemsSensorType type);

/* ─── L1: MEMS Transduction Mechanism ──────────────────────────────────────── */
typedef enum {
    MEMS_TRANSDUCE_CAPACITIVE     = 0,
    MEMS_TRANSDUCE_PIEZORESISTIVE = 1,
    MEMS_TRANSDUCE_PIEZOELECTRIC  = 2,
    MEMS_TRANSDUCE_THERMAL        = 3,
    MEMS_TRANSDUCE_OPTICAL        = 4,
    MEMS_TRANSDUCE_TUNNELING      = 5,
    MEMS_TRANSDUCE_RESONANT       = 6,
    MEMS_TRANSDUCE_MAGNETIC       = 7
} MemsTransduceType;

const char *mems_transduce_type_name(MemsTransduceType type);

/* ─── L1: Sensor Axis Enumeration ──────────────────────────────────────────── */
typedef enum {
    MEMS_AXIS_X     = 0,
    MEMS_AXIS_Y     = 1,
    MEMS_AXIS_Z     = 2,
    MEMS_AXIS_ROLL  = 0,
    MEMS_AXIS_PITCH = 1,
    MEMS_AXIS_YAW   = 2
} MemsAxis;

/* ─── L1: Core 3-Axis Vector ──────────────────────────────────────────────────
 * All MEMS inertial and magnetic measurements are 3-axis vectors.
 * Units: accelerometer [m/s²] or [g], gyroscope [rad/s] or [°/s],
 *        magnetometer [μT] or [Gauss]
 */
typedef struct {
    double x, y, z;
} MemsVector3;

/* ─── L1: Quaternion for orientation representation ────────────────────────── */
typedef struct {
    double w, x, y, z;
} MemsQuaternion;

/* ─── L1: Rotation Matrix 3x3 ───────────────────────────────────────────────── */
typedef struct {
    double m[3][3];
} MemsMatrix33;

/* ─── L1: Core MEMS Sensor Specification ──────────────────────────────────────
 * Standard datasheet parameters per IEEE 2700-2017
 */
typedef struct {
    /* Identity */
    MemsSensorType   type;
    MemsTransduceType transduce;
    char             model[32];
    char             manufacturer[32];
    char             serial[32];

    /* Full-scale range (sensor units) */
    double           fsr;
    double           fsr_min;
    double           fsr_max;

    /* Sensitivity / Scale Factor */
    double           sensitivity;        /* LSB/unit or V/unit */
    double           scale_factor;       /* nominal SF */
    double           scale_factor_error; /* % of FS */
    double           scale_factor_tempco;/* ppm/°C */

    /* Bias */
    double           bias_initial;       /* initial bias (sensor units) */
    double           bias_instability;   /* in-run bias stability */
    double           bias_tempco;        /* ppm/°C or mg/°C */
    double           bias_vib_rect;      /* VRE coefficient */

    /* Noise */
    double           noise_density;      /* μg/√Hz or °/s/√Hz */
    double           random_walk;        /* velocity/angle random walk */
    double           rate_noise_density; /* for gyros: °/s/√Hz */

    /* Dynamic */
    double           bandwidth;          /* -3dB bandwidth [Hz] */
    double           resonant_freq;      /* mechanical resonance [Hz] */
    double           q_factor;           /* quality factor */
    double           group_delay;        /* [ms] */

    /* Axis alignment */
    double           cross_axis_sens;    /* %  */
    double           non_orthogonality;  /* degrees */
    double           misalignment;       /* degrees */
    double           g_sensitivity;      /* °/s/g for gyros */

    /* Environmental */
    double           temp_range_min;     /* °C */
    double           temp_range_max;     /* °C */
    double           supply_voltage;     /* V */
    double           power_consumption;  /* mW */
    double           shock_limit;        /* g */

    /* Output */
    uint8_t          resolution_bits;
    double           odr;                /* output data rate [Hz] */
    const char       *interface_type;    /* SPI, I2C, UART */

} MemsSpec;

/* ─── L1: MEMS Sensor Reading (single timestamped sample) ─────────────────── */
typedef struct {
    MemsSensorType type;
    MemsVector3    value;           /* 3-axis measurement */
    MemsVector3    raw;             /* raw ADC counts */
    MemsQuaternion orientation;     /* optional orientation */
    double         temperature;     /* on-chip temperature [°C] */
    double         timestamp;       /* seconds */
    double         dt;              /* time since last sample [s] */
    uint32_t       sequence;
    uint8_t        status;
} MemsSample;

/* ─── L1: MEMS IMU Configuration (combined accel+gyro) ────────────────────── */
typedef struct {
    MemsSpec accel;
    MemsSpec gyro;
    MemsSpec mag;           /* optional magnetometer */
    double   sample_rate;   /* Hz, common rate */
    uint8_t  accel_enabled;
    uint8_t  gyro_enabled;
    uint8_t  mag_enabled;
    uint8_t  temp_enabled;
} MemsImuConfig;

/* ─── L1: MEMS Fabrication Process ────────────────────────────────────────────
 * L2 concept: overview of MEMS fabrication for datasheet literacy
 */
typedef enum {
    MEMS_FAB_BULK_MICROMACHINED  = 0,
    MEMS_FAB_SURFACE_MICROMACH   = 1,
    MEMS_FAB_LIGA                = 2,
    MEMS_FAB_DRIE                = 3,
    MEMS_FAB_SOI                 = 4,
    MEMS_FAB_CMOS_MEMS           = 5
} MemsFabricationType;

const char *mems_fab_type_name(MemsFabricationType fab);

/* ─── L2: Seismic Mass Parameters ─────────────────────────────────────────────
 * The fundamental mechanical element in MEMS inertial sensors.
 * Mass-spring-damper system governed by Newton's 2nd law.
 */
typedef struct {
    double mass;            /* proof mass [kg] */
    double spring_k;        /* spring constant [N/m] */
    double damping_c;       /* damping coefficient [N·s/m] */
    double resonant_freq;   /* ω_n = sqrt(k/m) / (2π) [Hz] */
    double q_factor;        /* Q = sqrt(k·m) / c */
    double damping_ratio;   /* ζ = c / (2·sqrt(k·m)) */
    double displacement_max; /* maximum displacement [m] */
    double gap_capacitive;  /* capacitive gap [m] */
    double area_electrode;  /* electrode area [m²] */
    double permittivity;    /* ε [F/m], typically 8.854e-12 */
} MemsSeismicMass;

/**
 * @brief Compute resonant frequency from mass and spring constant
 *        f_n = (1/2π) · sqrt(k/m)
 *        Reference: Thomson, "Theory of Vibration" (1996)
 */
double mems_resonant_freq(double mass_kg, double spring_k);

/**
 * @brief Compute quality factor from mass, spring constant, damping
 *        Q = sqrt(k·m) / c
 */
double mems_quality_factor(double mass_kg, double spring_k, double damping_c);

/**
 * @brief Compute damping ratio from mass, spring constant, damping
 *        ζ = c / (2·sqrt(k·m))
 */
double mems_damping_ratio(double mass_kg, double spring_k, double damping_c);

/**
 * @brief Compute capacitive sensitivity from displacement
 *        ΔC = ε·A·Δx / d²  (parallel plate, small displacement approx)
 */
double mems_capacitive_sensitivity(double epsilon, double area,
                                    double gap, double displacement);

/* ─── L1: Spec Initialization Functions ──────────────────────────────────── */

/**
 * @brief Configure MEMS spec for typical consumer accelerometer
 *        (Bosch BMI270 / ST LSM6DSO class)
 */
void mems_spec_typical_accel(MemsSpec *spec);

/**
 * @brief Configure MEMS spec for typical consumer gyroscope
 */
void mems_spec_typical_gyro(MemsSpec *spec);

/**
 * @brief Configure MEMS spec for typical magnetometer
 *        (ST LIS3MDL / AKM AK09915 class)
 */
void mems_spec_typical_mag(MemsSpec *spec);

/* ─── L1: Vector3 Operations ─────────────────────────────────────────────── */

MemsVector3 mems_vector3(double x, double y, double z);
double mems_vector3_mag(MemsVector3 v);
MemsVector3 mems_vector3_add(MemsVector3 a, MemsVector3 b);
MemsVector3 mems_vector3_sub(MemsVector3 a, MemsVector3 b);
MemsVector3 mems_vector3_scale(MemsVector3 v, double s);
double mems_vector3_dot(MemsVector3 a, MemsVector3 b);
MemsVector3 mems_vector3_cross(MemsVector3 a, MemsVector3 b);
MemsVector3 mems_vector3_normalize(MemsVector3 v);

/* ─── L1: Quaternion Utilities ───────────────────────────────────────────── */

MemsQuaternion mems_quat_identity_impl(void);
double mems_quat_norm(MemsQuaternion q);

/* ─── L1: Unit Conversions ───────────────────────────────────────────────── */

double mems_deg_to_rad(double deg);
double mems_rad_to_deg(double rad);
double mems_dps_to_radps(double dps);
double mems_ms2_to_g(double ms2);
double mems_g_to_ms2(double g);

#endif /* MEMS_BASICS_H */
