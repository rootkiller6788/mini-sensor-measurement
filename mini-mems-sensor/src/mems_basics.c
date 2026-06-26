/**
 * @file mems_basics.c
 * @brief L1-L2 MEMS Sensor Fundamentals — Type definitions, sensor specifications,
 *        transduction mechanisms, basic MEMS mechanical parameters
 */

#include "mems_basics.h"
#include <math.h>
#include <string.h>
#include <float.h>

/* ─── L1: Sensor Type Name Lookup ──────────────────────────────────────────── */

const char *mems_sensor_type_name(MemsSensorType type) {
    static const char *names[] = {
        [MEMS_TYPE_ACCELEROMETER] = "Accelerometer",
        [MEMS_TYPE_GYROSCOPE]     = "Gyroscope",
        [MEMS_TYPE_MAGNETOMETER]  = "Magnetometer",
        [MEMS_TYPE_PRESSURE]      = "Pressure Sensor",
        [MEMS_TYPE_MICROPHONE]    = "Microphone (MEMS)",
        [MEMS_TYPE_RESONATOR]     = "Resonator/Oscillator",
        [MEMS_TYPE_TEMPERATURE]   = "Temperature Sensor",
        [MEMS_TYPE_HUMIDITY]      = "Humidity Sensor",
        [MEMS_TYPE_GAS]           = "Gas Sensor",
        [MEMS_TYPE_FLOW]          = "Flow Sensor",
        [MEMS_TYPE_FORCE]         = "Force Sensor",
        [MEMS_TYPE_TILT]          = "Tilt/Inclinometer"
    };
    if (type <= MEMS_TYPE_TILT) return names[type];
    return "Unknown";
}

const char *mems_transduce_type_name(MemsTransduceType type) {
    static const char *names[] = {
        [MEMS_TRANSDUCE_CAPACITIVE]     = "Capacitive",
        [MEMS_TRANSDUCE_PIEZORESISTIVE] = "Piezoresistive",
        [MEMS_TRANSDUCE_PIEZOELECTRIC]  = "Piezoelectric",
        [MEMS_TRANSDUCE_THERMAL]        = "Thermal",
        [MEMS_TRANSDUCE_OPTICAL]        = "Optical",
        [MEMS_TRANSDUCE_TUNNELING]      = "Electron Tunneling",
        [MEMS_TRANSDUCE_RESONANT]       = "Resonant",
        [MEMS_TRANSDUCE_MAGNETIC]       = "Magnetic"
    };
    if (type <= MEMS_TRANSDUCE_MAGNETIC) return names[type];
    return "Unknown";
}

const char *mems_fab_type_name(MemsFabricationType fab) {
    static const char *names[] = {
        [MEMS_FAB_BULK_MICROMACHINED]  = "Bulk Micromachined",
        [MEMS_FAB_SURFACE_MICROMACH]   = "Surface Micromachined",
        [MEMS_FAB_LIGA]                = "LIGA",
        [MEMS_FAB_DRIE]                = "DRIE (Deep Reactive Ion Etch)",
        [MEMS_FAB_SOI]                 = "SOI (Silicon-On-Insulator)",
        [MEMS_FAB_CMOS_MEMS]           = "CMOS-MEMS"
    };
    if (fab <= MEMS_FAB_CMOS_MEMS) return names[fab];
    return "Unknown";
}

/* ─── L2: Seismic Mass Analysis ──────────────────────────────────────────────
 * The fundamental building block of MEMS inertial sensors.
 * A proof mass suspended by springs, whose displacement under acceleration
 * is measured (typically capacitively).
 *
 * Governing equation (Newton's 2nd law):
 *   m · ẍ + c · ẋ + k · x = m · a_ext
 *
 * Laplace domain transfer function:
 *   H(s) = X(s) / A_ext(s) = 1 / (s² + (c/m)·s + k/m)
 *        = 1 / (s² + 2·ζ·ω_n·s + ω_n²)
 *
 * where ω_n = sqrt(k/m), ζ = c / (2·sqrt(k·m))
 */

double mems_resonant_freq(double mass_kg, double spring_k) {
    if (mass_kg <= 0.0 || spring_k <= 0.0) return 0.0;
    return sqrt(spring_k / mass_kg) / (2.0 * M_PI);
}

double mems_quality_factor(double mass_kg, double spring_k, double damping_c) {
    if (mass_kg <= 0.0 || spring_k <= 0.0 || damping_c <= 0.0) return 0.0;
    return sqrt(spring_k * mass_kg) / damping_c;
}

double mems_damping_ratio(double mass_kg, double spring_k, double damping_c) {
    if (mass_kg <= 0.0 || spring_k <= 0.0) return 0.0;
    if (damping_c <= 0.0) return 0.0;
    return damping_c / (2.0 * sqrt(spring_k * mass_kg));
}

double mems_capacitive_sensitivity(double epsilon, double area,
                                    double gap, double displacement) {
    if (gap <= 0.0 || area <= 0.0) return 0.0;
    /* For parallel plate: C = ε·A/d, ΔC ≈ ε·A·Δx/d² (small displacement) */
    return epsilon * area * displacement / (gap * gap);
}

/* ─── L1: Sensor Spec Initialization ──────────────────────────────────────── */

static void mems_spec_init_defaults(MemsSpec *spec) {
    if (!spec) return;
    memset(spec, 0, sizeof(*spec));
    spec->q_factor = 1.0;
    spec->supply_voltage = 3.3;
    spec->resolution_bits = 16;
    spec->interface_type = "SPI";
}

/**
 * @brief Configure MEMS spec for a typical consumer accelerometer
 *        Modeled on Bosch BMI270 / ST LSM6DSO performance class.
 *        Reference: Bosch BMI270 Datasheet (2020), ST LSM6DSO Datasheet (2019)
 */
void mems_spec_typical_accel(MemsSpec *spec) {
    mems_spec_init_defaults(spec);
    spec->type = MEMS_TYPE_ACCELEROMETER;
    spec->transduce = MEMS_TRANSDUCE_CAPACITIVE;
    spec->fsr = 156.96;            /* m/s² (16g) */
    spec->fsr_min = 19.62;         /* m/s² (2g) */
    spec->fsr_max = 156.96;        /* m/s² (16g) */
    spec->sensitivity = 16384.0;   /* LSB/g at ±2g */
    spec->scale_factor = 1.0;
    spec->scale_factor_error = 0.5;      /* ±0.5% */
    spec->scale_factor_tempco = 0.01;    /* 0.01%/°C */
    spec->noise_density = 120.0;         /* μg/√Hz */
    spec->random_walk = 0.02;            /* m/s/√h velocity random walk */
    spec->bias_initial = 0.3924;         /* 40mg initial bias */
    spec->bias_instability = 0.03924;    /* 4mg in-run bias stability */
    spec->bias_tempco = 0.2;             /* mg/°C */
    spec->bandwidth = 200.0;             /* Hz */
    spec->resonant_freq = 2400.0;        /* Hz */
    spec->q_factor = 200.0;
    spec->cross_axis_sens = 1.0;         /* % */
    spec->non_orthogonality = 1.0;       /* degrees */
    spec->odr = 1600.0;                  /* Hz */
}

/**
 * @brief Configure for typical consumer gyroscope
 *        Modeled on Bosch BMI270 / ST LSM6DSO performance class.
 */
void mems_spec_typical_gyro(MemsSpec *spec) {
    mems_spec_init_defaults(spec);
    spec->type = MEMS_TYPE_GYROSCOPE;
    spec->transduce = MEMS_TRANSDUCE_CAPACITIVE;
    spec->fsr = 34.9066;           /* rad/s (2000°/s) */
    spec->fsr_max = 34.9066;
    spec->sensitivity = 262.0;     /* LSB/°/s at ±125°/s */
    spec->scale_factor = 1.0;
    spec->scale_factor_error = 0.5;
    spec->scale_factor_tempco = 0.03;    /* %/°C */
    spec->noise_density = 0.007;         /* °/s/√Hz = 0.000122 rad/s/√Hz */
    spec->random_walk = 0.002;           /* °/√h angle random walk */
    spec->bias_initial = 0.05236;        /* 3°/s initial bias */
    spec->bias_instability = 0.017453;   /* 1°/s in-run bias stability */
    spec->bias_tempco = 0.015;           /* °/s/°C */
    spec->bandwidth = 200.0;
    spec->resonant_freq = 15000.0;       /* Hz, typical drive mode */
    spec->q_factor = 500.0;
    spec->cross_axis_sens = 1.0;
    spec->non_orthogonality = 1.0;
    spec->g_sensitivity = 0.1;           /* °/s/g */
    spec->odr = 1600.0;
}

/**
 * @brief Configure for typical magnetometer
 *        Modeled on ST LIS3MDL / AKM AK09915 performance class.
 */
void mems_spec_typical_mag(MemsSpec *spec) {
    mems_spec_init_defaults(spec);
    spec->type = MEMS_TYPE_MAGNETOMETER;
    spec->transduce = MEMS_TRANSDUCE_MAGNETIC;
    spec->fsr = 1600.0;            /* μT (±16 Gauss) */
    spec->fsr_max = 1600.0;
    spec->sensitivity = 0.15;      /* μT/LSB */
    spec->scale_factor = 1.0;
    spec->scale_factor_error = 0.5;
    spec->noise_density = 0.003;   /* μT/√Hz */
    spec->bias_initial = 5.0;      /* μT */
    spec->bias_tempco = 0.01;      /* μT/°C */
    spec->bandwidth = 100.0;
    spec->cross_axis_sens = 2.0;
    spec->non_orthogonality = 2.0;
    spec->odr = 1000.0;
}

/* ─── L1: Utility Functions ───────────────────────────────────────────────── */

/**
 * @brief Create a 3D vector
 */
MemsVector3 mems_vector3(double x, double y, double z) {
    MemsVector3 v = {x, y, z};
    return v;
}

/**
 * @brief Vector magnitude
 */
double mems_vector3_mag(MemsVector3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/**
 * @brief Vector addition
 */
MemsVector3 mems_vector3_add(MemsVector3 a, MemsVector3 b) {
    MemsVector3 r = {a.x + b.x, a.y + b.y, a.z + b.z};
    return r;
}

/**
 * @brief Vector subtraction
 */
MemsVector3 mems_vector3_sub(MemsVector3 a, MemsVector3 b) {
    MemsVector3 r = {a.x - b.x, a.y - b.y, a.z - b.z};
    return r;
}

/**
 * @brief Vector scalar multiplication
 */
MemsVector3 mems_vector3_scale(MemsVector3 v, double s) {
    MemsVector3 r = {v.x * s, v.y * s, v.z * s};
    return r;
}

/**
 * @brief Vector dot product
 */
double mems_vector3_dot(MemsVector3 a, MemsVector3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/**
 * @brief Vector cross product
 */
MemsVector3 mems_vector3_cross(MemsVector3 a, MemsVector3 b) {
    MemsVector3 r = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
    return r;
}

/**
 * @brief Vector normalization
 */
MemsVector3 mems_vector3_normalize(MemsVector3 v) {
    double mag = mems_vector3_mag(v);
    if (mag < DBL_EPSILON) {
        MemsVector3 zero = {0.0, 0.0, 0.0};
        return zero;
    }
    return mems_vector3_scale(v, 1.0 / mag);
}

/**
 * @brief Initialize identity quaternion
 */
MemsQuaternion mems_quat_identity_impl(void) {
    MemsQuaternion q = {1.0, 0.0, 0.0, 0.0};
    return q;
}

/**
 * @brief Compute quaternion norm
 */
double mems_quat_norm(MemsQuaternion q) {
    return sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
}

/**
 * @brief Convert degrees to radians
 */
double mems_deg_to_rad(double deg) {
    return deg * M_PI / 180.0;
}

/**
 * @brief Convert radians to degrees
 */
double mems_rad_to_deg(double rad) {
    return rad * 180.0 / M_PI;
}

/**
 * @brief Convert angular rate units: °/s to rad/s
 */
double mems_dps_to_radps(double dps) {
    return dps * M_PI / 180.0;
}

/**
 * @brief Convert m/s² to g
 */
double mems_ms2_to_g(double ms2) {
    return ms2 / 9.80665;
}

/**
 * @brief Convert g to m/s²
 */
double mems_g_to_ms2(double g) {
    return g * 9.80665;
}
