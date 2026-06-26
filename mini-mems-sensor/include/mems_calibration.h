/**
 * @file mems_calibration.h
 * @brief L5 MEMS Calibration Methods — 6-position static, rate table,
 *        temperature compensation, axis alignment, scale factor correction
 *
 * Covers:
 *   L5 - Six-position calibration (accelerometer), rate table calibration (gyro),
 *        multi-position least-squares, temperature polynomial compensation,
 *        cross-axis sensitivity matrix, soft-iron / hard-iron correction
 *
 * Reference: IEEE Std 1554-2005 (Inertial Sensor Test Equipment)
 *            IEEE Std 2700-2017 (MEMS Performance Parameters)
 *            Titterton & Weston, "Strapdown Inertial Navigation Technology" (2004)
 *            Groves, "Principles of GNSS, Inertial, and Multisensor Navigation" (2013)
 * Courses: MIT 16.485 (GNSS/Inertial), Stanford EE267 (Virtual Reality),
 *          Berkeley EE123 (DSP), ETH 227-0436 (Communication/Estimation)
 */

#ifndef MEMS_CALIBRATION_H
#define MEMS_CALIBRATION_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include "mems_basics.h"

/* ─── L5: 6-Position Calibration Data ────────────────────────────────────────
 * Standard +X, -X, +Y, -Y, +Z, -Z orientations for accelerometer calibration.
 * Each position records the measured 3-axis accelerometer reading.
 */
typedef struct {
    MemsVector3 measured[6];  /* readings in each of 6 positions */
    MemsVector3 expected[6];  /* expected gravity vector */
    double      temperature[6]; /* temperature during each measurement */
    size_t      num_positions;  /* typically 6 */
} MemsSixPosData;

/* ─── L5: Accelerometer Calibration Parameters ────────────────────────────────
 * Classic model: a_calibrated = S · (a_raw - b)
 * where S is 3x3 scale+misalignment matrix and b is bias vector.
 */
typedef struct {
    MemsVector3 bias;           /* x, y, z bias [m/s²] or [g] */
    MemsMatrix33 scale_mis;     /* 3x3: diagonal = scale, off-diag = misalignment */
    MemsVector3 scale_factor;   /* per-axis scale factor */
    double      residuals[6];   /* per-position residual error */
    double      rms_error;      /* overall RMS calibration error */
} MemsAccelCalib;

/* ─── L5: Gyroscope Calibration Parameters ───────────────────────────────── */
typedef struct {
    MemsVector3 bias;           /* x, y, z rate bias [rad/s] or [°/s] */
    MemsMatrix33 scale_mis;     /* 3x3 scale+misalignment matrix */
    MemsVector3 scale_factor;   /* per-axis scale factor */
    MemsVector3 g_sensitivity;  /* g-sensitivity [°/s per g] */
    double      residuals[12];  /* residuals for multi-rate measurements */
    double      rms_error;
} MemsGyroCalib;

/* ─── L5: Magnetometer Calibration (Hard-Iron + Soft-Iron) ───────────────── */
typedef struct {
    MemsVector3 hard_iron;      /* additive offset from ferrous materials */
    MemsMatrix33 soft_iron;     /* 3x3 distortion from nearby metals */
    MemsVector3 scale_factor;
    double      field_magnitude;/* expected Earth field [μT or Gauss] */
    double      rms_error;
} MemsMagCalib;

/* ─── L5: Temperature Compensation Model ──────────────────────────────────── */

/**
 * @brief Temperature compensation polynomial coefficients
 *        Compensated = Raw / (1 + a₁·ΔT + a₂·ΔT² + a₃·ΔT³)
 *        where ΔT = T - T_ref
 */
typedef struct {
    double ref_temp;        /* T_ref [°C] */
    double a1;              /* linear tempco [ppm/°C] or [1/°C] */
    double a2;              /* quadratic tempco [ppm/°C²] */
    double a3;              /* cubic tempco [ppm/°C³] */
    double fit_r_squared;   /* goodness of fit */
} MemsTempComp;

/* ─── L5: Six-Position Accelerometer Calibration ───────────────────────────── */

/**
 * @brief Solve for bias and scale+misalignment matrix using 6-position data
 *        Minimizes ||a_expected - S·(a_measured - b)||²
 *        Reference: Titterton & Weston (2004) §11.3, Groves (2013) §4.3
 *
 * @param data   6-position measurement data
 * @param calib  Output calibration parameters
 * @return RMS error of fit, < 0 on error
 */
double mems_calib_accel_6pos(const MemsSixPosData *data, MemsAccelCalib *calib);

/**
 * @brief Apply accelerometer calibration: a_out = S · (a_in - bias)
 */
MemsVector3 mems_apply_accel_calib(MemsVector3 raw, const MemsAccelCalib *calib);

/* ─── L5: Gyroscope Rate Table Calibration ─────────────────────────────────── */

/**
 * @brief Gyro rate table calibration using multiple rotation rates
 *        Minimizes ||ω_known - S·(ω_meas - b - G·a)||²
 *        where G is g-sensitivity matrix.
 *        Reference: IEEE 2700-2017 §6.2.3
 *
 * @param rates        Array of known rotation rates (≤12 recommended)
 * @param measured     Array of corresponding gyro measurements
 * @param accelerations Array of applied accelerations (for g-sensitivity)
 * @param n            Number of rate points
 * @param calib        Output calibration parameters
 * @return RMS error of fit
 */
double mems_calib_gyro_rate(const double rates[][3], const double measured[][3],
                             const double accelerations[][3], size_t n,
                             MemsGyroCalib *calib);

/**
 * @brief Apply gyro calibration: ω_out = S · (ω_in - bias)
 */
MemsVector3 mems_apply_gyro_calib(MemsVector3 raw, const MemsGyroCalib *calib);

/* ─── L5: Magnetometer Hard/Soft Iron Calibration ───────────────────────── */

/**
 * @brief Ellipsoid fitting for magnetometer calibration
 *        Using iterative least-squares to fit:
 *        (x - hx)²/a² + (y - hy)²/b² + (z - hz)²/c² = R²
 *        Hard iron = (hx, hy, hz), Soft iron = from eigen-decomposition
 *        Reference: Gebre-Egziabher et al., "Calibration of Strapdown
 *        Magnetometers" (2006)
 *
 * @param measurements Array of 3D magnetic measurements (≥ 9 points)
 * @param n            Number of measurements
 * @param expected_mag Expected field magnitude [μT]
 * @param calib        Output calibration parameters
 * @return RMS error
 */
double mems_calib_mag_ellipsoid(const double measurements[][3], size_t n,
                                 double expected_mag, MemsMagCalib *calib);

/**
 * @brief Apply magnetometer calibration: m_out = S · (m_in - hard_iron)
 */
MemsVector3 mems_apply_mag_calib(MemsVector3 raw, const MemsMagCalib *calib);

/* ─── L5: Least-Squares Axis Alignment ─────────────────────────────────────── */

/**
 * @brief Compute cross-axis sensitivity matrix
 *        S = (A_meas · A_ref^T) · (A_ref · A_ref^T)⁻¹
 *        Reference: Gebre-Egziabher (2006)
 *
 * @param measured  N×3 matrix of measured vectors
 * @param reference N×3 matrix of reference vectors
 * @param n         Number of vector pairs
 * @param S         Output 3×3 alignment matrix
 * @return Condition number of A_ref·A_ref^T (quality indicator), < 0 on error
 */
double mems_alignment_matrix(const double measured[][3],
                              const double reference[][3],
                              size_t n, MemsMatrix33 *S);

/**
 * @brief Apply alignment matrix to a 3D vector
 */
MemsVector3 mems_apply_matrix33(MemsVector3 v, const MemsMatrix33 *M);

/* ─── L5: Temperature Compensation ──────────────────────────────────────── */

/**
 * @brief Fit polynomial temperature compensation model
 *        Uses ordinary least-squares on (T-T_ref) powers.
 *        Reference: IEEE 2700-2017 §6.4 (Temperature sensitivity)
 *
 * @param temperatures Array of temperature measurements [°C]
 * @param errors       Array of corresponding sensor errors at each temp
 * @param n            Number of (T, error) pairs
 * @param order        Polynomial order (1 to 4)
 * @param comp         Output temperature compensation model
 * @return R² goodness of fit, < 0 on error
 */
double mems_fit_temp_comp(const double *temperatures, const double *errors,
                           size_t n, int order, MemsTempComp *comp);

/**
 * @brief Apply temperature compensation
 *        out = raw / (1 + a1·ΔT + a2·ΔT² + a3·ΔT³)
 */
double mems_apply_temp_comp(double raw_value, double temperature,
                             const MemsTempComp *comp);

/**
 * @brief Compute temperature delta ΔT = T - T_ref
 */
double mems_temp_delta(double temperature, const MemsTempComp *comp);

#endif /* MEMS_CALIBRATION_H */
