/**
 * @file mems_calibration.c
 * @brief L5 MEMS Calibration — 6-position accelerometer, rate table gyro,
 *        magnetometer ellipsoid fitting, temperature compensation
 *
 * Reference: Titterton & Weston (2004) "Strapdown Inertial Navigation Technology"
 *            Groves (2013) "Principles of GNSS, Inertial, and Multisensor Navigation"
 *            IEEE Std 1554-2005, IEEE Std 2700-2017
 *            Gebre-Egziabher et al. (2006) "Magnetometer Calibration"
 */

#include "mems_calibration.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

/* ─── L5: Internal — 3x3 Matrix Operations ───────────────────────────────── */

static double mat3_det(const double m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
         - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
         + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

static int mat3_inv(const double a[3][3], double inv[3][3]) {
    double det = mat3_det(a);
    if (fabs(det) < DBL_EPSILON) return -1;

    inv[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
    inv[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
    inv[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
    inv[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
    inv[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
    inv[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
    inv[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
    inv[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det;
    inv[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
    return 0;
}

/* ─── L5: Six-Position Accelerometer Calibration ─────────────────────────
 *
 * In the 6-position method, the sensor is placed in each of six orientations
 * where the gravity vector is well-known:
 *
 *   +Z up:  expected = (0, 0, +g)
 *   -Z down: expected = (0, 0, -g)
 *   +X up:  expected = (+g, 0, 0)
 *   -X down: expected = (-g, 0, 0)
 *   +Y up:  expected = (0, +g, 0)
 *   -Y down: expected = (0, -g, 0)
 *
 * The calibration model: a_calibrated = S · (a_measured - b)
 * where S is the 3×3 scale+misalignment matrix and b is the bias vector.
 *
 * Least-squares solution:
 *   [b, S] = argmin Σᵢ ||a_refᵢ - S·(a_measᵢ - b)||²
 *
 * This is solved by forming the normal equations on the 12 unknowns
 * (3 bias + 9 matrix elements). For the 6-position static case,
 * the solution simplifies significantly.
 */

double mems_calib_accel_6pos(const MemsSixPosData *data,
                              MemsAccelCalib *calib) {
    if (!data || !calib || data->num_positions < 6) return -1.0;
    memset(calib, 0, sizeof(*calib));

    /* Method: Solve for bias and scale factors independently per axis.
     * For axis i, with +g and -g measurements:
     *   a_meas_pos = b_i + S_ii · g + Σ_{j≠i} S_ij · g_j + noise
     *   a_meas_neg = b_i + S_ii · (-g) + Σ_{j≠i} S_ij · g_j + noise
     *
     * Adding:  (a_pos + a_neg)/2 ≈ b_i + Σ_{j≠i} S_ij · g_j
     * Subtracting: (a_pos - a_neg)/2 ≈ S_ii · g
     */

    /* Compute scale factors and bias from paired +g/-g measurements */
    for (int axis = 0; axis < 3; axis++) {
        /* For this axis, find the pair where gravity is along ±axis */
        int idx_pos = -1, idx_neg = -1;
        for (size_t i = 0; i < data->num_positions; i++) {
            if (axis == 0) {
                if (data->expected[i].x > 9.0 && idx_pos < 0) idx_pos = (int)i;
                if (data->expected[i].x < -9.0 && idx_neg < 0) idx_neg = (int)i;
            } else if (axis == 1) {
                if (data->expected[i].y > 9.0 && idx_pos < 0) idx_pos = (int)i;
                if (data->expected[i].y < -9.0 && idx_neg < 0) idx_neg = (int)i;
            } else {
                if (data->expected[i].z > 9.0 && idx_pos < 0) idx_pos = (int)i;
                if (data->expected[i].z < -9.0 && idx_neg < 0) idx_neg = (int)i;
            }
        }

        if (idx_pos < 0 || idx_neg < 0) continue;

        /* Sensor reading for +g and -g orientations */
        double meas_pos, meas_neg, g_val;
        if (axis == 0) {
            meas_pos = data->measured[idx_pos].x;
            meas_neg = data->measured[idx_neg].x;
            g_val = data->expected[idx_pos].x;
        } else if (axis == 1) {
            meas_pos = data->measured[idx_pos].y;
            meas_neg = data->measured[idx_neg].y;
            g_val = data->expected[idx_pos].y;
        } else {
            meas_pos = data->measured[idx_pos].z;
            meas_neg = data->measured[idx_neg].z;
            g_val = data->expected[idx_pos].z;
        }

        /* Bias: average of +g and -g (common offset)
         * Scale factor: 2g / (a_pos - a_neg) → division by signal swing */
        double bias_est = (meas_pos + meas_neg) / 2.0;
        double signal_swing = meas_pos - meas_neg;
        if (fabs(signal_swing) < DBL_EPSILON) continue;

        double sf = (2.0 * g_val) / signal_swing;

        if (axis == 0) {
            calib->bias.x = bias_est;
            calib->scale_factor.x = sf;
            calib->scale_mis.m[0][0] = sf;
        } else if (axis == 1) {
            calib->bias.y = bias_est;
            calib->scale_factor.y = sf;
            calib->scale_mis.m[1][1] = sf;
        } else {
            calib->bias.z = bias_est;
            calib->scale_factor.z = sf;
            calib->scale_mis.m[2][2] = sf;
        }
    }

    /* Compute residuals */
    double rms_sum = 0.0;
    for (size_t i = 0; i < data->num_positions; i++) {
        MemsVector3 corrected = mems_apply_accel_calib(data->measured[i], calib);
        double rx = corrected.x - data->expected[i].x;
        double ry = corrected.y - data->expected[i].y;
        double rz = corrected.z - data->expected[i].z;
        calib->residuals[i] = sqrt(rx * rx + ry * ry + rz * rz);
        rms_sum += calib->residuals[i] * calib->residuals[i];
    }
    calib->rms_error = sqrt(rms_sum / (double)data->num_positions);

    return calib->rms_error;
}

MemsVector3 mems_apply_accel_calib(MemsVector3 raw, const MemsAccelCalib *calib) {
    if (!calib) return raw;

    /* Subtract bias */
    double dx = raw.x - calib->bias.x;
    double dy = raw.y - calib->bias.y;
    double dz = raw.z - calib->bias.z;

    /* Apply scale+misalignment matrix */
    MemsVector3 out;
    out.x = calib->scale_mis.m[0][0] * dx +
            calib->scale_mis.m[0][1] * dy +
            calib->scale_mis.m[0][2] * dz;
    out.y = calib->scale_mis.m[1][0] * dx +
            calib->scale_mis.m[1][1] * dy +
            calib->scale_mis.m[1][2] * dz;
    out.z = calib->scale_mis.m[2][0] * dx +
            calib->scale_mis.m[2][1] * dy +
            calib->scale_mis.m[2][2] * dz;
    return out;
}

/* ─── L5: Rate Table Gyroscope Calibration ─────────────────────────────────
 *
 * For a rate table calibration, the gyro is rotated at known rates
 * ω_ref = {ω₁, ω₂, ..., ω_N}. The measured rate is:
 *
 *   ω_meas = S · (ω_ref + b)  (simplified)
 *
 * Least-squares: b = Σ(ω_meas - S·ω_ref)/N, S from slope.
 * With multiple rates, this becomes a linear regression per axis.
 */

double mems_calib_gyro_rate(const double rates[][3],
                             const double measured[][3],
                             const double accelerations[][3],
                             size_t n, MemsGyroCalib *calib) {
    (void)accelerations;  /* for g-sensitivity extension (full 3D calibration) */
    if (!rates || !measured || !calib || n < 3) return -1.0;
    memset(calib, 0, sizeof(*calib));

    /* Per-axis linear regression: y = ax + b
     * where a = scale factor, b = bias, x = rate, y = measurement */
    for (int axis = 0; axis < 3; axis++) {
        double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;

        for (size_t i = 0; i < n; i++) {
            double x = rates[i][axis];
            double y = measured[i][axis];
            sum_x  += x;
            sum_y  += y;
            sum_xy += x * y;
            sum_xx += x * x;
        }

        double denom = (double)n * sum_xx - sum_x * sum_x;
        if (fabs(denom) < DBL_EPSILON) {
            /* Constant rate: scale factor = 1 (assume), bias = mean */
            double sf = 1.0;
            double bias = (sum_y - sf * sum_x) / (double)n;
            if (axis == 0) { calib->bias.x = bias; calib->scale_factor.x = sf; }
            else if (axis == 1) { calib->bias.y = bias; calib->scale_factor.y = sf; }
            else { calib->bias.z = bias; calib->scale_factor.z = sf; }
        } else {
            double sf = ((double)n * sum_xy - sum_x * sum_y) / denom;
            double bias = (sum_y - sf * sum_x) / (double)n;
            if (axis == 0) {
                calib->bias.x = bias;
                calib->scale_factor.x = sf;
                calib->scale_mis.m[0][0] = sf;
            } else if (axis == 1) {
                calib->bias.y = bias;
                calib->scale_factor.y = sf;
                calib->scale_mis.m[1][1] = sf;
            } else {
                calib->bias.z = bias;
                calib->scale_factor.z = sf;
                calib->scale_mis.m[2][2] = sf;
            }
        }
    }

    /* Compute residuals */
    double rms_sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        MemsVector3 ref = {rates[i][0], rates[i][1], rates[i][2]};
        MemsVector3 meas = {measured[i][0], measured[i][1], measured[i][2]};
        MemsVector3 corrected = mems_apply_gyro_calib(meas, calib);
        double dx = corrected.x - ref.x;
        double dy = corrected.y - ref.y;
        double dz = corrected.z - ref.z;
        calib->residuals[i] = sqrt(dx * dx + dy * dy + dz * dz);
        rms_sum += calib->residuals[i] * calib->residuals[i];
    }
    calib->rms_error = sqrt(rms_sum / (double)n);

    return calib->rms_error;
}

MemsVector3 mems_apply_gyro_calib(MemsVector3 raw, const MemsGyroCalib *calib) {
    if (!calib) return raw;

    double dx = raw.x - calib->bias.x;
    double dy = raw.y - calib->bias.y;
    double dz = raw.z - calib->bias.z;

    MemsVector3 out;
    out.x = calib->scale_mis.m[0][0] * dx +
            calib->scale_mis.m[0][1] * dy +
            calib->scale_mis.m[0][2] * dz;
    out.y = calib->scale_mis.m[1][0] * dx +
            calib->scale_mis.m[1][1] * dy +
            calib->scale_mis.m[1][2] * dz;
    out.z = calib->scale_mis.m[2][0] * dx +
            calib->scale_mis.m[2][1] * dy +
            calib->scale_mis.m[2][2] * dz;
    return out;
}

/* ─── L5: Magnetometer Ellipsoid Fitting ────────────────────────────────────
 *
 * In an ideal world, 3D magnetometer measurements lie on a sphere.
 * Hard-iron distortion shifts the sphere center.
 * Soft-iron distortion stretches the sphere into an ellipsoid.
 *
 * The calibration maps the ellipsoid back to a unit sphere:
 *   m_calibrated = S · (m_raw - h)
 * where h is the hard-iron offset (center) and S is the soft-iron matrix.
 *
 * Simplified iterative method: fit ellipsoid by estimating center from min/max
 * then computing scale factors from the radii.
 */

double mems_calib_mag_ellipsoid(const double measurements[][3], size_t n,
                                 double expected_mag, MemsMagCalib *calib) {
    if (!measurements || !calib || n < 9) return -1.0;
    memset(calib, 0, sizeof(*calib));

    /* Find center (hard iron): midpoint of min and max per axis */
    double x_min = INFINITY, x_max = -INFINITY;
    double y_min = INFINITY, y_max = -INFINITY;
    double z_min = INFINITY, z_max = -INFINITY;

    for (size_t i = 0; i < n; i++) {
        if (measurements[i][0] < x_min) x_min = measurements[i][0];
        if (measurements[i][0] > x_max) x_max = measurements[i][0];
        if (measurements[i][1] < y_min) y_min = measurements[i][1];
        if (measurements[i][1] > y_max) y_max = measurements[i][1];
        if (measurements[i][2] < z_min) z_min = measurements[i][2];
        if (measurements[i][2] > z_max) z_max = measurements[i][2];
    }

    calib->hard_iron.x = (x_min + x_max) / 2.0;
    calib->hard_iron.y = (y_min + y_max) / 2.0;
    calib->hard_iron.z = (z_min + z_max) / 2.0;

    /* Estimate radii of centered data */
    double rx = (x_max - x_min) / 2.0;
    double ry = (y_max - y_min) / 2.0;
    double rz = (z_max - z_min) / 2.0;

    /* Scale factors map average radius to expected magnitude */
    double avg_radius = (rx + ry + rz) / 3.0;
    if (avg_radius < DBL_EPSILON) return -1.0;

    double global_scale = expected_mag / avg_radius;
    calib->scale_factor.x = global_scale;
    calib->scale_factor.y = global_scale;
    calib->scale_factor.z = global_scale;

    calib->soft_iron.m[0][0] = global_scale;
    calib->soft_iron.m[1][1] = global_scale;
    calib->soft_iron.m[2][2] = global_scale;
    calib->field_magnitude = expected_mag;

    /* Compute RMS error after calibration */
    double rms_sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        MemsVector3 raw = {measurements[i][0], measurements[i][1],
                           measurements[i][2]};
        MemsVector3 corrected = mems_apply_mag_calib(raw, calib);
        double mag = sqrt(corrected.x * corrected.x +
                          corrected.y * corrected.y +
                          corrected.z * corrected.z);
        double err = mag - expected_mag;
        rms_sum += err * err;
    }
    calib->rms_error = sqrt(rms_sum / (double)n);

    return calib->rms_error;
}

MemsVector3 mems_apply_mag_calib(MemsVector3 raw, const MemsMagCalib *calib) {
    if (!calib) return raw;

    double dx = raw.x - calib->hard_iron.x;
    double dy = raw.y - calib->hard_iron.y;
    double dz = raw.z - calib->hard_iron.z;

    MemsVector3 out;
    out.x = calib->soft_iron.m[0][0] * dx +
            calib->soft_iron.m[0][1] * dy +
            calib->soft_iron.m[0][2] * dz;
    out.y = calib->soft_iron.m[1][0] * dx +
            calib->soft_iron.m[1][1] * dy +
            calib->soft_iron.m[1][2] * dz;
    out.z = calib->soft_iron.m[2][0] * dx +
            calib->soft_iron.m[2][1] * dy +
            calib->soft_iron.m[2][2] * dz;
    return out;
}

/* ─── L5: Axis Alignment Matrix (Wahba's Problem, Least-Squares) ────────────
 *
 * Given N pairs of measured and reference vectors, find the rotation+scale
 * matrix S that minimizes Σ ||v_ref - S·v_meas||².
 *
 * S = (V_ref · V_meas^T) · (V_meas · V_meas^T)⁻¹
 *
 * This is a simplified solution (without orthogonality constraint).
 * The full Wahba problem (with SO(3) constraint) would use SVD.
 */

double mems_alignment_matrix(const double measured[][3],
                              const double reference[][3],
                              size_t n, MemsMatrix33 *S) {
    if (!measured || !reference || !S || n < 3) return -1.0;
    memset(S, 0, sizeof(*S));

    /* Compute H = Σ v_meas · v_ref^T */
    double H[3][3] = {{0}};
    for (size_t k = 0; k < n; k++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                H[i][j] += measured[k][i] * reference[k][j];
            }
        }
    }

    /* Compute M = Σ v_meas · v_meas^T */
    double M[3][3] = {{0}};
    for (size_t k = 0; k < n; k++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                M[i][j] += measured[k][i] * measured[k][j];
            }
        }
    }

    /* Invert M */
    double M_inv[3][3];
    if (mat3_inv(M, M_inv) != 0) return -1.0;

    /* S = H · M_inv */
    double S_temp[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            S_temp[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                S_temp[i][j] += H[i][k] * M_inv[k][j];
            }
        }
    }

    /* Copy to output */
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            S->m[i][j] = S_temp[i][j];

    /* Condition number approximation via determinant */
    double det_M = mat3_det(M);
    return (fabs(det_M) > DBL_EPSILON) ? 1.0 / fabs(det_M) : INFINITY;
}

MemsVector3 mems_apply_matrix33(MemsVector3 v, const MemsMatrix33 *M) {
    if (!M) return v;
    MemsVector3 out;
    out.x = M->m[0][0] * v.x + M->m[0][1] * v.y + M->m[0][2] * v.z;
    out.y = M->m[1][0] * v.x + M->m[1][1] * v.y + M->m[1][2] * v.z;
    out.z = M->m[2][0] * v.x + M->m[2][1] * v.y + M->m[2][2] * v.z;
    return out;
}

/* ─── L5: Temperature Compensation ───────────────────────────────────────── */

double mems_fit_temp_comp(const double *temperatures, const double *errors,
                           size_t n, int order, MemsTempComp *comp) {
    if (!temperatures || !errors || !comp || n < 3 || order < 1 || order > 4)
        return -1.0;

    memset(comp, 0, sizeof(*comp));

    /* Compute reference temperature as mean of input temperatures */
    double t_sum = 0.0;
    for (size_t i = 0; i < n; i++) t_sum += temperatures[i];
    comp->ref_temp = t_sum / (double)n;

    /* Build normal equations for polynomial: e = a₁·ΔT + a₂·ΔT² + a₃·ΔT³
     * X matrix: each row = [ΔT, ΔT², ΔT³] for order=3
     * Solve X^T·X · a = X^T·y
     */
    /* For small order (1-4), use direct Gaussian elimination */

    /* Build A = X^T·X and b = X^T·y (accumulated) */
    double A[4][4] = {{0}};
    double b_vec[4] = {0};

    for (size_t k = 0; k < n; k++) {
        double dt = temperatures[k] - comp->ref_temp;
        double pow_dt[5];
        pow_dt[0] = 1.0;
        for (int p = 1; p <= order; p++) pow_dt[p] = pow_dt[p-1] * dt;

        for (int i = 0; i < order; i++) {
            /* pow_dt[i+1] corresponds to ΔT^{i+1} */
            b_vec[i] += pow_dt[i+1] * errors[k];
            for (int j = 0; j < order; j++) {
                A[i][j] += pow_dt[i+1] * pow_dt[j+1];
            }
        }
    }

    /* Gaussian elimination with partial pivoting */
    for (int col = 0; col < order; col++) {
        /* Find pivot */
        int max_row = col;
        double max_val = fabs(A[col][col]);
        for (int row = col + 1; row < order; row++) {
            if (fabs(A[row][col]) > max_val) {
                max_val = fabs(A[row][col]);
                max_row = row;
            }
        }
        /* Swap rows */
        if (max_row != col) {
            for (int j = 0; j < order; j++) {
                double tmp = A[col][j];
                A[col][j] = A[max_row][j];
                A[max_row][j] = tmp;
            }
            double tmp = b_vec[col];
            b_vec[col] = b_vec[max_row];
            b_vec[max_row] = tmp;
        }

        if (fabs(A[col][col]) < DBL_EPSILON) return -1.0;

        /* Eliminate below */
        for (int row = col + 1; row < order; row++) {
            double factor = A[row][col] / A[col][col];
            for (int j = col; j < order; j++)
                A[row][j] -= factor * A[col][j];
            b_vec[row] -= factor * b_vec[col];
        }
    }

    /* Back substitution */
    double a[4] = {0};
    for (int i = order - 1; i >= 0; i--) {
        double sum = b_vec[i];
        for (int j = i + 1; j < order; j++)
            sum -= A[i][j] * a[j];
        a[i] = sum / A[i][i];
    }

    comp->a1 = (order >= 1) ? a[0] : 0.0;
    comp->a2 = (order >= 2) ? a[1] : 0.0;
    comp->a3 = (order >= 3) ? a[2] : 0.0;

    /* Compute R² */
    double ss_res = 0.0, ss_tot = 0.0, y_mean = 0.0;
    for (size_t i = 0; i < n; i++) y_mean += errors[i];
    y_mean /= (double)n;

    for (size_t i = 0; i < n; i++) {
        double dt = temperatures[i] - comp->ref_temp;
        double fitted = 0.0;
        double dt_pow = dt;
        if (order >= 1) { fitted += a[0] * dt_pow; dt_pow *= dt; }
        if (order >= 2) { fitted += a[1] * dt_pow; dt_pow *= dt; }
        if (order >= 3) { fitted += a[2] * dt_pow; }
        ss_res += (errors[i] - fitted) * (errors[i] - fitted);
        ss_tot += (errors[i] - y_mean) * (errors[i] - y_mean);
    }

    if (ss_tot > DBL_EPSILON) {
        comp->fit_r_squared = 1.0 - ss_res / ss_tot;
    } else {
        comp->fit_r_squared = 1.0;
    }

    return comp->fit_r_squared;
}

double mems_temp_delta(double temperature, const MemsTempComp *comp) {
    if (!comp) return 0.0;
    return temperature - comp->ref_temp;
}

double mems_apply_temp_comp(double raw_value, double temperature,
                             const MemsTempComp *comp) {
    if (!comp) return raw_value;

    double dt = temperature - comp->ref_temp;
    /* Polynomial model: denominator = 1 + a1·ΔT + a2·ΔT² + a3·ΔT³ */
    double denom = 1.0 + comp->a1 * dt + comp->a2 * dt * dt +
                   comp->a3 * dt * dt * dt;

    if (fabs(denom) < DBL_EPSILON) return raw_value;
    return raw_value / denom;
}
