/**
 * @file    sensor_applications.c
 * @brief   Smart sensor real-world application implementations
 *
 * Knowledge Coverage:
 *   L7 (Applications): Industrial IoT sensor node (ISO 13374 condition
 *                      monitoring), automotive TPMS (tire pressure
 *                      monitoring), medical pulse oximeter signal chain,
 *                      environmental air quality monitoring (PM2.5/CO2),
 *                      structural health monitoring (strain + temperature)
 *
 * Reference:
 *   ISO 13374 — Condition monitoring and diagnostics of machines
 *   ISO 21750:2006 — Tyre Pressure Monitoring Systems (TPMS)
 *   Webster, J.G. (2009), "Medical Instrumentation: Application and
 *     Design", 4th ed., Wiley
 *   Fraden, J. (2016), "Handbook of Modern Sensors", 5th ed., Springer
 *
 * Course Mapping:
 *   MIT 6.003 — Signals and Systems (sensor applications)
 *   Michigan EECS 411 — Microwave Circuits (automotive radar/sensors)
 *   Georgia Tech ECE 6601 — Communications (IoT sensor networks)
 */

#include "smart_sensor.h"
#include "sensor_conditioning.h"
#include "sensor_calibration.h"
#include "sensor_fusion.h"
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =========================================================================
 * L7 Application 1: Industrial IoT Vibration Sensor Node
 *
 * Implements ISO 13374 condition monitoring for rotating machinery.
 * A smart vibration sensor node that captures accelerometer data,
 * computes RMS velocity (ISO 10816 vibration severity), and performs
 * FFT-based spectral analysis for bearing fault detection.
 *
 * ISO 10816-3 vibration severity zones (mm/s RMS, 10-1000 Hz):
 *   Zone A (new):  < 1.4 mm/s
 *   Zone B (good): 1.4-2.8 mm/s
 *   Zone C (satisfactory): 2.8-4.5 mm/s
 *   Zone D (alarm): > 4.5 mm/s
 *
 * Bearing fault frequencies (as multiples of shaft speed):
 *   BPFO (outer race) ≈ (N/2) * (1 - Bd/Pd * cos(phi)) * RPM/60
 *   BPFI (inner race) ≈ (N/2) * (1 + Bd/Pd * cos(phi)) * RPM/60
 *   BSF (ball spin)    ≈ (Pd/(2*Bd)) * (1 - (Bd/Pd)^2 * cos^2(phi)) * RPM/60
 *   FTF (cage)         ≈ (1/2) * (1 - Bd/Pd * cos(phi)) * RPM/60
 *
 * where N=number of balls, Bd=ball diameter, Pd=pitch diameter, phi=contact angle
 * ========================================================================= */

/**
 * @brief Compute ISO 10816 vibration severity from acceleration RMS
 *
 * Given RMS acceleration (in g), converts to velocity RMS (mm/s)
 * assuming a dominant frequency. For broadband vibration, use the
 * standard 10-1000 Hz band after bandpass filtering.
 *
 * Conversion: v_rms = a_rms * 9.81 * 1000 / (2*pi*f)
 * where a_rms is in g, output is in mm/s, f is dominant frequency in Hz.
 *
 * @param accel_rms_g       RMS acceleration in g (9.81 m/s^2)
 * @param dominant_freq_hz  Dominant vibration frequency (Hz)
 * @return                  RMS velocity in mm/s (ISO 10816)
 */
double ss_app_vibration_severity(double accel_rms_g, double dominant_freq_hz) {
    double accel_mm_s2 = accel_rms_g * 9810.0; /* g to mm/s^2 */
    double omega = 2.0 * M_PI * dominant_freq_hz;

    if (omega < 1e-6) return 0.0;

    /* Integration: acceleration -> velocity in frequency domain:
     * V(omega) = A(omega) / (j*omega)
     * |V| = |A| / omega */
    return accel_mm_s2 / omega;
}

/**
 * @brief ISO 10816 vibration severity zone classification
 *
 * @param velocity_rms_mm_s  RMS velocity in mm/s
 * @return                   0=Zone A(new), 1=Zone B(good),
 *                           2=Zone C(acceptable), 3=Zone D(alarm)
 */
int ss_app_vibration_zone(double velocity_rms_mm_s) {
    if (velocity_rms_mm_s < 1.4) return 0;
    if (velocity_rms_mm_s < 2.8) return 1;
    if (velocity_rms_mm_s < 4.5) return 2;
    return 3;
}

/**
 * @brief Calculate ball-pass frequency outer race (BPFO)
 *
 * BPFO = (N/2) * RPM/60 * (1 - (Bd/Pd) * cos(phi))
 *
 * @param n_balls       Number of rolling elements
 * @param rpm           Shaft speed in revolutions per minute
 * @param ball_diameter Ball diameter (mm)
 * @param pitch_diameter Pitch diameter (mm)
 * @param contact_angle Contact angle (radians, typically 0 for deep groove)
 * @return              BPFO frequency in Hz
 */
double ss_app_bpfo_frequency(int n_balls, double rpm,
                             double ball_diameter, double pitch_diameter,
                             double contact_angle) {
    double shaft_freq = rpm / 60.0;
    double ratio = ball_diameter / pitch_diameter;
    if (pitch_diameter <= 0.0) return 0.0;
    return (double)n_balls / 2.0 * shaft_freq * (1.0 - ratio * cos(contact_angle));
}

/**
 * @brief Calculate ball-pass frequency inner race (BPFI)
 *
 * BPFI = (N/2) * RPM/60 * (1 + (Bd/Pd) * cos(phi))
 */
double ss_app_bpfi_frequency(int n_balls, double rpm,
                             double ball_diameter, double pitch_diameter,
                             double contact_angle) {
    double shaft_freq = rpm / 60.0;
    double ratio = ball_diameter / pitch_diameter;
    if (pitch_diameter <= 0.0) return 0.0;
    return (double)n_balls / 2.0 * shaft_freq * (1.0 + ratio * cos(contact_angle));
}

/* =========================================================================
 * L7 Application 2: Automotive TPMS (Tire Pressure Monitoring)
 *
 * Direct TPMS uses pressure + temperature sensors mounted inside
 * each tire. The sensor measures absolute pressure and internal
 * temperature, transmitting data wirelessly (typically 315/433 MHz)
 * to the vehicle ECU.
 *
 * Key calculations:
 *   - Pressure compensation for temperature (ideal gas law)
 *   - Under-inflation detection threshold (typically 25% below placard)
 *   - Slow leak detection (rate of pressure loss)
 *
 * Reference: ISO 21750:2006, FMVSS 138 (US TPMS regulation)
 * ========================================================================= */

/**
 * @brief Compensate tire pressure to reference temperature
 *
 * Uses the ideal gas law (Gay-Lussac's law for constant volume):
 *   P_ref = P_meas * (T_ref / T_meas)
 * where temperatures are in Kelvin.
 *
 * This is critical because tire pressure changes approximately
 * 1 psi (6.9 kPa) per 10 degF (5.6 degC) temperature change.
 *
 * @param pressure_kpa     Measured absolute pressure (kPa)
 * @param temp_meas_c      Measured tire internal temperature (degC)
 * @param temp_ref_c       Reference temperature, typically 20 degC
 * @return                 Temperature-compensated pressure (kPa)
 */
double ss_app_tpms_compensate_pressure(double pressure_kpa,
                                       double temp_meas_c,
                                       double temp_ref_c) {
    double T_meas_K = temp_meas_c + 273.15;
    double T_ref_K  = temp_ref_c  + 273.15;

    if (T_meas_K <= 0.0) return pressure_kpa;
    return pressure_kpa * (T_ref_K / T_meas_K);
}

/**
 * @brief Compute tire pressure deviation from placard value
 *
 * Returns the deviation as a percentage. FMVSS 138 requires
 * a warning when any tire is 25% or more below placard pressure.
 *
 * @param measured_pressure_kpa  Current compensated pressure (kPa)
 * @param placard_pressure_kpa   Manufacturer recommended pressure (kPa)
 * @return                       Deviation in percent (negative = under-inflated)
 */
double ss_app_tpms_deviation_percent(double measured_pressure_kpa,
                                     double placard_pressure_kpa) {
    if (placard_pressure_kpa <= 0.0) return 0.0;
    return 100.0 * (measured_pressure_kpa - placard_pressure_kpa)
           / placard_pressure_kpa;
}

/**
 * @brief Detect slow tire leak from pressure trend
 *
 * Analyzes the rate of pressure change over time to detect
 * slow leaks (typically 0.5-2 psi/month or 3.4-13.8 kPa/month).
 *
 * A negative rate exceeding the threshold for a sustained period
 * indicates a slow leak.
 *
 * @param pressure_history  Array of pressure readings [n_samples]
 * @param time_hours        Time of each reading in hours [n_samples]
 * @param n_samples         Number of readings (>=2)
 * @param leak_rate_out     Output: estimated leak rate in kPa/hour (negative = leak)
 * @return                  1 if leak detected, 0 if normal
 */
int ss_app_tpms_detect_slow_leak(const double *pressure_history,
                                 const double *time_hours,
                                 size_t n_samples,
                                 double *leak_rate_out) {
    size_t i;
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    double slope, threshold = -0.001; /* -0.001 kPa/hr ≈ -0.02 kPa/day */

    if (!pressure_history || !time_hours || n_samples < 2) {
        if (leak_rate_out) *leak_rate_out = 0.0;
        return 0;
    }

    for (i = 0; i < n_samples; i++) {
        sx  += time_hours[i];
        sy  += pressure_history[i];
        sxx += time_hours[i] * time_hours[i];
        sxy += time_hours[i] * pressure_history[i];
    }

    {
        double denom = (double)n_samples * sxx - sx * sx;
        if (fabs(denom) < 1e-30) {
            if (leak_rate_out) *leak_rate_out = 0.0;
            return 0;
        }
        slope = ((double)n_samples * sxy - sx * sy) / denom;
    }

    if (leak_rate_out) *leak_rate_out = slope;
    return (slope < threshold) ? 1 : 0;
}

/* =========================================================================
 * L7 Application 3: Medical Pulse Oximeter Signal Processing
 *
 * Pulse oximetry measures arterial blood oxygen saturation (SpO2)
 * non-invasively using two wavelengths of light (typically 660 nm
 * red and 940 nm infrared) passed through a finger or earlobe.
 *
 * Principle: Beer-Lambert law for light absorption in tissue.
 * Oxygenated hemoglobin (HbO2) and deoxygenated hemoglobin (Hb)
 * have different absorption spectra at these wavelengths.
 *
 * The ratio of ratios R is computed from the AC and DC components
 * of the photoplethysmogram (PPG) at each wavelength:
 *   R = (AC_red / DC_red) / (AC_ir / DC_ir)
 *
 * Empirical calibration: SpO2 = a - b * R
 * (typical: a=110, b=25 for most commercial pulse oximeters)
 *
 * Reference: Webster, J.G. (2009), "Medical Instrumentation", 4th ed., Wiley: Ch. 11
 *            Mendelson, Y. (1992), Clin. Chem. 38(9):1601-1607
 * ========================================================================= */

/**
 * @brief Compute SpO2 from PPG ratio-of-ratios
 *
 * Uses the empirical linear calibration model:
 *   SpO2 = a - b * R
 *
 * where R = (AC_red/DC_red) / (AC_ir/DC_ir)
 *
 * Typical calibration coefficients:
 *   a = 110, b = 25 for transmission-mode finger probe
 *
 * @param ac_red   AC component (pulsatile) at 660 nm red
 * @param dc_red   DC component (non-pulsatile) at 660 nm red
 * @param ac_ir    AC component at 940 nm infrared
 * @param dc_ir    DC component at 940 nm infrared
 * @param cal_a    Calibration coefficient a (default 110)
 * @param cal_b    Calibration coefficient b (default 25)
 * @return         Estimated SpO2 in percent (0-100)
 */
double ss_app_spo2_estimate(double ac_red, double dc_red,
                            double ac_ir, double dc_ir,
                            double cal_a, double cal_b) {
    double R;

    if (dc_red <= 0.0 || dc_ir <= 0.0) return 0.0;

    /* Ratio of ratios */
    R = (ac_red / dc_red) / (ac_ir / dc_ir);

    /* Linear calibration */
    return cal_a - cal_b * R;
}

/**
 * @brief Estimate heart rate from PPG signal period
 *
 * Heart rate (BPM) = 60 / period_seconds
 *
 * The PPG signal period can be estimated from the time between
 * successive systolic peaks.
 *
 * @param peak_intervals_ms  Array of peak-to-peak intervals (ms)
 * @param n_intervals        Number of intervals (>=1)
 * @return                   Estimated heart rate in BPM
 */
double ss_app_heart_rate_from_ppg(const double *peak_intervals_ms,
                                  size_t n_intervals) {
    size_t i;
    double sum = 0.0;

    if (!peak_intervals_ms || n_intervals == 0) return 0.0;

    for (i = 0; i < n_intervals; i++) {
        sum += peak_intervals_ms[i];
    }

    /* Average interval in seconds */
    double avg_interval_s = (sum / (double)n_intervals) / 1000.0;

    if (avg_interval_s <= 0.0) return 0.0;

    return 60.0 / avg_interval_s;
}

/**
 * @brief Compute perfusion index from PPG signal
 *
 * Perfusion Index (PI) = AC / DC * 100%
 *
 * PI indicates the strength of peripheral blood flow.
 *   PI > 5%: strong perfusion (good signal)
 *   PI 1-5%: normal perfusion
 *   PI < 1%: weak perfusion (poor signal, may be unreliable)
 *
 * @param ac  Pulsatile component (peak-to-peak)
 * @param dc  Non-pulsatile baseline component
 * @return    Perfusion index in percent
 */
double ss_app_perfusion_index(double ac, double dc) {
    if (dc <= 0.0) return 0.0;
    return 100.0 * ac / dc;
}

/* =========================================================================
 * L7 Application 4: Environmental Air Quality Monitoring
 *
 * Combines multiple gas sensors (CO2, PM2.5, VOC, temperature,
 * humidity) to compute the Air Quality Index (AQI) and provide
 * actionable recommendations.
 *
 * PM2.5 AQI breakpoints (EPA, 24-hour average, ug/m^3):
 *   Good (0-50):        0.0 - 12.0
 *   Moderate (51-100):  12.1 - 35.4
 *   Unhealthy-S (101-150): 35.5 - 55.4
 *   Unhealthy (151-200):   55.5 - 150.4
 *   Very Unhealthy (201-300): 150.5 - 250.4
 *   Hazardous (301-500):     250.5 - 500.4
 *
 * CO2 concentration guidelines (ppm):
 *   400-600: excellent (outdoor equivalent)
 *   600-1000: good (well-ventilated indoor)
 *   1000-2000: moderate (stuffy — ventilate)
 *   2000-5000: poor (drowsiness, headaches possible)
 *   >5000: hazardous (OSHA 8-hour limit = 5000 ppm)
 * ========================================================================= */

/**
 * @brief Compute EPA PM2.5 AQI from 24-hour average concentration
 *
 * Uses piecewise-linear breakpoint formula per EPA AQI standard:
 *   AQI = (AQI_hi - AQI_lo)/(Conc_hi - Conc_lo) * (Conc - Conc_lo) + AQI_lo
 *
 * @param pm25_ug_m3  24-hour average PM2.5 concentration (ug/m^3)
 * @return            AQI value (0-500, capped)
 */
double ss_app_aqi_pm25(double pm25_ug_m3) {
    /* EPA PM2.5 AQI breakpoints */
    static const double conc_lo[] = {0.0, 12.1, 35.5, 55.5, 150.5, 250.5};
    static const double conc_hi[] = {12.0, 35.4, 55.4, 150.4, 250.4, 500.4};
    static const double aqi_lo[]  = {0, 51, 101, 151, 201, 301};
    static const double aqi_hi[]  = {50, 100, 150, 200, 300, 500};
    int i;

    if (pm25_ug_m3 <= 0.0) return 0.0;
    if (pm25_ug_m3 > 500.4) return 500.0;

    for (i = 0; i < 6; i++) {
        if (pm25_ug_m3 <= conc_hi[i]) {
            return (aqi_hi[i] - aqi_lo[i]) / (conc_hi[i] - conc_lo[i])
                   * (pm25_ug_m3 - conc_lo[i]) + aqi_lo[i];
        }
    }

    return 500.0;
}

/**
 * @brief CO2 concentration comfort level assessment
 *
 * @param co2_ppm  CO2 concentration in parts per million
 * @return         0=excellent, 1=good, 2=moderate, 3=poor, 4=hazardous
 */
int ss_app_co2_comfort_level(double co2_ppm) {
    if (co2_ppm < 600.0)  return 0; /* Excellent */
    if (co2_ppm < 1000.0) return 1; /* Good */
    if (co2_ppm < 2000.0) return 2; /* Moderate */
    if (co2_ppm < 5000.0) return 3; /* Poor */
    return 4;                        /* Hazardous */
}

/**
 * @brief Compute dew point from temperature and relative humidity
 *
 * Magnus formula (approximate, valid -45 to +60 degC):
 *   gamma = ln(RH/100) + (b*T) / (c+T)
 *   T_dew = (c * gamma) / (b - gamma)
 *
 * Constants (over water): b=17.27, c=237.7
 *
 * Reference: Sonntag, D. (1990), Z. Meteorol. 40(5):340-344
 *
 * @param temp_c       Air temperature in degC
 * @param rel_humidity Relative humidity in percent (0-100)
 * @return             Dew point temperature in degC
 */
double ss_app_dew_point(double temp_c, double rel_humidity) {
    double b = 17.27;
    double c = 237.7;
    double gamma;

    if (rel_humidity <= 0.0) return -273.15; /* No moisture */
    if (rel_humidity >= 100.0) return temp_c;  /* Saturated */

    gamma = log(rel_humidity / 100.0) + (b * temp_c) / (c + temp_c);
    return (c * gamma) / (b - gamma);
}

/**
 * @brief Compute heat index (apparent temperature) from temperature and humidity
 *
 * Rothfusz regression (NOAA/NWS, simplified for T >= 27 degC, RH >= 40%):
 *   HI = c1 + c2*T + c3*R + c4*T*R + c5*T^2 + c6*R^2 + c7*T^2*R + c8*T*R^2 + c9*T^2*R^2
 *
 * where T in degF, R in percent.
 *
 * @param temp_c       Air temperature in degC (converted internally to degF)
 * @param rel_humidity Relative humidity in percent
 * @return             Heat index in degC
 */
double ss_app_heat_index(double temp_c, double rel_humidity) {
    double T = temp_c * 9.0 / 5.0 + 32.0; /* Convert to degF */
    double R = rel_humidity;
    double HI;

    /* Simple formula for moderate conditions */
    if (T < 80.0 || R < 40.0) return temp_c;

    /* Rothfusz regression coefficients (NOAA) */
    HI = -42.379
         + 2.04901523 * T
         + 10.14333127 * R
         - 0.22475541 * T * R
         - 6.83783e-3 * T * T
         - 5.481717e-2 * R * R
         + 1.22874e-3 * T * T * R
         + 8.5282e-4 * T * R * R
         - 1.99e-6 * T * T * R * R;

    /* Convert back to degC */
    return (HI - 32.0) * 5.0 / 9.0;
}
