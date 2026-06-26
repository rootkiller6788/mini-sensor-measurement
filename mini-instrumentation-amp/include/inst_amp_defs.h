/**
 * inst_amp_defs.h - Core Type Definitions for Instrumentation Amplifiers
 *
 * L1 Definitions: amplifier topologies, gain structures, error types,
 * precision parameters, and bridge sensor specifications.
 *
 * Courses: Stanford EE247, Berkeley EE105, Michigan EECS 411,
 *          Georgia Tech ECE 6601, TU Munich High-Frequency Engineering
 * Reference:
 *   Kitchin & Counts, "A Designer's Guide to Instrumentation Amplifiers"
 *     (Analog Devices, 3rd ed., 2006)
 *   Sedra & Smith, "Microelectronic Circuits" (2020), Ch. 2, 9
 *   Franco, "Design with Operational Amplifiers and Analog ICs" (4th ed., 2015)
 *   Horowitz & Hill, "The Art of Electronics" (3rd ed., 2015), Ch. 5
 */

#ifndef INST_AMP_DEFS_H
#define INST_AMP_DEFS_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==================================================================
 * L1: Instrumentation Amplifier Topology Enumeration (7 topologies)
 *
 * The three-op-amp topology (classic IA) provides the best balance
 * of CMRR, input impedance, and gain adjustability. First described
 * in full by Tobey, Graeme, and Huelsman (1971).
 *
 * The two-op-amp topology saves one op-amp at the cost of asymmetric
 * common-mode path delay, limiting CMRR at high frequencies.
 *
 * Current-mode architectures (CCII-based) are preferred for low-voltage
 * operation below 2V supply rails (Toumazou, Lidgey, Haigh, 1990).
 * ================================================================== */

typedef enum {
    IA_TOPO_THREE_OPAMP      = 0,  /* Classic 3-op-amp: A1,A2 buffers + A3 diff amp */
    IA_TOPO_TWO_OPAMP        = 1,  /* 2-op-amp: differential + non-inverting gain stage */
    IA_TOPO_CURRENT_MODE     = 2,  /* Current-mode: transconductance + transimpedance */
    IA_TOPO_FLYING_CAP       = 3,  /* Flying-capacitor: auto-zero, CMRR > 140dB */
    IA_TOPO_INDIRECT_CURRENT = 4,  /* Indirect current feedback: wide bandwidth */
    IA_TOPO_DIFF_DIFF_AMP    = 5,  /* Differential-difference amplifier (DDA) */
    IA_TOPO_CCII_BASED       = 6   /* Second-generation current conveyor (CCII) */
} ia_topology_t;

/** Input coupling mode for sensor interfaces */
typedef enum {
    IA_COUPLING_DC           = 0,  /* Direct DC coupling, resistive sensors */
    IA_COUPLING_AC           = 1,  /* AC coupling through capacitor, biopotentials */
    IA_COUPLING_TRANSFORMER  = 2,  /* Transformer isolation, high CMRR at RF */
    IA_COUPLING_OPTO         = 3,  /* Opto-isolated, patient safety (IEC 60601) */
    IA_COUPLING_CHOPPER      = 4   /* Chopper-stabilized, nV-level DC precision */
} ia_coupling_t;

/** Gain programming method */
typedef enum {
    IA_GAIN_FIXED            = 0,  /* Fixed gain set by internal resistors */
    IA_GAIN_RESISTOR         = 1,  /* Single external resistor R_G sets gain */
    IA_GAIN_DIGITAL          = 2,  /* Digital SPI/I2C programmable gain */
    IA_GAIN_PIN_STRAP        = 3,  /* Pin-selectable gain: A0, A1, A2 levels */
    IA_GAIN_AUTORANGE        = 4   /* Automatic gain ranging based on input level */
} ia_gain_mode_t;

/* ==================================================================
 * L1: Precision Error Types - The 10 canonical DC/AC errors
 *
 * Each error source has distinct physical origin and temperature
 * behavior. Total error budget (RSS):
 *   V_err_total = sqrt(V_os^2 + (I_os*R_s)^2 + (V_cm/CMRR)^2 + ...)
 * Reference: Kitchin & Counts (2006), Ch. 3 - Error Budget Analysis
 * ================================================================== */

typedef enum {
    IA_ERR_INPUT_OFFSET_VOLTAGE  = 0,  /* V_OS: differential input offset voltage */
    IA_ERR_INPUT_OFFSET_CURRENT  = 1,  /* I_OS: I_B+ - I_B- mismatch */
    IA_ERR_INPUT_BIAS_CURRENT    = 2,  /* I_B: average of I_B+ and I_B- */
    IA_ERR_GAIN_ERROR            = 3,  /* Fractional deviation from nominal gain */
    IA_ERR_GAIN_NONLINEARITY     = 4,  /* Deviation from ideal linear Vout vs Vin */
    IA_ERR_CMRR_FINITE           = 5,  /* CMRR: common-mode couples to output */
    IA_ERR_PSRR_FINITE           = 6,  /* PSRR: supply ripple couples to output */
    IA_ERR_INPUT_IMPEDANCE       = 7,  /* Finite Z_in loads source, creates divider */
    IA_ERR_THERMAL_DRIFT         = 8,  /* Temperature drift of all above parameters */
    IA_ERR_NOISE_RTI             = 9,  /* Referred-to-input noise: thermal + 1/f */
    IA_ERR_COUNT                    = 10
} ia_error_type_t;

/* ==================================================================
 * L1: Instrumentation Amplifier Core Specification Structure
 *
 * Captures the complete datasheet parameters for a precision IA.
 * Modeled after AD620, AD8421, INA128, INA333 datasheets.
 * ================================================================== */

typedef struct {
    /* Gain parameters */
    double      gain;               /* Voltage gain Vout/Vin (V/V) */
    double      gain_min;           /* Minimum programmable gain */
    double      gain_max;           /* Maximum programmable gain */
    double      gain_error_pct;     /* Gain error at 25C (%) */
    double      gain_drift_ppm;     /* Gain drift (ppm/C) */
    double      gain_nonlinearity;  /* Integral nonlinearity (ppm of FS) */

    /* DC Precision */
    double      vos_uv;             /* Input offset voltage (uV) */
    double      vos_drift_nv_per_c; /* Offset drift (nV/C) */
    double      ib_pa;              /* Input bias current (pA) */
    double      ios_pa;             /* Input offset current (pA) */
    double      ib_drift_pa_per_c;  /* Bias current drift (pA/C) */

    /* AC Performance */
    double      bw_khz;             /* -3dB bandwidth at gain (kHz) */
    double      sr_v_per_us;        /* Slew rate (V/us) */
    double      settling_time_us;   /* Settling to 0.01% (us) */
    double      gbp_mhz;            /* Gain-bandwidth product (MHz) */

    /* Rejection */
    double      cmrr_db;            /* Common-mode rejection ratio at DC (dB) */
    double      cmrr_at_50hz_db;    /* CMRR at 50/60 Hz (dB), critical for ECG */
    double      psrr_db;            /* Power supply rejection ratio (dB) */

    /* Noise */
    double      en_rti_nv_per_rhz;  /* Input voltage noise density (nV/rtHz) */
    double      en_0p1_10hz_uvpp;   /* 0.1-10 Hz peak-to-peak noise (uVpp) */
    double      in_rti_fa_per_rhz;  /* Input current noise density (fA/rtHz) */
    double      en_corner_hz;        /* 1/f noise corner frequency (Hz) */

    /* Input/Output Range */
    double      zin_diff_mohm;      /* Differential input impedance (MOhm) */
    double      zin_cm_mohm;        /* Common-mode input impedance (MOhm) */
    double      input_range_v;      /* Input common-mode range (+/-V) */
    double      output_swing_v;     /* Output voltage swing (+/-V) */
    double      supply_min_v;       /* Minimum supply voltage */
    double      supply_max_v;       /* Maximum supply voltage */
    double      iq_ua;              /* Quiescent current (uA) */
} ia_spec_t;

/* ==================================================================
 * L1: Wheatstone Bridge Sensor Model
 *
 * The Wheatstone bridge is the fundamental sensor interface for
 * strain gauges, pressure sensors, RTDs, and load cells.
 *
 * Circuit topology (full bridge):
 *             V_exc
 *        R1 /    \ R2
 *          /      \
 *    V_a o-        -o V_b  (differential output = V_a - V_b)
 *          \      /
 *        R4 \    / R3
 *             |
 *            GND
 *
 * Exact transfer function:
 *   V_diff = V_exc * [R3/(R2+R3) - R4/(R1+R4)]
 *
 * For R1=R2=R3=R4=R and one active arm with change dR (dR << R):
 *   V_diff = V_exc * dR / (4R + 2*dR)  (exact)
 *   V_diff ~ V_exc * dR / (4R)          (linear approximation)
 *
 * Reference: Window & Holister, "Strain Gauge Technology" (1982)
 * ================================================================== */

typedef struct {
    double      r_nominal;          /* Nominal resistance (Ohm), typical 120/350/1000 */
    double      delta_r;            /* Resistance change due to measurand (Ohm) */
    double      v_excitation;       /* Bridge excitation voltage (V) */
    unsigned    active_arms;        /* 1=quarter, 2=half, 4=full bridge */
    double      tcr_ppm_per_c;      /* Temperature coefficient of resistance (ppm/C) */
    double      lead_resistance;    /* Lead wire resistance, 3W/4W config (Ohm) */
} bridge_sensor_t;

/** Bridge configuration types */
typedef enum {
    BRIDGE_QUARTER     = 1,  /* Single active arm, needs temperature compensation */
    BRIDGE_HALF        = 2,  /* Two active arms (adjacent), self-temp-compensated */
    BRIDGE_FULL        = 4,  /* Four active arms, maximum sensitivity, compensated */
    BRIDGE_QUARTER_3W  = 5,  /* Quarter bridge with 3-wire lead compensation */
    BRIDGE_HALF_4W     = 6   /* Half bridge with 4-wire Kelvin connection */
} bridge_config_t;

/* ==================================================================
 * L1: Sensor Types Library (14 sensor types commonly used with IAs)
 * ================================================================== */

typedef enum {
    SENSOR_STRAIN_GAUGE      = 0,  /* Foil/semiconductor strain gauge, 120/350 Ohm */
    SENSOR_LOAD_CELL         = 1,  /* Full-bridge load cell, mV/V output */
    SENSOR_PRESSURE_BRIDGE   = 2,  /* Piezoresistive pressure sensor bridge */
    SENSOR_THERMOCOUPLE      = 3,  /* Type K/J/T/E/R/S/B/N thermocouple, uV/C */
    SENSOR_RTD               = 4,  /* Pt100/Pt1000 RTD, alpha=0.00385 */
    SENSOR_THERMISTOR        = 5,  /* NTC/PTC thermistor, high sensitivity, nonlinear */
    SENSOR_PHOTODIODE        = 6,  /* Photodiode transimpedance front-end */
    SENSOR_HALL_EFFECT       = 7,  /* Hall effect sensor, mV/mT sensitivity */
    SENSOR_ACCELEROMETER     = 8,  /* MEMS capacitive accelerometer bridge */
    SENSOR_CURRENT_SHUNT     = 9,  /* Shunt resistor for current sensing */
    SENSOR_LVDT              = 10, /* Linear variable differential transformer */
    SENSOR_ECG_ELECTRODE     = 11, /* Biopotential electrode, ~1mV, 0.05-150Hz */
    SENSOR_EEG_ELECTRODE     = 12, /* EEG electrode, ~10-100uV, 0.5-50Hz */
    SENSOR_EMG_ELECTRODE     = 13  /* EMG electrode, ~0.1-5mV, 20-500Hz */
} sensor_type_t;

/* ==================================================================
 * L1: Operating State, Calibration Data, and Error Codes
 * ================================================================== */

typedef struct {
    bool        powered;            /* Amplifier is powered on */
    bool        calibrated;         /* Calibration has been performed */
    bool        overload;           /* Input over-range detected */
    bool        saturation;         /* Output saturated (clipped) */
    bool        thermal_warning;    /* Junction temperature above warning level */
    bool        open_input;         /* Input disconnected (open circuit detected) */
    double      junction_temp_c;    /* Estimated junction temperature (C) */
    double      power_mw;           /* Estimated power consumption (mW) */
} ia_state_t;

/** Calibration coefficients for 2-point calibration */
typedef struct {
    double      offset_v;           /* Measured offset voltage at zero input (V) */
    double      gain_slope;         /* Measured gain slope (V/V) */
    double      offset_drift;       /* Offset temperature drift coefficient (V/C) */
    double      gain_drift;         /* Gain temperature drift coefficient (1/C) */
    int         cal_date;           /* Calibration date (YYYYMMDD) */
    bool        valid;              /* Calibration data is valid */
} ia_calibration_t;

/** Error codes returned by library functions */
typedef enum {
    IA_OK                    = 0,   /* Operation successful */
    IA_ERR_NULL_POINTER      = 1,   /* NULL pointer argument */
    IA_ERR_INVALID_GAIN      = 2,   /* Gain out of valid range */
    IA_ERR_SATURATED         = 3,   /* Output is saturated */
    IA_ERR_OVERLOAD          = 4,   /* Input over-range */
    IA_ERR_CALIBRATION       = 5,   /* Calibration failed or expired */
    IA_ERR_THERMAL           = 6,   /* Thermal shutdown */
    IA_ERR_PARAMETER         = 7,   /* Invalid configuration parameter */
    IA_ERR_TIMEOUT           = 8,   /* Operation timed out */
    IA_ERR_NOT_SUPPORTED     = 9    /* Feature not supported in this topology */
} ia_error_code_t;

/* ==================================================================
 * L1: Thermocouple Type Definitions
 *
 * Seebeck coefficients at 25C (uV/C) per ANSI MC96.1 / IEC 60584.
 * Reference: NIST ITS-90 Thermocouple Database
 * ================================================================== */

typedef enum {
    TC_TYPE_K = 0,   /* Chromel-Alumel, -270 to 1372 C, 40.6 uV/C at 25C */
    TC_TYPE_J = 1,   /* Iron-Constantan, -210 to 1200 C, 51.9 uV/C */
    TC_TYPE_T = 2,   /* Copper-Constantan, -270 to 400 C, 40.6 uV/C */
    TC_TYPE_E = 3,   /* Chromel-Constantan, -270 to 1000 C, 60.9 uV/C */
    TC_TYPE_N = 4,   /* Nicrosil-Nisil, -270 to 1300 C, 26.0 uV/C */
    TC_TYPE_R = 5,   /* Pt13%Rh-Pt, -50 to 1768 C, 10.4 uV/C */
    TC_TYPE_S = 6,   /* Pt10%Rh-Pt, -50 to 1768 C, 9.4 uV/C */
    TC_TYPE_B = 7    /* Pt30%Rh-Pt6%Rh, 0 to 1820 C, 5.3 uV/C */
} thermocouple_type_t;

#ifdef __cplusplus
}
#endif

#endif /* INST_AMP_DEFS_H */
