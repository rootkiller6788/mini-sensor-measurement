#ifndef SMART_SENSOR_H
#define SMART_SENSOR_H
#include <stddef.h>
#include <stdint.h>
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1: Core Definitions — Smart Sensor Architecture & Sensor Types
 *
 * A smart sensor integrates: sensing element + signal conditioning +
 * ADC + microcontroller + digital communication interface per
 * IEEE 1451.0 standard for smart transducer interfaces.
 *
 * Reference: Frank, R. (2013), "Understanding Smart Sensors", 3rd ed., Artech House
 *            IEEE 1451.0-2007 Standard for Smart Transducer Interface
 * ========================================================================= */

/** Physical transduction principle classification */
typedef enum {
    SS_TRANSDUCER_RESISTIVE      = 0,   /* Strain gauge, thermistor, photoresistor */
    SS_TRANSDUCER_CAPACITIVE     = 1,   /* MEMS accelerometer, pressure, humidity */
    SS_TRANSDUCER_INDUCTIVE      = 2,   /* LVDT, eddy current proximity */
    SS_TRANSDUCER_PIEZOELECTRIC  = 3,   /* Vibration, ultrasound, dynamic pressure */
    SS_TRANSDUCER_THERMOELECTRIC = 4,   /* Thermocouple (Seebeck effect) */
    SS_TRANSDUCER_PHOTOELECTRIC  = 5,   /* Photodiode, phototransistor */
    SS_TRANSDUCER_HALL_EFFECT    = 6,   /* Magnetic field, current sensing */
    SS_TRANSDUCER_ELECTROCHEMICAL= 7,   /* Gas sensor, pH, biosensor */
    SS_TRANSDUCER_MEMS_THERMAL   = 8,   /* MEMS thermal (flow, vacuum) */
    SS_TRANSDUCER_MEMS_RESONANT  = 9,   /* MEMS resonant (pressure, mass) */
    SS_TRANSDUCER_OPTICAL        = 10,  /* Fiber optic, spectroscopy */
    SS_TRANSDUCER_ULTRASONIC     = 11,  /* Distance, flow, imaging */
    SS_TRANSDUCER_RADIATION      = 12,  /* Geiger-Müller, scintillation */
    SS_TRANSDUCER_MAGNETORESISTIVE=13   /* AMR, GMR, TMR sensors */
} ss_transducer_t;

/** Measurand (physical quantity being measured) categories */
typedef enum {
    SS_MEASURAND_TEMPERATURE     = 0,
    SS_MEASURAND_PRESSURE        = 1,
    SS_MEASURAND_HUMIDITY        = 2,
    SS_MEASURAND_ACCELERATION    = 3,
    SS_MEASURAND_ANGULAR_VELOCITY= 4,
    SS_MEASURAND_MAGNETIC_FIELD  = 5,
    SS_MEASURAND_DISPLACEMENT    = 6,
    SS_MEASURAND_STRAIN          = 7,
    SS_MEASURAND_FORCE_TORQUE    = 8,
    SS_MEASURAND_FLOW_RATE       = 9,
    SS_MEASURAND_LIGHT_INTENSITY = 10,
    SS_MEASURAND_GAS_CONCENTRATION=11,
    SS_MEASURAND_pH_CHEMICAL     = 12,
    SS_MEASURAND_PROXIMITY       = 13,
    SS_MEASURAND_SOUND_PRESSURE  = 14,
    SS_MEASURAND_VOLTAGE_CURRENT = 15,
    SS_MEASURAND_ANGLE_TILT      = 16
} ss_measurand_t;

/** Digital communication interface for smart sensors */
typedef enum {
    SS_INTERFACE_ANALOG_VOLTAGE  = 0,   /* 0-5V, 0-10V */
    SS_INTERFACE_ANALOG_CURRENT  = 1,   /* 4-20 mA current loop */
    SS_INTERFACE_PWM             = 2,   /* Pulse width modulation */
    SS_INTERFACE_FREQUENCY       = 3,   /* Frequency modulated output */
    SS_INTERFACE_I2C             = 10,  /* Inter-Integrated Circuit */
    SS_INTERFACE_SPI             = 11,  /* Serial Peripheral Interface */
    SS_INTERFACE_UART            = 12,  /* Universal Asynchronous Rx/Tx */
    SS_INTERFACE_ONEWIRE         = 13,  /* Dallas/Maxim 1-Wire */
    SS_INTERFACE_CAN             = 14,  /* Controller Area Network */
    SS_INTERFACE_MODBUS          = 15,  /* Modbus RTU/TCP */
    SS_INTERFACE_HART            = 16,  /* Highway Addressable Remote Transducer */
    SS_INTERFACE_WIRELESS_BLE    = 20,  /* Bluetooth Low Energy */
    SS_INTERFACE_WIRELESS_ZIGBEE = 21,  /* Zigbee/IEEE 802.15.4 */
    SS_INTERFACE_WIRELESS_LORA   = 22,  /* LoRa/LoRaWAN */
    SS_INTERFACE_WIRELESS_WIFI   = 23   /* WiFi 802.11 */
} ss_interface_t;

/** Sensor operating mode (IEEE 1451.0 operational states) */
typedef enum {
    SS_MODE_IDLE                 = 0,
    SS_MODE_NORMAL_SAMPLING      = 1,
    SS_MODE_BURST_SAMPLING       = 2,
    SS_MODE_SELF_TEST            = 3,
    SS_MODE_CALIBRATION          = 4,
    SS_MODE_SLEEP_LOW_POWER      = 5,
    SS_MODE_FAULT_DETECTED       = 6,
    SS_MODE_MAINTENANCE          = 7
} ss_operating_mode_t;

/** IEEE 1451.4 TEDS (Transducer Electronic Data Sheet) template types */
typedef enum {
    SS_TEDS_BASIC                = 0,
    SS_TEDS_STANDARD             = 1,
    SS_TEDS_CALIBRATION          = 2,
    SS_TEDS_USER_DEFINED         = 3,
    SS_TEDS_MANUFACTURER_EXT     = 4
} ss_teds_template_t;

/* =========================================================================
 * L1: Core Data Structures
 * ========================================================================= */

/**
 * @brief Sensor physical characteristics specification
 *
 * Captures the static transfer characteristic of a sensor element.
 * All values in SI units or derived units as specified.
 */
typedef struct {
    ss_transducer_t  transducer_type;
    ss_measurand_t   measurand;
    double           sensitivity;      /* Output change / input change (e.g., mV/degC) */
    double           offset;           /* Zero-input output value */
    double           full_scale_input; /* Maximum rated input (e.g., 100 degC, 10 bar) */
    double           full_scale_output;/* Output at full-scale input */
    double           resolution;       /* Smallest detectable input change */
    double           noise_density;    /* Noise spectral density (units/rtHz) */
    double           bandwidth;        /* -3 dB bandwidth in Hz */
    double           input_impedance;  /* Electrical input impedance (ohm) */
    double           output_impedance; /* Electrical output impedance (ohm) */
    double           supply_voltage_min;
    double           supply_voltage_max;
    double           supply_current_typ;
    double           temperature_min;   /* Operating temp range min (degC) */
    double           temperature_max;   /* Operating temp range max (degC) */
} ss_sensor_spec_t;

/**
 * @brief Accuracy and precision metrics for a sensor
 *
 * Statistical characterization per JCGM 100:2008 (GUM — Guide to
 * the Expression of Uncertainty in Measurement).
 */
typedef struct {
    double accuracy;          /* Max deviation from true value (pct FS or abs) */
    double precision;         /* Repeatability: std dev of repeated measurements */
    double linearity_error;   /* Max deviation from best-fit straight line (pct FS) */
    double hysteresis_error;  /* Max difference upscale vs. downscale (pct FS) */
    double drift_24h;         /* Short-term drift over 24 hours */
    double drift_long_term;   /* Long-term drift over 1 year */
    double temperature_coeff_zero;  /* Zero drift per degC */
    double temperature_coeff_span;  /* Sensitivity drift per degC */
    double response_time;     /* 10pct-90pct step response time (seconds) */
    double settling_time;     /* Time to settle within 1pct of final value */
} ss_accuracy_spec_t;

/**
 * @brief Smart sensor runtime state
 *
 * Maintains the operational context of a live smart sensor instance,
 * including current readings and diagnostic status.
 */
typedef struct {
    ss_operating_mode_t mode;
    uint64_t            sample_count;       /* Total samples acquired */
    uint64_t            fault_count;        /* Non-recoverable fault counter */
    uint64_t            recovery_count;     /* Self-recovery event counter */
    double              current_raw_value;  /* Most recent raw ADC reading */
    double              current_processed_value;  /* After conditioning */
    double              current_engineering_units; /* Final calibrated value */
    uint64_t            last_sample_timestamp_us;
    double              moving_average;     /* Short-term running average */
    double              moving_variance;    /* Short-term running variance */
    int                 self_test_passed;   /* 1 = passed, 0 = failed, -1 = not run */
    int                 calibration_valid;  /* 1 = valid, 0 = expired/invalid */
    double              supply_voltage_monitor; /* Current supply voltage reading */
    double              internal_temperature;   /* On-board temperature sensor */
} ss_sensor_state_t;

/**
 * @brief IEEE 1451.4 TEDS basic identification record
 *
 * Transducer Electronic Data Sheet — stores sensor identity,
 * calibration, and configuration in non-volatile memory per
 * IEEE 1451.4 standard for mixed-mode interface.
 */
typedef struct {
    ss_teds_template_t  template_type;
    char                manufacturer_name[32];
    char                model_number[32];
    char                serial_number[32];
    char                version_letter[8];
    uint16_t            manufacture_year;
    uint16_t            manufacture_month;
    uint32_t            teds_data_length;
    /* TEDS calibration data */
    double              cal_sensitivity;     /* Calibrated sensitivity */
    double              cal_offset;          /* Calibrated offset */
    double              cal_date;            /* Julian date of calibration */
    double              cal_interval_days;   /* Recommended recalibration interval */
    uint8_t             cal_polynomial_order;/* Order of calibration polynomial */
    double              cal_coeffs[8];       /* Polynomial coefficients (highest first) */
} ss_teds_t;

/**
 * @brief Measurement sample with timestamp and quality metadata
 *
 * Implements a time-stamped measurement record suitable for
 * time-series analysis and sensor data logging.
 */
typedef struct {
    uint64_t timestamp_us;       /* Microsecond timestamp */
    double   raw_value;          /* Unprocessed ADC reading */
    double   conditioned_value;  /* After analog/digital conditioning */
    double   calibrated_value;   /* Final engineering units */
    double   uncertainty_k2;     /* Expanded uncertainty (k=2, ~95pct CI) */
    uint32_t flags;              /* Quality flags (sensor_status bitfield) */
    double   internal_temp;      /* Internal temperature at time of reading */
    double   supply_voltage;     /* Supply voltage at time of reading */
} ss_measurement_t;

/* Quality flag bit definitions */
#define SS_FLAG_VALID            0x0001  /* Measurement is valid */
#define SS_FLAG_OVERRANGE        0x0002  /* Input exceeds sensor range */
#define SS_FLAG_UNDERRANGE       0x0004  /* Input below sensor range */
#define SS_FLAG_NOISE_DETECTED   0x0008  /* Excessive noise detected */
#define SS_FLAG_DRIFT_WARNING    0x0010  /* Calibration drift warning */
#define SS_FLAG_SELF_TEST_FAIL   0x0020  /* Self-test failed on this reading */
#define SS_FLAG_TRANSIENT        0x0040  /* Transient/settling condition */
#define SS_FLAG_COMPENSATED      0x0080  /* Temperature compensated */
#define SS_FLAG_FUSED            0x0100  /* Sensor fusion applied */
#define SS_FLAG_INTERPOLATED     0x0200  /* Missing sample — interpolated */
#define SS_FLAG_FAULT_MASKED     0x0400  /* Fault masked by redundancy */

/* =========================================================================
 * L2: Core Concepts — Signal Conditioning Configuration
 * ========================================================================= */

/**
 * @brief Analog front-end signal conditioning parameters
 *
 * Models the complete analog signal chain between the transducer
 * and the ADC: excitation source, amplifier, and anti-aliasing filter.
 */
typedef struct {
    /* Excitation / bias */
    double   excitation_voltage;    /* Bridge excitation voltage (V) */
    double   excitation_current;    /* Constant current excitation (mA) */
    uint8_t  excitation_type;       /* 0=voltage, 1=current, 2=AC, 3=none */

    /* Amplifier stage */
    double   amplifier_gain;        /* Voltage gain (V/V) */
    double   amplifier_offset;      /* Input-referred offset voltage (uV) */
    double   amplifier_noise;       /* Input-referred noise (nV/rtHz) */
    double   amplifier_bandwidth;   /* -3 dB bandwidth (Hz) */
    uint8_t  amplifier_topology;    /* 0=instr.amp, 1=diff, 2=charge, 3=TIA */

    /* Anti-aliasing filter */
    uint8_t  filter_order;          /* Filter order (1-8) */
    double   filter_cutoff_freq;    /* -3 dB cutoff frequency (Hz) */
    uint8_t  filter_type;           /* 0=Butterworth, 1=Bessel, 2=Chebyshev, 3=Elliptic */
    double   filter_ripple_db;      /* Passband ripple for Chebyshev/Elliptic (dB) */

    /* ADC parameters */
    uint8_t  adc_resolution_bits;   /* ADC resolution (8-24 bits) */
    double   adc_reference_voltage; /* ADC reference voltage (V) */
    double   adc_sample_rate_hz;    /* ADC sampling rate (samples/sec) */
    uint8_t  adc_oversampling_ratio;/* Oversampling ratio (OSR) for increased ENOB */
} ss_analog_frontend_t;

/**
 * @brief Digital post-processing configuration
 *
 * Digital signal processing chain applied after ADC conversion
 * to improve measurement quality.
 */
typedef struct {
    /* Digital filtering */
    uint8_t  moving_avg_window;      /* Moving average window length (samples) */
    uint8_t  median_filter_window;   /* Median filter window (1=disabled, 3/5/7) */
    double   lowpass_alpha;          /* 1st-order IIR lowpass coefficient (0-1) */
    double   highpass_alpha;         /* 1st-order IIR highpass coefficient (0-1) */

    /* Threshold detection */
    double   upper_threshold;        /* Upper alarm threshold (engineering units) */
    double   lower_threshold;        /* Lower alarm threshold (engineering units) */
    double   rate_of_change_limit;   /* Max allowed rate of change (units/sec) */
    uint32_t threshold_debounce_ms;  /* Debounce time for threshold crossing (ms) */

    /* Averaging / decimation */
    uint16_t decimation_factor;      /* Reduce sample rate by factor N */
    uint8_t  averaging_mode;         /* 0=simple, 1=running, 2=exponential, 3=boxcar */
} ss_digital_processing_t;

/**
 * @brief Smart sensor complete configuration
 *
 * Aggregates all configuration parameters for a complete smart sensor
 * channel including the sensing element, analog front-end, ADC,
 * digital processing, and communication interface.
 */
typedef struct {
    ss_sensor_spec_t       sensor_spec;
    ss_accuracy_spec_t     accuracy_spec;
    ss_analog_frontend_t   analog_frontend;
    ss_digital_processing_t digital_proc;
    ss_teds_t              teds;
    ss_interface_t         interface_type;
    uint8_t                interface_address;  /* I2C addr / SPI CS / UART ID */
    uint32_t               interface_baud_rate;/* Communication speed */
    uint16_t               sample_period_ms;    /* Sampling period in milliseconds */
    uint8_t                enable_temp_comp;    /* Enable temperature compensation */
    uint8_t                enable_self_test;    /* Enable periodic self-test */
    uint16_t               self_test_period_s;  /* Self-test period in seconds */
} ss_sensor_config_t;

/* =========================================================================
 * L3: Mathematical Structures — Sensor Transfer Function Models
 * ========================================================================= */

/**
 * @brief Linear sensor model: y = sensitivity * x + offset
 *
 * Implements the simplest sensor model characterizing a sensor
 * with a linear transfer function (most common model).
 */
typedef struct {
    double sensitivity;    /* Slope (gain) of the transfer function */
    double offset;         /* y-intercept (zero offset) */
    double r_squared;      /* Coefficient of determination from calibration */
    double residuals_rms;  /* RMS residual error from calibration fit */
} ss_linear_model_t;

/**
 * @brief Polynomial sensor model: y = a_n*x^n + ... + a_1*x + a_0
 *
 * Supports up to 7th-order polynomial calibration for sensors
 * with nonlinear transfer characteristics (e.g., thermistors,
 * thermocouples, nonlinear pressure sensors).
 */
typedef struct {
    uint8_t  order;           /* Polynomial order (1-7) */
    double   coeffs[8];       /* Coefficients [a_n, ..., a_1, a_0], highest first */
    double   x_min;           /* Minimum valid input for this calibration */
    double   x_max;           /* Maximum valid input for this calibration */
    double   r_squared;       /* Coefficient of determination */
    double   max_residual;    /* Maximum absolute residual error */
} ss_polynomial_model_t;

/**
 * @brief Steinhart-Hart thermistor model
 *
 * Accurate thermistor temperature model (typically +/-0.01 degC accuracy).
 * 1/T = A + B*ln(R) + C*[ln(R)]^3
 * where T is in Kelvin, R is the thermistor resistance in Ohms.
 *
 * Reference: Steinhart, J.S. & Hart, S.R. (1968), "Calibration curves
 *            for thermistors", Deep Sea Research, 15(4):497-503
 */
typedef struct {
    double coeff_A;          /* Steinhart-Hart A coefficient */
    double coeff_B;          /* Steinhart-Hart B coefficient */
    double coeff_C;          /* Steinhart-Hart C coefficient */
    double r_at_25c;         /* Reference resistance at 25 degC */
    double beta_value;       /* Beta value for simplified model */
} ss_steinhart_hart_t;

/**
 * @brief Thermocouple voltage-to-temperature model (ITS-90)
 *
 * Implements the NIST ITS-90 thermocouple inverse model.
 * EMF (mV) -> Temperature (degC) using polynomial coefficients
 * from the NIST Monograph 175 database.
 *
 * Reference: NIST ITS-90 Thermocouple Database (NIST Monograph 175)
 */
typedef struct {
    uint8_t  thermocouple_type;  /* 0=K, 1=J, 2=T, 3=E, 4=N, 5=R, 6=S, 7=B */
    double   cold_junction_temp; /* Reference junction temperature (degC) */
    double   coeffs[11];         /* NIST inverse polynomial coefficients */
    uint8_t  n_coeffs;           /* Number of valid coefficients */
    double   temp_min;           /* Valid temperature range minimum */
    double   temp_max;           /* Valid temperature range maximum */
} ss_thermocouple_model_t;

/**
 * @brief Allan variance statistical structure
 *
 * Characterizes sensor noise and stability as a function of
 * averaging time. Used for inertial sensors (gyroscopes,
 * accelerometers) to separate noise sources.
 *
 * Reference: IEEE Std 952-1997, "Standard Specification Format Guide
 *            and Test Procedure for Single-Axis Interferometric
 *            Fiber Optic Gyros"
 */
typedef struct {
    double   angle_random_walk;    /* ARW (units/rthr) — white noise */
    double   bias_instability;     /* Bias instability minimum (units) */
    double   rate_random_walk;     /* RRW (units/hr/rthr) — flicker noise */
    double   quantization_noise;   /* Quantization noise */
    double   rate_ramp;            /* Drift rate ramp */
    double   correlation_time;     /* Correlation time at bias instability minimum */
    uint32_t n_tau_points;         /* Number of averaging time points */
    double   tau_values[64];       /* Averaging times (seconds) */
    double   allan_deviation[64];  /* Allan deviation at each tau */
} ss_allan_variance_t;

/* =========================================================================
 * L2: Core Concepts — Sensor Data Structures for Aggregation
 * ========================================================================= */

/**
 * @brief Sensor data buffer for time-series storage
 */
typedef struct {
    ss_measurement_t *buffer;
    size_t            capacity;
    size_t            head;       /* Write position */
    size_t            count;      /* Number of valid entries */
    int               wrap_occurred;
} ss_measurement_buffer_t;

/**
 * @brief Multi-channel sensor aggregator
 *
 * Manages multiple sensor channels in a single smart sensor node
 * for sensor fusion and multi-parameter monitoring.
 */
#define SS_MAX_CHANNELS 16

typedef struct {
    ss_sensor_config_t      channels[SS_MAX_CHANNELS];
    ss_sensor_state_t       states[SS_MAX_CHANNELS];
    ss_measurement_buffer_t buffers[SS_MAX_CHANNELS];
    uint8_t             n_channels;
    uint8_t             active_channel_mask;   /* Bitmask of active channels */
    uint64_t            node_id;               /* Unique node identifier (e.g., EUI-64) */
    double              timestamp_offset_s;    /* Global time offset for synchronization */
} ss_node_t;

/* =========================================================================
 * API: Initialization and Configuration
 * ========================================================================= */

/**
 * @brief Initialize a sensor specification with sensible defaults
 * @param spec    Output sensor specification
 * @param ttype   Transducer type
 * @param meas    Measurand type
 * @param sens    Nominal sensitivity (output/input)
 * @param fsi     Full-scale input
 * @param fso     Full-scale output at FSI
 */
void ss_sensor_spec_init(ss_sensor_spec_t *spec, ss_transducer_t ttype,
                         ss_measurand_t meas, double sens,
                         double fsi, double fso);

/**
 * @brief Initialize accuracy specification with typical values
 * @param acc    Output accuracy spec
 * @param acc_percent  Accuracy as pct of full scale
 * @param precision_repeat   Repeatability (std dev)
 */
void ss_accuracy_spec_init(ss_accuracy_spec_t *acc, double acc_percent,
                           double precision_repeat);

/**
 * @brief Initialize analog front-end with default values
 * @param afe   Output analog front-end config
 */
void ss_analog_frontend_init(ss_analog_frontend_t *afe);

/**
 * @brief Initialize digital processing with default values
 * @param dp    Output digital processing config
 */
void ss_digital_processing_init(ss_digital_processing_t *dp);

/**
 * @brief Initialize a complete sensor config from sub-components
 * @param cfg   Output sensor config
 * @param spec  Sensor specification
 * @param acc   Accuracy specification
 * @param afe   Analog front-end
 * @param dp    Digital processing
 * @param iface Communication interface type
 */
void ss_sensor_config_init(ss_sensor_config_t *cfg,
                           const ss_sensor_spec_t *spec,
                           const ss_accuracy_spec_t *acc,
                           const ss_analog_frontend_t *afe,
                           const ss_digital_processing_t *dp,
                           ss_interface_t iface);

/**
 * @brief Initialize a sensor node with given number of channels
 * @param node       Output node structure
 * @param n_channels Number of sensor channels (1..SS_MAX_CHANNELS)
 * @param node_id    Unique node identifier
 * @return 0 on success, -1 on parameter error
 */
int ss_node_init(ss_node_t *node, uint8_t n_channels, uint64_t node_id);

/**
 * @brief Add a sensor channel to an existing node
 * @param node  Node to modify
 * @param cfg   Sensor configuration for the new channel
 * @return      Channel index (0..n_channels-1), or -1 on failure
 */
int ss_node_add_channel(ss_node_t *node, const ss_sensor_config_t *cfg);

/**
 * @brief Read sensor state fields
 * @param state Sensor state structure
 * @param mode  Output: current operating mode
 * @param raw   Output: current raw value
 * @param eng   Output: current engineering units value
 */
void ss_sensor_state_read(const ss_sensor_state_t *state,
                          ss_operating_mode_t *mode,
                          double *raw, double *eng);

/**
 * @brief Set sensor operating mode with transition checks
 * @param state  Sensor state to modify
 * @param mode   New operating mode
 * @return 0 on success, -1 if invalid transition
 */
int ss_sensor_state_set_mode(ss_sensor_state_t *state, ss_operating_mode_t mode);

/**
 * @brief Update sensor state with a new measurement
 * @param state  Sensor state to update
 * @param raw    New raw ADC reading
 * @param proc   Processed (conditioned) value
 * @param eng    Engineering units value
 * @param ts_us  Timestamp in microseconds
 */
void ss_sensor_state_update(ss_sensor_state_t *state, double raw,
                            double proc, double eng, uint64_t ts_us);

/**
 * @brief Update the running statistics (moving average and variance)
 *
 * Uses Welford's online algorithm for numerically stable
 * single-pass variance computation.
 *
 * Reference: Welford, B.P. (1962), "Note on a method for calculating
 *            corrected sums of squares and products", Technometrics 4(3):419-420
 *
 * @param state   Sensor state with running statistics
 * @param new_val New measurement value
 */
void ss_sensor_state_update_stats(ss_sensor_state_t *state, double new_val);

/* =========================================================================
 * API: Measurement Buffer Operations
 * ========================================================================= */

/**
 * @brief Initialize a measurement buffer
 * @param buf       Buffer structure to initialize
 * @param buffer    Pre-allocated measurement array
 * @param capacity  Number of measurements the buffer can hold
 */
void ss_measurement_buffer_init(ss_measurement_buffer_t *buf,
                                ss_measurement_t *buffer, size_t capacity);

/**
 * @brief Push a measurement into the ring buffer
 * @param buf    Ring buffer
 * @param meas   Measurement to push
 * @return       Number of valid entries after push
 */
size_t ss_measurement_buffer_push(ss_measurement_buffer_t *buf,
                                  const ss_measurement_t *meas);

/**
 * @brief Get a measurement by logical index (0 = oldest in buffer)
 * @param buf    Ring buffer
 * @param index  Logical index (0-based from oldest)
 * @param meas   Output measurement
 * @return       0 on success, -1 if index out of range
 */
int ss_measurement_buffer_get(const ss_measurement_buffer_t *buf,
                              size_t index, ss_measurement_t *meas);

/**
 * @brief Compute mean of all valid measurements in buffer
 * @param buf  Measurement buffer
 * @return     Mean value in engineering units, NaN if empty
 */
double ss_measurement_buffer_mean(const ss_measurement_buffer_t *buf);

/**
 * @brief Compute standard deviation of measurements in buffer
 * @param buf  Measurement buffer
 * @return     Standard deviation, NaN if fewer than 2 samples
 */
double ss_measurement_buffer_stddev(const ss_measurement_buffer_t *buf);

/**
 * @brief Find minimum and maximum values in buffer
 * @param buf  Measurement buffer
 * @param min  Output minimum
 * @param max  Output maximum
 * @return     0 on success, -1 if buffer empty
 */
int ss_measurement_buffer_range(const ss_measurement_buffer_t *buf,
                                double *min, double *max);

/**
 * @brief Clear all measurements from buffer
 * @param buf  Buffer to clear
 */
void ss_measurement_buffer_clear(ss_measurement_buffer_t *buf);

#ifdef __cplusplus
}
#endif
#endif /* SMART_SENSOR_H */
