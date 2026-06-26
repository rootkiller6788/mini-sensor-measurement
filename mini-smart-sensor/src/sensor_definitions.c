/**
 * @file    sensor_definitions.c
 * @brief   Smart sensor core definitions and initialization
 *
 * Knowledge Coverage:
 *   L1 (Definitions): Sensor types, measurand categories, interface types,
 *                     operating modes, TEDS structure, measurement records
 *   L2 (Core Concepts): Sensor configuration, measurement buffer management,
 *                       sensor state machine transitions
 *
 * Reference:
 *   IEEE 1451.0-2007 — Smart Transducer Interface Standard
 *   IEEE 1451.4-2004 — Mixed-Mode Communication Protocols and TEDS Formats
 *   Frank, R. (2013), "Understanding Smart Sensors", 3rd ed., Artech House
 *   JCGM 100:2008 — Evaluation of measurement data (GUM)
 *
 * Course Mapping:
 *   Berkeley EE16A/B — Designing Information Devices and Systems (sensors intro)
 *   MIT 6.003 — Signals and Systems (measurement systems)
 *   Michigan EECS 351 — DSP (sensor data processing)
 */

#include "smart_sensor.h"
#include <string.h>
#include <stdio.h>
#include <math.h>

/* =========================================================================
 * L1: Sensor Specification Initialization
 * ========================================================================= */

/**
 * @brief Initialize sensor specification with sensible defaults
 *
 * Sets up a sensor_spec structure with the core identifying parameters
 * and derives reasonable defaults for other fields based on the
 * transducer type. The sensitivity is the key parameter linking
 * the physical world to the electrical domain:
 *
 *   V_out = sensitivity * physical_input + offset
 *
 * @param spec   Output specification structure
 * @param ttype  Transducer classification
 * @param meas   Physical quantity being measured
 * @param sens   Sensitivity (output_units / input_units)
 * @param fsi    Full-scale input (max rated physical input)
 * @param fso    Full-scale output (output at full-scale input)
 */
void ss_sensor_spec_init(ss_sensor_spec_t *spec, ss_transducer_t ttype,
                         ss_measurand_t meas, double sens,
                         double fsi, double fso) {
    if (!spec) return;

    memset(spec, 0, sizeof(*spec));
    spec->transducer_type = ttype;
    spec->measurand = meas;
    spec->sensitivity = sens;
    spec->full_scale_input = fsi;
    spec->full_scale_output = fso;

    /* Derive offset assuming linear model: output = sens * input + offset
     * Default: zero output at zero input -> offset = 0 */
    spec->offset = 0.0;

    /* Estimate resolution: assume 12-bit effective ADC resolution
     * resolution = full_scale_input / 2^12 */
    if (fsi > 0.0) {
        spec->resolution = fsi / 4096.0;
    } else {
        spec->resolution = 1e-3; /* 1 mV/LSB default */
    }

    /* Default bandwidth based on transducer type */
    switch (ttype) {
    case SS_TRANSDUCER_PIEZOELECTRIC:
    case SS_TRANSDUCER_ULTRASONIC:
        spec->bandwidth = 100000.0; /* High bandwidth for dynamic sensors */
        break;
    case SS_TRANSDUCER_MEMS_THERMAL:
    case SS_TRANSDUCER_THERMOELECTRIC:
        spec->bandwidth = 1.0;  /* Slow thermal response */
        break;
    default:
        spec->bandwidth = 1000.0; /* General-purpose */
        break;
    }

    /* Input impedance defaults by transducer type */
    switch (ttype) {
    case SS_TRANSDUCER_PIEZOELECTRIC:
        spec->input_impedance = 1e9;  /* High Z for charge output */
        spec->output_impedance = 100.0;
        break;
    case SS_TRANSDUCER_RESISTIVE:
        spec->input_impedance = 350.0; /* Typical strain gauge */
        spec->output_impedance = 350.0;
        break;
    case SS_TRANSDUCER_THERMOELECTRIC:
        spec->input_impedance = 10.0;  /* Low Z for thermocouple */
        spec->output_impedance = 10.0;
        break;
    default:
        spec->input_impedance = 10000.0;
        spec->output_impedance = 1000.0;
        break;
    }

    /* Default noise density: thermal noise floor for the impedance */
    {
        double k = 1.380649e-23; /* Boltzmann constant */
        double T = 298.0;        /* Room temperature (K) */
        spec->noise_density = sqrt(4.0 * k * T * spec->output_impedance);
    }

    /* Supply voltage defaults (3.3V typical for digital sensors) */
    spec->supply_voltage_min = 2.7;
    spec->supply_voltage_max = 3.6;
    spec->supply_current_typ = 1.0; /* mA */

    /* Temperature range: industrial (-40 to +85 degC) */
    spec->temperature_min = -40.0;
    spec->temperature_max = 85.0;
}

/* =========================================================================
 * L1: Accuracy Specification Initialization
 * ========================================================================= */

/**
 * @brief Initialize accuracy specification
 *
 * Sets up accuracy metrics based on percent-of-full-scale accuracy
 * and repeatability. Derives secondary metrics from these primary ones.
 *
 * The relationship between the metrics:
 *   Total error <= |accuracy| + |linearity_error| + temperature_effects
 *
 * @param acc              Output accuracy specification
 * @param acc_percent      Accuracy as percentage of full scale
 * @param precision_repeat Repeatability standard deviation
 */
void ss_accuracy_spec_init(ss_accuracy_spec_t *acc, double acc_percent,
                           double precision_repeat) {
    if (!acc) return;

    memset(acc, 0, sizeof(*acc));

    acc->accuracy = acc_percent;
    acc->precision = precision_repeat;

    /* Linearity error: typically 1/3 to 1/2 of total accuracy */
    acc->linearity_error = acc_percent * 0.4;

    /* Hysteresis: typically smaller than linearity for electronic sensors */
    acc->hysteresis_error = acc_percent * 0.2;

    /* Drift estimates: 24-hour drift ~ repeatability */
    acc->drift_24h = precision_repeat * 2.0;

    /* Long-term drift: typically 2-5x the 24h drift per year */
    acc->drift_long_term = precision_repeat * 10.0;

    /* Temperature coefficients: 0.01% FS per degree C typical */
    acc->temperature_coeff_zero = acc_percent * 0.0001;
    acc->temperature_coeff_span = acc_percent * 0.0002;

    /* Response time defaults */
    acc->response_time = 0.001;  /* 1 ms typical for electronic sensors */
    acc->settling_time = 0.010;  /* 10 ms to settle within 1% */
}

/* =========================================================================
 * L2: Analog Front-End Initialization
 * ========================================================================= */

/**
 * @brief Initialize analog front-end with sensible defaults
 *
 * Configures a typical analog signal chain for a sensor:
 *   Sensor -> Instrumentation Amplifier -> Anti-aliasing Filter -> ADC
 *
 * Defaults target a general-purpose measurement setup with:
 *   - Bridge excitation for resistive sensors (3.3V)
 *   - Moderate gain (100x) for mV-level signals
 *   - 2nd-order Butterworth anti-aliasing at 1kHz
 *   - 12-bit ADC at 10kSPS
 *
 * @param afe  Output analog front-end configuration
 */
void ss_analog_frontend_init(ss_analog_frontend_t *afe) {
    if (!afe) return;

    memset(afe, 0, sizeof(*afe));

    /* Default: voltage excitation for resistive sensors */
    afe->excitation_type = 0;  /* Voltage excitation */
    afe->excitation_voltage = 3.3;
    afe->excitation_current = 1.0;

    /* Instrumentation amplifier with moderate gain */
    afe->amplifier_topology = 0;  /* Instrumentation amplifier */
    afe->amplifier_gain = 100.0;
    afe->amplifier_offset = 50.0;  /* 50 uV typical for good INA */
    afe->amplifier_noise = 10.0;   /* 10 nV/rtHz (e.g., AD8429) */
    afe->amplifier_bandwidth = 100000.0; /* 100 kHz GBW at G=100 */

    /* 2nd-order Butterworth anti-aliasing filter */
    afe->filter_type = 0;  /* Butterworth: maximally flat passband */
    afe->filter_order = 2;
    afe->filter_cutoff_freq = 1000.0; /* 1 kHz */
    afe->filter_ripple_db = 0.0;

    /* 12-bit ADC at 10 kSPS (general purpose) */
    afe->adc_resolution_bits = 12;
    afe->adc_reference_voltage = 3.3;
    afe->adc_sample_rate_hz = 10000.0;
    afe->adc_oversampling_ratio = 1;
}

/* =========================================================================
 * L2: Digital Processing Initialization
 * ========================================================================= */

/**
 * @brief Initialize digital processing with conservative defaults
 *
 * Sets up digital signal processing chain with:
 *   - 4-sample moving average for noise reduction
 *   - Disabled median filter (no impulse noise expected by default)
 *   - Low-pass IIR with alpha=0.1 for mild smoothing
 *   - No threshold detection enabled
 *   - No decimation
 *
 * @param dp  Output digital processing configuration
 */
void ss_digital_processing_init(ss_digital_processing_t *dp) {
    if (!dp) return;

    memset(dp, 0, sizeof(*dp));

    /* Light moving average: 4 samples, balanced noise/response */
    dp->moving_avg_window = 4;

    /* Median filter disabled by default */
    dp->median_filter_window = 1;

    /* Mild IIR lowpass: alpha=0.1, cutoff ~0.016*fs */
    dp->lowpass_alpha = 0.1;
    dp->highpass_alpha = 0.0;  /* No highpass filtering */

    /* Thresholds: disabled (set to extreme values) */
    dp->upper_threshold = 1e308;
    dp->lower_threshold = -1e308;
    dp->rate_of_change_limit = 1e308;
    dp->threshold_debounce_ms = 0;

    /* No decimation */
    dp->decimation_factor = 1;
    dp->averaging_mode = 1;  /* Running average */
}

/* =========================================================================
 * L1: Complete Sensor Configuration Initialization
 * ========================================================================= */

/**
 * @brief Initialize a complete sensor configuration
 *
 * Assembles all sub-configurations into a unified sensor config.
 * This represents a complete smart sensor channel ready for operation.
 *
 * All sub-structures are deep-copied so the caller retains ownership
 * of the input structures.
 *
 * @param cfg   Output complete sensor configuration
 * @param spec  Sensor specification (required)
 * @param acc   Accuracy specification (can be NULL for defaults)
 * @param afe   Analog front-end (can be NULL for defaults)
 * @param dp    Digital processing (can be NULL for defaults)
 * @param iface Communication interface type
 */
void ss_sensor_config_init(ss_sensor_config_t *cfg,
                           const ss_sensor_spec_t *spec,
                           const ss_accuracy_spec_t *acc,
                           const ss_analog_frontend_t *afe,
                           const ss_digital_processing_t *dp,
                           ss_interface_t iface) {
    if (!cfg || !spec) return;

    memset(cfg, 0, sizeof(*cfg));

    /* Copy sensor specification */
    memcpy(&cfg->sensor_spec, spec, sizeof(*spec));

    /* Copy or default accuracy */
    if (acc) {
        memcpy(&cfg->accuracy_spec, acc, sizeof(*acc));
    } else {
        ss_accuracy_spec_init(&cfg->accuracy_spec, 1.0, 0.01);
    }

    /* Copy or default analog front-end */
    if (afe) {
        memcpy(&cfg->analog_frontend, afe, sizeof(*afe));
    } else {
        ss_analog_frontend_init(&cfg->analog_frontend);
    }

    /* Copy or default digital processing */
    if (dp) {
        memcpy(&cfg->digital_proc, dp, sizeof(*dp));
    } else {
        ss_digital_processing_init(&cfg->digital_proc);
    }

    /* Set interface */
    cfg->interface_type = iface;
    cfg->interface_address = 0;
    cfg->interface_baud_rate = (iface == SS_INTERFACE_I2C) ? 400000 : 115200;

    /* Default sampling: 100 ms (10 Hz) */
    cfg->sample_period_ms = 100;

    /* Enable temperature compensation and self-test */
    cfg->enable_temp_comp = 1;
    cfg->enable_self_test = 1;
    cfg->self_test_period_s = 3600;  /* Every hour */

    /* Initialize TEDS with empty defaults */
    memset(&cfg->teds, 0, sizeof(cfg->teds));
    cfg->teds.template_type = SS_TEDS_BASIC;
    strncpy(cfg->teds.manufacturer_name, "Unknown", sizeof(cfg->teds.manufacturer_name) - 1);
    strncpy(cfg->teds.model_number, "SS-0001", sizeof(cfg->teds.model_number) - 1);
}

/* =========================================================================
 * L2: Sensor Node Initialization
 * ========================================================================= */

/**
 * @brief Initialize a multi-channel sensor node
 *
 * A sensor node aggregates multiple sensor channels into a single
 * smart sensor system. Each channel can be a different physical
 * quantity. The channels are not yet configured — use
 * ss_node_add_channel() to configure each one.
 *
 * This represents a typical IoT/industrial sensor node architecture
 * as specified in IEEE 1451.0.
 *
 * @param node       Output node structure (must be allocated by caller)
 * @param n_channels Number of sensor channels (1..SS_MAX_CHANNELS)
 * @param node_id    Unique 64-bit node identifier
 * @return          0 on success, -1 on parameter error
 */
int ss_node_init(ss_node_t *node, uint8_t n_channels, uint64_t node_id) {
    if (!node || n_channels == 0 || n_channels > SS_MAX_CHANNELS) {
        return -1;
    }

    memset(node, 0, sizeof(*node));
    node->n_channels = n_channels;
    node->node_id = node_id;
    node->active_channel_mask = 0;
    node->timestamp_offset_s = 0.0;

    /* Initialize all channel states to idle */
    {
        uint8_t i;
        for (i = 0; i < n_channels; i++) {
            memset(&node->states[i], 0, sizeof(node->states[i]));
            node->states[i].mode = SS_MODE_IDLE;
            node->states[i].self_test_passed = -1; /* Not run */
            node->states[i].calibration_valid = 0;
        }
    }

    return 0;
}

/**
 * @brief Add (configure) a sensor channel in the node
 *
 * Assigns a sensor configuration to the next available channel
 * and marks it as active. Each channel must be configured before
 * measurements can be taken.
 *
 * @param node  Sensor node (must be initialized via ss_node_init)
 * @param cfg   Sensor configuration for this channel
 * @return      Channel index on success, -1 if node is full or NULL
 */
int ss_node_add_channel(ss_node_t *node, const ss_sensor_config_t *cfg) {
    uint8_t i;
    if (!node || !cfg) return -1;

    /* Find first unused channel slot */
    for (i = 0; i < node->n_channels; i++) {
        if (!(node->active_channel_mask & (1u << i))) {
            memcpy(&node->channels[i], cfg, sizeof(*cfg));
            node->active_channel_mask |= (1u << i);
            node->states[i].mode = SS_MODE_NORMAL_SAMPLING;
            return (int)i;
        }
    }

    return -1; /* Node is full */
}

/* =========================================================================
 * L2: Sensor State Management
 * ========================================================================= */

/**
 * @brief Read key fields from a sensor state
 *
 * Extract operating mode and current values from the sensor state
 * structure. This is a convenience accessor that avoids exposing
 * the internal structure layout to callers.
 *
 * @param state  Sensor state structure
 * @param mode   Output: current operating mode (can be NULL)
 * @param raw    Output: raw ADC reading (can be NULL)
 * @param eng    Output: engineering units value (can be NULL)
 */
void ss_sensor_state_read(const ss_sensor_state_t *state,
                          ss_operating_mode_t *mode,
                          double *raw, double *eng) {
    if (!state) return;
    if (mode) *mode = state->mode;
    if (raw)  *raw  = state->current_raw_value;
    if (eng)  *eng  = state->current_engineering_units;
}

/**
 * @brief Set sensor operating mode with state transition validation
 *
 * Validates and enforces legal mode transitions. Some transitions
 * are restricted for safety/reliability:
 *   - Cannot enter CALIBRATION from FAULT_DETECTED
 *   - Cannot enter NORMAL_SAMPLING from FAULT_DETECTED
 *     (must go through MAINTENANCE first)
 *
 * Allowed transitions (from -> to):
 *   IDLE            -> NORMAL_SAMPLING, SELF_TEST, CALIBRATION, SLEEP, MAINTENANCE
 *   NORMAL_SAMPLING -> IDLE, BURST_SAMPLING, SELF_TEST, SLEEP, FAULT_DETECTED
 *   SELF_TEST       -> IDLE, NORMAL_SAMPLING, FAULT_DETECTED
 *   CALIBRATION     -> IDLE, NORMAL_SAMPLING
 *   FAULT_DETECTED  -> MAINTENANCE
 *   MAINTENANCE     -> IDLE, SELF_TEST, CALIBRATION
 *   SLEEP           -> IDLE, NORMAL_SAMPLING
 *
 * @param state  Sensor state to modify
 * @param mode   Desired new operating mode
 * @return       0 on success, -1 if invalid transition
 */
int ss_sensor_state_set_mode(ss_sensor_state_t *state,
                              ss_operating_mode_t mode) {
    if (!state) return -1;
    if (state->mode == mode) return 0; /* No change needed */

    /* Validate transition */
    switch (state->mode) {
    case SS_MODE_IDLE:
        /* Can go to any operational mode */
        if (mode == SS_MODE_FAULT_DETECTED) return -1;
        break;

    case SS_MODE_NORMAL_SAMPLING:
        /* Cannot directly go to calibration or maintenance */
        if (mode == SS_MODE_CALIBRATION || mode == SS_MODE_MAINTENANCE)
            return -1;
        break;

    case SS_MODE_SELF_TEST:
    case SS_MODE_CALIBRATION:
    case SS_MODE_BURST_SAMPLING:
        /* These modes can transition to idle or normal or fault */
        if (mode != SS_MODE_IDLE && mode != SS_MODE_NORMAL_SAMPLING
            && mode != SS_MODE_FAULT_DETECTED)
            return -1;
        break;

    case SS_MODE_FAULT_DETECTED:
        /* Must go through maintenance to recover */
        if (mode != SS_MODE_MAINTENANCE) return -1;
        break;

    case SS_MODE_MAINTENANCE:
        /* After maintenance, can go to idle, self-test, or calibration */
        if (mode != SS_MODE_IDLE && mode != SS_MODE_SELF_TEST
            && mode != SS_MODE_CALIBRATION)
            return -1;
        break;

    case SS_MODE_SLEEP_LOW_POWER:
        /* Wake up to idle or normal sampling */
        if (mode != SS_MODE_IDLE && mode != SS_MODE_NORMAL_SAMPLING)
            return -1;
        break;

    default:
        return -1;
    }

    /* Handle fault entry/exit */
    if (mode == SS_MODE_FAULT_DETECTED) {
        state->fault_count++;
    }
    if (state->mode == SS_MODE_FAULT_DETECTED
        && mode == SS_MODE_MAINTENANCE) {
        state->recovery_count++;
    }

    state->mode = mode;
    return 0;
}

/**
 * @brief Update sensor state with new measurement data
 *
 * Records a new measurement into the sensor state, updating
 * the raw, processed, and engineering values along with the
 * timestamp and sample counter.
 *
 * Also updates the moving average and variance statistics.
 *
 * @param state  Sensor state to update
 * @param raw    Raw ADC reading (unprocessed)
 * @param proc   Conditioned value (after analog/digital processing)
 * @param eng    Final calibrated value in engineering units
 * @param ts_us  Timestamp in microseconds
 */
void ss_sensor_state_update(ss_sensor_state_t *state, double raw,
                            double proc, double eng, uint64_t ts_us) {
    if (!state) return;

    state->current_raw_value = raw;
    state->current_processed_value = proc;
    state->current_engineering_units = eng;
    state->last_sample_timestamp_us = ts_us;
    state->sample_count++;

    /* Update running statistics with the new engineering value */
    ss_sensor_state_update_stats(state, eng);
}

/**
 * @brief Update running statistics using Welford's online algorithm
 *
 * Welford's algorithm computes mean and variance in a single pass
 * without storing all samples. It is numerically stable — avoids
 * catastrophic cancellation that affects the naive two-pass formula.
 *
 * Welford recurrence:
 *   delta = x_n - mean_{n-1}
 *   mean_n = mean_{n-1} + delta / n
 *   M2_n = M2_{n-1} + delta * (x_n - mean_n)
 *   variance_n = M2_n / (n - 1)   (sample variance)
 *
 * Reference: Welford, B.P. (1962), Technometrics 4(3):419-420
 *
 * Complexity: O(1) time, O(1) space
 *
 * @param state   Sensor state with running mean and variance
 * @param new_val New measurement value
 */
void ss_sensor_state_update_stats(ss_sensor_state_t *state, double new_val) {
    double delta, delta2;
    if (!state) return;

    /* Use sample_count (already incremented by _update) as n */
    double n = (double)state->sample_count;
    if (n < 1.0) n = 1.0;

    delta = new_val - state->moving_average;
    state->moving_average += delta / n;

    delta2 = new_val - state->moving_average;
    /* Accumulate M2 = sum of squared deviations */
    /* We store variance directly; update using cumulative formula */
    if (state->sample_count == 1) {
        state->moving_variance = 0.0;
    } else {
        /* Use the running variance update (Welford online algorithm):
         * M2_n = M2_{n-1} + delta * delta2
         * S_n^2 = M2_n / (n-1)
         * = S_{n-1}^2 * (n-2)/(n-1) + delta * delta2 / (n-1) */
        if (n > 1.0) {
            state->moving_variance = ((n - 2.0) / (n - 1.0)) * state->moving_variance
                                    + delta * delta2 / (n - 1.0);
        }
    }
}

/* =========================================================================
 * L2: Measurement Buffer Operations
 * ========================================================================= */

/**
 * @brief Initialize a ring buffer for measurement storage
 *
 * Associates a pre-allocated measurement array with the buffer
 * metadata structure. The buffer operates as a circular FIFO:
 * new measurements overwrite the oldest ones when full.
 *
 * @param buf       Buffer metadata structure
 * @param buffer    Pre-allocated array of measurements
 * @param capacity  Maximum number of measurements
 */
void ss_measurement_buffer_init(ss_measurement_buffer_t *buf,
                                ss_measurement_t *buffer, size_t capacity) {
    if (!buf || !buffer || capacity == 0) return;

    memset(buf, 0, sizeof(*buf));
    buf->buffer = buffer;
    buf->capacity = capacity;
    buf->head = 0;
    buf->count = 0;
    buf->wrap_occurred = 0;

    /* Initialize all buffer entries to zero */
    memset(buffer, 0, capacity * sizeof(ss_measurement_t));
}

/**
 * @brief Push a measurement onto the ring buffer
 *
 * Deep-copies the measurement into the next available slot.
 * When the buffer is full, the oldest measurement is overwritten
 * (circular FIFO behavior).
 *
 * @param buf   Ring buffer
 * @param meas  Measurement to store
 * @return      Number of valid entries in the buffer after push
 */
size_t ss_measurement_buffer_push(ss_measurement_buffer_t *buf,
                                  const ss_measurement_t *meas) {
    if (!buf || !buf->buffer || !meas) return 0;

    /* Copy measurement into the buffer at the head position */
    memcpy(&buf->buffer[buf->head], meas, sizeof(ss_measurement_t));

    /* Advance head pointer */
    buf->head = (buf->head + 1) % buf->capacity;

    /* Update count */
    if (buf->count < buf->capacity) {
        buf->count++;
    } else {
        buf->wrap_occurred = 1;
    }

    return buf->count;
}

/**
 * @brief Retrieve a measurement by logical index
 *
 * Index 0 = oldest measurement still in the buffer.
 * Index count-1 = most recent measurement.
 *
 * Time complexity: O(1)
 *
 * @param buf    Ring buffer
 * @param index  Logical index (0 to count-1)
 * @param meas   Output: copy of the requested measurement
 * @return       0 on success, -1 if index out of range
 */
int ss_measurement_buffer_get(const ss_measurement_buffer_t *buf,
                              size_t index, ss_measurement_t *meas) {
    size_t phys_index;

    if (!buf || !meas || index >= buf->count) {
        return -1;
    }

    /* Compute physical index accounting for wrap-around.
     * Oldest entry is at (head - count) mod capacity. */
    if (buf->count < buf->capacity || !buf->wrap_occurred) {
        /* Buffer hasn't wrapped — oldest is at index 0 */
        phys_index = index;
    } else {
        /* Buffer has wrapped — circular addressing */
        phys_index = (buf->head + buf->capacity - buf->count + index)
                     % buf->capacity;
    }

    memcpy(meas, &buf->buffer[phys_index], sizeof(ss_measurement_t));
    return 0;
}

/**
 * @brief Compute arithmetic mean of calibrated values in buffer
 *
 * Uses the engineering-units (calibrated_value) of each measurement.
 * Returns NaN for an empty buffer.
 *
 * Complexity: O(n) where n = buffer count
 *
 * @param buf  Measurement buffer
 * @return     Mean value, or NaN if empty
 */
double ss_measurement_buffer_mean(const ss_measurement_buffer_t *buf) {
    size_t i;
    double sum = 0.0;
    ss_measurement_t meas;

    if (!buf || buf->count == 0) {
        return nan("");
    }

    for (i = 0; i < buf->count; i++) {
        if (ss_measurement_buffer_get(buf, i, &meas) == 0) {
            sum += meas.calibrated_value;
        }
    }

    return sum / (double)buf->count;
}

/**
 * @brief Compute sample standard deviation of buffer values
 *
 * Uses the two-pass algorithm for numerical stability:
 *   1. Compute mean
 *   2. Compute sum of squared deviations
 *   3. sigma = sqrt(sum_squared_dev / (n-1))
 *
 * Complexity: O(n)
 *
 * @param buf  Measurement buffer
 * @return     Sample standard deviation, or NaN if < 2 samples
 */
double ss_measurement_buffer_stddev(const ss_measurement_buffer_t *buf) {
    size_t i;
    double mean, sum_sq_dev = 0.0;
    ss_measurement_t meas;

    if (!buf || buf->count < 2) {
        return nan("");
    }

    mean = ss_measurement_buffer_mean(buf);

    for (i = 0; i < buf->count; i++) {
        if (ss_measurement_buffer_get(buf, i, &meas) == 0) {
            double dev = meas.calibrated_value - mean;
            sum_sq_dev += dev * dev;
        }
    }

    return sqrt(sum_sq_dev / (double)(buf->count - 1));
}

/**
 * @brief Find min and max calibrated values in buffer
 *
 * Linear scan through all entries. For an empty buffer,
 * returns -1 and does not modify min/max outputs.
 *
 * Complexity: O(n)
 *
 * @param buf  Measurement buffer
 * @param min  Output: minimum value
 * @param max  Output: maximum value
 * @return     0 on success, -1 if buffer empty
 */
int ss_measurement_buffer_range(const ss_measurement_buffer_t *buf,
                                double *min, double *max) {
    size_t i;
    double vmin, vmax;
    ss_measurement_t meas;

    if (!buf || buf->count == 0 || !min || !max) {
        return -1;
    }

    /* Initialize with first value */
    if (ss_measurement_buffer_get(buf, 0, &meas) != 0) {
        return -1;
    }
    vmin = vmax = meas.calibrated_value;

    for (i = 1; i < buf->count; i++) {
        if (ss_measurement_buffer_get(buf, i, &meas) == 0) {
            double v = meas.calibrated_value;
            if (v < vmin) vmin = v;
            if (v > vmax) vmax = v;
        }
    }

    *min = vmin;
    *max = vmax;
    return 0;
}

/**
 * @brief Clear all measurements from the buffer
 *
 * Resets the buffer to empty state without deallocating memory.
 * All metadata counters are zeroed.
 *
 * @param buf  Buffer to clear
 */
void ss_measurement_buffer_clear(ss_measurement_buffer_t *buf) {
    if (!buf) return;

    buf->head = 0;
    buf->count = 0;
    buf->wrap_occurred = 0;

    if (buf->buffer) {
        memset(buf->buffer, 0, buf->capacity * sizeof(ss_measurement_t));
    }
}
