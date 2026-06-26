/**
 * @file adc_dac_core.h
 * @brief ADC/DAC core definitions, architecture enums, and fundamental structures.
 *
 * Knowledge Coverage:
 *   L1 - Definitions: ADC resolution, sampling rate, reference voltage,
 *        quantization step (LSB), conversion time, aperture jitter, input range
 *   L2 - Core Concepts: ADC architecture taxonomy, linearity metrics,
 *        Nyquist vs oversampling, single-ended vs differential
 *   L4 - Fundamental Laws: Nyquist-Shannon sampling theorem (stated),
 *        ideal SNR formula (stated)
 *
 * References:
 *   - Maloberti, "Data Converters" (2007), Ch. 1-2
 *   - Razavi, "Principles of Data Conversion System Design" (1995)
 *   - IEEE Std 1241-2010: Terminology and Test Methods for ADC
 *   - Sedra & Smith, "Microelectronic Circuits" (2020), Ch. 17
 *
 * Course Alignment:
 *   MIT 6.301 Solid-State Circuits, Stanford EE315 Data Converters,
 *   Berkeley EE140 Analog IC, Michigan EECS 511 Analog Circuits
 */

#ifndef ADC_DAC_CORE_H
#define ADC_DAC_CORE_H

#include <stdint.h>
#include <stddef.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1 Definitions — Core Quantities
 * ========================================================================= */

/** ADC Architecture enumerates the major converter topologies taught in
 *  graduate-level data converter courses (Stanford EE315, MIT 6.301). */
typedef enum {
    ADC_ARCH_FLASH,              /**< Full-parallel, fastest, N-bit needs 2^N-1 comparators */
    ADC_ARCH_SAR,                /**< Successive Approximation Register, medium speed/power */
    ADC_ARCH_PIPELINE,           /**< Multi-stage with inter-stage gain, high-speed/resolution */
    ADC_ARCH_SIGMA_DELTA,        /**< Oversampling + noise shaping, highest resolution */
    ADC_ARCH_DUAL_SLOPE,         /**< Integrating, very slow but high accuracy at DC */
    ADC_ARCH_FOLDING_INTERPOLATING, /**< Folding reduces comparator count vs flash */
    ADC_ARCH_TIME_INTERLEAVED,   /**< M parallel ADCs staggered in time to multiply sample rate */
    ADC_ARCH_SUBRANGING,         /**< Two-step flash, coarse + fine */
    ADC_ARCH_COUNT               /**< Sentinel */
} adc_architecture_t;

/** DAC Architecture enum. */
typedef enum {
    DAC_ARCH_BINARY_WEIGHTED,    /**< Binary-weighted resistors/currents */
    DAC_ARCH_R2R_LADDER,         /**< R-2R ladder, only two resistor values needed */
    DAC_ARCH_STRING,             /**< Resistor string / Kelvin divider, inherently monotonic */
    DAC_ARCH_CURRENT_STEERING,   /**< Switched current sources, high-speed */
    DAC_ARCH_SIGMA_DELTA,        /**< Oversampling DAC with noise shaping */
    DAC_ARCH_CHARGE_SCALING,     /**< Switched-capacitor DAC, CMOS compatible */
    DAC_ARCH_HYBRID,             /**< Combined architectures for segmentation */
    DAC_ARCH_COUNT
} dac_architecture_t;

/** Endianness of binary output representation. */
typedef enum {
    ADC_CODING_OFFSET_BINARY,    /**< 0...0 = -FS, 1...1 = +FS - 1LSB */
    ADC_CODING_TWOS_COMPLEMENT,  /**< MSB inverted, most common in DSP */
    ADC_CODING_GRAY,             /**< Single-bit transitions between adjacent codes */
    ADC_CODING_THERMOMETER,      /**< Flash ADC unary output */
    ADC_CODING_COUNT
} adc_coding_t;

/**
 * @brief Fundamental ADC specification struct.
 *
 * Captures the key parameters from a datasheet that determine
 * the converter's ideal and non-ideal behavior.
 *
 * Quantization step (V_LSB) = V_ref / (2^N) for full-scale V_ref.
 *
 * Ideal SNR (Nyquist): SNR[dB] = 6.02*N + 1.76  (L4 — derived from q²/12)
 */
typedef struct {
    uint32_t    resolution_bits;    /**< N, number of output bits */
    uint32_t    sample_rate_hz;     /**< f_s, samples per second */
    double      v_ref_volts;        /**< Reference voltage, full-scale = V_ref */
    double      v_lsb_volts;        /**< Derived: V_ref / 2^N */
    double      input_range_min_v;  /**< Minimum input voltage */
    double      input_range_max_v;  /**< Maximum input voltage */
    double      aperture_jitter_ps; /**< RMS aperture jitter in picoseconds */
    double      conversion_time_ns; /**< Conversion time in nanoseconds */
    adc_architecture_t architecture; /**< Converter topology */
    adc_coding_t coding;            /**< Output coding scheme */
    uint8_t     single_ended;       /**< 1 = single-ended, 0 = differential */
} adc_spec_t;

/**
 * @brief Fundamental DAC specification struct.
 */
typedef struct {
    uint32_t    resolution_bits;    /**< N, number of input bits */
    uint32_t    update_rate_hz;     /**< Maximum update rate */
    double      v_ref_volts;        /**< Reference voltage */
    double      v_lsb_volts;        /**< Derived: V_ref / 2^N */
    double      output_range_min_v; /**< Minimum output voltage */
    double      output_range_max_v; /**< Maximum output voltage */
    double      settling_time_ns;   /**< Settling time to within 0.5 LSB */
    double      glitch_energy_pVs;  /**< Glitch impulse energy [pV·s] */
    dac_architecture_t architecture;
    uint8_t     differential;       /**< 1 = differential output */
    double      output_impedance_ohm; /**< Output impedance */
} dac_spec_t;

/* =========================================================================
 * L1 Definitions — Conversion Result and Error Types
 * ========================================================================= */

/**
 * @brief Single ADC conversion result with metadata.
 *
 * Maps the analog input to a digital code and provides the
 * reconstructed analog value for error analysis.
 */
typedef struct {
    uint64_t    raw_code;           /**< Uncorrected output code */
    uint64_t    corrected_code;     /**< After gain/offset calibration */
    double      input_voltage;      /**< Actual analog input applied */
    double      reconstructed_volts;/**< code * V_LSB */
    double      quantization_error; /**< reconstructed - input */
    int64_t     timestamp_ns;       /**< Sample timestamp */
} adc_sample_t;

/**
 * @brief DAC output sample.
 */
typedef struct {
    uint64_t    input_code;         /**< Digital code applied */
    double      ideal_output;       /**< code * V_LSB, ideal output */
    double      actual_output;      /**< Measured output voltage */
    double      deviation;          /**< actual - ideal */
    double      settling_error;     /**< Residual error after settling time */
} dac_sample_t;

/* =========================================================================
 * L2 Core Concept — Non-Linearity Definitions
 * ========================================================================= */

/**
 * @brief DNL/INL result for one ADC/DAC code transition.
 *
 * DNL[k] = (V_transition[k+1] - V_transition[k]) / V_LSB - 1
 * INL[k] = (V_transition[k] - V_ideal[k]) / V_LSB
 *
 * DNL < -1 implies missing codes (non-monotonic for DAC).
 */
typedef struct {
    uint32_t    code;               /**< Digital code index k */
    double      dnl_lsb;            /**< Differential Non-Linearity [LSB] */
    double      inl_lsb;            /**< Integral Non-Linearity [LSB] */
    double      transition_voltage; /**< Measured code transition voltage */
    double      ideal_transition_v; /**< k * V_LSB (offset corrected) */
} linearity_point_t;

/**
 * @brief Full linearity characterization sweep.
 */
typedef struct {
    uint32_t        num_codes;      /**< Number of codes measured */
    linearity_point_t *points;      /**< Array of num_codes points */
    double           max_dnl;       /**< Peak |DNL| */
    double           max_inl;       /**< Peak |INL| */
    double           avg_dnl;       /**< Mean |DNL| */
    double           avg_inl;       /**< Mean |INL| */
    uint32_t         missing_codes; /**< Count of codes with DNL < -0.9 LSB */
} linearity_sweep_t;

/* =========================================================================
 * L3 Mathematical Structures — Transfer Characteristic
 * ========================================================================= */

/**
 * @brief Full transfer function representation of an ADC.
 *
 * V_out[k] = gain * V_in + offset + sum_{i=1}^{M} a_i * h_i(V_in)
 * where h_i are harmonic distortion terms.
 *
 * This captures: offset, gain error, and non-linear distortion
 * in a unified polynomial model.
 */
typedef struct {
    double          gain;           /**< Actual gain (ideal = 1.0) */
    double          offset_volts;   /**< DC offset voltage */
    uint32_t        harmonic_order; /**< Max harmonic modeled */
    double         *harmonic_coeffs;/**< Coefficients a_1..a_M, length = harmonic_order */
    double          thd_percent;    /**< Total Harmonic Distortion [%] */
} transfer_char_t;

/* =========================================================================
 * L2 — Sample-and-Hold (Track-and-Hold)
 * ========================================================================= */

/**
 * @brief Sample-and-hold amplifier model.
 *
 * Finite acquisition bandwidth and aperture uncertainty
 * are the dominant dynamic error sources in sampling ADCs.
 */
typedef struct {
    double          bandwidth_hz;   /**< -3 dB bandwidth of sampling switch + cap */
    double          aperture_rms_ps;/**< RMS aperture jitter (clock + sampling jitter) */
    double          droop_rate_vps; /**< Hold capacitor droop rate [V/s] */
    double          hold_cap_pF;    /**< Hold capacitance */
    double          acquisition_time_ns; /**< Time to settle to 0.1% */
    double          charge_injection_fC; /**< Switch charge injection [fC] */
} sample_hold_spec_t;

/* =========================================================================
 * L4 Fundamental Law — Ideal ADC SNR
 * ========================================================================= */

/**
 * @brief Compute ideal SNR for an N-bit Nyquist-rate ADC.
 *
 * Derivation:
 *   Quantization noise power = V_LSB² / 12  (uniform distribution [-q/2, q/2])
 *   Signal power (full-scale sine) = (2^{N-1} * V_LSB)² / 2
 *   SNR = 10*log10(Signal/Noise) = 6.02*N + 1.76 [dB]
 *
 * This is the **Shannon-theoretic limit** for uniform quantization
 * of a full-scale sinusoidal input.
 *
 * @param n_bits ADC resolution in bits.
 * @return Ideal SNR in dB.
 * 
 * Complexity: O(1)
 * Reference: IEEE Std 1241-2010, Eq. (36)
 */
double adc_ideal_snr_db(uint32_t n_bits);

/**
 * @brief Compute SNR from signal and noise powers.
 *
 * SNR = 10 * log10(P_signal / P_noise)
 *
 * @param p_signal Signal power (mean squared).
 * @param p_noise  Noise power (mean squared).
 * @return SNR in dB, NAN if inputs invalid.
 *
 * Complexity: O(1)
 */
double adc_snr_from_powers(double p_signal, double p_noise);

/* =========================================================================
 * L2 — Effective Number of Bits (ENOB)
 * ========================================================================= */

/**
 * @brief Compute ENOB from measured SINAD.
 *
 * ENOB = (SINAD_dB - 1.76) / 6.02
 *
 * ENOB accounts for all non-idealities: quantization noise,
 * thermal noise, jitter, and distortion.
 *
 * @param sinad_db Measured SINAD in dB.
 * @return ENOB in bits.
 *
 * Complexity: O(1)
 * Reference: Kester, "Data Conversion Handbook" (2005), Ch. 2
 */
double adc_enob_from_sinad(double sinad_db);

/**
 * @brief Compute ENOB from SNDR (redundant with SINAD, same formula).
 */
double adc_enob_from_sndr(double sndr_db);

/* =========================================================================
 * L4 — Jitter-Limited SNR
 * ========================================================================= */

/**
 * @brief Compute SNR limitation due to aperture jitter.
 *
 * For a sinusoidal input at frequency f_in:
 *   SNR_jitter = -20 * log10(2 * PI * f_in * t_jitter_rms)
 *
 * This is the fundamental limit beyond which increased resolution
 * yields no improvement (jitter noise floor).
 *
 * @param freq_hz       Input signal frequency [Hz].
 * @param jitter_rms_ps RMS aperture jitter [ps].
 * @return Jitter-limited SNR in dB, negative for invalid inputs.
 *
 * Complexity: O(1)
 * Reference: Maloberti (2007) §2.5, IEEE Std 1241-2010 Annex B
 */
double adc_jitter_snr_limit(double freq_hz, double jitter_rms_ps);

/* =========================================================================
 * L2 — Input-Referred Noise
 * ========================================================================= */

/**
 * @brief Compute total input-referred noise from all sources.
 *
 * v_noise_total = sqrt(v_thermal² + v_quant² + v_jitter²)
 *
 * Quantization noise: v_quant_rms = V_LSB / sqrt(12)
 * Thermal noise: k*T/C on sampling capacitor
 *
 * @param v_lsb         LSB voltage.
 * @param cap_pF        Sampling capacitance [pF].
 * @param temp_celsius  Temperature [°C].
 * @return Total input-referred RMS noise in volts.
 */
double adc_total_input_noise(double v_lsb, double cap_pF, double temp_celsius);

/* =========================================================================
 * L5 — Offset and Gain Calibration
 * ========================================================================= */

/**
 * @brief Two-point calibration: corrects offset and gain errors.
 *
 * Uses min/max code measurements to solve:
 *   corrected = (raw - offset_code) * gain_correction_factor
 *
 * @param raw_code           Uncalibrated output code.
 * @param code_at_min_input  Code with input at -FS (or min).
 * @param code_at_max_input  Code with input at +FS (or max).
 * @param n_bits             ADC resolution.
 * @return Calibrated code.
 */
uint64_t adc_two_point_calibrate(uint64_t raw_code,
                                  uint64_t code_at_min_input,
                                  uint64_t code_at_max_input,
                                  uint32_t n_bits);

/**
 * @brief Compute calibration parameters (gain factor and offset).
 *
 * @param[in]  code_at_min  Code at minimum input.
 * @param[in]  code_at_max  Code at maximum input.
 * @param[in]  n_bits       Resolution.
 * @param[out] offset_code  Computed offset in LSB.
 * @param[out] gain_factor  Gain correction factor (ideal = 1.0).
 */
void adc_calibration_params(uint64_t code_at_min,
                            uint64_t code_at_max,
                            uint32_t n_bits,
                            double *offset_code,
                            double *gain_factor);

/* =========================================================================
 * L6 — Classic Conversion Functions
 * ========================================================================= */

/**
 * @brief ADC code to voltage reconstruction.
 *
 * V = code * V_ref / 2^N
 *
 * @param code    Raw digital output.
 * @param v_ref   Reference voltage.
 * @param n_bits  ADC resolution.
 * @return Reconstructed analog voltage.
 */
double adc_code_to_voltage(uint64_t code, double v_ref, uint32_t n_bits);

/**
 * @brief Voltage to ideal ADC code (quantization mid-tread).
 *
 * Code = round(V_in / V_LSB), clamped to [0, 2^N - 1].
 *
 * @param v_in    Input voltage (must be within [0, V_ref]).
 * @param v_ref   Reference voltage.
 * @param n_bits  ADC resolution.
 * @return Quantized code.
 */
uint64_t adc_voltage_to_code(double v_in, double v_ref, uint32_t n_bits);

/* =========================================================================
 * L5 — DNL / INL Computation
 * ========================================================================= */

/**
 * @brief Compute DNL from a set of measured transition voltages.
 *
 * DNL[k] = (T[k+1] - T[k]) / V_LSB_avg - 1
 *
 * @param transitions  Array of measured transition voltages (length: 2^N + 1).
 * @param n_bits       ADC resolution.
 * @param dnl_out      Output DNL array (length: 2^N - 1).
 *
 * Complexity: O(2^N)
 * Reference: IEEE Std 1241-2010 §4.3.3
 */
void adc_compute_dnl(const double *transitions, uint32_t n_bits,
                     double *dnl_out);

/**
 * @brief Compute INL from DNL by cumulative sum.
 *
 * INL[k] = sum_{i=0}^{k-1} DNL[i],  INL[0] = 0.
 *
 * @param dnl      DNL array (length: 2^N - 1).
 * @param n_bits   ADC resolution.
 * @param inl_out  Output INL array (length: 2^N).
 *
 * Complexity: O(2^N)
 */
void adc_compute_inl_from_dnl(const double *dnl, uint32_t n_bits,
                               double *inl_out);

/**
 * @brief Compute THD from polynomial transfer characteristic model.
 *
 * Evaluates harmonic distortion through Chebyshev polynomial expansion
 * of the Taylor series representing the non-linear transfer function.
 *
 * @param tc         Transfer characteristic (coefficients + gain).
 * @param amplitude  Input sinusoid amplitude.
 * @return THD in percent.
 *
 * Complexity: O(M) where M is harmonic_order
 */
double transfer_char_thd(const transfer_char_t *tc, double amplitude);

/* =========================================================================
 * Utility
 * ========================================================================= */

/** Return string name for ADC architecture. */
const char *adc_arch_name(adc_architecture_t arch);

/** Return string name for DAC architecture. */
const char *dac_arch_name(dac_architecture_t arch);

/** Return string name for coding scheme. */
const char *adc_coding_name(adc_coding_t coding);

#ifdef __cplusplus
}
#endif

#endif /* ADC_DAC_CORE_H */
