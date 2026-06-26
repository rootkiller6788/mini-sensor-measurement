/**
 * @file sar_adc.h
 * @brief Successive Approximation Register (SAR) ADC models,
 *        binary search algorithm, charge-redistribution DAC,
 *        calibration and redundancy.
 *
 * Knowledge Coverage:
 *   L1 - Definitions: SAR, comparator, capacitive DAC (CDAC),
 *        sampling phase, conversion phase, bit-cycling
 *   L2 - Core Concepts: binary search successive approximation,
 *        charge redistribution, split capacitor array,
 *        monotonic switching, Vcm-based switching
 *   L3 - Mathematical Structures: Binary-weighted capacitor
 *        charge equations, comparator metastability model,
 *        transient settling model
 *   L4 - Fundamental Laws: SAR resolution limit from kT/C noise,
 *        DNL/INL bounds from capacitor mismatch
 *   L5 - Algorithms: SAR binary search algorithm, redundancy/
 *        sub-radix-2 search, digital calibration from
 *        capacitor mismatch measurement
 *   L6 - Canonical Problems: 12-bit SAR ADC design with
 *        capacitor mismatch calibration
 *
 * References:
 *   - McCreary & Gray, "All-MOS Charge Redistribution ADC"
 *     IEEE JSSC (1975)
 *   - Kuttner, "A 1.2-V 10-b 20-Msample/s Nonbinary SAR ADC"
 *     IEEE JSSC (2002)
 *   - Harpe, "Successive Approximation ADCs" (2018), Ch. 1-5
 *   - Baker, "CMOS Circuit Design, Layout, and Simulation" (2019)
 *
 * Course Alignment:
 *   MIT 6.775 CMOS Data Converters
 *   Stanford EE315B Advanced Data Converters
 *   Michigan EECS 511 Analog Circuits
 *   Berkeley EE247 Advanced Analog IC Design
 */

#ifndef SAR_ADC_H
#define SAR_ADC_H

#include <stdint.h>
#include "adc_dac_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1 — SAR ADC Specifications
 * ========================================================================= */

/** SAR ADC switching scheme enumeration. */
typedef enum {
    SAR_SWITCHING_TRADITIONAL,   /**< Traditional: V_ref, 0 switching */
    SAR_SWITCHING_MONOTONIC,     /**< Monotonic: only discharge */
    SAR_SWITCHING_VCM_BASED,     /**< V_cm-based switching: lower power */
    SAR_SWITCHING_SPLIT,         /**< Split capacitor for attenuation */
    SAR_SWITCHING_COUNT
} sar_switching_t;

/**
 * @brief SAR ADC specification including CDAC parameters.
 *
 * CDAC total capacitance C_total = 2^N * C_unit.
 * kT/C noise: v_n_rms = sqrt(kT / C_total).
 */
typedef struct {
    uint32_t        resolution_bits;    /**< N, number of bits */
    uint32_t        sample_rate_hz;     /**< Sampling frequency */
    double          v_ref_volts;        /**< Reference voltage */
    double          c_unit_fF;         /**< Unit capacitor [fF] */
    double          c_total_pF;        /**< Total DAC capacitance [pF] */
    sar_switching_t switching;          /**< Switching scheme */
    double          comparator_offset_mV; /**< Comparator offset [mV] */
    double          comparator_noise_uVrms; /**< Comparator RMS noise */
    double          comparator_delay_ns; /**< Comparator delay */
    double          logic_delay_ns;     /**< SAR logic delay per bit */
    uint32_t        redundancy_bits;    /**< Number of redundant bits (sub-radix-2) */
    double          radix;              /**< Actual radix (< 2.0 for redundancy) */
} sar_adc_spec_t;

/**
 * @brief SAR ADC runtime state (per conversion).
 */
typedef struct {
    uint32_t        current_bit;        /**< Bit index being tested (MSB=resolution-1) */
    uint64_t        dac_code;           /**< Current DAC code */
    double          dac_voltage;        /**< Current DAC output voltage */
    double          comparator_input;   /**< Vin - Vdac */
    int             comparator_result;  /**< 1, -1, or 0 (meta-stable) */
    uint32_t        cycles_remaining;   /**< Remaining bit trials */
    double          cdac_node_voltage;  /**< Voltage at top plate of CDAC */
    double          sample_voltage;     /**< Held input sample */
} sar_adc_state_t;

/* =========================================================================
 * L2 — CDAC Charge Model
 * ========================================================================= */

/**
 * @brief CDAC capacitor array model.
 *
 * Binary-weighted array: C_k = 2^k * C_unit for k = 0..N-1.
 * Terminal capacitor C_term = C_unit.
 */
typedef struct {
    uint32_t        num_bits;           /**< N */
    double          c_unit_fF;          /**< Unit capacitance [fF] */
    double         *capacitors_fF;      /**< Array: C[0]...C[N-1], C_term */
    double          c_total_fF;         /**< Sum of all capacitors */
    uint64_t        switch_state;       /**< Bitmask: 1=V_ref, 0=GND */
    double          top_plate_voltage;  /**< Voltage at CDAC top plate */
} cdac_array_t;

/**
 * @brief Initialize a binary-weighted CDAC array.
 *
 * @param num_bits  Resolution.
 * @param c_unit    Unit capacitance [fF].
 * @param[out] cdac Initialized CDAC struct (caller must call cdac_free).
 * @return 0 on success.
 */
int cdac_init(uint32_t num_bits, double c_unit, cdac_array_t *cdac);

/** Free CDAC memory. */
void cdac_free(cdac_array_t *cdac);

/**
 * @brief Compute CDAC output voltage for a given switch state.
 *
 * V_out = V_ref * sum(C_i * b_i) / C_total
 * where b_i = 1 if switch_i = V_ref, 0 if GND.
 *
 * @param cdac  CDAC array.
 * @param code  Switch state (bitmask).
 * @return Output voltage.
 *
 * Complexity: O(N)
 * Reference: McCreary & Gray (1975)
 */
double cdac_output_voltage(const cdac_array_t *cdac, uint64_t code);

/**
 * @brief Sample input voltage onto CDAC top plate.
 *
 * During sampling phase: top plate = V_in, bottom plates = V_cm.
 * During hold-to-conversion transition: bottom plates switch,
 * causing charge redistribution that produces V_out.
 *
 * @param cdac    CDAC array.
 * @param v_in    Input voltage to sample.
 * @param v_cm    Common-mode voltage for bottom plates during sampling.
 */
void cdac_sample(cdac_array_t *cdac, double v_in, double v_cm);

/**
 * @brief Set the switch state (code) for the CDAC.
 *
 * @param cdac  CDAC array.
 * @param code  Bitmask of switches to V_ref (LSB = C_unit).
 */
void cdac_set_code(cdac_array_t *cdac, uint64_t code);

/* =========================================================================
 * L5 — SAR Binary Search Algorithm
 * ========================================================================= */

/**
 * @brief Perform one bit-trial in SAR conversion.
 *
 * Standard binary search:
 *   1. Set test bit to 1
 *   2. If DAC voltage > V_in: clear the bit (V_in < DAC)
 *      Else: keep the bit (V_in ≥ DAC)
 *   3. Move to next bit (dividing step by 2)
 *
 * @param state           SAR state.
 * @param cdac            CDAC array.
 * @param v_in            Sampled input voltage.
 * @param comparator_offset Comparator offset (to be subtracted).
 *
 * Complexity: O(1) per bit trial
 */
void sar_bit_trial(sar_adc_state_t *state, cdac_array_t *cdac,
                    double v_in, double comparator_offset);

/**
 * @brief Run a complete SAR conversion (all bit trials).
 *
 * @param state      SAR state (initialized externally).
 * @param cdac       CDAC array.
 * @param v_in       Sampled input voltage.
 * @param spec       SAR ADC specification.
 * @param[out] code  Resulting digital code.
 *
 * Complexity: O(N^2) naive, O(N) with optimized CDAC
 */
void sar_full_conversion(sar_adc_state_t *state, cdac_array_t *cdac,
                          double v_in, const sar_adc_spec_t *spec,
                          uint64_t *code);

/* =========================================================================
 * L5 — Redundant SAR (Sub-Radix-2)
 * ========================================================================= */

/**
 * @brief Compute sub-radix-2 weights for redundant SAR.
 *
 * For radix r < 2, weights are w_k = r^k * w_unit.
 * Redundancy allows recovery from incorrect decisions.
 *
 * Total range = sum(r^k) for k=0..N-1, must cover 2^N LSB.
 *
 * @param resolution_bits  Target resolution.
 * @param radix            Actual radix (< 2.0).
 * @param[out] weights     Array of length resolution_bits + redundancy_bits.
 * @param[out] total_steps Number of steps (N + redundancy).
 *
 * Complexity: O(N)
 */
void sar_subradix_weights(uint32_t resolution_bits, double radix,
                           double *weights, uint32_t *total_steps);

/**
 * @brief Redundant SAR conversion with error correction.
 *
 * In sub-radix-2, the search range overlaps, so a "wrong" bit
 * decision can be compensated by later bits.
 *
 * @param state      SAR state.
 * @param cdac       CDAC array.
 * @param v_in       Input voltage.
 * @param weights    Sub-radix-2 weights.
 * @param n_steps    Number of steps.
 * @param[out] code  Output code.
 */
void sar_redundant_conversion(sar_adc_state_t *state, cdac_array_t *cdac,
                               double v_in,
                               const double *weights, uint32_t n_steps,
                               uint64_t *code);

/* =========================================================================
 * L4 — kT/C Noise Limit
 * ========================================================================= */

/**
 * @brief Compute kT/C noise RMS voltage.
 *
 * v_n_rms = sqrt(k * T / C)
 * where k = 1.380649e-23 J/K, T in Kelvin, C in Farads.
 *
 * For a 12-bit ADC with V_ref=3.3V, V_LSB=806μV, to keep noise < 0.5 LSB:
 * C > kT / (0.5*V_LSB)²
 *
 * @param temp_kelvin   Temperature in Kelvin.
 * @param capacitance_F Capacitance in Farads.
 * @return RMS noise voltage.
 *
 * Complexity: O(1)
 */
double ktc_noise_rms(double temp_kelvin, double capacitance_F);

/**
 * @brief Compute minimum capacitance for a given noise budget.
 *
 * C_min = kT / (v_noise_rms)²
 *
 * @param temp_kelvin       Temperature.
 * @param noise_target_vrms Target RMS noise.
 * @return Minimum capacitance [F].
 */
double ktc_min_capacitance(double temp_kelvin, double noise_target_vrms);

/* =========================================================================
 * L5 — Comparator Model
 * ========================================================================= */

/**
 * @brief Model comparator decision with offset and noise.
 *
 * Models the comparator as:
 *   decision = sign(V_diff + V_offset + V_noise)
 * where V_noise ~ N(0, σ²) and V_offset is systematic offset.
 *
 * @param v_diff    Differential input voltage.
 * @param offset    Systematic offset.
 * @param noise_rms RMS noise.
 * @param noise_val Random value from N(0,1) distribution.
 * @return 1 if (v_diff + offset + noise_rms*noise_val) > 0, else -1.
 */
int comparator_decision(double v_diff, double offset,
                         double noise_rms, double noise_val);

/**
 * @brief Model comparator metastability probability.
 *
 * For differential latch comparator with time constant τ:
 *   P(meta) = exp(-t_resolution / τ)
 *
 * @param decision_time_ns  Time allocated for decision.
 * @param tau_ns            Regeneration time constant.
 * @param v_diff            Input differential voltage.
 * @return Probability of metastable decision.
 */
double comparator_metastability_prob(double decision_time_ns,
                                      double tau_ns, double v_diff);

/* =========================================================================
 * L5 — Capacitor Mismatch Model
 * ========================================================================= */

/**
 * @brief Model capacitor mismatch with Pelgrom's law.
 *
 * σ(ΔC/C) = A_C / sqrt(W * L)
 * where A_C is the Pelgrom mismatch coefficient.
 *
 * @param c_nominal     Nominal capacitance [fF].
 * @param area_um2      Capacitor area [μm²].
 * @param pelgrom_coeff Pelgrom coefficient A_C [%·μm].
 * @param random_val    Random value from N(0,1).
 * @return Actual capacitance with mismatch [fF].
 */
double capacitor_with_mismatch(double c_nominal, double area_um2,
                                double pelgrom_coeff, double random_val);

/**
 * @brief Compute DNL from capacitor mismatch for binary-weighted CDAC.
 *
 * For binary-weighted array with independent mismatches σ_i:
 * DNL_max ≈ sqrt(sum(σ_i²)) * sqrt(2^{N-1}) in LSB.
 *
 * @param n_bits      Resolution.
 * @param mismatch_pct Capacitor mismatch standard deviation [%].
 * @return Predicted DNL [LSB].
 */
double sar_dnl_from_mismatch(uint32_t n_bits, double mismatch_pct);

/* =========================================================================
 * L6 — SAR ADC Complete Model
 * ========================================================================= */

/**
 * @brief Simulate a complete SAR ADC including non-idealities.
 *
 * Models: kT/C noise, comparator offset, comparator noise,
 * capacitor mismatch, sampling jitter.
 *
 * @param spec            SAR ADC specification.
 * @param v_in_array      Array of input voltages.
 * @param n_conversions   Number of conversions.
 * @param[out] codes      Output codes.
 * @param[out] noise_seed Seed for noise generator (updated).
 *
 * Complexity: O(n_conversions * N) where N is resolution
 */
void sar_simulate_conversions(const sar_adc_spec_t *spec,
                               const double *v_in_array,
                               uint32_t n_conversions,
                               uint64_t *codes,
                               uint64_t *noise_seed);

/* =========================================================================
 * L5 — SAR Calibration
 * ========================================================================= */

/**
 * @brief Compute calibration weights from measured bit weights.
 *
 * Measures actual weight of each bit and stores correction factors.
 * The calibrated code is computed as:
 *   code_cal = sum(b_i * w_measured_i) / w_ideal
 *
 * @param measured_weights Array of measured bit weights (length: N).
 * @param n_bits           Resolution.
 * @param[out] corrections Array of correction factors (length: N).
 */
void sar_calibrate_weights(const double *measured_weights,
                            uint32_t n_bits, double *corrections);

/**
 * @brief Apply calibration corrections to a raw SAR code.
 *
 * @param raw_code     Uncalibrated code.
 * @param corrections  Calibration correction factors.
 * @param n_bits       Resolution.
 * @return Calibrated code.
 */
uint64_t sar_apply_calibration(uint64_t raw_code,
                                const double *corrections,
                                uint32_t n_bits);

#ifdef __cplusplus
}
#endif

#endif /* SAR_ADC_H */
