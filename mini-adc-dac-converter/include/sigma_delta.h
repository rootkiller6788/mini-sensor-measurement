/**
 * @file sigma_delta.h
 * @brief Sigma-Delta (ΔΣ) modulator theory and simulation.
 *        Noise shaping, NTFs, stability, decimation filtering.
 *
 * Knowledge Coverage:
 *   L1 - Definitions: oversampling ratio (OSR), noise transfer function
 *        (NTF), signal transfer function (STF), in-band noise,
 *        idle tones, modulator order, quantizer levels
 *   L2 - Core Concepts: noise shaping principle, first-order and
 *        second-order ΔΣ modulation, multi-bit vs single-bit quantizers,
 *        dithering in ΔΣ, cascaded (MASH) architectures
 *   L3 - Mathematical Structures: linearized z-domain model,
 *        STF = z^{-1} (or unity in baseband), NTF = (1 - z^{-1})^L,
 *        in-band noise power: σ²_e_inband = σ²_e * ∫|NTF(f)|² df
 *   L4 - Fundamental Laws: in-band noise for L-th order ΔΣ:
 *        IBN ≈ σ²_e * π^{2L} / ((2L+1) * OSR^{2L+1})
 *        leading to DR = (3/2) * (2^B - 1)² * (2L+1)*OSR^{2L+1} / π^{2L}
 *   L5 - Algorithms: delta-sigma modulator simulation (time-domain),
 *        decimation filter design (CIC + compensation FIR),
 *        MASH 1-1-1, MASH 2-1 stability analysis
 *   L6 - Canonical Problems: design of a 2nd-order audio ΔΣ ADC,
 *        decimation chain for GSM baseband
 *
 * References:
 *   - Schreier & Temes, "Understanding Delta-Sigma Data Converters" (2005)
 *   - Norsworthy, Schreier, Temes, "Delta-Sigma Data Converters" (1997)
 *   - Candy & Temes, "Oversampling Delta-Sigma Data Converters" (1992)
 *
 * Course Alignment:
 *   Stanford EE315 Data Converters (ΔΣ theory and design)
 *   MIT 6.775 CMOS Data Converters
 *   Berkeley EE247 Advanced Analog IC (ΔΣ modulators)
 */

#ifndef SIGMA_DELTA_H
#define SIGMA_DELTA_H

#include <stdint.h>
#include "adc_dac_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1 — ΔΣ Modulator Specifications
 * ========================================================================= */

/** ΔΣ modulator architecture types. */
typedef enum {
    SDM_SINGLE_LOOP_1ST,        /**< First-order, single loop */
    SDM_SINGLE_LOOP_2ND,        /**< Second-order, single loop */
    SDM_SINGLE_LOOP_HIGHER,     /**< Higher-order single loop (with local feedback) */
    SDM_MASH_1_1,               /**< MASH 1-1 cascaded (2nd order equivalent) */
    SDM_MASH_2_1,               /**< MASH 2-1 cascaded (3rd order equivalent) */
    SDM_MASH_2_2,               /**< MASH 2-2 cascaded (4th order equivalent) */
    SDM_BANDPASS,               /**< Band-pass ΔΣ for IF sampling */
    SDM_COUNT
} sdm_architecture_t;

/**
 * @brief ΔΣ modulator design specification.
 */
typedef struct {
    uint32_t            order;          /**< L, noise shaping order */
    uint32_t            quantizer_bits; /**< B, internal quantizer resolution (1 = single-bit) */
    double              osr;            /**< Oversampling ratio */
    double              fs_hz;          /**< Sampling frequency (oversampled) */
    double              bandwidth_hz;   /**< Signal bandwidth (Nyquist baseband) */
    sdm_architecture_t  arch;           /**< Architecture type */
    double             *ntf_coeffs_num; /**< NTF numerator coefficients [order+1] */
    double             *ntf_coeffs_den; /**< NTF denominator coefficients [order+1] */
    double             *stf_coeffs_num; /**< STF numerator coefficients [order+1] */
    double             *stf_coeffs_den; /**< STF denominator coefficients [order+1] */
    double              v_ref;          /**< Reference voltage */
} sdm_spec_t;

/**
 * @brief ΔΣ modulator runtime state.
 */
typedef struct {
    sdm_spec_t          spec;           /**< Design specification */
    /* Integrator states: state[0..order-1] */
    double             *integrator_state;
    uint32_t            order;
    /* Error memory for noise shaping */
    double              prev_error;
    double              prev_error_2;
    /* Statistics */
    uint64_t            sample_count;
    double              accumulated_power;
} sdm_state_t;

/* =========================================================================
 * L4 — ΔΣ Performance Prediction
 * ========================================================================= */

/**
 * @brief Compute theoretical in-band noise power for L-th order ΔΣ.
 *
 * IBN = Δ²/12 * π^(2L) / ((2L+1) * OSR^(2L+1))
 * where Δ = V_ref / (2^B - 1) is the quantizer step size.
 *
 * This formula assumes: linear model holds (quantizer replaced by
 * additive white noise source), ideal NTF = (1 - z^{-1})^L,
 * and that the quantizer does not overload.
 *
 * @param order      L (noise shaping order).
 * @param osr        Oversampling ratio.
 * @param n_bits     Internal quantizer bits.
 * @param v_ref      Reference voltage.
 * @return In-band noise power (variance) in V².
 *
 * Complexity: O(1)
 * Reference: Schreier & Temes (2005), §2.4.2
 */
double sdm_inband_noise_power(uint32_t order, double osr,
                               uint32_t n_bits, double v_ref);

/**
 * @brief Compute ideal dynamic range of L-th order ΔΣ ADC.
 *
 * DR[dB] = 10*log10( P_signal_max / IBN )
 * with full-scale sine peak power P_signal_max = (V_ref/2)²/2.
 *
 * @param order      Noise shaping order.
 * @param osr        Oversampling ratio.
 * @param n_bits     Quantizer bits.
 * @param v_ref      Reference voltage.
 * @return Dynamic range in dB.
 */
double sdm_dynamic_range_db(uint32_t order, double osr,
                             uint32_t n_bits, double v_ref);

/**
 * @brief Compute equivalent resolution (ENOB) from ΔΣ parameters.
 *
 * ENOB = (DR_dB - 1.76) / 6.02
 */
double sdm_equivalent_enob(uint32_t order, double osr,
                            uint32_t n_bits, double v_ref);

/* =========================================================================
 * L5 — NTF Design
 * ========================================================================= */

/**
 * @brief Design an NTF with Butterworth high-pass characteristic.
 *
 * Creates NTF(z) coefficients for a high-pass noise shaping filter
 * with cutoff at signal bandwidth.
 *
 * @param order           Filter order L.
 * @param bandwidth_hz    Signal bandwidth.
 * @param fs_hz           Sampling frequency.
 * @param[out] ntf_num    Numerator coefficients [order+1].
 * @param[out] ntf_den    Denominator coefficients [order+1].
 *
 * Complexity: O(order²) for Butterworth design
 */
void sdm_design_ntf_butterworth(uint32_t order, double bandwidth_hz,
                                 double fs_hz,
                                 double *ntf_num, double *ntf_den);

/**
 * @brief Design pure differentiation NTF: NTF(z) = (1 - z^{-1})^L.
 *
 * @param order     L.
 * @param[out] ntf_num Numerator [order+1], e.g., [1, -1] for L=1.
 * @param[out] ntf_den Denominator [order+1], all 1.0 for FIR NTF.
 */
void sdm_design_ntf_pure_diff(uint32_t order, double *ntf_num, double *ntf_den);

/* =========================================================================
 * L5 — ΔΣ Modulator Time-Domain Simulation
 * ========================================================================= */

/**
 * @brief Initialize a ΔΣ modulator state.
 *
 * Allocates integrator states and initializes to zero.
 *
 * @param spec   Modulator specification.
 * @param[out] state Initialized state (caller must call sdm_free_state).
 * @return 0 on success, -1 on error.
 */
int sdm_init_state(const sdm_spec_t *spec, sdm_state_t *state);

/** Free state allocated by sdm_init_state. */
void sdm_free_state(sdm_state_t *state);

/**
 * @brief Process one sample through a first-order ΔΣ modulator.
 *
 * Implements:
 *   u = input_sample
 *   v = integrator + u
 *   y = quantize(v)    // single-bit: sign(v)
 *   e = y - v          // quantization error
 *   integrator = integrator + u - y   // discrete-time integrator
 *
 * NTF from input to output: (1 - z^{-1})
 * STF: z^{-1}
 *
 * @param input   Current analog input sample.
 * @param state   Modulator state (first-order).
 * @param v_ref   Reference voltage.
 * @return Modulator output code (-1, 0, 1 for three-level, etc.).
 *
 * Complexity: O(1)
 */
int64_t sdm_process_first_order(double input, sdm_state_t *state, double v_ref);

/**
 * @brief Process one sample through a second-order ΔΣ modulator.
 *
 * Boser-Wooley structure:
 *   i1[n] = i1[n-1] + input - y[n-1]
 *   i2[n] = i2[n-1] + i1[n]
 *   y[n] = quantize(i2[n])
 *
 * @param input   Current input.
 * @param state   Modulator state (2nd order).
 * @param v_ref   Reference voltage.
 * @return Quantized output code.
 *
 * Complexity: O(1)
 */
int64_t sdm_process_second_order(double input, sdm_state_t *state, double v_ref);

/**
 * @brief Process one sample through a ΔΣ modulator with general NTF.
 *
 * Uses a direct-form implementation of the loop filter:
 *   y[n] = Q( STF*u[n] + NTF*e[n-1] )
 *
 * @param input   Current input.
 * @param state   Modulator state (any order).
 * @param v_ref   Reference voltage.
 * @return Quantized output.
 *
 * Complexity: O(order)
 */
int64_t sdm_process_general(double input, sdm_state_t *state, double v_ref);

/* =========================================================================
 * L5 — MASH (Multi-stAge noise SHaping) Architecture
 * ========================================================================= */

/**
 * @brief Process one sample through a MASH 1-1 modulator.
 *
 * Two cascaded first-order stages with digital error cancellation.
 * Overall NTF = (1 - z^{-1})² (second-order shaping).
 *
 * @param input   Current input.
 * @param state_1 First-stage state.
 * @param state_2 Second-stage state.
 * @param v_ref   Reference voltage.
 * @param[out] mash_output Combined output code.
 *
 * Complexity: O(1)
 */
void sdm_mash_1_1_process(double input,
                           sdm_state_t *state_1, sdm_state_t *state_2,
                           double v_ref, int64_t *mash_output);

/* =========================================================================
 * L5 — Decimation Filter for ΔΣ Output
 * ========================================================================= */

/**
 * @brief Design coefficients for a CIC decimation filter.
 *
 * CIC filter of order N and rate R:
 * H(z) = [(1 - z^{-RM}) / (1 - z^{-1})]^N
 *
 * @param rate       Decimation rate R.
 * @param stages     Number of CIC stages N.
 * @param delay      Differential delay M (1 or 2).
 * @param[out] coeffs_len Number of output coefficients.
 * @return Array of coefficient values (caller must free).
 */
double *sdm_design_cic_filter(uint32_t rate, uint32_t stages,
                               uint32_t delay, uint32_t *coeffs_len);

/**
 * @brief Compensate CIC droop with a FIR filter.
 *
 * CIC filters have a sinc^N response that droops in the passband.
 * A compensation FIR can flatten the passband.
 *
 * @param cic_rate   CIC decimation rate.
 * @param cic_stages CIC stages.
 * @param cic_delay  CIC differential delay.
 * @param bandwidth_hz Signal bandwidth.
 * @param fs_out_hz  Output sample rate (after CIC).
 * @param[out] n_taps Number of compensation filter taps.
 * @return FIR coefficient array (caller must free).
 */
double *sdm_design_cic_compensation(uint32_t cic_rate, uint32_t cic_stages,
                                     uint32_t cic_delay,
                                     double bandwidth_hz, double fs_out_hz,
                                     uint32_t *n_taps);

/* =========================================================================
 * L6 — Complete ΔΣ ADC Simulation
 * ========================================================================= */

/**
 * @brief Simulate a complete ΔΣ ADC: modulator → decimation.
 *
 * Generates an input sine, processes through the ΔΣ modulator,
 * applies decimation filter, and returns baseband samples.
 *
 * @param spec           ΔΣ modulator specification.
 * @param f_in_hz        Input signal frequency.
 * @param amplitude      Input signal amplitude.
 * @param duration_s     Simulation duration in seconds.
 * @param[out] decimated_samples Decimated output (length ≈ duration * 2*BW).
 * @param[out] n_decimated Number of decimated samples.
 * @param[out] osr_samples Raw ΔΣ bitstream (length ≈ duration * fs).
 * @param[out] n_osr_samples Number of oversampled samples.
 * @return 0 on success.
 *
 * Complexity: O(fs * duration)
 */
int sdm_simulate_adc(const sdm_spec_t *spec,
                      double f_in_hz, double amplitude,
                      double duration_s,
                      double *decimated_samples, uint32_t *n_decimated,
                      int64_t *osr_samples, uint64_t *n_osr_samples);

/* =========================================================================
 * L4 — Stability Analysis
 * ========================================================================= */

/**
 * @brief Estimate the maximum stable input amplitude for a ΔΣ modulator.
 *
 * For single-bit quantizers, the modulator overloads when the
 * integrator output exceeds the quantizer input range.
 * Lee's rule: max |NTF| < 1.5 to 2.0 for stability.
 *
 * @param spec     Modulator specification.
 * @return Maximum stable input amplitude relative to V_ref (0 to 1).
 */
double sdm_stable_input_limit(const sdm_spec_t *spec);

/**
 * @brief Check Lee's criterion for NTF stability.
 *
 * max|NTF(e^{jω})| < 2.0 (for single-bit)
 * max|NTF(e^{jω})| < 1.5 (safe design margin)
 *
 * @param ntf_num NTF numerator coefficients.
 * @param ntf_den NTF denominator coefficients.
 * @param order   Filter order.
 * @param n_freq  Number of frequency points to evaluate.
 * @param[out] max_ntf_gain Max |NTF| over frequency.
 * @return 1 if stable, 0 if potentially unstable.
 */
int sdm_check_lee_criterion(const double *ntf_num, const double *ntf_den,
                             uint32_t order, uint32_t n_freq,
                             double *max_ntf_gain);

#ifdef __cplusplus
}
#endif

#endif /* SIGMA_DELTA_H */
