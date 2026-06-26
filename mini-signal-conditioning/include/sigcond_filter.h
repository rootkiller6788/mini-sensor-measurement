/**
 * sigcond_filter.h - Analog & Digital Filter Design
 *
 * Anti-aliasing filters, active filter design (Sallen-Key, MFB, Twin-T),
 * and digital filter coefficient calculation (FIR, IIR).
 *
 * L2: anti-alias, noise bandwidth, notch filtering
 * L3: Laplace/z-domain transfer functions
 * L4: Nyquist-Shannon sampling theorem
 * L5: Butterworth/Chebyshev/Bessel coefficient generation
 *
 * Reference: Zumbahlen (2008), Van Valkenburg (1982), Oppenheim & Schafer (2010)
 */
#ifndef SIGCOND_FILTER_H
#define SIGCOND_FILTER_H
#include "sigcond_defs.h"
#ifdef __cplusplus
extern "C" {
#endif

/** Minimum Butterworth order for anti-alias: n = ceil(log10((10^(As/10)-1)/(10^(Ap/10)-1))/(2*log10(fs/fc))) */
int filter_butterworth_order(double f_cutoff_hz, double f_stop_hz, double passband_ripple_db, double stopband_atten_db);

/** Verify anti-alias filter meets Nyquist: fc <= fs/(2*k_oversample) with k=2.5 for 12-bit */
bool filter_verify_antialias(double f_cutoff_hz, double f_sample_hz, unsigned adc_bits, double snr_target_db);

/** Bandwidth from settling time: BW = 0.35 / t_rise for 1st-order, adjusts for higher orders */
double filter_bandwidth_from_settling(double settling_time_us, unsigned order);

/** Butterworth poles on unit circle: p_k = -sin(theta_k) + j*cos(theta_k), theta_k = (2k-1)pi/(2n) */
unsigned butterworth_poles(unsigned order, double poles_real[], double poles_imag[], double *real_pole);

/** Butterworth biquads: Q = 1/(2*sin(theta/2)) for each conjugate pair */
void butterworth_biquads(unsigned order, double w0_norm[], double q_values[], unsigned *num_sections);

/** Butterworth denominator polynomial: B_n(s) = sum a_k*s^k */
void butterworth_denom_poly(unsigned order, double a_coeffs[]);

/** Chebyshev I poles on s-plane ellipse: p_k = -sinh(eta)*sin(theta_k) + j*cosh(eta)*cos(theta_k) */
unsigned chebyshev1_poles(unsigned order, double ripple_db, double poles_real[], double poles_imag[]);

/** Chebyshev I biquad parameters normalized to wc=1 */
void chebyshev1_biquads(unsigned order, double ripple_db, double w0_norm[], double q_values[], unsigned *num_sections);

/** Chebyshev I denominator polynomial coefficients */
void chebyshev1_denom_poly(unsigned order, double ripple_db, double a_coeffs[]);

/** Bessel poles for maximally flat group delay: B_n(s) = sum b_k*s^k, b_k = (2n-k)!/(2^(n-k)*k!*(n-k)!) */
unsigned bessel_poles(unsigned order, double poles_real[], double poles_imag[], double *real_pole);

/** Bessel biquad parameters: normalized Q values for each 2nd-order section */
void bessel_biquads(unsigned order, double w0_norm[], double q_values[], unsigned *num_sections);

/** Elliptic (Cauer) filter order estimation using complete elliptic integral ratio */
unsigned elliptic_order(double f_cutoff_hz, double f_stop_hz, double passband_ripple_db, double stopband_atten_db);

/** Sallen-Key LP design (equal-R): fc = 1/(2*pi*R*C), Q = 1/(3-G) where G=1+Rf/Rg. Returns Rf/Rg. */
double sallen_key_lp_design(double fc_hz, double q_target, double *r_ohm, double *c_farad, double *rf_by_rg);

/** Sallen-Key component sensitivity: S_R_fc = dfc/fc / dR/R */
void sallen_key_sensitivity(double r1, double r2, double c1, double c2, double sens_r1[], double sens_r2[], double sens_c1[], double sens_c2[]);

/** MFB LP design (equal-C): H(s) = -R5/R1 * 1/(s^2*R2*R5*C^2 + s*C*(R2+R5+R2*R5/R1) + 1) */
double mfb_lp_design(double fc_hz, double q_target, double dc_gain, double *r1, double *r2, double *r5, double *c3, double *c4);

/** Twin-T notch: H(s) = (s^2+w0^2)/(s^2+w0*s*(4-2*beta)+w0^2), w0=1/(RC) */
void twint_notch_design(double notch_freq_hz, double q_target, double *r_ohm, double *c_farad, double *beta_feedback);

/** Twin-T notch depth with component tolerances */
double twint_notch_depth(double nominal_r, double nominal_c, double tolerance_pct_rc);

/** Moving average FIR: y[n] = (1/N)*sum(x[n-k]), H(z) = (1/N)*(1-z^{-N})/(1-z^{-1}) */
void moving_average_coeffs(unsigned window_size, double coeffs[], unsigned *length);

/** First-order IIR LP alpha: y[n] = alpha*x[n] + (1-alpha)*y[n-1], alpha = 1-exp(-2*pi*fc/fs) */
double iir1_lp_alpha(double fc_hz, double fs_hz);

/** Bilinear transform: s = (2/T)*(1-z^{-1})/(1+z^{-1}), maps analog biquad to digital */
void bilinear_transform_sos(double fc_hz, double fs_hz, double a1_analog, double a0_analog, double b2, double b1, double b0, double *B0, double *B1, double *B2, double *A1, double *A2);

/** Kaiser FIR length estimate: N = (A-7.95)/(14.36*delta_f/fs), A=stopband_atten_dB */
unsigned kaiser_fir_length(double transition_width_hz, double fs_hz, double stopband_atten_db);

/** Kaiser window: w[n] = I0(beta*sqrt(1-(2n/(N-1)-1)^2))/I0(beta) */
void kaiser_window(unsigned length, double beta, double window[]);

/** FIR LP via window method: h[n] = sinc(2*fc/fs*(n-(N-1)/2)) * window[n] */
void fir_lp_window_design(double fc_hz, double fs_hz, unsigned length, const double window[], double fir_coeffs[]);

/** Group delay of linear-phase FIR: tau_g = (N-1)/2 samples */
double fir_group_delay_samples(unsigned length);

#ifdef __cplusplus
}
#endif
#endif /* SIGCOND_FILTER_H */
