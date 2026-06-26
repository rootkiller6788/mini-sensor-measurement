# Knowledge Graph — mini-measurement-error

## L1: Definitions — Complete
- ErrorCategory enum (16 values): systematic bias, random noise, hysteresis, nonlinearity, quantization, zero/span drift, loading, environmental, aging, crosstalk, common mode, ground loop, thermal EMF
- MeasurementType enum (6 values): direct, indirect, differential, ratiometric, null-balance, substitution
- MeasurementPoint struct: value, true_value, systematic_error, random_error, uncertainty, timestamp, category_mask
- MeasurementStats struct: mean, std_dev, variance, count, min/max, range, bias, precision, rms_error
- ErrorBudget struct: components[], total_systematic, total_random_rss, combined_uncertainty, expanded_uncertainty
- SensorSpec struct: full_scale_range, sensitivity, offset, nonlinearity_pct, hysteresis_pct, repeatability_pct, resolution, noise_rms, temp_coeff_zero/span, long_term_stability, bandwidth, input_impedance

## L2: Core Concepts — Complete
- Systematic vs random error classification
- Error computation functions (absolute, relative, percentage, full-scale)
- Accuracy class (IEC 60051): 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.5, 5.0
- Repeatability (ISO 5725-2): r = t_95 * sqrt(2) * s_r
- Reproducibility: R = t_95 * sqrt(2) * s_R (includes inter-lab variance)
- Hysteresis error computation
- Endpoint nonlinearity (maximum deviation from endpoint line)
- ADC quantization error (uniform distribution RMS = Q/sqrt(12))
- Loading error (voltage divider effect)
- Common-mode rejection error (CMRR-based)
- Thermal EMF error (Seebeck effect)
- Error budget construction and RSS combination

## L3: Mathematical Structures — Complete
- UncertaintyDistribution enum (6 types): normal, uniform, triangular, U-shaped, Student-t, trapezoid
- UncertaintyEvalType: Type A (statistical), Type B (other)
- UncertaintyComponent struct: name, eval_type, distribution, half_width, divisor, sensitivity_coeff, degrees_freedom, standard_uncertainty, contribution
- Gaussian PDF/CDF/quantile (Abramowitz & Stegun rational approximation)
- Uniform PDF
- Confidence interval computation (Student's t)
- Prediction interval
- Shapiro-Wilk W statistic (normality test)
- Kolmogorov-Smirnov test statistic
- t-test (one-sample, two-sample Welch)
- F-test for variance equality
- Chi-squared variance test
- Required sample size computation

## L4: Fundamental Laws — Complete
- GUM uncertainty framework (JCGM 100:2008): u_c = sqrt(SUM (c_i*u_i)^2)
- Welch-Satterthwaite equation for effective degrees of freedom
- Student's t rational approximation
- Error propagation law (first-order Taylor): u_c^2 = SUM (df/dX_i)^2 * u_i^2
- Propagation with correlation (covariance matrix)
- Common formulas: sum, product, quotient, power, exponential, logarithm
- Expanded uncertainty: U = k * u_c
- MA noise reduction: sigma_out = sigma_in / sqrt(N)
- Inverse-variance optimality (Gauss-Markov)
- Central Limit Theorem for measurement averages

## L5: Algorithms/Methods — Complete
- Linear least squares (direct 2x2 normal equations)
- Polynomial regression (Gaussian elimination with partial pivoting, up to 4th order)
- Weighted least squares (diagonal weight matrix)
- Exponential/power-law/logarithmic fit (log-space linearization)
- Calibration model evaluation and inverse (Newton-Raphson)
- Fit statistics: R-squared, RMSE, max residual, lack-of-fit F-test
- Prediction uncertainty (ISO 11095)
- Moving average filter (O(1) ring buffer)
- Exponential moving average (1st-order IIR)
- Median filter (quickselect O(n))
- Running median filter (sliding window)
- 1D Kalman filter (predict-update, Joseph form)
- Sensor fusion: weighted average, inverse-variance
- Gaussian MLE parameter estimation
- Grubbs outlier test (ISO 5725-2)
- Modified Z-score outlier detection (robust, median/MAD-based)

## L6: Canonical Problems — Complete
- Thermocouple cold-junction compensation
- Strain gauge temperature compensation (CTE mismatch)
- Wheatstone bridge nonlinearity correction (quarter/half bridge)
- Lead wire resistance compensation (3-wire, 4-wire Kelvin)
- Polynomial nonlinearity correction
- Piecewise linear LUT correction (binary search interpolation)
- Full compensation pipeline (zero → nonlinearity → temperature → gain)

## L7: Applications — Complete
- Thermocouple system error budget (Seebeck, CJC, amplifier, ADC, reference, EMI)
- Strain gauge bridge uncertainty budget (excitation, GF, imbalance, amplifier, ADC)
- ADC/DAQ error analysis: ENOB, SINAD, total noise, jitter SNR
- pH sensor: Nernst equation, temperature-compensated slope, drift budget

## L8: Advanced Topics — Complete
- Monte Carlo uncertainty propagation (GUM Supplement 1): Box-Muller, LCG, percentiles
- Bayesian sensor fault detection: Bayes rule update, prior/posterior probabilities
- Time-varying drift model (Wiener process): drift + diffusion, first-passage probability
- Allan variance (IEEE Std 1139): multi-tau analysis for noise characterization

## L9: Research Frontiers — Partial
- Heisenberg uncertainty principle (quantum measurement limit): documented
