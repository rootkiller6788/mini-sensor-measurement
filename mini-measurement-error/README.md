# mini-measurement-error — Measurement Error Analysis

## Module Status: COMPLETE &#x2705;

| Level | Name | Status |
|-------|------|--------|
| L1 | Definitions | **Complete** — 16 error categories, 6 measurement types, 5 core structs, SensorSpec |
| L2 | Core Concepts | **Complete** — 10 error computation functions, systematic vs random, hysteresis, loading |
| L3 | Math Structures | **Complete** — 6 distribution types, Gaussian PDF/CDF/quantile, confidence intervals |
| L4 | Fundamental Laws | **Complete** — GUM uncertainty framework, Welch-Satterthwaite, Student t, error propagation |
| L5 | Algorithms/Methods | **Complete** — OLS/WLS/polynomial regression, Kalman filter, moving avg, EMA, median filter |
| L6 | Canonical Problems | **Complete** — Thermocouple CJC, strain gauge bridge, ADC ENOB, pH calibration |
| L7 | Applications | **Complete** — Thermocouple system, strain gauge, ADC/DAQ, pH sensor error budgets |
| L8 | Advanced Topics | **Complete** — Monte Carlo uncertainty propagation, Bayesian fault detection, Wiener drift model |
| L9 | Research Frontiers | **Partial** — Heisenberg bound (documented) |

## Core Definitions (L1)

- **ErrorCategory**: 16 distinct error mechanisms (systematic bias, random noise, hysteresis, nonlinearity, quantization, zero/span drift, loading, environmental, aging, crosstalk, CM error, ground loop, thermal EMF)
- **MeasurementType**: Direct, indirect, differential, ratiometric, null-balance, substitution
- **MeasurementPoint**: Single measurement with full error attribution and timestamp
- **MeasurementStats**: Descriptive statistics with bias, precision, RMS error
- **ErrorBudget**: Multi-source error aggregation with RSS combination
- **SensorSpec**: Standard datasheet parameters (FSR, sensitivity, nonlinearity, noise, drift, etc.)

## Core Theorems (L4)

- **GUM Uncertainty Framework** (JCGM 100:2008): u_c = sqrt(SUM (c_i * u_i)^2)
- **Welch-Satterthwaite Equation**: nu_eff = u_c^4 / SUM((c_i*u_i)^4/nu_i)
- **Error Propagation Law**: u_c^2(y) = SUM (df/dX_i)^2 * u^2(x_i) (first-order Taylor)
- **Central Limit Theorem**: sigma_mean = sigma / sqrt(N) for measurement averages
- **Inverse-Variance Optimality**: Gauss-Markov theorem for sensor fusion
- **Moving Average Noise Reduction**: sigma_out = sigma_in / sqrt(N)
- **Heisenberg Uncertainty Principle**: dx * dp >= hbar/2 (L9 documented)

## Core Algorithms (L5)

- **Linear Least Squares**: O(N) direct solution of 2x2 normal equations
- **Polynomial Regression** (up to 4th order): Vandermonde normal equations with Gaussian elimination and partial pivoting
- **Weighted Least Squares**: Diagonal weight matrix for heteroscedastic data
- **Exponential/Power-law/Logarithmic Fit**: Log-space linearization with back-transformation
- **1D Kalman Filter**: Predict-update cycle with Joseph-form covariance
- **Moving Average**: O(1) ring buffer with running sum
- **EMA Filter**: 1st-order IIR, time-constant to alpha conversion
- **Median Filter**: Quickselect O(n) implementation
- **Inverse-Variance Sensor Fusion**: Optimal linear combination
- **Gaussian Elimination**: Partial pivoting, back substitution, O(n^3)
- **Monte Carlo MCM**: GUM Supplement 1, Box-Muller for normal samples
- **Grubbs/Modified Z-Score**: Outlier detection (ISO 5725-2)

## Canonical Problems (L6)

- **Thermocouple CJC**: Cold-junction compensation using Seebeck coefficient
- **Strain Gauge Bridge**: Quarter/half-bridge nonlinearity correction, lead wire compensation
- **ADC Analysis**: ENOB from SINAD, jitter SNR, INL/DNL noise budget
- **pH Sensor**: Nernst equation, temperature-compensated slope, drift modeling
- **Temperature Compensation**: Zero + span + quadratic drift model
- **Kelvin 4-Wire**: Lead resistance elimination for RTD
- **Piecewise LUT Correction**: Binary search interpolation for arbitrary nonlinearity

## Course Mapping

| School | Course | Key Topics |
|--------|--------|------------|
| MIT | 6.630 EM Waves & 6.450 Digital Comm | Sensor measurements, noise, SNR |
| Stanford | EE359 Wireless | Sensor fusion, estimation theory |
| Berkeley | EE117 Electromagnetics / EE123 DSP | Instrumentation, filtering |
| Illinois | ECE 310 DSP / ECE 451 EM | Measurement uncertainty |
| Michigan | EECS 351 DSP / EECS 411 Microwave | Sensor characterization |
| Georgia Tech | ECE 4270 DSP / ECE 6350 EM | Signal processing for sensors |
| TU Munich | Signal Processing / High-Frequency Eng | Messtechnik (measurement technology) |
| ETH | 227-0427 Signal Processing / 227-0455 EM | Precision measurements |
| Tsinghua | 信号与系统 / 电磁场 | Sensor principles |

## File Structure

```
mini-measurement-error/
├── Makefile
├── README.md
├── include/              (7 headers, 999 lines)
│   ├── measurement_error.h           — L1-L2: Types, categories, core functions
│   ├── measurement_uncertainty.h     — L3-L4: GUM, distributions, propagation
│   ├── measurement_calibration.h     — L2,L5: Regression, curve fitting
│   ├── measurement_statistics.h      — L3,L8: Hypothesis tests, Allan variance
│   ├── measurement_filtering.h       — L5: MA, EMA, Kalman, sensor fusion
│   ├── measurement_compensation.h    — L6: Temp comp, nonlinearity correction
│   └── measurement_applications.h    — L7-L8: TC, SG, ADC, pH, MC, Bayesian
├── src/                  (7 C files + 1 Lean, 3154 lines)
│   ├── measurement_error.c           — 329 lines
│   ├── measurement_uncertainty.c     — 322 lines
│   ├── measurement_calibration.c     — 579 lines
│   ├── measurement_statistics.c      — 455 lines
│   ├── measurement_filtering.c       — 300 lines
│   ├── measurement_compensation.c    — 260 lines
│   ├── measurement_applications.c    — 645 lines
│   └── measurement_error.lean        — 266 lines
├── tests/                (3 test suites, 44 tests total)
│   ├── test_error.c                  — 16 tests
│   ├── test_uncertainty.c            — 12 tests
│   └── test_calibration.c            — 16 tests
├── examples/             (4 end-to-end examples)
│   ├── example_thermocouple.c        — Thermocouple system error budget
│   ├── example_strain_gauge.c        — Strain gauge uncertainty analysis
│   ├── example_adc.c                 — ADC ENOB + Kalman + temp comp
│   └── example_sensor_fusion.c       — Multi-sensor fusion + outlier detection
├── demos/
├── benches/
└── docs/
```

## Build & Test

```bash
make          # Compile all source files
make test     # Build and run all tests (44/44 pass)
make examples # Build and run all examples
make clean    # Remove build artifacts
```

**include/ + src/ total: 4153 lines** &#x2705;

## References

- JCGM 100:2008 — Guide to the Expression of Uncertainty in Measurement (GUM)
- JCGM 101:2008 — GUM Supplement 1: Monte Carlo Method
- ISO 5725-1:1994 — Accuracy (trueness and precision)
- ISO 5725-2:2019 — Repeatability and reproducibility
- ISO 11095:1996 — Linear calibration using reference materials
- IEC 60051-1:2017 — Direct acting indicating analogue electrical measuring instruments
- IEEE Std 1139-2008 — Frequency stability (Allan variance)
- NIST/SEMATECH e-Handbook of Statistical Methods
- Abramowitz & Stegun — Handbook of Mathematical Functions (1964)
